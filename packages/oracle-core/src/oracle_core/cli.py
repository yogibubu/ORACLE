from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def find_repo_root(start: Path | None = None) -> Path:
    env_root = os.environ.get("ORACLE_HOME")
    if env_root:
        return Path(env_root).expanduser().resolve()

    search_from = Path.cwd() if start is None else Path(start).resolve()
    for candidate in (search_from, *search_from.parents):
        if (candidate / "packages").is_dir() and (candidate / "pyproject.toml").is_file():
            return candidate
    return Path.cwd().resolve()


def add_repo_packages_to_path(repo_root: Path | None = None) -> None:
    root = find_repo_root() if repo_root is None else Path(repo_root)
    packages = root / "packages"
    if not packages.is_dir():
        return
    for src in sorted(packages.glob("*/src")):
        text = str(src)
        if text not in sys.path:
            sys.path.insert(0, text)


def build_parser(*, repo_root: Path | None = None) -> argparse.ArgumentParser:
    root = find_repo_root() if repo_root is None else Path(repo_root)
    parser = argparse.ArgumentParser(prog="oracle", description="ORACLE workflow CLI")
    sub = parser.add_subparsers(dest="command")

    init = sub.add_parser("init", help="Create an ORACLE project workspace")
    init.add_argument("workdir", type=Path)

    validate = sub.add_parser("validate", help="Validate an enriched XYZ after preprocessing")
    validate.add_argument("xyzin", type=Path)
    validate.add_argument("--require-fragments", action="store_true")

    babel = sub.add_parser("babel", help="Run ORACLE-Babel preprocessing")
    babel_sub = babel.add_subparsers(dest="babel_command")
    preprocess = babel_sub.add_parser("preprocess", help="Import a source into enriched XYZ")
    preprocess.add_argument("source", type=Path)
    preprocess.add_argument("output", type=Path)
    preprocess.add_argument(
        "--source-kind",
        choices=("auto", "xyz", "enriched_xyz"),
        default="auto",
    )
    preprocess.add_argument("--symmetry-distance", type=float, default=1.0e-3)
    preprocess.add_argument("--symmetry-inertia", type=float, default=1.0e-3)
    preprocess.add_argument("--max-rotation-order", type=int, default=6)

    lcb25 = sub.add_parser("lcb25", help="Manage the local ORACLE LCB25 geometry cache")
    lcb25_sub = lcb25.add_subparsers(dest="lcb25_command")
    fetch = lcb25_sub.add_parser("fetch", help="Download/extract LCB25 geometries once")
    fetch.add_argument("--root", type=Path, default=root / "data" / "lcb25")
    fetch.add_argument("--dataset", action="append", help="PCS2, SE or HPCS2; repeatable")
    fetch.add_argument("--force", action="store_true")

    fragments = sub.add_parser("fragments", help="Manage topology-backed fragment workflows")
    fragments_sub = fragments.add_subparsers(dest="fragments_command")
    plan = fragments_sub.add_parser("plan", help="Write the initial #FRAGMENTS section")
    plan.add_argument("xyzin", type=Path)

    rovib = sub.add_parser("rovib", help="Standalone rovibrational xyzin utilities")
    rovib_sub = rovib.add_subparsers(dest="rovib_command")
    rovib_summary = rovib_sub.add_parser("summarize", help="Summarize rovib sections")
    rovib_summary.add_argument("xyzin", type=Path)

    gicforge = sub.add_parser("gicforge", help="Plan GICForge post-validation sections")
    gicforge_sub = gicforge.add_subparsers(dest="gicforge_command")
    gic_plan = gicforge_sub.add_parser("plan", help="Write planned #GIC/#SYCART sections")
    gic_plan.add_argument("xyzin", type=Path)
    gic_plan.add_argument("--symmetrize", action="store_true")
    gic_plan.add_argument("--sycart", action="store_true")
    gic_build = gicforge_sub.add_parser("build", help="Build frozen #GIC/#SYCART sections")
    gic_build.add_argument("xyzin", type=Path)
    gic_build.add_argument("--symmetrize", action="store_true")
    gic_build.add_argument("--sycart", action="store_true")
    corpus = gicforge_sub.add_parser(
        "corpus",
        help="List or summarize the demanding GIC regression corpus",
    )
    corpus.add_argument("--root", type=Path, help="Override the GIC corpus root directory")
    corpus.add_argument("--suffix", action="append", help="Filter by suffix, for example .inp or fchk")
    corpus.add_argument("--limit", type=int, help="Limit listed records")
    corpus.add_argument(
        "--format",
        choices=("summary", "paths", "json"),
        default="summary",
        help="Output format",
    )
    corpus_audit = gicforge_sub.add_parser(
        "corpus-audit",
        help="Audit geometry imports for the GIC regression corpus",
    )
    corpus_audit.add_argument("--root", type=Path, help="Override the GIC corpus root directory")
    corpus_audit.add_argument("--suffix", action="append", help="Filter by suffix; defaults to geometry inputs")
    corpus_audit.add_argument("--limit", type=int, help="Limit audited or listed records")
    corpus_audit.add_argument(
        "--format",
        choices=("summary", "failures", "json"),
        default="summary",
        help="Output format",
    )
    corpus_audit.add_argument(
        "--status",
        choices=("all", "pass", "fail"),
        default="all",
        help="Entry status filter for JSON output",
    )
    gaussian_input = gicforge_sub.add_parser(
        "gaussian-input",
        help="Write Gaussian input from validated #GIC state",
    )
    gaussian_input.add_argument("xyzin", type=Path)
    gaussian_input.add_argument("output", type=Path)
    gaussian_input.add_argument("--route", default="#p hf/sto-3g")
    gaussian_input.add_argument("--title")
    gaussian_input.add_argument("--charge", type=int)
    gaussian_input.add_argument("--multiplicity", type=int)

    return parser


def main(argv: list[str] | None = None, *, repo_root: Path | None = None) -> int:
    root = find_repo_root() if repo_root is None else Path(repo_root)
    add_repo_packages_to_path(root)
    parser = build_parser(repo_root=root)
    args = parser.parse_args(argv)
    if args.command == "init":
        from oracle_core.workspace import ensure_workspace

        ensure_workspace(args.workdir)
        print(f"Created ORACLE workspace: {args.workdir}")
        return 0
    if args.command == "validate":
        from oracle_chem import write_validation_section

        result = write_validation_section(args.xyzin, require_fragments=args.require_fragments)
        print(f"Validated ORACLE molecule: {args.xyzin} ({result.status})")
        return 0
    if args.command == "babel" and args.babel_command == "preprocess":
        from oracle_chem import SymmetryThresholds, preprocess_to_enriched_xyz

        result = preprocess_to_enriched_xyz(
            args.source,
            args.output,
            source_kind=args.source_kind,
            symmetry_thresholds=SymmetryThresholds(
                distance_angstrom=args.symmetry_distance,
                inertia_relative=args.symmetry_inertia,
                max_rotation_order=args.max_rotation_order,
            ),
        )
        print(
            "Preprocessed ORACLE-Babel molecule: "
            f"{result.path} ({result.geometry.natoms} atoms, "
            f"PG={result.point_group}, bonds={result.topology_bond_count}, "
            f"rings={result.ring_count})"
        )
        return 0
    if args.command == "lcb25" and args.lcb25_command == "fetch":
        from oracle_babel import sync_lcb25_library

        manifest = sync_lcb25_library(args.root, datasets=args.dataset, force=args.force)
        print(f"Synced LCB25 library: {manifest}")
        return 0
    if args.command == "fragments" and args.fragments_command == "plan":
        from oracle_fragments import write_fragment_plan_section

        write_fragment_plan_section(args.xyzin)
        print(f"Planned ORACLE fragment workflow: {args.xyzin}")
        return 0
    if args.command == "rovib" and args.rovib_command == "summarize":
        from oracle_rovib import rovib_summary_lines, summarize_xyzin

        print("\n".join(rovib_summary_lines(summarize_xyzin(args.xyzin))))
        return 0
    if args.command == "gicforge" and args.gicforge_command == "plan":
        from oracle_gicforge import write_gicforge_plan_sections

        write_gicforge_plan_sections(
            args.xyzin,
            symmetrize=args.symmetrize,
            sycart=args.sycart,
        )
        print(f"Planned GICForge workflow: {args.xyzin}")
        return 0
    if args.command == "gicforge" and args.gicforge_command == "build":
        from oracle_gicforge import write_gicforge_build_sections

        definition = write_gicforge_build_sections(
            args.xyzin,
            symmetrize=args.symmetrize,
            sycart=args.sycart,
        )
        print(
            "Built GICForge definition: "
            f"{args.xyzin} (GICs={len(definition.gics)}, rank={definition.rank})"
        )
        return 0
    if args.command == "gicforge" and args.gicforge_command == "corpus":
        from oracle_gicforge import (
            default_gic_corpus_root,
            format_gic_corpus_paths,
            format_gic_corpus_summary,
            gic_corpus_records,
            summarize_gic_corpus,
        )

        corpus_root = args.root or default_gic_corpus_root(root)
        summary = summarize_gic_corpus(corpus_root, suffixes=args.suffix)
        if args.format == "json":
            payload = {
                "root": str(summary.root),
                "total_files": summary.total_files,
                "suffix_counts": summary.suffix_counts,
                "role_counts": summary.role_counts,
                "entries": gic_corpus_records(summary, limit=args.limit),
            }
            print(json.dumps(payload, indent=2, sort_keys=True))
            return 0
        if args.format == "paths":
            print("\n".join(format_gic_corpus_paths(summary, limit=args.limit)))
            return 0
        print("\n".join(format_gic_corpus_summary(summary)))
        return 0
    if args.command == "gicforge" and args.gicforge_command == "corpus-audit":
        from oracle_gicforge import (
            audit_gic_corpus_geometry,
            default_gic_corpus_root,
            format_gic_corpus_geometry_audit_summary,
            format_gic_corpus_geometry_failures,
            gic_corpus_geometry_audit_records,
        )

        corpus_root = args.root or default_gic_corpus_root(root)
        audit = audit_gic_corpus_geometry(
            corpus_root,
            suffixes=args.suffix,
            limit=args.limit if args.format == "summary" else None,
        )
        if args.format == "json":
            payload = {
                "root": str(audit.root),
                "total_files": audit.total_files,
                "passed_files": audit.passed_files,
                "failed_files": audit.failed_files,
                "source_format_counts": audit.source_format_counts,
                "error_counts": audit.error_counts,
                "entries": gic_corpus_geometry_audit_records(
                    audit,
                    status=args.status,
                    limit=args.limit,
                ),
            }
            print(json.dumps(payload, indent=2, sort_keys=True))
            return 0
        if args.format == "failures":
            print("\n".join(format_gic_corpus_geometry_failures(audit, limit=args.limit)))
            return 0
        print("\n".join(format_gic_corpus_geometry_audit_summary(audit)))
        return 0
    if args.command == "gicforge" and args.gicforge_command == "gaussian-input":
        from oracle_gicforge import write_gicforge_gaussian_input

        output = write_gicforge_gaussian_input(
            args.xyzin,
            args.output,
            route=args.route,
            title=args.title,
            charge=args.charge,
            multiplicity=args.multiplicity,
        )
        print(f"Wrote Gaussian input: {output}")
        return 0
    parser.print_help()
    return 0
