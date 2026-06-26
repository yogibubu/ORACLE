from __future__ import annotations

import argparse
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
