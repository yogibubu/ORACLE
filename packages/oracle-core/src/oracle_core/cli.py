from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import TypeVar


_T = TypeVar("_T")


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

    gaussian = sub.add_parser("gaussian", help="Gaussian adapter and job utilities")
    gaussian_sub = gaussian.add_subparsers(dest="gaussian_command")
    gaussian_summary = gaussian_sub.add_parser("summary", help="Summarize a Gaussian log/out")
    gaussian_summary.add_argument("log", type=Path)
    gaussian_status = gaussian_sub.add_parser("status", help="Inspect Gaussian job state")
    gaussian_status.add_argument("workdir", type=Path)
    gaussian_run = gaussian_sub.add_parser("run", help="Run Gaussian from a work directory")
    gaussian_run.add_argument("workdir", type=Path)
    gaussian_run.add_argument("--executable")
    gaussian_run.add_argument("--input", type=Path)
    gaussian_run.add_argument("--background", action="store_true")
    gaussian_run.add_argument("--timeout", type=float)
    gaussian_formchk = gaussian_sub.add_parser("formchk", help="Run formchk on a checkpoint file")
    gaussian_formchk.add_argument("chk", type=Path)
    gaussian_formchk.add_argument("fchk", type=Path, nargs="?")
    gaussian_formchk.add_argument("--executable", default=None)
    gaussian_formchk.add_argument("--timeout", type=float)
    gaussian_fchk = gaussian_sub.add_parser("fchk-summary", help="Summarize FCHK/QFF blocks")
    gaussian_fchk.add_argument("fchk", type=Path)
    gaussian_promote_rovib = gaussian_sub.add_parser(
        "promote-rovib",
        help="Promote Gaussian rovibrational log data into an ORACLE xyzin",
    )
    gaussian_promote_rovib.add_argument("log", type=Path)
    gaussian_promote_rovib.add_argument("xyzin", type=Path)
    gaussian_promote_rovib.add_argument("--no-vibrational", action="store_true")
    gaussian_promote_rovib.add_argument("--no-rotational", action="store_true")
    gaussian_promote_rovib.add_argument("--no-deltabvib", action="store_true")
    gaussian_promote_rovib.add_argument("--no-invert-imaginary", action="store_true")
    gaussian_promote_rovib.add_argument(
        "--exclude-mode",
        type=int,
        action="append",
        default=[],
        help="Exclude a normal-mode index from alpha-derived DeltaBvib",
    )

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
    fragments_build = fragments_sub.add_parser("build", help="Build concrete #FRAGMENTS")
    fragments_build.add_argument("xyzin", type=Path)
    centers = fragments_sub.add_parser(
        "centers",
        help="Build virtual bond/ring interaction centers for GICForge",
    )
    centers.add_argument("xyzin", type=Path)

    rovib = sub.add_parser("rovib", help="Standalone rovibrational xyzin utilities")
    rovib_sub = rovib.add_subparsers(dest="rovib_command")
    rovib_summary = rovib_sub.add_parser("summarize", help="Summarize rovib sections")
    rovib_summary.add_argument("xyzin", type=Path)
    rovib_vibin = rovib_sub.add_parser("vibin", help="Build Merlino-compatible vibin from FCHK")
    rovib_vibin.add_argument("xyzin", type=Path)
    rovib_vibin.add_argument("--fchk", type=Path, required=True)
    rovib_vibin.add_argument("--workdir", type=Path)
    rovib_vibin.add_argument("--no-project-tr", action="store_true")
    rovib_vibin.add_argument("--no-update-vibrational", action="store_true")
    rovib_coriolis = rovib_sub.add_parser("coriolis", help="Compute sparse Coriolis terms")
    rovib_coriolis.add_argument("xyzin", type=Path)
    rovib_coriolis.add_argument("--vibin", type=Path)
    rovib_coriolis.add_argument("--threshold-cm1", type=float, default=1.0)
    rovib_coriolis.add_argument("--all-pairs", action="store_true")
    rovib_coriolis.add_argument("--append-vibin", action="store_true")
    rovib_coriolis.add_argument("--out", type=Path)
    rovib_qcent = rovib_sub.add_parser("qcent", help="Compute quartic centrifugal distortion")
    rovib_qcent.add_argument("xyzin", type=Path)
    rovib_qcent.add_argument("--vibin", type=Path)
    rovib_qcent.add_argument("--append-vibin", action="store_true")
    rovib_qcent.add_argument("--out", type=Path)
    rovib_dos = rovib_sub.add_parser("dos", help="Build direct vibrational DOS from #VIBRATIONAL")
    rovib_dos.add_argument("xyzin", type=Path)
    rovib_dos.add_argument("--vmax", type=int, default=6)
    rovib_dos.add_argument("--emax", type=float, default=8000.0)
    rovib_dos.add_argument("--bin-cm1", type=float, default=50.0)
    rovib_dos.add_argument("--ncap", type=float, default=10.0)
    rovib_dos.add_argument("--out", type=Path)
    rovib_dos_rovib = rovib_sub.add_parser("dos-rovib", help="Convolve vibrational and rotational DOS")
    rovib_dos_rovib.add_argument("xyzin", type=Path)
    rovib_dos_rovib.add_argument("--vib-dos", type=Path)
    rovib_dos_rovib.add_argument("--out", type=Path)
    rovib_dos_rovib.add_argument("--rot-out", type=Path)
    rovib_dos_rovib.add_argument("--q-out", type=Path)
    rovib_dos_rovib.add_argument("--emax-rot", type=float)
    rovib_dos_rovib.add_argument("--jmax", type=int)

    thermo = sub.add_parser("thermo", help="Run thermochemistry from an ORACLE xyzin")
    thermo.add_argument("xyzin", type=Path)
    thermo.add_argument("--out", type=Path, help="Write the readable thermo report here")
    thermo.add_argument("--no-report", action="store_true", help="Do not write thermo.report")
    thermo.add_argument("--no-write-section", action="store_true", help="Do not update #THERMO")
    thermo.add_argument("--cutoff-cm1", type=float, default=10.0)
    thermo.add_argument("--keep-low-positive", action="store_true")

    gf = sub.add_parser("gf", help="Run GF/PED analysis from a Cartesian Hessian")
    gf.add_argument("--fchk", type=Path, required=True)
    gf.add_argument("--xyzin", type=Path, help="Frozen ORACLE xyzin with a BUILT #GIC section")
    gf.add_argument("--out", type=Path, help="Write the GF/PED report")
    gf.add_argument("--csv-dir", type=Path, help="Write GF/PED CSV tables")
    gf.add_argument("--scale-file", type=Path)
    gf.add_argument("--scale", action="append", default=[])
    gf.add_argument("--local", action="store_true", help="Apply local force-field filtering")
    gf.add_argument("--symmetry-blocks", action="store_true", help="Solve separated irrep blocks")
    gf.add_argument("--force-threshold", type=float, help="Zero internal force constants below threshold")
    gf.add_argument(
        "--subtract-electrostatic",
        action="store_true",
        help="Subtract CM5/synthon electrostatic Hessian terms before GF",
    )
    gf.add_argument(
        "--subtract-uff-vdw",
        action="store_true",
        help="Subtract UFF van der Waals Hessian terms before GF",
    )
    gf.add_argument(
        "--nonbonded-14-scale",
        type=float,
        default=0.5,
        help="Scale factor for 1-4 electrostatic and UFF-vdW terms",
    )

    semiexp = sub.add_parser(
        "semiexp",
        help="Fit semiexperimental equilibrium geometry with ORACLE-MORPHEUS",
    )
    semiexp.add_argument(
        "--job",
        type=Path,
        help="ORACLE/Merlino semiexperimental job file or legacy MSR file",
    )
    semiexp.add_argument(
        "--xyz",
        "--geometry",
        dest="xyz",
        type=Path,
        help="Initial parent Cartesian geometry in XYZ or Gaussian .com/.gjf format",
    )
    semiexp.add_argument(
        "--observations",
        type=Path,
        help="CSV/JSON/TOML with isotopologue B0 constants and corrections, or legacy MSR file",
    )
    semiexp.add_argument(
        "--xyzin",
        type=Path,
        help="Canonical ORACLE xyzin container to create/update before SEfit",
    )
    semiexp.add_argument("--outdir", type=Path, required=True)
    semiexp.add_argument("--backend", choices=("python", "fortran77"), default="python")
    semiexp.add_argument(
        "--fixed",
        default="",
        help="Comma/semicolon-separated fixed GIC patterns or Gaussian-style constraints",
    )
    semiexp.add_argument("--fix-hydrogens", action="store_true")
    semiexp.add_argument("--max-iter", type=int, default=None)
    semiexp.add_argument("--step", type=float, default=1.0e-4)
    semiexp.add_argument("--damping", type=float, default=1.0e-8)
    semiexp.add_argument("--max-step", type=float, default=0.25)
    semiexp.add_argument("--prune-condition", type=float, default=0.0)
    semiexp.add_argument(
        "--robust-loss",
        choices=("none", "huber", "soft_l1", "cauchy"),
        default="none",
    )
    semiexp.add_argument("--robust-scale", type=float, default=0.0)
    semiexp.add_argument("--leave-one-out", action="store_true")
    semiexp.add_argument("--checkpoint", type=Path, default=None)
    semiexp.add_argument("--restart", type=Path, default=None)
    semiexp.add_argument(
        "--observable",
        choices=("moments", "rotational_constants", "auto"),
        default="moments",
    )
    semiexp.add_argument(
        "--coordinate-model",
        choices=("gic", "cartesian_symmetry"),
        default="gic",
    )
    semiexp.add_argument(
        "--rotational-components",
        choices=("auto", "ABC", "AB", "AC", "BC"),
        default="auto",
    )
    semiexp.add_argument(
        "--qm-predicate",
        action="append",
        default=[],
        help="QM prior as label_pattern:value:sigma[:source]; can be repeated",
    )
    semiexp.add_argument(
        "--parameter-class",
        action="append",
        default=[],
        help="Class constraint as name:shared|fixed:pattern[|pattern...]; can be repeated",
    )

    semiexp_ensemble = sub.add_parser(
        "semiexp-ensemble",
        help="Fit shared class corrections across multiple semiexperimental molecule jobs",
    )
    semiexp_ensemble.add_argument("--job", type=Path, required=True)
    semiexp_ensemble.add_argument("--outdir", type=Path, required=True)

    gicforge = sub.add_parser("gicforge", help="Plan GICForge post-validation sections")
    gicforge_sub = gicforge.add_subparsers(dest="gicforge_command")
    gic_plan = gicforge_sub.add_parser("plan", help="Write planned #GIC/#SYCART sections")
    gic_plan.add_argument("xyzin", type=Path)
    gic_plan.add_argument("--symmetrize", action="store_true")
    gic_plan.add_argument("--sycart", action="store_true")
    gic_plan.add_argument("--improper-dihedrals", action="store_true")
    gic_build = gicforge_sub.add_parser("build", help="Build frozen #GIC/#SYCART sections")
    gic_build.add_argument("xyzin", type=Path)
    gic_build.add_argument("--symmetrize", action="store_true")
    gic_build.add_argument("--sycart", action="store_true")
    gic_build.add_argument("--improper-dihedrals", action="store_true")
    bmatrix = gicforge_sub.add_parser("bmatrix", help="Evaluate the frozen GIC B matrix")
    bmatrix.add_argument("xyzin", type=Path)
    bmatrix.add_argument("output", type=Path, nargs="?")
    report = gicforge_sub.add_parser("report", help="Write a readable frozen-GIC report")
    report.add_argument("xyzin", type=Path)
    report.add_argument("output", type=Path, nargs="?")
    corpus = gicforge_sub.add_parser(
        "corpus",
        help="List or summarize the demanding GIC regression corpus",
    )
    corpus.add_argument("--root", type=Path, help="Override the GIC corpus root directory")
    corpus.add_argument(
        "--suffix",
        action="append",
        help="Filter by suffix, for example .inp or fchk",
    )
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
    corpus_audit.add_argument(
        "--suffix",
        action="append",
        help="Filter by suffix; defaults to geometry inputs",
    )
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
    if args.command == "gaussian" and args.gaussian_command == "summary":
        from oracle_gaussian import summarize_gaussian_log

        summary = summarize_gaussian_log(args.log)
        print(f"path: {summary.path}")
        print(f"normal_termination: {int(summary.normal_termination)}")
        print(f"scf_count: {len(summary.scf_energies_hartree)}")
        if summary.scf_energies_hartree:
            print(f"last_scf_hartree: {summary.scf_energies_hartree[-1]:.12g}")
        print(f"standard_orientations: {summary.standard_orientation_count}")
        print(f"input_orientations: {summary.input_orientation_count}")
        print(f"frequencies: {len(summary.frequencies_cm)}")
        return 0
    if args.command == "gaussian" and args.gaussian_command == "status":
        from oracle_gaussian import gaussian_job_status

        status = gaussian_job_status(args.workdir)
        print(f"status: {status.status}")
        print(f"workdir: {status.workdir}")
        print(f"log: {status.log_path}")
        if status.input_path is not None:
            print(f"input: {status.input_path}")
        if status.pid is not None:
            print(f"pid: {status.pid}")
        print(f"normal_termination: {int(status.normal_termination)}")
        print(f"error_termination: {int(status.error_termination)}")
        print(f"message: {status.message}")
        return 0
    if args.command == "gaussian" and args.gaussian_command == "run":
        from oracle_gaussian import run_gaussian_job

        result = run_gaussian_job(
            args.workdir,
            executable=args.executable,
            input_path=args.input,
            background=args.background,
            timeout=args.timeout,
        )
        print(f"gaussian_input: {result.input_path}")
        print(f"gaussian_log: {result.log_path}")
        if result.pid is not None:
            print(f"pid: {result.pid}")
        if result.exit_code is not None:
            print(f"exit_code: {result.exit_code}")
        if result.success is not None:
            print(f"success: {int(result.success)}")
        print(f"message: {result.message}")
        return 0
    if args.command == "gaussian" and args.gaussian_command == "formchk":
        from oracle_gaussian import FORMCHK_EXECUTABLE, formchk_checkpoint

        output = formchk_checkpoint(
            args.chk,
            args.fchk,
            executable=args.executable or FORMCHK_EXECUTABLE,
            timeout=args.timeout,
        )
        print(f"Wrote formatted checkpoint: {output}")
        return 0
    if args.command == "gaussian" and args.gaussian_command == "fchk-summary":
        from oracle_gaussian import read_gaussian_fchk_qff

        data = read_gaussian_fchk_qff(args.fchk)
        print(f"path: {args.fchk}")
        print(f"atoms: {len(data.atomic_numbers)}")
        print(f"hessian_lower: {len(data.cartesian_hessian_lower)}")
        print(f"harmonic_frequencies: {len(data.harmonic_frequencies_cm)}")
        print(f"anharmonic_frequencies: {len(data.anharmonic_frequencies_cm)}")
        print(f"anharmonic_e2_values: {len(data.anharmonic_e2)}")
        print(f"normal_mode_values: {len(data.normal_modes)}")
        return 0
    if args.command == "gaussian" and args.gaussian_command == "promote-rovib":
        from oracle_gaussian import promote_gaussian_rovib_to_xyzin

        result = promote_gaussian_rovib_to_xyzin(
            args.log,
            args.xyzin,
            write_vibrational=not args.no_vibrational,
            write_rotational=not args.no_rotational,
            write_deltabvib=not args.no_deltabvib,
            invert_imaginary_modes=not args.no_invert_imaginary,
            exclude_modes=tuple(args.exclude_mode),
        )
        print(f"Promoted Gaussian rovib data: {result.log_path} -> {result.xyzin}")
        print(f"wrote_vibrational: {int(result.wrote_vibrational)}")
        print(f"wrote_rotational: {int(result.wrote_rotational)}")
        print(f"wrote_deltabvib: {int(result.wrote_deltabvib)}")
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
    if args.command == "fragments" and args.fragments_command == "build":
        from oracle_fragments import write_fragment_build_section

        definition = write_fragment_build_section(args.xyzin)
        print(
            "Built ORACLE fragments: "
            f"{args.xyzin} (fragments={len(definition.fragments)}, "
            f"reference={definition.reference_fragment})"
        )
        return 0
    if args.command == "fragments" and args.fragments_command == "centers":
        from oracle_fragments import write_interaction_center_section

        definition = write_interaction_center_section(args.xyzin)
        print(
            "Built ORACLE interaction centers: "
            f"{args.xyzin} (centers={len(definition.centers)}, "
            f"interactions={len(definition.interactions)})"
        )
        return 0
    if args.command == "rovib" and args.rovib_command == "summarize":
        from oracle_rovib import rovib_summary_lines, summarize_xyzin

        print("\n".join(rovib_summary_lines(summarize_xyzin(args.xyzin))))
        return 0
    if args.command == "rovib" and args.rovib_command == "vibin":
        from oracle_rovib import vibin_from_xyzin_fchk

        result = vibin_from_xyzin_fchk(
            args.xyzin,
            args.fchk,
            workdir=args.workdir,
            project_TR=not args.no_project_tr,
            update_vibrational_section=not args.no_update_vibrational,
        )
        print(
            "Built rovib vibin: "
            f"{result.vibin} (nvib={result.data.nvib}, "
            f"imag_like={result.n_imag_like})"
        )
        return 0
    if args.command == "rovib" and args.rovib_command == "coriolis":
        from oracle_rovib import (
            append_coriolis_to_vibin,
            compute_coriolis_from_xyzin,
            coriolis_report_lines,
        )

        result = compute_coriolis_from_xyzin(
            args.xyzin,
            vibin=args.vibin,
            Geff_thr_cm1=args.threshold_cm1,
            only_upper=not args.all_pairs,
        )
        text = "\n".join(coriolis_report_lines(result))
        if args.out is not None:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(text + "\n", encoding="utf-8")
            print(f"Wrote Coriolis report: {args.out}")
        else:
            print(text)
        if args.append_vibin:
            vibin_path = args.vibin or (args.xyzin.parent / "vibin")
            append_coriolis_to_vibin(vibin_path, result)
            print(f"Appended Coriolis block: {vibin_path}")
        return 0
    if args.command == "rovib" and args.rovib_command == "qcent":
        from oracle_rovib import append_qcent_to_vibin, compute_qcent_from_xyzin, qcent_report_lines

        result = compute_qcent_from_xyzin(args.xyzin, vibin=args.vibin)
        text = "\n".join(qcent_report_lines(result))
        if args.out is not None:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(text + "\n", encoding="utf-8")
            print(f"Wrote QCENT report: {args.out}")
        else:
            print(text)
        if args.append_vibin:
            vibin_path = args.vibin or (args.xyzin.parent / "vibin")
            append_qcent_to_vibin(vibin_path, result)
            print(f"Appended QCENT block: {vibin_path}")
        return 0
    if args.command == "rovib" and args.rovib_command == "dos":
        from oracle_rovib import direct_vibrational_dos_from_xyzin

        out = args.out or (args.xyzin.parent / "dos_vib.dat")
        result = direct_vibrational_dos_from_xyzin(
            args.xyzin,
            vmax=args.vmax,
            emax_cm1=args.emax,
            bin_cm1=args.bin_cm1,
            ncap=args.ncap,
            out=out,
        )
        print(f"Wrote vibrational DOS: {result.path} (bins={len(result.bins_logg)})")
        return 0
    if args.command == "rovib" and args.rovib_command == "dos-rovib":
        from oracle_rovib import rovib_pipeline

        result = rovib_pipeline(
            args.xyzin,
            vib_dos=args.vib_dos,
            out=args.out,
            rot_out=args.rot_out,
            q_out=args.q_out,
            emax_rot=args.emax_rot,
            jmax=args.jmax,
        )
        print(f"Wrote rovibrational DOS: {result.dos_rovib}")
        if result.q_path is not None:
            print(f"Wrote rovib Q(T): {result.q_path} (Q={result.Q_rovib:.8g})")
        return 0
    if args.command == "thermo":
        from oracle_thermo import run_thermo_on_xyzin

        result = run_thermo_on_xyzin(
            args.xyzin,
            report=not args.no_report or args.out is not None,
            report_path=args.out,
            write_section=not args.no_write_section,
            cutoff_cm1=args.cutoff_cm1,
            keep_low_positive=args.keep_low_positive,
        )
        total = result.total
        q_text = "" if total is None or total.Q_dimless is None else f", Q={total.Q_dimless:.8g}"
        print(f"Ran ORACLE Thermo: {args.xyzin}{q_text}")
        if args.out is not None:
            print(f"Wrote thermo report: {args.out}")
        elif not args.no_report:
            print(f"Wrote thermo report: {args.xyzin.parent / 'thermo.report'}")
        if not args.no_write_section:
            print(f"Updated #THERMO: {args.xyzin}")
        return 0
    if args.command == "gf":
        from oracle_gf import (
            run_gf_report_from_fchk,
            run_xyzin_gf_report_from_fchk,
            write_csv_tables,
        )

        if args.xyzin is None:
            report = run_gf_report_from_fchk(args.fchk)
            prefix = "gf"
        else:
            report = run_xyzin_gf_report_from_fchk(
                args.fchk,
                args.xyzin,
                scale_path=args.scale_file,
                scale_records=tuple(args.scale),
                local=args.local,
                force_threshold=args.force_threshold,
                block_by_irrep=args.symmetry_blocks,
                subtract_electrostatic=args.subtract_electrostatic,
                subtract_uff_vdw=args.subtract_uff_vdw,
                nonbonded_14_scale=args.nonbonded_14_scale,
            )
            prefix = "gic_gf"
        if args.out is not None:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(report.text + "\n", encoding="utf-8")
            print(f"Wrote GF/PED report: {args.out}")
        else:
            print(report.text)
        if args.csv_dir is not None:
            written = write_csv_tables(report, args.csv_dir, prefix=prefix)
            print(f"Wrote GF/PED CSV tables: {len(written)} files in {args.csv_dir}")
        return 0
    if args.command == "semiexp":
        from oracle_morpheus import (
            DEFAULT_SEMIEXP_OBSERVABLE,
            DEFAULT_SEMIEXP_ROBUST_LOSS,
            DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS,
            HYDROGEN_PARAMETER_CONSTRAINT,
            ParameterClassConstraint,
            QMParameterPredicate,
            SemiexperimentalFitRequest,
            fit_semiexperimental_geometry,
            is_msr_legacy_file,
            prepare_semiexperimental_xyzin,
            read_observations,
            read_semiexperimental_job,
            semiexperimental_latex_tables,
            write_semiexperimental_html_report,
        )

        legacy_msr_job = bool(args.job and is_msr_legacy_file(args.job))
        job = None if legacy_msr_job or not args.job else read_semiexperimental_job(args.job)
        geometry_path = args.xyz or (job.path if job is not None else None)
        observations_inline = job.observations_inline if job is not None else ()
        observations_path = (
            args.observations
            or (None if observations_inline else (job.observations if job is not None else None))
        )
        if legacy_msr_job:
            geometry_path = args.xyz or args.job
            observations_path = args.observations or args.job
            observations_inline = ()
        if geometry_path is None:
            raise ValueError("semiexp needs --geometry or --job")
        if observations_path is None and not observations_inline:
            raise ValueError(
                "semiexp needs --observations, inline [[isotopologues]], "
                "or a [files].observations entry in --job"
            )
        preprocess = prepare_semiexperimental_xyzin(
            Path(geometry_path),
            observations_source=Path(observations_path) if observations_path is not None else None,
            observations_inline=observations_inline,
            xyzin_path=args.xyzin,
        )
        geometry_path = preprocess.xyzin
        observations = read_observations(preprocess.xyzin)
        print(f"semiexp_xyzin: {preprocess.xyzin}")
        if preprocess.created_or_updated_geometry:
            print("semiexp_xyzin_geometry: updated")
        if preprocess.updated_isotopologues:
            print("semiexp_xyzin_isotopologues: updated")

        fixed = _merge_unique(preprocess.source_fixed_parameters, job.fixed_parameters if job else ())
        fixed = _merge_unique(fixed, _parse_fixed_parameters(args.fixed))
        if args.fix_hydrogens:
            fixed = _merge_unique(fixed, (HYDROGEN_PARAMETER_CONSTRAINT,))
        observable = _job_default(
            args.observable,
            DEFAULT_SEMIEXP_OBSERVABLE,
            job.observable if job else None,
        )
        coordinate_model = _job_default(
            args.coordinate_model,
            "gic",
            job.coordinate_model if job else None,
        )
        rotational_components = _job_default(
            args.rotational_components,
            DEFAULT_SEMIEXP_ROTATIONAL_COMPONENTS,
            job.rotational_components if job else None,
        )
        qm_predicates = _merge_unique(
            job.qm_predicates if job else (),
            _parse_qm_predicates(args.qm_predicate, QMParameterPredicate),
        )
        parameter_classes = _merge_unique(
            job.parameter_classes if job else (),
            _parse_parameter_classes(args.parameter_class, ParameterClassConstraint),
        )
        backend = _job_default(args.backend, "python", job.backend if job else None)
        max_iter = args.max_iter if args.max_iter is not None else (job.max_iter if job else None)
        step = _job_default(args.step, 1.0e-4, job.step if job else None)
        damping = _job_default(args.damping, 1.0e-8, job.damping if job else None)
        max_step = _job_default(args.max_step, 0.25, job.max_step if job else None)
        prune_condition = _job_default(
            args.prune_condition,
            0.0,
            job.prune_condition if job else None,
        )
        robust_loss = _job_default(
            args.robust_loss,
            DEFAULT_SEMIEXP_ROBUST_LOSS,
            job.robust_loss if job else None,
        )
        robust_scale = _job_default(args.robust_scale, 0.0, job.robust_scale if job else None)
        leave_one_out = bool(args.leave_one_out or (job.leave_one_out if job else False))
        checkpoint = args.checkpoint if args.checkpoint is not None else (job.checkpoint if job else None)
        restart = args.restart if args.restart is not None else (job.restart if job else None)
        request = SemiexperimentalFitRequest(
            initial_geometry=geometry_path,
            observations=observations,
            fixed_parameters=fixed,
            observable=observable,
            rotational_components=rotational_components,
            qm_predicates=qm_predicates,
            parameter_classes=parameter_classes,
            coordinate_model=coordinate_model,
            robust_loss=robust_loss,
            robust_scale=robust_scale,
            leave_one_out=leave_one_out,
        )
        result = fit_semiexperimental_geometry(
            request,
            max_iter=max_iter,
            step=step,
            damping=damping,
            max_step=max_step,
            prune_condition=prune_condition,
            checkpoint=checkpoint,
            restart=restart,
            outdir=args.outdir,
        )
        report_path = write_semiexperimental_html_report(
            args.outdir / "semiexp_report.html",
            result,
            request,
        )
        tables_path = args.outdir / "semiexp_tables.tex"
        tables = semiexperimental_latex_tables(result)
        tables_path.write_text(
            "\n\n".join(f"% {name}\n{table}" for name, table in tables.items()),
            encoding="utf-8",
        )
        _append_manifest_output(args.outdir / "semiexp_manifest.json", "html_report", report_path)
        _append_manifest_output(args.outdir / "semiexp_manifest.json", "latex_tables", tables_path)
        print(f"manifest: {result.manifest}")
        print(f"report: {report_path}")
        rms_label = "rms_MHz" if result.diagnostics.observable == "rotational_constants" else "rms_observable"
        print(f"{rms_label}: {result.rms_MHz:.8g}")
        rot_diffs = [row.difference_MHz for row in result.rotational_constants]
        rotational_rms = math.sqrt(sum(diff * diff for diff in rot_diffs) / len(rot_diffs)) if rot_diffs else 0.0
        rotational_mse = sum(diff * diff for diff in rot_diffs) / len(rot_diffs) if rot_diffs else 0.0
        print(f"rotational_rms_MHz: {rotational_rms:.8g}")
        print(f"rotational_mean_square_MHz2: {rotational_mse:.8g}")
        print(f"rotational_mean_square_1e3_MHz2: {1000.0 * rotational_mse:.8g}")
        print(f"iterations: {result.iterations}")
        print(f"stationary_point: {result.stationary_point}")
        print(f"convergence: {result.diagnostics.convergence_reason}")
        print(f"rank: {result.diagnostics.rank}")
        print(f"condition_number: {result.diagnostics.condition_number:.8g}")
        print(f"observable: {result.diagnostics.observable}")
        print(f"components: {','.join(result.diagnostics.components)}")
        print(f"backend: {backend}")
        print(f"coordinate_model: {result.diagnostics.coordinate_model}")
        return 0
    if args.command == "semiexp-ensemble":
        from oracle_core.manifest import build_run_manifest
        from oracle_morpheus import fit_ensemble_job

        result = fit_ensemble_job(args.job, outdir=args.outdir)
        outputs = _ensemble_output_paths(args.outdir)
        build_run_manifest(
            workflow="semiexp_ensemble",
            status=result.acceptance.status,
            run_dir=args.outdir,
            inputs={"job": args.job},
            outputs=outputs,
            parameters={
                "classes": len(result.classes),
                "molecules": len(result.molecule_blocks),
                "rank": result.rank,
                "scaled_condition_number": result.condition_number,
                "weighted_rms_before": result.weighted_rms_before,
                "weighted_rms_after": result.weighted_rms_after,
                "accepted": result.acceptance.accepted,
            },
            backend={"solver": "python", "model": "linearized shared class corrections"},
            messages=list(result.acceptance.reasons) + list(result.acceptance.review_items),
        ).write(args.outdir / "run_manifest.json")
        print(f"report: {args.outdir / 'ensemble_class_corrections.txt'}")
        print(f"manifest: {args.outdir / 'run_manifest.json'}")
        print(f"classes: {len(result.classes)}")
        print(f"molecules: {len(result.molecule_blocks)}")
        print(f"rank: {result.rank}")
        print(f"scaled_condition_number: {result.condition_number:.8g}")
        print(f"acceptance_status: {result.acceptance.status}")
        if result.acceptance.reasons:
            print("acceptance_failures: " + " | ".join(result.acceptance.reasons))
        if result.acceptance.review_items:
            print("acceptance_review: " + " | ".join(result.acceptance.review_items))
        print(f"weighted_rms_before: {result.weighted_rms_before:.8g}")
        print(f"weighted_rms_after: {result.weighted_rms_after:.8g}")
        for item in result.classes:
            print(
                f"class:{item.name}: correction={result.corrections[item.name]:.10g} "
                f"sigma={result.sigma[item.name]:.4g}"
            )
        return 0
    if args.command == "gicforge" and args.gicforge_command == "plan":
        from oracle_gicforge import write_gicforge_plan_sections

        plan_kwargs = {
            "symmetrize": args.symmetrize,
            "sycart": args.sycart,
        }
        if args.improper_dihedrals:
            plan_kwargs["improper_dihedrals"] = True
        write_gicforge_plan_sections(
            args.xyzin,
            **plan_kwargs,
        )
        print(f"Planned GICForge workflow: {args.xyzin}")
        return 0
    if args.command == "gicforge" and args.gicforge_command == "build":
        from oracle_gicforge import write_gicforge_build_sections

        build_kwargs = {
            "symmetrize": args.symmetrize,
            "sycart": args.sycart,
        }
        if args.improper_dihedrals:
            build_kwargs["improper_dihedrals"] = True
        definition = write_gicforge_build_sections(
            args.xyzin,
            **build_kwargs,
        )
        print(
            "Built GICForge definition: "
            f"{args.xyzin} (GICs={len(definition.gics)}, rank={definition.rank})"
        )
        return 0
    if args.command == "gicforge" and args.gicforge_command == "bmatrix":
        from oracle_gicforge import (
            build_gic_b_matrix_from_xyzin,
            gic_b_matrix_lines,
            write_gic_b_matrix,
        )

        if args.output is None:
            matrix = build_gic_b_matrix_from_xyzin(args.xyzin)
            print("\n".join(gic_b_matrix_lines(matrix)))
            return 0
        matrix = write_gic_b_matrix(args.xyzin, args.output)
        print(
            "Wrote GIC B matrix: "
            f"{args.output} (rows={len(matrix.rows)}, "
            f"columns={len(matrix.cartesian_columns)})"
        )
        return 0
    if args.command == "gicforge" and args.gicforge_command == "report":
        from oracle_gicforge import gic_report_from_xyzin, write_gic_report

        if args.output is None:
            print("\n".join(gic_report_from_xyzin(args.xyzin)))
            return 0
        output = write_gic_report(args.xyzin, args.output)
        print(f"Wrote GICForge report: {output}")
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


def _parse_fixed_parameters(raw: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in _split_top_level(raw, separators=",;") if part.strip())


def _split_top_level(raw: str, *, separators: str) -> list[str]:
    parts: list[str] = []
    current: list[str] = []
    round_depth = 0
    square_depth = 0
    brace_depth = 0
    for char in str(raw):
        if char == "(":
            round_depth += 1
        elif char == ")" and round_depth > 0:
            round_depth -= 1
        elif char == "[":
            square_depth += 1
        elif char == "]" and square_depth > 0:
            square_depth -= 1
        elif char == "{":
            brace_depth += 1
        elif char == "}" and brace_depth > 0:
            brace_depth -= 1
        if char in separators and round_depth == 0 and square_depth == 0 and brace_depth == 0:
            parts.append("".join(current))
            current = []
            continue
        current.append(char)
    parts.append("".join(current))
    return parts


def _parse_qm_predicates(items: list[str], predicate_type: type) -> tuple:
    predicates = []
    for item in items:
        parts = item.split(":")
        if len(parts) not in {3, 4}:
            raise ValueError("--qm-predicate must be label_pattern:value:sigma[:source]")
        source = parts[3] if len(parts) == 4 else "qm"
        predicates.append(predicate_type(parts[0], float(parts[1]), float(parts[2]), source=source))
    return tuple(predicates)


def _parse_parameter_classes(items: list[str], class_type: type) -> tuple:
    constraints = []
    for item in items:
        parts = item.split(":", 2)
        if len(parts) != 3:
            raise ValueError("--parameter-class must be name:shared|fixed:pattern[|pattern...]")
        patterns = tuple(part.strip() for part in parts[2].split("|") if part.strip())
        constraints.append(class_type(parts[0].strip(), patterns, parts[1].strip()))
    return tuple(constraints)


def _merge_unique(left: tuple[_T, ...], right: tuple[_T, ...]) -> tuple[_T, ...]:
    result: list[_T] = []
    for item in (*left, *right):
        if item not in result:
            result.append(item)
    return tuple(result)


def _job_default(value: _T, default: _T, job_value: _T | None) -> _T:
    if job_value is not None and value == default:
        return job_value
    return value


def _append_manifest_output(manifest_path: Path, name: str, path: Path) -> None:
    if not manifest_path.exists():
        return
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    data.setdefault("outputs", {})[name] = str(path)
    if path.is_file():
        from oracle_core.manifest import sha256_file

        data.setdefault("output_sha256", {})[name] = sha256_file(path)
    manifest_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _ensemble_output_paths(outdir: Path) -> dict[str, Path]:
    root = Path(outdir)
    return {
        "text_report": root / "ensemble_class_corrections.txt",
        "class_corrections_csv": root / "ensemble_class_corrections.csv",
        "class_report_csv": root / "ensemble_class_report.csv",
        "molecule_blocks_csv": root / "ensemble_molecule_blocks.csv",
        "scientific_manifest": root / "ensemble_manifest.json",
        "covariance_csv": root / "ensemble_covariance.csv",
        "correlation_csv": root / "ensemble_correlation.csv",
    }
