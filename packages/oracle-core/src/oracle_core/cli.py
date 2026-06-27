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


def build_parser(
    *,
    repo_root: Path | None = None,
    prog: str = "oracle",
) -> argparse.ArgumentParser:
    root = find_repo_root() if repo_root is None else Path(repo_root)
    parser = argparse.ArgumentParser(prog=prog, description="ORACLE/MATRIX workflow CLI")
    sub = parser.add_subparsers(dest="command")

    init = sub.add_parser("init", help="Create an ORACLE project workspace")
    init.add_argument("workdir", type=Path)

    validate = sub.add_parser("validate", help="Validate an enriched XYZ after preprocessing")
    validate.add_argument("xyzin", type=Path)
    validate.add_argument("--require-fragments", action="store_true")

    contracts = sub.add_parser("contracts", help="List standalone xyzin tool contracts")
    contracts.add_argument("--tool", help="Show one tool contract by key or planned name")
    contracts.add_argument("--framework", action="store_true", help="Show the planned MATRIX name")
    contracts.add_argument(
        "--check-xyzin",
        type=Path,
        help="Check whether an xyzin contains the required sections for selected contracts",
    )
    contracts.add_argument(
        "--format",
        choices=("text", "json", "markdown"),
        default="text",
        help="Output format",
    )
    contracts.add_argument("--no-gui", action="store_true", help="Omit the GUI orchestrator")

    babel = sub.add_parser("babel", help="Run ORACLE-Babel preprocessing")
    babel_sub = babel.add_subparsers(dest="babel_command")
    preprocess = babel_sub.add_parser("preprocess", help="Import a source into enriched XYZ")
    preprocess.add_argument("source", type=Path)
    preprocess.add_argument("output", type=Path)
    preprocess.add_argument(
        "--source-kind",
        choices=("auto", "xyz", "enriched_xyz", "gaussian", "molpro", "mrcc"),
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
    gaussian_promote_fchk = gaussian_sub.add_parser(
        "promote-fchk",
        help="Promote Gaussian FCHK Hessian/normal-mode/QFF data into an ORACLE xyzin",
    )
    gaussian_promote_fchk.add_argument("fchk", type=Path)
    gaussian_promote_fchk.add_argument("xyzin", type=Path)
    gaussian_promote_fchk.add_argument("--no-cartesian-hessian", action="store_true")
    gaussian_promote_fchk.add_argument("--no-normal-modes", action="store_true")
    gaussian_promote_fchk.add_argument("--no-qff", action="store_true")
    gaussian_promote_fchk.add_argument("--no-electronic", action="store_true")
    gaussian_promote_fchk.add_argument("--no-orbitals", action="store_true")
    gaussian_promote_electronic = gaussian_sub.add_parser(
        "promote-electronic",
        help="Promote Gaussian electronic states/transitions into an ORACLE xyzin",
    )
    gaussian_promote_electronic.add_argument("log", type=Path)
    gaussian_promote_electronic.add_argument("xyzin", type=Path)
    gaussian_promote_electronic.add_argument("--no-electronic", action="store_true")
    gaussian_promote_electronic.add_argument("--no-transitions", action="store_true")
    gaussian_promote_electronic.add_argument(
        "--orbital-file",
        type=Path,
        action="append",
        default=[],
        help="Register an external Molden/Cube/FCHK orbital or density file in #ORBITALS",
    )
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

    molpro = sub.add_parser("molpro", help="Molpro output adapter utilities")
    molpro_sub = molpro.add_subparsers(dest="molpro_command")
    molpro_summary = molpro_sub.add_parser("summary", help="Summarize a Molpro output")
    molpro_summary.add_argument("output", type=Path)
    molpro_promote = molpro_sub.add_parser(
        "promote",
        help="Preprocess a Molpro output into an ORACLE xyzin",
    )
    molpro_promote.add_argument("output", type=Path)
    molpro_promote.add_argument("xyzin", type=Path)
    molpro_promote.add_argument("--symmetry-distance", type=float, default=1.0e-3)
    molpro_promote.add_argument("--symmetry-inertia", type=float, default=1.0e-3)
    molpro_promote.add_argument("--max-rotation-order", type=int, default=6)

    mrcc = sub.add_parser("mrcc", help="MRCC output adapter utilities")
    mrcc_sub = mrcc.add_subparsers(dest="mrcc_command")
    mrcc_summary = mrcc_sub.add_parser("summary", help="Summarize an MRCC output")
    mrcc_summary.add_argument("output", type=Path)
    mrcc_promote = mrcc_sub.add_parser(
        "promote",
        help="Preprocess an MRCC output into an ORACLE xyzin",
    )
    mrcc_promote.add_argument("output", type=Path)
    mrcc_promote.add_argument("xyzin", type=Path)
    mrcc_promote.add_argument("--symmetry-distance", type=float, default=1.0e-3)
    mrcc_promote.add_argument("--symmetry-inertia", type=float, default=1.0e-3)
    mrcc_promote.add_argument("--max-rotation-order", type=int, default=6)

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
    rovib_wmsrot = rovib_sub.add_parser(
        "wmsrot-input",
        help="Export a WMS-Rot browser input file from normalized xyzin sections",
    )
    rovib_wmsrot.add_argument("xyzin", type=Path)
    rovib_wmsrot.add_argument("--out", type=Path)
    rovib_wmsrot.add_argument("--j-min", type=int, default=0)
    rovib_wmsrot.add_argument("--j-max", type=int, default=30)
    rovib_wmsrot.add_argument("--auto-estimate-j-range", action="store_true")
    rovib_wmsrot.add_argument("--reduction", choices=("A", "S"))
    rovib_wmsrot_run = rovib_sub.add_parser(
        "wmsrot-run",
        help="Run the vendored WMS-Rot Hamiltonian engine on normalized xyzin data",
    )
    rovib_wmsrot_run.add_argument("xyzin", type=Path)
    rovib_wmsrot_run.add_argument("--out", type=Path, required=True)
    rovib_wmsrot_run.add_argument("--j-min", type=int, default=0)
    rovib_wmsrot_run.add_argument("--j-max", type=int, default=30)
    rovib_wmsrot_run.add_argument("--intensity-cut", type=float, default=1.0e-20)
    rovib_wmsrot_run.add_argument("--reduction", choices=("A", "S"))
    rovib_wmsrot_run.add_argument("--no-a-type", action="store_true")
    rovib_wmsrot_run.add_argument("--no-b-type", action="store_true")
    rovib_wmsrot_run.add_argument("--no-c-type", action="store_true")
    rovib_vib_spectrum = rovib_sub.add_parser(
        "vib-spectrum",
        help="Build a broadened IR/Raman/VCD/ROA spectrum from #VIBRATIONAL",
    )
    rovib_vib_spectrum.add_argument("xyzin", type=Path)
    rovib_vib_spectrum.add_argument(
        "--observable",
        choices=("IR", "RAMAN", "VCD", "ROA"),
        default="IR",
    )
    rovib_vib_spectrum.add_argument(
        "--source",
        choices=("harmonic", "anharmonic", "hybrid"),
        default="harmonic",
    )
    rovib_vib_spectrum.add_argument(
        "--level2-xyzin",
        type=Path,
        help="Level-2 xyzin used for hybrid harmonic(level1)+anharmonic correction(level2)",
    )
    rovib_vib_spectrum.add_argument("--csv", type=Path, required=True)
    rovib_vib_spectrum.add_argument("--plot", type=Path)
    rovib_vib_spectrum.add_argument("--peaks", type=Path)
    rovib_vib_spectrum.add_argument("--mode-match-csv", type=Path)
    rovib_vib_spectrum.add_argument("--min-mode-overlap", type=float, default=0.70)
    rovib_vib_spectrum.add_argument("--fwhm-cm1", type=float, default=10.0)
    rovib_vib_spectrum.add_argument("--step-cm1", type=float, default=1.0)
    rovib_vib_spectrum.add_argument(
        "--lineshape",
        choices=("gaussian", "lorentzian"),
        default="gaussian",
    )
    rovib_vib_spectrum.add_argument("--no-normalize", action="store_true")
    rovib_vib_compare = rovib_sub.add_parser(
        "vib-compare",
        help="Compare two vibrational spectra with mirror plotting for IR/Raman",
    )
    rovib_vib_compare.add_argument("xyzin", type=Path)
    rovib_vib_compare.add_argument(
        "second_xyzin",
        type=Path,
        nargs="?",
        help="Optional second xyzin file; defaults to the first file",
    )
    rovib_vib_compare.add_argument(
        "--observable",
        choices=("IR", "RAMAN", "VCD", "ROA"),
        default="IR",
    )
    rovib_vib_compare.add_argument(
        "--first-source",
        choices=("harmonic", "anharmonic", "hybrid"),
        default="harmonic",
    )
    rovib_vib_compare.add_argument(
        "--second-source",
        choices=("harmonic", "anharmonic", "hybrid"),
        default="anharmonic",
    )
    rovib_vib_compare.add_argument("--csv", type=Path, required=True)
    rovib_vib_compare.add_argument("--plot", type=Path)
    rovib_vib_compare.add_argument("--mode-match-csv", type=Path)
    rovib_vib_compare.add_argument("--min-mode-overlap", type=float, default=0.70)
    rovib_vib_compare.add_argument("--fwhm-cm1", type=float, default=10.0)
    rovib_vib_compare.add_argument("--step-cm1", type=float, default=1.0)
    rovib_vib_compare.add_argument(
        "--lineshape",
        choices=("gaussian", "lorentzian"),
        default="gaussian",
    )
    rovib_vib_compare.add_argument("--no-normalize", action="store_true")
    rovib_vib_compare.add_argument(
        "--no-mirror-second",
        action="store_true",
        help="Do not mirror the second spectrum even for IR/Raman",
    )
    rovib_nist_ir = rovib_sub.add_parser(
        "nist-ir",
        help="Download a NIST gas-phase IR JCAMP spectrum and convert it to CSV",
    )
    rovib_nist_ir.add_argument(
        "identifier",
        help="NIST ID, CAS registry number or molecule name",
    )
    rovib_nist_ir.add_argument("--out", type=Path, required=True)
    rovib_nist_ir.add_argument("--index", type=int, default=1)
    rovib_nist_ir.add_argument("--timeout", type=float, default=20.0)
    rovib_dos = rovib_sub.add_parser("dos", help="Build direct vibrational DOS from #VIBRATIONAL")
    rovib_dos.add_argument("xyzin", type=Path)
    rovib_dos.add_argument("--vmax", type=int, default=6)
    rovib_dos.add_argument("--emax", type=float, default=8000.0)
    rovib_dos.add_argument("--bin-cm1", type=float, default=50.0)
    rovib_dos.add_argument("--ncap", type=float, default=10.0)
    rovib_dos.add_argument("--out", type=Path)
    rovib_dos_rovib = rovib_sub.add_parser(
        "dos-rovib", help="Convolve vibrational and rotational DOS"
    )
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
    gf.add_argument("--fchk", type=Path)
    gf.add_argument("--xyzin", type=Path, help="Frozen ORACLE xyzin with a BUILT #GIC section")
    gf.add_argument("--out", type=Path, help="Write the GF/PED report")
    gf.add_argument("--csv-dir", type=Path, help="Write GF/PED CSV tables")
    gf.add_argument("--scale-file", type=Path)
    gf.add_argument("--scale", action="append", default=[])
    gf.add_argument("--local", action="store_true", help="Apply local force-field filtering")
    gf.add_argument("--symmetry-blocks", action="store_true", help="Solve separated irrep blocks")
    gf.add_argument(
        "--force-threshold", type=float, help="Zero internal force constants below threshold"
    )
    gf.add_argument("--no-write-section", action="store_true", help="Do not update #GF_PED")
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

    vpt2_vci = sub.add_parser("vpt2-vci", help="Run VPT2/VCI from normalized ORACLE QFF data")
    vpt2_source = vpt2_vci.add_mutually_exclusive_group(required=True)
    vpt2_source.add_argument("--xyzin", type=Path, help="ORACLE xyzin containing #QFF")
    vpt2_source.add_argument("--fchk", type=Path, help="Gaussian FCHK adapter input")
    vpt2_source.add_argument("--qff-file", type=Path, help="Indexed QFF text adapter input")
    vpt2_source.add_argument(
        "--collect", type=Path, help="Collect post-run VPT2/VCI outputs into #VPT2_VCI"
    )
    vpt2_vci.add_argument("--max-quanta", type=int, default=2)
    vpt2_vci.add_argument("--roots", type=int, default=10)
    vpt2_vci.add_argument("--vci-method", choices=("dense", "davidson"), default="dense")
    vpt2_vci.add_argument("--run-dir", type=Path, help="Write report, CSV tables and manifest here")
    vpt2_vci.add_argument("--out", type=Path, help="Write the readable VPT2/VCI report")
    vpt2_vci.add_argument("--csv-dir", type=Path, help="Write VPT2/VCI CSV tables")
    vpt2_vci.add_argument(
        "--no-write", action="store_true", help="With --collect, do not update #VPT2_VCI"
    )

    dvr = sub.add_parser("dvr", help="Prepare scan/path DVR workflows")
    dvr_sub = dvr.add_subparsers(dest="dvr_command")
    dvr_prepare = dvr_sub.add_parser(
        "prepare", help="Prepare DVR manifest from a Gaussian scan log"
    )
    dvr_prepare.add_argument("log", type=Path)
    dvr_prepare.add_argument("--outdir", type=Path, required=True)
    dvr_prepare.add_argument("--figdir", type=Path)
    dvr_prepare.add_argument("--prefix", default="puckering_dvr")
    dvr_prepare.add_argument("--boundary", default="periodic")
    dvr_prepare.add_argument(
        "--solver",
        choices=("fourier", "sinc-dvr", "fortran-sinc-dvr", "fortran-gaussian"),
        default="fourier",
    )
    dvr_prepare.add_argument("--no-rotconst", action="store_true")
    dvr_prepare.add_argument("--no-cremer-pople", action="store_true")
    dvr_prepare.add_argument("--check-only", action="store_true")
    dvr_prepare.add_argument("--xyzin", type=Path, help="Update this ORACLE xyzin with #DVR")
    dvr_run = dvr_sub.add_parser("run", help="Run DVR directly and update #DVR outputs")
    dvr_run.add_argument("log", type=Path, nargs="?")
    dvr_run.add_argument("--outdir", type=Path)
    dvr_run.add_argument("--figdir", type=Path)
    dvr_run.add_argument("--prefix", default="puckering_dvr")
    dvr_run.add_argument("--boundary", default="periodic")
    dvr_run.add_argument(
        "--solver",
        choices=("fourier", "sinc-dvr", "fortran-sinc-dvr", "fortran-gaussian"),
        default="fourier",
    )
    dvr_run.add_argument("--no-rotconst", action="store_true")
    dvr_run.add_argument("--no-cremer-pople", action="store_true")
    dvr_run.add_argument("--check-only", action="store_true")
    dvr_run.add_argument(
        "--xyzin",
        type=Path,
        help="Read #DVR when LOG is omitted; otherwise update this ORACLE xyzin",
    )
    dvr_run.add_argument("--timeout", type=float)
    dvr_collect = dvr_sub.add_parser("collect", help="Collect post-run DVR outputs into #DVR")
    dvr_collect.add_argument("xyzin", type=Path)
    dvr_collect.add_argument(
        "--no-write",
        action="store_true",
        help="Read and summarize outputs without updating #DVR",
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
    semiexp.add_argument("--no-write-section", action="store_true", help="Do not update #MORPHEUS")
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

    semiexp_ensemble_paper = sub.add_parser(
        "semiexp-ensemble-paper",
        help="Run ensemble paper comparisons and write JPCL-ready artifacts",
    )
    semiexp_ensemble_paper.add_argument("--job", type=Path, required=True)
    semiexp_ensemble_paper.add_argument("--paper-dir", type=Path, required=True)
    semiexp_ensemble_paper.add_argument("--outdir", type=Path)
    semiexp_ensemble_paper.add_argument("--soft-prior-sigma", type=float, default=1.0e-3)

    semiexp_ensemble_prior_scan = sub.add_parser(
        "semiexp-ensemble-prior-scan",
        help="Scan ensemble soft-prior sigma values",
    )
    semiexp_ensemble_prior_scan.add_argument("--job", type=Path, required=True)
    semiexp_ensemble_prior_scan.add_argument("--outdir", type=Path, required=True)
    semiexp_ensemble_prior_scan.add_argument("--sigma", type=float, action="append", default=[])

    semiexp_ensemble_synthon_scan = sub.add_parser(
        "semiexp-ensemble-synthon-scan",
        help="Scan Zeff synthon thresholds for an ensemble job",
    )
    semiexp_ensemble_synthon_scan.add_argument("--job", type=Path, required=True)
    semiexp_ensemble_synthon_scan.add_argument("--outdir", type=Path, required=True)
    semiexp_ensemble_synthon_scan.add_argument(
        "--threshold", type=float, action="append", default=[]
    )

    semiexp_benchmark = sub.add_parser(
        "semiexp-benchmark",
        help="Generate MORPHEUS benchmark and paper tables from a regression snapshot",
    )
    semiexp_benchmark.add_argument("--snapshot", type=Path)
    semiexp_benchmark.add_argument("--outdir", type=Path)
    semiexp_benchmark.add_argument("--no-refresh", action="store_true")
    semiexp_benchmark.add_argument("--update-snapshot", action="store_true")

    trinity = sub.add_parser(
        "trinity",
        help="Prepare TRINITY external energy/gradient geometry optimization state",
    )
    trinity_sub = trinity.add_subparsers(dest="trinity_command")
    trinity_prepare = trinity_sub.add_parser("prepare", help="Write a prepared #TRINITY section")
    trinity_prepare.add_argument("xyzin", type=Path)
    trinity_prepare.add_argument("--run-dir", type=Path, required=True)
    trinity_prepare.add_argument("--engine-command", required=True)
    trinity_prepare.add_argument("--coordinate-model", choices=("gic", "cartesian"), default="gic")
    trinity_prepare.add_argument("--active-space", default="total_symmetric")
    trinity_prepare.add_argument("--max-steps", type=int, default=50)
    trinity_prepare.add_argument("--trust-radius", type=float, default=0.2)
    trinity_prepare.add_argument("--gradient-tolerance", type=float, default=1.0e-5)
    trinity_prepare.add_argument("--step-tolerance", type=float, default=1.0e-5)
    trinity_prepare.add_argument("--energy-tolerance", type=float, default=1.0e-8)
    trinity_prepare.add_argument("--energy-unit", default="hartree")
    trinity_prepare.add_argument("--gradient-unit", default="hartree/bohr")
    trinity_prepare.add_argument("--external-protocol", default="xyz-energy-gradient-v1")
    trinity_status = trinity_sub.add_parser("status", help="Summarize #TRINITY state")
    trinity_status.add_argument("xyzin", type=Path)

    reference_search = sub.add_parser(
        "multistructure-reference-search",
        help="Search the local semiexperimental geometry reference library",
    )
    reference_search.add_argument("--query-xyz", type=Path, required=True)
    reference_search.add_argument("--library-root", type=Path)
    reference_search.add_argument("--outdir", type=Path, required=True)
    reference_search.add_argument("--top-k", type=int, default=10)
    reference_search.add_argument("--ring-weight", type=float, default=0.25)
    reference_search.add_argument("--no-ring-comparison", action="store_true")

    reference_build = sub.add_parser(
        "multistructure-build-reference-geometry",
        help="Build a reference-assisted geometry from local semiexperimental fragments",
    )
    reference_build.add_argument("--query-xyz", type=Path, required=True)
    reference_build.add_argument("--library-root", type=Path)
    reference_build.add_argument("--outdir", type=Path, required=True)
    reference_build.add_argument("--top-library-matches", type=int, default=25)
    reference_build.add_argument("--max-fragment-matches", type=int, default=8)
    reference_build.add_argument("--min-fragment-support", type=int, default=1)
    reference_build.add_argument("--zeff-threshold", type=float, default=0.08)
    reference_build.add_argument("--apply-kind", action="append", default=[])
    reference_build.add_argument("--ring-weight", type=float, default=0.25)
    reference_build.add_argument("--no-ring-comparison", action="store_true")

    gicforge = sub.add_parser(
        "gicforge",
        aliases=("neo",),
        help="Plan NEO/GICForge post-validation sections",
    )
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
    fortran_audit = gicforge_sub.add_parser(
        "fortran-audit",
        help="Compare ORACLE GIC/B rows against the vendored Merlino Fortran backend",
    )
    fortran_audit.add_argument("--root", type=Path, help="Override the GIC corpus root directory")
    fortran_audit.add_argument("--workdir", type=Path, help="Keep audit work directories here")
    fortran_audit.add_argument(
        "--molecule",
        action="append",
        help="Corpus-relative molecule path to audit; repeatable",
    )
    fortran_audit.add_argument("--limit", type=int)
    fortran_audit.add_argument("--tolerance", type=float, default=2.0e-8)
    fortran_audit.add_argument(
        "--format",
        choices=("summary", "cases", "failures", "json"),
        default="summary",
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


def main(
    argv: list[str] | None = None,
    *,
    repo_root: Path | None = None,
    prog: str = "oracle",
) -> int:
    root = find_repo_root() if repo_root is None else Path(repo_root)
    add_repo_packages_to_path(root)
    parser = build_parser(repo_root=root, prog=prog)
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
    if args.command == "contracts":
        from oracle_core import (
            PLANNED_FRAMEWORK_EXPANSION,
            PLANNED_FRAMEWORK_NAME,
            tool_contract,
            tool_contract_lines,
            tool_contract_markdown_table,
            tool_contract_readinesses,
            tool_contracts,
            tool_contracts_json,
            tool_readiness_json,
            tool_readiness_lines,
            tool_readiness_markdown_table,
        )

        if args.framework:
            if args.format == "json":
                print(
                    json.dumps(
                        {
                            "planned_name": PLANNED_FRAMEWORK_NAME,
                            "expanded_name": PLANNED_FRAMEWORK_EXPANSION,
                        },
                        indent=2,
                        sort_keys=True,
                    )
                )
            elif args.format == "markdown":
                print("| Planned name | Expanded name |")
                print("| --- | --- |")
                print(f"| {PLANNED_FRAMEWORK_NAME} | {PLANNED_FRAMEWORK_EXPANSION} |")
            else:
                print(f"framework: {PLANNED_FRAMEWORK_NAME}")
                print(f"  expanded_name: {PLANNED_FRAMEWORK_EXPANSION}")
            return 0

        rows = (
            (tool_contract(args.tool),)
            if args.tool
            else tool_contracts(include_gui=not args.no_gui)
        )
        if args.check_xyzin is not None:
            readinesses = tool_contract_readinesses(args.check_xyzin, rows)
            if args.format == "json":
                print(tool_readiness_json(readinesses))
            elif args.format == "markdown":
                print(tool_readiness_markdown_table(readinesses))
            else:
                print("\n".join(tool_readiness_lines(readinesses)))
            return 0 if all(readiness.ready for readiness in readinesses) else 2
        if args.format == "json":
            print(tool_contracts_json(rows))
        elif args.format == "markdown":
            print(tool_contract_markdown_table(rows))
        else:
            print("\n".join(tool_contract_lines(rows)))
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
    if args.command == "gaussian" and args.gaussian_command == "promote-fchk":
        from oracle_gaussian import promote_gaussian_fchk_to_xyzin

        result = promote_gaussian_fchk_to_xyzin(
            args.fchk,
            args.xyzin,
            write_cartesian_hessian=not args.no_cartesian_hessian,
            write_normal_modes=not args.no_normal_modes,
            write_qff=not args.no_qff,
            write_electronic=not args.no_electronic,
            write_orbitals=not args.no_orbitals,
        )
        print(f"Promoted Gaussian FCHK data: {result.fchk_path} -> {result.xyzin}")
        print(f"wrote_cartesian_hessian: {int(result.wrote_cartesian_hessian)}")
        print(f"wrote_normal_modes: {int(result.wrote_normal_modes)}")
        print(f"wrote_qff: {int(result.wrote_qff)}")
        print(f"wrote_electronic: {int(result.wrote_electronic)}")
        print(f"wrote_orbitals: {int(result.wrote_orbitals)}")
        return 0
    if args.command == "gaussian" and args.gaussian_command == "promote-electronic":
        from oracle_gaussian import promote_gaussian_electronic_log_to_xyzin

        result = promote_gaussian_electronic_log_to_xyzin(
            args.log,
            args.xyzin,
            write_electronic=not args.no_electronic,
            write_transitions=not args.no_transitions,
            orbital_files=tuple(args.orbital_file),
        )
        print(f"Promoted Gaussian electronic data: {result.log_path} -> {result.xyzin}")
        print(f"wrote_electronic: {int(result.wrote_electronic)}")
        print(f"wrote_transitions: {int(result.wrote_transitions)}")
        print(f"wrote_orbitals: {int(result.wrote_orbitals)}")
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
    if args.command == "molpro" and args.molpro_command == "summary":
        from oracle_molpro import summarize_molpro_output

        summary = summarize_molpro_output(args.output)
        print(f"path: {summary.path}")
        print(f"atoms: {summary.geometry.natoms}")
        print(f"charge: {summary.charge}")
        print(f"multiplicity: {summary.multiplicity}")
        print(f"atomic_coordinate_blocks: {summary.atomic_coordinate_blocks}")
        return 0
    if args.command == "molpro" and args.molpro_command == "promote":
        from oracle_molpro import promote_molpro_output_to_xyzin

        result = promote_molpro_output_to_xyzin(
            args.output,
            args.xyzin,
            symmetry_distance=args.symmetry_distance,
            symmetry_inertia=args.symmetry_inertia,
            max_rotation_order=args.max_rotation_order,
        )
        print(
            "Promoted Molpro output: "
            f"{args.output} -> {result.path} ({result.geometry.natoms} atoms, "
            f"PG={result.point_group}, bonds={result.topology_bond_count}, "
            f"rings={result.ring_count})"
        )
        return 0
    if args.command == "mrcc" and args.mrcc_command == "summary":
        from oracle_mrcc import summarize_mrcc_output

        summary = summarize_mrcc_output(args.output)
        print(f"path: {summary.path}")
        print(f"atoms: {summary.geometry.natoms}")
        print(f"charge: {summary.charge}")
        print(f"multiplicity: {summary.multiplicity}")
        print(f"cartesian_coordinate_blocks: {summary.cartesian_coordinate_blocks}")
        return 0
    if args.command == "mrcc" and args.mrcc_command == "promote":
        from oracle_mrcc import promote_mrcc_output_to_xyzin

        result = promote_mrcc_output_to_xyzin(
            args.output,
            args.xyzin,
            symmetry_distance=args.symmetry_distance,
            symmetry_inertia=args.symmetry_inertia,
            max_rotation_order=args.max_rotation_order,
        )
        print(
            "Promoted MRCC output: "
            f"{args.output} -> {result.path} ({result.geometry.natoms} atoms, "
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
    if args.command == "rovib" and args.rovib_command == "wmsrot-input":
        from oracle_rovib import (
            WMSRotInputOptions,
            wmsrot_input_text_from_xyzin,
            write_wmsrot_input,
        )

        options = WMSRotInputOptions(
            j_min=args.j_min,
            j_max=args.j_max,
            auto_estimate_j_range=args.auto_estimate_j_range,
            reduction=args.reduction,
        )
        if args.out is not None:
            out = write_wmsrot_input(args.xyzin, args.out, options=options)
            print(f"Wrote WMS-Rot input: {out}")
        else:
            print(wmsrot_input_text_from_xyzin(args.xyzin, options=options), end="")
        return 0
    if args.command == "rovib" and args.rovib_command == "wmsrot-run":
        from oracle_rovib import (
            WMSRotEngineUnavailable,
            WMSRotSimulationOptions,
            write_wmsrot_spectrum_csv,
        )

        options = WMSRotSimulationOptions(
            j_min=args.j_min,
            j_max=args.j_max,
            intensity_cut=args.intensity_cut,
            reduction=args.reduction,
            a_type=not args.no_a_type,
            b_type=not args.no_b_type,
            c_type=not args.no_c_type,
        )
        try:
            out = write_wmsrot_spectrum_csv(args.xyzin, args.out, options=options)
        except WMSRotEngineUnavailable as exc:
            print(str(exc), file=sys.stderr)
            return 2
        print(f"Wrote WMS-Rot spectrum line list: {out}")
        return 0
    if args.command == "rovib" and args.rovib_command == "vib-spectrum":
        from oracle_rovib import VibrationalSpectrumOptions, write_vibrational_spectrum_outputs

        if args.source == "hybrid" and args.level2_xyzin is None:
            print("rovib vib-spectrum --source hybrid requires --level2-xyzin", file=sys.stderr)
            return 2
        options = VibrationalSpectrumOptions(
            fwhm_cm1=args.fwhm_cm1,
            step_cm1=args.step_cm1,
            lineshape=args.lineshape,
            normalize=not args.no_normalize,
        )
        spectrum = write_vibrational_spectrum_outputs(
            args.xyzin,
            csv_path=args.csv,
            plot_path=args.plot,
            peaks_path=args.peaks,
            level2_xyzin=args.level2_xyzin,
            mode_match_csv_path=args.mode_match_csv,
            observable=args.observable,
            source=args.source,
            min_mode_overlap=args.min_mode_overlap,
            options=options,
        )
        outputs = [f"CSV: {args.csv}"]
        if args.plot is not None:
            outputs.append(f"plot: {args.plot}")
        if args.peaks is not None:
            outputs.append(f"peaks: {args.peaks}")
        if args.mode_match_csv is not None:
            outputs.append(f"mode matches: {args.mode_match_csv}")
        print(
            f"Wrote {spectrum.source} {spectrum.observable} spectrum "
            f"({len(spectrum.peaks)} peaks, {len(spectrum.x_cm1)} points): " + ", ".join(outputs)
        )
        return 0
    if args.command == "rovib" and args.rovib_command == "vib-compare":
        from oracle_rovib import (
            VibrationalSpectrumOptions,
            write_vibrational_spectrum_comparison_outputs,
        )

        if "hybrid" in {args.first_source, args.second_source} and args.second_xyzin is None:
            print(
                "rovib vib-compare with source=hybrid requires a second xyzin file", file=sys.stderr
            )
            return 2
        options = VibrationalSpectrumOptions(
            fwhm_cm1=args.fwhm_cm1,
            step_cm1=args.step_cm1,
            lineshape=args.lineshape,
            normalize=not args.no_normalize,
        )
        comparison = write_vibrational_spectrum_comparison_outputs(
            args.xyzin,
            csv_path=args.csv,
            plot_path=args.plot,
            second_xyzin=args.second_xyzin,
            observable=args.observable,
            first_source=args.first_source,
            second_source=args.second_source,
            options=options,
            mirror_second=False if args.no_mirror_second else None,
            min_mode_overlap=args.min_mode_overlap,
            mode_match_csv_path=args.mode_match_csv,
        )
        outputs = [f"CSV: {args.csv}"]
        if args.plot is not None:
            outputs.append(f"plot: {args.plot}")
        if args.mode_match_csv is not None:
            outputs.append(f"mode matches: {args.mode_match_csv}")
        mirror = "mirrored" if comparison.mirror_second else "not mirrored"
        second_file = args.second_xyzin if args.second_xyzin is not None else args.xyzin
        print(
            f"Wrote {args.observable} comparison "
            f"({args.xyzin} {args.first_source} vs {second_file} {args.second_source}, "
            f"second {mirror}): " + ", ".join(outputs)
        )
        return 0
    if args.command == "rovib" and args.rovib_command == "nist-ir":
        from oracle_rovib import fetch_nist_ir_gas_phase_csv

        result = fetch_nist_ir_gas_phase_csv(
            args.identifier,
            args.out,
            index=args.index,
            timeout=args.timeout,
        )
        if result.status != "downloaded":
            print(result.message, file=sys.stderr)
            return 3
        print(f"{result.message}: {result.csv_path}")
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
            run_xyzin_gf_report_from_xyzin,
            write_csv_tables,
            write_gf_ped_section_from_report,
        )

        if args.xyzin is None:
            if args.fchk is None:
                raise ValueError("gf needs --fchk, or --xyzin containing #CARTESIAN_HESSIAN")
            report = run_gf_report_from_fchk(args.fchk)
            prefix = "gf"
            source_kind = "fchk"
            source_path = args.fchk
        elif args.fchk is None:
            report = run_xyzin_gf_report_from_xyzin(
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
            source_kind = "xyzin"
            source_path = args.xyzin
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
            source_kind = "fchk"
            source_path = args.fchk
        if args.out is not None:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(report.text + "\n", encoding="utf-8")
            print(f"Wrote GF/PED report: {args.out}")
        else:
            print(report.text)
        if args.csv_dir is not None:
            written = write_csv_tables(report, args.csv_dir, prefix=prefix)
            print(f"Wrote GF/PED CSV tables: {len(written)} files in {args.csv_dir}")
        if args.xyzin is not None and not args.no_write_section:
            write_gf_ped_section_from_report(
                args.xyzin,
                report,
                source_kind=source_kind,
                source_path=source_path,
                report_path=args.out,
                csv_dir=args.csv_dir,
            )
            print(f"Updated #GF_PED: {args.xyzin}")
        return 0
    if args.command == "vpt2-vci":
        from oracle_vpt2_vci import (
            VCIOptions,
            collect_vpt2_vci_outputs_from_xyzin,
            refresh_vpt2_vci_section,
            load_force_field,
            run_vpt2_vci_report,
            vpt2_vci_output_summary_lines,
            vpt2_vci_section_from_run,
            write_csv_tables,
            write_vpt2_vci_manifest,
            write_vpt2_vci_section,
        )

        if args.collect is not None:
            snapshot = (
                collect_vpt2_vci_outputs_from_xyzin(args.collect)
                if args.no_write
                else refresh_vpt2_vci_section(args.collect)
            )
            print("\n".join(vpt2_vci_output_summary_lines(snapshot)))
            if not args.no_write:
                print(f"Updated #VPT2_VCI: {args.collect}")
            return 0

        qff = load_force_field(fchk_path=args.fchk, qff_path=args.qff_file, xyzin_path=args.xyzin)
        report = run_vpt2_vci_report(
            qff,
            max_quanta=args.max_quanta,
            roots=args.roots,
            options=VCIOptions(),
            vci_method=args.vci_method,
        )
        report_path = args.out
        csv_dir = args.csv_dir
        run_dir = args.run_dir
        if run_dir is not None:
            run_dir.mkdir(parents=True, exist_ok=True)
            report_path = report_path or (run_dir / "vpt2_vci.report")
            csv_dir = csv_dir or run_dir
        if report_path is not None:
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_path.write_text(report.text + "\n", encoding="utf-8")
            print(f"Wrote VPT2/VCI report: {report_path}")
        else:
            print(report.text)
        written_csv: dict[str, Path] = {}
        if csv_dir is not None:
            written_csv = write_csv_tables(report, csv_dir)
            print(f"Wrote VPT2/VCI CSV tables: {len(written_csv)} files in {csv_dir}")
        outputs = {}
        if report_path is not None:
            outputs["report"] = report_path
        outputs.update(
            {
                key.removesuffix(".csv").replace("mode_contributions", "mode_contributions"): path
                for key, path in written_csv.items()
            }
        )
        manifest_path = None
        source_kind, source_path = _vpt2_vci_source(args)
        if run_dir is not None:
            manifest_path = write_vpt2_vci_manifest(
                run_dir=run_dir,
                inputs={source_kind: source_path} if source_path is not None else {},
                outputs=outputs,
                max_quanta=args.max_quanta,
                roots=args.roots,
                vci_method=args.vci_method,
                source_kind=source_kind,
                status="complete" if "comparison" in outputs and "report" in outputs else "partial",
            )
            outputs["manifest"] = manifest_path
            print(f"manifest: {manifest_path}")
        if args.xyzin is not None and args.xyzin.exists() and outputs:
            status = "complete" if "comparison" in outputs and "report" in outputs else "partial"
            write_vpt2_vci_section(
                args.xyzin,
                vpt2_vci_section_from_run(
                    source_kind=source_kind,
                    source_path=source_path,
                    run_dir=run_dir,
                    report_path=report_path,
                    csv_dir=csv_dir,
                    manifest_path=manifest_path,
                    max_quanta=args.max_quanta,
                    roots=args.roots,
                    vci_method=args.vci_method,
                    outputs=outputs,
                    status=status,
                ),
            )
            print(f"Updated #VPT2_VCI: {args.xyzin}")
        return 0
    if args.command == "dvr" and args.dvr_command == "prepare":
        import shlex

        from oracle_dvr import (
            DVRRequest,
            build_fortran_bridge_args,
            build_fortran_shell_command,
            build_path_analysis_args,
            dvr_section_from_request,
            is_fortran_solver,
            resolve_dvr_executable,
            write_dvr_manifest,
            write_dvr_section,
        )

        request = DVRRequest(
            repo_root=root,
            log_path=args.log,
            outdir=args.outdir,
            figdir=args.figdir or (args.outdir / "figures"),
            prefix=args.prefix,
            boundary=args.boundary,
            solver=args.solver,
            compute_rotconst=not args.no_rotconst,
            label_cremer_pople=not args.no_cremer_pople,
            check_only=args.check_only,
        )
        python_args = build_path_analysis_args(request)
        command = f"{shlex.quote(request.python_executable)} {shlex.join(python_args)}"
        manifest_args = python_args
        if is_fortran_solver(request.solver):
            bridge_args = build_fortran_bridge_args(request, resolve_dvr_executable(root))
            command = build_fortran_shell_command(request, python_args, bridge_args)
            manifest_args = [command]
        manifest = write_dvr_manifest(request, manifest_args)
        if args.xyzin is not None:
            write_dvr_section(
                args.xyzin,
                dvr_section_from_request(request, manifest_path=manifest),
            )
            print(f"Updated #DVR: {args.xyzin}")
        print(f"manifest: {manifest}")
        print(f"command: {command}")
        return 0
    if args.command == "dvr" and args.dvr_command == "run":
        from oracle_dvr import (
            DVRRequest,
            dvr_output_summary_lines,
            dvr_request_from_section,
            dvr_section_from_request,
            read_dvr_section,
            refresh_dvr_section,
            run_dvr_request,
            write_dvr_section,
        )

        if args.log is None:
            if args.xyzin is None:
                raise ValueError("dvr run needs LOG --outdir, or --xyzin containing #DVR")
            request = dvr_request_from_section(read_dvr_section(args.xyzin), repo_root=root)
        else:
            if args.outdir is None:
                raise ValueError("dvr run with LOG needs --outdir")
            request = DVRRequest(
                repo_root=root,
                log_path=args.log,
                outdir=args.outdir,
                figdir=args.figdir or (args.outdir / "figures"),
                prefix=args.prefix,
                boundary=args.boundary,
                solver=args.solver,
                compute_rotconst=not args.no_rotconst,
                label_cremer_pople=not args.no_cremer_pople,
                check_only=args.check_only,
            )

        result = run_dvr_request(request, timeout=args.timeout)
        print(f"manifest: {result.manifest_path}")
        print(f"status: {result.status}")
        if args.xyzin is not None:
            write_dvr_section(
                args.xyzin,
                dvr_section_from_request(
                    request,
                    manifest_path=result.manifest_path,
                    status=result.status,
                ),
            )
            snapshot = refresh_dvr_section(args.xyzin)
            print("\n".join(dvr_output_summary_lines(snapshot)))
            print(f"Updated #DVR: {args.xyzin}")
        return 0
    if args.command == "dvr" and args.dvr_command == "collect":
        from oracle_dvr import (
            collect_dvr_outputs_from_xyzin,
            dvr_output_summary_lines,
            refresh_dvr_section,
        )

        snapshot = (
            collect_dvr_outputs_from_xyzin(args.xyzin)
            if args.no_write
            else refresh_dvr_section(args.xyzin)
        )
        print("\n".join(dvr_output_summary_lines(snapshot)))
        if not args.no_write:
            print(f"Updated #DVR: {args.xyzin}")
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
            write_morpheus_section_from_result,
            write_semiexperimental_html_report,
        )

        legacy_msr_job = bool(args.job and is_msr_legacy_file(args.job))
        job = None if legacy_msr_job or not args.job else read_semiexperimental_job(args.job)
        geometry_path = args.xyz or (job.path if job is not None else None)
        observations_inline = job.observations_inline if job is not None else ()
        observations_path = args.observations or (
            None if observations_inline else (job.observations if job is not None else None)
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

        fixed = _merge_unique(
            preprocess.source_fixed_parameters, job.fixed_parameters if job else ()
        )
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
        checkpoint = (
            args.checkpoint if args.checkpoint is not None else (job.checkpoint if job else None)
        )
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
        if args.xyzin is not None and not args.no_write_section:
            write_morpheus_section_from_result(
                preprocess.xyzin,
                result,
                outdir=args.outdir,
                backend=backend,
                source_path=args.job or args.xyz or observations_path,
                html_report_path=report_path,
                latex_tables_path=tables_path,
            )
            print(f"updated_morpheus_section: {preprocess.xyzin}")
        print(f"manifest: {result.manifest}")
        print(f"report: {report_path}")
        rms_label = (
            "rms_MHz"
            if result.diagnostics.observable == "rotational_constants"
            else "rms_observable"
        )
        print(f"{rms_label}: {result.rms_MHz:.8g}")
        rot_diffs = [row.difference_MHz for row in result.rotational_constants]
        rotational_rms = (
            math.sqrt(sum(diff * diff for diff in rot_diffs) / len(rot_diffs)) if rot_diffs else 0.0
        )
        rotational_mse = (
            sum(diff * diff for diff in rot_diffs) / len(rot_diffs) if rot_diffs else 0.0
        )
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
    if args.command == "semiexp-ensemble-paper":
        from oracle_morpheus import write_ensemble_jpcl_artifacts

        artifacts = write_ensemble_jpcl_artifacts(
            args.job,
            args.paper_dir,
            outdir=args.outdir,
            soft_prior_sigma=args.soft_prior_sigma,
        )
        print(f"paper_dir: {args.paper_dir}")
        if args.outdir is not None:
            print(f"analysis_dir: {args.outdir}")
        for name, path in sorted(artifacts.items()):
            print(f"{name}: {path}")
        return 0
    if args.command == "semiexp-ensemble-prior-scan":
        from oracle_morpheus import run_ensemble_prior_scan

        kwargs = {}
        if args.sigma:
            kwargs["sigmas"] = tuple(args.sigma)
        rows = run_ensemble_prior_scan(args.job, args.outdir, **kwargs)
        print(f"rows: {len(rows)}")
        print(f"csv: {args.outdir / 'prior_sigma_scan.csv'}")
        print(f"json: {args.outdir / 'prior_sigma_scan.json'}")
        return 0
    if args.command == "semiexp-ensemble-synthon-scan":
        from oracle_morpheus import run_ensemble_synthon_threshold_scan

        kwargs = {}
        if args.threshold:
            kwargs["thresholds"] = tuple(args.threshold)
        rows = run_ensemble_synthon_threshold_scan(args.job, args.outdir, **kwargs)
        print(f"rows: {len(rows)}")
        print(f"csv: {args.outdir / 'synthon_threshold_scan.csv'}")
        print(f"json: {args.outdir / 'synthon_threshold_scan.json'}")
        return 0
    if args.command == "semiexp-benchmark":
        from oracle_morpheus import generate_paper_benchmark_artifacts

        snapshot, artifacts = generate_paper_benchmark_artifacts(
            snapshot_path=args.snapshot,
            outdir=args.outdir,
            refresh_from_outputs=not args.no_refresh,
            update_snapshot=args.update_snapshot,
        )
        print(f"cases: {len(snapshot.get('cases', {}))}")
        print(f"planar_diagnostics: {len(snapshot.get('planar_pair_diagnostics', {}))}")
        for name, path in sorted(artifacts.items()):
            print(f"{name}: {path}")
        return 0
    if args.command == "trinity" and args.trinity_command == "prepare":
        from oracle_trinity import prepare_trinity_section

        section = prepare_trinity_section(
            args.xyzin,
            run_dir=args.run_dir,
            engine_command=args.engine_command,
            coordinate_model=args.coordinate_model,
            active_space=args.active_space,
            max_steps=args.max_steps,
            trust_radius=args.trust_radius,
            gradient_tolerance=args.gradient_tolerance,
            step_tolerance=args.step_tolerance,
            energy_tolerance=args.energy_tolerance,
            energy_unit=args.energy_unit,
            gradient_unit=args.gradient_unit,
            external_protocol=args.external_protocol,
        )
        print(f"Updated #TRINITY: {args.xyzin}")
        print(f"manifest: {section.manifest_path}")
        print(f"run_dir: {section.run_dir}")
        print(f"coordinate_model: {section.coordinate_model}")
        print(f"engine_command: {section.engine_command}")
        return 0
    if args.command == "trinity" and args.trinity_command == "status":
        from oracle_trinity import read_trinity_section, trinity_section_summary_lines

        print("\n".join(trinity_section_summary_lines(read_trinity_section(args.xyzin))))
        return 0
    if args.command == "multistructure-reference-search":
        from oracle_morpheus import search_reference_library

        result = search_reference_library(
            args.query_xyz,
            library_root=args.library_root,
            top_k=args.top_k,
            include_ring_comparison=not args.no_ring_comparison,
            ring_weight=args.ring_weight,
            outdir=args.outdir,
        )
        print(f"reference_library: {result.library_root}")
        print(f"matches: {len(result.matches)}")
        print(f"skipped: {len(result.skipped)}")
        if result.matches:
            top = result.matches[0]
            print(f"top_match: {top.slug} similarity={top.similarity_combined:.8g}")
        print(f"outputs: {args.outdir}")
        return 0
    if args.command == "multistructure-build-reference-geometry":
        from oracle_morpheus import build_reference_assisted_geometry

        apply_kinds = tuple(args.apply_kind) if args.apply_kind else None
        kwargs = {}
        if apply_kinds is not None:
            kwargs["apply_kinds"] = apply_kinds
        result = build_reference_assisted_geometry(
            args.query_xyz,
            library_root=args.library_root,
            top_library_matches=args.top_library_matches,
            max_fragment_matches=args.max_fragment_matches,
            min_fragment_support=args.min_fragment_support,
            zeff_threshold=args.zeff_threshold,
            include_ring_comparison=not args.no_ring_comparison,
            ring_weight=args.ring_weight,
            outdir=args.outdir,
            **kwargs,
        )
        print(f"reference_library: {result.library_root}")
        print(f"targets: {len(result.targets)}")
        print(f"unmatched: {len(result.unmatched)}")
        print(f"iterations: {result.iterations}")
        print(f"rms_target_residual_final: {result.rms_target_residual_final:.8g}")
        print(f"outputs: {args.outdir}")
        return 0
    if _is_neo_command(args) and args.gicforge_command == "plan":
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
    if _is_neo_command(args) and args.gicforge_command == "build":
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
    if _is_neo_command(args) and args.gicforge_command == "bmatrix":
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
    if _is_neo_command(args) and args.gicforge_command == "report":
        from oracle_gicforge import gic_report_from_xyzin, write_gic_report

        if args.output is None:
            print("\n".join(gic_report_from_xyzin(args.xyzin)))
            return 0
        output = write_gic_report(args.xyzin, args.output)
        print(f"Wrote GICForge report: {output}")
        return 0
    if _is_neo_command(args) and args.gicforge_command == "corpus":
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
    if _is_neo_command(args) and args.gicforge_command == "corpus-audit":
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
    if _is_neo_command(args) and args.gicforge_command == "fortran-audit":
        from oracle_gicforge import (
            audit_gicforge_fortran_corpus,
            default_gic_corpus_root,
            format_gicforge_fortran_audit_cases,
            format_gicforge_fortran_audit_summary,
            gicforge_fortran_audit_records,
        )

        corpus_root = args.root or default_gic_corpus_root(root)
        audit = audit_gicforge_fortran_corpus(
            root=corpus_root,
            molecules=args.molecule,
            workdir=args.workdir,
            repo_root=root,
            limit=args.limit,
            tolerance=args.tolerance,
        )
        if args.format == "json":
            payload = {
                "root": str(audit.root),
                "workdir": None if audit.workdir is None else str(audit.workdir),
                "tolerance": audit.tolerance,
                "cases": len(audit.results),
                "passed": audit.passed,
                "failed": audit.failed,
                "errored": audit.errored,
                "skipped": audit.skipped,
                "max_row_space_residual": audit.max_row_space_residual,
                "results": gicforge_fortran_audit_records(audit),
            }
            print(json.dumps(payload, indent=2, sort_keys=True))
            return 0
        if args.format == "failures":
            lines = format_gicforge_fortran_audit_cases(audit, status="fail")
            lines.extend(format_gicforge_fortran_audit_cases(audit, status="error"))
            print("\n".join(lines))
            return 0
        if args.format == "cases":
            print("\n".join(format_gicforge_fortran_audit_cases(audit)))
            return 0
        print("\n".join(format_gicforge_fortran_audit_summary(audit)))
        return 0
    if _is_neo_command(args) and args.gicforge_command == "gaussian-input":
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


def matrix_main(argv: list[str] | None = None) -> int:
    """Console-script alias for the MATRIX framework CLI."""
    return main(argv, prog="matrix")


def neo_main(argv: list[str] | None = None) -> int:
    """Console-script alias for the NEO/GICForge coordinate tool."""
    command_args = sys.argv[1:] if argv is None else argv
    return main(["neo", *command_args], prog="neo")


def _is_neo_command(args: argparse.Namespace) -> bool:
    return args.command in {"gicforge", "neo"}


def _parse_fixed_parameters(raw: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in _split_top_level(raw, separators=",;") if part.strip())


def _vpt2_vci_source(args) -> tuple[str, Path | None]:
    if getattr(args, "xyzin", None) is not None:
        return "xyzin", args.xyzin
    if getattr(args, "fchk", None) is not None:
        return "fchk", args.fchk
    if getattr(args, "qff_file", None) is not None:
        return "qff_file", args.qff_file
    return "unknown", None


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
