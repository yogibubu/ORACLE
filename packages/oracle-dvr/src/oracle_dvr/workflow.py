from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

from oracle_core import build_run_manifest
from oracle_engines import dvr_executable, puckering_dvr_layout


FORTRAN_SOLVERS = {"fortran-sinc-dvr", "fortran-gaussian"}


@dataclass(frozen=True)
class DVRRequest:
    repo_root: Path
    log_path: Path
    outdir: Path
    figdir: Path
    prefix: str = "puckering_dvr"
    boundary: str = "periodic"
    solver: str = "fourier"
    compute_rotconst: bool = True
    label_cremer_pople: bool = True
    check_only: bool = False
    python_executable: str = sys.executable

    @property
    def dvr_root(self) -> Path:
        return puckering_dvr_layout(self.repo_root).engine_dir

    @property
    def script(self) -> Path:
        return puckering_dvr_layout(self.repo_root).path_analysis_script

    @property
    def fortran_bridge(self) -> Path:
        return puckering_dvr_layout(self.repo_root).fortran_bridge_script

    @property
    def normalized_prefix(self) -> str:
        return self.prefix or "puckering_dvr"

    @property
    def grid_csv(self) -> Path:
        return self.outdir / f"{self.normalized_prefix}_grid.csv"


def is_fortran_solver(solver: str) -> bool:
    return solver in FORTRAN_SOLVERS


def effective_python_solver(solver: str) -> str:
    return "sinc-dvr" if is_fortran_solver(solver) else solver


def resolve_dvr_executable(repo_root: Path) -> Path:
    return dvr_executable(repo_root)


def build_path_analysis_args(request: DVRRequest) -> list[str]:
    args = [
        str(request.script),
        "--gaussian-log",
        str(request.log_path),
        "--log-selection",
        "last-per-link",
        "--boundary",
        request.boundary,
        "--solver",
        effective_python_solver(request.solver),
        "--outdir",
        str(request.outdir),
        "--figdir",
        str(request.figdir),
        "--prefix",
        request.normalized_prefix,
    ]
    if request.compute_rotconst:
        args.append("--compute-rotconst")
    if request.label_cremer_pople:
        args.append("--label-cremer-pople")
    if request.check_only:
        args.append("--check-only")
    return args


def build_fortran_bridge_args(request: DVRRequest, fortran_exe: Path) -> list[str]:
    args = [
        str(request.fortran_bridge),
        "--grid-csv",
        str(request.grid_csv),
        "--exe",
        str(fortran_exe),
        "--outdir",
        str(request.outdir),
        "--prefix",
        request.normalized_prefix,
        "--boundary",
        request.boundary,
    ]
    if request.solver == "fortran-gaussian":
        args.extend(["--mode", "gaussian"])
    return args


def build_fortran_shell_command(
    request: DVRRequest,
    python_args: list[str],
    bridge_args: list[str],
) -> str:
    return (
        f"{_quote(request.python_executable)} {' '.join(_quote(arg) for arg in python_args)}"
        f" && {_quote(request.python_executable)} {' '.join(_quote(arg) for arg in bridge_args)}"
    )


def _quote(value: str | Path) -> str:
    text = str(value)
    return "'" + text.replace("'", "'\"'\"'") + "'"


def write_dvr_manifest(request: DVRRequest, args: list[str], *, status: str = "prepared") -> Path:
    """Write a standard ORACLE run manifest for a DVR request."""
    outputs = {
        "grid_csv": request.grid_csv,
        "summary": request.outdir / f"{request.normalized_prefix}_summary.txt",
        "levels": request.outdir / f"{request.normalized_prefix}_levels.csv",
    }
    manifest = build_run_manifest(
        workflow="dvr",
        status=status,
        run_dir=request.outdir,
        inputs={"gaussian_log": request.log_path},
        outputs={name: path for name, path in outputs.items() if path.exists()},
        parameters={
            "prefix": request.normalized_prefix,
            "boundary": request.boundary,
            "solver": request.solver,
            "compute_rotconst": request.compute_rotconst,
            "label_cremer_pople": request.label_cremer_pople,
            "check_only": request.check_only,
        },
        backend={"python_executable": request.python_executable, "args": args},
    )
    return manifest.write(request.outdir / f"{request.normalized_prefix}_manifest.json")
