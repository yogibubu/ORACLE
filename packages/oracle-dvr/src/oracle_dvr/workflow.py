from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

from oracle_core import (
    build_run_manifest,
    key_value_section_lines,
    parse_key_value_section,
    read_sectioned_lines,
    replace_section,
    section_content,
)
from oracle_engines import dvr_executable, puckering_dvr_layout


FORTRAN_SOLVERS = {"fortran-sinc-dvr", "fortran-gaussian"}
ORACLE_XYZ_DVR_SCHEMA = "oracle.xyz.dvr.v1"


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


@dataclass(frozen=True)
class DVRSection:
    log_path: Path
    outdir: Path
    figdir: Path
    prefix: str = "puckering_dvr"
    boundary: str = "periodic"
    solver: str = "fourier"
    compute_rotconst: bool = True
    label_cremer_pople: bool = True
    check_only: bool = False
    manifest_path: Path | None = None
    grid_csv: Path | None = None
    summary: Path | None = None
    levels: Path | None = None
    status: str = "prepared"
    source: str = "gaussian-log"
    schema: str = ORACLE_XYZ_DVR_SCHEMA

    def __post_init__(self) -> None:
        object.__setattr__(self, "log_path", Path(self.log_path))
        object.__setattr__(self, "outdir", Path(self.outdir))
        object.__setattr__(self, "figdir", Path(self.figdir))
        for attr in ("manifest_path", "grid_csv", "summary", "levels"):
            value = getattr(self, attr)
            if value is not None:
                object.__setattr__(self, attr, Path(value))
        object.__setattr__(self, "prefix", self.prefix or "puckering_dvr")
        object.__setattr__(self, "compute_rotconst", bool(self.compute_rotconst))
        object.__setattr__(self, "label_cremer_pople", bool(self.label_cremer_pople))
        object.__setattr__(self, "check_only", bool(self.check_only))


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


def dvr_section_from_request(
    request: DVRRequest,
    *,
    manifest_path: Path | None = None,
    status: str = "prepared",
) -> DVRSection:
    return DVRSection(
        log_path=request.log_path,
        outdir=request.outdir,
        figdir=request.figdir,
        prefix=request.normalized_prefix,
        boundary=request.boundary,
        solver=request.solver,
        compute_rotconst=request.compute_rotconst,
        label_cremer_pople=request.label_cremer_pople,
        check_only=request.check_only,
        manifest_path=manifest_path,
        grid_csv=request.grid_csv,
        summary=request.outdir / f"{request.normalized_prefix}_summary.txt",
        levels=request.outdir / f"{request.normalized_prefix}_levels.csv",
        status=status,
    )


def dvr_section_lines(section: DVRSection) -> list[str]:
    return key_value_section_lines(
        ORACLE_XYZ_DVR_SCHEMA,
        {
            "SOURCE": section.source,
            "STATUS": section.status,
            "LOG_PATH": section.log_path,
            "OUTDIR": section.outdir,
            "FIGDIR": section.figdir,
            "PREFIX": section.prefix,
            "BOUNDARY": section.boundary,
            "SOLVER": section.solver,
            "COMPUTE_ROTCONST": int(section.compute_rotconst),
            "LABEL_CREMER_POPLE": int(section.label_cremer_pople),
            "CHECK_ONLY": int(section.check_only),
            "MANIFEST": section.manifest_path,
            "GRID_CSV": section.grid_csv,
            "SUMMARY": section.summary,
            "LEVELS": section.levels,
        },
        key_order=(
            "SOURCE",
            "STATUS",
            "LOG_PATH",
            "OUTDIR",
            "FIGDIR",
            "PREFIX",
            "BOUNDARY",
            "SOLVER",
            "COMPUTE_ROTCONST",
            "LABEL_CREMER_POPLE",
            "CHECK_ONLY",
            "MANIFEST",
            "GRID_CSV",
            "SUMMARY",
            "LEVELS",
        ),
    )


def parse_dvr_section(lines: list[str] | tuple[str, ...]) -> DVRSection:
    values = parse_key_value_section(lines)
    schema = values.get("SCHEMA", ORACLE_XYZ_DVR_SCHEMA)
    if schema != ORACLE_XYZ_DVR_SCHEMA:
        raise ValueError(f"unsupported DVR schema: {schema}")
    outdir = _path_value(values, "OUTDIR", Path("."))
    prefix = values.get("PREFIX", "puckering_dvr") or "puckering_dvr"
    return DVRSection(
        log_path=_path_value(values, "LOG_PATH", Path("")),
        outdir=outdir,
        figdir=_path_value(values, "FIGDIR", outdir / "figures"),
        prefix=prefix,
        boundary=values.get("BOUNDARY", "periodic"),
        solver=values.get("SOLVER", "fourier"),
        compute_rotconst=_bool_value(values, "COMPUTE_ROTCONST", True),
        label_cremer_pople=_bool_value(values, "LABEL_CREMER_POPLE", True),
        check_only=_bool_value(values, "CHECK_ONLY", False),
        manifest_path=_optional_path(values.get("MANIFEST")),
        grid_csv=_optional_path(values.get("GRID_CSV")) or outdir / f"{prefix}_grid.csv",
        summary=_optional_path(values.get("SUMMARY")) or outdir / f"{prefix}_summary.txt",
        levels=_optional_path(values.get("LEVELS")) or outdir / f"{prefix}_levels.csv",
        status=values.get("STATUS", "prepared"),
        source=values.get("SOURCE", "gaussian-log"),
        schema=schema,
    )


def read_dvr_section(path: Path | str) -> DVRSection:
    content = section_content(read_sectioned_lines(Path(path)), "DVR")
    if not content:
        raise ValueError("missing #DVR section")
    return parse_dvr_section(content)


def write_dvr_section(path: Path | str, section: DVRSection) -> None:
    replace_section(Path(path), "DVR", dvr_section_lines(section))


def _optional_path(raw: str | None) -> Path | None:
    if raw is None or not raw.strip():
        return None
    return Path(raw)


def _path_value(values: dict[str, str], key: str, default: Path) -> Path:
    return _optional_path(values.get(key)) or default


def _bool_value(values: dict[str, str], key: str, default: bool) -> bool:
    raw = values.get(key)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "t", "yes", "y", "on"}
