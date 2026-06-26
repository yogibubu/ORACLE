from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from typing import Iterable

import numpy as np

from oracle_chem import preprocess_to_enriched_xyz, read_enriched_xyz, write_validation_section
from oracle_engines import gicforge_fortran_layout, run_legacy_gicforge

from .corpus import default_gic_corpus_root
from .definition import build_gic_b_matrix, write_gicforge_build_sections


DEFAULT_FORTRAN_AUDIT_MOLECULES = ("naphtalene.inp", "phenantrene.inp", "pyrene.inp")


@dataclass(frozen=True)
class GICForgeFortranAuditResult:
    molecule: str
    source: Path
    status: str
    oracle_rank: int = 0
    fortran_rank: int = 0
    oracle_shape: tuple[int, int] = (0, 0)
    fortran_shape: tuple[int, int] = (0, 0)
    oracle_row_rank: int = 0
    fortran_row_rank: int = 0
    row_space_residual: float | None = None
    oracle_ring_pucker_components: int = 0
    fortran_label_prefixes: tuple[str, ...] = ()
    message: str = ""
    workdir: Path | None = None

    @property
    def passed(self) -> bool:
        return self.status == "PASS"

    def to_record(self, *, root: Path | None = None) -> dict[str, object]:
        data = asdict(self)
        data["source"] = _display_path(self.source, root=root)
        data["workdir"] = None if self.workdir is None else str(self.workdir)
        return data


@dataclass(frozen=True)
class GICForgeFortranAudit:
    root: Path
    workdir: Path | None
    tolerance: float
    results: tuple[GICForgeFortranAuditResult, ...]

    @property
    def passed(self) -> int:
        return sum(1 for result in self.results if result.status == "PASS")

    @property
    def failed(self) -> int:
        return sum(1 for result in self.results if result.status == "FAIL")

    @property
    def skipped(self) -> int:
        return sum(1 for result in self.results if result.status == "SKIP")

    @property
    def errored(self) -> int:
        return sum(1 for result in self.results if result.status == "ERROR")

    @property
    def max_row_space_residual(self) -> float | None:
        values = [
            float(result.row_space_residual)
            for result in self.results
            if result.row_space_residual is not None and np.isfinite(result.row_space_residual)
        ]
        return max(values) if values else None


def audit_gicforge_fortran_corpus(
    *,
    root: Path | None = None,
    molecules: Iterable[str | Path] | None = None,
    workdir: Path | None = None,
    repo_root: Path | None = None,
    limit: int | None = None,
    tolerance: float = 2.0e-8,
) -> GICForgeFortranAudit:
    corpus_root = Path(root) if root is not None else default_gic_corpus_root(repo_root)
    selected = tuple(molecules or DEFAULT_FORTRAN_AUDIT_MOLECULES)
    if limit is not None:
        if limit < 0:
            raise ValueError("fortran audit limit cannot be negative")
        selected = selected[:limit]

    layout = gicforge_fortran_layout(repo_root)
    if shutil.which("gfortran") is None and not layout.legacy_executable.is_file():
        results = tuple(
            GICForgeFortranAuditResult(
                molecule=str(item),
                source=_resolve_corpus_source(corpus_root, item),
                status="SKIP",
                message="gfortran is not available and the legacy GICForge executable is not built",
            )
            for item in selected
        )
        return GICForgeFortranAudit(
            root=corpus_root,
            workdir=Path(workdir) if workdir is not None else None,
            tolerance=float(tolerance),
            results=results,
        )

    if workdir is None:
        with TemporaryDirectory(prefix="oracle_gicforge_fortran_audit_") as tmp:
            results = _audit_many(
                corpus_root,
                selected,
                workdir=Path(tmp),
                repo_root=repo_root,
                tolerance=float(tolerance),
            )
        return GICForgeFortranAudit(corpus_root, None, float(tolerance), results)

    target = Path(workdir)
    target.mkdir(parents=True, exist_ok=True)
    results = _audit_many(
        corpus_root,
        selected,
        workdir=target,
        repo_root=repo_root,
        tolerance=float(tolerance),
    )
    return GICForgeFortranAudit(corpus_root, target, float(tolerance), results)


def gicforge_fortran_audit_records(
    audit: GICForgeFortranAudit,
) -> list[dict[str, object]]:
    return [result.to_record(root=audit.root) for result in audit.results]


def format_gicforge_fortran_audit_summary(audit: GICForgeFortranAudit) -> list[str]:
    max_residual = audit.max_row_space_residual
    lines = [
        f"ROOT {audit.root}",
        f"WORKDIR {audit.workdir or 'TEMPORARY'}",
        f"TOLERANCE {audit.tolerance:.12g}",
        f"CASES {len(audit.results)}",
        f"PASS {audit.passed}",
        f"FAIL {audit.failed}",
        f"ERROR {audit.errored}",
        f"SKIP {audit.skipped}",
        "MAX_ROW_SPACE_RESIDUAL "
        + ("NONE" if max_residual is None else f"{max_residual:.12g}"),
    ]
    lines.extend(format_gicforge_fortran_audit_cases(audit))
    return lines


def format_gicforge_fortran_audit_cases(
    audit: GICForgeFortranAudit,
    *,
    status: str = "all",
) -> list[str]:
    normalized = status.strip().upper()
    if normalized not in {"ALL", "PASS", "FAIL", "ERROR", "SKIP"}:
        raise ValueError(f"unsupported fortran audit status filter: {status}")
    lines: list[str] = []
    for result in audit.results:
        if normalized != "ALL" and result.status != normalized:
            continue
        residual = "NONE" if result.row_space_residual is None else f"{result.row_space_residual:.12g}"
        lines.append(
            "CASE "
            f"{result.status} "
            f"{_display_path(result.source, root=audit.root)} "
            f"oracle_rank={result.oracle_rank} "
            f"fortran_rank={result.fortran_rank} "
            f"oracle_row_rank={result.oracle_row_rank} "
            f"fortran_row_rank={result.fortran_row_rank} "
            f"row_space_residual={residual}"
            + (f" message={result.message}" if result.message else "")
        )
    return lines


def _audit_many(
    corpus_root: Path,
    selected: tuple[str | Path, ...],
    *,
    workdir: Path,
    repo_root: Path | None,
    tolerance: float,
) -> tuple[GICForgeFortranAuditResult, ...]:
    return tuple(
        _audit_one(
            _resolve_corpus_source(corpus_root, item),
            workdir=workdir / _audit_slug(item),
            repo_root=repo_root,
            tolerance=tolerance,
        )
        for item in selected
    )


def _audit_one(
    source: Path,
    *,
    workdir: Path,
    repo_root: Path | None,
    tolerance: float,
) -> GICForgeFortranAuditResult:
    molecule = source.name
    try:
        case_dir = Path(workdir)
        case_dir.mkdir(parents=True, exist_ok=True)
        xyzin = case_dir / f"{source.stem}.xyzin"
        preprocess_to_enriched_xyz(source, xyzin)
        write_validation_section(xyzin)
        definition = write_gicforge_build_sections(xyzin)
        geometry = read_enriched_xyz(xyzin)
        oracle_b = np.asarray(build_gic_b_matrix(definition).rows, dtype=float)
        legacy = run_legacy_gicforge(
            case_dir / "legacy",
            atoms=geometry.atoms,
            coordinates_angstrom=geometry.coordinates_angstrom,
            point_group="C1",
            title=source.stem,
            keywords=("GNIC", "BMAT"),
            repo_root=repo_root,
        )
        fortran_b = np.asarray(legacy.b_matrix_rows, dtype=float)
        oracle_row_rank = _row_space_rank(oracle_b)
        fortran_row_rank = _row_space_rank(fortran_b)
        residual = _row_space_residual(oracle_b, fortran_b)
        oracle_rank = int(definition.rank)
        fortran_rank = int(legacy.final_counts[-1])
        passed = (
            oracle_rank == fortran_rank
            and oracle_row_rank == fortran_row_rank == oracle_rank
            and oracle_b.shape == fortran_b.shape
            and residual <= tolerance
        )
        return GICForgeFortranAuditResult(
            molecule=molecule,
            source=source,
            status="PASS" if passed else "FAIL",
            oracle_rank=oracle_rank,
            fortran_rank=fortran_rank,
            oracle_shape=tuple(int(value) for value in oracle_b.shape),
            fortran_shape=tuple(int(value) for value in fortran_b.shape),
            oracle_row_rank=oracle_row_rank,
            fortran_row_rank=fortran_row_rank,
            row_space_residual=float(residual),
            oracle_ring_pucker_components=sum(
                1 for gic in definition.gics if gic.family == "RING_PUCKER_COMPONENT"
            ),
            fortran_label_prefixes=tuple(sorted({label[:4] for label in legacy.gic_labels})),
            message="" if passed else "rank, shape or row-space residual mismatch",
            workdir=case_dir,
        )
    except Exception as exc:
        return GICForgeFortranAuditResult(
            molecule=molecule,
            source=source,
            status="ERROR",
            message=f"{type(exc).__name__}: {exc}",
            workdir=workdir,
        )


def _row_space_basis(matrix: np.ndarray) -> tuple[np.ndarray, int]:
    if matrix.size == 0:
        return np.zeros((0, matrix.shape[1] if matrix.ndim == 2 else 0)), 0
    _u, singular_values, vh = np.linalg.svd(matrix, full_matrices=False)
    if not singular_values.size:
        return np.zeros((0, matrix.shape[1])), 0
    tolerance = max(matrix.shape) * float(singular_values[0]) * np.finfo(float).eps * 10.0
    rank = int(np.sum(singular_values > tolerance))
    return vh[:rank], rank


def _row_space_rank(matrix: np.ndarray) -> int:
    return _row_space_basis(matrix)[1]


def _row_space_residual(left: np.ndarray, right: np.ndarray) -> float:
    left_basis, left_rank = _row_space_basis(left)
    right_basis, right_rank = _row_space_basis(right)
    if left_rank != right_rank:
        return float("inf")
    if left_rank == 0:
        return 0.0
    projected = left_basis @ right_basis.T @ right_basis
    return float(np.linalg.norm(left_basis - projected) / np.sqrt(float(left_rank)))


def _resolve_corpus_source(root: Path, item: str | Path) -> Path:
    path = Path(item)
    if path.is_absolute() or path.exists():
        return path
    return Path(root) / path


def _audit_slug(item: str | Path) -> str:
    stem = Path(item).stem if Path(item).suffix else str(item)
    return "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in stem)


def _display_path(path: Path, *, root: Path | None) -> str:
    target = Path(path)
    if root is not None:
        try:
            return str(target.relative_to(root))
        except ValueError:
            pass
    return str(target)
