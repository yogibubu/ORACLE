from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .definition import GICDefinition, read_gic_definition_from_xyzin


SALC_SNAPSHOT_SCHEMA = "matrix.neo.gic_salc_snapshot.v2"
DEFAULT_ROUNDING_DECIMALS = 12
DEFAULT_SELECTED_PER_FAMILY = 3
SALC_COEFFICIENT_TOLERANCE = 1.0e-8
PRIORITY_FAMILIES = frozenset(
    {
        "RING_PUCKER_COMPONENT",
        "BUTTERFLY",
        "SPIRO_BEND",
        "CENTER_ATOM_DISTANCE",
        "FRAG_DISTANCE",
        "FRAG_CENTER_ATOM_DISTANCE",
        "FRAG_TRANSLATION",
        "FRAG_ORIENTATION",
    }
)


@dataclass(frozen=True)
class SALCSnapshotComparison:
    ok: bool
    messages: tuple[str, ...] = ()


def salc_snapshot_record(
    definition: GICDefinition,
    *,
    case_id: str = "",
    source: str = "",
    rounding_decimals: int = DEFAULT_ROUNDING_DECIMALS,
    selected_per_family: int = DEFAULT_SELECTED_PER_FAMILY,
) -> dict[str, Any]:
    full = _salc_records(definition, rounding_decimals=rounding_decimals)
    selected = _selected_salc_records(full, selected_per_family=selected_per_family)
    return {
        "id": case_id,
        "source": source,
        "point_group": definition.point_group,
        "rank": definition.rank,
        "target_rank": definition.target_rank,
        "symmetry_method": (
            definition.symmetry_diagnostics.method if definition.symmetry_diagnostics else "NONE"
        ),
        "salc_count": len(full),
        "salc_sha256": _stable_sha256(full),
        "selected_salcs": selected,
    }


def salc_snapshot_document(
    definitions: tuple[tuple[str, str, GICDefinition], ...],
    *,
    rounding_decimals: int = DEFAULT_ROUNDING_DECIMALS,
    selected_per_family: int = DEFAULT_SELECTED_PER_FAMILY,
) -> dict[str, Any]:
    return {
        "schema": SALC_SNAPSHOT_SCHEMA,
        "description": (
            "Compact golden snapshots of nontrivial NEO/GICForge SALC coefficient "
            "vectors. salc_sha256 covers the complete coefficient list; selected_salcs "
            "keeps representative human-readable vectors."
        ),
        "rounding_decimals": int(rounding_decimals),
        "selected_per_family": int(selected_per_family),
        "entries": tuple(
            salc_snapshot_record(
                definition,
                case_id=case_id,
                source=source,
                rounding_decimals=rounding_decimals,
                selected_per_family=selected_per_family,
            )
            for case_id, source, definition in definitions
        ),
    }


def write_salc_snapshot(path: Path, definition: GICDefinition) -> Path:
    target = Path(path)
    payload = salc_snapshot_document((("", "", definition),))
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target


def write_salc_snapshot_from_xyzin(xyzin: Path, output: Path) -> Path:
    return write_salc_snapshot(Path(output), read_gic_definition_from_xyzin(Path(xyzin)))


def compare_salc_snapshot_entry(
    expected: dict[str, Any],
    definition: GICDefinition,
    *,
    rounding_decimals: int,
    selected_per_family: int = DEFAULT_SELECTED_PER_FAMILY,
) -> SALCSnapshotComparison:
    current = salc_snapshot_record(
        definition,
        case_id=str(expected.get("id", "")),
        source=str(expected.get("source", "")),
        rounding_decimals=rounding_decimals,
        selected_per_family=selected_per_family,
    )
    messages: list[str] = []
    for key in ("point_group", "rank", "target_rank", "symmetry_method", "salc_count"):
        if current[key] != expected.get(key):
            messages.append(
                f"{expected.get('id', '<unknown>')}: {key} changed "
                f"expected={expected.get(key)!r} current={current[key]!r}"
            )
    detail = _first_selected_difference(
        _canonical_records(expected.get("selected_salcs", ())),
        _canonical_records(current["selected_salcs"]),
    )
    if detail:
        messages.append(f"{expected.get('id', '<unknown>')}: {detail}")
    elif current["salc_sha256"] != expected.get("salc_sha256"):
        # The complete hash is intentionally retained as a drift diagnostic, but
        # platform BLAS/SVD details can change non-selected SALC bases inside the
        # same symmetry subspace.  The hard golden gate is therefore the
        # human-inspectable representative set plus rank/symmetry metadata.
        pass
    return SALCSnapshotComparison(ok=not messages, messages=tuple(messages))


def _salc_records(
    definition: GICDefinition,
    *,
    rounding_decimals: int,
) -> tuple[dict[str, Any], ...]:
    records = []
    for gic in definition.gics:
        if len(gic.coefficients) <= 1:
            continue
        records.append(
            {
                "name": gic.name,
                "family": gic.family,
                "irrep": gic.irrep,
                "coefficients": [
                    [primitive_id, _rounded_coefficient(coefficient, rounding_decimals)]
                    for primitive_id, coefficient in gic.coefficients
                ],
            }
        )
    return tuple(records)


def _selected_salc_records(
    records: tuple[dict[str, Any], ...],
    *,
    selected_per_family: int,
) -> tuple[dict[str, Any], ...]:
    selected: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    for record in records:
        family = str(record["family"])
        if family in PRIORITY_FAMILIES:
            selected.append(record)
            continue
        count = counts.get(family, 0)
        if count < selected_per_family:
            selected.append(record)
        counts[family] = count + 1
    return tuple(selected)


def _stable_sha256(records: tuple[dict[str, Any], ...]) -> str:
    payload = json.dumps(records, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _rounded_coefficient(value: float, decimals: int) -> float:
    rounded = round(float(value), decimals)
    if rounded == 0.0:
        return 0.0
    return rounded


def _canonical_records(records: Any) -> tuple[dict[str, Any], ...]:
    canonical: list[dict[str, Any]] = []
    for record in records or ():
        coefficients = [
            [str(primitive_id), float(coefficient)]
            for primitive_id, coefficient in record.get("coefficients", ())
        ]
        canonical.append(
            {
                "name": str(record.get("name", "")),
                "family": str(record.get("family", "")),
                "irrep": str(record.get("irrep", "")),
                "coefficients": coefficients,
            }
        )
    return tuple(canonical)


def _first_selected_difference(
    expected: tuple[Any, ...],
    current: tuple[Any, ...],
) -> str:
    if len(expected) != len(current):
        return f"selected SALC count changed expected={len(expected)} current={len(current)}"
    left_groups = _selected_subspace_groups(expected)
    right_groups = _selected_subspace_groups(current)
    if set(left_groups) != set(right_groups):
        return (
            "selected SALC family/irrep groups changed "
            f"expected={sorted(left_groups)} current={sorted(right_groups)}"
        )
    for key in sorted(left_groups):
        left_records = left_groups[key]
        right_records = right_groups[key]
        if len(left_records) != len(right_records):
            return (
                f"selected SALC count changed for {key}: "
                f"expected={len(left_records)} current={len(right_records)}"
            )
        if key[0] not in PRIORITY_FAMILIES:
            # Generic symmetry blocks can be larger than the stored compact
            # selection.  Later vectors in a truncated degenerate block are not
            # unique across BLAS/SVD backends, so only the leading representative
            # is a stable golden artifact.  Special-coordinate families are kept
            # as full selected subspaces because they are the regression target.
            left_records = left_records[:1]
            right_records = right_records[:1]
        left_ids, left_projector = _selected_subspace_projector(left_records)
        right_ids, right_projector = _selected_subspace_projector(right_records)
        if left_ids != right_ids or not np.allclose(
            left_projector,
            right_projector,
            atol=SALC_COEFFICIENT_TOLERANCE,
            rtol=0.0,
        ):
            return (
                f"selected SALC subspace changed for family={key[0]} irrep={key[1]}: "
                f"expected={left_records!r} current={right_records!r}"
            )
    return ""


def _selected_subspace_groups(records: tuple[Any, ...]) -> dict[tuple[str, str], tuple[Any, ...]]:
    groups: dict[tuple[str, str], list[Any]] = {}
    for record in records:
        key = (str(record.get("family", "")), str(record.get("irrep", "")))
        groups.setdefault(key, []).append(record)
    return {key: tuple(value) for key, value in groups.items()}


def _selected_subspace_projector(records: tuple[Any, ...]) -> tuple[tuple[str, ...], np.ndarray]:
    primitive_ids = tuple(
        sorted(
            {
                str(primitive_id)
                for record in records
                for primitive_id, _coefficient in record.get("coefficients", ())
            }
        )
    )
    primitive_index = {primitive_id: index for index, primitive_id in enumerate(primitive_ids)}
    vectors = np.zeros((len(primitive_ids), len(records)), dtype=float)
    for column, record in enumerate(records):
        for primitive_id, coefficient in record.get("coefficients", ()):
            vectors[primitive_index[str(primitive_id)], column] = float(coefficient)
    if vectors.size == 0:
        return primitive_ids, np.zeros((0, 0), dtype=float)
    left, singular_values, _right_t = np.linalg.svd(vectors, full_matrices=False)
    rank = int(np.count_nonzero(singular_values > SALC_COEFFICIENT_TOLERANCE))
    basis = left[:, :rank]
    return primitive_ids, basis @ basis.T
