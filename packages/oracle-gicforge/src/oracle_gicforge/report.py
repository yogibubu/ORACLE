"""Human-readable GICForge reports."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from .definition import GICDefinition, read_gic_definition_from_xyzin
from .policy import (
    REDUCTION_POLICY,
    SPECIAL_REDUCTION_CLASS,
    primitive_reduction_class,
)
from .symmetry import gic_symmetry_source_blocks


def gic_report_lines(definition: GICDefinition) -> list[str]:
    family_counts = Counter(primitive.family for primitive in definition.primitives)
    protected_count = sum(
        1
        for primitive in definition.primitives
        if primitive.reduction_class == SPECIAL_REDUCTION_CLASS
    )
    diagnostics = definition.reduction_diagnostics
    skipped_singular = diagnostics.skipped_singular if diagnostics else ()
    skipped_dependent = diagnostics.skipped_dependent if diagnostics else ()
    rank_method = diagnostics.rank_method if diagnostics else "UNKNOWN"
    reduction_policy = diagnostics.reduction_policy if diagnostics else REDUCTION_POLICY

    lines = [
        "ORACLE GICForge Report",
        "======================",
        "",
        f"Backend: {definition.backend}",
        f"Point group: {definition.point_group}",
        f"Symmetrize requested: {definition.symmetrize}",
        f"Target rank: {definition.target_rank}",
        f"Final rank: {definition.rank}",
        f"Candidate count: {definition.candidate_count}",
        f"Selected primitive count: {len(definition.primitives)}",
        f"Frozen GIC count: {len(definition.gics)}",
        f"Protected selected count: {protected_count}",
        f"Rank method: {rank_method}",
        f"Reduction policy: {reduction_policy}",
        f"Skipped singular/zero rows: {len(skipped_singular)}",
        f"Skipped dependent rows: {len(skipped_dependent)}",
        "",
        "Selected Families",
        "-----------------",
    ]
    if family_counts:
        for family in sorted(family_counts):
            lines.append(
                f"{family}: {family_counts[family]} "
                f"({primitive_reduction_class(family)})"
            )
    else:
        lines.append("NONE")

    lines.extend(["", "Symmetry Source Blocks", "----------------------"])
    blocks = gic_symmetry_source_blocks(definition)
    if blocks:
        for block in blocks:
            lines.append(
                f"{block.block}: family={block.family} "
                f"class={block.reduction_class} count={len(block.gic_names)}"
            )
    else:
        lines.append("NONE")

    lines.extend(["", "Reduction Diagnostics", "---------------------"])
    lines.append("Selected: " + _list_or_none(diagnostics.selected if diagnostics else ()))
    lines.append("Skipped singular: " + _list_or_none(skipped_singular))
    lines.append("Skipped dependent: " + _list_or_none(skipped_dependent))

    lines.extend(["", "Frozen GICs", "-----------"])
    if definition.gics:
        primitive_by_id = {primitive.identifier: primitive for primitive in definition.primitives}
        for gic in definition.gics:
            primitive = primitive_by_id.get(gic.primitive_id)
            reduction_class = primitive.reduction_class if primitive else "UNKNOWN"
            lines.append(
                f"{gic.identifier} {gic.name} family={gic.family} "
                f"class={reduction_class} irrep={gic.irrep} primitive={gic.primitive_id}"
            )
    else:
        lines.append("NONE")
    return lines


def gic_report_from_xyzin(path: Path) -> list[str]:
    return gic_report_lines(read_gic_definition_from_xyzin(Path(path)))


def write_gic_report(path: Path, output: Path) -> Path:
    target = Path(output)
    target.write_text("\n".join(gic_report_from_xyzin(Path(path))) + "\n", encoding="utf-8")
    return target


def _list_or_none(values: tuple[str, ...]) -> str:
    return ", ".join(values) if values else "NONE"
