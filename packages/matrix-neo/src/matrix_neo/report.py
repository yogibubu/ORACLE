"""Human-readable GICForge reports."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from .definition import (
    GICDefinition,
    read_gic_definition_from_xyzin,
    total_symmetric_gic_names,
)
from .policy import (
    FRAGMENT_MODE_PSEUDO_BONDS,
    FRAGMENT_MODE_SPECIAL_COORDINATES,
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
    selected = diagnostics.selected if diagnostics else tuple(
        primitive.identifier for primitive in definition.primitives
    )
    selected_by_family = (
        diagnostics.selected_by_family
        if diagnostics and diagnostics.selected_by_family
        else _family_count_tokens(family_counts)
    )
    primitive_by_id = {primitive.identifier: primitive for primitive in definition.primitives}

    lines = [
        "ORACLE GICForge Report",
        "======================",
        "",
        f"Backend: {definition.backend}",
        f"Point group: {definition.point_group}",
        f"Symmetry group: {definition.point_group}",
        f"Totally symmetric GICs: {_list_or_none(total_symmetric_gic_names(definition))}",
        f"Symmetrize requested: {definition.symmetrize}",
        f"Target rank: {definition.target_rank}",
        f"Target rank rationale: {_target_rank_rationale(definition)}",
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
        "Fragment Mode Policy",
        "--------------------",
        *_fragment_policy_lines(definition),
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

    lines.extend(["", "Protected Coordinates", "---------------------"])
    protected = [
        primitive
        for primitive in definition.primitives
        if primitive.reduction_class == SPECIAL_REDUCTION_CLASS
    ]
    if protected:
        for primitive in protected:
            lines.append(_primitive_summary(primitive))
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

    symmetry = definition.symmetry_diagnostics
    lines.extend(["", "Symmetrization Diagnostics", "---------------------------"])
    if symmetry:
        lines.append(f"Method: {symmetry.method}")
        lines.append(f"Policy: {symmetry.policy}")
        lines.append(f"Status: {symmetry.status}")
        lines.append(f"Symmetry group: {symmetry.symmetry_group}")
        lines.append(f"Total irrep: {symmetry.total_symmetric_irrep}")
        lines.append("Total GICs: " + _list_or_none(symmetry.total_symmetric_gics))
        lines.append(f"Groups: {len(symmetry.groups)}")
        for group in symmetry.groups:
            lines.append(
                f"{group.block} {group.signature}: "
                f"{','.join(group.source_gics)} -> {','.join(group.output_gics)}"
            )
    else:
        lines.append("NONE")

    lines.extend(["", "SALC Coefficients", "-----------------"])
    coefficient_lines = _salc_coefficient_lines(definition)
    lines.extend(coefficient_lines or ["NONE"])

    lines.extend(["", "Reduction Diagnostics", "---------------------"])
    lines.append("Selected: " + _list_or_none(selected))
    lines.append("Selected by family: " + _list_or_none(selected_by_family))
    lines.append("Skipped singular: " + _list_or_none(skipped_singular))
    lines.append("Skipped dependent: " + _list_or_none(skipped_dependent))
    lines.append(
        "Skipped singular details: "
        + _list_or_none(diagnostics.skipped_singular_details if diagnostics else ())
    )
    lines.append(
        "Skipped dependent details: "
        + _list_or_none(diagnostics.skipped_dependent_details if diagnostics else ())
    )

    lines.extend(["", "Frozen GICs", "-----------"])
    if definition.gics:
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


def _family_count_tokens(family_counts: Counter[str]) -> tuple[str, ...]:
    return tuple(f"{family}:{family_counts[family]}" for family in sorted(family_counts))


def _target_rank_rationale(definition: GICDefinition) -> str:
    natoms = len(definition.reference_coordinates_angstrom)
    nonlinear = 3 * natoms - 6
    linear = 3 * natoms - 5
    if definition.target_rank == nonlinear:
        return f"3N-6 non-linear vibrational rank for N={natoms}"
    if definition.target_rank == linear:
        return f"3N-5 linear vibrational rank for N={natoms}"
    return f"contract-specified rank {definition.target_rank} for N={natoms}"


def _fragment_policy_lines(definition: GICDefinition) -> list[str]:
    mode = definition.fragment_mode
    lines = [f"Mode: {mode}"]
    if mode == FRAGMENT_MODE_SPECIAL_COORDINATES:
        lines.extend(
            [
                "Policy: automatic weak-complex default; keep fragments as protected bodies.",
                "Coordinates: fragment-center distances, center-atom distances, translations and orientations.",
                "Rationale: preserve inter-fragment rigid-body motion before ordinary valence pruning.",
            ]
        )
    elif mode == FRAGMENT_MODE_PSEUDO_BONDS:
        kinds = definition.pseudo_bond_kinds or (
            ("INTERFRAGMENT_CLOSEST",) * len(definition.pseudo_bonds)
        )
        contacts = tuple(
            f"{left}-{right}:{kind}"
            for (left, right), kind in zip(definition.pseudo_bonds, kinds)
        )
        lines.extend(
            [
                "Policy: explicit graph-joining mode; do not build protected fragment coordinates.",
                "Selection: Merlino/BDPCS3 H-bonds first, closest inter-fragment fallback otherwise.",
                "Pseudo-bonds: " + _list_or_none(contacts),
                "Rationale: use ordinary stretch/bend/torsion coordinates on the augmented construction graph.",
            ]
        )
    else:
        lines.append("Policy: no built fragments were consumed by this GIC definition.")
    return lines


def _primitive_summary(primitive) -> str:
    atoms = ",".join(str(atom) for atom in primitive.atoms) or "NONE"
    refs = ",".join(primitive.refs) or "NONE"
    return (
        f"{primitive.identifier} {primitive.name} family={primitive.family} "
        f"class={primitive.reduction_class} atoms={atoms} refs={refs}"
    )


def _salc_coefficient_lines(definition: GICDefinition) -> list[str]:
    lines: list[str] = []
    for gic in definition.gics:
        if len(gic.coefficients) <= 1:
            continue
        terms = ",".join(
            f"{primitive_id}:{coefficient:+.12g}"
            for primitive_id, coefficient in gic.coefficients
        )
        norm2 = sum(float(coefficient) ** 2 for _primitive_id, coefficient in gic.coefficients)
        lines.append(
            f"{gic.name} irrep={gic.irrep} family={gic.family} "
            f"norm2={norm2:.12g} coeffs={terms}"
        )
    return lines
