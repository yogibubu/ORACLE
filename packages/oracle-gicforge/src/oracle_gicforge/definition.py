from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np

from oracle_chem import read_enriched_xyz
from oracle_core import read_sectioned_lines, replace_section, section_content

from .contracts import (
    ORACLE_XYZ_GIC_SCHEMA,
    ORACLE_XYZ_SYCART_SCHEMA,
    GICForgeContractError,
    validate_gicforge_prerequisites,
)
from .policy import (
    B_MATRIX_BACKEND,
    DIAGNOSTIC_FINITE_DIFFERENCE_STEP,
    GIC_BACKEND,
    LINEAR_ANGLE_DEGREES,
    LOCAL_SYMMETRIZATION_METHOD,
    ORDINARY_REDUCTION_CLASS,
    POINT_GROUP_PROJECTOR_METHOD,
    PRIMITIVE_FAMILY_ORDER,
    PROJECTOR_SYMMETRIZATION_POLICY,
    RANK_METHOD,
    RANK_TOLERANCE,
    REDUCTION_POLICY,
    SPECIAL_REDUCTION_CLASS,
    SYMMETRIZATION_POLICY,
    SYCART_BACKEND,
    primitive_prefix,
    primitive_reduction_class,
    primitive_symmetry_block,
)
from .symmetry_labels import (
    irrep_characters_for_operations,
    irrep_name_prefix,
    is_total_symmetric_irrep,
    non_total_irrep_sequence,
    total_symmetric_irrep,
)


@dataclass(frozen=True)
class GICPrimitive:
    identifier: str
    name: str
    family: str
    function: str
    atoms: tuple[int, ...]
    mode: int = 0
    ref_atoms: tuple[int, ...] = ()
    refs: tuple[str, ...] = ()
    frame_atoms: tuple[int, ...] = ()
    ref_frame_atoms: tuple[int, ...] = ()

    def gaussian_expression(self) -> str:
        if not self.is_gaussian_native:
            return "NONE"
        atoms = ",".join(str(atom) for atom in self.atoms)
        if self.function == "L":
            return f"L({atoms},0,{self.mode})"
        return f"{self.function}({atoms})"

    @property
    def is_gaussian_native(self) -> bool:
        return self.function in {"R", "A", "L", "D", "U"}

    @property
    def reduction_class(self) -> str:
        return primitive_reduction_class(self.family)


@dataclass(frozen=True)
class FrozenGIC:
    identifier: str
    name: str
    family: str
    irrep: str
    primitive_id: str
    gaussian_expression: str
    coefficients: tuple[tuple[str, float], ...] = ()


@dataclass(frozen=True)
class GICReductionDiagnostics:
    rank_method: str
    reduction_policy: str
    selected: tuple[str, ...] = ()
    skipped_singular: tuple[str, ...] = ()
    skipped_dependent: tuple[str, ...] = ()


@dataclass(frozen=True)
class GICSymmetrizedGroup:
    block: str
    family: str
    signature: str
    source_gics: tuple[str, ...]
    output_gics: tuple[str, ...]


@dataclass(frozen=True)
class GICSymmetrizationDiagnostics:
    method: str
    policy: str
    status: str
    point_group: str
    symmetry_group: str
    total_symmetric_irrep: str
    total_symmetric_gics: tuple[str, ...] = ()
    groups: tuple[GICSymmetrizedGroup, ...] = ()


@dataclass(frozen=True)
class GICPointGroupOperation:
    label: str
    rotation: tuple[tuple[float, float, float], ...]
    permutation: tuple[int, ...]


@dataclass(frozen=True)
class GICDefinition:
    backend: str
    point_group: str
    symmetrize: bool
    target_rank: int
    rank: int
    candidate_count: int
    reference_coordinates_angstrom: tuple[tuple[float, float, float], ...]
    primitives: tuple[GICPrimitive, ...]
    gics: tuple[FrozenGIC, ...]
    reduction_diagnostics: GICReductionDiagnostics | None = None
    symmetry_diagnostics: GICSymmetrizationDiagnostics | None = None


@dataclass(frozen=True)
class GICBMatrix:
    backend: str
    coordinate_labels: tuple[str, ...]
    coordinate_names: tuple[str, ...]
    irreps: tuple[str, ...]
    cartesian_columns: tuple[str, ...]
    rows: tuple[tuple[float, ...], ...]


@dataclass(frozen=True)
class SYCartDefinition:
    backend: str
    point_group: str
    target_rank: int
    vectors: tuple[tuple[float, ...], ...]


def build_gic_definition_from_xyzin(
    path: Path,
    *,
    symmetrize: bool = False,
    rank_tolerance: float = RANK_TOLERANCE,
) -> GICDefinition:
    """Build a frozen ORACLE GIC definition from saved xyzin state."""
    definition, atom_symbols, operations = construct_gic_definition_from_xyzin(
        path,
        rank_tolerance=rank_tolerance,
    )
    if not symmetrize:
        return definition
    return symmetrize_gic_definition(
        definition,
        atom_symbols=atom_symbols,
        symmetry_operations=operations,
    )


def construct_gic_definition_from_xyzin(
    path: Path,
    *,
    rank_tolerance: float = RANK_TOLERANCE,
) -> tuple[GICDefinition, tuple[str, ...], tuple[GICPointGroupOperation, ...]]:
    """Construct and reduce GICs without applying symmetry adaptation."""
    target = Path(path)
    validate_gicforge_prerequisites(target)
    lines = read_sectioned_lines(target)
    geometry = read_enriched_xyz(target)
    coords = np.asarray(geometry.coordinates_angstrom, dtype=float)
    point_group = _point_group(lines)
    symmetry_operations = _symmetry_operations(lines)

    bonds = _topology_bonds(lines, natoms=geometry.natoms)
    rings = _topology_rings(lines, natoms=geometry.natoms)
    candidates = _primitive_candidates(
        bonds,
        rings=rings,
        coords=coords,
        natoms=geometry.natoms,
        fragment_records=_fragment_records(target),
        interaction_centers=_interaction_center_definition(target),
    )
    target_rank = _vibrational_rank(coords)
    selected, rank, reduction_diagnostics = _select_ranked_primitives_with_diagnostics(
        candidates,
        coords,
        target_rank=target_rank,
        rank_tolerance=rank_tolerance,
    )
    if rank != target_rank:
        raise GICForgeContractError(
            "insufficient independent primitive coordinates: "
            f"need {target_rank}, selected rank {rank} from {len(candidates)} candidates"
        )

    gics = tuple(
        FrozenGIC(
            identifier=f"GIC{idx:03d}",
            name=primitive.name,
            family=primitive.family,
            irrep="A" if point_group.upper() == "C1" else "UNASSIGNED",
            primitive_id=primitive.identifier,
            gaussian_expression=primitive.gaussian_expression(),
            coefficients=((primitive.identifier, 1.0),),
        )
        for idx, primitive in enumerate(selected, start=1)
    )
    definition = GICDefinition(
        backend=GIC_BACKEND,
        point_group=point_group,
        symmetrize=False,
        target_rank=target_rank,
        rank=rank,
        candidate_count=len(candidates),
        reference_coordinates_angstrom=tuple(
            tuple(float(value) for value in row) for row in coords
        ),
        primitives=tuple(selected),
        gics=gics,
        reduction_diagnostics=reduction_diagnostics,
        symmetry_diagnostics=_empty_symmetry_diagnostics(point_group, requested=False),
    )
    return definition, tuple(geometry.atoms), symmetry_operations


def symmetrize_gic_definition(
    definition: GICDefinition,
    *,
    atom_symbols: tuple[str, ...],
    symmetry_operations: tuple[GICPointGroupOperation, ...] = (),
) -> GICDefinition:
    """Apply the frozen GIC symmetrization utility to a reduced definition."""
    gics, symmetry_diagnostics = _apply_local_symmetrization(
        definition.gics,
        definition.primitives,
        atom_symbols=tuple(atom_symbols),
        point_group=definition.point_group,
        requested=True,
        symmetry_operations=tuple(symmetry_operations),
    )
    return GICDefinition(
        backend=definition.backend,
        point_group=definition.point_group,
        symmetrize=True,
        target_rank=definition.target_rank,
        rank=definition.rank,
        candidate_count=definition.candidate_count,
        reference_coordinates_angstrom=definition.reference_coordinates_angstrom,
        primitives=definition.primitives,
        gics=gics,
        reduction_diagnostics=definition.reduction_diagnostics,
        symmetry_diagnostics=symmetry_diagnostics,
    )


def build_sycart_definition_from_xyzin(path: Path) -> SYCartDefinition:
    """Build an external-mode-free Cartesian basis for #SYCART."""
    target = Path(path)
    validate_gicforge_prerequisites(target)
    lines = read_sectioned_lines(target)
    geometry = read_enriched_xyz(target)
    coords = np.asarray(geometry.coordinates_angstrom, dtype=float)
    target_rank = _vibrational_rank(coords)
    basis = _cartesian_vibrational_basis(coords, target_rank=target_rank)
    return SYCartDefinition(
        backend=SYCART_BACKEND,
        point_group=_point_group(lines),
        target_rank=target_rank,
        vectors=tuple(tuple(float(value) for value in row) for row in basis.T),
    )


def gic_definition_section_lines(definition: GICDefinition) -> list[str]:
    lines = [
        f"SCHEMA {ORACLE_XYZ_GIC_SCHEMA}",
        "STATUS BUILT",
        "DEPENDENCIES VALIDATION=oracle.xyz.validation.v1 "
        "TOPOLOGY=oracle.xyz.topology.v1 SYNTHONS=oracle.xyz.synthons.v1 "
        "SYMMETRY=oracle.xyz.symmetry.v1",
        "INDEXING ATOMS=ONE_BASED",
        f"BACKEND {definition.backend}",
        f"POINT_GROUP {definition.point_group}",
        f"SYMMETRY_GROUP {definition.point_group}",
        f"TOTAL_SYMMETRIC_IRREP {total_symmetric_irrep(definition.point_group)}",
        f"TOTAL_SYMMETRIC_GIC_COUNT {len(total_symmetric_gics(definition))}",
        f"TOTAL_SYMMETRIC_GICS {_csv_or_none(total_symmetric_gic_names(definition))}",
        f"SYMMETRIZE {_bool_text(definition.symmetrize)}",
        f"SYMMETRY_MODE {_symmetry_mode(definition)}",
        f"TARGET_RANK {definition.target_rank}",
        f"RANK {definition.rank}",
        f"CANDIDATE_COUNT {definition.candidate_count}",
        f"PRIMITIVE_COUNT {len(definition.primitives)}",
        f"GIC_COUNT {len(definition.gics)}",
        f"PROTECTED_GIC_COUNT {_protected_gic_count(definition.primitives)}",
        f"SKIPPED_SINGULAR_COUNT {_skipped_singular_count(definition)}",
        f"SKIPPED_DEPENDENT_COUNT {_skipped_dependent_count(definition)}",
        f"RANK_METHOD {RANK_METHOD}",
        f"REDUCTION_POLICY {REDUCTION_POLICY}",
        "B_MATRIX_DERIVATIVE_MODE ANALYTIC",
        f"RANK_TOLERANCE {RANK_TOLERANCE:.12g}",
        "[PRIMITIVES]",
    ]
    if definition.primitives:
        lines.extend(_primitive_line(primitive) for primitive in definition.primitives)
    else:
        lines.append("NONE")
    lines.append("[FROZEN_GICS]")
    if definition.gics:
        lines.extend(_frozen_gic_line(gic) for gic in definition.gics)
    else:
        lines.append("NONE")
    lines.append("[REDUCTION_DIAGNOSTICS]")
    lines.extend(_reduction_diagnostics_lines(definition))
    lines.append("[SYMMETRY_DIAGNOSTICS]")
    lines.extend(_symmetry_diagnostics_lines(definition))
    lines.append("[GAUSSIAN_GIC]")
    gaussian_gics = _gaussian_gic_block_lines(definition)
    if gaussian_gics:
        lines.extend(gaussian_gics)
    else:
        lines.append("NONE")
    return lines


def sycart_definition_section_lines(definition: SYCartDefinition) -> list[str]:
    lines = [
        f"SCHEMA {ORACLE_XYZ_SYCART_SCHEMA}",
        "STATUS BUILT",
        "DEPENDENCIES VALIDATION=oracle.xyz.validation.v1 GIC=oracle.xyz.gic.v1 "
        "SYMMETRY=oracle.xyz.symmetry.v1",
        "INDEXING ATOMS=ONE_BASED",
        f"BACKEND {definition.backend}",
        f"POINT_GROUP {definition.point_group}",
        (
            "SYMMETRY_MODE IDENTITY_C1"
            if definition.point_group.upper() == "C1"
            else "SYMMETRY_MODE UNSYMMETRIZED"
        ),
        f"TARGET_RANK {definition.target_rank}",
        f"COORD_COUNT {len(definition.vectors)}",
        "[SYCART]",
    ]
    if definition.vectors:
        lines.extend(
            f"SYC{idx:03d} IRREP=A " + _sycart_components(vector)
            for idx, vector in enumerate(definition.vectors, start=1)
        )
    else:
        lines.append("NONE")
    return lines


def write_gicforge_build_sections(
    path: Path,
    *,
    symmetrize: bool = False,
    sycart: bool = False,
) -> GICDefinition:
    target = Path(path)
    definition = build_gic_definition_from_xyzin(target, symmetrize=symmetrize)
    replace_section(target, "GIC", gic_definition_section_lines(definition))
    if sycart:
        sycart_definition = build_sycart_definition_from_xyzin(target)
        replace_section(target, "SYCART", sycart_definition_section_lines(sycart_definition))
    return definition


def build_gic_b_matrix(
    definition: GICDefinition,
    *,
    coordinates_angstrom: tuple[tuple[float, float, float], ...] | np.ndarray | None = None,
) -> GICBMatrix:
    """Evaluate the Wilson B matrix for a frozen GIC definition."""
    coords = (
        np.asarray(definition.reference_coordinates_angstrom, dtype=float)
        if coordinates_angstrom is None
        else np.asarray(coordinates_angstrom, dtype=float)
    )
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise GICForgeContractError("B-matrix coordinates must have shape (natoms, 3)")
    primitive_by_id = {primitive.identifier: primitive for primitive in definition.primitives}
    rows: list[tuple[float, ...]] = []
    for gic in definition.gics:
        row = np.zeros(coords.size, dtype=float)
        coefficients = gic.coefficients or ((gic.primitive_id, 1.0),)
        for primitive_id, coefficient in coefficients:
            primitive = primitive_by_id.get(primitive_id)
            if primitive is None:
                raise GICForgeContractError(
                    f"unknown primitive {primitive_id!r} in frozen GIC {gic.identifier}"
                )
            row += float(coefficient) * _analytic_b_row(primitive, coords)
        if not np.all(np.isfinite(row)):
            raise GICForgeContractError(
                f"non-finite B-matrix row for frozen GIC {gic.identifier}"
            )
        rows.append(tuple(float(value) for value in row))
    return GICBMatrix(
        backend=B_MATRIX_BACKEND,
        coordinate_labels=tuple(gic.identifier for gic in definition.gics),
        coordinate_names=tuple(gic.name for gic in definition.gics),
        irreps=tuple(gic.irrep for gic in definition.gics),
        cartesian_columns=_cartesian_column_labels(coords.shape[0]),
        rows=tuple(rows),
    )


def build_gic_b_matrix_from_xyzin(path: Path) -> GICBMatrix:
    """Evaluate the B matrix from the frozen #GIC section of an enriched XYZ."""
    target = Path(path)
    definition = read_gic_definition_from_xyzin(target)
    geometry = read_enriched_xyz(target)
    return build_gic_b_matrix(
        definition,
        coordinates_angstrom=geometry.coordinates_angstrom,
    )


def total_symmetric_gics(definition: GICDefinition) -> tuple[FrozenGIC, ...]:
    """Return frozen GICs active in symmetry-preserving optimization/fitting."""
    return tuple(
        gic
        for gic in definition.gics
        if is_total_symmetric_irrep(definition.point_group, gic.irrep)
    )


def total_symmetric_gic_names(definition: GICDefinition) -> tuple[str, ...]:
    return tuple(gic.name for gic in total_symmetric_gics(definition))


def read_gic_definition_from_xyzin(path: Path) -> GICDefinition:
    """Read a frozen ORACLE GIC definition without regenerating coordinates."""
    target = Path(path)
    lines = read_sectioned_lines(target)
    geometry = read_enriched_xyz(target)
    section = section_content(lines, "GIC")
    if not section:
        raise GICForgeContractError("missing #GIC section")
    if section[0].strip() != f"SCHEMA {ORACLE_XYZ_GIC_SCHEMA}":
        raise GICForgeContractError("invalid #GIC schema")
    status = _section_value(section, "STATUS")
    if (status or "").upper() != "BUILT":
        raise GICForgeContractError(f"#GIC status must be BUILT; found {status or 'UNKNOWN'}")
    primitives = tuple(
        _parse_primitive_line(line)
        for line in _subsection(section, "PRIMITIVES")
        if line.strip() and line.strip().upper() != "NONE"
    )
    gics = tuple(
        _parse_frozen_gic_line(line)
        for line in _subsection(section, "FROZEN_GICS")
        if line.strip() and line.strip().upper() != "NONE"
    )
    diagnostics = _parse_reduction_diagnostics(
        section,
        selected=tuple(p.identifier for p in primitives),
    )
    symmetry_diagnostics = _parse_symmetry_diagnostics(section)
    return GICDefinition(
        backend=_section_value(section, "BACKEND") or GIC_BACKEND,
        point_group=_section_value(section, "POINT_GROUP") or _point_group(lines),
        symmetrize=_parse_bool(_section_value(section, "SYMMETRIZE")),
        target_rank=_parse_int(_section_value(section, "TARGET_RANK")),
        rank=_parse_int(_section_value(section, "RANK")),
        candidate_count=_parse_int(_section_value(section, "CANDIDATE_COUNT")),
        reference_coordinates_angstrom=tuple(
            tuple(float(value) for value in row)
            for row in geometry.coordinates_angstrom
        ),
        primitives=primitives,
        gics=gics,
        reduction_diagnostics=diagnostics,
        symmetry_diagnostics=symmetry_diagnostics,
    )


def gic_b_matrix_lines(matrix: GICBMatrix) -> list[str]:
    """Serialize a GIC B matrix in a compact machine-readable text format."""
    lines = [
        "SCHEMA oracle.gic.bmatrix.v1",
        f"BACKEND {matrix.backend}",
        "UNITS MIXED_GIC_PER_ANGSTROM",
        "DERIVATIVE_MODE ANALYTIC",
        f"ROW_COUNT {len(matrix.rows)}",
        f"COLUMN_COUNT {len(matrix.cartesian_columns)}",
        "[COLUMNS]",
        " ".join(matrix.cartesian_columns) if matrix.cartesian_columns else "NONE",
        "[ROWS]",
    ]
    if not matrix.rows:
        lines.append("NONE")
        return lines
    for label, name, irrep, row in zip(
        matrix.coordinate_labels,
        matrix.coordinate_names,
        matrix.irreps,
        matrix.rows,
    ):
        values = ",".join(f"{value:.12g}" for value in row)
        lines.append(f"{label} NAME={name} IRREP={irrep} VALUES={values}")
    return lines


def write_gic_b_matrix(path: Path, output: Path) -> GICBMatrix:
    matrix = build_gic_b_matrix_from_xyzin(path)
    target = Path(output)
    target.write_text("\n".join(gic_b_matrix_lines(matrix)) + "\n", encoding="utf-8")
    return matrix


def gaussian_gic_lines_from_xyzin(path: Path) -> list[str]:
    gic = section_content(read_sectioned_lines(Path(path)), "GIC")
    block = _subsection(gic, "GAUSSIAN_GIC")
    return [line for line in block if line.strip() and line.strip().upper() != "NONE"]


def _topology_bonds(lines: list[str], *, natoms: int) -> tuple[tuple[int, int], ...]:
    topology = section_content(lines, "TOPOLOGY")
    expected = "SCHEMA oracle.xyz.topology.v1"
    if not topology or topology[0].strip() != expected:
        raise GICForgeContractError("missing valid #TOPOLOGY section")
    bond_lines = _subsection(topology, "BONDS")
    if not bond_lines or any(line.strip().upper() == "NONE" for line in bond_lines):
        raise GICForgeContractError("#TOPOLOGY contains no bonds")
    bonds: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for line in bond_lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            i, j = int(parts[0]), int(parts[1])
        except ValueError as exc:
            raise GICForgeContractError(f"invalid #TOPOLOGY bond line: {line}") from exc
        if i == j or i < 1 or j < 1 or i > natoms or j > natoms:
            raise GICForgeContractError(f"invalid #TOPOLOGY bond indexes: {line}")
        bond = tuple(sorted((i, j)))
        if bond not in seen:
            seen.add(bond)
            bonds.append(bond)
    return tuple(sorted(bonds))


def _topology_rings(lines: list[str], *, natoms: int) -> tuple[tuple[int, tuple[int, ...]], ...]:
    topology = section_content(lines, "TOPOLOGY")
    expected = "SCHEMA oracle.xyz.topology.v1"
    if not topology or topology[0].strip() != expected:
        raise GICForgeContractError("missing valid #TOPOLOGY section")
    ring_lines = _subsection(topology, "RINGS")
    rings: list[tuple[int, tuple[int, ...]]] = []
    for line in ring_lines:
        if line.strip().upper() == "NONE":
            continue
        parts = line.replace(",", " ").replace("[", " ").replace("]", " ").split()
        if not parts:
            continue
        try:
            ring_index = int(parts[0])
        except ValueError:
            continue
        atoms: list[int] = []
        reading_atoms = False
        for part in parts[1:]:
            token = part.strip()
            if token.upper().startswith("ATOMS="):
                reading_atoms = True
                token = token.split("=", 1)[1]
            elif "=" in token and reading_atoms:
                break
            if not reading_atoms or not token:
                continue
            try:
                atoms.append(int(token))
            except ValueError as exc:
                raise GICForgeContractError(f"invalid #TOPOLOGY ring line: {line}") from exc
        if len(atoms) < 3:
            continue
        if any(atom < 1 or atom > natoms for atom in atoms):
            raise GICForgeContractError(f"invalid #TOPOLOGY ring atom indexes: {line}")
        rings.append((ring_index, tuple(dict.fromkeys(atoms))))
    return tuple(rings)


def _primitive_candidates(
    bonds: tuple[tuple[int, int], ...],
    *,
    rings: tuple[tuple[int, tuple[int, ...]], ...] = (),
    coords: np.ndarray,
    natoms: int,
    fragment_records: tuple[object, ...] = (),
    interaction_centers: object | None = None,
) -> tuple[GICPrimitive, ...]:
    adjacency = _adjacency(bonds, natoms=natoms)
    counters: dict[str, int] = {family: 0 for family in PRIMITIVE_FAMILY_ORDER}
    candidates: list[GICPrimitive] = []

    for i, j in bonds:
        candidates.append(_make_primitive("STRETCH", "R", (i, j), counters))

    candidates.extend(
        _fragment_primitive_candidates(fragment_records, coords=coords, counters=counters)
    )
    candidates.extend(
        _interaction_center_primitive_candidates(interaction_centers, counters=counters)
    )

    for center in range(1, natoms + 1):
        for i, k in combinations(sorted(adjacency[center]), 2):
            angle = _angle_value(coords, (i, center, k))
            if np.degrees(angle) >= LINEAR_ANGLE_DEGREES:
                candidates.append(
                    _make_primitive("LINEAR_BEND", "L", (i, center, k), counters, mode=-1)
                )
                candidates.append(
                    _make_primitive("LINEAR_BEND", "L", (i, center, k), counters, mode=-2)
                )
            else:
                family = (
                    "CYCLIC_BEND"
                    if _ring_index_for_atoms((i, center, k), rings) is not None
                    else "BEND"
                )
                candidates.append(_make_primitive(family, "A", (i, center, k), counters))

    seen_torsions: set[tuple[int, int, int, int]] = set()
    butterfly_torsions: list[tuple[str, tuple[int, int, int, int]]] = []
    condensed_torsions: list[tuple[str, tuple[int, int, int, int]]] = []
    ordinary_torsions: list[tuple[str, tuple[int, int, int, int]]] = []
    for j, k in bonds:
        for i in sorted(adjacency[j] - {k}):
            for l in sorted(adjacency[k] - {j}):
                torsion = (i, j, k, l)
                canonical = min(torsion, tuple(reversed(torsion)))
                if canonical in seen_torsions:
                    continue
                seen_torsions.add(canonical)
                family = _torsion_family(canonical, rings)
                if family == "CYCLIC_TORSION":
                    continue
                if family == "BUTTERFLY":
                    butterfly_torsions.append((family, canonical))
                elif family == "CONDENSED_RING_TORSION":
                    condensed_torsions.append((family, canonical))
                else:
                    ordinary_torsions.append((family, canonical))

    for family, torsion in butterfly_torsions:
        candidates.append(_make_primitive(family, "D", torsion, counters))
    candidates.extend(_ring_pucker_component_candidates(rings, counters=counters))
    for family, torsion in condensed_torsions:
        candidates.append(_make_primitive(family, "D", torsion, counters))
    for family, torsion in ordinary_torsions:
        candidates.append(_make_primitive(family, "D", torsion, counters))

    for center in range(1, natoms + 1):
        neighbors = sorted(adjacency[center])
        if len(neighbors) < 3:
            continue
        for n1, n2, n3 in combinations(neighbors, 3):
            candidates.append(
                _make_primitive("OUT_OF_PLANE", "U", (center, n1, n2, n3), counters)
            )

    return tuple(candidates)


def _torsion_family(
    atoms: tuple[int, int, int, int],
    rings: tuple[tuple[int, tuple[int, ...]], ...],
) -> str:
    if _butterfly_torsion(atoms, rings):
        return "BUTTERFLY"
    if _ring_index_for_atoms(atoms, rings) is not None:
        return "CYCLIC_TORSION"
    component = _ring_component_for_atoms(atoms, rings)
    if component is not None and len(component) > 1:
        return "CONDENSED_RING_TORSION"
    return "TORSION"


def _ring_pucker_component_candidates(
    rings: tuple[tuple[int, tuple[int, ...]], ...],
    *,
    counters: dict[str, int],
) -> tuple[GICPrimitive, ...]:
    candidates: list[GICPrimitive] = []
    for _ring_index, ring_atoms in rings:
        for terms in _ring_pucker_component_terms(ring_atoms):
            candidates.append(
                _make_primitive(
                    "RING_PUCKER_COMPONENT",
                    "RPCK",
                    tuple(ring_atoms),
                    counters,
                    refs=tuple(
                        _encode_ring_pucker_term(coefficient, atoms)
                        for coefficient, atoms in terms
                    ),
                )
            )
    return tuple(candidates)


def _ring_pucker_component_terms(
    ring_atoms: tuple[int, ...],
) -> tuple[tuple[tuple[float, tuple[int, int, int, int]], ...], ...]:
    """Return Merlino-style RPck linear combinations for one ordered ring."""
    ncyc = len(ring_atoms)
    if ncyc <= 3:
        return ()
    vnorm = float(np.sqrt(2.0 / float(ncyc)))
    vnorm1 = float(np.sqrt(1.0 / float(ncyc)))
    istart = ncyc
    components: list[tuple[tuple[float, tuple[int, int, int, int]], ...]] = []
    for ivar in range(1, ncyc - 2):
        even = ivar == 2 * (ivar // 2)
        terms: list[tuple[float, tuple[int, int, int, int]]] = []
        for iterm in range(1, ncyc + 1):
            iang1 = _cyclic_index(iterm + istart - 1, ncyc)
            iang2 = _cyclic_index(iterm + istart, ncyc)
            iang3 = _cyclic_index(iterm + istart + 1, ncyc)
            iang4 = _cyclic_index(iterm + istart + 2, ncyc)
            ivar1 = ivar
            if ivar == 1:
                ivar1 = 2
            elif ivar == 4:
                ivar1 = 3
            elif ivar in {5, 6}:
                ivar1 = 4
            value = np.pi * (2.0 * float(ivar1) * float(iterm - 1)) / float(ncyc)
            if even:
                coefficient = vnorm * float(np.sin(value))
            elif ivar < ncyc - 3:
                coefficient = vnorm * float(np.cos(value))
            else:
                coefficient = vnorm1 * float(np.cos(float(iterm - 1) * np.pi))
            if abs(coefficient) <= 1.0e-14:
                coefficient = 0.0
            terms.append(
                (
                    coefficient,
                    (
                        ring_atoms[iang1],
                        ring_atoms[iang2],
                        ring_atoms[iang3],
                        ring_atoms[iang4],
                    ),
                )
            )
        components.append(tuple(terms))
    return tuple(components)


def _cyclic_index(index_1based: int, ncyc: int) -> int:
    while index_1based > ncyc:
        index_1based -= ncyc
    while index_1based <= 0:
        index_1based += ncyc
    return index_1based - 1


def _encode_ring_pucker_term(
    coefficient: float,
    atoms: tuple[int, int, int, int],
) -> str:
    atom_text = "-".join(str(atom) for atom in atoms)
    return f"{float(coefficient):.17g}:{atom_text}"


def _ring_pucker_terms_from_refs(
    primitive: GICPrimitive,
) -> tuple[tuple[float, tuple[int, int, int, int]], ...]:
    terms: list[tuple[float, tuple[int, int, int, int]]] = []
    for ref in primitive.refs:
        if ":" not in ref:
            raise GICForgeContractError(
                f"invalid RPck term {ref!r} in primitive {primitive.identifier}"
            )
        coefficient_text, atom_text = ref.split(":", 1)
        try:
            coefficient = float(coefficient_text)
            atoms = tuple(int(atom) for atom in atom_text.split("-") if atom)
        except ValueError as exc:
            raise GICForgeContractError(
                f"invalid RPck term {ref!r} in primitive {primitive.identifier}"
            ) from exc
        if len(atoms) != 4:
            raise GICForgeContractError(
                f"invalid RPck dihedral term {ref!r} in primitive {primitive.identifier}"
            )
        terms.append((coefficient, atoms))
    return tuple(terms)


def _ring_index_for_atoms(
    atoms: tuple[int, ...],
    rings: tuple[tuple[int, tuple[int, ...]], ...],
) -> int | None:
    atom_set = set(atoms)
    best_index: int | None = None
    best_size: int | None = None
    for ring_index, ring_atoms in rings:
        ring_set = set(ring_atoms)
        if not atom_set.issubset(ring_set):
            continue
        ring_size = len(ring_atoms)
        if best_size is None or ring_size < best_size or (
            ring_size == best_size and ring_index < (best_index or ring_index)
        ):
            best_index = ring_index
            best_size = ring_size
    return best_index


def _ring_component_for_atoms(
    atoms: tuple[int, ...],
    rings: tuple[tuple[int, tuple[int, ...]], ...],
) -> tuple[int, ...] | None:
    atom_set = set(atoms)
    for component in _ring_components(rings):
        component_atoms: set[int] = set()
        for ring_index in component:
            component_atoms.update(dict(rings)[ring_index])
        if atom_set.issubset(component_atoms):
            return component
    return None


def _ring_components(
    rings: tuple[tuple[int, tuple[int, ...]], ...],
) -> tuple[tuple[int, ...], ...]:
    if not rings:
        return ()
    ring_ids = tuple(ring_index for ring_index, _atoms in rings)
    bond_to_rings = _ring_bond_to_rings(rings)
    neighbors: dict[int, set[int]] = {ring_index: set() for ring_index in ring_ids}
    for ring_indices in bond_to_rings.values():
        if len(ring_indices) < 2:
            continue
        for left, right in combinations(ring_indices, 2):
            neighbors[left].add(right)
            neighbors[right].add(left)
    components: list[tuple[int, ...]] = []
    seen: set[int] = set()
    for ring_index in ring_ids:
        if ring_index in seen:
            continue
        stack = [ring_index]
        seen.add(ring_index)
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in sorted(neighbors[current]):
                if neighbor in seen:
                    continue
                seen.add(neighbor)
                stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return tuple(components)


def _butterfly_torsion(
    atoms: tuple[int, int, int, int],
    rings: tuple[tuple[int, tuple[int, ...]], ...],
) -> bool:
    if not rings:
        return False
    i, j, k, l = atoms
    central_bond = tuple(sorted((j, k)))
    ring_by_index = {ring_index: set(ring_atoms) for ring_index, ring_atoms in rings}
    sharing_rings = _ring_bond_to_rings(rings).get(central_bond, ())
    if len(sharing_rings) < 2:
        return False
    for left_index, right_index in combinations(sharing_rings, 2):
        left_only = ring_by_index[left_index] - ring_by_index[right_index]
        right_only = ring_by_index[right_index] - ring_by_index[left_index]
        if (i in left_only and l in right_only) or (i in right_only and l in left_only):
            return True
    return False


def _ring_bond_to_rings(
    rings: tuple[tuple[int, tuple[int, ...]], ...],
) -> dict[tuple[int, int], tuple[int, ...]]:
    mapping: dict[tuple[int, int], list[int]] = {}
    for ring_index, ring_atoms in rings:
        if len(ring_atoms) < 2:
            continue
        for left, right in zip(ring_atoms, ring_atoms[1:] + ring_atoms[:1]):
            bond = tuple(sorted((left, right)))
            mapping.setdefault(bond, []).append(ring_index)
    return {bond: tuple(indices) for bond, indices in mapping.items()}


def _make_primitive(
    family: str,
    function: str,
    atoms: tuple[int, ...],
    counters: dict[str, int],
    *,
    mode: int = 0,
    ref_atoms: tuple[int, ...] = (),
    refs: tuple[str, ...] = (),
    frame_atoms: tuple[int, ...] = (),
    ref_frame_atoms: tuple[int, ...] = (),
) -> GICPrimitive:
    counters[family] += 1
    prefix = primitive_prefix(family)
    serial = sum(counters.values())
    return GICPrimitive(
        identifier=f"P{serial:03d}",
        name=f"{prefix}{counters[family]:04d}",
        family=family,
        function=function,
        atoms=tuple(int(atom) for atom in atoms),
        mode=int(mode),
        ref_atoms=tuple(int(atom) for atom in ref_atoms),
        refs=tuple(str(ref) for ref in refs),
        frame_atoms=tuple(int(atom) for atom in frame_atoms),
        ref_frame_atoms=tuple(int(atom) for atom in ref_frame_atoms),
    )


def _select_ranked_primitives(
    candidates: tuple[GICPrimitive, ...],
    coords: np.ndarray,
    *,
    target_rank: int,
    rank_tolerance: float,
) -> tuple[tuple[GICPrimitive, ...], int]:
    selected, rank, _ = _select_ranked_primitives_with_diagnostics(
        candidates,
        coords,
        target_rank=target_rank,
        rank_tolerance=rank_tolerance,
    )
    return selected, rank


def _select_ranked_primitives_with_diagnostics(
    candidates: tuple[GICPrimitive, ...],
    coords: np.ndarray,
    *,
    target_rank: int,
    rank_tolerance: float,
) -> tuple[tuple[GICPrimitive, ...], int, GICReductionDiagnostics]:
    skipped_singular: list[str] = []
    skipped_dependent: list[str] = []
    if target_rank == 0:
        return (
            (),
            0,
            GICReductionDiagnostics(
                rank_method=RANK_METHOD,
                reduction_policy=REDUCTION_POLICY,
            ),
        )
    selected: list[GICPrimitive] = []
    basis: list[np.ndarray] = []
    rank = 0

    special_candidates = [
        primitive for primitive in candidates if _is_special_primitive(primitive)
    ]
    ordinary_candidates = [
        primitive for primitive in candidates if not _is_special_primitive(primitive)
    ]
    for index, primitive in enumerate(special_candidates):
        if rank == target_rank:
            singular, dependent = _raise_if_remaining_special_independent(
                tuple(special_candidates[index:]),
                coords,
                basis,
                rank,
                rank_tolerance=rank_tolerance,
            )
            skipped_singular.extend(singular)
            skipped_dependent.extend(dependent)
            return (
                tuple(selected),
                rank,
                _make_reduction_diagnostics(
                    selected,
                    skipped_singular=skipped_singular,
                    skipped_dependent=skipped_dependent,
                ),
            )
        rank, status = _try_select_ranked_primitive(
            primitive,
            coords,
            selected,
            basis,
            rank,
            rank_tolerance=rank_tolerance,
        )
        _record_skip(primitive, status, skipped_singular, skipped_dependent)

    for primitive in ordinary_candidates:
        if rank == target_rank:
            break
        rank, status = _try_select_ranked_primitive(
            primitive,
            coords,
            selected,
            basis,
            rank,
            rank_tolerance=rank_tolerance,
        )
        _record_skip(primitive, status, skipped_singular, skipped_dependent)
    return (
        tuple(selected),
        rank,
        _make_reduction_diagnostics(
            selected,
            skipped_singular=skipped_singular,
            skipped_dependent=skipped_dependent,
        ),
    )


def _is_special_primitive(primitive: GICPrimitive) -> bool:
    return primitive.reduction_class == SPECIAL_REDUCTION_CLASS


def _protected_gic_count(primitives: tuple[GICPrimitive, ...]) -> int:
    return sum(1 for primitive in primitives if _is_special_primitive(primitive))


def _skipped_singular_count(definition: GICDefinition) -> int:
    if definition.reduction_diagnostics is None:
        return 0
    return len(definition.reduction_diagnostics.skipped_singular)


def _skipped_dependent_count(definition: GICDefinition) -> int:
    if definition.reduction_diagnostics is None:
        return 0
    return len(definition.reduction_diagnostics.skipped_dependent)


def _reduction_diagnostics_lines(definition: GICDefinition) -> list[str]:
    diagnostics = definition.reduction_diagnostics or GICReductionDiagnostics(
        rank_method=RANK_METHOD,
        reduction_policy=REDUCTION_POLICY,
        selected=tuple(primitive.identifier for primitive in definition.primitives),
    )
    return [
        f"RANK_METHOD {diagnostics.rank_method}",
        f"REDUCTION_POLICY {diagnostics.reduction_policy}",
        f"SELECTED {_csv_or_none(diagnostics.selected)}",
        f"SKIPPED_SINGULAR {_csv_or_none(diagnostics.skipped_singular)}",
        f"SKIPPED_DEPENDENT {_csv_or_none(diagnostics.skipped_dependent)}",
    ]


def _symmetry_diagnostics_lines(definition: GICDefinition) -> list[str]:
    diagnostics = definition.symmetry_diagnostics or _empty_symmetry_diagnostics(
        definition.point_group,
        requested=definition.symmetrize,
    )
    lines = [
        f"METHOD {diagnostics.method}",
        f"POLICY {diagnostics.policy}",
        f"STATUS {diagnostics.status}",
        f"POINT_GROUP {diagnostics.point_group}",
        f"SYMMETRY_GROUP {diagnostics.symmetry_group}",
        f"TOTAL_SYMMETRIC_IRREP {diagnostics.total_symmetric_irrep}",
        f"TOTAL_SYMMETRIC_GICS {_csv_or_none(diagnostics.total_symmetric_gics)}",
        f"GROUP_COUNT {len(diagnostics.groups)}",
        f"SYMMETRIZED_GIC_COUNT {_symmetrized_gic_count(diagnostics)}",
    ]
    for index, group in enumerate(diagnostics.groups, start=1):
        lines.append(
            f"GROUP {index} BLOCK={group.block} FAMILY={group.family} "
            f"SIGNATURE={group.signature} "
            f"SOURCES={_csv_or_none(group.source_gics)} "
            f"OUTPUTS={_csv_or_none(group.output_gics)}"
        )
    return lines


def _empty_symmetry_diagnostics(
    point_group: str,
    *,
    requested: bool,
) -> GICSymmetrizationDiagnostics:
    return GICSymmetrizationDiagnostics(
        method="NONE",
        policy=SYMMETRIZATION_POLICY,
        status="NOT_REQUESTED" if not requested else "NO_ELIGIBLE_GROUPS",
        point_group=point_group,
        symmetry_group=point_group,
        total_symmetric_irrep=total_symmetric_irrep(point_group),
        total_symmetric_gics=(),
    )


def _symmetrized_gic_count(diagnostics: GICSymmetrizationDiagnostics) -> int:
    return sum(len(group.output_gics) for group in diagnostics.groups)


def _make_reduction_diagnostics(
    selected: list[GICPrimitive],
    *,
    skipped_singular: list[str],
    skipped_dependent: list[str],
) -> GICReductionDiagnostics:
    return GICReductionDiagnostics(
        rank_method=RANK_METHOD,
        reduction_policy=REDUCTION_POLICY,
        selected=tuple(primitive.identifier for primitive in selected),
        skipped_singular=tuple(skipped_singular),
        skipped_dependent=tuple(skipped_dependent),
    )


def _record_skip(
    primitive: GICPrimitive,
    status: str,
    skipped_singular: list[str],
    skipped_dependent: list[str],
) -> None:
    if status == "singular":
        skipped_singular.append(primitive.identifier)
    elif status == "dependent":
        skipped_dependent.append(primitive.identifier)


def _apply_local_symmetrization(
    gics: tuple[FrozenGIC, ...],
    primitives: tuple[GICPrimitive, ...],
    *,
    atom_symbols: tuple[str, ...],
    point_group: str,
    requested: bool,
    symmetry_operations: tuple[GICPointGroupOperation, ...] = (),
) -> tuple[tuple[FrozenGIC, ...], GICSymmetrizationDiagnostics]:
    if not requested:
        return gics, _empty_symmetry_diagnostics(point_group, requested=False)

    projected = _apply_point_group_projector(
        gics,
        primitives,
        point_group=point_group,
        symmetry_operations=symmetry_operations,
    )
    if projected is not None:
        return projected

    source_groups = _local_symmetry_groups(gics, primitives, atom_symbols=atom_symbols)
    if not source_groups:
        prefixed_gics = _prefix_symmetrized_singletons(gics, point_group=point_group)
        return (
            prefixed_gics,
            GICSymmetrizationDiagnostics(
                method=LOCAL_SYMMETRIZATION_METHOD,
                policy=SYMMETRIZATION_POLICY,
                status="NO_ELIGIBLE_GROUPS",
                point_group=point_group,
                symmetry_group=point_group,
                total_symmetric_irrep=total_symmetric_irrep(point_group),
                total_symmetric_gics=tuple(
                    gic.name
                    for gic in prefixed_gics
                    if is_total_symmetric_irrep(point_group, gic.irrep)
                ),
            ),
        )

    groups_by_first = {
        group[0].identifier: (key, group)
        for key, group in source_groups.items()
    }
    grouped_ids = {
        gic.identifier
        for group in source_groups.values()
        for gic in group
    }
    name_counters: dict[tuple[str, str, str], int] = {}
    output: list[FrozenGIC] = []
    diagnostics: list[GICSymmetrizedGroup] = []

    for gic in gics:
        if gic.identifier in groups_by_first:
            key, group = groups_by_first[gic.identifier]
            new_gics = _symmetrized_group_gics(
                key,
                group,
                first_index=len(output) + 1,
                name_counters=name_counters,
                point_group=point_group,
            )
            output.extend(new_gics)
            diagnostics.append(
                GICSymmetrizedGroup(
                    block=key[0],
                    family=key[1],
                    signature=key[2],
                    source_gics=tuple(source.name for source in group),
                    output_gics=tuple(new_gic.name for new_gic in new_gics),
                )
            )
            continue
        if gic.identifier in grouped_ids:
            continue
        output.append(
            _renumber_frozen_gic(
                _prefix_symmetrized_gic(gic, point_group=point_group),
                len(output) + 1,
            )
        )

    output_tuple = tuple(output)
    return (
        output_tuple,
        GICSymmetrizationDiagnostics(
            method=LOCAL_SYMMETRIZATION_METHOD,
            policy=SYMMETRIZATION_POLICY,
            status="APPLIED",
            point_group=point_group,
            symmetry_group=point_group,
            total_symmetric_irrep=total_symmetric_irrep(point_group),
            total_symmetric_gics=tuple(
                gic.name
                for gic in output_tuple
                if is_total_symmetric_irrep(point_group, gic.irrep)
            ),
            groups=tuple(diagnostics),
        ),
    )


def _apply_point_group_projector(
    gics: tuple[FrozenGIC, ...],
    primitives: tuple[GICPrimitive, ...],
    *,
    point_group: str,
    symmetry_operations: tuple[GICPointGroupOperation, ...],
) -> tuple[tuple[FrozenGIC, ...], GICSymmetrizationDiagnostics] | None:
    operations = _valid_projector_operations(symmetry_operations)
    if len(operations) <= 1 or point_group.upper() in {"C1", "UNKNOWN"}:
        return None

    primitive_by_id = {primitive.identifier: primitive for primitive in primitives}
    blocks: dict[tuple[str, str], list[FrozenGIC]] = {}
    for gic in gics:
        key = (primitive_symmetry_block(gic.family), gic.family)
        blocks.setdefault(key, []).append(gic)

    output: list[FrozenGIC] = []
    diagnostics: list[GICSymmetrizedGroup] = []
    name_counters: dict[tuple[str, str], int] = {}
    for key, block_gics in blocks.items():
        projected = _project_gic_block(
            key,
            tuple(block_gics),
            primitive_by_id=primitive_by_id,
            operations=operations,
            point_group=point_group,
            first_index=len(output) + 1,
            name_counters=name_counters,
        )
        if projected is None:
            return None
        block_output, block_diagnostics = projected
        output.extend(block_output)
        diagnostics.append(block_diagnostics)

    output_tuple = tuple(output)
    return (
        output_tuple,
        GICSymmetrizationDiagnostics(
            method=POINT_GROUP_PROJECTOR_METHOD,
            policy=PROJECTOR_SYMMETRIZATION_POLICY,
            status="APPLIED",
            point_group=point_group,
            symmetry_group=point_group,
            total_symmetric_irrep=total_symmetric_irrep(point_group),
            total_symmetric_gics=tuple(
                gic.name
                for gic in output_tuple
                if is_total_symmetric_irrep(point_group, gic.irrep)
            ),
            groups=tuple(diagnostics),
        ),
    )


def _valid_projector_operations(
    operations: tuple[GICPointGroupOperation, ...],
) -> tuple[GICPointGroupOperation, ...]:
    if not operations:
        return ()
    natoms = len(operations[0].permutation)
    expected = tuple(range(1, natoms + 1))
    if natoms == 0:
        return ()
    validated: list[GICPointGroupOperation] = []
    seen: set[tuple[tuple[int, ...], tuple[float, ...]]] = set()
    for operation in operations:
        if len(operation.permutation) != natoms:
            return ()
        if tuple(sorted(operation.permutation)) != expected:
            return ()
        unique_key = (
            operation.permutation,
            tuple(
                round(float(value), 10)
                for row in operation.rotation
                for value in row
            ),
        )
        if unique_key in seen:
            continue
        seen.add(unique_key)
        validated.append(operation)
    identity_index = next(
        (
            idx
            for idx, operation in enumerate(validated)
            if operation.label == "E" and operation.permutation == expected
        ),
        None,
    )
    if identity_index is None:
        return ()
    if identity_index:
        identity = validated.pop(identity_index)
        validated.insert(0, identity)
    return tuple(validated)


def _project_gic_block(
    key: tuple[str, str],
    gics: tuple[FrozenGIC, ...],
    *,
    primitive_by_id: dict[str, GICPrimitive],
    operations: tuple[GICPointGroupOperation, ...],
    point_group: str,
    first_index: int,
    name_counters: dict[tuple[str, str], int],
) -> tuple[tuple[FrozenGIC, ...], GICSymmetrizedGroup] | None:
    block, family = key
    block_primitives = _block_primitives_for_gics(
        gics,
        primitive_by_id=primitive_by_id,
        key=key,
    )
    if block_primitives is None:
        return None

    primitive_index = {
        primitive.identifier: idx
        for idx, primitive in enumerate(block_primitives)
    }
    source_vectors = tuple(
        _gic_coefficient_vector(gic, primitive_index=primitive_index)
        for gic in gics
    )
    if any(vector is None for vector in source_vectors):
        return None

    primitive_key_index = _primitive_projector_key_index(block_primitives)
    if primitive_key_index is None:
        return None
    transforms = tuple(
        _operation_primitive_transform(
            block_primitives,
            operation=operation,
            primitive_key_index=primitive_key_index,
        )
        for operation in operations
    )
    if any(transform is None for transform in transforms):
        return None

    projected_vectors: list[tuple[str, np.ndarray]] = []
    basis: list[np.ndarray] = []
    operation_labels = tuple(operation.label for operation in operations)
    operation_matrices = tuple(operation.rotation for operation in operations)
    for irrep, characters in irrep_characters_for_operations(
        operation_labels,
        point_group,
        operation_matrices=operation_matrices,
    ):
        if len(characters) != len(operations):
            return None
        if all(abs(character) <= 1.0e-14 for character in characters):
            continue
        for source_vector in source_vectors:
            assert source_vector is not None
            projected = _project_vector_for_irrep(
                source_vector,
                characters=characters,
                transforms=transforms,
            )
            normalized = _normalized_coefficient_vector_or_none(projected)
            if normalized is None:
                continue
            independent = _orthonormal_coefficient_residual_or_none(basis, normalized)
            if independent is None:
                continue
            basis.append(independent)
            projected_vectors.append((irrep, independent))
            if len(projected_vectors) == len(block_primitives):
                break
        if len(projected_vectors) == len(block_primitives):
            break
    if len(projected_vectors) != len(block_primitives):
        return None

    output: list[FrozenGIC] = []
    for offset, (irrep, vector) in enumerate(projected_vectors):
        coefficients = _coefficients_from_vector(block_primitives, vector)
        if not coefficients:
            return None
        output.append(
            FrozenGIC(
                identifier=f"GIC{first_index + offset:03d}",
                name=_next_projected_name(family, irrep, name_counters),
                family=family,
                irrep=irrep,
                primitive_id=coefficients[0][0],
                gaussian_expression="LINEAR_COMBINATION",
                coefficients=coefficients,
            )
        )

    return (
        tuple(output),
        GICSymmetrizedGroup(
            block=block,
            family=family,
            signature="OPS=" + ",".join(operation_labels),
            source_gics=tuple(gic.name for gic in gics),
            output_gics=tuple(gic.name for gic in output),
        ),
    )


def _block_primitives_for_gics(
    gics: tuple[FrozenGIC, ...],
    *,
    primitive_by_id: dict[str, GICPrimitive],
    key: tuple[str, str],
) -> tuple[GICPrimitive, ...] | None:
    block, family = key
    ordered: list[GICPrimitive] = []
    seen: set[str] = set()
    for gic in gics:
        coefficients = gic.coefficients or ((gic.primitive_id, 1.0),)
        for primitive_id, _coefficient in coefficients:
            primitive = primitive_by_id.get(primitive_id)
            if primitive is None:
                return None
            primitive_key = (primitive_symmetry_block(primitive.family), primitive.family)
            if primitive_key != (block, family):
                return None
            if primitive.identifier in seen:
                continue
            seen.add(primitive.identifier)
            ordered.append(primitive)
    return tuple(ordered)


def _gic_coefficient_vector(
    gic: FrozenGIC,
    *,
    primitive_index: dict[str, int],
) -> np.ndarray | None:
    vector = np.zeros(len(primitive_index), dtype=float)
    coefficients = gic.coefficients or ((gic.primitive_id, 1.0),)
    for primitive_id, coefficient in coefficients:
        idx = primitive_index.get(primitive_id)
        if idx is None:
            return None
        vector[idx] += float(coefficient)
    return vector


def _primitive_projector_key_index(
    primitives: tuple[GICPrimitive, ...],
) -> dict[tuple[object, ...], int] | None:
    out: dict[tuple[object, ...], int] = {}
    for idx, primitive in enumerate(primitives):
        key = _primitive_projector_key(primitive)
        if key is None or key in out:
            return None
        out[key] = idx
    return out


def _operation_primitive_transform(
    primitives: tuple[GICPrimitive, ...],
    *,
    operation: GICPointGroupOperation,
    primitive_key_index: dict[tuple[object, ...], int],
) -> np.ndarray | None:
    matrix = np.zeros((len(primitives), len(primitives)), dtype=float)
    for source_index, primitive in enumerate(primitives):
        terms = _mapped_primitive_projector_terms(primitive, operation)
        if terms is None:
            return None
        for target_key, coefficient in terms:
            if abs(float(coefficient)) <= 1.0e-12:
                continue
            target_index = primitive_key_index.get(target_key)
            if target_index is None:
                return None
            matrix[target_index, source_index] += float(coefficient)
    return matrix


def _project_vector_for_irrep(
    vector: np.ndarray,
    *,
    characters: tuple[float, ...],
    transforms: tuple[np.ndarray | None, ...],
) -> np.ndarray:
    projected = np.zeros_like(vector, dtype=float)
    for character, transform in zip(characters, transforms):
        assert transform is not None
        projected += float(character) * (transform @ vector)
    return projected / float(len(transforms))


def _normalized_coefficient_vector_or_none(vector: np.ndarray) -> np.ndarray | None:
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm <= 1.0e-10:
        return None
    return vector / norm


def _orthonormal_coefficient_residual_or_none(
    basis: list[np.ndarray],
    normalized: np.ndarray,
) -> np.ndarray | None:
    residual = np.array(normalized, dtype=float, copy=True)
    for vector in basis:
        residual -= float(np.dot(residual, vector)) * vector
    norm = float(np.linalg.norm(residual))
    if not np.isfinite(norm) or norm <= 1.0e-10:
        return None
    return residual / norm


def _coefficients_from_vector(
    primitives: tuple[GICPrimitive, ...],
    vector: np.ndarray,
) -> tuple[tuple[str, float], ...]:
    return tuple(
        (primitive.identifier, float(value))
        for primitive, value in zip(primitives, vector)
        if abs(float(value)) > 1.0e-12
    )


def _next_projected_name(
    family: str,
    irrep: str,
    counters: dict[tuple[str, str], int],
) -> str:
    key = (family, irrep)
    counters[key] = counters.get(key, 0) + 1
    return f"{irrep_name_prefix(irrep)}{primitive_prefix(family)}{counters[key]:03d}"


def _primitive_projector_key(primitive: GICPrimitive) -> tuple[object, ...] | None:
    if primitive.family == "STRETCH" and len(primitive.atoms) == 2:
        return ("STRETCH", tuple(sorted(primitive.atoms)))
    if primitive.family == "BEND" and len(primitive.atoms) == 3:
        return (
            "BEND",
            primitive.atoms[1],
            tuple(sorted((primitive.atoms[0], primitive.atoms[2]))),
        )
    if primitive.family == "LINEAR_BEND" and len(primitive.atoms) == 3:
        return (
            "LINEAR_BEND",
            primitive.atoms[1],
            tuple(sorted((primitive.atoms[0], primitive.atoms[2]))),
            primitive.mode,
        )
    if primitive.family == "TORSION" and len(primitive.atoms) == 4:
        canonical, _sign = _canonical_torsion_key_and_sign(primitive.atoms)
        return ("TORSION", canonical)
    if primitive.family == "RING_PUCKER_COMPONENT" and primitive.function == "RPCK":
        signature = _ring_pucker_projector_signature(primitive)
        if signature is None:
            return None
        key, _sign = signature
        return ("RING_PUCKER_COMPONENT", key)
    if primitive.family == "OUT_OF_PLANE" and len(primitive.atoms) == 4:
        return (
            "OUT_OF_PLANE",
            primitive.atoms[0],
            tuple(sorted(primitive.atoms[1:])),
        )
    if primitive.family == "FRAG_DISTANCE":
        pair = tuple(
            sorted((_atom_set_key(primitive.atoms), _atom_set_key(primitive.ref_atoms)))
        )
        return ("FRAG_DISTANCE", pair)
    if primitive.family == "FRAG_CENTER_ATOM_DISTANCE":
        return (
            "FRAG_CENTER_ATOM_DISTANCE",
            _atom_set_key(primitive.atoms),
            _atom_set_key(primitive.ref_atoms),
        )
    if primitive.family == "CENTER_ATOM_DISTANCE":
        return (
            "CENTER_ATOM_DISTANCE",
            _atom_set_key(primitive.atoms),
            _atom_set_key(primitive.ref_atoms),
        )
    if primitive.family == "FRAG_TRANSLATION":
        return (
            "FRAG_TRANSLATION",
            primitive.mode,
            _atom_set_key(primitive.atoms),
            _atom_set_key(primitive.ref_atoms),
        )
    if primitive.family == "FRAG_ORIENTATION":
        return (
            "FRAG_ORIENTATION",
            primitive.mode,
            _atom_set_key(primitive.atoms),
            _atom_set_key(primitive.ref_atoms),
            _atom_set_key(primitive.frame_atoms),
            _atom_set_key(primitive.ref_frame_atoms),
        )
    return None


def _mapped_primitive_projector_terms(
    primitive: GICPrimitive,
    operation: GICPointGroupOperation,
) -> tuple[tuple[tuple[object, ...], float], ...] | None:
    mapped_atoms = tuple(_mapped_atom(operation, atom) for atom in primitive.atoms)
    mapped_refs = tuple(_mapped_atom(operation, atom) for atom in primitive.ref_atoms)
    mapped_frame = tuple(_mapped_atom(operation, atom) for atom in primitive.frame_atoms)
    mapped_ref_frame = tuple(
        _mapped_atom(operation, atom) for atom in primitive.ref_frame_atoms
    )
    if any(atom < 1 for atom in mapped_atoms + mapped_refs + mapped_frame + mapped_ref_frame):
        return None

    if primitive.family == "STRETCH" and len(mapped_atoms) == 2:
        return ((("STRETCH", tuple(sorted(mapped_atoms))), 1.0),)
    if primitive.family == "BEND" and len(mapped_atoms) == 3:
        return (
            (
                (
                    "BEND",
                    mapped_atoms[1],
                    tuple(sorted((mapped_atoms[0], mapped_atoms[2]))),
                ),
                1.0,
            ),
        )
    if primitive.family == "LINEAR_BEND":
        if not _is_identity_operation(operation):
            return None
        key = _primitive_projector_key(primitive)
        return ((key, 1.0),) if key is not None else None
    if primitive.family == "TORSION" and len(mapped_atoms) == 4:
        canonical, sign = _canonical_torsion_key_and_sign(mapped_atoms)
        return ((("TORSION", canonical), sign),)
    if primitive.family == "RING_PUCKER_COMPONENT" and primitive.function == "RPCK":
        mapped_terms = []
        for coefficient, atoms in _ring_pucker_terms_from_refs(primitive):
            mapped_term_atoms = tuple(_mapped_atom(operation, atom) for atom in atoms)
            if any(atom < 1 for atom in mapped_term_atoms):
                return None
            mapped_terms.append((coefficient, mapped_term_atoms))
        signature = _ring_pucker_projector_signature_from_terms(tuple(mapped_terms))
        if signature is None:
            return None
        key, sign = signature
        return ((("RING_PUCKER_COMPONENT", key), sign),)
    if primitive.family == "OUT_OF_PLANE" and len(mapped_atoms) == 4:
        center = mapped_atoms[0]
        substituents = mapped_atoms[1:]
        sorted_substituents = tuple(sorted(substituents))
        return (
            (
                (
                    "OUT_OF_PLANE",
                    center,
                    sorted_substituents,
                ),
                _permutation_parity_sign(substituents, sorted_substituents),
            ),
        )
    if primitive.family == "FRAG_DISTANCE":
        pair = tuple(sorted((_atom_set_key(mapped_atoms), _atom_set_key(mapped_refs))))
        return ((("FRAG_DISTANCE", pair), 1.0),)
    if primitive.family == "FRAG_CENTER_ATOM_DISTANCE":
        return (
            (
                (
                    "FRAG_CENTER_ATOM_DISTANCE",
                    _atom_set_key(mapped_atoms),
                    _atom_set_key(mapped_refs),
                ),
                1.0,
            ),
        )
    if primitive.family == "CENTER_ATOM_DISTANCE":
        return (
            (
                (
                    "CENTER_ATOM_DISTANCE",
                    _atom_set_key(mapped_atoms),
                    _atom_set_key(mapped_refs),
                ),
                1.0,
            ),
        )
    if primitive.family == "FRAG_TRANSLATION":
        return tuple(
            (
                (
                    "FRAG_TRANSLATION",
                    target_mode,
                    _atom_set_key(mapped_atoms),
                    _atom_set_key(mapped_refs),
                ),
                coefficient,
            )
            for target_mode, coefficient in _vector_component_terms(
                operation,
                source_mode=primitive.mode,
                axial=False,
            )
        )
    if primitive.family == "FRAG_ORIENTATION":
        return tuple(
            (
                (
                    "FRAG_ORIENTATION",
                    target_mode,
                    _atom_set_key(mapped_atoms),
                    _atom_set_key(mapped_refs),
                    _atom_set_key(mapped_frame),
                    _atom_set_key(mapped_ref_frame),
                ),
                coefficient,
            )
            for target_mode, coefficient in _vector_component_terms(
                operation,
                source_mode=primitive.mode,
                axial=True,
            )
        )
    return None


def _mapped_atom(operation: GICPointGroupOperation, atom: int) -> int:
    if 1 <= atom <= len(operation.permutation):
        return operation.permutation[atom - 1]
    return -1


def _is_identity_operation(operation: GICPointGroupOperation) -> bool:
    return operation.permutation == tuple(range(1, len(operation.permutation) + 1))


def _vector_component_terms(
    operation: GICPointGroupOperation,
    *,
    source_mode: int,
    axial: bool,
) -> tuple[tuple[int, float], ...]:
    if source_mode not in {0, 1, 2}:
        return ()
    rotation = np.asarray(operation.rotation, dtype=float)
    if rotation.shape != (3, 3):
        return ()
    if axial:
        rotation = float(np.linalg.det(rotation)) * rotation
    return tuple(
        (target_mode, float(rotation[source_mode, target_mode]))
        for target_mode in range(3)
        if abs(float(rotation[source_mode, target_mode])) > 1.0e-10
    )


def _atom_set_key(atoms: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(int(atom) for atom in atoms))


def _canonical_torsion_key_and_sign(atoms: tuple[int, ...]) -> tuple[tuple[int, ...], float]:
    forward = tuple(atoms)
    backward = tuple(reversed(atoms))
    if backward < forward:
        return backward, -1.0
    return forward, 1.0


def _ring_pucker_projector_signature(
    primitive: GICPrimitive,
) -> tuple[tuple[tuple[tuple[int, ...], float], ...], float] | None:
    return _ring_pucker_projector_signature_from_terms(
        _ring_pucker_terms_from_refs(primitive)
    )


def _ring_pucker_projector_signature_from_terms(
    terms: tuple[tuple[float, tuple[int, ...]], ...],
) -> tuple[tuple[tuple[tuple[int, ...], float], ...], float] | None:
    by_torsion: dict[tuple[int, ...], float] = {}
    for coefficient, atoms in terms:
        if len(atoms) != 4:
            return None
        canonical, sign = _canonical_torsion_key_and_sign(atoms)
        by_torsion[canonical] = by_torsion.get(canonical, 0.0) + float(coefficient) * sign
    compact = {
        atoms: coefficient
        for atoms, coefficient in by_torsion.items()
        if abs(float(coefficient)) > 1.0e-12
    }
    if not compact:
        return None
    dominant_atoms, dominant_coefficient = max(
        compact.items(),
        key=lambda item: (abs(float(item[1])), tuple(-atom for atom in item[0])),
    )
    del dominant_atoms
    overall_sign = -1.0 if dominant_coefficient < 0.0 else 1.0
    key = tuple(
        (atoms, round(float(coefficient) * overall_sign, 12))
        for atoms, coefficient in sorted(compact.items())
    )
    return key, overall_sign


def _permutation_parity_sign(
    order: tuple[int, ...],
    target_order: tuple[int, ...],
) -> float:
    positions = {atom: idx for idx, atom in enumerate(target_order)}
    try:
        indexes = [positions[atom] for atom in order]
    except KeyError:
        return 1.0
    inversions = 0
    for left in range(len(indexes)):
        for right in range(left + 1, len(indexes)):
            if indexes[left] > indexes[right]:
                inversions += 1
    return -1.0 if inversions % 2 else 1.0


def _local_symmetry_groups(
    gics: tuple[FrozenGIC, ...],
    primitives: tuple[GICPrimitive, ...],
    *,
    atom_symbols: tuple[str, ...],
) -> dict[tuple[str, str, str], tuple[FrozenGIC, ...]]:
    primitive_by_id = {primitive.identifier: primitive for primitive in primitives}
    grouped: dict[tuple[str, str, str], list[FrozenGIC]] = {}
    for gic in gics:
        primitive = _single_source_primitive(gic, primitive_by_id)
        if primitive is None:
            continue
        signature = _local_symmetry_signature(primitive, atom_symbols=atom_symbols)
        if not signature:
            continue
        key = (primitive_symmetry_block(primitive.family), primitive.family, signature)
        grouped.setdefault(key, []).append(gic)
    return {
        key: tuple(group)
        for key, group in grouped.items()
        if len(group) > 1
    }


def _single_source_primitive(
    gic: FrozenGIC,
    primitive_by_id: dict[str, GICPrimitive],
) -> GICPrimitive | None:
    coefficients = gic.coefficients or ((gic.primitive_id, 1.0),)
    if len(coefficients) != 1:
        return None
    primitive_id, coefficient = coefficients[0]
    if abs(float(coefficient) - 1.0) > 1.0e-12:
        return None
    return primitive_by_id.get(primitive_id)


def _local_symmetry_signature(
    primitive: GICPrimitive,
    *,
    atom_symbols: tuple[str, ...],
) -> str | None:
    if primitive.family == "STRETCH" and len(primitive.atoms) == 2:
        return "R:" + "-".join(_sorted_atom_symbols(primitive.atoms, atom_symbols))
    if primitive.family == "BEND" and len(primitive.atoms) == 3:
        end_symbols = _sorted_atom_symbols(
            (primitive.atoms[0], primitive.atoms[2]),
            atom_symbols,
        )
        return f"A:{_atom_symbol(primitive.atoms[1], atom_symbols)}:{'-'.join(end_symbols)}"
    if primitive.family == "LINEAR_BEND" and len(primitive.atoms) == 3:
        end_symbols = _sorted_atom_symbols(
            (primitive.atoms[0], primitive.atoms[2]),
            atom_symbols,
        )
        return (
            f"L:{primitive.mode}:{_atom_symbol(primitive.atoms[1], atom_symbols)}:"
            f"{'-'.join(end_symbols)}"
        )
    if primitive.family == "TORSION" and len(primitive.atoms) == 4:
        return "D:" + "-".join(_atom_symbol(atom, atom_symbols) for atom in primitive.atoms)
    if primitive.family == "OUT_OF_PLANE" and len(primitive.atoms) == 4:
        substituents = _sorted_atom_symbols(primitive.atoms[1:], atom_symbols)
        return (
            f"U:{_atom_symbol(primitive.atoms[0], atom_symbols)}:"
            f"{'-'.join(substituents)}"
        )
    if primitive.family == "FRAG_DISTANCE":
        left = _atom_multiset_signature(primitive.atoms, atom_symbols)
        right = _atom_multiset_signature(primitive.ref_atoms, atom_symbols)
        pair = tuple(sorted((left, right)))
        return f"FC_DIST:{pair[0]}:{pair[1]}"
    if primitive.family == "FRAG_CENTER_ATOM_DISTANCE":
        return (
            "FCA_DIST:"
            f"{_atom_multiset_signature(primitive.atoms, atom_symbols)}:"
            f"{_atom_multiset_signature(primitive.ref_atoms, atom_symbols)}"
        )
    if primitive.family == "FRAG_TRANSLATION":
        return (
            f"FTRANS:{primitive.mode}:"
            f"{_atom_multiset_signature(primitive.atoms, atom_symbols)}:"
            f"{_atom_multiset_signature(primitive.ref_atoms, atom_symbols)}"
        )
    if primitive.family == "FRAG_ORIENTATION":
        return (
            f"FROT:{primitive.mode}:"
            f"{_atom_multiset_signature(primitive.atoms, atom_symbols)}:"
            f"{_atom_multiset_signature(primitive.ref_atoms, atom_symbols)}:"
            f"{_atom_multiset_signature(primitive.frame_atoms, atom_symbols)}:"
            f"{_atom_multiset_signature(primitive.ref_frame_atoms, atom_symbols)}"
        )
    if primitive.family == "CENTER_ATOM_DISTANCE":
        return (
            "CENTER_ATOM_DIST:"
            f"{_atom_multiset_signature(primitive.atoms, atom_symbols)}:"
            f"{_atom_multiset_signature(primitive.ref_atoms, atom_symbols)}"
        )
    return None


def _symmetrized_group_gics(
    key: tuple[str, str, str],
    group: tuple[FrozenGIC, ...],
    *,
    first_index: int,
    name_counters: dict[tuple[str, str, str], int],
    point_group: str,
) -> tuple[FrozenGIC, ...]:
    _block, family, _signature = key
    size = len(group)
    output: list[FrozenGIC] = []
    symmetric_weight = 1.0 / np.sqrt(float(size))
    symmetric_irrep = _local_symmetry_irrep(point_group=point_group, kind="S", index=0)
    output.append(
        FrozenGIC(
            identifier=f"GIC{first_index:03d}",
            name=_next_symmetrized_name(family, "S", name_counters, irrep=symmetric_irrep),
            family=family,
            irrep=symmetric_irrep,
            primitive_id=group[0].primitive_id,
            gaussian_expression="LINEAR_COMBINATION",
            coefficients=_combine_gic_coefficients(
                group,
                tuple(symmetric_weight for _idx in range(size)),
            ),
        )
    )
    for group_index in range(1, size):
        weights = [0.0 for _idx in group]
        weights[group_index - 1] = 1.0 / np.sqrt(2.0)
        weights[group_index] = -1.0 / np.sqrt(2.0)
        difference_irrep = _local_symmetry_irrep(
            point_group=point_group,
            kind="D",
            index=group_index - 1,
        )
        output.append(
            FrozenGIC(
                identifier=f"GIC{first_index + group_index:03d}",
                name=_next_symmetrized_name(
                    family,
                    "D",
                    name_counters,
                    irrep=difference_irrep,
                ),
                family=family,
                irrep=difference_irrep,
                primitive_id=group[group_index - 1].primitive_id,
                gaussian_expression="LINEAR_COMBINATION",
                coefficients=_combine_gic_coefficients(group, tuple(weights)),
            )
        )
    return tuple(output)


def _combine_gic_coefficients(
    gics: tuple[FrozenGIC, ...],
    weights: tuple[float, ...],
) -> tuple[tuple[str, float], ...]:
    totals: dict[str, float] = {}
    order: list[str] = []
    for gic, weight in zip(gics, weights):
        if abs(weight) <= 1.0e-14:
            continue
        coefficients = gic.coefficients or ((gic.primitive_id, 1.0),)
        for primitive_id, coefficient in coefficients:
            if primitive_id not in totals:
                order.append(primitive_id)
                totals[primitive_id] = 0.0
            totals[primitive_id] += float(weight) * float(coefficient)
    return tuple(
        (primitive_id, totals[primitive_id])
        for primitive_id in order
        if abs(totals[primitive_id]) > 1.0e-14
    )


def _renumber_frozen_gic(gic: FrozenGIC, index: int) -> FrozenGIC:
    return FrozenGIC(
        identifier=f"GIC{index:03d}",
        name=gic.name,
        family=gic.family,
        irrep=gic.irrep,
        primitive_id=gic.primitive_id,
        gaussian_expression=gic.gaussian_expression,
        coefficients=gic.coefficients,
    )


def _prefix_symmetrized_singletons(
    gics: tuple[FrozenGIC, ...],
    *,
    point_group: str,
) -> tuple[FrozenGIC, ...]:
    return tuple(_prefix_symmetrized_gic(gic, point_group=point_group) for gic in gics)


def _prefix_symmetrized_gic(gic: FrozenGIC, *, point_group: str) -> FrozenGIC:
    irrep = (
        gic.irrep
        if gic.irrep and gic.irrep != "UNASSIGNED"
        else total_symmetric_irrep(point_group)
    )
    prefix = irrep_name_prefix(irrep)
    name = gic.name if gic.name.startswith(prefix) else f"{prefix}{gic.name}"
    return FrozenGIC(
        identifier=gic.identifier,
        name=name,
        family=gic.family,
        irrep=irrep,
        primitive_id=gic.primitive_id,
        gaussian_expression=gic.gaussian_expression,
        coefficients=gic.coefficients,
    )


def _next_symmetrized_name(
    family: str,
    kind: str,
    counters: dict[tuple[str, str, str], int],
    *,
    irrep: str,
) -> str:
    key = (family, kind, irrep)
    counters[key] = counters.get(key, 0) + 1
    return f"{irrep_name_prefix(irrep)}{primitive_prefix(family)}{kind}{counters[key]:03d}"


def _local_symmetry_irrep(
    *,
    point_group: str,
    kind: str,
    index: int,
) -> str:
    if kind == "S":
        return total_symmetric_irrep(point_group)
    non_total = non_total_irrep_sequence(point_group)
    if not non_total:
        return total_symmetric_irrep(point_group)
    return non_total[index % len(non_total)]


def _sorted_atom_symbols(
    atoms: tuple[int, ...],
    atom_symbols: tuple[str, ...],
) -> tuple[str, ...]:
    return tuple(sorted(_atom_symbol(atom, atom_symbols) for atom in atoms))


def _atom_multiset_signature(
    atoms: tuple[int, ...],
    atom_symbols: tuple[str, ...],
) -> str:
    if not atoms:
        return "NONE"
    return ".".join(_sorted_atom_symbols(atoms, atom_symbols))


def _atom_symbol(atom: int, atom_symbols: tuple[str, ...]) -> str:
    if 1 <= atom <= len(atom_symbols):
        return atom_symbols[atom - 1].upper()
    return f"A{atom}"


def _try_select_ranked_primitive(
    primitive: GICPrimitive,
    coords: np.ndarray,
    selected: list[GICPrimitive],
    basis: list[np.ndarray],
    rank: int,
    *,
    rank_tolerance: float,
) -> tuple[int, str]:
    normalized = _normalized_b_row_or_none(
        primitive,
        coords,
        rank_tolerance=rank_tolerance,
    )
    if normalized is None:
        return rank, "singular"
    orthonormal = _orthonormal_residual_or_none(
        basis,
        normalized,
        rank_tolerance=rank_tolerance,
    )
    if orthonormal is None:
        return rank, "dependent"
    selected.append(primitive)
    basis.append(orthonormal)
    return rank + 1, "selected"


def _raise_if_remaining_special_independent(
    primitives: tuple[GICPrimitive, ...],
    coords: np.ndarray,
    basis: list[np.ndarray],
    rank: int,
    *,
    rank_tolerance: float,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    skipped_singular: list[str] = []
    skipped_dependent: list[str] = []
    for primitive in primitives:
        normalized = _normalized_b_row_or_none(
            primitive,
            coords,
            rank_tolerance=rank_tolerance,
        )
        if normalized is None:
            skipped_singular.append(primitive.identifier)
            continue
        if _orthonormal_residual_or_none(
            basis,
            normalized,
            rank_tolerance=rank_tolerance,
        ) is not None:
            raise GICForgeContractError(
                "protected special primitive set exceeds the vibrational rank: "
                f"{primitive.identifier} {primitive.name} would add an independent "
                "row after the target rank was reached"
            )
        skipped_dependent.append(primitive.identifier)
    return tuple(skipped_singular), tuple(skipped_dependent)


def _normalized_b_row_or_none(
    primitive: GICPrimitive,
    coords: np.ndarray,
    *,
    rank_tolerance: float,
) -> np.ndarray | None:
    try:
        row = _analytic_b_row(primitive, coords)
    except FloatingPointError:
        return None
    norm = float(np.linalg.norm(row))
    if not np.isfinite(norm) or norm <= rank_tolerance:
        return None
    return row / norm


def _orthonormal_residual_or_none(
    basis: list[np.ndarray],
    normalized: np.ndarray,
    *,
    rank_tolerance: float,
) -> np.ndarray | None:
    residual = np.array(normalized, dtype=float, copy=True)
    for vector in basis:
        residual -= float(np.dot(residual, vector)) * vector
    norm = float(np.linalg.norm(residual))
    if not np.isfinite(norm) or norm <= rank_tolerance:
        return None
    return residual / norm


def _analytic_b_row(primitive: GICPrimitive, coords: np.ndarray) -> np.ndarray:
    coords = np.asarray(coords, dtype=float)
    if primitive.function == "R":
        return _distance_b_row(coords, primitive.atoms)
    if primitive.function == "A":
        return _angle_b_row(coords, primitive.atoms)
    if primitive.function == "FC_DIST":
        return _fragment_center_distance_b_row(coords, primitive.atoms, primitive.ref_atoms)
    if primitive.function == "FCA_DIST":
        return _fragment_center_atom_distance_b_row(coords, primitive.atoms, primitive.ref_atoms)
    if primitive.function == "CENTER_ATOM_DIST":
        return _fragment_center_atom_distance_b_row(coords, primitive.atoms, primitive.ref_atoms)
    if primitive.function == "FTRANS":
        return _fragment_translation_b_row(
            coords,
            primitive.atoms,
            primitive.ref_atoms,
            mode=primitive.mode,
        )
    if primitive.function == "RPCK":
        return _ring_pucker_component_b_row(primitive, coords)
    return _dual_b_row(primitive, coords)


def _distance_b_row(coords: np.ndarray, atoms: tuple[int, ...]) -> np.ndarray:
    i, j = (atom - 1 for atom in atoms)
    delta = coords[i] - coords[j]
    distance = float(np.linalg.norm(delta))
    if distance <= RANK_TOLERANCE:
        raise FloatingPointError("zero-length distance coordinate")
    row = np.zeros(coords.size, dtype=float)
    unit = delta / distance
    row[3 * i : 3 * i + 3] = unit
    row[3 * j : 3 * j + 3] = -unit
    return row


def _angle_b_row(coords: np.ndarray, atoms: tuple[int, ...]) -> np.ndarray:
    i, j, k = (atom - 1 for atom in atoms)
    rji = coords[i] - coords[j]
    rjk = coords[k] - coords[j]
    dji = float(np.linalg.norm(rji))
    djk = float(np.linalg.norm(rjk))
    if dji <= RANK_TOLERANCE or djk <= RANK_TOLERANCE:
        raise FloatingPointError("zero-length angle arm")
    eji = rji / dji
    ejk = rjk / djk
    cosine = float(np.clip(np.dot(eji, ejk), -1.0, 1.0))
    sine = float(np.sqrt(max(1.0 - cosine * cosine, 0.0)))
    if sine <= RANK_TOLERANCE:
        raise FloatingPointError("linear angle has no ordinary bend derivative")
    gi = (cosine * eji - ejk) / (dji * sine)
    gk = (cosine * ejk - eji) / (djk * sine)
    row = np.zeros(coords.size, dtype=float)
    row[3 * i : 3 * i + 3] = gi
    row[3 * k : 3 * k + 3] = gk
    row[3 * j : 3 * j + 3] = -(gi + gk)
    return row


def _fragment_center_distance_b_row(
    coords: np.ndarray,
    atoms: tuple[int, ...],
    ref_atoms: tuple[int, ...],
) -> np.ndarray:
    center = _fragment_center(coords, atoms)
    ref_center = _fragment_center(coords, ref_atoms)
    delta = center - ref_center
    distance = float(np.linalg.norm(delta))
    if distance <= RANK_TOLERANCE:
        raise FloatingPointError("coincident fragment centers")
    unit = delta / distance
    row = np.zeros(coords.size, dtype=float)
    _accumulate_center_gradient(row, atoms, unit / len(atoms))
    _accumulate_center_gradient(row, ref_atoms, -unit / len(ref_atoms))
    return row


def _fragment_center_atom_distance_b_row(
    coords: np.ndarray,
    atoms: tuple[int, ...],
    ref_atoms: tuple[int, ...],
) -> np.ndarray:
    if len(ref_atoms) != 1:
        raise FloatingPointError("center-atom distance needs exactly one reference atom")
    atom = ref_atoms[0]
    delta = _fragment_center(coords, atoms) - coords[atom - 1]
    distance = float(np.linalg.norm(delta))
    if distance <= RANK_TOLERANCE:
        raise FloatingPointError("fragment center and atom are coincident")
    unit = delta / distance
    row = np.zeros(coords.size, dtype=float)
    _accumulate_center_gradient(row, atoms, unit / len(atoms))
    row[3 * (atom - 1) : 3 * atom] -= unit
    return row


def _fragment_translation_b_row(
    coords: np.ndarray,
    atoms: tuple[int, ...],
    ref_atoms: tuple[int, ...],
    *,
    mode: int,
) -> np.ndarray:
    row = np.zeros(coords.size, dtype=float)
    axis = np.zeros(3, dtype=float)
    axis[mode] = 1.0
    _accumulate_center_gradient(row, atoms, axis / len(atoms))
    _accumulate_center_gradient(row, ref_atoms, -axis / len(ref_atoms))
    return row


def _ring_pucker_component_b_row(
    primitive: GICPrimitive,
    coords: np.ndarray,
) -> np.ndarray:
    row = np.zeros(coords.size, dtype=float)
    for coefficient, atoms in _ring_pucker_terms_from_refs(primitive):
        if coefficient == 0.0:
            continue
        term = GICPrimitive(
            identifier=f"{primitive.identifier}_D",
            name="RPckD",
            family="TORSION",
            function="D",
            atoms=atoms,
        )
        row += coefficient * _dual_b_row(term, coords)
    return row


def _accumulate_center_gradient(
    row: np.ndarray,
    atoms: tuple[int, ...],
    gradient: np.ndarray,
) -> None:
    for atom in atoms:
        start = 3 * (atom - 1)
        row[start : start + 3] += gradient


def _dual_b_row(primitive: GICPrimitive, coords: np.ndarray) -> np.ndarray:
    dcoords = _dual_coordinates(coords)
    value = _dual_primitive_value(primitive, dcoords, coords)
    if not np.isfinite(value.val):
        raise FloatingPointError("non-finite analytic derivative value")
    return np.asarray(value.der, dtype=float)


def _finite_difference_b_row(
    primitive: GICPrimitive,
    coords: np.ndarray,
    *,
    step_angstrom: float = DIAGNOSTIC_FINITE_DIFFERENCE_STEP,
) -> np.ndarray:
    flat = np.asarray(coords, dtype=float).reshape(-1)
    base = _primitive_value(primitive, coords)
    row = np.zeros_like(flat)
    for idx in range(flat.size):
        plus = flat.copy()
        minus = flat.copy()
        plus[idx] += step_angstrom
        minus[idx] -= step_angstrom
        try:
            value_plus = _primitive_value(primitive, plus.reshape(coords.shape))
            value_minus = _primitive_value(primitive, minus.reshape(coords.shape))
        except FloatingPointError:
            row[idx] = np.nan
            continue
        if primitive.function == "RPCK":
            delta = _ring_pucker_component_periodic_delta(
                primitive,
                plus.reshape(coords.shape),
                minus.reshape(coords.shape),
            )
        elif primitive.function == "D":
            delta = _periodic_delta(value_plus, value_minus)
        else:
            delta = value_plus - value_minus
        if not np.isfinite(delta) or not np.isfinite(base):
            row[idx] = np.nan
        else:
            row[idx] = delta / (2.0 * step_angstrom)
    return row


def _ring_pucker_component_periodic_delta(
    primitive: GICPrimitive,
    plus_coords: np.ndarray,
    minus_coords: np.ndarray,
) -> float:
    delta = 0.0
    for coefficient, atoms in _ring_pucker_terms_from_refs(primitive):
        delta += coefficient * _periodic_delta(
            _dihedral_value(plus_coords, atoms),
            _dihedral_value(minus_coords, atoms),
        )
    return float(delta)


class _Dual:
    __slots__ = ("val", "der")

    def __init__(self, val: float, der: np.ndarray):
        self.val = float(val)
        self.der = np.asarray(der, dtype=float)

    def _coerce(self, other: object) -> "_Dual":
        if isinstance(other, _Dual):
            return other
        return _Dual(float(other), np.zeros_like(self.der))

    def __add__(self, other: object) -> "_Dual":
        rhs = self._coerce(other)
        return _Dual(self.val + rhs.val, self.der + rhs.der)

    def __radd__(self, other: object) -> "_Dual":
        return self.__add__(other)

    def __sub__(self, other: object) -> "_Dual":
        rhs = self._coerce(other)
        return _Dual(self.val - rhs.val, self.der - rhs.der)

    def __rsub__(self, other: object) -> "_Dual":
        lhs = self._coerce(other)
        return _Dual(lhs.val - self.val, lhs.der - self.der)

    def __mul__(self, other: object) -> "_Dual":
        rhs = self._coerce(other)
        return _Dual(self.val * rhs.val, self.val * rhs.der + rhs.val * self.der)

    def __rmul__(self, other: object) -> "_Dual":
        return self.__mul__(other)

    def __truediv__(self, other: object) -> "_Dual":
        rhs = self._coerce(other)
        inv = 1.0 / rhs.val
        return _Dual(self.val * inv, (self.der - self.val * rhs.der * inv) * inv)

    def __rtruediv__(self, other: object) -> "_Dual":
        lhs = self._coerce(other)
        inv = 1.0 / self.val
        return _Dual(lhs.val * inv, (lhs.der - lhs.val * self.der * inv) * inv)

    def __neg__(self) -> "_Dual":
        return _Dual(-self.val, -self.der)


def _dual_coordinates(coords: np.ndarray) -> list[list[_Dual]]:
    flat = np.asarray(coords, dtype=float).reshape(-1)
    dim = flat.size
    out: list[list[_Dual]] = []
    for atom in range(coords.shape[0]):
        row = []
        for axis in range(3):
            idx = 3 * atom + axis
            der = np.zeros(dim, dtype=float)
            der[idx] = 1.0
            row.append(_Dual(float(coords[atom, axis]), der))
        out.append(row)
    return out


def _dual_primitive_value(
    primitive: GICPrimitive,
    dcoords: list[list[_Dual]],
    coords: np.ndarray,
) -> _Dual:
    if primitive.function == "L":
        return _dual_linear_bend_value(dcoords, primitive.atoms, mode=primitive.mode)
    if primitive.function == "D":
        return _dual_dihedral_value(dcoords, primitive.atoms)
    if primitive.function == "U":
        return _dual_out_of_plane_value(dcoords, primitive.atoms)
    if primitive.function == "FROT":
        return _dual_fragment_rotation_value(
            dcoords,
            coords,
            primitive.atoms,
            primitive.ref_atoms,
            mode=primitive.mode,
            frame_atoms=primitive.frame_atoms,
            ref_frame_atoms=primitive.ref_frame_atoms,
        )
    raise GICForgeContractError(
        f"analytic B row is not implemented for function {primitive.function}"
    )


def _dual_linear_bend_value(
    dcoords: list[list[_Dual]],
    atoms: tuple[int, ...],
    *,
    mode: int,
) -> _Dual:
    i, j, k = (atom - 1 for atom in atoms)
    left = _d_unit(_d_vec_sub(dcoords[i], dcoords[j]))
    right = _d_unit(_d_vec_sub(dcoords[k], dcoords[j]))
    axis = _d_unit(_d_vec_sub(right, left))
    e1, e2 = _d_orthogonal_frame(axis)
    bend = _d_vec_add(left, right)
    return _d_dot(bend, e1 if mode == -1 else e2)


def _dual_dihedral_value(dcoords: list[list[_Dual]], atoms: tuple[int, ...]) -> _Dual:
    i, j, k, l = (atom - 1 for atom in atoms)
    p0, p1, p2, p3 = dcoords[i], dcoords[j], dcoords[k], dcoords[l]
    b0 = _d_vec_neg(_d_vec_sub(p1, p0))
    b1 = _d_vec_sub(p2, p1)
    b2 = _d_vec_sub(p3, p2)
    b1 = _d_unit(b1)
    v = _d_vec_sub(b0, _d_vec_scale(b1, _d_dot(b0, b1)))
    w = _d_vec_sub(b2, _d_vec_scale(b1, _d_dot(b2, b1)))
    x = _d_dot(v, w)
    y = _d_dot(_d_cross(b1, v), w)
    return _d_atan2(y, x)


def _dual_out_of_plane_value(
    dcoords: list[list[_Dual]],
    atoms: tuple[int, ...],
) -> _Dual:
    center, n1, n2, n3 = (atom - 1 for atom in atoms)
    r1 = _d_vec_sub(dcoords[n1], dcoords[center])
    r2 = _d_vec_sub(dcoords[n2], dcoords[center])
    r3 = _d_vec_sub(dcoords[n3], dcoords[center])
    normal = _d_unit(_d_cross(r2, r3))
    return _d_asin(_d_dot(_d_unit(r1), normal))


def _dual_fragment_rotation_value(
    dcoords: list[list[_Dual]],
    coords: np.ndarray,
    atoms: tuple[int, ...],
    ref_atoms: tuple[int, ...],
    *,
    mode: int,
    frame_atoms: tuple[int, ...] = (),
    ref_frame_atoms: tuple[int, ...] = (),
) -> _Dual:
    frame_atoms = frame_atoms or _fragment_frame_anchor_atoms(atoms, coords=coords)
    ref_frame_atoms = ref_frame_atoms or _fragment_frame_anchor_atoms(ref_atoms, coords=coords)
    frame_frag = _d_fragment_frame(dcoords, atoms, frame_atoms=frame_atoms)
    frame_ref = _d_fragment_frame(dcoords, ref_atoms, frame_atoms=ref_frame_atoms)
    rotation = [
        [_d_dot(frame_frag[left], frame_ref[right]) for right in range(3)]
        for left in range(3)
    ]
    return _d_rotation_vector(rotation)[mode]


def _d_fragment_frame(
    dcoords: list[list[_Dual]],
    atoms: tuple[int, ...],
    *,
    frame_atoms: tuple[int, ...],
) -> list[list[_Dual]]:
    p_atom, q_atom = frame_atoms
    center = _d_fragment_center(dcoords, atoms)
    p_axis = _d_unit(_d_vec_sub(dcoords[p_atom - 1], center))
    q_raw = _d_cross(p_axis, _d_vec_sub(dcoords[q_atom - 1], center))
    q_axis = _d_unit(q_raw)
    s_axis = _d_unit(_d_cross(p_axis, q_axis))
    return [p_axis, q_axis, s_axis]


def _d_fragment_center(
    dcoords: list[list[_Dual]],
    atoms: tuple[int, ...],
) -> list[_Dual]:
    if not atoms:
        raise FloatingPointError("fragment has no atoms")
    total = [_d_zero_like(dcoords[atoms[0] - 1][0]) for _axis in range(3)]
    for atom in atoms:
        total = _d_vec_add(total, dcoords[atom - 1])
    return _d_vec_scale(total, 1.0 / len(atoms))


def _d_quaternion_vector(rotation: list[list[_Dual]]) -> tuple[_Dual, _Dual, _Dual]:
    trace = rotation[0][0] + rotation[1][1] + rotation[2][2]
    if trace.val <= -1.0 + RANK_TOLERANCE:
        raise FloatingPointError("fragment quaternion is singular near 180 degrees")
    kw = 0.5 * _d_sqrt(trace + 1.0)
    denom = 4.0 * kw
    return (
        (rotation[1][2] - rotation[2][1]) / denom,
        (rotation[2][0] - rotation[0][2]) / denom,
        (rotation[0][1] - rotation[1][0]) / denom,
    )


def _d_rotation_vector(rotation: list[list[_Dual]]) -> tuple[_Dual, _Dual, _Dual]:
    trace = rotation[0][0] + rotation[1][1] + rotation[2][2]
    if trace.val <= -1.0 + RANK_TOLERANCE:
        raise FloatingPointError("fragment exponential map is singular near 180 degrees")
    kw = 0.5 * _d_sqrt(trace + 1.0)
    kx, ky, kz = _d_quaternion_vector(rotation)
    kn2 = kx * kx + ky * ky + kz * kz
    if kn2.val <= RANK_TOLERANCE * RANK_TOLERANCE:
        return 2.0 * kx, 2.0 * ky, 2.0 * kz
    kn = _d_sqrt(kn2)
    factor = (2.0 * _d_atan2(kn, kw)) / kn
    return factor * kx, factor * ky, factor * kz


def _d_zero_like(value: _Dual) -> _Dual:
    return _Dual(0.0, np.zeros_like(value.der))


def _d_vec_add(left: list[_Dual], right: list[_Dual]) -> list[_Dual]:
    return [left[idx] + right[idx] for idx in range(3)]


def _d_vec_sub(left: list[_Dual], right: list[_Dual]) -> list[_Dual]:
    return [left[idx] - right[idx] for idx in range(3)]


def _d_vec_neg(vector: list[_Dual]) -> list[_Dual]:
    return [-item for item in vector]


def _d_vec_scale(vector: list[_Dual], scale: float | _Dual) -> list[_Dual]:
    return [scale * item for item in vector]


def _d_dot(left: list[_Dual], right: list[_Dual]) -> _Dual:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _d_cross(left: list[_Dual], right: list[_Dual]) -> list[_Dual]:
    return [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]


def _d_norm(vector: list[_Dual]) -> _Dual:
    return _d_sqrt(_d_dot(vector, vector))


def _d_unit(vector: list[_Dual]) -> list[_Dual]:
    norm = _d_norm(vector)
    if norm.val <= RANK_TOLERANCE:
        raise FloatingPointError("zero-length vector")
    return [item / norm for item in vector]


def _d_orthogonal_frame(axis: list[_Dual]) -> tuple[list[_Dual], list[_Dual]]:
    trial = [1.0, 0.0, 0.0]
    if abs(axis[0].val) > 0.9:
        trial = [0.0, 1.0, 0.0]
    e1 = _d_unit(_d_cross(axis, [_d_constant(axis[0], value) for value in trial]))
    e2 = _d_unit(_d_cross(axis, e1))
    return e1, e2


def _d_constant(template: _Dual, value: float) -> _Dual:
    return _Dual(value, np.zeros_like(template.der))


def _d_sqrt(value: _Dual) -> _Dual:
    if value.val <= 0.0:
        raise FloatingPointError("square root of non-positive value")
    root = float(np.sqrt(value.val))
    return _Dual(root, value.der / (2.0 * root))


def _d_acos(value: _Dual) -> _Dual:
    clipped = float(np.clip(value.val, -1.0, 1.0))
    denom = float(np.sqrt(max(1.0 - clipped * clipped, 0.0)))
    if denom <= RANK_TOLERANCE:
        raise FloatingPointError("acos derivative is singular")
    return _Dual(float(np.arccos(clipped)), -value.der / denom)


def _d_asin(value: _Dual) -> _Dual:
    clipped = float(np.clip(value.val, -1.0, 1.0))
    denom = float(np.sqrt(max(1.0 - clipped * clipped, 0.0)))
    if denom <= RANK_TOLERANCE:
        raise FloatingPointError("asin derivative is singular")
    return _Dual(float(np.arcsin(clipped)), value.der / denom)


def _d_atan2(y_value: _Dual, x_value: _Dual) -> _Dual:
    denom = x_value.val * x_value.val + y_value.val * y_value.val
    if denom <= RANK_TOLERANCE:
        raise FloatingPointError("atan2 derivative is singular")
    der = (x_value.val * y_value.der - y_value.val * x_value.der) / denom
    return _Dual(float(np.arctan2(y_value.val, x_value.val)), der)


def _primitive_value(primitive: GICPrimitive, coords: np.ndarray) -> float:
    if primitive.function == "R":
        return _distance_value(coords, primitive.atoms)
    if primitive.function == "A":
        return _angle_value(coords, primitive.atoms)
    if primitive.function == "L":
        return _linear_bend_value(coords, primitive.atoms, mode=primitive.mode)
    if primitive.function == "D":
        return _dihedral_value(coords, primitive.atoms)
    if primitive.function == "U":
        return _out_of_plane_value(coords, primitive.atoms)
    if primitive.function == "FC_DIST":
        return _fragment_center_distance_value(coords, primitive.atoms, primitive.ref_atoms)
    if primitive.function == "FCA_DIST":
        return _fragment_center_atom_distance_value(coords, primitive.atoms, primitive.ref_atoms)
    if primitive.function == "CENTER_ATOM_DIST":
        return _fragment_center_atom_distance_value(coords, primitive.atoms, primitive.ref_atoms)
    if primitive.function == "FTRANS":
        return _fragment_translation_value(
            coords,
            primitive.atoms,
            primitive.ref_atoms,
            mode=primitive.mode,
        )
    if primitive.function == "RPCK":
        return _ring_pucker_component_value(primitive, coords)
    if primitive.function == "FROT":
        return _fragment_rotation_value(
            coords,
            primitive.atoms,
            primitive.ref_atoms,
            mode=primitive.mode,
            frame_atoms=primitive.frame_atoms,
            ref_frame_atoms=primitive.ref_frame_atoms,
        )
    raise GICForgeContractError(f"unsupported primitive function: {primitive.function}")


def _ring_pucker_component_value(primitive: GICPrimitive, coords: np.ndarray) -> float:
    value = 0.0
    for coefficient, atoms in _ring_pucker_terms_from_refs(primitive):
        value += coefficient * _dihedral_value(coords, atoms)
    return float(value)


def _distance_value(coords: np.ndarray, atoms: tuple[int, ...]) -> float:
    i, j = (atom - 1 for atom in atoms)
    return float(np.linalg.norm(coords[i] - coords[j]))


def _angle_value(coords: np.ndarray, atoms: tuple[int, ...]) -> float:
    i, j, k = (atom - 1 for atom in atoms)
    u = coords[i] - coords[j]
    v = coords[k] - coords[j]
    return float(np.arccos(np.clip(_dot_unit(u, v), -1.0, 1.0)))


def _linear_bend_value(coords: np.ndarray, atoms: tuple[int, ...], *, mode: int) -> float:
    i, j, k = (atom - 1 for atom in atoms)
    left = _unit(coords[i] - coords[j])
    right = _unit(coords[k] - coords[j])
    axis = _unit(right - left)
    e1, e2 = _orthogonal_frame(axis)
    bend = left + right
    return float(np.dot(bend, e1 if mode == -1 else e2))


def _dihedral_value(coords: np.ndarray, atoms: tuple[int, ...]) -> float:
    i, j, k, l = (atom - 1 for atom in atoms)
    p0, p1, p2, p3 = coords[i], coords[j], coords[k], coords[l]
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = _unit(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.arctan2(y, x))


def _out_of_plane_value(coords: np.ndarray, atoms: tuple[int, ...]) -> float:
    center, n1, n2, n3 = (atom - 1 for atom in atoms)
    r1 = coords[n1] - coords[center]
    r2 = coords[n2] - coords[center]
    r3 = coords[n3] - coords[center]
    normal = _unit(np.cross(r2, r3))
    return float(np.arcsin(np.clip(np.dot(_unit(r1), normal), -1.0, 1.0)))


def _fragment_center_distance_value(
    coords: np.ndarray,
    atoms: tuple[int, ...],
    ref_atoms: tuple[int, ...],
) -> float:
    delta = _fragment_center(coords, atoms) - _fragment_center(coords, ref_atoms)
    return float(np.linalg.norm(delta))


def _fragment_center_atom_distance_value(
    coords: np.ndarray,
    atoms: tuple[int, ...],
    ref_atoms: tuple[int, ...],
) -> float:
    if len(ref_atoms) != 1:
        raise FloatingPointError("center-atom distance needs exactly one reference atom")
    return float(np.linalg.norm(_fragment_center(coords, atoms) - coords[ref_atoms[0] - 1]))


def _fragment_translation_value(
    coords: np.ndarray,
    atoms: tuple[int, ...],
    ref_atoms: tuple[int, ...],
    *,
    mode: int,
) -> float:
    delta = _fragment_center(coords, atoms) - _fragment_center(coords, ref_atoms)
    return float(delta[mode])


def _fragment_rotation_value(
    coords: np.ndarray,
    atoms: tuple[int, ...],
    ref_atoms: tuple[int, ...],
    *,
    mode: int,
    frame_atoms: tuple[int, ...] = (),
    ref_frame_atoms: tuple[int, ...] = (),
) -> float:
    frame_frag = _fragment_frame(coords, atoms, frame_atoms=frame_atoms)
    frame_ref = _fragment_frame(coords, ref_atoms, frame_atoms=ref_frame_atoms)
    rotation = frame_frag.T @ frame_ref
    return float(_rotation_vector(rotation)[mode])


def _cartesian_vibrational_basis(coords: np.ndarray, *, target_rank: int) -> np.ndarray:
    flat_dim = int(coords.shape[0] * 3)
    if target_rank == 0:
        return np.zeros((flat_dim, 0), dtype=float)
    centered = coords - np.mean(coords, axis=0)
    rows: list[np.ndarray] = []
    for axis in range(3):
        row = np.zeros((coords.shape[0], 3), dtype=float)
        row[:, axis] = 1.0
        rows.append(row.reshape(-1))
    for axis_vector in np.eye(3):
        row = np.cross(axis_vector, centered).reshape(-1)
        if np.linalg.norm(row) > RANK_TOLERANCE:
            rows.append(row)
    external = np.asarray(rows, dtype=float)
    _, singular_values, vh = np.linalg.svd(external, full_matrices=True)
    rank = int(np.sum(singular_values > RANK_TOLERANCE))
    basis = vh[rank:].T
    if basis.shape[1] < target_rank:
        raise GICForgeContractError(
            f"cannot build SYCART basis: need {target_rank}, got {basis.shape[1]}"
        )
    return basis[:, :target_rank]


def _fragment_primitive_candidates(
    fragment_records: tuple[object, ...],
    *,
    coords: np.ndarray,
    counters: dict[str, int],
) -> tuple[GICPrimitive, ...]:
    if len(fragment_records) <= 1:
        return ()
    records = sorted(fragment_records, key=lambda item: getattr(item, "identifier"))
    reference = max(
        records,
        key=lambda item: (len(getattr(item, "atoms")), -_fragment_number(item)),
    )
    candidates: list[GICPrimitive] = []

    for left, right in combinations(records, 2):
        candidates.append(
            _make_primitive(
                "FRAG_DISTANCE",
                "FC_DIST",
                tuple(getattr(left, "atoms")),
                counters,
                ref_atoms=tuple(getattr(right, "atoms")),
                refs=(getattr(left, "identifier"), getattr(right, "identifier")),
            )
        )

    anchor_atoms = _fragment_anchor_atoms(reference, coords=coords)
    for record in records:
        if getattr(record, "identifier") == getattr(reference, "identifier"):
            continue
        atoms = tuple(getattr(record, "atoms"))
        ref_atoms = tuple(getattr(reference, "atoms"))
        for atom in anchor_atoms:
            candidates.append(
                _make_primitive(
                    "FRAG_CENTER_ATOM_DISTANCE",
                    "FCA_DIST",
                    atoms,
                    counters,
                    ref_atoms=(atom,),
                    refs=(getattr(record, "identifier"), f"A{atom}"),
                )
            )
        for axis in range(3):
            candidates.append(
                _make_primitive(
                    "FRAG_TRANSLATION",
                    "FTRANS",
                    atoms,
                    counters,
                    mode=axis,
                    ref_atoms=ref_atoms,
                    refs=(getattr(record, "identifier"), getattr(reference, "identifier")),
                )
            )
        if (
            _fragment_frame_rank(coords, atoms) >= 2
            and _fragment_frame_rank(coords, ref_atoms) >= 2
        ):
            frame_atoms = _fragment_frame_anchor_atoms(atoms, coords=coords)
            ref_frame_atoms = _fragment_frame_anchor_atoms(ref_atoms, coords=coords)
            for axis in range(3):
                candidates.append(
                    _make_primitive(
                        "FRAG_ORIENTATION",
                        "FROT",
                        atoms,
                        counters,
                        mode=axis,
                        ref_atoms=ref_atoms,
                        refs=(getattr(record, "identifier"), getattr(reference, "identifier")),
                        frame_atoms=frame_atoms,
                        ref_frame_atoms=ref_frame_atoms,
                    )
                )
    return tuple(candidates)


def _interaction_center_primitive_candidates(
    definition: object | None,
    *,
    counters: dict[str, int],
) -> tuple[GICPrimitive, ...]:
    if definition is None:
        return ()
    centers = {
        getattr(center, "identifier"): center
        for center in getattr(definition, "centers", ())
    }
    candidates: list[GICPrimitive] = []
    for interaction in getattr(definition, "interactions", ()):
        center = centers.get(getattr(interaction, "center_id"))
        if center is None:
            continue
        atom = int(getattr(interaction, "atom"))
        atoms = tuple(int(item) for item in getattr(center, "atoms"))
        candidates.append(
            _make_primitive(
                "CENTER_ATOM_DISTANCE",
                "CENTER_ATOM_DIST",
                atoms,
                counters,
                ref_atoms=(atom,),
                refs=(getattr(center, "identifier"), f"A{atom}"),
            )
        )
    return tuple(candidates)


def _fragment_records(path: Path) -> tuple[object, ...]:
    try:
        from oracle_fragments import read_fragment_records
    except ImportError:
        return ()
    return tuple(read_fragment_records(Path(path)))


def _interaction_center_definition(path: Path) -> object | None:
    try:
        from oracle_fragments import read_interaction_center_definition
    except ImportError:
        return None
    return read_interaction_center_definition(Path(path))


def _fragment_number(record: object) -> int:
    identifier = str(getattr(record, "identifier"))
    try:
        return int(identifier[1:])
    except ValueError:
        return 0


def _fragment_anchor_atoms(
    record: object,
    *,
    coords: np.ndarray,
    limit: int = 3,
) -> tuple[int, ...]:
    atoms = tuple(getattr(record, "atoms"))
    if len(atoms) <= limit:
        return atoms
    center = _fragment_center(coords, atoms)
    ranked = sorted(
        atoms,
        key=lambda atom: (-float(np.linalg.norm(coords[atom - 1] - center)), atom),
    )
    return tuple(sorted(ranked[:limit]))


def _fragment_center(coords: np.ndarray, atoms: tuple[int, ...]) -> np.ndarray:
    if not atoms:
        raise FloatingPointError("fragment has no atoms")
    return np.mean(coords[[atom - 1 for atom in atoms]], axis=0)


def _fragment_frame_rank(coords: np.ndarray, atoms: tuple[int, ...]) -> int:
    if len(atoms) < 2:
        return 0
    centered = coords[[atom - 1 for atom in atoms]] - _fragment_center(coords, atoms)
    singular_values = np.linalg.svd(centered, compute_uv=False)
    return int(np.sum(singular_values > 1.0e-8))


def _fragment_frame(
    coords: np.ndarray,
    atoms: tuple[int, ...],
    *,
    frame_atoms: tuple[int, ...] = (),
) -> np.ndarray:
    if _fragment_frame_rank(coords, atoms) < 2:
        raise FloatingPointError("fragment orientation is underdefined")
    p_atom, q_atom = (
        tuple(frame_atoms)
        if frame_atoms
        else _fragment_frame_anchor_atoms(atoms, coords=coords)
    )
    center = _fragment_center(coords, atoms)
    p_axis = _unit(coords[p_atom - 1] - center)
    q_raw = np.cross(p_axis, coords[q_atom - 1] - center)
    q_axis = _unit(q_raw)
    s_axis = _unit(np.cross(p_axis, q_axis))
    return np.column_stack([p_axis, q_axis, s_axis])


def _fragment_frame_anchor_atoms(
    atoms: tuple[int, ...],
    *,
    coords: np.ndarray,
) -> tuple[int, int]:
    center = _fragment_center(coords, atoms)
    ranked = sorted(
        atoms,
        key=lambda atom: (-float(np.linalg.norm(coords[atom - 1] - center)), atom),
    )
    p_atom = ranked[0]
    p_axis = _unit(coords[p_atom - 1] - center)
    q_candidates = []
    for atom in atoms:
        if atom == p_atom:
            continue
        vector = coords[atom - 1] - center
        norm = float(np.linalg.norm(vector))
        if norm <= RANK_TOLERANCE:
            continue
        dot = abs(float(np.dot(p_axis, vector / norm)))
        q_candidates.append((dot, -norm, atom))
    if not q_candidates:
        raise FloatingPointError("fragment has no second orientation anchor")
    _dot, _norm, q_atom = min(q_candidates)
    return p_atom, q_atom


def _canonicalize_frame(frame: np.ndarray) -> np.ndarray:
    out = np.array(frame, dtype=float, copy=True)
    for axis in range(3):
        column = out[:, axis]
        pivot = int(np.argmax(np.abs(column)))
        if column[pivot] < 0.0:
            out[:, axis] *= -1.0
    if np.linalg.det(out) < 0.0:
        out[:, -1] *= -1.0
    return out


def _rotation_vector(rotation: np.ndarray) -> np.ndarray:
    trace = float(np.trace(rotation))
    if trace > 0.0:
        scale = 0.5 / np.sqrt(trace + 1.0)
        qw = 0.25 / scale
        qx = (rotation[1, 2] - rotation[2, 1]) * scale
        qy = (rotation[2, 0] - rotation[0, 2]) * scale
        qz = (rotation[0, 1] - rotation[1, 0]) * scale
    else:
        qw, qx, qy, qz = _fallback_quaternion(rotation)
    quat = np.array([qw, qx, qy, qz], dtype=float)
    if quat[0] < 0.0:
        quat = -quat
    norm = float(np.linalg.norm(quat))
    if norm <= 1.0e-12:
        return np.zeros(3, dtype=float)
    quat /= norm
    vector = quat[1:]
    vector_norm = float(np.linalg.norm(vector))
    if vector_norm < 1.0e-12:
        return np.zeros(3, dtype=float)
    angle = 2.0 * np.arctan2(vector_norm, quat[0])
    return vector / vector_norm * angle


def _quaternion_vector(rotation: np.ndarray) -> np.ndarray:
    trace = float(np.trace(rotation))
    if trace > 0.0:
        scale = 0.5 / np.sqrt(trace + 1.0)
        qw = 0.25 / scale
        qx = (rotation[1, 2] - rotation[2, 1]) * scale
        qy = (rotation[2, 0] - rotation[0, 2]) * scale
        qz = (rotation[0, 1] - rotation[1, 0]) * scale
    else:
        qw, qx, qy, qz = _fallback_quaternion(rotation)
    quat = np.array([qw, qx, qy, qz], dtype=float)
    if quat[0] < 0.0:
        quat = -quat
    norm = float(np.linalg.norm(quat))
    if norm <= 1.0e-12:
        return np.zeros(3, dtype=float)
    quat /= norm
    return quat[1:]


def _fallback_quaternion(rotation: np.ndarray) -> tuple[float, float, float, float]:
    if rotation[0, 0] > rotation[1, 1] and rotation[0, 0] > rotation[2, 2]:
        scale = 2.0 * np.sqrt(max(0.0, 1.0 + rotation[0, 0] - rotation[1, 1] - rotation[2, 2]))
        return (
            (rotation[2, 1] - rotation[1, 2]) / scale,
            0.25 * scale,
            (rotation[0, 1] + rotation[1, 0]) / scale,
            (rotation[0, 2] + rotation[2, 0]) / scale,
        )
    if rotation[1, 1] > rotation[2, 2]:
        scale = 2.0 * np.sqrt(max(0.0, 1.0 + rotation[1, 1] - rotation[0, 0] - rotation[2, 2]))
        return (
            (rotation[0, 2] - rotation[2, 0]) / scale,
            (rotation[0, 1] + rotation[1, 0]) / scale,
            0.25 * scale,
            (rotation[1, 2] + rotation[2, 1]) / scale,
        )
    scale = 2.0 * np.sqrt(max(0.0, 1.0 + rotation[2, 2] - rotation[0, 0] - rotation[1, 1]))
    return (
        (rotation[1, 0] - rotation[0, 1]) / scale,
        (rotation[0, 2] + rotation[2, 0]) / scale,
        (rotation[1, 2] + rotation[2, 1]) / scale,
        0.25 * scale,
    )


def _vibrational_rank(coords: np.ndarray) -> int:
    natoms = int(coords.shape[0])
    if natoms <= 1:
        return 0
    return max(0, 3 * natoms - (5 if _is_linear(coords) else 6))


def _is_linear(coords: np.ndarray) -> bool:
    if coords.shape[0] <= 2:
        return True
    centered = coords - np.mean(coords, axis=0)
    singular_values = np.linalg.svd(centered, compute_uv=False)
    if singular_values.size < 2 or singular_values[0] <= RANK_TOLERANCE:
        return False
    return bool(singular_values[1] / singular_values[0] <= 1.0e-6)


def _adjacency(bonds: tuple[tuple[int, int], ...], *, natoms: int) -> dict[int, set[int]]:
    graph = {idx: set() for idx in range(1, natoms + 1)}
    for i, j in bonds:
        graph[i].add(j)
        graph[j].add(i)
    return graph


def _point_group(lines: list[str]) -> str:
    for line in section_content(lines, "SYMMETRY"):
        parts = line.split()
        if len(parts) >= 2 and parts[0].upper() == "POINT_GROUP":
            return parts[1]
    return "UNKNOWN"


def _symmetry_operations(lines: list[str]) -> tuple[GICPointGroupOperation, ...]:
    symmetry = section_content(lines, "SYMMETRY")
    if not symmetry:
        return ()
    operation_lines = _subsection(symmetry, "OPERATIONS")
    if not operation_lines:
        return ()
    operations: list[GICPointGroupOperation] = []
    for line in operation_lines:
        text = line.strip()
        if not text or text.upper() == "NONE":
            continue
        parts = text.split()
        fields = _key_values(parts[1:])
        try:
            matrix_values = _parse_float_list(fields["MATRIX"])
            if len(matrix_values) != 9:
                raise ValueError("operation matrix must have 9 values")
            permutation = _parse_atom_list(fields["PERMUTATION"])
            operations.append(
                GICPointGroupOperation(
                    label=fields["LABEL"],
                    rotation=tuple(
                        tuple(float(value) for value in matrix_values[start : start + 3])
                        for start in (0, 3, 6)
                    ),
                    permutation=permutation,
                )
            )
        except KeyError as exc:
            raise GICForgeContractError(
                f"invalid #SYMMETRY operation line: {line}"
            ) from exc
        except ValueError as exc:
            raise GICForgeContractError(
                f"invalid #SYMMETRY operation numeric field: {line}"
            ) from exc
    return tuple(operations)


def _subsection(section_lines: list[str], name: str) -> list[str]:
    header = f"[{name.upper()}]"
    start = None
    for idx, line in enumerate(section_lines):
        if line.strip().upper() == header:
            start = idx + 1
            break
    if start is None:
        return []
    end = len(section_lines)
    for idx in range(start, len(section_lines)):
        text = section_lines[idx].strip()
        if text.startswith("[") and text.endswith("]"):
            end = idx
            break
    return list(section_lines[start:end])


def _section_value(section_lines: list[str], key: str) -> str | None:
    key_upper = key.upper()
    for line in section_lines:
        parts = line.split(maxsplit=1)
        if len(parts) == 2 and parts[0].upper() == key_upper:
            return parts[1].strip()
    return None


def _parse_primitive_line(line: str) -> GICPrimitive:
    parts = line.split()
    if not parts:
        raise GICForgeContractError("empty primitive line")
    fields = _key_values(parts[1:])
    try:
        return GICPrimitive(
            identifier=parts[0],
            name=fields["NAME"],
            family=fields["FAMILY"],
            function=fields["FUNCTION"],
            atoms=_parse_atom_list(fields["ATOMS"]),
            mode=int(fields.get("MODE", "0")),
            ref_atoms=_parse_atom_list(fields.get("REF_ATOMS", "")),
            refs=_parse_text_list(fields.get("REFS", "")),
            frame_atoms=_parse_atom_list(fields.get("FRAME_ATOMS", "")),
            ref_frame_atoms=_parse_atom_list(fields.get("REF_FRAME_ATOMS", "")),
        )
    except KeyError as exc:
        raise GICForgeContractError(f"invalid primitive line: {line}") from exc
    except ValueError as exc:
        raise GICForgeContractError(f"invalid primitive numeric field: {line}") from exc


def _parse_frozen_gic_line(line: str) -> FrozenGIC:
    parts = line.split()
    if not parts:
        raise GICForgeContractError("empty frozen GIC line")
    fields = _key_values(parts[1:])
    coefficients = _parse_coefficients(fields.get("COEFFS", ""))
    if not coefficients:
        primitive_id = fields.get("PRIMITIVE")
        if not primitive_id:
            raise GICForgeContractError(f"invalid frozen GIC coefficients: {line}")
        coefficients = ((primitive_id, 1.0),)
    try:
        return FrozenGIC(
            identifier=parts[0],
            name=fields["NAME"],
            family=fields["FAMILY"],
            irrep=fields["IRREP"],
            primitive_id=coefficients[0][0],
            gaussian_expression=fields.get("GAUSSIAN", "NONE"),
            coefficients=coefficients,
        )
    except KeyError as exc:
        raise GICForgeContractError(f"invalid frozen GIC line: {line}") from exc


def _parse_reduction_diagnostics(
    section_lines: list[str],
    *,
    selected: tuple[str, ...],
) -> GICReductionDiagnostics:
    lines = _subsection(section_lines, "REDUCTION_DIAGNOSTICS")
    if not lines:
        return GICReductionDiagnostics(
            rank_method=_section_value(section_lines, "RANK_METHOD") or RANK_METHOD,
            reduction_policy=_section_value(section_lines, "REDUCTION_POLICY")
            or REDUCTION_POLICY,
            selected=selected,
        )
    return GICReductionDiagnostics(
        rank_method=_section_value(lines, "RANK_METHOD") or RANK_METHOD,
        reduction_policy=_section_value(lines, "REDUCTION_POLICY") or REDUCTION_POLICY,
        selected=_parse_text_list(_section_value(lines, "SELECTED") or "") or selected,
        skipped_singular=_parse_text_list(_section_value(lines, "SKIPPED_SINGULAR") or ""),
        skipped_dependent=_parse_text_list(_section_value(lines, "SKIPPED_DEPENDENT") or ""),
    )


def _parse_symmetry_diagnostics(
    section_lines: list[str],
) -> GICSymmetrizationDiagnostics | None:
    lines = _subsection(section_lines, "SYMMETRY_DIAGNOSTICS")
    if not lines:
        return None
    groups: list[GICSymmetrizedGroup] = []
    for line in lines:
        parts = line.split()
        if not parts or parts[0].upper() != "GROUP":
            continue
        fields = _key_values(parts[2:] if len(parts) > 1 else parts[1:])
        try:
            groups.append(
                GICSymmetrizedGroup(
                    block=fields["BLOCK"],
                    family=fields["FAMILY"],
                    signature=fields["SIGNATURE"],
                    source_gics=_parse_text_list(fields.get("SOURCES", "")),
                    output_gics=_parse_text_list(fields.get("OUTPUTS", "")),
                )
            )
        except KeyError as exc:
            raise GICForgeContractError(
                f"invalid symmetry diagnostic group line: {line}"
            ) from exc
    return GICSymmetrizationDiagnostics(
        method=_section_value(lines, "METHOD") or "UNKNOWN",
        policy=_section_value(lines, "POLICY") or SYMMETRIZATION_POLICY,
        status=_section_value(lines, "STATUS") or "UNKNOWN",
        point_group=_section_value(lines, "POINT_GROUP")
        or _section_value(section_lines, "POINT_GROUP")
        or "UNKNOWN",
        symmetry_group=_section_value(lines, "SYMMETRY_GROUP")
        or _section_value(section_lines, "SYMMETRY_GROUP")
        or _section_value(lines, "POINT_GROUP")
        or _section_value(section_lines, "POINT_GROUP")
        or "UNKNOWN",
        total_symmetric_irrep=_section_value(lines, "TOTAL_SYMMETRIC_IRREP")
        or _section_value(section_lines, "TOTAL_SYMMETRIC_IRREP")
        or total_symmetric_irrep(_section_value(section_lines, "POINT_GROUP")),
        total_symmetric_gics=_parse_text_list(
            _section_value(lines, "TOTAL_SYMMETRIC_GICS")
            or _section_value(section_lines, "TOTAL_SYMMETRIC_GICS")
            or ""
        ),
        groups=tuple(groups),
    )


def _key_values(parts: list[str]) -> dict[str, str]:
    fields: dict[str, str] = {}
    for part in parts:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        fields[key.upper()] = value
    return fields


def _parse_atom_list(text: str) -> tuple[int, ...]:
    if not text:
        return ()
    try:
        return tuple(int(item) for item in text.split(",") if item)
    except ValueError as exc:
        raise GICForgeContractError(f"invalid atom list: {text}") from exc


def _parse_float_list(text: str) -> tuple[float, ...]:
    if not text:
        return ()
    try:
        return tuple(float(item) for item in text.split(",") if item)
    except ValueError as exc:
        raise GICForgeContractError(f"invalid float list: {text}") from exc


def _parse_text_list(text: str) -> tuple[str, ...]:
    if not text or text.upper() == "NONE":
        return ()
    return tuple(item for item in text.split(",") if item)


def _csv_or_none(values: tuple[str, ...]) -> str:
    return ",".join(values) if values else "NONE"


def _parse_coefficients(text: str) -> tuple[tuple[str, float], ...]:
    if not text:
        return ()
    coefficients: list[tuple[str, float]] = []
    for item in text.replace(";", ",").split(","):
        if not item:
            continue
        if ":" not in item:
            raise GICForgeContractError(f"invalid GIC coefficient: {item}")
        primitive_id, value = item.split(":", 1)
        try:
            coefficients.append((primitive_id, float(value)))
        except ValueError as exc:
            raise GICForgeContractError(f"invalid GIC coefficient value: {item}") from exc
    return tuple(coefficients)


def _parse_bool(text: str | None) -> bool:
    return bool((text or "").strip().upper() == "TRUE")


def _parse_int(text: str | None) -> int:
    if text is None:
        return 0
    try:
        return int(text)
    except ValueError as exc:
        raise GICForgeContractError(f"invalid integer field: {text}") from exc


def _cartesian_column_labels(natoms: int) -> tuple[str, ...]:
    axes = ("X", "Y", "Z")
    return tuple(f"{atom}:{axis}" for atom in range(1, natoms + 1) for axis in axes)


def _primitive_line(primitive: GICPrimitive) -> str:
    atoms = ",".join(str(atom) for atom in primitive.atoms)
    mode = f" MODE={primitive.mode}" if primitive.function == "L" else ""
    if primitive.function in {"FTRANS", "FROT"}:
        mode = f" MODE={primitive.mode}"
    ref_atoms = (
        " REF_ATOMS=" + ",".join(str(atom) for atom in primitive.ref_atoms)
        if primitive.ref_atoms
        else ""
    )
    refs = " REFS=" + ",".join(primitive.refs) if primitive.refs else ""
    frame_atoms = (
        " FRAME_ATOMS=" + ",".join(str(atom) for atom in primitive.frame_atoms)
        if primitive.frame_atoms
        else ""
    )
    ref_frame_atoms = (
        " REF_FRAME_ATOMS=" + ",".join(str(atom) for atom in primitive.ref_frame_atoms)
        if primitive.ref_frame_atoms
        else ""
    )
    return (
        f"{primitive.identifier} NAME={primitive.name} FAMILY={primitive.family} "
        f"CLASS={primitive.reduction_class} FUNCTION={primitive.function} "
        f"ATOMS={atoms}{ref_atoms}{refs}"
        f"{frame_atoms}{ref_frame_atoms}{mode} "
        f"GAUSSIAN={primitive.gaussian_expression()}"
    )


def _gaussian_gic_block_lines(definition: GICDefinition) -> list[str]:
    coords = np.asarray(definition.reference_coordinates_angstrom, dtype=float)
    lines: list[str] = []
    virtual_centers = _gaussian_virtual_center_atoms(definition.primitives)
    if virtual_centers:
        for center_id, atoms in sorted(virtual_centers.items()):
            lines.extend(_gaussian_virtual_center_lines(center_id, atoms))
    fragment_atoms = _gaussian_fragment_atoms(definition.primitives)
    if fragment_atoms:
        for fragment_id, atoms in sorted(fragment_atoms.items()):
            lines.append(f"{fragment_id}=Fragment({_atom_interval(atoms)})")
        for fragment_id in sorted(fragment_atoms):
            lines.extend(_gaussian_center_lines(fragment_id))
        frot_pairs = _gaussian_frot_pairs(definition.primitives)
        frame_atoms = _gaussian_fragment_frame_atoms(definition.primitives)
        frame_fragments = sorted({fragment for pair in frot_pairs for fragment in pair})
        for fragment_id in frame_fragments:
            atoms = fragment_atoms[fragment_id]
            lines.extend(
                _gaussian_frame_lines(
                    fragment_id,
                    atoms,
                    coords=coords,
                    frame_atoms=frame_atoms.get(fragment_id, ()),
                )
            )
        for frag_id, ref_id in sorted(frot_pairs):
            lines.extend(_gaussian_quaternion_lines(frag_id, ref_id))

    for gic in definition.gics:
        expression = _gaussian_expression_for_gic(definition, gic)
        if expression:
            label = _gaussian_label_for_gic(definition, gic)
            if gic.family == "RING_PUCKER_COMPONENT":
                label = f"{label}(Inactive)"
            lines.append(f"{label} = {expression}")
    lines.extend(_gaussian_ring_puckering_function_lines(definition))
    return lines


def _gaussian_label_for_gic(definition: GICDefinition, gic: FrozenGIC) -> str:
    if gic.family == "RING_PUCKER_COMPONENT":
        return gic.name
    return gic.name if definition.symmetrize else gic.identifier


def _gaussian_ring_puckering_function_lines(definition: GICDefinition) -> list[str]:
    if definition.symmetrize:
        return _gaussian_symmetrized_ring_puckering_function_lines(definition)
    primitive_by_id = {primitive.identifier: primitive for primitive in definition.primitives}
    groups: dict[tuple[int, ...], list[tuple[str, GICPrimitive]]] = {}
    for gic in definition.gics:
        if gic.family != "RING_PUCKER_COMPONENT":
            continue
        primitive = _single_source_primitive(gic, primitive_by_id)
        if primitive is None or primitive.function != "RPCK":
            continue
        groups.setdefault(primitive.atoms, []).append(
            (_gaussian_label_for_gic(definition, gic), primitive)
        )

    lines: list[str] = []
    pair_index = 0
    for components in groups.values():
        if len(components) == 1 and len(components[0][1].atoms) == 4:
            label = components[0][0]
            pair_index += 1
            lines.append(f"QPck{pair_index:04d} = SQRT({label}*{label})")
            lines.append(f"PhiP{pair_index:04d} = {label}")
            continue
        component_index = 0
        while component_index + 1 < len(components):
            left = components[component_index][0]
            right = components[component_index + 1][0]
            pair_index += 1
            lines.append(
                f"QPck{pair_index:04d} = SQRT({left}*{left}+{right}*{right})"
            )
            lines.append(f"PhiP{pair_index:04d} = ATAN2({right},{left})")
            component_index += 2
    return lines


def _gaussian_symmetrized_ring_puckering_function_lines(
    definition: GICDefinition,
) -> list[str]:
    groups: dict[str, list[str]] = {}
    for gic in definition.gics:
        if gic.family != "RING_PUCKER_COMPONENT":
            continue
        groups.setdefault(gic.irrep, []).append(_gaussian_label_for_gic(definition, gic))
    lines: list[str] = []
    pair_index = 0
    for labels in groups.values():
        component_index = 0
        while component_index + 1 < len(labels):
            left = labels[component_index]
            right = labels[component_index + 1]
            pair_index += 1
            lines.append(
                f"QPck{pair_index:04d} = SQRT({left}*{left}+{right}*{right})"
            )
            lines.append(f"PhiP{pair_index:04d} = ATAN2({right},{left})")
            component_index += 2
    return lines


def _gaussian_expression_for_gic(
    definition: GICDefinition,
    gic: FrozenGIC,
) -> str | None:
    coefficients = gic.coefficients or ((gic.primitive_id, 1.0),)
    primitive_by_id = {primitive.identifier: primitive for primitive in definition.primitives}
    if len(coefficients) == 1 and abs(coefficients[0][1] - 1.0) <= 1.0e-12:
        primitive = primitive_by_id.get(coefficients[0][0])
        if primitive is None:
            raise GICForgeContractError(
                f"unknown primitive {coefficients[0][0]!r} in frozen GIC {gic.identifier}"
            )
        return _gaussian_expression_for_primitive(primitive)

    terms: list[str] = []
    for primitive_id, coefficient in coefficients:
        primitive = primitive_by_id.get(primitive_id)
        if primitive is None:
            raise GICForgeContractError(
                f"unknown primitive {primitive_id!r} in frozen GIC {gic.identifier}"
            )
        expression = _gaussian_expression_for_primitive(primitive)
        if expression is None:
            return None
        terms.append(_gaussian_linear_term(coefficient, expression, first=not terms))
    return "".join(terms) if terms else None


def _gaussian_linear_term(coefficient: float, expression: str, *, first: bool) -> str:
    sign = "-" if coefficient < 0.0 else "+"
    magnitude = abs(float(coefficient))
    coefficient_text = f"{magnitude:.12g}"
    term = f"{coefficient_text}*({_strip_gaussian_outer_parentheses(expression)})"
    if first:
        return f"-{term}" if sign == "-" else term
    return f"{sign}{term}"


def _strip_gaussian_outer_parentheses(expression: str) -> str:
    text = expression.strip()
    if text.startswith("(") and text.endswith(")"):
        return text[1:-1]
    return text


def _gaussian_expression_for_primitive(primitive: GICPrimitive) -> str | None:
    if primitive.is_gaussian_native:
        return primitive.gaussian_expression()
    if primitive.function == "FC_DIST":
        frag_id, ref_id = primitive.refs
        return (
            f"SQRT((Cx{frag_id}-Cx{ref_id})**2+"
            f"(Cy{frag_id}-Cy{ref_id})**2+"
            f"(Cz{frag_id}-Cz{ref_id})**2)"
        )
    if primitive.function == "FCA_DIST":
        frag_id, atom_ref = primitive.refs
        atom = atom_ref[1:]
        return (
            f"SQRT((Cx{frag_id}-X({atom}))**2+"
            f"(Cy{frag_id}-Y({atom}))**2+"
            f"(Cz{frag_id}-Z({atom}))**2)"
        )
    if primitive.function == "CENTER_ATOM_DIST":
        center_id, atom_ref = primitive.refs
        atom = atom_ref[1:]
        return (
            f"SQRT((Cx{center_id}-X({atom}))**2+"
            f"(Cy{center_id}-Y({atom}))**2+"
            f"(Cz{center_id}-Z({atom}))**2)"
        )
    if primitive.function == "FTRANS":
        frag_id, ref_id = primitive.refs
        axis = ("x", "y", "z")[primitive.mode]
        return f"C{axis}{frag_id}-C{axis}{ref_id}"
    if primitive.function == "FROT":
        frag_id, ref_id = primitive.refs
        axis = ("x", "y", "z")[primitive.mode]
        return f"E{axis}{frag_id}{ref_id}"
    if primitive.function == "RPCK":
        terms: list[str] = []
        for coefficient, atoms in _ring_pucker_terms_from_refs(primitive):
            expression = "D(" + ",".join(str(atom) for atom in atoms) + ")"
            terms.append(_gaussian_linear_term(coefficient, expression, first=not terms))
        return "".join(terms) if terms else None
    return None


def _gaussian_fragment_atoms(
    primitives: tuple[GICPrimitive, ...],
) -> dict[str, tuple[int, ...]]:
    fragments: dict[str, tuple[int, ...]] = {}
    for primitive in primitives:
        if not primitive.refs:
            continue
        first = primitive.refs[0]
        if first.startswith("F"):
            fragments.setdefault(first, primitive.atoms)
        if len(primitive.refs) < 2:
            continue
        second = primitive.refs[1]
        if second.startswith("F"):
            fragments.setdefault(second, primitive.ref_atoms)
    return fragments


def _gaussian_virtual_center_atoms(
    primitives: tuple[GICPrimitive, ...],
) -> dict[str, tuple[int, ...]]:
    centers: dict[str, tuple[int, ...]] = {}
    for primitive in primitives:
        if primitive.function != "CENTER_ATOM_DIST" or not primitive.refs:
            continue
        center_id = primitive.refs[0]
        if center_id.startswith("C"):
            centers.setdefault(center_id, primitive.atoms)
    return centers


def _gaussian_frot_pairs(primitives: tuple[GICPrimitive, ...]) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for primitive in primitives:
        if primitive.function == "FROT" and len(primitive.refs) >= 2:
            pairs.add((primitive.refs[0], primitive.refs[1]))
    return pairs


def _gaussian_fragment_frame_atoms(
    primitives: tuple[GICPrimitive, ...],
) -> dict[str, tuple[int, ...]]:
    frame_atoms: dict[str, tuple[int, ...]] = {}
    for primitive in primitives:
        if primitive.function != "FROT" or len(primitive.refs) < 2:
            continue
        if primitive.frame_atoms:
            frame_atoms.setdefault(primitive.refs[0], primitive.frame_atoms)
        if primitive.ref_frame_atoms:
            frame_atoms.setdefault(primitive.refs[1], primitive.ref_frame_atoms)
    return frame_atoms


def _gaussian_center_lines(fragment_id: str) -> list[str]:
    return [
        f"Cx{fragment_id}(Inactive)=XCntr({fragment_id})",
        f"Cy{fragment_id}(Inactive)=YCntr({fragment_id})",
        f"Cz{fragment_id}(Inactive)=ZCntr({fragment_id})",
    ]


def _gaussian_virtual_center_lines(center_id: str, atoms: tuple[int, ...]) -> list[str]:
    if not atoms:
        return []
    denominator = len(atoms)
    return [
        f"Cx{center_id}(Inactive)=({_gaussian_axis_sum('X', atoms)})/{denominator}",
        f"Cy{center_id}(Inactive)=({_gaussian_axis_sum('Y', atoms)})/{denominator}",
        f"Cz{center_id}(Inactive)=({_gaussian_axis_sum('Z', atoms)})/{denominator}",
    ]


def _gaussian_axis_sum(axis: str, atoms: tuple[int, ...]) -> str:
    return "+".join(f"{axis}({atom})" for atom in atoms)


def _gaussian_frame_lines(
    fragment_id: str,
    atoms: tuple[int, ...],
    *,
    coords: np.ndarray,
    frame_atoms: tuple[int, ...] = (),
) -> list[str]:
    p_atom, q_atom = (
        tuple(frame_atoms)
        if frame_atoms
        else _fragment_frame_anchor_atoms(atoms, coords=coords)
    )
    return [
        f"RP{fragment_id}(Inactive)=SQRT((X({p_atom})-Cx{fragment_id})**2+"
        f"(Y({p_atom})-Cy{fragment_id})**2+(Z({p_atom})-Cz{fragment_id})**2)",
        f"Px{fragment_id}(Inactive)=[X({p_atom})-Cx{fragment_id}]/RP{fragment_id}",
        f"Py{fragment_id}(Inactive)=[Y({p_atom})-Cy{fragment_id}]/RP{fragment_id}",
        f"Pz{fragment_id}(Inactive)=[Z({p_atom})-Cz{fragment_id}]/RP{fragment_id}",
        f"QQx{fragment_id}(Inactive)=Py{fragment_id}*[Z({q_atom})-Cz{fragment_id}]-"
        f"Pz{fragment_id}*[Y({q_atom})-Cy{fragment_id}]",
        f"QQy{fragment_id}(Inactive)=Pz{fragment_id}*[X({q_atom})-Cx{fragment_id}]-"
        f"Px{fragment_id}*[Z({q_atom})-Cz{fragment_id}]",
        f"QQz{fragment_id}(Inactive)=Px{fragment_id}*[Y({q_atom})-Cy{fragment_id}]-"
        f"Py{fragment_id}*[X({q_atom})-Cx{fragment_id}]",
        f"RQ{fragment_id}(Inactive)=SQRT(QQx{fragment_id}**2+"
        f"QQy{fragment_id}**2+QQz{fragment_id}**2)",
        f"Qx{fragment_id}(Inactive)=QQx{fragment_id}/RQ{fragment_id}",
        f"Qy{fragment_id}(Inactive)=QQy{fragment_id}/RQ{fragment_id}",
        f"Qz{fragment_id}(Inactive)=QQz{fragment_id}/RQ{fragment_id}",
        f"Sx{fragment_id}(Inactive)=Py{fragment_id}*Qz{fragment_id}-"
        f"Pz{fragment_id}*Qy{fragment_id}",
        f"Sy{fragment_id}(Inactive)=Pz{fragment_id}*Qx{fragment_id}-"
        f"Px{fragment_id}*Qz{fragment_id}",
        f"Sz{fragment_id}(Inactive)=Px{fragment_id}*Qy{fragment_id}-"
        f"Py{fragment_id}*Qx{fragment_id}",
    ]


def _gaussian_quaternion_lines(frag_id: str, ref_id: str) -> list[str]:
    pair = f"{frag_id}{ref_id}"
    rows = []
    for left_axis, left_prefix in (("1", "P"), ("2", "Q"), ("3", "S")):
        for right_axis, right_prefix in (("1", "P"), ("2", "Q"), ("3", "S")):
            rows.append(
                f"R{left_axis}{right_axis}{pair}(Inactive)="
                f"{left_prefix}x{frag_id}*{right_prefix}x{ref_id}+"
                f"{left_prefix}y{frag_id}*{right_prefix}y{ref_id}+"
                f"{left_prefix}z{frag_id}*{right_prefix}z{ref_id}"
            )
    return [
        *rows,
        f"Kw{pair}(Inactive)=0.5*SQRT(R11{pair}+R22{pair}+R33{pair}+1)",
        f"Kx{pair}(Inactive)=(R23{pair}-R32{pair})/(4*Kw{pair})",
        f"Ky{pair}(Inactive)=(R31{pair}-R13{pair})/(4*Kw{pair})",
        f"Kz{pair}(Inactive)=(R12{pair}-R21{pair})/(4*Kw{pair})",
        f"Kn{pair}(Inactive)=SQRT(Kx{pair}**2+Ky{pair}**2+Kz{pair}**2+1.0D-24)",
        f"Th{pair}(Inactive)=2*ATAN(Kn{pair}/Kw{pair})",
        f"Ex{pair}(Inactive)=Th{pair}*Kx{pair}/Kn{pair}",
        f"Ey{pair}(Inactive)=Th{pair}*Ky{pair}/Kn{pair}",
        f"Ez{pair}(Inactive)=Th{pair}*Kz{pair}/Kn{pair}",
    ]


def _primitive_by_id(
    primitives: tuple[GICPrimitive, ...],
    identifier: str,
) -> GICPrimitive:
    for primitive in primitives:
        if primitive.identifier == identifier:
            return primitive
    raise GICForgeContractError(f"unknown primitive id in frozen GIC: {identifier}")


def _atom_interval(atoms: tuple[int, ...]) -> str:
    sorted_atoms = sorted(atoms)
    if not sorted_atoms:
        return ""
    ranges: list[str] = []
    start = previous = sorted_atoms[0]
    for atom in sorted_atoms[1:]:
        if atom == previous + 1:
            previous = atom
            continue
        ranges.append(f"{start}-{previous}" if start != previous else str(start))
        start = previous = atom
    ranges.append(f"{start}-{previous}" if start != previous else str(start))
    return ",".join(ranges)


def _frozen_gic_line(gic: FrozenGIC) -> str:
    coefficients = gic.coefficients or ((gic.primitive_id, 1.0),)
    coeffs = ",".join(
        f"{primitive_id}:{coefficient:.12f}"
        for primitive_id, coefficient in coefficients
    )
    return (
        f"{gic.identifier} NAME={gic.name} FAMILY={gic.family} IRREP={gic.irrep} "
        f"COEFFS={coeffs}"
    )


def _sycart_components(vector: tuple[float, ...]) -> str:
    parts = []
    axes = ("X", "Y", "Z")
    for idx, value in enumerate(vector):
        if abs(value) <= 1.0e-12:
            continue
        atom = idx // 3 + 1
        axis = axes[idx % 3]
        parts.append(f"{atom}:{axis}={value:.12g}")
    return "COMPONENTS=" + (";".join(parts) if parts else "NONE")


def _symmetry_mode(definition: GICDefinition) -> str:
    diagnostics = definition.symmetry_diagnostics
    if (
        definition.symmetrize
        and diagnostics is not None
        and diagnostics.method == POINT_GROUP_PROJECTOR_METHOD
        and diagnostics.groups
    ):
        return "POINT_GROUP_PROJECTOR"
    if (
        definition.symmetrize
        and diagnostics is not None
        and diagnostics.method == LOCAL_SYMMETRIZATION_METHOD
        and diagnostics.groups
    ):
        return "LOCAL_BLOCK_C1" if definition.point_group.upper() == "C1" else "LOCAL_BLOCK"
    if definition.point_group.upper() == "C1":
        return "IDENTITY_C1" if definition.symmetrize else "UNSYMMETRIZED_C1"
    return "SYMMETRIZED" if definition.symmetrize else "UNSYMMETRIZED"


def _bool_text(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def _unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm <= RANK_TOLERANCE:
        raise FloatingPointError("zero-length vector")
    return vector / norm


def _dot_unit(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.dot(_unit(left), _unit(right)))


def _orthogonal_frame(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    trial = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(axis, trial))) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])
    e1 = _unit(np.cross(axis, trial))
    e2 = _unit(np.cross(axis, e1))
    return e1, e2


def _periodic_delta(value_plus: float, value_minus: float) -> float:
    delta = float(value_plus - value_minus)
    while delta > np.pi:
        delta -= 2.0 * np.pi
    while delta < -np.pi:
        delta += 2.0 * np.pi
    return delta
