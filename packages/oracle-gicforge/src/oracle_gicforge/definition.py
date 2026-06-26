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


GIC_BACKEND = "oracle-native-primitive.v1"
SYCART_BACKEND = "oracle-native-cartesian-nullspace.v1"
RANK_TOLERANCE = 1.0e-7
FINITE_DIFFERENCE_STEP = 1.0e-5
LINEAR_ANGLE_DEGREES = 175.0


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


@dataclass(frozen=True)
class FrozenGIC:
    identifier: str
    name: str
    family: str
    irrep: str
    primitive_id: str
    gaussian_expression: str


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
    target = Path(path)
    validate_gicforge_prerequisites(target)
    lines = read_sectioned_lines(target)
    geometry = read_enriched_xyz(target)
    coords = np.asarray(geometry.coordinates_angstrom, dtype=float)
    point_group = _point_group(lines)
    if symmetrize and point_group.upper() != "C1":
        raise GICForgeContractError(
            "non-C1 GIC symmetrization is not available in the native backend yet"
        )

    bonds = _topology_bonds(lines, natoms=geometry.natoms)
    candidates = _primitive_candidates(
        bonds,
        coords=coords,
        natoms=geometry.natoms,
        fragment_records=_fragment_records(target),
    )
    target_rank = _vibrational_rank(coords)
    selected, rank = _select_ranked_primitives(
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
        )
        for idx, primitive in enumerate(selected, start=1)
    )
    return GICDefinition(
        backend=GIC_BACKEND,
        point_group=point_group,
        symmetrize=bool(symmetrize),
        target_rank=target_rank,
        rank=rank,
        candidate_count=len(candidates),
        reference_coordinates_angstrom=tuple(
            tuple(float(value) for value in row) for row in coords
        ),
        primitives=tuple(selected),
        gics=gics,
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
        f"SYMMETRIZE {_bool_text(definition.symmetrize)}",
        f"SYMMETRY_MODE {_symmetry_mode(definition)}",
        f"TARGET_RANK {definition.target_rank}",
        f"RANK {definition.rank}",
        f"CANDIDATE_COUNT {definition.candidate_count}",
        f"PRIMITIVE_COUNT {len(definition.primitives)}",
        f"GIC_COUNT {len(definition.gics)}",
        "RANK_METHOD finite_difference_b_matrix_greedy",
        f"FINITE_DIFFERENCE_STEP_ANGSTROM {FINITE_DIFFERENCE_STEP:.12g}",
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


def _primitive_candidates(
    bonds: tuple[tuple[int, int], ...],
    *,
    coords: np.ndarray,
    natoms: int,
    fragment_records: tuple[object, ...] = (),
) -> tuple[GICPrimitive, ...]:
    adjacency = _adjacency(bonds, natoms=natoms)
    counters: dict[str, int] = {
        "STRETCH": 0,
        "BEND": 0,
        "LINEAR_BEND": 0,
        "TORSION": 0,
        "OUT_OF_PLANE": 0,
        "FRAG_DISTANCE": 0,
        "FRAG_CENTER_ATOM_DISTANCE": 0,
        "FRAG_TRANSLATION": 0,
        "FRAG_ORIENTATION": 0,
    }
    candidates: list[GICPrimitive] = []

    for i, j in bonds:
        candidates.append(_make_primitive("STRETCH", "R", (i, j), counters))

    candidates.extend(
        _fragment_primitive_candidates(fragment_records, coords=coords, counters=counters)
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
                candidates.append(_make_primitive("BEND", "A", (i, center, k), counters))

    seen_torsions: set[tuple[int, int, int, int]] = set()
    for j, k in bonds:
        for i in sorted(adjacency[j] - {k}):
            for l in sorted(adjacency[k] - {j}):
                torsion = (i, j, k, l)
                canonical = min(torsion, tuple(reversed(torsion)))
                if canonical in seen_torsions:
                    continue
                seen_torsions.add(canonical)
                candidates.append(_make_primitive("TORSION", "D", canonical, counters))

    for center in range(1, natoms + 1):
        neighbors = sorted(adjacency[center])
        if len(neighbors) < 3:
            continue
        for n1, n2, n3 in combinations(neighbors, 3):
            candidates.append(
                _make_primitive("OUT_OF_PLANE", "U", (center, n1, n2, n3), counters)
            )

    return tuple(candidates)


def _make_primitive(
    family: str,
    function: str,
    atoms: tuple[int, ...],
    counters: dict[str, int],
    *,
    mode: int = 0,
    ref_atoms: tuple[int, ...] = (),
    refs: tuple[str, ...] = (),
) -> GICPrimitive:
    counters[family] += 1
    prefix = {
        "STRETCH": "Str",
        "BEND": "Bend",
        "LINEAR_BEND": "LinB",
        "TORSION": "Tors",
        "OUT_OF_PLANE": "OuPl",
        "FRAG_DISTANCE": "FCDi",
        "FRAG_CENTER_ATOM_DISTANCE": "FCAt",
        "FRAG_TRANSLATION": "FTrn",
        "FRAG_ORIENTATION": "FRot",
    }[family]
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
    )


def _select_ranked_primitives(
    candidates: tuple[GICPrimitive, ...],
    coords: np.ndarray,
    *,
    target_rank: int,
    rank_tolerance: float,
) -> tuple[tuple[GICPrimitive, ...], int]:
    if target_rank == 0:
        return (), 0
    selected: list[GICPrimitive] = []
    rows: list[np.ndarray] = []
    rank = 0
    for primitive in candidates:
        row = _finite_difference_b_row(primitive, coords)
        norm = float(np.linalg.norm(row))
        if not np.isfinite(norm) or norm <= rank_tolerance:
            continue
        normalized = row / norm
        trial = np.vstack([*rows, normalized]) if rows else normalized.reshape(1, -1)
        trial_rank = int(np.linalg.matrix_rank(trial, tol=rank_tolerance))
        if trial_rank <= rank:
            continue
        selected.append(primitive)
        rows.append(normalized)
        rank = trial_rank
        if rank == target_rank:
            break
    return tuple(selected), rank


def _finite_difference_b_row(primitive: GICPrimitive, coords: np.ndarray) -> np.ndarray:
    flat = np.asarray(coords, dtype=float).reshape(-1)
    base = _primitive_value(primitive, coords)
    row = np.zeros_like(flat)
    for idx in range(flat.size):
        plus = flat.copy()
        minus = flat.copy()
        plus[idx] += FINITE_DIFFERENCE_STEP
        minus[idx] -= FINITE_DIFFERENCE_STEP
        try:
            value_plus = _primitive_value(primitive, plus.reshape(coords.shape))
            value_minus = _primitive_value(primitive, minus.reshape(coords.shape))
        except FloatingPointError:
            row[idx] = np.nan
            continue
        if primitive.function == "D":
            delta = _periodic_delta(value_plus, value_minus)
        else:
            delta = value_plus - value_minus
        if not np.isfinite(delta) or not np.isfinite(base):
            row[idx] = np.nan
        else:
            row[idx] = delta / (2.0 * FINITE_DIFFERENCE_STEP)
    return row


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
    if primitive.function == "FTRANS":
        return _fragment_translation_value(
            coords,
            primitive.atoms,
            primitive.ref_atoms,
            mode=primitive.mode,
        )
    if primitive.function == "FROT":
        return _fragment_rotation_value(
            coords,
            primitive.atoms,
            primitive.ref_atoms,
            mode=primitive.mode,
        )
    raise GICForgeContractError(f"unsupported primitive function: {primitive.function}")


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
) -> float:
    frame_ref = _fragment_frame(coords, ref_atoms)
    frame_frag = _fragment_frame(coords, atoms)
    rotation = frame_ref.T @ frame_frag
    return float(_quaternion_vector(rotation)[mode])


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
                    )
                )
    return tuple(candidates)


def _fragment_records(path: Path) -> tuple[object, ...]:
    try:
        from oracle_fragments import read_fragment_records
    except ImportError:
        return ()
    return tuple(read_fragment_records(Path(path)))


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


def _fragment_frame(coords: np.ndarray, atoms: tuple[int, ...]) -> np.ndarray:
    if _fragment_frame_rank(coords, atoms) < 2:
        raise FloatingPointError("fragment orientation is underdefined")
    p_atom, q_atom = _fragment_frame_anchor_atoms(atoms, coords=coords)
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
        qx = (rotation[2, 1] - rotation[1, 2]) * scale
        qy = (rotation[0, 2] - rotation[2, 0]) * scale
        qz = (rotation[1, 0] - rotation[0, 1]) * scale
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
        qx = (rotation[2, 1] - rotation[1, 2]) * scale
        qy = (rotation[0, 2] - rotation[2, 0]) * scale
        qz = (rotation[1, 0] - rotation[0, 1]) * scale
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
    return (
        f"{primitive.identifier} NAME={primitive.name} FAMILY={primitive.family} "
        f"FUNCTION={primitive.function} ATOMS={atoms}{ref_atoms}{refs}{mode} "
        f"GAUSSIAN={primitive.gaussian_expression()}"
    )


def _gaussian_gic_block_lines(definition: GICDefinition) -> list[str]:
    coords = np.asarray(definition.reference_coordinates_angstrom, dtype=float)
    lines: list[str] = []
    fragment_atoms = _gaussian_fragment_atoms(definition.primitives)
    if fragment_atoms:
        for fragment_id, atoms in sorted(fragment_atoms.items()):
            lines.append(f"{fragment_id}=Fragment({_atom_interval(atoms)})")
        for fragment_id in sorted(fragment_atoms):
            lines.extend(_gaussian_center_lines(fragment_id))
        frot_pairs = _gaussian_frot_pairs(definition.primitives)
        frame_fragments = sorted({fragment for pair in frot_pairs for fragment in pair})
        for fragment_id in frame_fragments:
            atoms = fragment_atoms[fragment_id]
            lines.extend(_gaussian_frame_lines(fragment_id, atoms, coords=coords))
        for frag_id, ref_id in sorted(frot_pairs):
            lines.extend(_gaussian_quaternion_lines(frag_id, ref_id))

    for gic in definition.gics:
        primitive = _primitive_by_id(definition.primitives, gic.primitive_id)
        expression = _gaussian_expression_for_primitive(primitive)
        if expression:
            lines.append(f"{gic.identifier} = {expression}")
    return lines


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
    if primitive.function == "FTRANS":
        frag_id, ref_id = primitive.refs
        axis = ("x", "y", "z")[primitive.mode]
        return f"C{axis}{frag_id}-C{axis}{ref_id}"
    if primitive.function == "FROT":
        frag_id, ref_id = primitive.refs
        axis = ("x", "y", "z")[primitive.mode]
        return f"K{axis}{frag_id}{ref_id}"
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


def _gaussian_frot_pairs(primitives: tuple[GICPrimitive, ...]) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for primitive in primitives:
        if primitive.function == "FROT" and len(primitive.refs) >= 2:
            pairs.add((primitive.refs[0], primitive.refs[1]))
    return pairs


def _gaussian_center_lines(fragment_id: str) -> list[str]:
    return [
        f"Cx{fragment_id}(Inactive)=XCntr({fragment_id})",
        f"Cy{fragment_id}(Inactive)=YCntr({fragment_id})",
        f"Cz{fragment_id}(Inactive)=ZCntr({fragment_id})",
    ]


def _gaussian_frame_lines(
    fragment_id: str,
    atoms: tuple[int, ...],
    *,
    coords: np.ndarray,
) -> list[str]:
    p_atom, q_atom = _fragment_frame_anchor_atoms(atoms, coords=coords)
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
    return (
        f"{gic.identifier} NAME={gic.name} FAMILY={gic.family} IRREP={gic.irrep} "
        f"COEFFS={gic.primitive_id}:1.000000000000"
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
