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

    def gaussian_expression(self) -> str:
        atoms = ",".join(str(atom) for atom in self.atoms)
        if self.function == "L":
            return f"L({atoms},0,{self.mode})"
        return f"{self.function}({atoms})"


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
    if definition.gics:
        lines.extend(f"{gic.identifier} = {gic.gaussian_expression}" for gic in definition.gics)
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
) -> tuple[GICPrimitive, ...]:
    adjacency = _adjacency(bonds, natoms=natoms)
    counters: dict[str, int] = {
        "STRETCH": 0,
        "BEND": 0,
        "LINEAR_BEND": 0,
        "TORSION": 0,
        "OUT_OF_PLANE": 0,
    }
    candidates: list[GICPrimitive] = []

    for i, j in bonds:
        candidates.append(_make_primitive("STRETCH", "R", (i, j), counters))

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
) -> GICPrimitive:
    counters[family] += 1
    prefix = {
        "STRETCH": "Str",
        "BEND": "Bend",
        "LINEAR_BEND": "LinB",
        "TORSION": "Tors",
        "OUT_OF_PLANE": "OuPl",
    }[family]
    serial = sum(counters.values())
    return GICPrimitive(
        identifier=f"P{serial:03d}",
        name=f"{prefix}{counters[family]:04d}",
        family=family,
        function=function,
        atoms=tuple(int(atom) for atom in atoms),
        mode=int(mode),
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
    return (
        f"{primitive.identifier} NAME={primitive.name} FAMILY={primitive.family} "
        f"FUNCTION={primitive.function} ATOMS={atoms}{mode} "
        f"GAUSSIAN={primitive.gaussian_expression()}"
    )


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
