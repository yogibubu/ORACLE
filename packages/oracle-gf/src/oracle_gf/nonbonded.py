from __future__ import annotations

from pathlib import Path
from collections import deque

import numpy as np

from oracle_chem.topology.vdw_radii import uff_well_depth_kcal, vdw_radius
from oracle_core import read_sectioned_lines, section_content


BOHR_TO_ANGSTROM = 0.52917721092
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
KCAL_MOL_TO_HARTREE = 1.0 / 627.5094740631


def synthon_charges_from_xyzin(path: Path, natoms: int) -> tuple[np.ndarray, str]:
    """Read the canonical GF charge vector from #SYNTHONS.

    The values are CM5 when Gaussian CM5 was imported during preprocessing;
    otherwise they are the ORACLE synthon charges from electronegativity.
    """
    content = section_content(read_sectioned_lines(Path(path)), "SYNTHONS")
    charges = np.full(int(natoms), np.nan, dtype=float)
    source = "Synthons electronegativity model"
    for raw in content:
        text = raw.strip()
        if not text:
            continue
        upper = text.upper()
        if upper.startswith("CHARGE_SOURCE"):
            source = text.split(None, 1)[1].strip() if len(text.split(None, 1)) == 2 else source
            continue
        if upper.startswith(("SCHEMA ", "INDEXING ", "BOND_ORDER_SOURCE", "COLUMNS ")):
            continue
        parts = text.split()
        if len(parts) < 5:
            continue
        try:
            idx = int(parts[0]) - 1
            charge = float(parts[3])
        except ValueError:
            continue
        if 0 <= idx < len(charges):
            charges[idx] = charge
    if np.any(~np.isfinite(charges)):
        raise ValueError("#SYNTHONS does not contain a complete charge column")
    return charges, source


def nonbonded_cartesian_hessian_correction(
    coordinates_bohr: np.ndarray,
    atomic_numbers: np.ndarray,
    topology_bonds: tuple[tuple[int, int], ...],
    *,
    charges: np.ndarray | None = None,
    electrostatic: bool = False,
    uff_vdw: bool = False,
    one_four_scale: float = 0.5,
) -> np.ndarray:
    """Build the Cartesian Hessian contribution to subtract before GF."""
    coords = np.asarray(coordinates_bohr, dtype=float)
    atomic_numbers = np.asarray(atomic_numbers, dtype=int)
    natoms = len(atomic_numbers)
    correction = np.zeros((3 * natoms, 3 * natoms), dtype=float)
    graph = _zero_based_graph(natoms, topology_bonds)
    distances = _topological_distances(graph, natoms)
    scale_14 = float(one_four_scale)
    if scale_14 < 0.0:
        raise ValueError("1-4 non-bonded scale factor must be non-negative")

    if electrostatic:
        if charges is None:
            raise ValueError("Electrostatic Hessian correction requires atomic charges")
        correction += _electrostatic_hessian(
            coords,
            np.asarray(charges, dtype=float),
            distances,
            one_four_scale=scale_14,
        )
    if uff_vdw:
        correction += _uff_vdw_hessian(coords, atomic_numbers, distances, one_four_scale=scale_14)
    return 0.5 * (correction + correction.T)


def _electrostatic_hessian(
    coordinates_bohr: np.ndarray,
    charges: np.ndarray,
    distances: dict[tuple[int, int], int],
    *,
    one_four_scale: float,
) -> np.ndarray:
    natoms = len(charges)
    hessian = np.zeros((3 * natoms, 3 * natoms), dtype=float)
    for i in range(natoms):
        for j in range(i + 1, natoms):
            pair_scale = _nonbonded_pair_scale(distances.get((i, j), 999999), one_four_scale)
            if pair_scale == 0.0:
                continue
            qiqj = float(charges[i]) * float(charges[j])
            if abs(qiqj) <= 1.0e-16:
                continue
            rvec = coordinates_bohr[i] - coordinates_bohr[j]
            r = float(np.linalg.norm(rvec))
            if r <= 1.0e-12:
                continue
            first = pair_scale * (-qiqj / r**2)
            second = pair_scale * (2.0 * qiqj / r**3)
            _accumulate_pair_hessian(hessian, i, j, rvec, first, second)
    return hessian


def _uff_vdw_hessian(
    coordinates_bohr: np.ndarray,
    atomic_numbers: np.ndarray,
    distances: dict[tuple[int, int], int],
    *,
    one_four_scale: float,
) -> np.ndarray:
    natoms = len(atomic_numbers)
    hessian = np.zeros((3 * natoms, 3 * natoms), dtype=float)
    for i in range(natoms):
        for j in range(i + 1, natoms):
            pair_scale = _nonbonded_pair_scale(distances.get((i, j), 999999), one_four_scale)
            if pair_scale == 0.0:
                continue
            xij_bohr, dij_hartree = _uff_pair_parameters(int(atomic_numbers[i]), int(atomic_numbers[j]))
            rvec = coordinates_bohr[i] - coordinates_bohr[j]
            r = float(np.linalg.norm(rvec))
            if r <= 1.0e-12:
                continue
            x6 = xij_bohr**6
            x12 = x6 * x6
            r7 = r**7
            r8 = r**8
            r13 = r**13
            r14 = r**14
            first = pair_scale * dij_hartree * (-12.0 * x12 / r13 + 12.0 * x6 / r7)
            second = pair_scale * dij_hartree * (156.0 * x12 / r14 - 84.0 * x6 / r8)
            _accumulate_pair_hessian(hessian, i, j, rvec, first, second)
    return hessian


def _uff_pair_parameters(zi: int, zj: int) -> tuple[float, float]:
    ri = vdw_radius(zi, scheme="uff")
    rj = vdw_radius(zj, scheme="uff")
    di = uff_well_depth_kcal(zi)
    dj = uff_well_depth_kcal(zj)
    if ri is None or rj is None or di is None or dj is None:
        raise ValueError(f"UFF vdW parameters are not available for pair Z={zi}, Z={zj}")
    xij_angstrom = 2.0 * float(np.sqrt(float(ri) * float(rj)))
    dij_hartree = float(np.sqrt(float(di) * float(dj))) * KCAL_MOL_TO_HARTREE
    return xij_angstrom * ANGSTROM_TO_BOHR, dij_hartree


def _nonbonded_pair_scale(topological_distance: int, one_four_scale: float) -> float:
    if topological_distance <= 2:
        return 0.0
    if topological_distance == 3:
        return float(one_four_scale)
    return 1.0


def _accumulate_pair_hessian(
    hessian: np.ndarray,
    i: int,
    j: int,
    rvec: np.ndarray,
    first_derivative: float,
    second_derivative: float,
) -> None:
    r = float(np.linalg.norm(rvec))
    unit = rvec / r
    block = (second_derivative - first_derivative / r) * np.outer(unit, unit)
    block += (first_derivative / r) * np.eye(3)
    si = slice(3 * i, 3 * i + 3)
    sj = slice(3 * j, 3 * j + 3)
    hessian[si, si] += block
    hessian[sj, sj] += block
    hessian[si, sj] -= block
    hessian[sj, si] -= block


def _zero_based_graph(natoms: int, bonds: tuple[tuple[int, int], ...]) -> dict[int, set[int]]:
    graph = {idx: set() for idx in range(natoms)}
    for left, right in bonds:
        i = int(left) - 1
        j = int(right) - 1
        if 0 <= i < natoms and 0 <= j < natoms and i != j:
            graph[i].add(j)
            graph[j].add(i)
    return graph


def _topological_distances(graph: dict[int, set[int]], natoms: int) -> dict[tuple[int, int], int]:
    distances: dict[tuple[int, int], int] = {}
    for start in range(natoms):
        seen = {start: 0}
        queue: deque[int] = deque([start])
        while queue:
            current = queue.popleft()
            for neighbor in graph.get(current, set()):
                if neighbor in seen:
                    continue
                seen[neighbor] = seen[current] + 1
                queue.append(neighbor)
        for end, distance in seen.items():
            if start < end:
                distances[(start, end)] = distance
    return distances
