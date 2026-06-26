from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from merlino_fit.survibfit.geometry import angle
from merlino_fit.survibfit.pipeline import b_matrix_analytic
from merlino_fit.survibfit.primitives import Primitive
from merlino_fit.topology.covalent_radii import covalent_radius
from merlino_fit.topology.pipeline import build_topology_objects
from merlino_fit.topology.ringset import RingSet
from topology.elements import atomic_number

from .model import (
    GICDefinition,
    _definition_coordinate_kind_counts,
    _gicforge_cartesian_from_gauin,
    _primitive_signature,
    define_gics_from_cartesian,
)


LINEAR_THRESHOLD_RAD = np.deg2rad(170.0)
COLLAPSED_BOND_THRESHOLD_ANGSTROM = 0.2


@dataclass(frozen=True)
class GICForgePythonCoordinate:
    name: str
    block: str
    terms: tuple[tuple[float, Primitive], ...]
    type_index: int = 0

    @property
    def dominant_kind(self) -> str:
        return self.terms[0][1].kind if self.terms else "unknown"


@dataclass(frozen=True)
class GICForgePythonModel:
    atom_symbols: tuple[str, ...]
    atomic_numbers: tuple[int, ...]
    coordinates_angstrom: tuple[tuple[float, float, float], ...]
    primitive_candidates: tuple[GICForgePythonCoordinate, ...]
    coordinates: tuple[GICForgePythonCoordinate, ...]
    target_rank: int
    primitive_fallback: bool
    diagnostics: dict[str, object]

    def to_definition(self, *, workdir: Path | None = None) -> GICDefinition:
        primitive_basis = _primitive_basis(self.coordinates)
        row_index = {primitive: index for index, primitive in enumerate(primitive_basis)}
        u_matrix = np.zeros((len(primitive_basis), len(self.coordinates)), dtype=float)
        labels: list[str] = []
        names: list[str] = []
        for column, coordinate in enumerate(self.coordinates):
            names.append(coordinate.name)
            for coefficient, primitive in coordinate.terms:
                u_matrix[row_index[primitive], column] += float(coefficient)
            labels.append(
                f"GIC{column + 1:03d} GICForgePython {coordinate.name} "
                f"irrep=UNK {_format_terms(coordinate.terms)}"
            )
        definition = GICDefinition(
            atom_symbols=self.atom_symbols,
            atomic_numbers=self.atomic_numbers,
            reference_coordinates_angstrom=self.coordinates_angstrom,
            primitives=primitive_basis,
            u_matrix=u_matrix,
            labels=tuple(labels),
            names=tuple(names),
            irreps=tuple("UNK" for _ in names),
            point_group="UNKNOWN",
            symmetrized=False,
            symmetry_source="none",
            gaussian_input="\n".join(_format_readgic(name, coord.terms) for name, coord in zip(names, self.coordinates))
            + "\n",
            source="gicforge-python",
            generation_workdir=str(workdir) if workdir is not None else None,
            provenance={
                "backend": "gicforge-python",
                "target_vibrational_rank": str(self.target_rank),
                "primitive_fallback": str(self.primitive_fallback).lower(),
                "svd_local": str(self.diagnostics.get("svd_local", False)).lower(),
            },
        )
        if workdir is not None:
            Path(workdir).mkdir(parents=True, exist_ok=True)
            (Path(workdir) / "gicforge_python_diagnostics.json").write_text(
                json.dumps(self.diagnostics, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        return definition


def build_gicforge_python_model(
    atom_symbols: Iterable[str],
    coordinates_angstrom: np.ndarray,
    *,
    impdih: bool = False,
    onedih: bool = True,
    svd_local: bool = False,
    max_linear_angle_pairs_per_center: int = 3,
    linear_threshold: float = LINEAR_THRESHOLD_RAD,
    primitive_fallback: bool = True,
) -> GICForgePythonModel:
    atoms = tuple(str(atom).strip() for atom in atom_symbols)
    coords = np.asarray(coordinates_angstrom, dtype=float)
    if coords.shape != (len(atoms), 3):
        raise ValueError(f"Expected coordinate shape ({len(atoms)}, 3), got {coords.shape}")
    atomic_numbers = tuple(atomic_number(atom) for atom in atoms)
    _cg, graph, ringset, _synthons, _aromaticity = build_topology_objects(coords, np.asarray(atomic_numbers))
    _validate_no_spurious_hh_contacts(coords, atomic_numbers, graph.bonds)
    if _remove_collapsed_bonds(graph, coords):
        ringset = RingSet(graph, coords=coords)
    primitive_blocks = _fortran_like_primitive_blocks(
        graph,
        coords,
        atomic_numbers=atomic_numbers,
        ringset=ringset,
        impdih=impdih,
        onedih=onedih,
        svd_local=svd_local,
        max_linear_angle_pairs_per_center=max_linear_angle_pairs_per_center,
        linear_threshold=linear_threshold,
    )
    primitive_candidates = tuple(coord for block in primitive_blocks for coord in block)
    target = _target_rank(coords, graph)
    candidates = primitive_candidates
    if primitive_fallback and len(candidates) < target:
        primitive_blocks = _primitive_fallback_blocks(
            graph,
            coords,
            atomic_numbers=atomic_numbers,
            ringset=ringset,
            impdih=impdih,
            linear_threshold=linear_threshold,
        )
        candidates = tuple(coord for block in primitive_blocks for coord in block)
    if not primitive_fallback and len(candidates) < target:
        raise ValueError(f"GICForge Python candidates below vibrational rank ({len(candidates)} < {target})")
    if len(candidates) < target:
        raise ValueError(f"Primitive candidates below vibrational rank ({len(candidates)} < {target})")
    coordinates = _prune_type_local(candidates, coords, target_rank=target, block_pruning=svd_local)
    if primitive_fallback and (
        len(coordinates) < target or (svd_local and _coordinate_b_rank(coordinates, coords) < target)
    ):
        primitive_blocks = _primitive_fallback_blocks(
            graph,
            coords,
            atomic_numbers=atomic_numbers,
            ringset=ringset,
            impdih=impdih,
            linear_threshold=linear_threshold,
        )
        candidates = tuple(coord for block in primitive_blocks for coord in block)
        if len(candidates) < target:
            raise ValueError(f"Primitive candidates below vibrational rank ({len(candidates)} < {target})")
        coordinates = _prune_type_local(candidates, coords, target_rank=target, block_pruning=svd_local)
        if len(coordinates) < target or _coordinate_b_rank(coordinates, coords) < target:
            raise ValueError(
                f"GICForge Python primitive fallback did not reach vibrational rank "
                f"({len(coordinates)} coordinates for target {target})"
            )
    diagnostics = _python_model_diagnostics(
        candidates,
        coordinates,
        coords,
        target_rank=target,
        svd_local=svd_local,
        onedih=onedih,
        max_linear_angle_pairs_per_center=max_linear_angle_pairs_per_center,
    )
    return GICForgePythonModel(
        atom_symbols=atoms,
        atomic_numbers=atomic_numbers,
        coordinates_angstrom=tuple(tuple(float(value) for value in row) for row in coords),
        primitive_candidates=primitive_candidates,
        coordinates=coordinates,
        target_rank=target,
        primitive_fallback=primitive_fallback,
        diagnostics=diagnostics,
    )


def _validate_no_spurious_hh_contacts(
    coords: np.ndarray,
    atomic_numbers: tuple[int, ...],
    bonds: Iterable[tuple[int, int]],
) -> None:
    bonded = {tuple(sorted((int(i), int(j)))) for i, j in bonds}
    contacts: list[str] = []
    for i, zi in enumerate(atomic_numbers):
        if zi != 1:
            continue
        ri = covalent_radius(zi)
        if ri is None:
            continue
        for j in range(i + 1, len(atomic_numbers)):
            if atomic_numbers[j] != 1 or (i, j) in bonded:
                continue
            rj = covalent_radius(atomic_numbers[j])
            if rj is None:
                continue
            distance = float(np.linalg.norm(coords[i] - coords[j]))
            if distance <= 1.25 * (float(ri) + float(rj)):
                contacts.append(f"{i + 1}-{j + 1} ({distance:.3f} A)")
    if contacts:
        preview = ", ".join(contacts[:8])
        extra = f"; {len(contacts) - 8} additional H-H contacts" if len(contacts) > 8 else ""
        raise ValueError(f"GICForge Python input topology validation failed: spurious nonbonded H-H contact {preview}{extra}")


def compare_gicforge_python_to_fortran(
    atom_symbols: Iterable[str],
    coordinates_angstrom: np.ndarray,
    *,
    workdir: Path,
    executable: Path | None = None,
    impdih: bool = False,
    onedih: bool = True,
    svd_local: bool = False,
) -> dict[str, object]:
    workdir = Path(workdir)
    fortran_dir = workdir / "fortran"
    python_model = build_gicforge_python_model(
        atom_symbols,
        coordinates_angstrom,
        impdih=impdih,
        onedih=onedih,
        svd_local=svd_local,
    )
    extra_keywords = []
    if impdih:
        extra_keywords.append("IMPDIH")
    if not onedih:
        extra_keywords.append("NOONEDIH")
    if svd_local:
        extra_keywords.append("LOCSVD")
    fortran_definition = define_gics_from_cartesian(
        tuple(atom_symbols),
        np.asarray(coordinates_angstrom, dtype=float),
        workdir=fortran_dir,
        executable=executable,
        symmetrize=False,
        extra_keywords=tuple(extra_keywords),
    )
    raw_coords = _gicforge_cartesian_from_gauin(fortran_dir / "gauin", len(fortran_definition.atom_symbols))
    fortran_signatures = tuple(_primitive_signature(primitive) for primitive in fortran_definition.primitives)
    python_candidates: list[tuple[str, GICDefinition]] = [
        ("input", python_model.to_definition(workdir=workdir / "python-input"))
    ]
    try:
        python_candidates.append(
            (
                "fortran-frame",
                build_gicforge_python_model(
                    atom_symbols,
                    raw_coords,
                    impdih=impdih,
                    onedih=onedih,
                    svd_local=svd_local,
                ).to_definition(workdir=workdir / "python-fortran-frame"),
            )
        )
    except Exception:
        pass
    selected_frame, python_definition = python_candidates[0]
    for candidate_frame, candidate_definition in python_candidates:
        candidate_signatures = tuple(_primitive_signature(primitive) for primitive in candidate_definition.primitives)
        if (
            len(candidate_definition.names) == len(fortran_definition.names)
            and _definition_coordinate_kind_counts(candidate_definition)
            == _definition_coordinate_kind_counts(fortran_definition)
            and candidate_signatures == fortran_signatures
        ):
            selected_frame = candidate_frame
            python_definition = candidate_definition
            break
    python_signatures = tuple(_primitive_signature(primitive) for primitive in python_definition.primitives)
    same_ordered_primitives = fortran_signatures == python_signatures
    b_max_abs_diff = None
    if same_ordered_primitives and fortran_definition.u_matrix.shape == python_definition.u_matrix.shape:
        fortran_b = fortran_definition.u_matrix.T @ b_matrix_analytic(fortran_definition.primitives, raw_coords)
        python_b = python_definition.u_matrix.T @ b_matrix_analytic(
            python_definition.primitives,
            raw_coords,
        )
        b_max_abs_diff = float(np.max(np.abs(python_b - fortran_b))) if python_b.size else 0.0
    return {
        "passed": (
            len(python_definition.names) == len(fortran_definition.names)
            and _definition_coordinate_kind_counts(python_definition)
            == _definition_coordinate_kind_counts(fortran_definition)
            and same_ordered_primitives
            and (b_max_abs_diff is None or b_max_abs_diff <= 1.0e-7)
        ),
        "target_rank": python_model.target_rank,
        "python_gic_count": len(python_definition.names),
        "fortran_gic_count": len(fortran_definition.names),
        "python_kind_counts": _definition_coordinate_kind_counts(python_definition),
        "fortran_kind_counts": _definition_coordinate_kind_counts(fortran_definition),
        "python_primitive_count": len(python_definition.primitives),
        "fortran_primitive_count": len(fortran_definition.primitives),
        "same_ordered_primitives": same_ordered_primitives,
        "b_max_abs_diff": b_max_abs_diff,
        "python_names": list(python_definition.names),
        "fortran_names": list(fortran_definition.names),
        "python_comparison_frame": selected_frame,
        "python_workdir": str(workdir / f"python-{selected_frame}"),
        "fortran_workdir": str(fortran_dir),
    }


def _fortran_like_primitive_blocks(
    graph,
    coords: np.ndarray,
    *,
    atomic_numbers: tuple[int, ...],
    ringset,
    impdih: bool,
    onedih: bool,
    svd_local: bool,
    max_linear_angle_pairs_per_center: int,
    linear_threshold: float,
):
    bonds: list[GICForgePythonCoordinate] = []
    bends: list[GICForgePythonCoordinate] = []
    linears: list[GICForgePythonCoordinate] = []
    torsions: list[GICForgePythonCoordinate] = []
    oops: list[GICForgePythonCoordinate] = []
    neighbors = [sorted(graph.adjacency[index]) for index in range(graph.natoms)]
    effective_atomic_numbers = _effective_atomic_numbers(graph, coords, atomic_numbers, neighbors)
    selected_rings = _minimum_cycle_basis(
        graph,
        ringset,
        effective_atomic_numbers=effective_atomic_numbers,
        neighbors=neighbors,
    )
    atom_ring = _atom_ring_map_from_rings(selected_rings, graph.natoms)
    ring_counts = _atom_selected_ring_counts(selected_rings, graph.natoms)
    ring_bonds = {
        tuple(sorted((ring[i], ring[(i + 1) % len(ring)])))
        for ring in selected_rings
        for i in range(len(ring))
    }
    bridge_bonds = _bridge_bonds(selected_rings)

    for center in range(graph.natoms):
        neigh = neighbors[center]
        for first in neigh:
            if first < center:
                continue
            bonds.append(_primitive_coordinate("Stre", len(bonds) + 1, Primitive("bond", (center, first))))
        if svd_local and len(neigh) > 1:
            exo_primitives, exo_linears = _exocyclic_angle_primitives(
                center,
                neigh,
                selected_rings=selected_rings,
                coords=coords,
                linear_threshold=linear_threshold,
            )
            if exo_primitives:
                bends.extend(
                    _svd_local_coordinates(
                        exo_primitives,
                        coords=coords,
                        prefix="XAng",
                        start=len(bends) + 1,
                        kind_type_index=0,
                    )
                )
            exo_linears = exo_linears[:max(0, max_linear_angle_pairs_per_center)] if len(neigh) == 2 else []
            for primitive in exo_linears:
                linears.append(_primitive_coordinate("LAng", len(linears) + 1, primitive))
                linears.append(
                    _primitive_coordinate(
                        "LAng",
                        len(linears) + 1,
                        Primitive("linear_bend", primitive.atoms, mode=-2),
                    )
                )
        elif len(neigh) == 3:
            bends.extend(
                _c2v3_angle_coordinates(
                    center,
                    neigh,
                    atomic_numbers=atomic_numbers,
                    effective_atomic_numbers=effective_atomic_numbers,
                    neighbors=neighbors,
                    atom_ring=atom_ring,
                    start=len(bends) + 1,
                )
            )
        elif len(neigh) > 1:
            if len(neigh) == 2 and atom_ring[center] != 0:
                continue
            if len(neigh) == 4 and not _has_linear_pair(center, neigh, coords, linear_threshold):
                bends.extend(
                    _four_atom_angle_coordinates(
                        center,
                        neigh,
                        effective_atomic_numbers=effective_atomic_numbers,
                        atom_ring=atom_ring,
                        coords=coords,
                        start=len(bends) + 1,
                    )
                )
                continue
            if len(neigh) > 4:
                high_bends, high_linears = _high_coord_angle_coordinates(
                    center,
                    neigh,
                    coords=coords,
                    linear_threshold=linear_threshold,
                    angle_start=len(bends) + 1,
                    linear_start=len(linears) + 1,
                )
                bends.extend(high_bends)
                linears.extend(high_linears)
                continue
            for ib, first_angle in enumerate(neigh[:-1]):
                for second_angle in neigh[ib + 1 :]:
                    value = angle(first_angle, center, second_angle, coords)
                    left, right = sorted((first_angle, second_angle))
                    if value < linear_threshold:
                        bends.append(
                            _primitive_coordinate("Bend", len(bends) + 1, Primitive("angle", (left, center, right)))
                        )
                    else:
                        linears.append(
                            _primitive_coordinate(
                                "LAng",
                                len(linears) + 1,
                                Primitive("linear_bend", (left, center, right), mode=-1),
                            )
                        )
                        linears.append(
                            _primitive_coordinate(
                                "LAng",
                                len(linears) + 1,
                                Primitive("linear_bend", (left, center, right), mode=-2),
                            )
                        )

    for ring in selected_rings:
        if _all_atoms_in_three_selected_rings(ring, selected_rings):
            continue
        if svd_local:
            bends.extend(
                _cyclic_svd_coordinates(
                    ring,
                    valence_angle=True,
                    coords=coords,
                    prefix="RDef",
                    start=len(bends) + 1,
                )
            )
        else:
            bends.extend(_cyclic_coordinates(ring, valence_angle=True, prefix="RDef", start=len(bends) + 1))

    for bond in bonds:
        _coef, primitive = bond.terms[0]
        center, right = primitive.atoms
        if tuple(sorted((center, right))) in bridge_bonds:
            if ring_counts[center] < 3 and ring_counts[right] < 3:
                butterfly = _butterfly_coordinate(
                    center,
                    right,
                    neighbors=neighbors,
                    atom_ring=atom_ring,
                    selected_rings=selected_rings,
                    coords=coords,
                    linear_threshold=linear_threshold,
                    index=len(torsions) + 1,
                )
                if butterfly is not None:
                    torsions.append(butterfly)
            continue
        if tuple(sorted((center, right))) in ring_bonds:
            continue
        if len(neighbors[center]) == 1 or len(neighbors[right]) == 1:
            continue
        torsion_factory = _onedih_torsion_coordinate if onedih else _torsion_coordinate
        torsion = torsion_factory(
            center,
            right,
            neighbors=neighbors,
            atomic_numbers=atomic_numbers,
            effective_atomic_numbers=effective_atomic_numbers,
            atom_ring=atom_ring,
            ring_counts=ring_counts,
            coords=coords,
            linear_threshold=linear_threshold,
            index=len(torsions) + 1,
        )
        if torsion is not None:
            torsions.append(torsion)

    for ring in selected_rings:
        if _all_atoms_in_three_selected_rings(ring, selected_rings):
            continue
        if svd_local:
            torsions.extend(
                _cyclic_svd_coordinates(
                    ring,
                    valence_angle=False,
                    coords=coords,
                    prefix="RPck",
                    start=len(torsions) + 1,
                )
            )
        else:
            torsions.extend(_cyclic_coordinates(ring, valence_angle=False, prefix="RPck", start=len(torsions) + 1))

    oop_prefix = "ImpD" if impdih else "OuPl"
    for center in range(graph.natoms):
        neigh = neighbors[center]
        if len(neigh) != 3:
            continue
        first, second, third = neigh
        if all(atom_ring[atom] != 0 for atom in (center, first, second, third)):
            continue
        if impdih:
            primitive = Primitive("dihedral", (first, center, third, second))
        else:
            primitive = Primitive("out_of_plane", (center, first, second, third))
        oops.append(_primitive_coordinate(oop_prefix, len(oops) + 1, primitive))

    return bonds, bends, linears, torsions, oops


def _primitive_fallback_blocks(
    graph,
    coords: np.ndarray,
    *,
    atomic_numbers: tuple[int, ...],
    ringset,
    impdih: bool,
    linear_threshold: float,
):
    bonds: list[GICForgePythonCoordinate] = []
    bends: list[GICForgePythonCoordinate] = []
    linears: list[GICForgePythonCoordinate] = []
    torsions: list[GICForgePythonCoordinate] = []
    oops: list[GICForgePythonCoordinate] = []
    neighbors = [sorted(graph.adjacency[index]) for index in range(graph.natoms)]
    effective_atomic_numbers = _effective_atomic_numbers(graph, coords, atomic_numbers, neighbors)
    selected_rings = _minimum_cycle_basis(
        graph,
        ringset,
        effective_atomic_numbers=effective_atomic_numbers,
        neighbors=neighbors,
    )
    atom_ring = _atom_ring_map_from_rings(selected_rings, graph.natoms)

    for center in range(graph.natoms):
        neigh = neighbors[center]
        for first in neigh:
            if first < center:
                continue
            bonds.append(_primitive_coordinate("Stre", len(bonds) + 1, Primitive("bond", (center, first))))
        for ib, first_angle in enumerate(neigh[:-1]):
            for second_angle in neigh[ib + 1 :]:
                left, right = sorted((first_angle, second_angle))
                primitive = Primitive("angle", (left, center, right))
                if angle(first_angle, center, second_angle, coords) < linear_threshold:
                    bends.append(_primitive_coordinate("Bend", len(bends) + 1, primitive))
                else:
                    linears.append(
                        _primitive_coordinate(
                            "LAng",
                            len(linears) + 1,
                            Primitive("linear_bend", primitive.atoms, mode=-1),
                        )
                    )
                    linears.append(
                        _primitive_coordinate(
                            "LAng",
                            len(linears) + 1,
                            Primitive("linear_bend", primitive.atoms, mode=-2),
                        )
                    )

    for bond in bonds:
        _coef, primitive = bond.terms[0]
        center, right = primitive.atoms
        if len(neighbors[center]) == 1 or len(neighbors[right]) == 1:
            continue
        for left in neighbors[center]:
            if left == right:
                continue
            if angle(left, center, right, coords) > linear_threshold:
                continue
            for far in neighbors[right]:
                if far == center or far == left:
                    continue
                if angle(center, right, far, coords) > linear_threshold:
                    continue
                torsions.append(
                    _primitive_coordinate(
                        "Dihe",
                        len(torsions) + 1,
                        Primitive("dihedral", (left, center, right, far)),
                    )
                )

    oop_prefix = "ImpD" if impdih else "OuPl"
    for center in range(graph.natoms):
        neigh = neighbors[center]
        if len(neigh) != 3:
            continue
        first, second, third = neigh
        if all(atom_ring[atom] != 0 for atom in (center, first, second, third)):
            continue
        if impdih:
            primitive = Primitive("dihedral", (first, center, third, second))
        else:
            primitive = Primitive("out_of_plane", (center, first, second, third))
        oops.append(_primitive_coordinate(oop_prefix, len(oops) + 1, primitive))

    return bonds, bends, linears, torsions, oops


def _has_linear_pair(center: int, neigh: list[int], coords: np.ndarray, linear_threshold: float) -> bool:
    for ib, first in enumerate(neigh[:-1]):
        for second in neigh[ib + 1 :]:
            if angle(first, center, second, coords) >= linear_threshold:
                return True
    return False


def _remove_collapsed_bonds(graph, coords: np.ndarray) -> bool:
    removed = False
    kept_bonds = []
    for first, second in graph.bonds:
        distance = float(np.linalg.norm(coords[first] - coords[second]))
        if distance < COLLAPSED_BOND_THRESHOLD_ANGSTROM:
            graph.adjacency[first].discard(second)
            graph.adjacency[second].discard(first)
            removed = True
        else:
            kept_bonds.append((first, second))
    if removed:
        graph.bonds = kept_bonds
    return removed


def _minimum_cycle_basis(
    graph,
    ringset,
    *,
    effective_atomic_numbers: tuple[float, ...],
    neighbors: list[list[int]],
) -> list[tuple[int, ...]]:
    if ringset is None:
        return []
    edges = sorted(tuple(sorted(edge)) for edge in graph.bonds)
    edge_index = {edge: index for index, edge in enumerate(edges)}
    target = len(edges) - graph.natoms + len(_connected_components(graph))
    if target <= 0:
        return []
    basis: list[set[int]] = []
    pivots: list[int] = []
    selected: list[tuple[int, ...]] = []

    def reduce_vector(vector: set[int]) -> set[int]:
        reduced = set(vector)
        for pivot, row in zip(pivots, basis):
            if pivot in reduced:
                reduced ^= row
        return reduced

    def add_vector(vector: set[int]) -> bool:
        reduced = reduce_vector(vector)
        if not reduced:
            return False
        pivot = min(reduced)
        index = 0
        while index < len(pivots) and pivots[index] < pivot:
            index += 1
        pivots.insert(index, pivot)
        basis.insert(index, reduced)
        return True

    for ring in sorted(ringset.rings, key=lambda item: (len(item.atoms), item.index)):
        atoms = tuple(int(atom) for atom in ring.atoms)
        vector = {
            edge_index[tuple(sorted((atoms[index], atoms[(index + 1) % len(atoms)])))]
            for index in range(len(atoms))
        }
        if add_vector(vector):
            selected.append(
                _orient_ring_for_gicforge(
                    atoms,
                    effective_atomic_numbers=effective_atomic_numbers,
                    neighbors=neighbors,
                )
            )
        if len(selected) >= target:
            break
    return selected


def _orient_ring_for_gicforge(
    atoms: tuple[int, ...],
    *,
    effective_atomic_numbers: tuple[float, ...],
    neighbors: list[list[int]],
) -> tuple[int, ...]:
    if len(atoms) <= 3:
        return atoms

    traversal = _connected_ring_traversal(atoms, neighbors)
    candidates = []
    for seq in (traversal, tuple(reversed(traversal))):
        for index in range(len(seq)):
            candidates.append(seq[index:] + seq[:index])
    return max(
        candidates,
        key=lambda seq: _ring_priority_key(
            seq,
            effective_atomic_numbers=effective_atomic_numbers,
            neighbors=neighbors,
        ),
    )


def _connected_ring_traversal(atoms: tuple[int, ...], neighbors: list[list[int]]) -> tuple[int, ...]:
    remaining = list(atoms)
    first = min(remaining)
    start = remaining.index(first)
    remaining[0], remaining[start] = remaining[start], remaining[0]
    for index in range(len(remaining) - 1):
        current = remaining[index]
        next_index = index + 1
        for candidate_index in range(index + 1, len(remaining)):
            if remaining[candidate_index] in neighbors[current]:
                next_index = candidate_index
                break
        remaining[index + 1], remaining[next_index] = remaining[next_index], remaining[index + 1]
    return tuple(remaining)


def _ring_priority_key(
    ring: tuple[int, ...],
    *,
    effective_atomic_numbers: tuple[float, ...],
    neighbors: list[list[int]],
) -> tuple[tuple[float, int, tuple[float, ...], int], ...]:
    return tuple(
        (
            round(effective_atomic_numbers[atom], 12),
            len(neighbors[atom]),
            tuple(_exocyclic_neighbor_priorities(atom, ring, effective_atomic_numbers, neighbors)),
            -atom,
        )
        for atom in ring
    )


def _exocyclic_neighbor_priorities(
    atom: int,
    ring: tuple[int, ...],
    effective_atomic_numbers: tuple[float, ...],
    neighbors: list[list[int]],
) -> tuple[float, ...]:
    ring_atoms = set(ring)
    return tuple(
        sorted(
            (round(effective_atomic_numbers[neighbor], 12) for neighbor in neighbors[atom] if neighbor not in ring_atoms),
            reverse=True,
        )
    )


def _effective_atomic_numbers(
    graph,
    coords: np.ndarray,
    atomic_numbers: tuple[int, ...],
    neighbors: list[list[int]],
) -> tuple[float, ...]:
    synthons: list[float] = []
    for center in range(graph.natoms):
        neigh = neighbors[center]
        if not neigh:
            synthons.append(0.0)
            continue
        z_center = atomic_numbers[center]
        nval = z_center
        if z_center > 2:
            nval -= 2
        if z_center > 10:
            nval -= 8
        nmax = 8 - nval
        neff = nmax - len(neigh) + 1
        if nval == 1:
            t0 = 180.0
        else:
            t0_values = {1: 109.47, 2: 120.0, 3: 180.0}
            t0 = t0_values.get(neff, 109.47)

        delocalization = 1.0
        coordination = 0.0
        rigidity = 0.0
        angle_count = 0
        for pos, neighbor in enumerate(neigh):
            distance = float(np.linalg.norm(coords[center] - coords[neighbor]))
            radius_center = covalent_radius(z_center) or 0.0
            radius_neighbor = covalent_radius(atomic_numbers[neighbor]) or 0.0
            bond_order = float(np.exp(((radius_center + radius_neighbor) - distance) / 0.3))
            delocalization *= bond_order
            coordination += float(atomic_numbers[neighbor]) * bond_order
            if len(neigh) == 1:
                continue
            for other in neigh[pos:]:
                if other == neighbor:
                    continue
                angle_count += 1
                rigidity += abs(np.sin(angle(neighbor, center, other, coords) - np.deg2rad(t0)))
        if angle_count == 0:
            angle_count = 1
        synthons.append(coordination / len(neigh) + delocalization / len(neigh) + rigidity / angle_count)

    synmax = max(synthons) if synthons else 0.0
    denominator = synmax + 0.1
    return tuple(float(z) - 0.495 + synthon / denominator for z, synthon in zip(atomic_numbers, synthons))


def _all_atoms_in_three_selected_rings(ring: tuple[int, ...], rings: list[tuple[int, ...]]) -> bool:
    counts = {atom: 0 for atom in ring}
    for candidate in rings:
        for atom in ring:
            if atom in candidate:
                counts[atom] += 1
    return all(count >= 3 for count in counts.values())


def _bridge_bonds(rings: list[tuple[int, ...]]) -> set[tuple[int, int]]:
    counts: dict[tuple[int, int], int] = {}
    for ring in rings:
        for index, atom in enumerate(ring):
            edge = tuple(sorted((atom, ring[(index + 1) % len(ring)])))
            counts[edge] = counts.get(edge, 0) + 1
    return {edge for edge, count in counts.items() if count >= 2}


def _butterfly_coordinate(
    center: int,
    right: int,
    *,
    neighbors: list[list[int]],
    atom_ring: list[int],
    selected_rings: list[tuple[int, ...]],
    coords: np.ndarray,
    linear_threshold: float,
    index: int,
) -> GICForgePythonCoordinate | None:
    terms: list[tuple[float, Primitive]] = []
    for left in neighbors[center]:
        if left == right:
            continue
        if angle(left, center, right, coords) > linear_threshold:
            continue
        if atom_ring[left] == 0:
            continue
        for far in neighbors[right]:
            if far == center or far == left:
                continue
            if angle(center, right, far, coords) > linear_threshold:
                continue
            if atom_ring[far] == 0:
                continue
            if _atoms_share_selected_ring(left, far, selected_rings):
                continue
            coefficient = 1.0 if not terms else -1.0
            terms.append((coefficient, Primitive("dihedral", (left, center, right, far))))
    if not terms:
        return None
    norm = np.sqrt(float(len(terms)))
    return GICForgePythonCoordinate(
        name=f"BtFl{index:04d}",
        block="BtFl",
        type_index=2,
        terms=tuple((coefficient / norm, primitive) for coefficient, primitive in terms),
    )


def _atoms_share_selected_ring(first: int, second: int, rings: list[tuple[int, ...]]) -> bool:
    for ring in rings:
        if first in ring and second in ring:
            return True
    return False


def _exocyclic_angle_primitives(
    center: int,
    neigh: list[int],
    *,
    selected_rings: list[tuple[int, ...]],
    coords: np.ndarray,
    linear_threshold: float,
) -> tuple[list[Primitive], list[Primitive]]:
    angles: list[Primitive] = []
    linears: list[Primitive] = []
    for index, first in enumerate(neigh[:-1]):
        for second in neigh[index + 1 :]:
            left, right = sorted((first, second))
            if _is_endocyclic_angle(left, center, right, selected_rings):
                continue
            primitive = Primitive("angle", (left, center, right))
            if angle(left, center, right, coords) < linear_threshold:
                angles.append(primitive)
            else:
                linears.append(Primitive("linear_bend", (left, center, right), mode=-1))
    return angles, linears


def _is_endocyclic_angle(left: int, center: int, right: int, rings: list[tuple[int, ...]]) -> bool:
    for ring in rings:
        size = len(ring)
        for index, atom in enumerate(ring):
            if atom != center:
                continue
            if {left, right} == {ring[(index - 1) % size], ring[(index + 1) % size]}:
                return True
    return False


def _torsion_coordinate(
    center: int,
    right: int,
    *,
    neighbors: list[list[int]],
    atomic_numbers: tuple[int, ...],
    effective_atomic_numbers: tuple[float, ...],
    atom_ring: list[int],
    ring_counts: list[int],
    coords: np.ndarray,
    linear_threshold: float,
    index: int,
) -> GICForgePythonCoordinate | None:
    candidates: list[tuple[int, int]] = []
    for left in neighbors[center]:
        if left == right:
            continue
        if angle(left, center, right, coords) > linear_threshold:
            continue
        for far in neighbors[right]:
            if far == center or far == left:
                continue
            if angle(center, right, far, coords) > linear_threshold:
                continue
            if ring_counts[left] >= 2 or ring_counts[far] >= 2:
                continue
            candidates.append((left, far))
    if not candidates:
        return None
    coefficient = 1.0 / np.sqrt(float(len(candidates)))
    return GICForgePythonCoordinate(
        name=f"Tors{index:04d}",
        block="Tors",
        type_index=-1,
        terms=tuple((coefficient, Primitive("dihedral", (left, center, right, far))) for left, far in candidates),
    )


def _onedih_torsion_coordinate(
    center: int,
    right: int,
    *,
    neighbors: list[list[int]],
    atomic_numbers: tuple[int, ...],
    effective_atomic_numbers: tuple[float, ...],
    atom_ring: list[int],
    ring_counts: list[int],
    coords: np.ndarray,
    linear_threshold: float,
    index: int,
) -> GICForgePythonCoordinate | None:
    candidates: list[tuple[int, int]] = []
    for left in neighbors[center]:
        if left == right:
            continue
        if angle(left, center, right, coords) > linear_threshold:
            continue
        if ring_counts[left] >= 2:
            continue
        for far in neighbors[right]:
            if far == center or far == left:
                continue
            if angle(center, right, far, coords) > linear_threshold:
                continue
            if ring_counts[far] >= 2:
                continue
            candidates.append((left, far))
    if not candidates:
        return None

    threshold = 5.0e-4
    selected_left, selected_far = max(
        candidates,
        key=lambda pair: (
            round(effective_atomic_numbers[pair[0]] / threshold) * threshold,
            round(effective_atomic_numbers[pair[1]] / threshold) * threshold,
            len(neighbors[center]),
            len(neighbors[right]),
            -pair[0],
            -pair[1],
        ),
    )
    orbit = [
        (left, far)
        for left, far in candidates
        if atom_ring[left] == atom_ring[selected_left]
        and atomic_numbers[left] == atomic_numbers[selected_left]
        and len(neighbors[left]) == len(neighbors[selected_left])
        and atom_ring[far] == atom_ring[selected_far]
        and atomic_numbers[far] == atomic_numbers[selected_far]
        and len(neighbors[far]) == len(neighbors[selected_far])
    ]
    if not orbit:
        orbit = [(selected_left, selected_far)]
    coefficient = 1.0 / np.sqrt(float(len(orbit)))
    return GICForgePythonCoordinate(
        name=f"Tors{index:04d}",
        block="Tors",
        type_index=-1,
        terms=tuple((coefficient, Primitive("dihedral", (left, center, right, far))) for left, far in orbit),
    )


def _cyclic_coordinates(
    ring: tuple[int, ...],
    *,
    valence_angle: bool,
    prefix: str,
    start: int,
) -> list[GICForgePythonCoordinate]:
    ncyc = len(ring)
    if ncyc == 3:
        return []
    istart = 2 if valence_angle else ncyc
    if valence_angle:
        if ncyc == 6:
            istart = 2
        elif ncyc == 7:
            istart = 4
        elif ncyc == 8:
            istart = 2
        else:
            istart = 3
    vnorm = np.sqrt(2.0 / float(ncyc))
    vnorm1 = np.sqrt(1.0 / float(ncyc))
    coordinates: list[GICForgePythonCoordinate] = []
    for ivar in range(1, ncyc - 2):
        even = ivar == 2 * (ivar // 2)
        terms = []
        for iterm in range(1, ncyc + 1):
            iang1 = _cyclic_index(iterm + istart - 1, ncyc)
            iang2 = _cyclic_index(iterm + istart, ncyc)
            iang3 = _cyclic_index(iterm + istart + 1, ncyc)
            iang4 = _cyclic_index(iterm + istart + 2, ncyc)
            ivar1 = ivar
            if ivar == 1:
                ivar1 = ivar + 1
            if ivar == 4:
                ivar1 = ivar - 1
            if ivar == 5:
                ivar1 = ivar - 1
            if ivar == 6:
                ivar1 = ivar - 2
            snum = float(2 * ivar1 * (iterm - 1))
            value = np.pi * snum / float(ncyc)
            if even:
                coefficient = vnorm * np.sin(value)
            elif ivar < ncyc - 3:
                coefficient = vnorm * np.cos(value)
            else:
                coefficient = vnorm1 * np.cos(float(iterm - 1) * np.pi)
            if abs(coefficient) < 1.0e-14:
                coefficient = 0.0
            if valence_angle:
                primitive = Primitive("angle", (ring[iang1], ring[iang2], ring[iang3]))
            else:
                primitive = Primitive("dihedral", (ring[iang1], ring[iang2], ring[iang3], ring[iang4]))
            terms.append((float(coefficient), primitive))
        coordinates.append(
            GICForgePythonCoordinate(
                name=f"{prefix}{start + len(coordinates):04d}",
                block=prefix,
                type_index=14 if valence_angle else 1,
                terms=tuple(terms),
            )
        )
    return coordinates


def _cyclic_svd_coordinates(
    ring: tuple[int, ...],
    *,
    valence_angle: bool,
    coords: np.ndarray,
    prefix: str,
    start: int,
) -> list[GICForgePythonCoordinate]:
    primitives = _cyclic_primitives_legacy_order(ring, valence_angle=valence_angle)
    if not primitives:
        return []
    reference = _cyclic_reference_coefficients(len(ring), valence_angle=valence_angle)
    primitive_b = b_matrix_analytic(tuple(primitives), coords)
    u_matrix, singular_values, _vh = np.linalg.svd(primitive_b, full_matrices=False)
    rank = min(_svd_rank(singular_values), max(0, len(ring) - 3))
    if rank == 0:
        return []
    coefficients = _align_svd_modes_to_reference(u_matrix[:, :rank], reference[:, :rank])
    coordinates: list[GICForgePythonCoordinate] = []
    for mode in range(rank):
        coeffs = coefficients[:, mode].astype(float)
        coeffs[np.abs(coeffs) < 1.0e-14] = 0.0
        terms = tuple(
            (float(coefficient), primitive)
            for coefficient, primitive in zip(coeffs, primitives)
            if abs(float(coefficient)) > 1.0e-12
        )
        if not terms:
            continue
        coordinates.append(
            GICForgePythonCoordinate(
                name=f"{prefix}{start + len(coordinates):04d}",
                block=prefix,
                type_index=14 if valence_angle else 1,
                terms=terms,
            )
        )
    return coordinates


def _cyclic_primitives(ring: tuple[int, ...], *, valence_angle: bool) -> list[Primitive]:
    ncyc = len(ring)
    if ncyc == 3:
        return []
    primitives: list[Primitive] = []
    for term in range(ncyc):
        if valence_angle:
            primitives.append(
                Primitive(
                    "angle",
                    (ring[(term - 1) % ncyc], ring[term], ring[(term + 1) % ncyc]),
                )
            )
        else:
            primitives.append(
                Primitive(
                    "dihedral",
                    (ring[(term - 1) % ncyc], ring[term], ring[(term + 1) % ncyc], ring[(term + 2) % ncyc]),
                )
            )
    return primitives


def _cyclic_primitives_legacy_order(ring: tuple[int, ...], *, valence_angle: bool) -> list[Primitive]:
    ncyc = len(ring)
    if ncyc == 3:
        return []
    istart = _cyclic_legacy_start(ncyc, valence_angle=valence_angle)
    primitives: list[Primitive] = []
    for iterm in range(1, ncyc + 1):
        iang1 = _cyclic_index(iterm + istart - 1, ncyc)
        iang2 = _cyclic_index(iterm + istart, ncyc)
        iang3 = _cyclic_index(iterm + istart + 1, ncyc)
        iang4 = _cyclic_index(iterm + istart + 2, ncyc)
        if valence_angle:
            primitives.append(Primitive("angle", (ring[iang1], ring[iang2], ring[iang3])))
        else:
            primitives.append(Primitive("dihedral", (ring[iang1], ring[iang2], ring[iang3], ring[iang4])))
    return primitives


def _cyclic_reference_coefficients(ncyc: int, *, valence_angle: bool) -> np.ndarray:
    reference = np.zeros((ncyc, max(0, ncyc - 3)), dtype=float)
    if ncyc <= 3:
        return reference
    vnorm = np.sqrt(2.0 / float(ncyc))
    vnorm1 = np.sqrt(1.0 / float(ncyc))
    for ivar in range(1, ncyc - 2):
        even = ivar == 2 * (ivar // 2)
        ivar1 = ivar
        if ivar == 1:
            ivar1 = ivar + 1
        if ivar == 4:
            ivar1 = ivar - 1
        if ivar == 5:
            ivar1 = ivar - 1
        if ivar == 6:
            ivar1 = ivar - 2
        for iterm in range(1, ncyc + 1):
            snum = float(2 * ivar1 * (iterm - 1))
            value = np.pi * snum / float(ncyc)
            if even:
                coefficient = vnorm * np.sin(value)
            elif ivar < ncyc - 3:
                coefficient = vnorm * np.cos(value)
            else:
                coefficient = vnorm1 * np.cos(float(iterm - 1) * np.pi)
            if abs(coefficient) < 1.0e-14:
                coefficient = 0.0
            reference[iterm - 1, ivar - 1] = float(coefficient)
    return reference


def _align_svd_modes_to_reference(u_matrix: np.ndarray, reference: np.ndarray) -> np.ndarray:
    aligned = np.zeros_like(reference)
    used: set[int] = set()
    for mode in range(reference.shape[1]):
        best_index = -1
        best_dot = 0.0
        best_score = -1.0
        for candidate in range(u_matrix.shape[1]):
            if candidate in used:
                continue
            dot = float(np.dot(reference[:, mode], u_matrix[:, candidate]))
            score = abs(dot)
            if score > best_score:
                best_index = candidate
                best_dot = dot
                best_score = score
        if best_index < 0:
            continue
        used.add(best_index)
        sign = -1.0 if best_dot < 0.0 else 1.0
        aligned[:, mode] = sign * u_matrix[:, best_index]
    return aligned


def _cyclic_legacy_start(ncyc: int, *, valence_angle: bool) -> int:
    if not valence_angle:
        return ncyc
    if ncyc == 6:
        return 2
    if ncyc == 7:
        return 4
    if ncyc == 8:
        return 2
    return 3


def _svd_local_coordinates(
    primitives: list[Primitive],
    *,
    coords: np.ndarray,
    prefix: str,
    start: int,
    kind_type_index: int,
    max_modes: int | None = None,
) -> list[GICForgePythonCoordinate]:
    if not primitives:
        return []
    primitive_b = b_matrix_analytic(tuple(primitives), coords)
    u_matrix, singular_values, _vh = np.linalg.svd(primitive_b, full_matrices=False)
    rank = _svd_rank(singular_values)
    if max_modes is not None:
        rank = min(rank, max_modes)
    coordinates: list[GICForgePythonCoordinate] = []
    for mode in range(rank):
        coeffs = u_matrix[:, mode].astype(float)
        coeffs = _canonical_svd_coefficients(coeffs)
        terms = tuple(
            (float(coefficient), primitive)
            for coefficient, primitive in zip(coeffs, primitives)
            if abs(float(coefficient)) > 1.0e-12
        )
        if not terms:
            continue
        coordinates.append(
            GICForgePythonCoordinate(
                name=f"{prefix}{start + len(coordinates):04d}",
                block=prefix,
                type_index=kind_type_index,
                terms=terms,
            )
        )
    return coordinates


def _svd_rank(singular_values: np.ndarray) -> int:
    if singular_values.size == 0:
        return 0
    tolerance = max(1.0e-10, 1.0e-8 * float(singular_values[0]))
    return int(np.sum(singular_values > tolerance))


def _canonical_svd_coefficients(coefficients: np.ndarray) -> np.ndarray:
    if coefficients.size == 0:
        return coefficients
    dominant = int(np.argmax(np.abs(coefficients)))
    if coefficients[dominant] < 0.0:
        coefficients = -coefficients
    coefficients[np.abs(coefficients) < 1.0e-14] = 0.0
    return coefficients


def _cyclic_index(index_1based: int, ncyc: int) -> int:
    while index_1based > ncyc:
        index_1based -= ncyc
    while index_1based <= 0:
        index_1based += ncyc
    return index_1based - 1


def _c2v3_angle_coordinates(
    center: int,
    neigh: list[int],
    *,
    atomic_numbers: tuple[int, ...],
    effective_atomic_numbers: tuple[float, ...],
    neighbors: list[list[int]],
    atom_ring: list[int],
    start: int,
) -> list[GICForgePythonCoordinate]:
    first, second, third = neigh
    threshold = 5.0e-4
    eq12 = abs(effective_atomic_numbers[first] - effective_atomic_numbers[second]) < threshold
    eq13 = abs(effective_atomic_numbers[first] - effective_atomic_numbers[third]) < threshold
    eq23 = abs(effective_atomic_numbers[second] - effective_atomic_numbers[third]) < threshold
    if not eq12 and not eq13 and not eq23:
        different = first
        if atomic_numbers[second] == 1:
            different = second
        elif atomic_numbers[third] == 1:
            different = third
        elif len(neighbors[second]) == 1:
            different = second
        elif len(neighbors[third]) == 1:
            different = third
    elif eq12 and eq13:
        different = first
    elif eq12 and not eq13:
        different = third
    elif eq13 and not eq12:
        different = second
    else:
        different = first

    if atom_ring[center] != 0:
        if atom_ring[first] != 0 and atom_ring[second] != 0:
            different = third
        if atom_ring[first] != 0 and atom_ring[third] != 0:
            different = second
        if atom_ring[second] != 0 and atom_ring[third] != 0:
            different = first

    if different == first:
        jat, kat, lat = first, second, third
    elif different == second:
        jat, kat, lat = second, first, third
    else:
        jat, kat, lat = third, second, first

    if atom_ring[center] != 0 and atom_ring[jat] != 0:
        return []

    coords: list[GICForgePythonCoordinate] = []
    if atom_ring[center] == 0:
        den = np.sqrt(6.0)
        coords.append(
            GICForgePythonCoordinate(
                name=f"SymD{start:04d}",
                block="SymD",
                type_index=1,
                terms=(
                    (2.0 / den, Primitive("angle", (kat, center, lat))),
                    (-1.0 / den, Primitive("angle", (jat, center, kat))),
                    (-1.0 / den, Primitive("angle", (jat, center, lat))),
                ),
            )
        )
        start += 1
    den = np.sqrt(2.0)
    coords.append(
        GICForgePythonCoordinate(
            name=f"Rock{start:04d}",
            block="Rock",
            type_index=2,
            terms=(
                (1.0 / den, Primitive("angle", (jat, center, kat))),
                (-1.0 / den, Primitive("angle", (jat, center, lat))),
            ),
        )
    )
    return coords


def _four_atom_angle_coordinates(
    center: int,
    neigh: list[int],
    *,
    effective_atomic_numbers: tuple[float, ...],
    atom_ring: list[int],
    coords: np.ndarray,
    start: int,
) -> list[GICForgePythonCoordinate]:
    first, second, third, fourth = _order_four_atom_neighbors(
        tuple(neigh),
        center=center,
        effective_atomic_numbers=effective_atomic_numbers,
        atom_ring=atom_ring,
    )
    frozen = {
        atom: atom_ring[atom] != 0 and atom_ring[center] != 0
        for atom in (first, second, third, fourth)
    }
    equal_count = _four_atom_equal_count((first, second, third, fourth), effective_atomic_numbers)
    pivot_count = sum(1 for atom in (first, second, third, fourth) if frozen[atom])
    if equal_count == 4:
        return _td_four_atom_coordinates(
            center,
            first,
            second,
            third,
            fourth,
            frozen=frozen,
            start=start,
        )
    if equal_count == 3 or pivot_count in {1, 3}:
        return _wxy3_coordinates(
            center,
            first,
            second,
            third,
            fourth,
            frozen=frozen,
            start=start,
        )
    return _w2xy2_coordinates(
        center,
        first,
        second,
        third,
        fourth,
        equal_count=equal_count,
        frozen=frozen,
        start=start,
    )


def _order_four_atom_neighbors(
    atoms: tuple[int, int, int, int],
    *,
    center: int,
    effective_atomic_numbers: tuple[float, ...],
    atom_ring: list[int],
) -> tuple[int, int, int, int]:
    jat, kat, lat, mat = atoms
    j1, k1, l1, m1 = atoms
    threshold = 5.0e-4

    def equivalent(first: int, second: int) -> bool:
        return abs(effective_atomic_numbers[first] - effective_atomic_numbers[second]) < threshold

    if equivalent(j1, k1):
        if equivalent(j1, l1):
            if not equivalent(j1, m1):
                pass
        elif equivalent(j1, m1):
            lat, mat = m1, l1
        elif equivalent(l1, m1):
            pass
    elif equivalent(j1, l1):
        jat, kat, lat, mat = j1, l1, k1, m1
        if equivalent(j1, m1):
            jat, kat, lat, mat = j1, l1, m1, k1
        elif equivalent(k1, m1):
            jat, kat, lat, mat = j1, l1, k1, m1
    elif equivalent(j1, m1):
        jat, kat, lat, mat = j1, m1, k1, l1
    elif equivalent(k1, l1):
        jat, kat, lat, mat = k1, l1, j1, m1
        if equivalent(l1, m1):
            jat, kat, lat, mat = k1, l1, m1, j1
    elif equivalent(k1, m1):
        jat, kat, lat, mat = k1, m1, j1, l1
    elif equivalent(l1, m1):
        jat, kat, lat, mat = l1, m1, j1, k1

    pivot_count = 0
    if atom_ring[center] != 0:
        pivot_count = sum(1 for atom in (jat, kat, lat, mat) if atom_ring[atom] != 0)
    if pivot_count == 1:
        if atom_ring[jat] != 0:
            return jat, kat, lat, mat
        if atom_ring[kat] != 0:
            return kat, jat, lat, mat
        if atom_ring[lat] != 0:
            return lat, kat, jat, mat
        if atom_ring[mat] != 0:
            return mat, kat, lat, jat
    if pivot_count == 3:
        if atom_ring[jat] == 0:
            ordered = (jat, kat, lat, mat)
        elif atom_ring[kat] == 0:
            ordered = (kat, jat, lat, mat)
        elif atom_ring[lat] == 0:
            ordered = (lat, kat, jat, mat)
        else:
            ordered = (mat, kat, lat, jat)
        jat, kat, lat, mat = ordered
        if abs(effective_atomic_numbers[kat] - effective_atomic_numbers[lat]) < threshold:
            if abs(effective_atomic_numbers[lat] - effective_atomic_numbers[mat]) >= threshold:
                kat, mat = mat, kat
        elif abs(effective_atomic_numbers[lat] - effective_atomic_numbers[mat]) >= threshold:
            kat, lat = lat, kat
    return jat, kat, lat, mat


def _four_atom_equal_count(atoms: tuple[int, int, int, int], effective_atomic_numbers: tuple[float, ...]) -> int:
    jat, kat, lat, mat = atoms
    threshold = 5.0e-4

    def equivalent(first: int, second: int) -> bool:
        return abs(effective_atomic_numbers[first] - effective_atomic_numbers[second]) < threshold

    neq = 0
    if equivalent(jat, kat):
        neq += 1
        if equivalent(jat, lat):
            neq += 2
            if equivalent(jat, mat):
                neq += 1
        elif equivalent(jat, mat):
            neq += 2
        elif equivalent(lat, mat):
            neq += 1
    elif equivalent(jat, lat):
        neq += 1
        if equivalent(jat, mat):
            neq += 2
        elif equivalent(kat, mat):
            neq += 1
    elif equivalent(jat, mat):
        neq += 1
        if equivalent(kat, lat):
            neq += 1
    elif equivalent(kat, lat):
        neq += 1
        if equivalent(lat, mat):
            neq += 2
    elif equivalent(kat, mat):
        neq += 1
    elif equivalent(lat, mat):
        neq += 1
    return neq


def _w2xy2_coordinates(
    center: int,
    jat: int,
    kat: int,
    lat: int,
    mat: int,
    *,
    equal_count: int,
    frozen: dict[int, bool],
    start: int,
) -> list[GICForgePythonCoordinate]:
    if all(frozen[atom] for atom in (jat, kat, lat, mat)):
        return []
    inot1, inot2 = jat, kat
    iyes1, iyes2 = lat, mat
    if frozen[jat] and frozen[lat]:
        inot1, inot2 = jat, lat
        iyes1, iyes2 = kat, mat
    elif frozen[jat] and frozen[mat]:
        inot1, inot2 = jat, mat
        iyes1, iyes2 = kat, lat
    elif frozen[kat] and frozen[lat]:
        inot1, inot2 = kat, lat
        iyes1, iyes2 = jat, mat
    elif frozen[kat] and frozen[mat]:
        inot1, inot2 = kat, mat
        iyes1, iyes2 = jat, lat
    elif frozen[lat] and frozen[mat]:
        inot1, inot2 = lat, mat
        iyes1, iyes2 = jat, kat

    den_sym = np.sqrt(6.0)
    den_rock = np.sqrt(2.0)
    coordinates = [
        GICForgePythonCoordinate(
            name=f"SymD{start:04d}",
            block="SymD",
            type_index=1,
            terms=(
                (2.0 / den_sym, Primitive("angle", (iyes1, center, iyes2))),
                (-1.0 / den_sym, Primitive("angle", (inot1, center, iyes1))),
                (-1.0 / den_sym, Primitive("angle", (inot1, center, iyes2))),
            ),
        ),
        GICForgePythonCoordinate(
            name=f"Rock{start + 1:04d}",
            block="Rock",
            type_index=2,
            terms=(
                (1.0 / den_rock, Primitive("angle", (inot1, center, iyes1))),
                (-1.0 / den_rock, Primitive("angle", (inot1, center, iyes2))),
            ),
        ),
        GICForgePythonCoordinate(
            name=f"SymD{start + 2:04d}",
            block="SymD",
            type_index=1,
            terms=(
                (2.0 / den_sym, Primitive("angle", (iyes1, center, iyes2))),
                (-1.0 / den_sym, Primitive("angle", (inot2, center, iyes1))),
                (-1.0 / den_sym, Primitive("angle", (inot2, center, iyes2))),
            ),
        ),
        GICForgePythonCoordinate(
            name=f"Rock{start + 3:04d}",
            block="Rock",
            type_index=2,
            terms=(
                (1.0 / den_rock, Primitive("angle", (inot2, center, iyes1))),
                (-1.0 / den_rock, Primitive("angle", (inot2, center, iyes2))),
            ),
        ),
    ]
    if not any(frozen.values()):
        coordinates.append(
            _primitive_coordinate("Bend" if equal_count != 2 else "Bend", start + 4, Primitive("angle", (inot1, center, inot2)))
        )
    return coordinates


def _wxy3_coordinates(
    center: int,
    jat: int,
    kat: int,
    lat: int,
    mat: int,
    *,
    frozen: dict[int, bool],
    start: int,
) -> list[GICForgePythonCoordinate]:
    if all(frozen[atom] for atom in (jat, kat, lat, mat)):
        return []
    den = np.sqrt(2.0)
    coordinates = [
        GICForgePythonCoordinate(
            name=f"Rock{start:04d}",
            block="Rock",
            type_index=2,
            terms=(
                (0.5, Primitive("angle", (jat, center, kat))),
                (-0.25, Primitive("angle", (jat, center, lat))),
                (-0.25, Primitive("angle", (jat, center, mat))),
            ),
        ),
        GICForgePythonCoordinate(
            name=f"Rock{start + 1:04d}",
            block="Rock",
            type_index=2,
            terms=(
                (1.0 / den, Primitive("angle", (jat, center, lat))),
                (-1.0 / den, Primitive("angle", (jat, center, mat))),
            ),
        ),
    ]
    if sum(1 for value in frozen.values() if value) == 3:
        return coordinates
    coordinates.append(_primitive_coordinate("Bend", start + 2, Primitive("angle", (lat, center, mat))))
    if sum(1 for value in frozen.values() if value) == 2:
        return coordinates
    coordinates.append(_primitive_coordinate("Bend", start + 3, Primitive("angle", (kat, center, lat))))
    coordinates.append(_primitive_coordinate("Bend", start + 4, Primitive("angle", (kat, center, mat))))
    return coordinates


def _td_four_atom_coordinates(
    center: int,
    jat: int,
    kat: int,
    lat: int,
    mat: int,
    *,
    frozen: dict[int, bool],
    start: int,
) -> list[GICForgePythonCoordinate]:
    den_ea = np.sqrt(12.0)
    den_eb = 2.0
    den_t2 = np.sqrt(2.0)
    return [
        GICForgePythonCoordinate(
            name=f"EEee{start:04d}",
            block="EEee",
            type_index=8,
            terms=(
                (2.0 / den_ea, Primitive("angle", (jat, center, kat))),
                (-1.0 / den_ea, Primitive("angle", (jat, center, lat))),
                (-1.0 / den_ea, Primitive("angle", (jat, center, mat))),
                (-1.0 / den_ea, Primitive("angle", (kat, center, lat))),
                (-1.0 / den_ea, Primitive("angle", (kat, center, mat))),
                (2.0 / den_ea, Primitive("angle", (lat, center, mat))),
            ),
        ),
        GICForgePythonCoordinate(
            name=f"EEee{start + 1:04d}",
            block="EEee",
            type_index=8,
            terms=(
                (1.0 / den_eb, Primitive("angle", (jat, center, lat))),
                (-1.0 / den_eb, Primitive("angle", (jat, center, mat))),
                (-1.0 / den_eb, Primitive("angle", (kat, center, lat))),
                (1.0 / den_eb, Primitive("angle", (kat, center, mat))),
            ),
        ),
        GICForgePythonCoordinate(
            name=f"T2xx{start + 2:04d}",
            block="T2xx",
            type_index=9,
            terms=(
                (1.0 / den_t2, Primitive("angle", (jat, center, lat))),
                (-1.0 / den_t2, Primitive("angle", (kat, center, mat))),
            ),
        ),
        GICForgePythonCoordinate(
            name=f"T2yy{start + 3:04d}",
            block="T2yy",
            type_index=10,
            terms=(
                (1.0 / den_t2, Primitive("angle", (kat, center, lat))),
                (-1.0 / den_t2, Primitive("angle", (jat, center, mat))),
            ),
        ),
        GICForgePythonCoordinate(
            name=f"T2zz{start + 4:04d}",
            block="T2zz",
            type_index=11,
            terms=(
                (1.0 / den_t2, Primitive("angle", (jat, center, kat))),
                (-1.0 / den_t2, Primitive("angle", (lat, center, mat))),
            ),
        ),
    ]


def _high_coord_angle_coordinates(
    center: int,
    neigh: list[int],
    *,
    coords: np.ndarray,
    linear_threshold: float,
    angle_start: int,
    linear_start: int,
) -> tuple[list[GICForgePythonCoordinate], list[GICForgePythonCoordinate]]:
    coordinates: list[GICForgePythonCoordinate] = []
    linears: list[GICForgePythonCoordinate] = []
    for ib, first in enumerate(neigh[:-1]):
        for second in neigh[ib + 1 :]:
            left, right = sorted((first, second))
            value = angle(left, center, right, coords)
            coordinates.append(
                _primitive_coordinate("HCAn", angle_start + len(coordinates), Primitive("angle", (left, center, right)))
            )
            if value >= linear_threshold:
                linears.append(
                    _primitive_coordinate(
                        "LAng",
                        linear_start + len(linears),
                        Primitive("linear_bend", (left, center, right), mode=-1),
                    )
                )
                linears.append(
                    _primitive_coordinate(
                        "LAng",
                        linear_start + len(linears),
                        Primitive("linear_bend", (left, center, right), mode=-2),
                    )
                )
    return coordinates, linears


def _atom_ring_map_from_rings(rings: list[tuple[int, ...]], natoms: int) -> list[int]:
    atom_ring = [0 for _ in range(natoms)]
    for index, ring in enumerate(rings, start=1):
        for atom in ring:
            atom_ring[int(atom)] = index
    return atom_ring


def _atom_selected_ring_counts(rings: list[tuple[int, ...]], natoms: int) -> list[int]:
    counts = [0 for _ in range(natoms)]
    for ring in rings:
        for atom in ring:
            counts[int(atom)] += 1
    return counts


def _primitive_coordinate(prefix: str, index: int, primitive: Primitive) -> GICForgePythonCoordinate:
    return GICForgePythonCoordinate(
        name=f"{prefix}{index:04d}",
        block=prefix,
        terms=((1.0, primitive),),
    )


def _python_model_diagnostics(
    candidates: tuple[GICForgePythonCoordinate, ...],
    coordinates: tuple[GICForgePythonCoordinate, ...],
    coords: np.ndarray,
    *,
    target_rank: int,
    svd_local: bool,
    onedih: bool,
    max_linear_angle_pairs_per_center: int,
) -> dict[str, object]:
    candidate_counts = _coordinate_block_counts(candidates)
    kept_counts = _coordinate_block_counts(coordinates)
    removed_counts = {
        block: int(candidate_counts.get(block, 0) - kept_counts.get(block, 0))
        for block in sorted(set(candidate_counts) | set(kept_counts))
    }
    candidate_rank = _coordinate_b_rank(candidates, coords)
    final_rank = _coordinate_b_rank(coordinates, coords)
    return {
        "backend": "gicforge-python",
        "svd_local": bool(svd_local),
        "onedih": bool(onedih),
        "target_rank": int(target_rank),
        "candidate_count": int(len(candidates)),
        "final_count": int(len(coordinates)),
        "candidate_rank": int(candidate_rank),
        "final_rank": int(final_rank),
        "rank_complete": bool(final_rank == target_rank),
        "count_complete": bool(len(coordinates) == target_rank),
        "max_linear_angle_pairs_per_center": int(max_linear_angle_pairs_per_center),
        "candidate_counts_by_block": candidate_counts,
        "final_counts_by_block": kept_counts,
        "removed_counts_by_block": removed_counts,
    }


def _coordinate_block_counts(coordinates: tuple[GICForgePythonCoordinate, ...]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for coordinate in coordinates:
        counts[coordinate.block] = counts.get(coordinate.block, 0) + 1
    return dict(sorted(counts.items()))


def _coordinate_b_rank(coordinates: tuple[GICForgePythonCoordinate, ...], coords: np.ndarray) -> int:
    if not coordinates:
        return 0
    primitive_basis = _primitive_basis(coordinates)
    row_index = {primitive: index for index, primitive in enumerate(primitive_basis)}
    primitive_b = b_matrix_analytic(primitive_basis, coords)
    b_rows = []
    for coordinate in coordinates:
        row = np.zeros(primitive_b.shape[1], dtype=float)
        for coefficient, primitive in coordinate.terms:
            row += coefficient * primitive_b[row_index[primitive]]
        b_rows.append(row)
    singular_values = np.linalg.svd(np.vstack(b_rows), compute_uv=False)
    return _svd_rank(singular_values)


def _prune_type_local(
    coordinates: tuple[GICForgePythonCoordinate, ...],
    coords: np.ndarray,
    *,
    target_rank: int,
    block_pruning: bool = False,
) -> tuple[GICForgePythonCoordinate, ...]:
    if block_pruning:
        return _prune_block_local(coordinates, coords, target_rank=target_rank)
    by_kind = {
        "bond": [coord for coord in coordinates if coord.dominant_kind == "bond"],
        "angle": [coord for coord in coordinates if coord.dominant_kind == "angle"],
        "linear_bend": [coord for coord in coordinates if coord.dominant_kind == "linear_bend"],
        "dihedral": [coord for coord in coordinates if coord.dominant_kind == "dihedral"],
        "out_of_plane": [coord for coord in coordinates if coord.dominant_kind == "out_of_plane"],
    }
    linear_bends = by_kind["linear_bend"]
    forced_linear_atoms: tuple[int, ...] | None = None
    if any(coord.block == "HCAn" for coord in by_kind["angle"]):
        linear_bends = _high_coord_linear_pruning_order(linear_bends)
        if linear_bends:
            forced_linear_atoms = linear_bends[0].terms[0][1].atoms
    ordered = by_kind["bond"] + by_kind["angle"] + linear_bends + by_kind["dihedral"] + by_kind["out_of_plane"]
    if len(ordered) <= target_rank:
        return tuple(ordered)
    primitive_basis = _primitive_basis(ordered)
    row_index = {primitive: index for index, primitive in enumerate(primitive_basis)}
    primitive_b = b_matrix_analytic(primitive_basis, coords)
    b_rows = []
    for coordinate in ordered:
        row = np.zeros(primitive_b.shape[1], dtype=float)
        for coefficient, primitive in coordinate.terms:
            row += coefficient * primitive_b[row_index[primitive]]
        b_rows.append(row)

    basis: list[np.ndarray] = []
    keep: list[GICForgePythonCoordinate] = []
    for index, coordinate in enumerate(ordered):
        if coordinate.dominant_kind == "bond":
            _seed_basis_row(b_rows[index], basis)
            keep.append(coordinate)
    for kind in ("angle", "linear_bend", "dihedral", "out_of_plane"):
        for index, coordinate in enumerate(ordered):
            if coordinate.dominant_kind != kind:
                continue
            if len(basis) >= target_rank or len(keep) >= target_rank:
                continue
            force_keep = (
                kind == "linear_bend"
                and forced_linear_atoms is not None
                and coordinate.terms[0][1].atoms == forced_linear_atoms
            )
            if _seed_basis_row(b_rows[index], basis):
                keep.append(coordinate)
            elif force_keep:
                keep.append(coordinate)
    return tuple(keep)


def _prune_block_local(
    coordinates: tuple[GICForgePythonCoordinate, ...],
    coords: np.ndarray,
    *,
    target_rank: int,
) -> tuple[GICForgePythonCoordinate, ...]:
    ordered = sorted(coordinates, key=lambda coord: (_block_pruning_priority(coord), coord.name))
    if len(ordered) <= target_rank:
        return tuple(ordered)
    primitive_basis = _primitive_basis(ordered)
    row_index = {primitive: index for index, primitive in enumerate(primitive_basis)}
    primitive_b = b_matrix_analytic(primitive_basis, coords)
    b_rows = []
    for coordinate in ordered:
        row = np.zeros(primitive_b.shape[1], dtype=float)
        for coefficient, primitive in coordinate.terms:
            row += coefficient * primitive_b[row_index[primitive]]
        b_rows.append(row)

    basis_rows: list[np.ndarray] = []
    keep: list[GICForgePythonCoordinate] = []
    for index, coordinate in enumerate(ordered):
        if len(basis_rows) >= target_rank or len(keep) >= target_rank:
            break
        if _seed_basis_row_by_svd(b_rows[index], basis_rows):
            keep.append(coordinate)
    return tuple(sorted(keep, key=lambda coord: coordinates.index(coord)))


def _block_pruning_priority(coordinate: GICForgePythonCoordinate) -> int:
    if coordinate.dominant_kind == "bond":
        return 0
    if coordinate.block == "XAng":
        return 1
    if coordinate.dominant_kind == "linear_bend":
        return 2
    if coordinate.block in {"Dihe", "Tors"}:
        return 3
    if coordinate.block == "BtFl":
        return 4
    if coordinate.block == "RDef":
        return 5
    if coordinate.block == "RPck":
        return 6
    if coordinate.dominant_kind == "out_of_plane":
        return 7
    return 8


def _high_coord_linear_pruning_order(
    coordinates: list[GICForgePythonCoordinate],
) -> list[GICForgePythonCoordinate]:
    groups: list[list[GICForgePythonCoordinate]] = []
    index_by_atoms: dict[tuple[int, ...], int] = {}
    for coordinate in coordinates:
        _coefficient, primitive = coordinate.terms[0]
        atoms = primitive.atoms
        if atoms not in index_by_atoms:
            index_by_atoms[atoms] = len(groups)
            groups.append([])
        groups[index_by_atoms[atoms]].append(coordinate)
    if not groups:
        return []

    ordered: list[GICForgePythonCoordinate] = []
    ordered.extend(groups[0])
    for group in groups[1:]:
        ordered.extend(coord for coord in group if coord.terms[0][1].mode == -2)
    for group in groups[1:]:
        ordered.extend(coord for coord in group if coord.terms[0][1].mode != -2)
    seen: set[GICForgePythonCoordinate] = set(ordered)
    ordered.extend(coord for coord in coordinates if coord not in seen)
    return ordered


def _seed_basis_row(row: np.ndarray, basis: list[np.ndarray], *, t_abs: float = 1.0e-10, t_rel: float = 1.0e-8) -> bool:
    candidate = np.asarray(row, dtype=float).copy()
    norm0 = float(np.linalg.norm(candidate))
    if norm0 <= t_abs:
        return False
    for existing in basis:
        candidate -= float(np.dot(candidate, existing)) * existing
    norm = float(np.linalg.norm(candidate))
    if norm > t_abs and norm > t_rel * norm0:
        basis.append(candidate / norm)
        return True
    return False


def _seed_basis_row_by_svd(row: np.ndarray, basis_rows: list[np.ndarray]) -> bool:
    candidate = np.asarray(row, dtype=float).copy()
    norm = float(np.linalg.norm(candidate))
    if norm <= 1.0e-10:
        return False
    if not basis_rows:
        basis_rows.append(candidate)
        return True
    old = np.vstack(basis_rows)
    new = np.vstack([old, candidate])
    old_rank = _svd_rank(np.linalg.svd(old, compute_uv=False))
    new_rank = _svd_rank(np.linalg.svd(new, compute_uv=False))
    if new_rank > old_rank:
        basis_rows.append(candidate)
        return True
    return False


def _target_rank(coords: np.ndarray, graph) -> int:
    components = _connected_components(graph)
    rank = 0
    for component in components:
        rank += 3 * len(component) - (5 if _is_linear(coords[list(component)]) else 6)
    return rank


def _connected_components(graph) -> list[tuple[int, ...]]:
    seen: set[int] = set()
    components: list[tuple[int, ...]] = []
    for start in range(graph.natoms):
        if start in seen:
            continue
        stack = [start]
        component = []
        seen.add(start)
        while stack:
            atom = stack.pop()
            component.append(atom)
            for neighbor in graph.adjacency[atom]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return components


def _is_linear(coords: np.ndarray) -> bool:
    if coords.shape[0] <= 2:
        return True
    centered = coords - coords.mean(axis=0)
    singular_values = np.linalg.svd(centered, compute_uv=False)
    return bool(singular_values[1] <= max(1.0e-8, 1.0e-6 * singular_values[0]))


def _primitive_basis(coordinates: Iterable[GICForgePythonCoordinate]) -> tuple[Primitive, ...]:
    basis: list[Primitive] = []
    seen: set[Primitive] = set()
    for coordinate in coordinates:
        for _coefficient, primitive in coordinate.terms:
            if primitive in seen:
                continue
            seen.add(primitive)
            basis.append(primitive)
    return tuple(basis)


def _format_readgic(name: str, terms: tuple[tuple[float, Primitive], ...]) -> str:
    return f"{name}={_format_terms(terms)}"


def _format_terms(terms: tuple[tuple[float, Primitive], ...]) -> str:
    chunks = []
    for coefficient, primitive in terms:
        atom_text = ",".join(str(atom + 1) for atom in primitive.atoms)
        symbol = {
            "bond": "R",
            "angle": "A",
            "linear_bend": "L",
            "dihedral": "D",
            "out_of_plane": "U",
        }[primitive.kind]
        if primitive.kind == "linear_bend":
            atom_text = f"{atom_text},{primitive.mode}"
        chunks.append(f"{coefficient:.10g}*{symbol}({atom_text})")
    return "+".join(chunks)
