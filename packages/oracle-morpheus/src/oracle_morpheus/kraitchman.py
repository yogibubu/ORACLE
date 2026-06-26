from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from oracle_chem.inertia import center_of_mass, inertia_tensor
from oracle_chem.isotopes_table import get_default_isotope, get_isotope
from oracle_chem.physical_constants import Phy, get_physical_constants
from oracle_chem.structure import Structure
from oracle_chem.topology.elements import atomic_number

from .contracts import IsotopologueObservation


ROTCONST_TO_MOMENT = (
    get_physical_constants()[Phy.PLANCK]
    / (8.0 * np.pi**2 * get_physical_constants()[Phy.TO_KG] * (1.0e-10) ** 2)
    * 1.0e-6
)


@dataclass(frozen=True)
class KraitchmanComparison:
    isotopologue: str
    atom_index: int
    atom: str
    isotope_mass_number: int
    coordinate: str
    kraitchman_abs_angstrom: float
    fitted_abs_angstrom: float
    difference_angstrom: float
    signed_kraitchman_angstrom: float = 0.0
    signed_reference_angstrom: float = 0.0
    substitution_mass_amu: float = 0.0


@dataclass(frozen=True)
class KraitchmanSeedResult:
    coordinates_angstrom: np.ndarray
    rows: tuple[KraitchmanComparison, ...]
    method: str
    fitted_atom_indices: tuple[int, ...]
    rms_atom_displacement_angstrom: float


def kraitchman_comparison(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
) -> tuple[KraitchmanComparison, ...]:
    parent = next((obs for obs in observations if not obs.substitutions), None)
    if parent is None:
        return ()
    axis_coords = principal_axis_coordinates(atoms, coords)
    parent_moments = np.array(constants_to_moments(parent.corrected.as_tuple()), dtype=float)
    parent_mass = parent_total_mass(atoms, coords)
    rows: list[KraitchmanComparison] = []
    for obs in observations:
        if len(obs.substitutions) != 1:
            continue
        atom_index, isotope_a = next(iter(obs.substitutions.items()))
        atom_pos = atom_index - 1
        if atom_pos < 0 or atom_pos >= len(atoms):
            continue
        substitution_mass = kraitchman_substitution_mass(atoms[atom_pos], int(isotope_a), parent_mass)
        if substitution_mass <= 0.0:
            continue
        moments = np.array(constants_to_moments(obs.corrected.as_tuple()), dtype=float)
        delta = moments - parent_moments
        squared = (
            (delta[1] + delta[2] - delta[0]) / (2.0 * substitution_mass),
            (delta[0] + delta[2] - delta[1]) / (2.0 * substitution_mass),
            (delta[0] + delta[1] - delta[2]) / (2.0 * substitution_mass),
        )
        for axis, value, fitted in zip(("a", "b", "c"), squared, axis_coords[atom_pos]):
            kraitchman_abs = float(np.sqrt(max(value, 0.0)))
            signed = float(np.copysign(kraitchman_abs, fitted)) if kraitchman_abs else 0.0
            fitted_abs = float(abs(fitted))
            rows.append(
                KraitchmanComparison(
                    isotopologue=obs.label,
                    atom_index=atom_index,
                    atom=str(atoms[atom_pos]),
                    isotope_mass_number=int(isotope_a),
                    coordinate=axis,
                    kraitchman_abs_angstrom=kraitchman_abs,
                    fitted_abs_angstrom=fitted_abs,
                    difference_angstrom=kraitchman_abs - fitted_abs,
                    signed_kraitchman_angstrom=signed,
                    signed_reference_angstrom=float(fitted),
                    substitution_mass_amu=float(substitution_mass),
                )
            )
    return tuple(rows)


def kraitchman_seed_geometry(
    atoms: list[str] | tuple[str, ...],
    coords: np.ndarray,
    observations: tuple[IsotopologueObservation, ...],
    rows: tuple[KraitchmanComparison, ...] | None = None,
) -> KraitchmanSeedResult | None:
    kraitchman_rows = rows if rows is not None else kraitchman_comparison(atoms, coords, observations)
    atom_targets = _atom_targets(kraitchman_rows)
    if not atom_targets:
        return None
    axis_coords = principal_axis_coordinates(atoms, coords)
    seeded = axis_coords.copy()
    atom_indices = tuple(sorted(atom_targets))
    reference = np.array([axis_coords[idx - 1] for idx in atom_indices], dtype=float)
    target = np.array([atom_targets[idx] for idx in atom_indices], dtype=float)
    method = "direct_substitution"
    if len(atom_indices) >= 3 and _rank3(reference) >= 2 and _rank3(target) >= 2:
        rotated, ok = _rigid_kabsch_update(axis_coords, reference, target)
        if ok:
            seeded = rotated
            method = "rigid_kabsch_plus_exact_substitution"
    for atom_index, target_coord in atom_targets.items():
        seeded[atom_index - 1] = target_coord
    displacements = seeded[[idx - 1 for idx in atom_indices]] - reference
    rms = float(np.sqrt(np.mean(np.sum(displacements * displacements, axis=1)))) if atom_indices else 0.0
    return KraitchmanSeedResult(
        coordinates_angstrom=seeded,
        rows=kraitchman_rows,
        method=method,
        fitted_atom_indices=atom_indices,
        rms_atom_displacement_angstrom=rms,
    )


def constants_to_moments(constants: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(ROTCONST_TO_MOMENT / value if value > 0.0 else 0.0 for value in constants)


def principal_axis_coordinates(atoms: list[str] | tuple[str, ...], coords: np.ndarray) -> np.ndarray:
    structure = Structure.from_atoms_coords(list(atoms), [tuple(row) for row in np.asarray(coords, dtype=float)])
    inertia = inertia_tensor(structure, isotopic=True)
    eigvals, eigvecs = np.linalg.eigh(inertia)
    order = np.argsort(eigvals)
    eigvecs = eigvecs[:, order]
    if np.linalg.det(eigvecs) < 0.0:
        eigvecs[:, -1] *= -1.0
    centered = np.asarray(coords, dtype=float) - center_of_mass(structure, isotopic=True)
    return centered @ eigvecs


def parent_total_mass(atoms: list[str] | tuple[str, ...], coords: np.ndarray) -> float:
    structure = Structure.from_atoms_coords(list(atoms), [tuple(row) for row in np.asarray(coords, dtype=float)])
    return float(sum(structure.mass_isotope))


def kraitchman_substitution_mass(atom: str, isotope_a: int, parent_mass: float) -> float:
    z_number = atomic_number(atom)
    if z_number is None:
        return 0.0
    default_iso = get_default_isotope(int(z_number))
    substituted_iso = get_isotope(int(z_number), int(isotope_a))
    if default_iso is None or substituted_iso is None:
        return 0.0
    delta_mass = float(substituted_iso.mass - default_iso.mass)
    if delta_mass <= 0.0 or parent_mass <= 0.0:
        return 0.0
    return delta_mass * parent_mass / (parent_mass + delta_mass)


def _atom_targets(rows: tuple[KraitchmanComparison, ...]) -> dict[int, np.ndarray]:
    grouped: dict[int, dict[str, float]] = {}
    for row in rows:
        grouped.setdefault(row.atom_index, {})[row.coordinate] = row.signed_kraitchman_angstrom
    targets = {}
    for atom_index, axes in grouped.items():
        if {"a", "b", "c"}.issubset(axes):
            targets[atom_index] = np.array([axes["a"], axes["b"], axes["c"]], dtype=float)
    return targets


def _rank3(points: np.ndarray) -> int:
    centered = np.asarray(points, dtype=float) - np.mean(points, axis=0)
    singular = np.linalg.svd(centered, compute_uv=False)
    if not singular.size:
        return 0
    tol = max(centered.shape) * np.finfo(float).eps * max(float(singular[0]), 1.0)
    return int(np.sum(singular > tol))


def _rigid_kabsch_update(coords: np.ndarray, reference: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, bool]:
    ref_center = np.mean(reference, axis=0)
    target_center = np.mean(target, axis=0)
    ref0 = reference - ref_center
    target0 = target - target_center
    cov = ref0.T @ target0
    try:
        u_matrix, _singular, vt_matrix = np.linalg.svd(cov)
    except np.linalg.LinAlgError:
        return coords.copy(), False
    rotation = vt_matrix.T @ u_matrix.T
    if np.linalg.det(rotation) < 0.0:
        vt_matrix[-1, :] *= -1.0
        rotation = vt_matrix.T @ u_matrix.T
    if not np.all(np.isfinite(rotation)):
        return coords.copy(), False
    if abs(np.linalg.det(rotation) - 1.0) > 1.0e-6:
        return coords.copy(), False
    return (np.asarray(coords, dtype=float) - ref_center) @ rotation.T + target_center, True
