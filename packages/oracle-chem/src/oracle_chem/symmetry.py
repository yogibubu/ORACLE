"""Molecular point-group detection and serialized symmetry operations."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import gcd
import re

import numpy as np

from .geometry import MolecularGeometry
from .topology.elements import atomic_number


@dataclass(frozen=True)
class SymmetryOperation:
    label: str
    rotation: tuple[tuple[float, float, float], ...]
    permutation: tuple[int, ...]
    max_deviation: float


@dataclass(frozen=True)
class MolecularSymmetry:
    point_group: str
    operations: tuple[SymmetryOperation, ...]
    atom_classes: tuple[tuple[int, ...], ...]
    max_deviation: float
    mean_deviation: float


def analyze_molecular_symmetry(
    geometry: MolecularGeometry,
    *,
    distance_tolerance: float,
    inertia_tolerance: float,
    max_rotation_order: int,
) -> MolecularSymmetry:
    symbols = list(geometry.atoms)
    weights = np.array([atomic_number(symbol) or 1 for symbol in symbols], dtype=float)
    oriented = orient_coords(geometry.coordinates_angstrom, weights=weights)
    elements, atom_classes, permutations = symmetry_elements_from_geometry(
        symbols,
        oriented,
        tol=distance_tolerance,
        max_n=max_rotation_order,
        tol_H=distance_tolerance,
        ignore_isotopes=True,
        auto_max_n=True,
        inertia_tol=inertia_tolerance,
    )
    if not elements:
        elements = [("E", np.eye(3), 0.0)]
        permutations = [tuple(range(len(symbols)))]
        atom_classes = tuple((idx,) for idx in range(len(symbols)))
    operations = tuple(
        SymmetryOperation(
            label=str(label),
            rotation=tuple(tuple(float(value) for value in row) for row in rotation),
            permutation=tuple(int(item) + 1 for item in permutation),
            max_deviation=float(max_deviation),
        )
        for (label, rotation, max_deviation), permutation in zip(elements, permutations)
    )
    return MolecularSymmetry(
        point_group=group_label(elements, linear=is_linear(oriented, tol=distance_tolerance)),
        operations=operations,
        atom_classes=tuple(tuple(int(atom) + 1 for atom in cls) for cls in atom_classes),
        max_deviation=float(max((op.max_deviation for op in operations), default=0.0)),
        mean_deviation=(
            float(np.mean([op.max_deviation for op in operations])) if operations else 0.0
        ),
    )


def symmetry_section_lines(symmetry: MolecularSymmetry, *, thresholds) -> list[str]:
    lines = [
        "SCHEMA oracle.xyz.symmetry.v1",
        f"POINT_GROUP {symmetry.point_group}",
        f"OPERATION_COUNT {len(symmetry.operations)}",
        f"MAX_OPERATION_DEVIATION_ANGSTROM {symmetry.max_deviation:.12g}",
        f"MEAN_OPERATION_DEVIATION_ANGSTROM {symmetry.mean_deviation:.12g}",
        f"THRESHOLD_DISTANCE_ANGSTROM {thresholds.distance_angstrom:.12g}",
        f"THRESHOLD_INERTIA_RELATIVE {thresholds.inertia_relative:.12g}",
        f"MAX_ROTATION_ORDER {thresholds.max_rotation_order}",
        "[OPERATIONS]",
    ]
    for idx, operation in enumerate(symmetry.operations, start=1):
        matrix = ",".join(
            f"{value:.12g}" for row in operation.rotation for value in row
        )
        permutation = ",".join(str(atom) for atom in operation.permutation)
        lines.append(
            f"{idx} LABEL={operation.label} "
            f"MAX_DEVIATION={operation.max_deviation:.12g} "
            f"PERMUTATION={permutation} MATRIX={matrix}"
        )
    lines.append("[ATOM_CLASSES]")
    if symmetry.atom_classes:
        for idx, atoms in enumerate(symmetry.atom_classes, start=1):
            lines.append(f"{idx} ATOMS=" + ",".join(str(atom) for atom in atoms))
    else:
        lines.append("NONE")
    return lines


def orient_coords(coords, weights=None):
    x = np.array(coords, dtype=float)
    w = np.ones(len(x), dtype=float) if weights is None else np.array(weights, dtype=float)
    center = np.sum(x * w[:, None], axis=0) / max(float(np.sum(w)), 1.0e-12)
    x = x - center
    inertia = np.zeros((3, 3), dtype=float)
    for idx, vec in enumerate(x):
        inertia += w[idx] * ((np.dot(vec, vec) * np.eye(3)) - np.outer(vec, vec))
    evals, evecs = np.linalg.eigh(inertia)
    frame = evecs[:, np.argsort(evals)]
    if np.linalg.det(frame) < 0.0:
        frame[:, -1] *= -1.0
    return x @ frame


def is_linear(coords, tol=1.0e-3):
    x = np.array(coords, dtype=float)
    inertia = np.zeros((3, 3), dtype=float)
    for vec in x:
        inertia += np.dot(vec, vec) * np.eye(3) - np.outer(vec, vec)
    return bool(np.linalg.eigvalsh(inertia)[0] < tol)


def symmetry_elements_from_geometry(
    symbols,
    coords_oriented,
    tol=1.0e-3,
    max_n=6,
    tol_H=None,
    ignore_isotopes=False,
    auto_max_n=False,
    inertia_tol=1.0e-3,
):
    coords = np.asarray(coords_oriented, dtype=float)
    radii = np.linalg.norm(coords, axis=1)
    max_radius = float(np.max(radii)) if len(radii) else 1.0
    if max_radius <= 0.0:
        max_radius = 1.0
    scaled = coords / max_radius
    sym_use = [symbol[0] for symbol in symbols] if ignore_isotopes else list(symbols)
    if auto_max_n:
        inertia = np.zeros((3, 3), dtype=float)
        for vec in coords:
            inertia += np.dot(vec, vec) * np.eye(3) - np.outer(vec, vec)
        evals = np.linalg.eigvalsh(inertia)
        max_inertia = float(np.max(evals)) if len(evals) else 0.0
        if max_inertia > 0.0:
            d01 = abs(evals[0] - evals[1]) / max_inertia
            d12 = abs(evals[1] - evals[2]) / max_inertia
            if d01 > inertia_tol and d12 > inertia_tol:
                max_n = min(max_n, 2)
    elements = []
    permutations = []
    seen: set[tuple[tuple[int, ...], tuple[float, ...]]] = set()
    for label, rotation in candidate_ops(max_n=max_n):
        mapped, max_dev = _match_with_map(
            sym_use,
            scaled,
            scaled @ rotation.T,
            tol,
            tol_H=tol_H,
        )
        if mapped is not None:
            unique_key = (
                tuple(int(item) for item in mapped),
                tuple(round(float(value), 10) for value in rotation.reshape(-1)),
            )
            if unique_key in seen:
                continue
            seen.add(unique_key)
            elements.append((label, rotation, float(max_dev)))
            permutations.append(tuple(mapped))
    return elements, _atom_classes(len(symbols), permutations), permutations


def group_label(elements, linear=False):
    labels = [item[0] for item in elements]
    nmax, axis = _highest_cn_axis(labels)
    has_i = "i" in labels
    has_sigma = any(label.startswith("sigma") for label in labels)
    has_c2 = any(label.startswith("C2") for label in labels)
    if linear:
        return "Dinfh" if has_i else "Cinfv"
    if nmax >= 2:
        sigma_h = {"x": "sigma_yz", "y": "sigma_xz", "z": "sigma_xy"}.get(axis or "z")
        has_sigma_h = sigma_h in labels
        has_sigma_v = has_sigma and not has_sigma_h
        if has_sigma_h and has_c2:
            return f"D{nmax}h"
        if has_sigma_h:
            return f"C{nmax}h"
        if has_sigma_v:
            return f"C{nmax}v"
        if has_c2:
            return f"D{nmax}"
        return f"C{nmax}"
    if has_i:
        return "Ci"
    if has_sigma:
        return "Cs"
    return "C1"


def _match_with_map(symbols, coords1, coords2, tol, tol_H=None):
    used = np.zeros(len(coords2), dtype=bool)
    mapping = [-1] * len(coords1)
    by_symbol: dict[str, list[int]] = {}
    for idx, symbol in enumerate(symbols):
        by_symbol.setdefault(symbol, []).append(idx)
    radii2 = np.linalg.norm(coords2, axis=1)
    max_dev = 0.0
    for idx in sorted(range(len(coords1)), key=lambda item: (len(by_symbol[symbols[item]]), item)):
        eff_tol = tol_H if tol_H is not None and symbols[idx] == "H" else tol
        radius = float(np.linalg.norm(coords1[idx]))
        candidates = [
            cand
            for cand in by_symbol.get(symbols[idx], ())
            if not used[cand] and abs(radii2[cand] - radius) <= eff_tol
        ]
        candidates.sort(key=lambda cand: abs(radii2[cand] - radius))
        for cand in candidates:
            deviation = float(np.linalg.norm(coords1[idx] - coords2[cand]))
            if deviation < eff_tol:
                mapping[idx] = cand
                used[cand] = True
                max_dev = max(max_dev, deviation)
                break
        if mapping[idx] < 0:
            return None, None
    return tuple(mapping), max_dev


def _atom_classes(natoms: int, permutations) -> tuple[tuple[int, ...], ...]:
    parent = list(range(natoms))

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left, right):
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for permutation in permutations:
        for left, right in enumerate(permutation):
            union(left, int(right))
    classes: dict[int, list[int]] = {}
    for atom in range(natoms):
        classes.setdefault(find(atom), []).append(atom)
    return tuple(tuple(values) for values in classes.values())


def _highest_cn_axis(labels):
    best = 1
    axis = None
    for label in labels:
        match = re.match(r"C(\d+)([xyz])", str(label))
        if match and int(match.group(1)) > best:
            best = int(match.group(1))
            axis = match.group(2)
    return best, axis


@lru_cache(maxsize=16)
def candidate_ops(max_n=6):
    ops = [("E", np.eye(3)), ("i", -np.eye(3))]
    for axis, name in [(0, "sigma_yz"), (1, "sigma_xz"), (2, "sigma_xy")]:
        matrix = np.eye(3)
        matrix[axis, axis] = -1.0
        ops.append((name, matrix))
    for n in range(2, max_n + 1):
        for power in range(1, n):
            if gcd(n, power) != 1:
                continue
            theta = 2.0 * np.pi * power / n
            ops.append((f"C{n}z^{power}", _rotation_matrix((0, 0, 1), theta)))
            ops.append((f"C{n}x^{power}", _rotation_matrix((1, 0, 0), theta)))
            ops.append((f"C{n}y^{power}", _rotation_matrix((0, 1, 0), theta)))
    for n in range(2, max_n + 1):
        for k in range(n):
            theta = np.pi * k / n
            ops.append(
                (
                    f"C2_xy_{n}_{k}",
                    _rotation_matrix((np.cos(theta), np.sin(theta), 0), np.pi),
                )
            )
    for n in range(3, max_n + 1):
        for k in range(n):
            theta = np.pi * k / n
            normal = np.array((np.cos(theta), np.sin(theta), 0.0), dtype=float)
            normal /= np.linalg.norm(normal)
            ops.append((f"sigma_v_{n}_{k}", np.eye(3) - 2.0 * np.outer(normal, normal)))
    return ops


def _rotation_matrix(axis, theta):
    axis = np.array(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(theta)
    s = np.sin(theta)
    one_c = 1.0 - c
    return np.array(
        [
            [c + x * x * one_c, x * y * one_c - z * s, x * z * one_c + y * s],
            [y * x * one_c + z * s, c + y * y * one_c, y * z * one_c - x * s],
            [z * x * one_c - y * s, z * y * one_c + x * s, c + z * z * one_c],
        ],
        dtype=float,
    )
