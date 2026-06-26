#!/usr/bin/env python3
"""Path Hamiltonians for one- and two-dimensional large-amplitude motions.

The one-dimensional workflow reads a sequence of Cartesian geometries, removes
translation and rotation by local mass-weighted alignment, computes the
mass-weighted arc-length coordinate, and solves the resulting path Hamiltonian.
The two-dimensional workflow reads a rectangular potential grid and solves a
reduced product-basis Hamiltonian after pruning the basis with two one-
dimensional reference problems.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator, UnivariateSpline
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import eigh
from scipy.optimize import least_squares, minimize, minimize_scalar


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "gaussian_outputs"
DEFAULT_FIGS = ROOT / "figs"

HARTREE_TO_CM = 219474.6313705
KJMOL_TO_CM = 83.59347225
KCALMOL_TO_CM = 349.75508885
KB_CM_PER_K = 0.695034800
AMU_TO_ME = 1822.888486209
BOHR_TO_ANGSTROM = 0.529177210903
MW_ANGSTROM_TO_AU = math.sqrt(AMU_TO_ME) / BOHR_TO_ANGSTROM
DEBYE_PER_AU = 2.541746473
TORSION_PHASE_OFFSET_DEG = 234.0

PLANCK = 6.62607015e-34
AMU_ANGSTROM2_TO_KG_M2 = 1.66053906660e-47

ELEMENTS: dict[str, tuple[int, float]] = {
    "H": (1, 1.00782503223),
    "D": (1, 2.01410177812),
    "He": (2, 4.00260325413),
    "Li": (3, 7.0160034366),
    "Be": (4, 9.012183065),
    "B": (5, 11.00930536),
    "C": (6, 12.0),
    "N": (7, 14.00307400443),
    "O": (8, 15.99491461957),
    "F": (9, 18.99840316273),
    "Ne": (10, 19.9924401762),
    "Na": (11, 22.9897692820),
    "Mg": (12, 23.985041697),
    "Al": (13, 26.98153853),
    "Si": (14, 27.97692653465),
    "P": (15, 30.97376199842),
    "S": (16, 31.9720711744),
    "Cl": (17, 34.968852682),
    "Ar": (18, 39.9623831237),
    "K": (19, 38.9637064864),
    "Ca": (20, 39.962590863),
    "Br": (35, 78.9183376),
    "I": (53, 126.9044719),
}

ATOMIC_NUMBER_TO_SYMBOL = {
    number: symbol for symbol, (number, _) in ELEMENTS.items() if symbol != "D"
}
ATOMIC_NUMBER_TO_MASS = {
    number: mass for symbol, (number, mass) in ELEMENTS.items() if symbol != "D"
}

ENERGY_KEYS = (
    "energy_hartree",
    "hartree",
    "e_h",
    "eh",
    "hf",
    "scf",
    "energy_cm-1",
    "energy_cm1",
    "v_cm-1",
    "v_cm1",
    "relative_energy_cm-1",
    "relative_energy_cm1",
    "energy_kjmol",
    "energy_kj_mol",
    "energy_kcalmol",
    "energy_kcal_mol",
    "e",
    "energy",
)

METRIC_KEYS = ("g11", "g22", "g12")
COVARIANT_METRIC_KEYS = ("cov11", "cov22", "cov12")

@dataclass
class Structure:
    atoms: np.ndarray
    symbols: list[str]
    coords_angstrom: np.ndarray
    props: dict[str, float] = field(default_factory=dict)
    comment: str = ""


@dataclass
class SourceData:
    structures: list[Structure]
    path: Path
    charge: int = 0
    multiplicity: int = 1


@dataclass
class GicToCremerPopleBridge:
    mode_index: int
    cp_mode: int
    matrix: np.ndarray


@dataclass
class Grid1DModel:
    grid_au: np.ndarray
    potential_cm: np.ndarray
    properties: dict[str, np.ndarray]
    expectation_properties: dict[str, np.ndarray]
    model_points_au: np.ndarray
    model_potential_cm: np.ndarray
    model_properties: dict[str, np.ndarray]
    info: dict[str, float | int | str] = field(default_factory=dict)


@dataclass
class PropertyDerivativeSpec:
    name: str
    value: float
    derivatives: dict[int, float]
    origin: str | float = "zero"
    parity: str = "auto"


@dataclass
class AnharmonicModeDerivatives:
    requested_mode: int
    force_mode: int
    mode_order: str
    f2_cm: float
    f3_cm: float
    f4_cm: float
    frequency_cm: float | None = None
    n_modes: int = 0
    source: str = "gaussian"
    notes: tuple[str, ...] = ()


@dataclass
class AnharmonicDerivativePotential:
    q: np.ndarray
    potential_cm: np.ndarray
    model_potential_cm: np.ndarray
    taylor_potential_cm: np.ndarray
    info: dict[str, float | int | str] = field(default_factory=dict)


@dataclass
class GaussianNormalMode:
    mode: int
    frequency_cm: float
    reduced_mass_amu: float
    force_constant_mdyne_A: float
    atoms: np.ndarray
    symbols: list[str]
    displacements: np.ndarray


@dataclass
class NormalModeTailSide:
    side: str
    mode: int
    force_mode: int
    frequency_cm: float
    force_constant_mdyne_A: float
    reduced_mass_amu: float
    overlap: float
    abs_overlap: float
    curvature_cm_au2: float
    diagonal_f2_cm: float
    diagonal_f3_cm: float
    diagonal_f4_cm: float
    vpt2_xe_cm: float
    vpt2_zpe_cm: float
    vpt2_fundamental_cm: float
    log_min_rmsd_angstrom: float
    scan_endpoint_index: int
    scan_neighbor_index: int
    tangent_steps: int
    notes: tuple[str, ...] = ()


@dataclass
class NormalModeTailData:
    sides: dict[str, NormalModeTailSide]
    info: dict[str, float | int | str] = field(default_factory=dict)


@dataclass
class AnharmonicVPT2Comparison:
    harmonic_levels_cm: np.ndarray
    vpt2_levels_cm: np.ndarray
    variational_levels_cm: np.ndarray
    info: dict[str, float | int | str] = field(default_factory=dict)


@dataclass
class PropertyVPT2Result:
    comparisons: dict[str, dict[str, float | str]]
    info: dict[str, float | int | str] = field(default_factory=dict)


@dataclass
class AxisBasis2D:
    grid: np.ndarray
    boundary: str
    solver: str
    reference_potential: np.ndarray
    energies: np.ndarray
    vectors: np.ndarray
    kinetic: np.ndarray
    derivative: np.ndarray
    info: dict[str, float | int | str] = field(default_factory=dict)


@dataclass
class Grid2DData:
    path: Path
    q1_key: str
    q2_key: str
    energy_key: str
    energy_unit: str
    q1: np.ndarray
    q2: np.ndarray
    potential_cm: np.ndarray
    properties: dict[str, np.ndarray] = field(default_factory=dict)
    metric: dict[str, np.ndarray] = field(default_factory=dict)
    metric_note: str = "constant"
    metric_diagnostics: dict[str, float | int | str] = field(default_factory=dict)


@dataclass
class Result2D:
    data: Grid2DData
    basis1: AxisBasis2D
    basis2: AxisBasis2D
    reference_note: str
    levels: np.ndarray
    coeff: np.ndarray
    expectations: list[dict[str, float]]
    product_basis_initial_size: int
    product_basis_final_size: int
    product_pruning_note: str


def normalized_key(key: str) -> str:
    return key.strip().lower().replace("/", "_").replace(" ", "_")


def parse_float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def parse_properties_from_comment(comment: str) -> dict[str, float]:
    props: dict[str, float] = {}
    pattern = re.compile(
        r"([A-Za-z_][A-Za-z0-9_.:/-]*)\s*[=:]\s*"
        r"([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?|[-+]?\.\d+(?:[EeDd][-+]?\d+)?)"
    )
    for key, value in pattern.findall(comment):
        props[key] = parse_float(value)
    return props


def format_float(value: float) -> str:
    return f"{value:.10f}"


def atom_from_token(token: str) -> tuple[int, str, float]:
    if token.isdigit():
        number = int(token)
        if number not in ATOMIC_NUMBER_TO_MASS:
            raise ValueError(f"No mass is defined for atomic number {number}")
        return number, ATOMIC_NUMBER_TO_SYMBOL.get(number, str(number)), ATOMIC_NUMBER_TO_MASS[number]

    symbol = token[0].upper() + token[1:].lower()
    if token == "D":
        symbol = "D"
    if symbol not in ELEMENTS:
        raise ValueError(f"No mass is defined for element {token}")
    number, mass = ELEMENTS[symbol]
    return number, symbol, mass


def read_xyz(path: Path) -> list[Structure]:
    lines = path.read_text().splitlines()
    structures: list[Structure] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        natoms = int(lines[i].strip())
        if i + natoms + 1 >= len(lines):
            raise ValueError(f"Incomplete XYZ frame starting at line {i + 1}")
        comment = lines[i + 1].strip()
        atoms: list[int] = []
        symbols: list[str] = []
        coords: list[list[float]] = []
        for line in lines[i + 2 : i + 2 + natoms]:
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"Malformed XYZ atom line: {line}")
            number, symbol, _ = atom_from_token(fields[0])
            atoms.append(number)
            symbols.append(symbol)
            coords.append([parse_float(fields[1]), parse_float(fields[2]), parse_float(fields[3])])
        structures.append(
            Structure(
                atoms=np.array(atoms, dtype=int),
                symbols=symbols,
                coords_angstrom=np.array(coords, dtype=float),
                props=parse_properties_from_comment(comment),
                comment=comment,
            )
        )
        i += natoms + 2
    return structures


def parse_orientation_block(lines: list[str], start: int) -> tuple[np.ndarray, list[str], np.ndarray, int]:
    i = start + 1
    dash_count = 0
    while i < len(lines):
        if set(lines[i].strip()) == {"-"}:
            dash_count += 1
            if dash_count == 2:
                i += 1
                break
        i += 1

    atoms: list[int] = []
    symbols: list[str] = []
    coords: list[list[float]] = []
    while i < len(lines):
        if set(lines[i].strip()) == {"-"}:
            break
        fields = lines[i].split()
        if len(fields) >= 6:
            number = int(fields[1])
            atoms.append(number)
            symbols.append(ATOMIC_NUMBER_TO_SYMBOL.get(number, str(number)))
            coords.append([float(fields[3]), float(fields[4]), float(fields[5])])
        i += 1
    return np.array(atoms, dtype=int), symbols, np.array(coords, dtype=float), i


def read_gaussian_charge_multiplicity(path: Path) -> tuple[int, int]:
    charge = 0
    multiplicity = 1
    pattern = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
    for line in path.read_text(errors="ignore").splitlines():
        match = pattern.search(line)
        if match:
            charge = int(match.group(1))
            multiplicity = int(match.group(2))
            break
    return charge, multiplicity


def read_gaussian_log(
    path: Path,
    selection: str = "all",
    energy_source: str = "auto",
) -> list[Structure]:
    lines = path.read_text(errors="ignore").splitlines()
    structures: list[Structure] = []
    scan_structures: list[Structure] = []
    scan_blocks: list[list[Structure]] = []
    scan_link_indices: list[int] = []
    scan_link_last_indices: list[int] = []
    link_last_indices: list[int] = []
    link_indices: list[int] = []
    last_atoms: np.ndarray | None = None
    last_symbols: list[str] | None = None
    last_coords: np.ndarray | None = None
    preferred_scan_atoms: np.ndarray | None = None
    preferred_scan_symbols: list[str] | None = None
    preferred_scan_coords: np.ndarray | None = None
    last_scf_energy: float | None = None
    last_post_scf_energy: float | None = None
    last_rot_constants_ghz: tuple[float, float, float] | None = None
    last_dipole_debye: tuple[float, float, float] | None = None
    last_gic_values: dict[str, float] = {}
    pending_optimized = False
    pending_scf_energy: float | None = None
    pending_post_scf_energy: float | None = None

    if energy_source not in {"auto", "scf", "post-scf"}:
        raise ValueError(f"Unsupported Gaussian energy source: {energy_source}")

    gaussian_number = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[DEde][-+]?\d+)?"
    scf_re = re.compile(rf"SCF Done:\s+E\([^)]*\)\s*=\s*({gaussian_number})")
    eparen_re = re.compile(
        r"(?<![A-Za-z0-9])E\(([^)\n]+)\)\s*=\s*"
        rf"({gaussian_number})"
    )
    emp2_re = re.compile(
        rf"\bE[UR]?MP2\s*=\s*({gaussian_number})",
        re.IGNORECASE,
    )
    rot_re = re.compile(
        r"Rotational constants \(GHZ\):\s+"
        rf"({gaussian_number})\s+({gaussian_number})\s+({gaussian_number})",
        re.IGNORECASE,
    )
    dipole_debye_re = re.compile(
        rf"X=\s*({gaussian_number})\s+Y=\s*({gaussian_number})\s+"
        rf"Z=\s*({gaussian_number})\s+Tot=\s*({gaussian_number})",
        re.IGNORECASE,
    )
    dipole_au_re = re.compile(r"\s*Dipole\s*=\s*(.*)$", re.IGNORECASE)
    scan_step_re = re.compile(
        r"Step number\s+(\d+)\b.*?\bon scan point\s+(\d+)"
        r"(?:\s+out of\s+(\d+))?",
        re.IGNORECASE,
    )
    gic_table_re = re.compile(
        rf"!\s*([A-Za-z][A-Za-z0-9_]*?)\s+GIC-\d+\s+({gaussian_number})\b"
    )
    gic_step_re = re.compile(
        rf"^\s*([A-Za-z][A-Za-z0-9_]*?)\s+"
        rf"({gaussian_number})(?:\s+{gaussian_number}){{4}}\s+({gaussian_number})\s*$"
    )

    def gic_prop_key(name: str) -> str:
        return f"gic_{name}"

    def store_gic_value(name: str, value: float) -> None:
        last_gic_values[gic_prop_key(name)] = value
        # If Gaussian prints the final GIC table after the scan-point record,
        # keep the selected scan structure synchronized with the latest values.
        if scan_link_indices and preferred_scan_atoms is None:
            scan_structures[scan_link_indices[-1]].props[gic_prop_key(name)] = value

    def selected_energy(
        scf_energy: float | None,
        post_scf_energy: float | None,
    ) -> tuple[float | None, str]:
        if energy_source == "scf":
            return scf_energy, "scf"
        if energy_source == "post-scf":
            return post_scf_energy, "post-scf"
        if post_scf_energy is not None:
            return post_scf_energy, "post-scf"
        return scf_energy, "scf"

    def structure_props(
        scf_energy: float | None,
        post_scf_energy: float | None,
        rot_constants_ghz: tuple[float, float, float] | None = None,
        dipole_debye: tuple[float, float, float] | None = None,
    ) -> dict[str, float]:
        props: dict[str, float] = {}
        if scf_energy is not None:
            props["energy_scf_hartree"] = scf_energy
        if post_scf_energy is not None:
            props["energy_post_scf_hartree"] = post_scf_energy
        energy, _source = selected_energy(scf_energy, post_scf_energy)
        if energy is not None:
            props["energy_hartree"] = energy
        if rot_constants_ghz is not None:
            a, b, c = rot_constants_ghz
            props["A_GHz"] = a
            props["B_GHz"] = b
            props["C_GHz"] = c
            props["A_MHz"] = 1000.0 * a
            props["B_MHz"] = 1000.0 * b
            props["C_MHz"] = 1000.0 * c
        if dipole_debye is not None:
            mux, muy, muz = dipole_debye
            props["dipole_x_debye"] = mux
            props["dipole_y_debye"] = muy
            props["dipole_z_debye"] = muz
            props["dipole_debye"] = math.sqrt(mux * mux + muy * muy + muz * muz)
        props.update(last_gic_values)
        return props

    def scan_geometry() -> tuple[np.ndarray, list[str], np.ndarray]:
        if preferred_scan_atoms is not None and preferred_scan_symbols is not None and preferred_scan_coords is not None:
            return preferred_scan_atoms, preferred_scan_symbols, preferred_scan_coords
        if last_atoms is not None and last_symbols is not None and last_coords is not None:
            return last_atoms, last_symbols, last_coords
        raise ValueError("Found a scan-point record without a preceding orientation block")

    def append_pending_optimized(rot_constants_ghz: tuple[float, float, float] | None = None) -> None:
        nonlocal pending_optimized, pending_scf_energy, pending_post_scf_energy
        if not pending_optimized:
            return
        if last_atoms is None or last_symbols is None or last_coords is None:
            raise ValueError("Found an optimized structure without a final orientation block")

        props = structure_props(
            pending_scf_energy,
            pending_post_scf_energy,
            rot_constants_ghz,
            last_dipole_debye,
        )
        structures.append(
            Structure(
                atoms=last_atoms.copy(),
                symbols=list(last_symbols),
                coords_angstrom=last_coords.copy(),
                props=props,
                comment="Optimization completed",
            )
        )
        link_indices.append(len(structures) - 1)
        pending_optimized = False
        pending_scf_energy = None
        pending_post_scf_energy = None

    def append_scan_step(step: int, point: int, total_points: int) -> None:
        nonlocal preferred_scan_atoms, preferred_scan_symbols, preferred_scan_coords
        if scan_link_indices:
            previous = int(scan_structures[scan_link_indices[-1]].props["scan_point"])
            if point < previous:
                finalize_scan_block()
        atoms, symbols, coords = scan_geometry()
        props = structure_props(
            last_scf_energy,
            last_post_scf_energy,
            last_rot_constants_ghz,
            last_dipole_debye,
        )
        props["scan_point"] = float(point)
        props["scan_step"] = float(step)
        props["scan_points_total"] = float(total_points)
        scan_structures.append(
            Structure(
                atoms=atoms.copy(),
                symbols=list(symbols),
                coords_angstrom=coords.copy(),
                props=props,
                comment=f"scan_point={point} scan_step={step}",
            )
        )
        scan_link_indices.append(len(scan_structures) - 1)
        preferred_scan_atoms = None
        preferred_scan_symbols = None
        preferred_scan_coords = None

    def finalize_scan_block() -> None:
        nonlocal scan_link_indices
        if not scan_link_indices:
            return
        by_point: dict[int, Structure] = {}
        for index in scan_link_indices:
            structure = scan_structures[index]
            by_point[int(structure.props["scan_point"])] = structure
        scan_blocks.append([by_point[point] for point in sorted(by_point)])
        scan_link_last_indices.append(scan_link_indices[-1])
        scan_link_indices = []

    i = 0
    while i < len(lines):
        line = lines[i]
        match = scf_re.search(line)
        if match:
            last_scf_energy = parse_float(match.group(1))
            last_post_scf_energy = None
            last_dipole_debye = None

        match = gic_table_re.search(line)
        if match:
            store_gic_value(match.group(1), parse_float(match.group(2)))

        match = gic_step_re.match(line)
        if match:
            store_gic_value(match.group(1), parse_float(match.group(3)))

        for match in eparen_re.finditer(line):
            label = match.group(1).strip().lower()
            if label in {"corr", "correlation"}:
                continue
            last_post_scf_energy = parse_float(match.group(2))

        match = emp2_re.search(line)
        if match:
            last_post_scf_energy = parse_float(match.group(1))

        if "Input orientation:" in line or "Z-Matrix orientation:" in line:
            last_atoms, last_symbols, last_coords, i = parse_orientation_block(lines, i)
            preferred_scan_atoms = last_atoms
            preferred_scan_symbols = last_symbols
            preferred_scan_coords = last_coords
        elif "Standard orientation:" in line:
            last_atoms, last_symbols, last_coords, i = parse_orientation_block(lines, i)
            if preferred_scan_atoms is None:
                preferred_scan_atoms = last_atoms
                preferred_scan_symbols = last_symbols
                preferred_scan_coords = last_coords

        if "Optimization completed." in line:
            pending_optimized = True
            pending_scf_energy = last_scf_energy
            pending_post_scf_energy = last_post_scf_energy

        match = rot_re.search(line)
        if match:
            last_rot_constants_ghz = tuple(parse_float(match.group(j)) for j in (1, 2, 3))
            if pending_optimized:
                append_pending_optimized(last_rot_constants_ghz)

        if "Dipole moment" in line and "Debye" in line and i + 1 < len(lines):
            match = dipole_debye_re.search(lines[i + 1])
            if match:
                last_dipole_debye = tuple(parse_float(match.group(j)) for j in (1, 2, 3))

        match = dipole_au_re.match(line)
        if match:
            fields = re.findall(gaussian_number, match.group(1))
            if len(fields) >= 3:
                last_dipole_debye = tuple(DEBYE_PER_AU * parse_float(fields[j]) for j in range(3))

        match = scan_step_re.search(line)
        if match:
            total_text = match.group(3)
            append_scan_step(
                step=int(match.group(1)),
                point=int(match.group(2)),
                total_points=int(total_text) if total_text is not None else 0,
            )

        if "Normal termination of Gaussian" in line:
            finalize_scan_block()
            append_pending_optimized()
            if link_indices:
                link_last_indices.append(link_indices[-1])
                link_indices = []
        i += 1

    append_pending_optimized()
    finalize_scan_block()
    if link_indices:
        link_last_indices.append(link_indices[-1])

    if scan_structures:
        scan_points = [structure for block in scan_blocks for structure in block]
        if energy_source == "post-scf" and any("energy_post_scf_hartree" not in s.props for s in scan_points):
            raise ValueError("Requested --gaussian-energy post-scf, but at least one scan point lacks a post-SCF energy")
        if energy_source == "scf" and any("energy_scf_hartree" not in s.props for s in scan_points):
            raise ValueError("Requested --gaussian-energy scf, but at least one scan point lacks an SCF energy")
        if selection == "all":
            return scan_points
        if selection == "last":
            return scan_points[-1:] if scan_points else []
        if selection == "last-per-link":
            seen: set[int] = set()
            selected = []
            for index in scan_link_last_indices:
                if index not in seen:
                    selected.append(scan_structures[index])
                    seen.add(index)
            return selected
        raise ValueError(f"Unsupported Gaussian log selection mode: {selection}")

    if energy_source == "post-scf" and any("energy_post_scf_hartree" not in s.props for s in structures):
        raise ValueError("Requested --gaussian-energy post-scf, but at least one structure lacks a post-SCF energy")
    if energy_source == "scf" and any("energy_scf_hartree" not in s.props for s in structures):
        raise ValueError("Requested --gaussian-energy scf, but at least one structure lacks an SCF energy")

    if selection == "all":
        return structures
    if selection == "last":
        return structures[-1:] if structures else []
    if selection == "last-per-link":
        seen: set[int] = set()
        selected = []
        for index in link_last_indices:
            if index not in seen:
                selected.append(structures[index])
                seen.add(index)
        return selected
    raise ValueError(f"Unsupported Gaussian log selection mode: {selection}")


def numeric_fields(line: str) -> list[float]:
    pattern = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[DEde][-+]?\d+)?"
    return [parse_float(item) for item in re.findall(pattern, line)]


def integer_text(text: str) -> bool:
    return bool(re.fullmatch(r"[-+]?\d+", text))


def parse_inline_anharmonic_derivatives(
    lines: list[str],
    requested_mode: int,
) -> AnharmonicModeDerivatives | None:
    compact_pattern = re.compile(
        r"ANHARMONIC_MODE_DERIVATIVES\s+"
        r"mode\s*=\s*(\d+)\s+"
        r"d2\s*=\s*([-+]?\S+)\s+"
        r"d3\s*=\s*([-+]?\S+)\s+"
        r"d4\s*=\s*([-+]?\S+)",
        re.IGNORECASE,
    )
    verbose_pattern = re.compile(
        r"Anharmonic\s+mode\s+(\d+)\s+derivatives\s+cm-?1\s*:\s*"
        r"D2\s*=\s*([-+]?\S+)\s+"
        r"D3\s*=\s*([-+]?\S+)\s+"
        r"D4\s*=\s*([-+]?\S+)",
        re.IGNORECASE,
    )
    for line in lines:
        match = compact_pattern.search(line) or verbose_pattern.search(line)
        if not match:
            continue
        mode = int(match.group(1))
        if mode != requested_mode:
            continue
        return AnharmonicModeDerivatives(
            requested_mode=requested_mode,
            force_mode=mode,
            mode_order="force-table",
            f2_cm=parse_float(match.group(2)),
            f3_cm=parse_float(match.group(3)),
            f4_cm=parse_float(match.group(4)),
            frequency_cm=parse_float(match.group(2)),
            n_modes=max(mode, requested_mode),
            source="inline",
        )
    return None


def parse_gaussian_frequency_blocks(lines: list[str], n_modes: int) -> list[float]:
    limit = len(lines)
    for i, line in enumerate(lines):
        if "QUADRATIC FORCE CONSTANTS IN NORMAL MODES" in line:
            limit = i
            break

    frequencies: list[float] = []
    for line in lines[:limit]:
        if "Frequencies --" not in line:
            continue
        _, values_text = line.split("--", 1)
        frequencies.extend(numeric_fields(values_text))
    if len(frequencies) < n_modes:
        return []
    return frequencies[-n_modes:]


def parse_final_principal_axis_geometry(lines: list[str]) -> Structure | None:
    start: int | None = None
    for i, line in enumerate(lines):
        if "Principal axis orientation:" in line:
            start = i
    if start is None:
        return None

    atoms: list[int] = []
    symbols: list[str] = []
    coords: list[list[float]] = []
    in_table = False
    for line in lines[start:]:
        if set(line.strip()) == {"-"}:
            if in_table and coords:
                break
            in_table = True
            continue
        if not in_table:
            continue
        fields = line.split()
        if len(fields) == 5 and integer_text(fields[0]) and integer_text(fields[1]):
            number = int(fields[1])
            if number not in ATOMIC_NUMBER_TO_MASS:
                continue
            atoms.append(number)
            symbols.append(ATOMIC_NUMBER_TO_SYMBOL.get(number, str(number)))
            coords.append([parse_float(fields[2]), parse_float(fields[3]), parse_float(fields[4])])
    if not atoms:
        return None
    return Structure(
        atoms=np.array(atoms, dtype=int),
        symbols=symbols,
        coords_angstrom=np.array(coords, dtype=float),
    )


def parse_last_orientation_geometry(lines: list[str]) -> Structure | None:
    parsed: Structure | None = None
    i = 0
    while i < len(lines):
        if (
            "Input orientation:" in lines[i]
            or "Z-Matrix orientation:" in lines[i]
            or "Standard orientation:" in lines[i]
        ):
            atoms, symbols, coords, i = parse_orientation_block(lines, i)
            parsed = Structure(atoms=atoms, symbols=symbols, coords_angstrom=coords)
        i += 1
    return parsed


def parse_gaussian_normal_modes(path: Path) -> tuple[Structure, list[GaussianNormalMode]]:
    lines = path.read_text(errors="ignore").splitlines()
    stop = len(lines)
    for i, line in enumerate(lines):
        if "Second-order Perturbative Anharmonic Analysis" in line:
            stop = i
            break

    starts = [i for i, line in enumerate(lines[:stop]) if "Harmonic frequencies" in line]
    if not starts:
        raise ValueError(f"Could not find a Gaussian harmonic frequency block in {path}")
    start = starts[-1]

    modes: dict[int, GaussianNormalMode] = {}
    i = start
    while i < stop:
        if "Frequencies --" not in lines[i]:
            i += 1
            continue

        mode_numbers = [int(value) for value in lines[i - 2].split() if integer_text(value)]
        freqs = numeric_fields(lines[i].split("--", 1)[1])
        red_masses = numeric_fields(lines[i + 1].split("--", 1)[1])
        frc_consts = numeric_fields(lines[i + 2].split("--", 1)[1])
        n_block = len(mode_numbers)
        if not (len(freqs) >= n_block and len(red_masses) >= n_block and len(frc_consts) >= n_block):
            i += 1
            continue

        j = i + 1
        while j < stop and not lines[j].strip().startswith("Atom  AN"):
            j += 1
        if j >= stop:
            break

        block_atoms: list[int] = []
        block_symbols: list[str] = []
        block_disp = {mode: [] for mode in mode_numbers}
        row = j + 1
        while row < stop:
            fields = lines[row].split()
            if len(fields) < 2 + 3 * n_block or not (integer_text(fields[0]) and integer_text(fields[1])):
                break
            number = int(fields[1])
            block_atoms.append(number)
            block_symbols.append(ATOMIC_NUMBER_TO_SYMBOL.get(number, str(number)))
            values = [parse_float(value) for value in fields[2 : 2 + 3 * n_block]]
            for col, mode in enumerate(mode_numbers):
                block_disp[mode].append(values[3 * col : 3 * col + 3])
            row += 1

        atoms_array = np.array(block_atoms, dtype=int)
        for col, mode in enumerate(mode_numbers):
            modes[mode] = GaussianNormalMode(
                mode=mode,
                frequency_cm=float(freqs[col]),
                reduced_mass_amu=float(red_masses[col]),
                force_constant_mdyne_A=float(frc_consts[col]),
                atoms=atoms_array.copy(),
                symbols=list(block_symbols),
                displacements=np.array(block_disp[mode], dtype=float),
            )
        i = row

    if not modes:
        raise ValueError(f"Could not parse Gaussian normal-coordinate displacements in {path}")

    geometry = parse_final_principal_axis_geometry(lines) or parse_last_orientation_geometry(lines)
    if geometry is None:
        raise ValueError(f"Could not find a minimum geometry in {path}")
    return geometry, [modes[index] for index in sorted(modes)]


def parse_gaussian_anharmonic_force_tables(
    lines: list[str],
) -> tuple[dict[int, float], dict[int, float], dict[int, float]]:
    f2: dict[int, float] = {}
    f3: dict[int, float] = {}
    f4: dict[int, float] = {}
    section: str | None = None

    for line in lines:
        if "QUADRATIC FORCE CONSTANTS IN NORMAL MODES" in line:
            section = "f2"
            continue
        if "CUBIC FORCE CONSTANTS IN NORMAL MODES" in line:
            section = "f3"
            continue
        if "QUARTIC FORCE CONSTANTS IN NORMAL MODES" in line:
            section = "f4"
            continue
        if section and re.search(r"Num\. of\s+[234](?:nd|rd|th)\s+derivatives", line):
            section = None
            continue
        if section is None:
            continue

        fields = line.split()
        if section == "f2":
            if len(fields) < 5 or not (integer_text(fields[0]) and integer_text(fields[1])):
                continue
            i, j = int(fields[0]), int(fields[1])
            if i == j:
                f2[i] = parse_float(fields[2])
        elif section == "f3":
            if len(fields) < 6 or not all(integer_text(item) for item in fields[:3]):
                continue
            i, j, k = (int(item) for item in fields[:3])
            if i == j == k:
                f3[i] = parse_float(fields[3])
        elif section == "f4":
            if len(fields) < 7 or not all(integer_text(item) for item in fields[:4]):
                continue
            i, j, k, l = (int(item) for item in fields[:4])
            if i == j == k == l:
                f4[i] = parse_float(fields[4])
    return f2, f3, f4


def resolve_anharmonic_force_mode(
    requested_mode: int,
    mode_order: str,
    f2_by_force_mode: dict[int, float],
    standard_frequencies: list[float],
) -> tuple[int, float | None, tuple[str, ...]]:
    notes: list[str] = []
    n_modes = len(f2_by_force_mode)
    if requested_mode < 1 or requested_mode > n_modes:
        raise ValueError(f"--anharmonic-mode must be between 1 and {n_modes}")
    if mode_order == "force-table":
        return requested_mode, f2_by_force_mode.get(requested_mode), tuple(notes)
    if mode_order != "frequency-ascending":
        raise ValueError(f"Unsupported --anharmonic-mode-order: {mode_order}")
    if len(standard_frequencies) != n_modes:
        notes.append("standard_frequency_block_not_found")
        return requested_mode, f2_by_force_mode.get(requested_mode), tuple(notes)

    target_frequency = standard_frequencies[requested_mode - 1]
    force_mode = min(
        f2_by_force_mode,
        key=lambda idx: abs(abs(f2_by_force_mode[idx]) - abs(target_frequency)),
    )
    mismatch = abs(abs(f2_by_force_mode[force_mode]) - abs(target_frequency))
    tolerance = max(1.0e-3, 1.0e-5 * max(abs(target_frequency), 1.0))
    if mismatch > tolerance:
        notes.append(f"frequency_force_table_mismatch_cm-1={mismatch:.6g}")
    return force_mode, target_frequency, tuple(notes)


def parse_gaussian_anharmonic_derivatives(
    path: Path,
    requested_mode: int,
    mode_order: str,
) -> AnharmonicModeDerivatives:
    lines = path.read_text(errors="ignore").splitlines()
    f2_by_force_mode, f3_by_force_mode, f4_by_force_mode = parse_gaussian_anharmonic_force_tables(lines)
    if not f2_by_force_mode:
        inline = parse_inline_anharmonic_derivatives(lines, requested_mode)
        if inline is not None:
            return inline
        raise ValueError(
            "Could not find Gaussian 'QUADRATIC FORCE CONSTANTS IN NORMAL MODES' "
            "or inline ANHARMONIC_MODE_DERIVATIVES data"
        )

    n_modes = len(f2_by_force_mode)
    standard_frequencies = parse_gaussian_frequency_blocks(lines, n_modes)
    force_mode, frequency_cm, notes = resolve_anharmonic_force_mode(
        requested_mode,
        mode_order,
        f2_by_force_mode,
        standard_frequencies,
    )
    notes_list = list(notes)
    f3_missing = force_mode not in f3_by_force_mode
    f4_missing = force_mode not in f4_by_force_mode
    if f3_missing:
        notes_list.append("diagonal_F3_not_printed_assumed_zero")
    if f4_missing:
        notes_list.append("diagonal_F4_not_printed_assumed_zero")
    return AnharmonicModeDerivatives(
        requested_mode=requested_mode,
        force_mode=force_mode,
        mode_order=mode_order,
        f2_cm=f2_by_force_mode[force_mode],
        f3_cm=f3_by_force_mode.get(force_mode, 0.0),
        f4_cm=f4_by_force_mode.get(force_mode, 0.0),
        frequency_cm=frequency_cm,
        n_modes=n_modes,
        source="gaussian",
        notes=tuple(notes_list),
    )


def merge_property_csv(structures: list[Structure], path: Path) -> None:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return

    for row_index, row in enumerate(rows):
        point_text = row.get("point") or row.get("index") or row.get("frame")
        if point_text:
            target = int(float(point_text)) - 1
        else:
            target = row_index
        if target < 0 or target >= len(structures):
            raise ValueError(f"Property row {row_index + 1} targets missing point {target + 1}")

        for key, value in row.items():
            if key is None or normalized_key(key) in {"point", "index", "frame"}:
                continue
            if value is None or not value.strip():
                continue
            structures[target].props[key] = parse_float(value)


def masses_for_atoms(atoms: np.ndarray, symbols: Iterable[str]) -> np.ndarray:
    masses: list[float] = []
    for number, symbol in zip(atoms, symbols):
        if symbol == "D":
            masses.append(ELEMENTS["D"][1])
        elif int(number) in ATOMIC_NUMBER_TO_MASS:
            masses.append(ATOMIC_NUMBER_TO_MASS[int(number)])
        else:
            raise ValueError(f"No mass is defined for atom {symbol} with atomic number {number}")
    return np.array(masses, dtype=float)


def check_atom_order(structures: list[Structure]) -> None:
    if not structures:
        raise ValueError("No structures were read")
    atoms = structures[0].atoms
    natoms = len(atoms)
    for i, structure in enumerate(structures, start=1):
        if len(structure.atoms) != natoms or not np.array_equal(structure.atoms, atoms):
            raise ValueError(f"Atom order changes at point {i}")


def center_of_mass(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
    return np.sum(coords * masses[:, None], axis=0) / np.sum(masses)


def remove_translation(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
    return coords - center_of_mass(coords, masses)


def mass_weighted_kabsch(
    mobile: np.ndarray, reference: np.ndarray, masses: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Rotate mobile onto reference under the local Eckart condition.

    The mass-weighted least-squares rotation minimizes the step length after
    translation removal. For an infinitesimal path this is equivalent to
    removing the rotational component of the displacement, i.e. minimizing the
    angular momentum associated with the step.
    """
    x = remove_translation(mobile, masses)
    y = remove_translation(reference, masses)
    covariance = x.T @ (masses[:, None] * y)
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(u @ vt) < 0.0:
        correction[-1, -1] = -1.0
    rotation = u @ correction @ vt
    return x @ rotation, rotation


def orient_path(
    structures: list[Structure], masses: np.ndarray
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray]:
    """Return locally Eckart-oriented geometries and path diagnostics."""
    oriented = [remove_translation(structures[0].coords_angstrom, masses)]
    step_lengths = np.zeros(len(structures), dtype=float)
    angular_residuals = np.zeros(len(structures), dtype=float)

    for i in range(1, len(structures)):
        aligned, _ = mass_weighted_kabsch(structures[i].coords_angstrom, oriented[i - 1], masses)
        displacement = aligned - oriented[i - 1]
        step_lengths[i] = math.sqrt(float(np.sum(masses[:, None] * displacement * displacement)))

        # Discrete Eckart residual for the reoriented step. It should be small;
        # otherwise the path length still contains rotational contamination.
        midpoint = 0.5 * (aligned + oriented[i - 1])
        angular = np.sum(masses[:, None] * np.cross(midpoint, displacement), axis=0)
        angular_residuals[i] = float(np.linalg.norm(angular))
        oriented.append(aligned)

    s_sqrtamu_angstrom = np.cumsum(step_lengths)
    return oriented, s_sqrtamu_angstrom, angular_residuals


def normalized_vector(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0:
        raise ValueError("Cannot normalize a zero vector")
    return np.asarray(vector, dtype=float) / norm


def harmonic_curvature_cm_au2(frequency_cm: float) -> float:
    return float(frequency_cm) * float(frequency_cm) / HARTREE_TO_CM


def endpoint_tangent_from_steps(
    oriented: list[np.ndarray],
    s_au: np.ndarray,
    endpoint_index: int,
    *,
    side: str,
    steps: int,
) -> tuple[np.ndarray, int, int]:
    if len(oriented) != len(s_au):
        raise ValueError("Path geometries and path coordinates have different lengths")
    n_steps = max(1, min(int(steps), len(oriented) - 1))
    if side == "right":
        start = len(oriented) - n_steps - 1
        stop = len(oriented)
        direction = 1.0
    elif side == "left":
        start = 0
        stop = n_steps + 1
        direction = -1.0
    else:
        raise ValueError(f"Unsupported endpoint side: {side}")

    x = np.asarray(s_au[start:stop], dtype=float) - float(s_au[endpoint_index])
    coords = np.asarray(oriented[start:stop], dtype=float)
    if len(x) < 2 or np.ptp(x) <= 0.0:
        raise ValueError("Endpoint tangent needs at least two distinct scan points")
    design = np.column_stack((np.ones_like(x), x))
    coeff, *_ = np.linalg.lstsq(design, coords.reshape(len(x), -1), rcond=None)
    tangent = direction * coeff[1].reshape(coords.shape[1], 3)
    neighbor_index = start if side == "right" else stop - 1
    return tangent, neighbor_index, n_steps


def selected_mode_diagonal_derivatives(
    modes: list[GaussianNormalMode],
    requested_mode: int,
    frequency_cm: float,
    f2_by_force_mode: dict[int, float],
    f3_by_force_mode: dict[int, float],
    f4_by_force_mode: dict[int, float],
) -> tuple[int, float, float, float, tuple[str, ...]]:
    notes: list[str] = []
    force_mode = requested_mode
    if f2_by_force_mode:
        standard_frequencies = [mode.frequency_cm for mode in modes]
        force_mode, _, resolve_notes = resolve_anharmonic_force_mode(
            requested_mode,
            "frequency-ascending",
            f2_by_force_mode,
            standard_frequencies,
        )
        notes.extend(resolve_notes)
    else:
        notes.append("diagonal_force_tables_not_found")

    f2_cm = float(f2_by_force_mode.get(force_mode, abs(frequency_cm)))
    f3_cm = float(f3_by_force_mode.get(force_mode, 0.0))
    f4_cm = float(f4_by_force_mode.get(force_mode, 0.0))
    if force_mode not in f3_by_force_mode:
        notes.append("diagonal_F3_not_printed_assumed_zero")
    if force_mode not in f4_by_force_mode:
        notes.append("diagonal_F4_not_printed_assumed_zero")
    return force_mode, f2_cm, f3_cm, f4_cm, tuple(notes)


def one_mode_vpt2_summary_from_diagonal_force_constants(
    f2_cm: float,
    f3_cm: float,
    f4_cm: float,
) -> tuple[float, float, float]:
    harmonic, vpt2, xe = one_mode_vpt2_levels(f2_cm, f3_cm, f4_cm, 2)
    _ = harmonic
    return float(xe), float(vpt2[0]), float(vpt2[1] - vpt2[0])


def select_normal_mode_for_endpoint(
    endpoint_side: str,
    endpoint_index: int,
    oriented: list[np.ndarray],
    s_au: np.ndarray,
    masses: np.ndarray,
    log_geometry: Structure,
    modes: list[GaussianNormalMode],
    f2_by_force_mode: dict[int, float],
    f3_by_force_mode: dict[int, float],
    f4_by_force_mode: dict[int, float],
    *,
    max_frequency_cm: float,
    use_heavy_atoms: bool,
    min_overlap: float,
    tangent_steps: int,
) -> NormalModeTailSide:
    if not np.array_equal(log_geometry.atoms, modes[0].atoms):
        raise ValueError("Normal-mode atom order does not match the Gaussian minimum geometry")
    if len(log_geometry.atoms) != len(masses):
        raise ValueError("Normal-mode atom count does not match the DVR path")

    endpoint = oriented[endpoint_index]
    aligned_min, rotation = mass_weighted_kabsch(log_geometry.coords_angstrom, endpoint, masses)
    rmsd = math.sqrt(float(np.sum(masses[:, None] * (aligned_min - endpoint) ** 2) / np.sum(masses)))

    tangent, neighbor_index, used_steps = endpoint_tangent_from_steps(
        oriented,
        s_au,
        endpoint_index,
        side=endpoint_side,
        steps=tangent_steps,
    )
    tangent = tangent - center_of_mass(tangent, masses)
    mask = np.array([True] * len(masses), dtype=bool)
    if use_heavy_atoms:
        mask = np.array([symbol != "H" and symbol != "D" for symbol in log_geometry.symbols], dtype=bool)
        if int(np.sum(mask)) < 2:
            raise ValueError("Heavy-atom normal-mode overlap needs at least two non-H atoms")

    sqrt_masses = np.sqrt(masses)[:, None]
    tangent_mw = normalized_vector((sqrt_masses[mask] * tangent[mask]).reshape(-1))

    candidates = [mode for mode in modes if abs(mode.frequency_cm) <= max_frequency_cm]
    if not candidates:
        raise ValueError(
            f"No Gaussian normal modes found below --tail-normal-mode-max-frequency {max_frequency_cm:g} cm^-1"
        )

    best: NormalModeTailSide | None = None
    for mode in candidates:
        if not np.array_equal(mode.atoms, log_geometry.atoms):
            raise ValueError("Normal-mode atom order changes between Gaussian frequency blocks")
        rotated_displacement = mode.displacements @ rotation
        mode_mw = normalized_vector((sqrt_masses[mask] * rotated_displacement[mask]).reshape(-1))
        overlap = float(np.dot(tangent_mw, mode_mw))
        force_mode, f2_cm, f3_cm, f4_cm, notes = selected_mode_diagonal_derivatives(
            modes,
            mode.mode,
            mode.frequency_cm,
            f2_by_force_mode,
            f3_by_force_mode,
            f4_by_force_mode,
        )
        try:
            vpt2_xe, vpt2_zpe, vpt2_fundamental = one_mode_vpt2_summary_from_diagonal_force_constants(
                f2_cm,
                f3_cm,
                f4_cm,
            )
        except ValueError:
            vpt2_xe, vpt2_zpe, vpt2_fundamental = 0.0, 0.0, 0.0
            notes = (*notes, "vpt2_not_available_nonpositive_F2")
        side = NormalModeTailSide(
            side=endpoint_side,
            mode=mode.mode,
            force_mode=force_mode,
            frequency_cm=mode.frequency_cm,
            force_constant_mdyne_A=mode.force_constant_mdyne_A,
            reduced_mass_amu=mode.reduced_mass_amu,
            overlap=overlap,
            abs_overlap=abs(overlap),
            curvature_cm_au2=harmonic_curvature_cm_au2(abs(mode.frequency_cm)),
            diagonal_f2_cm=f2_cm,
            diagonal_f3_cm=f3_cm,
            diagonal_f4_cm=f4_cm,
            vpt2_xe_cm=vpt2_xe,
            vpt2_zpe_cm=vpt2_zpe,
            vpt2_fundamental_cm=vpt2_fundamental,
            log_min_rmsd_angstrom=rmsd,
            scan_endpoint_index=endpoint_index,
            scan_neighbor_index=neighbor_index,
            tangent_steps=used_steps,
            notes=notes,
        )
        if best is None or side.abs_overlap > best.abs_overlap:
            best = side

    if best is None:
        raise ValueError("Could not select a normal mode for the endpoint tail")
    if best.abs_overlap < min_overlap:
        raise ValueError(
            f"Best endpoint normal-mode overlap is {best.abs_overlap:.4f}, below "
            f"--tail-normal-mode-min-overlap {min_overlap:.4f}"
        )
    return best


def build_normal_mode_tail_data(
    log_path: Path,
    oriented: list[np.ndarray],
    s_au: np.ndarray,
    rel_energy_cm: np.ndarray,
    masses: np.ndarray,
    args: argparse.Namespace,
) -> NormalModeTailData:
    if len(oriented) < 2:
        raise ValueError("Normal-mode tail selection needs at least two path geometries")
    lines = log_path.read_text(errors="ignore").splitlines()
    log_geometry, modes = parse_gaussian_normal_modes(log_path)
    f2_by_force_mode, f3_by_force_mode, f4_by_force_mode = parse_gaussian_anharmonic_force_tables(lines)
    if len(log_geometry.atoms) != len(masses) or not np.array_equal(log_geometry.atoms, modes[0].atoms):
        raise ValueError("Gaussian normal-mode log is incompatible with the path atom order")

    requested_side = args.tail_normal_mode_side
    sides: list[str]
    if requested_side == "minimum":
        min_index = int(np.argmin(rel_energy_cm))
        if min_index == 0:
            sides = ["left"]
        elif min_index == len(oriented) - 1:
            sides = ["right"]
        else:
            left_gap = abs(float(rel_energy_cm[0] - rel_energy_cm[min_index]))
            right_gap = abs(float(rel_energy_cm[-1] - rel_energy_cm[min_index]))
            sides = ["left" if left_gap <= right_gap else "right"]
    elif requested_side == "both":
        sides = ["left", "right"]
    else:
        sides = [requested_side]

    selected: dict[str, NormalModeTailSide] = {}
    for side in sides:
        if side == "left":
            endpoint_index = 0
        elif side == "right":
            endpoint_index = len(oriented) - 1
        else:
            raise ValueError(f"Unsupported normal-mode tail side: {side}")
        selected[side] = select_normal_mode_for_endpoint(
            side,
            endpoint_index,
            oriented,
            s_au,
            masses,
            log_geometry,
            modes,
            f2_by_force_mode,
            f3_by_force_mode,
            f4_by_force_mode,
            max_frequency_cm=float(args.tail_normal_mode_max_frequency),
            use_heavy_atoms=bool(args.tail_normal_mode_heavy_atoms),
            min_overlap=float(args.tail_normal_mode_min_overlap),
            tangent_steps=int(args.tail_normal_mode_steps),
        )

    info: dict[str, float | int | str] = {
        "normal_mode_tail_log": str(log_path),
        "normal_mode_tail_side_requested": requested_side,
        "normal_mode_tail_max_frequency_cm-1": float(args.tail_normal_mode_max_frequency),
        "normal_mode_tail_heavy_atoms": int(bool(args.tail_normal_mode_heavy_atoms)),
        "normal_mode_tail_steps": int(args.tail_normal_mode_steps),
    }
    for side, selected_side in selected.items():
        prefix = f"{side}_normal_mode_tail"
        info.update(
            {
                f"{prefix}_mode": int(selected_side.mode),
                f"{prefix}_force_mode": int(selected_side.force_mode),
                f"{prefix}_frequency_cm-1": float(selected_side.frequency_cm),
                f"{prefix}_force_constant_mdyne_A": float(selected_side.force_constant_mdyne_A),
                f"{prefix}_reduced_mass_amu": float(selected_side.reduced_mass_amu),
                f"{prefix}_overlap": float(selected_side.overlap),
                f"{prefix}_abs_overlap": float(selected_side.abs_overlap),
                f"{prefix}_curvature_cm_au2": float(selected_side.curvature_cm_au2),
                f"{prefix}_diagonal_F2_cm-1": float(selected_side.diagonal_f2_cm),
                f"{prefix}_diagonal_F3_cm-1": float(selected_side.diagonal_f3_cm),
                f"{prefix}_diagonal_F4_cm-1": float(selected_side.diagonal_f4_cm),
                f"{prefix}_vpt2_xe_cm-1": float(selected_side.vpt2_xe_cm),
                f"{prefix}_vpt2_zpe_cm-1": float(selected_side.vpt2_zpe_cm),
                f"{prefix}_vpt2_fundamental_cm-1": float(selected_side.vpt2_fundamental_cm),
                f"{prefix}_log_min_rmsd_angstrom": float(selected_side.log_min_rmsd_angstrom),
                f"{prefix}_tangent_steps": int(selected_side.tangent_steps),
                f"{prefix}_notes": ";".join(selected_side.notes),
            }
        )
    return NormalModeTailData(sides=selected, info=info)


def rotational_constants_mhz(coords: np.ndarray, masses: np.ndarray) -> tuple[float, float, float]:
    centered = remove_translation(coords, masses)
    inertia = np.zeros((3, 3), dtype=float)
    for mass, (x, y, z) in zip(masses, centered):
        r2 = x * x + y * y + z * z
        inertia += mass * (r2 * np.eye(3) - np.outer([x, y, z], [x, y, z]))

    moments = np.linalg.eigvalsh(inertia)
    moments = np.sort(np.maximum(moments, 0.0))
    constants: list[float] = []
    for moment in moments:
        if moment < 1.0e-12:
            constants.append(float("inf"))
        else:
            constants.append(PLANCK / (8.0 * math.pi**2 * moment * AMU_ANGSTROM2_TO_KG_M2) / 1.0e6)
    return constants[0], constants[1], constants[2]


def add_rotational_constants(structures: list[Structure], masses: np.ndarray) -> None:
    for structure in structures:
        a, b, c = rotational_constants_mhz(structure.coords_angstrom, masses)
        structure.props["A_MHz_calc"] = a
        structure.props["B_MHz_calc"] = b
        structure.props["C_MHz_calc"] = c


def dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.atan2(y, x)


def canonical_ring_indices(ring_indices: list[int]) -> list[int]:
    """Return ORACLE's deterministic cyclic ring numbering."""
    indices = list(ring_indices)
    if len(indices) < 4:
        raise ValueError("A canonical ring needs at least 4 atoms")
    if len(set(indices)) != len(indices):
        raise ValueError("A canonical ring cannot contain duplicated atoms")
    n = len(indices)
    start = min(range(n), key=lambda i: indices[i])
    forward = [indices[(start + i) % n] for i in range(n)]
    backward = [indices[(start - i) % n] for i in range(n)]
    return forward if forward[1] <= backward[1] else backward


def five_ring_endocyclic_torsions(coords: np.ndarray, ring_indices: list[int]) -> np.ndarray:
    ring_indices = canonical_ring_indices(ring_indices)
    ring = coords[ring_indices]
    return np.array(
        [
            dihedral(ring[0], ring[1], ring[2], ring[3]),
            dihedral(ring[1], ring[2], ring[3], ring[4]),
            dihedral(ring[2], ring[3], ring[4], ring[0]),
            dihedral(ring[3], ring[4], ring[0], ring[1]),
            dihedral(ring[4], ring[0], ring[1], ring[2]),
        ],
        dtype=float,
    )


def four_ring_pucker_torsion(coords: np.ndarray, ring_indices: list[int]) -> float:
    ring_indices = canonical_ring_indices(ring_indices)
    ring = coords[ring_indices]
    return dihedral(ring[1], ring[2], ring[3], ring[0])


def ring_coordinate_components(
    coords: np.ndarray, ring_indices: list[int]
) -> tuple[float, float, float, float]:
    if len(ring_indices) == 5:
        return five_ring_components(five_ring_endocyclic_torsions(coords, ring_indices))
    if len(ring_indices) == 4:
        puck = four_ring_pucker_torsion(coords, ring_indices)
        return puck, 0.0, abs(puck), math.degrees(puck)
    raise ValueError("Only four- and five-membered rings are supported")


def five_ring_component_coefficients() -> tuple[np.ndarray, np.ndarray]:
    scale = math.sqrt(2.0 / 5.0)
    qc0 = scale * np.array(
        [1.0, -0.8090169944, 0.3090169944, 0.3090169944, -0.8090169944],
        dtype=float,
    )
    qs0 = scale * np.array(
        [0.0, 0.5877852523, -0.9510565163, 0.9510565163, -0.5877852523],
        dtype=float,
    )
    phase = math.radians(TORSION_PHASE_OFFSET_DEG)
    qc = math.cos(phase) * qc0 + math.sin(phase) * qs0
    qs = -math.sin(phase) * qc0 + math.cos(phase) * qs0
    return qc, qs


def five_ring_components(torsions: np.ndarray) -> tuple[float, float, float, float]:
    qc_coeff, qs_coeff = five_ring_component_coefficients()
    qc = float(qc_coeff @ torsions)
    qs = float(qs_coeff @ torsions)
    q = math.hypot(qc, qs)
    phi = math.degrees(math.atan2(qs, qc)) % 360.0
    return qc, qs, q, phi


def cremer_pople_five_membered(
    coords: np.ndarray, ring_indices: list[int]
) -> tuple[float, float, float, float]:
    """Return five-membered-ring Cremer-Pople labels from Cartesian geometry.

    The values are diagnostics only. They are not used to generate constraints
    or to build the one-dimensional Hamiltonian.
    """
    ring = coords[ring_indices]
    center = np.mean(ring, axis=0)
    j = np.arange(5, dtype=float)
    angle1 = 2.0 * math.pi * j / 5.0
    plane_x = np.sum(ring * np.cos(angle1)[:, None], axis=0)
    plane_y = np.sum(ring * np.sin(angle1)[:, None], axis=0)
    normal = np.cross(plane_x, plane_y)
    norm = float(np.linalg.norm(normal))
    if norm < 1.0e-12:
        raise ValueError("Cannot define Cremer-Pople plane for the selected five-membered ring")
    normal /= norm

    z = (ring - center) @ normal
    angle2 = 4.0 * math.pi * j / 5.0
    cp_x = math.sqrt(2.0 / 5.0) * float(np.sum(z * np.cos(angle2)))
    cp_y = -math.sqrt(2.0 / 5.0) * float(np.sum(z * np.sin(angle2)))
    q2 = math.hypot(cp_x, cp_y)
    phi2 = math.degrees(math.atan2(cp_y, cp_x)) % 360.0
    return q2, phi2, cp_x, cp_y


def cremer_pople_mode_values(coords: np.ndarray, ring_indices: list[int]) -> dict[str, float]:
    """Return Cremer-Pople out-of-plane components for an arbitrary ring.

    The coordinate frame is defined by the first harmonic plane of the current
    ring geometry. Odd rings have paired modes m=2..(N-1)/2; even rings have
    paired modes m=2..N/2-1 plus the special unpaired m=N/2 coordinate.
    """
    ring_indices = canonical_ring_indices(ring_indices)
    ring_size = len(ring_indices)
    if ring_size < 4:
        raise ValueError("Cremer-Pople labels require a ring with at least four atoms")

    ring = coords[ring_indices]
    center = np.mean(ring, axis=0)
    j = np.arange(ring_size, dtype=float)
    angle1 = 2.0 * math.pi * j / float(ring_size)
    plane_x = np.sum(ring * np.cos(angle1)[:, None], axis=0)
    plane_y = np.sum(ring * np.sin(angle1)[:, None], axis=0)
    normal = np.cross(plane_x, plane_y)
    norm = float(np.linalg.norm(normal))
    if norm < 1.0e-12:
        raise ValueError("Cannot define Cremer-Pople plane for the selected ring")
    normal /= norm

    z = (ring - center) @ normal
    values: dict[str, float] = {}
    for mode in range(2, (ring_size - 1) // 2 + 1):
        angle = 2.0 * math.pi * float(mode) * j / float(ring_size)
        cp_x = math.sqrt(2.0 / float(ring_size)) * float(np.sum(z * np.cos(angle)))
        cp_y = -math.sqrt(2.0 / float(ring_size)) * float(np.sum(z * np.sin(angle)))
        q = math.hypot(cp_x, cp_y)
        phi = math.degrees(math.atan2(cp_y, cp_x)) % 360.0
        values[f"CP_m{mode}_q_angstrom"] = q
        values[f"CP_m{mode}_phi_deg"] = phi
        values[f"CP_m{mode}_x_angstrom"] = cp_x
        values[f"CP_m{mode}_y_angstrom"] = cp_y

    if ring_size % 2 == 0:
        mode = ring_size // 2
        q = math.sqrt(1.0 / float(ring_size)) * float(np.sum(((-1.0) ** j) * z))
        values[f"CP_m{mode}_q_angstrom"] = q
    return values


def gaussian_puckering_mode_values(props: dict[str, float]) -> dict[int, dict[str, float]]:
    modes: dict[int, dict[str, float]] = {}
    for key, value in props.items():
        match = re.fullmatch(r"gic_QPck0*(\d+)", key)
        if match:
            modes.setdefault(int(match.group(1)), {})["q"] = float(value)
            continue
        match = re.fullmatch(r"gic_PhiP0*(\d+)", key)
        if match:
            phi_rad = float(value)
            modes.setdefault(int(match.group(1)), {})["phi_rad"] = phi_rad
            modes[int(match.group(1))]["phi_deg"] = math.degrees(phi_rad)
    for values in modes.values():
        if "q" in values and "phi_rad" in values:
            q = values["q"]
            phi = values["phi_rad"]
            values["x"] = q * math.cos(phi)
            values["y"] = q * math.sin(phi)
    return modes


def cp_paired_modes_for_ring(ring_size: int) -> list[int]:
    if ring_size < 5:
        return []
    return list(range(2, (ring_size - 1) // 2 + 1))


def fit_gic_to_cremer_pople_bridges(
    structures: list[Structure],
    cp_ring_indices: list[int] | None,
) -> dict[int, GicToCremerPopleBridge]:
    if cp_ring_indices is None:
        return {}
    cp_modes = cp_paired_modes_for_ring(len(cp_ring_indices))
    bridges: dict[int, GicToCremerPopleBridge] = {}
    for mode_index, cp_mode in enumerate(cp_modes, start=1):
        gic_rows = []
        cp_rows = []
        for structure in structures:
            gic = gaussian_puckering_mode_values(structure.props).get(mode_index)
            if not gic or "x" not in gic or "y" not in gic:
                continue
            cp = cremer_pople_mode_values(structure.coords_angstrom, cp_ring_indices)
            x_key = f"CP_m{cp_mode}_x_angstrom"
            y_key = f"CP_m{cp_mode}_y_angstrom"
            if x_key not in cp or y_key not in cp:
                continue
            gic_rows.append([gic["x"], gic["y"]])
            cp_rows.append([cp[x_key], cp[y_key]])
        if len(gic_rows) >= 2:
            matrix, *_ = np.linalg.lstsq(np.asarray(gic_rows), np.asarray(cp_rows), rcond=None)
            bridges[mode_index] = GicToCremerPopleBridge(
                mode_index=mode_index,
                cp_mode=cp_mode,
                matrix=matrix,
            )
    return bridges


def bridge_gic_to_cremer_pople_values(
    props: dict[str, float],
    bridges: dict[int, GicToCremerPopleBridge],
) -> dict[str, float]:
    values: dict[str, float] = {}
    gic_modes = gaussian_puckering_mode_values(props)
    for mode_index, bridge in bridges.items():
        gic = gic_modes.get(mode_index)
        if not gic or "x" not in gic or "y" not in gic:
            continue
        cp_x, cp_y = np.asarray([gic["x"], gic["y"]], dtype=float) @ bridge.matrix
        q = math.hypot(float(cp_x), float(cp_y))
        phi = math.degrees(math.atan2(float(cp_y), float(cp_x))) % 360.0
        prefix = f"CP_from_GIC_m{bridge.cp_mode}"
        values[f"{prefix}_q_angstrom"] = q
        values[f"{prefix}_phi_deg"] = phi
        values[f"{prefix}_x_angstrom"] = float(cp_x)
        values[f"{prefix}_y_angstrom"] = float(cp_y)
    return values


def parse_ring_indices(text: str, natoms: int) -> list[int]:
    indices = [int(item.strip()) - 1 for item in text.split(",") if item.strip()]
    if len(indices) < 4:
        raise ValueError("--ring must contain at least four comma-separated one-based atom indices")
    if min(indices) < 0 or max(indices) >= natoms:
        raise ValueError("--ring contains an atom index outside the molecular geometry")
    if len(set(indices)) != len(indices):
        raise ValueError("--ring contains duplicated atom indices")
    return canonical_ring_indices(indices)


def find_energy_key(structures: list[Structure], requested: str | None) -> str:
    if requested:
        for key in structures[0].props:
            if normalized_key(key) == normalized_key(requested):
                return key
        raise ValueError(f"Requested energy key {requested!r} was not found")

    available = {normalized_key(key): key for key in structures[0].props}
    for candidate in ENERGY_KEYS:
        if candidate in available:
            return available[candidate]
    raise ValueError(
        "No energy column was found. Use energy_hartree, energy_cm-1, "
        "energy_kjmol, or pass --energy-key."
    )


def infer_energy_unit(key: str, requested: str) -> str:
    if requested != "auto":
        return requested
    nkey = normalized_key(key)
    if "cm" in nkey:
        return "cm-1"
    if "kj" in nkey:
        return "kjmol"
    if "kcal" in nkey:
        return "kcalmol"
    if nkey in {"v", "v_cm1"}:
        return "cm-1"
    return "hartree"


def potential_cm(
    structures: list[Structure], energy_key: str, energy_unit: str
) -> tuple[np.ndarray, np.ndarray]:
    raw = np.array([structure.props[energy_key] for structure in structures], dtype=float)
    unit = infer_energy_unit(energy_key, energy_unit)
    values = energy_values_cm(raw, unit)
    return raw, values - np.min(values)


def energy_values_cm(raw: np.ndarray, unit: str) -> np.ndarray:
    if unit == "hartree":
        values = raw * HARTREE_TO_CM
    elif unit == "cm-1":
        values = raw
    elif unit == "kjmol":
        values = raw * KJMOL_TO_CM
    elif unit == "kcalmol":
        values = raw * KCALMOL_TO_CM
    else:
        raise ValueError(f"Unsupported energy unit: {unit}")
    return values


def scalar_property_keys(
    structures: list[Structure],
    energy_key: str,
    requested: list[str] | None,
    derivative_property_names: Iterable[str] = (),
) -> list[str]:
    derivative_by_norm = {normalized_key(key): key for key in derivative_property_names}
    if requested:
        by_norm = {normalized_key(key): key for key in structures[0].props}
        keys = []
        for item in requested:
            if normalized_key(item) == "all":
                numeric_keys = []
                energy_norm = normalized_key(energy_key)
                for key in structures[0].props:
                    nkey = normalized_key(key)
                    if nkey == energy_norm or nkey in ENERGY_KEYS:
                        continue
                    if key.startswith("gic_"):
                        continue
                    try:
                        [float(structure.props[key]) for structure in structures]
                    except (KeyError, TypeError, ValueError):
                        continue
                    numeric_keys.append(key)
                keys.extend(key for key in numeric_keys if key not in keys)
                keys.extend(key for key in derivative_property_names if key not in keys)
                continue
            nitem = normalized_key(item)
            if nitem in by_norm:
                keys.append(by_norm[nitem])
            elif nitem in derivative_by_norm:
                keys.append(derivative_by_norm[nitem])
            else:
                raise ValueError(f"Requested property {item!r} was not found")
        return keys

    auto_candidates = ["A_MHz", "B_MHz", "C_MHz", "dipole_debye"]
    keys: list[str] = []
    by_norm = {normalized_key(key): key for key in structures[0].props}
    for item in auto_candidates:
        key = by_norm.get(normalized_key(item))
        if key is None:
            continue
        try:
            [float(structure.props[key]) for structure in structures]
        except (KeyError, TypeError, ValueError):
            continue
        keys.append(key)
    keys.extend(key for key in derivative_property_names if key not in keys)
    return keys


def derivative_order_from_key(key: str) -> int | None:
    normalized = normalized_key(key).replace("-", "_")
    match = re.fullmatch(r"(?:d|deriv|derivative|order)_?(\d+)", normalized)
    if match:
        return int(match.group(1))
    match = re.fullmatch(r"(?:d|deriv|derivative|order)(\d+)", normalized)
    if match:
        return int(match.group(1))
    return None


def read_property_derivative_specs(path: Path) -> dict[str, PropertyDerivativeSpec]:
    specs: dict[str, PropertyDerivativeSpec] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row_number, row in enumerate(reader, start=2):
            name = (
                row.get("property")
                or row.get("name")
                or row.get("key")
                or row.get("label")
                or ""
            ).strip()
            if not name:
                raise ValueError(f"Missing property name in {path} row {row_number}")
            value_text = row.get("value") or row.get("d0") or row.get("f0") or row.get("p0")
            if value_text is None or not value_text.strip():
                raise ValueError(f"Missing central value for property {name!r} in {path} row {row_number}")
            derivatives: dict[int, float] = {}
            for key, text in row.items():
                if key is None or text is None or not text.strip():
                    continue
                order = derivative_order_from_key(key)
                if order is None or order == 0:
                    continue
                derivatives[order] = parse_float(text)
            origin_text = (row.get("origin") or "zero").strip()
            try:
                origin: str | float = parse_float(origin_text)
            except ValueError:
                origin = origin_text.lower()
            parity = (row.get("parity") or row.get("symmetry") or "auto").strip().lower()
            if parity not in {"auto", "even", "odd", "full"}:
                raise ValueError(
                    f"Unsupported parity {parity!r} for property {name!r}; use auto, even, odd, or full"
                )
            if name in specs:
                raise ValueError(f"Duplicate property derivative specification for {name!r}")
            specs[name] = PropertyDerivativeSpec(
                name=name,
                value=parse_float(value_text),
                derivatives=derivatives,
                origin=origin,
                parity=parity,
            )
    return specs


def property_derivative_degree(spec: PropertyDerivativeSpec) -> int:
    nonzero_orders = [order for order, value in spec.derivatives.items() if abs(value) > 0.0]
    return max(nonzero_orders, default=0)


def resolve_property_derivative_origin(
    origin: str | float,
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
) -> float:
    if isinstance(origin, (int, float)):
        return float(origin)
    origin_key = origin.lower()
    if origin_key in {"zero", "center", "centre", "origin"}:
        return 0.0
    if origin_key == "first":
        return float(coordinate_au[0])
    if origin_key == "last":
        return float(coordinate_au[-1])
    if origin_key in {"minimum", "min"}:
        return float(coordinate_au[int(np.argmin(potential_cm))])
    raise ValueError(f"Unsupported property derivative origin: {origin!r}")


def evaluate_property_derivative_spec(
    spec: PropertyDerivativeSpec,
    coordinate_au: np.ndarray,
    *,
    origin_au: float,
    symmetric_potential: bool,
) -> tuple[np.ndarray, str, int]:
    degree = property_derivative_degree(spec)
    parity = spec.parity
    if parity == "auto":
        if symmetric_potential:
            nonzero_orders = [
                order for order, value in spec.derivatives.items() if abs(value) > 0.0
            ]
            if not nonzero_orders or all(order % 2 == 0 for order in nonzero_orders):
                parity = "even"
            elif all(order % 2 == 1 for order in nonzero_orders):
                parity = "odd"
            else:
                raise ValueError(
                    f"Property {spec.name!r} mixes even and odd derivatives; use explicit parity=full"
                )
        else:
            parity = "full"

    q = np.asarray(coordinate_au, dtype=float) - float(origin_au)
    values = np.full_like(q, float(spec.value), dtype=float)
    for order, derivative in sorted(spec.derivatives.items()):
        if symmetric_potential and parity == "even" and order % 2 == 1:
            continue
        if symmetric_potential and parity == "odd" and order % 2 == 0:
            continue
        values += float(derivative) * np.power(q, order) / math.factorial(order)
    return values, parity, degree


def periodic_adjust(values: np.ndarray, mode: str) -> np.ndarray:
    adjusted = np.array(values, dtype=float).copy()
    if mode == "none" or len(adjusted) < 2:
        return adjusted
    if mode == "average":
        endpoint = 0.5 * (adjusted[0] + adjusted[-1])
        adjusted[0] = endpoint
        adjusted[-1] = endpoint
    elif mode == "first":
        adjusted[-1] = adjusted[0]
    elif mode == "last":
        adjusted[0] = adjusted[-1]
    else:
        raise ValueError(f"Unsupported periodic endpoint mode: {mode}")
    return adjusted


def interpolate_on_dvr_grid(
    s_au: np.ndarray,
    values: np.ndarray,
    *,
    boundary: str,
    n_grid: int,
    repeat: int,
    endpoint_mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    x = s_au - s_au[0]
    if np.any(np.diff(x) <= 0.0):
        raise ValueError("The mass-weighted path coordinate is not strictly increasing")

    y = np.array(values, dtype=float)
    if boundary == "periodic":
        y = periodic_adjust(y, endpoint_mode)
        period = float(x[-1])
        length = period * repeat
        grid = np.linspace(0.0, length, n_grid, endpoint=False)
        phase = np.mod(grid, period)
        interpolator = PchipInterpolator(x, y)
        return grid, interpolator(phase)

    grid = np.linspace(float(x[0]), float(x[-1]), n_grid, endpoint=True)
    interpolator = PchipInterpolator(x, y)
    return grid, interpolator(grid)


def as_2d_float_array(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim == 1:
        return array[:, None]
    if array.ndim != 2:
        raise ValueError("Expected a one- or two-dimensional numeric array")
    return array


def sorted_unique_average(
    coordinate: np.ndarray, values: np.ndarray, *, atol: float = 1.0e-10
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(coordinate, dtype=float)
    y = as_2d_float_array(values)
    if len(x) != len(y):
        raise ValueError("Coordinate and value arrays must have the same length")
    if len(x) == 0:
        raise ValueError("At least one coordinate point is required")
    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]

    unique_x: list[float] = []
    unique_y: list[np.ndarray] = []
    start = 0
    while start < len(x_sorted):
        stop = start + 1
        while stop < len(x_sorted) and abs(x_sorted[stop] - x_sorted[start]) <= atol:
            stop += 1
        unique_x.append(float(np.mean(x_sorted[start:stop])))
        unique_y.append(np.mean(y_sorted[start:stop], axis=0))
        start = stop

    return np.array(unique_x, dtype=float), np.vstack(unique_y)


def symmetrize_1d_samples(
    coordinate_au: np.ndarray,
    values: np.ndarray,
    *,
    mode: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int | str]]:
    x = np.asarray(coordinate_au, dtype=float)
    y = as_2d_float_array(values)
    if len(x) != len(y):
        raise ValueError("Coordinate and value arrays must have the same length")
    if np.any(np.diff(x) <= 0.0):
        raise ValueError("The one-dimensional path coordinate must be strictly increasing")

    x = x - x[0]
    info: dict[str, float | int | str] = {"path_symmetry": mode}
    if mode == "none":
        x_out, y_out = sorted_unique_average(x, y)
        return x_out, y_out, info

    if mode in {"half-even", "half-even-origin"}:
        if abs(x[0]) > 1.0e-10:
            raise ValueError("Internal error: shifted half path does not start at zero")
        x_full = np.concatenate((-x[:0:-1], x))
        y_full = np.concatenate((y[:0:-1], y), axis=0)
        x_out, y_out = sorted_unique_average(x_full, y_full)
        info["symmetry_origin_au"] = 0.0
        info["symmetry_reference_point"] = "first"
        info["symmetry_positive_extent_au"] = float(np.max(x))
        return x_out, y_out, info

    if mode == "half-even-last":
        x_end = float(x[-1])
        if x_end <= 0.0:
            raise ValueError("--path-symmetry half-even-last needs a nonzero path extent")
        signed = x - x_end
        x_full = np.concatenate((signed, -signed[-2::-1]))
        y_full = np.concatenate((y, y[-2::-1]), axis=0)
        x_out, y_out = sorted_unique_average(x_full, y_full)
        info["symmetry_origin_au"] = x_end
        info["symmetry_reference_point"] = "last"
        info["symmetry_positive_extent_au"] = x_end
        return x_out, y_out, info

    if mode == "center-even":
        center = 0.5 * (x[0] + x[-1])
        signed = x - center
        positive_extent = min(abs(float(signed[0])), abs(float(signed[-1])))
        if positive_extent <= 0.0:
            raise ValueError("--path-symmetry center-even needs a nonzero path extent")
        keep = np.abs(signed) <= positive_extent + 1.0e-10
        radii = np.unique(np.round(np.abs(signed[keep]), 12))
        radii = np.array(sorted(float(value) for value in radii), dtype=float)
        interpolator = PchipInterpolator(signed, y, axis=0)
        y_pos = interpolator(radii)
        y_neg = interpolator(-radii)
        y_half = 0.5 * (y_pos + y_neg)
        x_full = np.concatenate((-radii[:0:-1], radii))
        y_full = np.concatenate((y_half[:0:-1], y_half), axis=0)
        x_out, y_out = sorted_unique_average(x_full, y_full)
        info["symmetry_origin_au"] = float(center)
        info["symmetry_reference_point"] = "center"
        info["symmetry_positive_extent_au"] = float(positive_extent)
        return x_out, y_out, info

    raise ValueError(f"Unsupported path symmetry mode: {mode}")


def endpoint_quadratic_tail(
    x: np.ndarray, potential_cm: np.ndarray, *, side: str, fit_points: int
) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    potential_cm = np.asarray(potential_cm, dtype=float)
    if len(x) != len(potential_cm):
        raise ValueError("Coordinate and potential arrays must have the same length")
    n_fit = max(3, min(int(fit_points), len(x)))
    if side == "left":
        x0 = float(x[0])
        local_x = x[:n_fit] - x0
        local_y = potential_cm[:n_fit]
        outward_sign = -1.0
    elif side == "right":
        x0 = float(x[-1])
        local_x = x[-n_fit:] - x0
        local_y = potential_cm[-n_fit:]
        outward_sign = 1.0
    else:
        raise ValueError(f"Unsupported tail side: {side}")

    scale = float(np.max(np.abs(local_x)))
    if scale <= 0.0:
        return 0.0, 0.0
    u = local_x / scale
    delta_y = local_y - local_y[0 if side == "left" else -1]
    design = np.column_stack((u, 0.5 * u * u))
    coeff, *_ = np.linalg.lstsq(design, delta_y, rcond=None)
    derivative = outward_sign * float(coeff[0] / scale)
    curvature = float(coeff[1] / (scale * scale))
    return max(derivative, 0.0), max(curvature, 0.0)


def canonical_extension_mode(mode: str) -> str:
    aliases = {
        "repulsive-quartic": "repulsive-polynomial",
        "morse-quartic": "morse-polynomial",
    }
    return aliases.get(mode, mode)


def repulsive_polynomial_tail(
    endpoint_value_cm: float,
    slope_cm_au: float,
    curvature_cm_au2: float,
    distance_au: np.ndarray,
    *,
    length_au: float,
    target_cm: float,
    degree: int,
) -> tuple[np.ndarray, float]:
    if degree not in (6, 8):
        raise ValueError("Repulsive tail degree must be 6 or 8")
    smooth_target = (
        endpoint_value_cm
        + slope_cm_au * length_au
        + 0.5 * curvature_cm_au2 * length_au**2
    )
    coefficient = max((target_cm - smooth_target) / length_au**degree, 0.0)
    values = (
        endpoint_value_cm
        + slope_cm_au * distance_au
        + 0.5 * curvature_cm_au2 * distance_au**2
        + coefficient * distance_au**degree
    )
    return values, coefficient


def normal_mode_polynomial_tail(
    endpoint_value_cm: float,
    slope_cm_au: float,
    curvature_cm_au2: float,
    distance_au: np.ndarray,
    *,
    length_au: float,
    target_cm: float,
    degree: int,
) -> tuple[np.ndarray, float]:
    return repulsive_polynomial_tail(
        endpoint_value_cm,
        slope_cm_au,
        curvature_cm_au2,
        distance_au,
        length_au=length_au,
        target_cm=target_cm,
        degree=degree,
    )


def tail_shape_diagnostics(endpoint_value_cm: float, tail_values_cm: np.ndarray) -> tuple[float, float]:
    outward_values = np.concatenate(([endpoint_value_cm], np.asarray(tail_values_cm, dtype=float)))
    first_differences = np.diff(outward_values)
    second_differences = np.diff(outward_values, n=2)
    min_step = float(np.min(first_differences)) if len(first_differences) else 0.0
    min_second = float(np.min(second_differences)) if len(second_differences) else 0.0
    return min_step, min_second


def fit_double_morse(
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
) -> tuple[tuple[float, float, float, float], float]:
    radius, radial_values = sorted_unique_average(np.abs(coordinate_au), potential_cm)
    y = radial_values[:, 0]
    r_max = float(radius[-1])
    y_min = float(np.min(y))
    y_range = max(float(np.ptp(y)), 1.0)
    r_min = float(radius[int(np.argmin(y))])
    a0 = 1.0 / max(r_max - r_min, 1.0)
    initial = np.array([y_range, a0, r_min, y_min], dtype=float)
    lower = np.array([1.0e-8, 1.0e-8, 0.0, y_min - 2.0 * y_range], dtype=float)
    upper = np.array([100.0 * y_range, 10.0 / max(r_max, 1.0), r_max, y_min + y_range], dtype=float)

    def model(params: np.ndarray, r: np.ndarray) -> np.ndarray:
        depth, alpha, minimum, offset = params
        return offset + depth * (1.0 - np.exp(-alpha * (r - minimum))) ** 2

    scale = max(y_range, 1.0)

    def residual(params: np.ndarray) -> np.ndarray:
        return (model(params, radius) - y) / scale

    result = least_squares(
        residual,
        initial,
        bounds=(lower, upper),
        loss="soft_l1",
        max_nfev=5000,
    )
    rms = float(math.sqrt(np.mean((model(result.x, radius) - y) ** 2)))
    return (float(result.x[0]), float(result.x[1]), float(result.x[2]), float(result.x[3])), rms


def double_morse_values(params: tuple[float, float, float, float], radius_au: np.ndarray) -> np.ndarray:
    depth, alpha, minimum, offset = params
    return offset + depth * (1.0 - np.exp(-alpha * (radius_au - minimum))) ** 2


def double_morse_derivatives(params: tuple[float, float, float, float], radius_au: float) -> tuple[float, float]:
    depth, alpha, minimum, _offset = params
    exponential = math.exp(-alpha * (radius_au - minimum))
    first = 2.0 * depth * alpha * exponential * (1.0 - exponential)
    second = 2.0 * depth * alpha**2 * (2.0 * exponential**2 - exponential)
    return float(first), float(second)


def infer_well_type(
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
    *,
    requested: str,
) -> tuple[str, dict[str, float | int | str]]:
    x = np.asarray(coordinate_au, dtype=float)
    y = np.asarray(potential_cm, dtype=float)
    if len(x) != len(y):
        raise ValueError("Coordinate and potential arrays must have the same length")
    y = y - float(np.min(y))
    info: dict[str, float | int | str] = {"well_type_requested": requested}
    if requested != "auto":
        info["well_type_used"] = requested
        info["well_type_reason"] = "user-specified"
        return requested, info

    local_minima = [
        index
        for index in range(1, len(y) - 1)
        if y[index] <= y[index - 1] and y[index] <= y[index + 1]
    ]
    spans_origin = np.min(x) < 0.0 < np.max(x)
    if spans_origin:
        radius, radial_values = sorted_unique_average(np.abs(x), y[:, None])
        radial_y = radial_values[:, 0]
        radial_minimum = int(np.argmin(radial_y))
        if radial_minimum > 0 and radial_y[0] > radial_y[radial_minimum] + 1.0e-6:
            info.update(
                {
                    "well_type_used": "double",
                    "well_type_reason": "symmetric double minimum detected",
                    "well_type_local_minima": int(max(2, len(local_minima))),
                }
            )
            return "double", info

    used = "double" if len(local_minima) >= 2 else "single"
    reason = "multiple local minima detected" if used == "double" else "single minimum detected"
    info.update(
        {
            "well_type_used": used,
            "well_type_reason": reason,
            "well_type_local_minima": int(len(local_minima)),
        }
    )
    return used, info


def double_well_criterion_info(
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
    *,
    fit_points: int,
) -> dict[str, float | int | str]:
    """Estimate the Flanigan-de la Vega symmetric double-well shape criterion."""

    def skipped(reason: str) -> dict[str, float | int | str]:
        return {"double_well_criterion_status": f"skipped: {reason}"}

    x = np.asarray(coordinate_au, dtype=float)
    y_input = np.asarray(potential_cm, dtype=float)
    if len(x) != len(y_input):
        raise ValueError("Coordinate and potential arrays must have the same length")
    if len(x) < 5:
        return skipped("fewer than five core points")
    if np.min(x) >= 0.0 or np.max(x) <= 0.0:
        return skipped("symmetrized path does not span the origin")

    radius, radial_values = sorted_unique_average(np.abs(x), y_input[:, None])
    y = np.asarray(radial_values[:, 0], dtype=float)
    finite = np.isfinite(radius) & np.isfinite(y)
    radius = radius[finite]
    y = y[finite]
    if len(radius) < 5:
        return skipped("fewer than five finite radial points")
    if radius[0] > 1.0e-8 * max(1.0, float(radius[-1])):
        return skipped("no point at the symmetry origin")

    y = y - float(np.min(y))
    origin_energy = float(y[0])
    minimum_index = int(np.argmin(y))
    if minimum_index == 0:
        return skipped("minimum is at the symmetry origin")
    if origin_energy <= float(y[minimum_index]) + 1.0e-8:
        return skipped("barrier at the symmetry origin is not positive")

    n_fit = max(5, min(int(fit_points), len(radius)))
    selected = np.argsort(np.abs(radius - radius[minimum_index]))[:n_fit]
    selected = np.array(sorted(selected), dtype=int)
    if minimum_index == selected[0] or minimum_index == selected[-1]:
        return skipped("minimum is at the edge of the local fit window")

    local_r = radius[selected]
    local_y = y[selected]
    r_guess = float(radius[minimum_index])
    scale = float(np.max(np.abs(local_r - r_guess)))
    if scale <= 0.0:
        return skipped("zero radial span in local fit")

    u = (local_r - r_guess) / scale
    degree = min(4, len(local_r) - 1)
    if degree < 2:
        return skipped("local fit degree below quadratic")

    coeff = np.polyfit(u, local_y, degree)
    polynomial = np.poly1d(coeff)
    lower = float(u[0])
    upper = float(u[-1])
    result = minimize_scalar(lambda value: float(polynomial(value)), bounds=(lower, upper), method="bounded")
    if not result.success:
        return skipped("local polynomial minimization failed")
    u_min = float(result.x)
    edge_tol = 1.0e-5 * max(1.0, upper - lower)
    if u_min <= lower + edge_tol or u_min >= upper - edge_tol:
        return skipped("local polynomial minimum is at the fit-window edge")

    curvature = float(np.polyder(polynomial, 2)(u_min) / (scale * scale))
    minimum_radius = r_guess + u_min * scale
    minimum_energy = float(polynomial(u_min))
    barrier = origin_energy - minimum_energy
    if not np.isfinite(curvature) or curvature <= 0.0:
        return skipped("nonpositive fitted curvature at the minimum")
    if not np.isfinite(minimum_radius) or minimum_radius <= 0.0:
        return skipped("nonpositive fitted minimum radius")
    if not np.isfinite(barrier) or barrier <= 0.0:
        return skipped("nonpositive fitted barrier")

    distance = 2.0 * minimum_radius
    threshold = 32.0 * barrier / (distance * distance)
    criterion = 0.25 * distance * math.sqrt(curvature / (2.0 * barrier))
    recommended = "parabola+Gaussian" if criterion < 1.0 else "double-Morse"

    info: dict[str, float | int | str] = {
        "double_well_criterion_status": "ok",
        "double_well_criterion_reference": "Flanigan-de_la_Vega_1978",
        "double_well_fit_points": int(n_fit),
        "double_well_fit_degree": int(degree),
        "double_well_barrier_cm-1": float(barrier),
        "double_well_minimum_au": float(minimum_radius),
        "double_well_distance_au": float(distance),
        "double_well_curvature_cm_au2": float(curvature),
        "double_well_threshold_curvature_cm_au2": float(threshold),
        "double_well_criterion": float(criterion),
        "double_well_recommended_model": recommended,
    }
    if 0.95 <= criterion <= 1.05:
        info["double_well_criterion_note"] = "near threshold; compare both analytic forms"
    return info


def solve_exponential_tail_rate(
    endpoint_value_cm: float,
    slope_cm_au: float,
    curvature_cm_au2: float,
    *,
    length_au: float,
    target_cm: float,
) -> tuple[float, float]:
    def exp_minus_one_minus_z_scalar(z: float) -> float:
        if abs(z) < 1.0e-5:
            return 0.5 * z * z + z**3 / 6.0 + z**4 / 24.0
        return math.expm1(z) - z

    residual_target = target_cm - endpoint_value_cm - slope_cm_au * length_au
    curvature = max(float(curvature_cm_au2), 1.0e-10)
    if residual_target <= 0.0:
        return 1.0 / max(length_au, 1.0), curvature

    quadratic_residual = 0.5 * curvature * length_au * length_au
    if residual_target <= quadratic_residual:
        curvature = 2.0 * residual_target / (length_au * length_au)
        return 1.0e-8 / max(length_au, 1.0), max(curvature, 1.0e-10)

    def residual(rate: float) -> float:
        z = rate * length_au
        if z > 700.0:
            return math.inf
        return curvature * exp_minus_one_minus_z_scalar(z) / (rate * rate) - residual_target

    low = 1.0e-8 / max(length_au, 1.0)
    high = 1.0 / max(length_au, 1.0)
    while residual(high) < 0.0 and high * length_au < 200.0:
        high *= 2.0
    if residual(high) < 0.0:
        return high, curvature
    result = minimize_scalar(
        lambda value: abs(residual(value)),
        bounds=(low, high),
        method="bounded",
        options={"xatol": 1.0e-12},
    )
    return float(result.x), curvature


def repulsive_exponential_tail(
    endpoint_value_cm: float,
    slope_cm_au: float,
    curvature_cm_au2: float,
    distance_au: np.ndarray,
    *,
    length_au: float,
    target_cm: float,
) -> tuple[np.ndarray, float, float]:
    rate, curvature = solve_exponential_tail_rate(
        endpoint_value_cm,
        slope_cm_au,
        curvature_cm_au2,
        length_au=length_au,
        target_cm=target_cm,
    )
    z = np.clip(rate * distance_au, None, 700.0)
    exponential_shape = np.where(
        np.abs(z) < 1.0e-5,
        0.5 * z * z + z**3 / 6.0 + z**4 / 24.0,
        np.expm1(z) - z,
    )
    values = (
        endpoint_value_cm
        + slope_cm_au * distance_au
        + curvature * exponential_shape / (rate * rate)
    )
    return values, float(rate), float(curvature)


def fit_single_morse(
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
) -> tuple[tuple[float, float, float, float, float], float]:
    x = np.asarray(coordinate_au, dtype=float)
    y = np.asarray(potential_cm, dtype=float)
    y = y - float(np.min(y))
    y_range = max(float(np.ptp(y)), 1.0)
    x_span = max(float(np.ptp(x)), 1.0)
    x_min = float(x[int(np.argmin(y))])

    best_params: tuple[float, float, float, float, float] | None = None
    best_rms = math.inf
    for direction in (1.0, -1.0):
        transformed = direction * x
        minimum_initial = direction * x_min
        initial = np.array([y_range, 1.0 / x_span, minimum_initial, 0.0], dtype=float)
        lower = np.array([1.0e-8, 1.0e-8, float(np.min(transformed)) - x_span, -2.0 * y_range], dtype=float)
        upper = np.array([100.0 * y_range, 20.0 / x_span, float(np.max(transformed)) + x_span, y_range], dtype=float)

        def model(params: np.ndarray, values: np.ndarray) -> np.ndarray:
            depth, alpha, minimum, offset = params
            z = np.clip(-alpha * (values - minimum), -200.0, 200.0)
            return offset + depth * (1.0 - np.exp(z)) ** 2

        def residual(params: np.ndarray) -> np.ndarray:
            return (model(params, transformed) - y) / y_range

        result = least_squares(
            residual,
            initial,
            bounds=(lower, upper),
            loss="soft_l1",
            max_nfev=10000,
        )
        fit = model(result.x, transformed)
        rms = float(math.sqrt(np.mean((fit - y) ** 2)))
        if rms < best_rms:
            best_rms = rms
            best_params = (
                float(result.x[0]),
                float(result.x[1]),
                float(result.x[2]),
                float(result.x[3]),
                float(direction),
            )

    if best_params is None:
        raise ValueError("Single-Morse fit failed")
    return best_params, best_rms


def single_morse_values(params: tuple[float, float, float, float, float], coordinate_au: np.ndarray) -> np.ndarray:
    depth, alpha, minimum, offset, direction = params
    transformed = direction * np.asarray(coordinate_au, dtype=float)
    z = np.clip(-alpha * (transformed - minimum), -200.0, 200.0)
    return offset + depth * (1.0 - np.exp(z)) ** 2


def fit_single_inverse_power(
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
    *,
    extension_length_au: float,
) -> tuple[tuple[float, float, float, float], float]:
    x = np.asarray(coordinate_au, dtype=float)
    y = np.asarray(potential_cm, dtype=float)
    y = y - float(np.min(y))
    y_range = max(float(np.ptp(y)), 1.0)
    span = max(float(np.ptp(x)), 1.0)
    left_extension = max(float(extension_length_au), 0.0)
    shift = float(np.min(x) - left_extension - 0.05 * span - 1.0e-6)
    r = x - shift
    if np.any(r <= 0.0):
        raise ValueError("Internal error: nonpositive radius in inverse-power fit")
    r_min_initial = max(float(r[int(np.argmin(y))]), 1.0e-6)
    initial = np.array([y_range, r_min_initial, 0.0], dtype=float)
    lower = np.array([1.0e-8, 1.0e-8, -2.0 * y_range], dtype=float)
    upper = np.array([100.0 * y_range, float(np.max(r)) + span, y_range], dtype=float)

    def model(params: np.ndarray, radius: np.ndarray) -> np.ndarray:
        depth, minimum_radius, offset = params
        return offset + depth * (1.0 - minimum_radius / radius) ** 2

    def residual(params: np.ndarray) -> np.ndarray:
        return (model(params, r) - y) / y_range

    result = least_squares(
        residual,
        initial,
        bounds=(lower, upper),
        loss="soft_l1",
        max_nfev=10000,
    )
    fit = model(result.x, r)
    rms = float(math.sqrt(np.mean((fit - y) ** 2)))
    params = (float(result.x[0]), float(result.x[1]), float(result.x[2]), shift)
    return params, rms


def single_inverse_power_values(params: tuple[float, float, float, float], coordinate_au: np.ndarray) -> np.ndarray:
    depth, minimum_radius, offset, shift = params
    radius = np.asarray(coordinate_au, dtype=float) - shift
    if np.any(radius <= 0.0):
        raise ValueError("single-inverse-power tail crossed the fitted radial origin")
    return offset + depth * (1.0 - minimum_radius / radius) ** 2


def asymmetric_parabola_gaussian_values(params: tuple[float, float, float, float, float, float], coordinate_au: np.ndarray) -> np.ndarray:
    offset, quadratic, parabola_center, gaussian_height, gaussian_alpha, gaussian_center = params
    x = np.asarray(coordinate_au, dtype=float)
    exponent = -gaussian_alpha * (x - gaussian_center) ** 2
    return offset + quadratic * (x - parabola_center) ** 2 + gaussian_height * np.exp(exponent)


def fit_asymmetric_parabola_gaussian(
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
) -> tuple[tuple[float, float, float, float, float, float], float, float]:
    x = np.asarray(coordinate_au, dtype=float)
    y = np.asarray(potential_cm, dtype=float)
    if len(x) != len(y):
        raise ValueError("Coordinate and potential arrays must have the same length")
    if len(x) < 5:
        raise ValueError("The asymmetric parabola+Gaussian core fit needs at least five points")
    if np.any(np.diff(x) <= 0.0):
        raise ValueError("The asymmetric parabola+Gaussian core fit needs sorted unique coordinates")

    y = y - float(np.min(y))
    span = max(float(np.ptp(x)), 1.0)
    y_range = max(float(np.ptp(y)), 1.0)
    midpoint = 0.5 * (float(x[0]) + float(x[-1]))
    maximum_x = float(x[int(np.argmax(y))])
    minimum_x = float(x[int(np.argmin(y))])
    local_minima = [
        index
        for index in range(1, len(y) - 1)
        if y[index] <= y[index - 1] and y[index] <= y[index + 1]
    ]
    local_maxima = [
        index
        for index in range(1, len(y) - 1)
        if y[index] >= y[index - 1] and y[index] >= y[index + 1]
    ]
    if len(local_minima) >= 2:
        two_minima = sorted(local_minima, key=lambda index: y[index])[:2]
        minima_midpoint = float(np.mean(x[two_minima]))
    else:
        minima_midpoint = midpoint
    if local_maxima:
        local_barrier_x = float(x[max(local_maxima, key=lambda index: y[index])])
    else:
        local_barrier_x = maximum_x
    endpoint_rise = max(float(y[0]), float(y[-1]), y_range)
    quadratic0 = max(endpoint_rise / (0.5 * span) ** 2, 1.0e-8)
    lower = np.array(
        [
            -2.0 * y_range,
            1.0e-12,
            float(x[0]),
            1.0e-12,
            1.0e-12,
            float(x[0]),
        ],
        dtype=float,
    )
    upper = np.array(
        [
            2.0 * y_range,
            1000.0 * y_range / (span * span),
            float(x[-1]),
            1000.0 * y_range,
            1000.0 / (span * span),
            float(x[-1]),
        ],
        dtype=float,
    )

    def residual(raw_params: np.ndarray) -> np.ndarray:
        params = tuple(float(value) for value in raw_params)
        fit = asymmetric_parabola_gaussian_values(params, x)
        return (fit - y) / y_range

    center_guesses = sorted({minimum_x, midpoint, minima_midpoint})
    gaussian_center_guesses = sorted({maximum_x, local_barrier_x, midpoint})
    alpha_guesses = [4.0 / (span * span), 16.0 / (span * span), 40.0 / (span * span)]
    height_guesses = [y_range, 1.5 * y_range]
    offset_guesses = [-0.1 * y_range, 0.0]
    starts: list[np.ndarray] = []
    seen_starts: set[tuple[float, ...]] = set()
    for offset0 in offset_guesses:
        for parabola_center0 in center_guesses:
            for gaussian_height0 in height_guesses:
                for gaussian_alpha0 in alpha_guesses:
                    for gaussian_center0 in gaussian_center_guesses:
                        start = np.array(
                            [
                                offset0,
                                quadratic0,
                                parabola_center0,
                                gaussian_height0,
                                gaussian_alpha0,
                                gaussian_center0,
                            ],
                            dtype=float,
                        )
                        key = tuple(np.round(start, 12))
                        if key not in seen_starts:
                            starts.append(start)
                            seen_starts.add(key)

    best_start = None
    best_score = math.inf
    for initial in starts:
        result = least_squares(
            residual,
            initial,
            bounds=(lower, upper),
            loss="linear",
            max_nfev=300,
        )
        score = float(np.dot(result.fun, result.fun))
        if score < best_score:
            best_score = score
            best_start = result.x
        if score < 1.0e-16:
            break
    if best_start is None:
        raise ValueError("Asymmetric parabola+Gaussian fit failed")
    best_result = least_squares(
        residual,
        best_start,
        bounds=(lower, upper),
        loss="linear",
        max_nfev=10000,
    )
    params = tuple(float(value) for value in best_result.x)
    fit = asymmetric_parabola_gaussian_values(params, x)
    residuals = fit - y
    rms = float(math.sqrt(np.mean(residuals**2)))
    max_abs = float(np.max(np.abs(residuals)))
    return params, rms, max_abs


def apply_1d_core_model(
    coordinate_au: np.ndarray,
    values: np.ndarray,
    args: argparse.Namespace,
    *,
    well_type: str,
) -> tuple[np.ndarray, dict[str, float | int | str]]:
    x = np.asarray(coordinate_au, dtype=float)
    y = as_2d_float_array(values).copy()
    requested = args.core_model
    info: dict[str, float | int | str] = {"core_model_requested": requested}

    use_model = requested
    if requested == "auto":
        sparse = len(x) <= int(args.core_auto_max_points)
        if well_type == "double" and args.path_symmetry == "none" and sparse and len(x) >= 5:
            use_model = "asymmetric-parabola-gaussian"
            info["core_model_auto_reason"] = "sparse asymmetric double-minimum scan"
        else:
            use_model = "sampled"
            info["core_model_auto_reason"] = "sampled profile is safer for this scan"

    info["core_model_used"] = use_model
    if use_model == "sampled":
        return y, info
    if use_model == "asymmetric-parabola-gaussian":
        if well_type != "double":
            raise ValueError("--core-model asymmetric-parabola-gaussian requires --well-type double or auto")
        params, rms, max_abs = fit_asymmetric_parabola_gaussian(x, y[:, 0])
        y[:, 0] = asymmetric_parabola_gaussian_values(params, x)
        info.update(
            {
                "asym_pg_offset_cm-1": params[0],
                "asym_pg_quadratic_cm_au2": params[1],
                "asym_pg_parabola_center_au": params[2],
                "asym_pg_gaussian_height_cm-1": params[3],
                "asym_pg_gaussian_alpha_au-2": params[4],
                "asym_pg_gaussian_center_au": params[5],
                "asym_pg_fit_rms_cm-1": rms,
                "asym_pg_fit_max_abs_cm-1": max_abs,
            }
        )
        return y, info
    raise ValueError(f"Unsupported core model: {requested}")


def extend_1d_samples(
    coordinate_au: np.ndarray,
    values: np.ndarray,
    args: argparse.Namespace,
    normal_mode_tail: NormalModeTailData | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int | str]]:
    x = np.asarray(coordinate_au, dtype=float)
    y = as_2d_float_array(values)
    if len(x) != len(y):
        raise ValueError("Coordinate and value arrays must have the same length")
    extension_mode = canonical_extension_mode(args.potential_extension)
    if extension_mode == "none":
        return x, y, {"potential_extension": "none"}
    if len(x) < 4:
        raise ValueError("Potential extension needs at least four path points")
    if args.extension_length_au <= 0.0:
        raise ValueError("--extension-length-au must be positive when extension is enabled")
    if args.extension_points < 1:
        raise ValueError("--extension-points must be positive when extension is enabled")

    length = float(args.extension_length_au)
    degree = int(args.extension_degree)
    distances = np.linspace(length / args.extension_points, length, args.extension_points)
    info: dict[str, float | int | str] = {
        "potential_extension": extension_mode,
        "extension_degree": degree,
        "extension_length_au": length,
        "extension_points": int(args.extension_points),
        "extension_target_cm-1": float(args.extension_target_cm),
    }

    if extension_mode == "repulsive-polynomial":
        left_slope, left_curvature = endpoint_quadratic_tail(
            x, y[:, 0], side="left", fit_points=args.extension_fit_points
        )
        right_slope, right_curvature = endpoint_quadratic_tail(
            x, y[:, 0], side="right", fit_points=args.extension_fit_points
        )
        left_distances = distances[::-1]
        left_potential, left_coefficient = repulsive_polynomial_tail(
            float(y[0, 0]),
            left_slope,
            left_curvature,
            left_distances,
            length_au=length,
            target_cm=float(args.extension_target_cm),
            degree=degree,
        )
        right_potential, right_coefficient = repulsive_polynomial_tail(
            float(y[-1, 0]),
            right_slope,
            right_curvature,
            distances,
            length_au=length,
            target_cm=float(args.extension_target_cm),
            degree=degree,
        )
        left_x = x[0] - left_distances
        right_x = x[-1] + distances
        info.update(
            {
                "left_tail_slope_cm_au": left_slope,
                "right_tail_slope_cm_au": right_slope,
                "left_tail_curvature_cm_au2": left_curvature,
                "right_tail_curvature_cm_au2": right_curvature,
                "left_tail_polynomial_coefficient": left_coefficient,
                "right_tail_polynomial_coefficient": right_coefficient,
            }
        )
    elif extension_mode == "normal-mode-polynomial":
        if normal_mode_tail is None:
            raise ValueError(
                "--potential-extension normal-mode-polynomial requires --tail-normal-mode-log"
            )
        left_slope_data, left_curvature_data = endpoint_quadratic_tail(
            x, y[:, 0], side="left", fit_points=args.extension_fit_points
        )
        right_slope_data, right_curvature_data = endpoint_quadratic_tail(
            x, y[:, 0], side="right", fit_points=args.extension_fit_points
        )
        left_slope = 0.0 if args.tail_normal_mode_slope == "zero" else left_slope_data
        right_slope = 0.0 if args.tail_normal_mode_slope == "zero" else right_slope_data
        left_curvature = left_curvature_data
        right_curvature = right_curvature_data
        if "left" in normal_mode_tail.sides:
            left_curvature = normal_mode_tail.sides["left"].curvature_cm_au2
        if "right" in normal_mode_tail.sides:
            right_curvature = normal_mode_tail.sides["right"].curvature_cm_au2
        left_distances = distances[::-1]
        left_potential, left_coefficient = normal_mode_polynomial_tail(
            float(y[0, 0]),
            left_slope,
            left_curvature,
            left_distances,
            length_au=length,
            target_cm=float(args.extension_target_cm),
            degree=degree,
        )
        right_potential, right_coefficient = normal_mode_polynomial_tail(
            float(y[-1, 0]),
            right_slope,
            right_curvature,
            distances,
            length_au=length,
            target_cm=float(args.extension_target_cm),
            degree=degree,
        )
        left_x = x[0] - left_distances
        right_x = x[-1] + distances
        info.update(normal_mode_tail.info)
        info.update(
            {
                "normal_mode_tail_slope_mode": args.tail_normal_mode_slope,
                "left_tail_slope_cm_au": left_slope,
                "right_tail_slope_cm_au": right_slope,
                "left_tail_data_slope_cm_au": left_slope_data,
                "right_tail_data_slope_cm_au": right_slope_data,
                "left_tail_curvature_cm_au2": left_curvature,
                "right_tail_curvature_cm_au2": right_curvature,
                "left_tail_data_curvature_cm_au2": left_curvature_data,
                "right_tail_data_curvature_cm_au2": right_curvature_data,
                "left_tail_polynomial_coefficient": left_coefficient,
                "right_tail_polynomial_coefficient": right_coefficient,
            }
        )
    elif extension_mode == "repulsive-exponential":
        left_slope, left_curvature = endpoint_quadratic_tail(
            x, y[:, 0], side="left", fit_points=args.extension_fit_points
        )
        right_slope, right_curvature = endpoint_quadratic_tail(
            x, y[:, 0], side="right", fit_points=args.extension_fit_points
        )
        left_distances = distances[::-1]
        left_potential, left_rate, left_curvature_used = repulsive_exponential_tail(
            float(y[0, 0]),
            left_slope,
            left_curvature,
            left_distances,
            length_au=length,
            target_cm=float(args.extension_target_cm),
        )
        right_potential, right_rate, right_curvature_used = repulsive_exponential_tail(
            float(y[-1, 0]),
            right_slope,
            right_curvature,
            distances,
            length_au=length,
            target_cm=float(args.extension_target_cm),
        )
        left_x = x[0] - left_distances
        right_x = x[-1] + distances
        info.update(
            {
                "left_tail_slope_cm_au": left_slope,
                "right_tail_slope_cm_au": right_slope,
                "left_tail_curvature_cm_au2": left_curvature,
                "right_tail_curvature_cm_au2": right_curvature,
                "left_tail_exponential_rate_au-1": left_rate,
                "right_tail_exponential_rate_au-1": right_rate,
                "left_tail_exponential_curvature_cm_au2": left_curvature_used,
                "right_tail_exponential_curvature_cm_au2": right_curvature_used,
            }
        )
    elif extension_mode == "morse-polynomial":
        if args.path_symmetry == "none":
            raise ValueError("--potential-extension morse-polynomial requires --path-symmetry")
        params, rms = fit_double_morse(x, y[:, 0])
        r_max = max(abs(float(x[0])), abs(float(x[-1])))
        end_value = float(np.mean(y[np.isclose(np.abs(x), r_max), 0]))
        end_fit = float(double_morse_values(params, np.array([r_max]))[0])
        shift = end_value - end_fit
        morse_slope, morse_curvature = double_morse_derivatives(params, r_max)
        data_slope, data_curvature = endpoint_quadratic_tail(
            np.abs(x), y[:, 0], side="right", fit_points=args.extension_fit_points
        )
        slope = max(morse_slope, data_slope, 0.0)
        curvature = max(morse_curvature, data_curvature, 0.0)
        r_new = r_max + distances
        tail_potential, coefficient = repulsive_polynomial_tail(
            end_value,
            slope,
            curvature,
            distances,
            length_au=length,
            target_cm=float(args.extension_target_cm),
            degree=degree,
        )
        left_x = -r_new[::-1]
        right_x = r_new
        left_potential = tail_potential[::-1]
        right_potential = tail_potential
        info.update(
            {
                "morse_depth_cm-1": params[0],
                "morse_alpha_au-1": params[1],
                "morse_minimum_au": params[2],
                "morse_offset_cm-1": params[3],
                "morse_fit_rms_cm-1": rms,
                "morse_endpoint_shift_cm-1": shift,
                "morse_tail_data_slope_cm_au": data_slope,
                "morse_tail_data_curvature_cm_au2": data_curvature,
                "morse_tail_slope_cm_au": slope,
                "morse_tail_curvature_cm_au2": curvature,
                "morse_tail_polynomial_coefficient": coefficient,
            }
        )
    elif extension_mode == "single-morse":
        params, rms = fit_single_morse(x, y[:, 0])
        left_x = x[0] - distances[::-1]
        right_x = x[-1] + distances
        left_potential = single_morse_values(params, left_x)
        right_potential = single_morse_values(params, right_x)
        direction_label = "right-long-range" if params[4] > 0.0 else "left-long-range"
        info.update(
            {
                "single_morse_depth_cm-1": params[0],
                "single_morse_alpha_au-1": params[1],
                "single_morse_minimum_transformed_au": params[2],
                "single_morse_offset_cm-1": params[3],
                "single_morse_direction": direction_label,
                "single_morse_fit_rms_cm-1": rms,
            }
        )
    elif extension_mode == "single-inverse-power":
        params, rms = fit_single_inverse_power(
            x,
            y[:, 0],
            extension_length_au=length,
        )
        left_x = x[0] - distances[::-1]
        right_x = x[-1] + distances
        left_potential = single_inverse_power_values(params, left_x)
        right_potential = single_inverse_power_values(params, right_x)
        info.update(
            {
                "single_inverse_power_depth_cm-1": params[0],
                "single_inverse_power_minimum_radius_au": params[1],
                "single_inverse_power_offset_cm-1": params[2],
                "single_inverse_power_shift_au": params[3],
                "single_inverse_power_fit_rms_cm-1": rms,
            }
        )
    else:
        raise ValueError(f"Unsupported potential extension: {args.potential_extension}")

    left_min_step, left_min_second = tail_shape_diagnostics(float(y[0, 0]), left_potential[::-1])
    right_min_step, right_min_second = tail_shape_diagnostics(float(y[-1, 0]), right_potential)
    info.update(
        {
            "left_tail_min_outward_step_cm-1": left_min_step,
            "right_tail_min_outward_step_cm-1": right_min_step,
            "left_tail_min_second_difference_cm-1": left_min_second,
            "right_tail_min_second_difference_cm-1": right_min_second,
        }
    )

    if y.shape[1] > 1:
        left_properties = np.repeat(y[[0], 1:], args.extension_points, axis=0)
        right_properties = np.repeat(y[[-1], 1:], args.extension_points, axis=0)
        left_values = np.column_stack((left_potential, left_properties))
        right_values = np.column_stack((right_potential, right_properties))
    else:
        left_values = left_potential[:, None]
        right_values = right_potential[:, None]

    x_out = np.concatenate((left_x, x, right_x))
    y_out = np.concatenate((left_values, y, right_values), axis=0)
    x_unique, y_unique = sorted_unique_average(x_out, y_out)
    return x_unique, y_unique, info


def spline_or_pchip_values(
    coordinate_au: np.ndarray,
    values: np.ndarray,
    grid_au: np.ndarray,
    *,
    smoothing: str,
    smoothing_factor: float,
) -> np.ndarray:
    if smoothing == "spline":
        degree = min(3, len(coordinate_au) - 1)
        spline = UnivariateSpline(coordinate_au, values, k=degree, s=smoothing_factor)
        return np.asarray(spline(grid_au), dtype=float)
    if smoothing == "none":
        interpolator = PchipInterpolator(coordinate_au, values)
        return np.asarray(interpolator(grid_au), dtype=float)
    raise ValueError(f"Unsupported potential smoothing mode: {smoothing}")


def add_derivative_property_grids(
    property_grid: dict[str, np.ndarray],
    expectation_properties: dict[str, np.ndarray],
    model_properties: dict[str, np.ndarray],
    *,
    model_x: np.ndarray,
    model_potential_cm: np.ndarray,
    grid_x_for_derivatives: np.ndarray,
    property_derivative_specs: dict[str, PropertyDerivativeSpec],
    symmetric_potential: bool,
    info: dict[str, float | int | str],
) -> None:
    for name, spec in property_derivative_specs.items():
        origin_au = resolve_property_derivative_origin(
            spec.origin,
            model_x,
            model_potential_cm,
        )
        model_values, parity, degree = evaluate_property_derivative_spec(
            spec,
            model_x,
            origin_au=origin_au,
            symmetric_potential=symmetric_potential,
        )
        grid_values, _, _ = evaluate_property_derivative_spec(
            spec,
            grid_x_for_derivatives,
            origin_au=origin_au,
            symmetric_potential=symmetric_potential,
        )
        if symmetric_potential and parity == "odd":
            expectation_grid = np.full_like(grid_values, float(spec.value), dtype=float)
        else:
            expectation_grid = grid_values

        property_grid[name] = grid_values
        expectation_properties[name] = expectation_grid
        model_properties[name] = model_values

        tag = normalized_key(name)
        info[f"property_{tag}_source"] = "derivatives"
        info[f"property_{tag}_degree"] = int(degree)
        info[f"property_{tag}_parity"] = parity
        info[f"property_{tag}_origin_au"] = float(origin_au)
        if symmetric_potential and parity == "odd":
            info[f"property_{tag}_expectation_rule"] = "central_value_for_symmetric_potential"


def build_1d_grid_model(
    s_au: np.ndarray,
    rel_energy_cm: np.ndarray,
    property_samples: dict[str, np.ndarray],
    property_derivative_specs: dict[str, PropertyDerivativeSpec],
    args: argparse.Namespace,
    *,
    n_grid: int,
    repeat: int,
    normal_mode_tail: NormalModeTailData | None = None,
) -> Grid1DModel:
    columns = [np.asarray(rel_energy_cm, dtype=float)]
    property_names = list(property_samples)
    for key in property_names:
        columns.append(np.asarray(property_samples[key], dtype=float))
    samples = np.column_stack(columns)

    model_x, model_values, info = symmetrize_1d_samples(
        s_au,
        samples,
        mode=args.path_symmetry,
    )
    well_type, well_info = infer_well_type(
        model_x,
        model_values[:, 0],
        requested=args.well_type,
    )
    info.update(well_info)
    extension_mode = canonical_extension_mode(args.potential_extension)
    symmetric_potential = args.path_symmetry != "none"

    if args.boundary == "periodic":
        if extension_mode != "none":
            raise ValueError("--potential-extension is supported only for nonperiodic 1D paths")
        if args.potential_smoothing != "none":
            raise ValueError("--potential-smoothing is supported only for nonperiodic 1D paths")
        if args.core_model == "asymmetric-parabola-gaussian":
            raise ValueError("--core-model asymmetric-parabola-gaussian is supported only for nonperiodic 1D paths")
        if args.double_well_criterion == "auto" and well_type == "double":
            info.update(
                double_well_criterion_info(
                    model_x,
                    model_values[:, 0],
                    fit_points=args.double_well_fit_points,
                )
            )
        elif args.double_well_criterion == "auto":
            info["double_well_criterion_status"] = "skipped: well_type is single"
        info.update(
            {
                "core_model_requested": args.core_model,
                "core_model_used": "sampled",
                "potential_extension": "none",
                "potential_smoothing": "none",
                "potential_spline_smoothing": 0.0,
            }
        )
        grid_au, potential_grid_cm = interpolate_on_dvr_grid(
            model_x,
            model_values[:, 0],
            boundary=args.boundary,
            n_grid=n_grid,
            repeat=repeat,
            endpoint_mode=args.periodic_endpoints,
        )
        potential_grid_cm = potential_grid_cm - float(np.min(potential_grid_cm))
        model_values[:, 0] = model_values[:, 0] - float(np.min(model_values[:, 0]))
        property_grid: dict[str, np.ndarray] = {}
        expectation_properties: dict[str, np.ndarray] = {}
        model_properties: dict[str, np.ndarray] = {}
        for i, key in enumerate(property_names, start=1):
            _, grid_values = interpolate_on_dvr_grid(
                model_x,
                model_values[:, i],
                boundary=args.boundary,
                n_grid=n_grid,
                repeat=repeat,
                endpoint_mode=args.periodic_endpoints,
            )
            property_grid[key] = grid_values
            expectation_properties[key] = grid_values
            model_properties[key] = model_values[:, i].copy()
        if property_derivative_specs:
            period = float(model_x[-1] - model_x[0])
            if period <= 0.0:
                raise ValueError("Periodic derivative properties need a nonzero model period")
            grid_x_for_derivatives = np.mod(grid_au, period) + float(model_x[0])
            add_derivative_property_grids(
                property_grid,
                expectation_properties,
                model_properties,
                model_x=model_x,
                model_potential_cm=model_values[:, 0],
                grid_x_for_derivatives=grid_x_for_derivatives,
                property_derivative_specs=property_derivative_specs,
                symmetric_potential=symmetric_potential,
                info=info,
            )
        info.update(
            {
                "model_points": int(len(model_x)),
                "model_min_au": float(model_x[0]),
                "model_max_au": float(model_x[-1]),
                "model_potential_max_cm-1": float(np.max(model_values[:, 0])),
                "grid_potential_max_cm-1": float(np.max(potential_grid_cm)),
            }
        )
        return Grid1DModel(
            grid_au=grid_au,
            potential_cm=potential_grid_cm,
            properties=property_grid,
            expectation_properties=expectation_properties,
            model_points_au=model_x,
            model_potential_cm=model_values[:, 0],
            model_properties=model_properties,
            info=info,
        )

    if extension_mode in {"single-morse", "single-inverse-power"} and well_type != "single":
        raise ValueError(f"--potential-extension {extension_mode} requires a single-minimum core")
    if extension_mode == "morse-polynomial" and well_type != "double":
        raise ValueError("--potential-extension morse-polynomial requires a double-minimum core")
    if args.double_well_criterion == "auto" and well_type == "double":
        info.update(
            double_well_criterion_info(
                model_x,
                model_values[:, 0],
                fit_points=args.double_well_fit_points,
            )
        )
    elif args.double_well_criterion == "auto":
        info["double_well_criterion_status"] = "skipped: well_type is single"
    model_values, core_info = apply_1d_core_model(
        model_x,
        model_values,
        args,
        well_type=well_type,
    )
    info.update(core_info)
    model_x, model_values, extension_info = extend_1d_samples(
        model_x,
        model_values,
        args,
        normal_mode_tail=normal_mode_tail,
    )
    info.update(extension_info)
    model_values[:, 0] = model_values[:, 0] - float(np.min(model_values[:, 0]))

    grid_au = np.linspace(float(model_x[0]), float(model_x[-1]), n_grid, endpoint=True)
    potential_grid_cm = spline_or_pchip_values(
        model_x,
        model_values[:, 0],
        grid_au,
        smoothing=args.potential_smoothing,
        smoothing_factor=args.potential_spline_smoothing,
    )
    potential_grid_cm = potential_grid_cm - float(np.min(potential_grid_cm))

    property_grid = {}
    expectation_properties = {}
    model_properties = {}
    for i, key in enumerate(property_names, start=1):
        interpolator = PchipInterpolator(model_x, model_values[:, i])
        grid_values = np.asarray(interpolator(grid_au), dtype=float)
        property_grid[key] = grid_values
        expectation_properties[key] = grid_values
        model_properties[key] = model_values[:, i].copy()
    if property_derivative_specs:
        add_derivative_property_grids(
            property_grid,
            expectation_properties,
            model_properties,
            model_x=model_x,
            model_potential_cm=model_values[:, 0],
            grid_x_for_derivatives=grid_au,
            property_derivative_specs=property_derivative_specs,
            symmetric_potential=symmetric_potential,
            info=info,
        )

    info.update(
        {
            "potential_smoothing": args.potential_smoothing,
            "potential_spline_smoothing": float(args.potential_spline_smoothing),
            "model_points": int(len(model_x)),
            "model_min_au": float(model_x[0]),
            "model_max_au": float(model_x[-1]),
            "model_potential_max_cm-1": float(np.max(model_values[:, 0])),
            "grid_potential_max_cm-1": float(np.max(potential_grid_cm)),
        }
    )
    return Grid1DModel(
        grid_au=grid_au,
        potential_cm=potential_grid_cm,
        properties=property_grid,
        expectation_properties=expectation_properties,
        model_points_au=model_x,
        model_potential_cm=model_values[:, 0],
        model_properties=model_properties,
        info=info,
    )


def periodic_fourier_kinetic(n_grid: int, length_au: float) -> np.ndarray:
    dx = length_au / n_grid
    wavevectors = 2.0 * math.pi * np.fft.fftfreq(n_grid, d=dx)
    identity = np.eye(n_grid)
    kinetic = np.fft.ifft(
        (0.5 * wavevectors * wavevectors)[:, None] * np.fft.fft(identity, axis=0),
        axis=0,
    )
    return HARTREE_TO_CM * kinetic.real


def sinc_kinetic(n_grid: int, dx_au: float) -> np.ndarray:
    indices = np.arange(n_grid)
    diff = indices[:, None] - indices[None, :]
    kinetic = np.empty((n_grid, n_grid), dtype=float)
    diagonal = diff == 0
    kinetic[diagonal] = math.pi**2 / (6.0 * dx_au * dx_au)
    off = ~diagonal
    kinetic[off] = (
        (-1.0) ** diff[off] / (dx_au * dx_au * diff[off].astype(float) ** 2)
    )
    return HARTREE_TO_CM * kinetic


def sinc_kinetic_with_prefactor(n_grid: int, dx: float, prefactor_cm: float) -> np.ndarray:
    indices = np.arange(n_grid)
    diff = indices[:, None] - indices[None, :]
    kinetic = np.empty((n_grid, n_grid), dtype=float)
    diagonal = diff == 0
    kinetic[diagonal] = math.pi**2 / (6.0 * dx * dx)
    off = ~diagonal
    kinetic[off] = (-1.0) ** diff[off] / (dx * dx * diff[off].astype(float) ** 2)
    return prefactor_cm * kinetic


def anharmonic_taylor_potential(q: np.ndarray, f2: float, f3: float, f4: float) -> np.ndarray:
    return 0.5 * f2 * q * q + (f3 / 6.0) * q**3 + (f4 / 24.0) * q**4


def gaussian_derivative_potential(
    q: np.ndarray,
    f2: float,
    f4: float,
) -> tuple[np.ndarray, dict[str, float | int | str]]:
    if abs(f2) < 1.0e-14:
        potential = (f4 / 24.0) * q**4
        return potential, {"gaussian_status": "f2_zero_quartic_fallback"}
    if abs(f4) < 1.0e-14:
        potential = 0.5 * f2 * q * q
        return potential, {"gaussian_status": "f4_zero_harmonic_fallback"}

    exponent_sign = 1.0 if f4 / f2 >= 0.0 else -1.0
    exponent_scale = abs(f4 / f2) / 6.0
    amplitude = f2 / (2.0 * exponent_sign * exponent_scale)
    argument = np.clip(exponent_sign * exponent_scale * q * q, -100.0, 100.0)
    potential = amplitude * (np.exp(argument) - 1.0)
    info = {
        "gaussian_amplitude_cm-1": float(amplitude),
        "gaussian_exponent_scale": float(exponent_scale),
        "gaussian_exponent_sign": int(exponent_sign),
        "gaussian_status": "matched_f2_f4",
    }
    return potential, info


def morse_derivative_potential(
    q: np.ndarray,
    f2: float,
    f3: float,
) -> tuple[np.ndarray, dict[str, float | int | str]]:
    if f2 <= 0.0:
        raise ValueError("Morse derivative model requires positive second derivative")
    if abs(f3) < 1.0e-14:
        raise ValueError("Morse derivative model requires non-zero third derivative")
    exponent_scale = abs(f3) / (3.0 * f2)
    direction = -1.0 if f3 > 0.0 else 1.0
    dissociation = f2 / (2.0 * exponent_scale * exponent_scale)
    argument = np.clip(-exponent_scale * direction * q, -100.0, 100.0)
    potential = dissociation * (1.0 - np.exp(argument)) ** 2
    predicted_f4 = 7.0 * f2 * exponent_scale * exponent_scale
    return potential, {
        "morse_dissociation_cm-1": float(dissociation),
        "morse_exponent_scale": float(exponent_scale),
        "morse_direction": int(direction),
        "morse_predicted_f4_cm-1": float(predicted_f4),
    }


def handy_morse_derivative_potential(
    q: np.ndarray,
    f2: float,
    f3: float,
    f4: float,
) -> tuple[np.ndarray, dict[str, float | int | str]]:
    if f2 <= 0.0:
        raise ValueError("Handy Morse-like coordinate requires positive second derivative")
    if abs(f3) < 1.0e-14:
        raise ValueError("Handy Morse-like coordinate requires non-zero third derivative")
    alpha = -f3 / (3.0 * f2)
    y = 1.0 - np.exp(np.clip(-alpha * q, -100.0, 100.0))
    c2 = 0.5 * f2 / (alpha * alpha)
    c3 = 0.5 * f2 / (alpha * alpha) + f3 / (6.0 * alpha**3)
    c4 = (
        11.0 * f2 / (24.0 * alpha * alpha)
        + f3 / (4.0 * alpha**3)
        + f4 / (24.0 * alpha**4)
    )
    potential = c2 * y * y + c3 * y**3 + c4 * y**4
    info = {
        "handy_morse_alpha": float(alpha),
        "handy_morse_y2_coeff_cm-1": float(c2),
        "handy_morse_y3_coeff_cm-1": float(c3),
        "handy_morse_y4_coeff_cm-1": float(c4),
        "handy_morse_reference": "Burcl-Carter-Handy CPL 373, Eqs. 6 and 10-13",
    }
    return potential, info


def handy_gaussian_derivative_potential(
    q: np.ndarray,
    f2: float,
    f4: float,
    beta: float,
) -> tuple[np.ndarray, dict[str, float | int | str]]:
    if f2 <= 0.0:
        raise ValueError("Handy Gauss-like coordinate requires positive second derivative")
    if beta <= 0.0:
        raise ValueError("Handy Gauss-like coordinate requires positive beta")
    z2 = 1.0 - np.exp(np.clip(-beta * q * q, -100.0, 100.0))
    c2 = 0.5 * f2 / beta
    c4 = 0.25 * f2 / beta + f4 / (24.0 * beta * beta)
    potential = c2 * z2 + c4 * z2 * z2
    info = {
        "handy_gaussian_beta": float(beta),
        "handy_gaussian_z2_coeff_cm-1": float(c2),
        "handy_gaussian_z4_coeff_cm-1": float(c4),
        "handy_gaussian_reference": "Burcl-Carter-Handy CPL 373, Eqs. 19-22",
    }
    return potential, info


def inverse_power_coefficients(
    r0: float,
    direction: float,
    f2: float,
    f3: float,
) -> tuple[np.ndarray, float]:
    matrix = np.array(
        [
            [
                direction * (-1.0 / r0**2),
                direction * (-2.0 / r0**3),
                direction * (-3.0 / r0**4),
            ],
            [
                2.0 / r0**3,
                6.0 / r0**4,
                12.0 / r0**5,
            ],
            [
                direction**3 * (-6.0 / r0**4),
                direction**3 * (-24.0 / r0**5),
                direction**3 * (-60.0 / r0**6),
            ],
        ],
        dtype=float,
    )
    rhs = np.array([0.0, f2, f3], dtype=float)
    coeff = np.linalg.solve(matrix, rhs)
    predicted_f4 = float(
        24.0 * coeff[0] / r0**5
        + 120.0 * coeff[1] / r0**6
        + 360.0 * coeff[2] / r0**7
    )
    return coeff, predicted_f4


def fit_inverse_power_derivatives(
    f2: float,
    f3: float,
    f4: float,
) -> dict[str, float | int | str]:
    if abs(f3) < 1.0e-14:
        raise ValueError("Inverse-power derivative model requires non-zero third derivative")

    candidates: list[dict[str, float | int | str]] = []
    discriminant = f3 * f3 - f2 * f4
    if abs(f4) < 1.0e-14:
        for direction in (-1.0, 1.0):
            r0 = -3.0 * direction * f2 / f3
            if r0 > 0.0 and math.isfinite(r0):
                candidates.append({"inverse_power_r0": float(r0), "inverse_power_direction": int(direction)})
    elif discriminant >= -1.0e-10 * max(f3 * f3, abs(f2 * f4), 1.0):
        root = math.sqrt(max(discriminant, 0.0))
        for direction in (-1.0, 1.0):
            for sign in (-1.0, 1.0):
                r0 = 6.0 * (-f3 + sign * root) / (direction * f4)
                if r0 > 0.0 and math.isfinite(r0):
                    candidates.append({"inverse_power_r0": float(r0), "inverse_power_direction": int(direction)})
    else:
        raise ValueError(
            "No real inverse-power derivative fit exists because F3^2 - F2*F4 is negative"
        )

    fitted: list[dict[str, float | int | str]] = []
    for candidate in candidates:
        r0 = float(candidate["inverse_power_r0"])
        direction = float(candidate["inverse_power_direction"])
        try:
            coeff, predicted_f4 = inverse_power_coefficients(r0, direction, f2, f3)
        except np.linalg.LinAlgError:
            continue
        if not np.all(np.isfinite(coeff)):
            continue
        c0 = -float(coeff[0] / r0 + coeff[1] / r0**2 + coeff[2] / r0**3)
        fitted.append(
            {
                "inverse_power_r0": float(r0),
                "inverse_power_direction": int(direction),
                "inverse_power_c0": c0,
                "inverse_power_c1": float(coeff[0]),
                "inverse_power_c2": float(coeff[1]),
                "inverse_power_c3": float(coeff[2]),
                "inverse_power_predicted_f4_cm-1": float(predicted_f4),
                "inverse_power_f4_residual_cm-1": float(predicted_f4 - f4),
                "inverse_power_discriminant": float(discriminant),
            }
        )
    if not fitted:
        raise ValueError("Could not fit inverse-power derivative model")
    return max(fitted, key=lambda item: float(item["inverse_power_r0"]))


def inverse_power_derivative_potential(
    q: np.ndarray,
    params: dict[str, float | int | str],
) -> np.ndarray:
    r0 = float(params["inverse_power_r0"])
    direction = float(params["inverse_power_direction"])
    c0 = float(params["inverse_power_c0"])
    c1 = float(params["inverse_power_c1"])
    c2 = float(params["inverse_power_c2"])
    c3 = float(params["inverse_power_c3"])
    r = r0 + direction * q
    if np.any(r <= 0.0):
        raise ValueError("Inverse-power grid crosses the fitted repulsive singularity")
    return c0 + c1 / r + c2 / (r * r) + c3 / (r * r * r)


def inverse_power_grid_bounds(r0: float, direction: float, half_width: float) -> tuple[float, float, str]:
    left = -half_width
    right = half_width
    note = "symmetric"
    singular_q = -r0 / direction
    guard = max(0.02 * r0, 0.02)
    if left < singular_q < 0.0:
        left = singular_q + guard
        note = "left_limited_by_inverse_power_singularity"
    elif 0.0 < singular_q < right:
        right = singular_q - guard
        note = "right_limited_by_inverse_power_singularity"
    if left >= 0.0 or right <= 0.0 or right <= left:
        raise ValueError("Inverse-power singularity leaves no usable grid around Q=0")
    return left, right, note


def choose_anharmonic_well_type(args: argparse.Namespace, derivatives: AnharmonicModeDerivatives) -> str:
    if args.well_type in {"single", "double"}:
        return args.well_type
    return "double" if derivatives.f2_cm < 0.0 else "single"


def build_anharmonic_derivative_potential(
    derivatives: AnharmonicModeDerivatives,
    args: argparse.Namespace,
) -> AnharmonicDerivativePotential:
    f2 = derivatives.f2_cm
    f3 = derivatives.f3_cm
    f4 = derivatives.f4_cm
    half_width = float(args.anharmonic_grid_half_width)
    cubic_threshold = float(args.anharmonic_cubic_threshold)
    resolved_well = choose_anharmonic_well_type(args, derivatives)
    requested_model = args.anharmonic_derivative_model
    model = requested_model
    if model == "auto":
        if resolved_well == "single" and f2 > 0.0 and abs(f3) <= cubic_threshold:
            model = "handy-gaussian"
        elif resolved_well == "single" and f2 > 0.0:
            model = "handy-morse"
        else:
            model = "inverse-power"

    info: dict[str, float | int | str] = {
        "resolved_well_type": resolved_well,
        "reference_point": "TS" if resolved_well == "double" else "minimum",
        "requested_derivative_model": requested_model,
        "derivative_model": model,
    }

    left, right = -half_width, half_width
    if model == "inverse-power":
        try:
            inverse_params = fit_inverse_power_derivatives(f2, f3, f4)
            left, right, grid_note = inverse_power_grid_bounds(
                float(inverse_params["inverse_power_r0"]),
                float(inverse_params["inverse_power_direction"]),
                half_width,
            )
            info.update(inverse_params)
            info["grid_note"] = grid_note
        except ValueError:
            if requested_model == "inverse-power":
                raise
            model = "handy-morse" if f2 > 0.0 and abs(f3) > cubic_threshold else "taylor"
            info["derivative_model"] = model
            info["auto_fallback"] = "inverse_power_failed"

    q = np.linspace(left, right, args.grid)
    taylor = anharmonic_taylor_potential(q, f2, f3, f4)
    if model == "gaussian":
        model_potential, model_info = gaussian_derivative_potential(q, f2, f4)
        info.update(model_info)
    elif model == "handy-gaussian":
        if args.anharmonic_handy_beta is not None:
            beta = float(args.anharmonic_handy_beta)
        elif abs(f4) > 1.0e-14:
            beta = abs(f4) / (6.0 * max(abs(f2), 1.0e-14))
        else:
            beta = 1.0 / max(half_width * half_width, 1.0e-12)
        model_potential, model_info = handy_gaussian_derivative_potential(q, f2, f4, beta)
        info.update(model_info)
    elif model == "morse":
        model_potential, model_info = morse_derivative_potential(q, f2, f3)
        morse_f4 = float(model_info["morse_predicted_f4_cm-1"])
        info.update(model_info)
        info["morse_f4_residual_cm-1"] = morse_f4 - f4
    elif model == "handy-morse":
        model_potential, model_info = handy_morse_derivative_potential(q, f2, f3, f4)
        info.update(model_info)
    elif model == "inverse-power":
        model_potential = inverse_power_derivative_potential(q, info)
    elif model == "taylor":
        model_potential = taylor.copy()
    else:
        raise ValueError(f"Unsupported anharmonic derivative model: {model}")

    wall = np.zeros_like(q)
    if args.anharmonic_wall_height_cm > 0.0:
        wall = float(args.anharmonic_wall_height_cm) * (
            np.abs(q) / max(half_width, 1.0e-12)
        ) ** int(args.anharmonic_wall_degree)
        info["wall_height_cm-1"] = float(args.anharmonic_wall_height_cm)
        info["wall_degree"] = int(args.anharmonic_wall_degree)
    potential = model_potential + wall
    info["grid_left_Q"] = float(left)
    info["grid_right_Q"] = float(right)
    info["grid_points"] = int(args.grid)
    return AnharmonicDerivativePotential(
        q=q,
        potential_cm=potential,
        model_potential_cm=model_potential,
        taylor_potential_cm=taylor,
        info=info,
    )


def solve_anharmonic_derivative_dvr(
    q: np.ndarray,
    potential_cm: np.ndarray,
    kinetic_frequency_cm: float,
    n_levels: int,
) -> tuple[np.ndarray, np.ndarray]:
    dx = float(q[1] - q[0])
    kinetic = sinc_kinetic_with_prefactor(len(q), dx, kinetic_frequency_cm)
    hamiltonian = kinetic + np.diag(potential_cm)
    max_level = min(n_levels, len(q)) - 1
    values, vectors = eigh(hamiltonian, subset_by_index=[0, max_level])
    return values, vectors


def one_mode_vpt2_levels(
    f2_cm: float,
    f3_cm: float,
    f4_cm: float,
    n_levels: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    if f2_cm <= 0.0:
        raise ValueError("One-mode VPT2 comparison requires positive F2")
    omega = float(f2_cm)
    states = np.arange(n_levels, dtype=float)
    quanta = states + 0.5
    harmonic = omega * quanta
    morse_anharmonicity = 5.0 * f3_cm * f3_cm / (48.0 * omega) - f4_cm / 16.0
    constant_shift = f4_cm / 64.0 - 7.0 * f3_cm * f3_cm / (576.0 * omega)
    return harmonic, harmonic - morse_anharmonicity * quanta * quanta + constant_shift, morse_anharmonicity


def build_one_mode_vpt2_comparison(
    derivatives: AnharmonicModeDerivatives,
    potential: AnharmonicDerivativePotential,
    variational_levels_cm: np.ndarray,
) -> AnharmonicVPT2Comparison | None:
    if potential.info.get("resolved_well_type") == "double":
        return None
    if derivatives.f2_cm <= 0.0:
        return None
    harmonic, vpt2, morse_anharmonicity = one_mode_vpt2_levels(
        derivatives.f2_cm,
        derivatives.f3_cm,
        derivatives.f4_cm,
        len(variational_levels_cm),
    )
    model = str(potential.info.get("derivative_model", "unknown"))
    info = {
        "comparison": "standard normal-coordinate VPT2 vs variational levels on the selected analytic potential",
        "formula": "E_n=w(n+1/2)-xe(n+1/2)^2+constant; xe=5F3^2/(48w)-F4/16",
        "variational_model": model,
        "handy_transform": "yes" if model in {"handy-morse", "handy-gaussian"} else "no",
        "omega_cm-1": float(derivatives.f2_cm),
        "vpt2_xe_cm-1": float(morse_anharmonicity),
        "fundamental_harmonic_cm-1": float(harmonic[1] - harmonic[0]) if len(harmonic) > 1 else 0.0,
        "fundamental_vpt2_cm-1": float(vpt2[1] - vpt2[0]) if len(vpt2) > 1 else 0.0,
        "fundamental_variational_cm-1": (
            float(variational_levels_cm[1] - variational_levels_cm[0])
            if len(variational_levels_cm) > 1
            else 0.0
        ),
    }
    if len(variational_levels_cm) > 1:
        info["fundamental_variational_minus_vpt2_cm-1"] = (
            float(variational_levels_cm[1] - variational_levels_cm[0])
            - float(vpt2[1] - vpt2[0])
        )
    return AnharmonicVPT2Comparison(
        harmonic_levels_cm=harmonic,
        vpt2_levels_cm=vpt2,
        variational_levels_cm=np.array(variational_levels_cm, dtype=float),
        info=info,
    )


def solve_dvr(
    grid_au: np.ndarray,
    potential_grid_cm: np.ndarray,
    *,
    boundary: str,
    n_levels: int,
) -> tuple[np.ndarray, np.ndarray]:
    if boundary == "periodic":
        length = float(grid_au[-1] + (grid_au[1] - grid_au[0]) - grid_au[0])
        kinetic = periodic_fourier_kinetic(len(grid_au), length)
    else:
        dx = float(grid_au[1] - grid_au[0])
        kinetic = sinc_kinetic(len(grid_au), dx)

    hamiltonian = kinetic + np.diag(potential_grid_cm)
    max_level = min(n_levels, len(grid_au)) - 1
    values, vectors = eigh(hamiltonian, subset_by_index=[0, max_level])
    return values, vectors


def gaussian_overlap_kinetic(
    centers_au: np.ndarray,
    alpha: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ai = alpha[:, None]
    aj = alpha[None, :]
    xi = centers_au[:, None]
    xj = centers_au[None, :]
    asum = ai + aj
    dx2 = (xi - xj) ** 2
    overlap = np.sqrt(2.0) * (ai * aj) ** 0.25 / np.sqrt(asum) * np.exp(
        -ai * aj * dx2 / asum
    )
    kinetic_hartree = overlap * (ai * aj / asum) * (1.0 - 2.0 * ai * aj * dx2 / asum)
    return overlap, HARTREE_TO_CM * kinetic_hartree


def symmetrize(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + matrix.T)


def gaussian_potential_matrix(
    centers_au: np.ndarray,
    alpha: np.ndarray,
    overlap: np.ndarray,
    coordinate_au: np.ndarray,
    values: np.ndarray,
    *,
    quadrature_order: int,
) -> np.ndarray:
    interpolator = PchipInterpolator(coordinate_au, values)
    nodes, weights = hermgauss(quadrature_order)
    ai = alpha[:, None]
    aj = alpha[None, :]
    xi = centers_au[:, None]
    xj = centers_au[None, :]
    asum = ai + aj
    product_center = (ai * xi + aj * xj) / asum
    lower = float(coordinate_au[0])
    upper = float(coordinate_au[-1])
    points = product_center[..., None] + nodes[None, None, :] / np.sqrt(asum[..., None])
    clipped = np.clip(points, lower, upper)
    sampled = interpolator(clipped.reshape(-1)).reshape(clipped.shape)
    average = np.tensordot(sampled, weights, axes=([-1], [0])) / math.sqrt(math.pi)
    return symmetrize(overlap * average)


def gaussian_basis_values(
    grid_au: np.ndarray,
    centers_au: np.ndarray,
    alpha: np.ndarray,
) -> np.ndarray:
    prefactor = (2.0 * alpha / math.pi) ** 0.25
    diff = grid_au[:, None] - centers_au[None, :]
    return prefactor[None, :] * np.exp(-alpha[None, :] * diff * diff)


def initial_gaussian_basis(
    lower: float,
    upper: float,
    basis_size: int,
    width_scale: float,
) -> tuple[np.ndarray, np.ndarray]:
    nbasis = max(2, int(basis_size))
    centers = np.linspace(lower, upper, nbasis)
    spacing = (upper - lower) / max(nbasis - 1, 1)
    alpha0 = max(float(width_scale), 1.0e-12) / (spacing * spacing)
    alpha = np.full(nbasis, alpha0)
    return centers, alpha


def gaussian_hamiltonian(
    centers_au: np.ndarray,
    alpha: np.ndarray,
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
    *,
    quadrature_order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    overlap, kinetic = gaussian_overlap_kinetic(centers_au, alpha)
    potential = gaussian_potential_matrix(
        centers_au,
        alpha,
        overlap,
        coordinate_au,
        potential_cm,
        quadrature_order=quadrature_order,
    )
    return symmetrize(overlap), symmetrize(kinetic + potential), potential


def solve_stable_generalized(
    hamiltonian: np.ndarray,
    overlap: np.ndarray,
    *,
    n_levels: int,
    overlap_threshold: float,
    condition_limit: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    s_eval, s_vec = eigh(symmetrize(overlap))
    scale = max(float(np.max(np.abs(s_eval))), 1.0)
    condition_limit = max(float(condition_limit), 1.0)
    keep = s_eval > float(overlap_threshold) * scale
    if not np.any(keep):
        raise ValueError("Gaussian overlap matrix is numerically singular")
    while int(np.count_nonzero(keep)) > max(1, int(n_levels)):
        kept = s_eval[keep]
        condition = float(np.max(kept) / np.min(kept))
        if condition <= float(condition_limit):
            break
        first = np.flatnonzero(keep)[0]
        keep[first] = False
    x = s_vec[:, keep] / np.sqrt(s_eval[keep])[None, :]
    h_orth = symmetrize(x.T @ symmetrize(hamiltonian) @ x)
    n_roots = min(int(n_levels), h_orth.shape[0])
    values, y = eigh(h_orth, subset_by_index=[0, n_roots - 1])
    coeff = x @ y
    kept = s_eval[keep]
    info = {
        "basis_size": float(overlap.shape[0]),
        "active_basis_size": float(len(kept)),
        "overlap_condition": float(np.max(kept) / np.min(kept)),
        "overlap_min_eigenvalue": float(np.min(kept)),
        "overlap_max_eigenvalue": float(np.max(kept)),
    }
    return values, coeff, info


def pack_gaussian_parameters(
    centers_au: np.ndarray,
    alpha: np.ndarray,
    mode: str,
) -> np.ndarray:
    parts: list[np.ndarray] = []
    if mode == "centers-widths" and len(centers_au) > 2:
        parts.append(centers_au[1:-1])
    if mode in {"widths", "centers-widths"}:
        parts.append(np.log(alpha))
    if not parts:
        return np.array([], dtype=float)
    return np.concatenate(parts)


def unpack_gaussian_parameters(
    params: np.ndarray,
    centers_au: np.ndarray,
    mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    cursor = 0
    centers = np.array(centers_au, dtype=float)
    if mode == "centers-widths" and len(centers) > 2:
        n_inner = len(centers) - 2
        centers[1:-1] = np.sort(params[cursor : cursor + n_inner])
        cursor += n_inner
    log_alpha = params[cursor : cursor + len(centers)]
    alpha = np.exp(log_alpha)
    return centers, alpha


def optimize_gaussian_basis(
    centers_au: np.ndarray,
    alpha: np.ndarray,
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
    *,
    mode: str,
    objective_levels: int,
    quadrature_order: int,
    overlap_threshold: float,
    condition_limit: float,
    maxiter: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int | str]]:
    if mode == "none" or maxiter <= 0:
        return centers_au, alpha, {
            "gaussian_optimization": "none",
            "gaussian_optimization_success": "not_requested",
            "gaussian_optimization_iterations": 0,
            "gaussian_optimization_evaluations": 0,
        }

    lower = float(coordinate_au[0])
    upper = float(coordinate_au[-1])
    length = upper - lower
    if length <= 0.0:
        raise ValueError("Gaussian basis optimization needs an increasing coordinate")

    params0 = pack_gaussian_parameters(centers_au, alpha, mode)
    if len(params0) == 0:
        return centers_au, alpha, {
            "gaussian_optimization": "none",
            "gaussian_optimization_success": "not_applicable",
            "gaussian_optimization_iterations": 0,
            "gaussian_optimization_evaluations": 0,
        }

    alpha_ref = float(np.median(alpha))
    log_alpha_min = math.log(alpha_ref * 1.0e-2)
    log_alpha_max = math.log(alpha_ref * 1.0e2)
    bounds: list[tuple[float, float]] = []
    if mode == "centers-widths" and len(centers_au) > 2:
        bounds.extend([(lower, upper)] * (len(centers_au) - 2))
    bounds.extend([(log_alpha_min, log_alpha_max)] * len(centers_au))

    min_gap = length / max(200.0 * len(centers_au), 1.0)
    penalty_scale = max(float(np.ptp(potential_cm)), 1.0e3)
    target_levels = max(1, min(int(objective_levels), len(centers_au)))

    def objective(params: np.ndarray) -> float:
        centers, trial_alpha = unpack_gaussian_parameters(params, centers_au, mode)
        gaps = np.diff(centers)
        gap_penalty = float(np.sum(np.maximum(0.0, min_gap - gaps) ** 2))
        try:
            overlap, hamiltonian, _ = gaussian_hamiltonian(
                centers,
                trial_alpha,
                coordinate_au,
                potential_cm,
                quadrature_order=quadrature_order,
            )
            values, _, info = solve_stable_generalized(
                hamiltonian,
                overlap,
                n_levels=target_levels,
                overlap_threshold=overlap_threshold,
                condition_limit=condition_limit,
            )
        except (FloatingPointError, ValueError, np.linalg.LinAlgError):
            return 1.0e30
        if float(values[0]) < float(np.min(potential_cm)) - 1.0e-6:
            return 1.0e30
        inactive = len(centers) - int(info["active_basis_size"])
        linear_penalty = inactive * penalty_scale * penalty_scale
        condition_penalty = max(0.0, float(info["overlap_condition"]) / condition_limit - 1.0)
        return float(
            np.sum(values)
            + penalty_scale * gap_penalty
            + linear_penalty
            + penalty_scale * penalty_scale * condition_penalty
        )

    initial_objective = objective(params0)
    result = minimize(
        objective,
        params0,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": int(maxiter), "ftol": 1.0e-9, "maxls": 30},
    )
    use_optimized = bool(np.isfinite(result.fun) and result.fun < initial_objective)
    best = result.x if use_optimized else params0
    centers, optimized_alpha = unpack_gaussian_parameters(best, centers_au, mode)
    return centers, optimized_alpha, {
        "gaussian_optimization": mode,
        "gaussian_optimization_success": str(bool(result.success)),
        "gaussian_optimization_used": str(use_optimized),
        "gaussian_optimization_initial_objective_cm-1": float(initial_objective),
        "gaussian_optimization_iterations": int(result.nit),
        "gaussian_optimization_evaluations": int(result.nfev),
        "gaussian_optimization_objective_cm-1": float(result.fun),
        "gaussian_optimization_message": str(result.message),
    }


def solve_gaussian_basis(
    grid_au: np.ndarray,
    potential_grid_cm: np.ndarray,
    property_grid: dict[str, np.ndarray],
    *,
    n_levels: int,
    basis_size: int,
    width_scale: float,
    quadrature_order: int,
    optimize_mode: str,
    optimization_levels: int,
    optimization_maxiter: int,
    overlap_threshold: float,
    condition_limit: float,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, float]], dict[str, float | int | str]]:
    if len(grid_au) < 2:
        raise ValueError("Gaussian basis solver needs at least two grid points")
    lower = float(grid_au[0])
    upper = float(grid_au[-1])
    if upper <= lower:
        raise ValueError("Gaussian basis solver needs an increasing coordinate")
    nbasis = max(2, min(int(basis_size), max(2, len(grid_au) - 1)))
    centers, alpha = initial_gaussian_basis(lower, upper, nbasis, width_scale)
    objective_levels = optimization_levels or min(n_levels, nbasis)
    centers, alpha, opt_info = optimize_gaussian_basis(
        centers,
        alpha,
        grid_au,
        potential_grid_cm,
        mode=optimize_mode,
        objective_levels=objective_levels,
        quadrature_order=quadrature_order,
        overlap_threshold=overlap_threshold,
        condition_limit=condition_limit,
        maxiter=optimization_maxiter,
    )
    overlap, hamiltonian, _ = gaussian_hamiltonian(
        centers,
        alpha,
        grid_au,
        potential_grid_cm,
        quadrature_order=quadrature_order,
    )
    values, coeff, eig_info = solve_stable_generalized(
        hamiltonian,
        overlap,
        n_levels=n_levels,
        overlap_threshold=overlap_threshold,
        condition_limit=condition_limit,
    )

    basis_grid = gaussian_basis_values(grid_au, centers, alpha)
    wavefunctions = basis_grid @ coeff
    dx = float(grid_au[1] - grid_au[0])
    for state in range(wavefunctions.shape[1]):
        norm = math.sqrt(float(np.trapezoid(wavefunctions[:, state] ** 2, grid_au)))
        if norm > 0.0:
            wavefunctions[:, state] /= norm
    vectors_for_output = wavefunctions * math.sqrt(dx)

    property_matrices = {
        key: gaussian_potential_matrix(
            centers,
            alpha,
            overlap,
            grid_au,
            values_grid,
            quadrature_order=quadrature_order,
        )
        for key, values_grid in property_grid.items()
    }
    expectations: list[dict[str, float]] = []
    for state in range(coeff.shape[1]):
        c = coeff[:, state]
        row = {key: float(c @ matrix @ c) for key, matrix in property_matrices.items()}
        expectations.append(row)
    info: dict[str, float | int | str] = {
        **opt_info,
        **eig_info,
        "gaussian_width_min": float(np.min(alpha)),
        "gaussian_width_max": float(np.max(alpha)),
        "gaussian_center_min_au": float(np.min(centers)),
        "gaussian_center_max_au": float(np.max(centers)),
    }
    return values, vectors_for_output, expectations, info


def solve_path_hamiltonian(
    grid_au: np.ndarray,
    potential_grid_cm: np.ndarray,
    property_grid: dict[str, np.ndarray],
    args: argparse.Namespace,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, float]], str]:
    solver = resolve_axis_solver(args.boundary, args.solver)
    if solver == "fourier":
        if args.boundary != "periodic":
            raise ValueError("--solver fourier requires --boundary periodic")
        levels, vectors = solve_dvr(
            grid_au,
            potential_grid_cm,
            boundary="periodic",
            n_levels=args.levels,
        )
        return levels, vectors, expectation_values(vectors, property_grid), "fourier"
    if solver == "sinc-dvr":
        levels, vectors = solve_dvr(
            grid_au,
            potential_grid_cm,
            boundary="nonperiodic",
            n_levels=args.levels,
        )
        return levels, vectors, expectation_values(vectors, property_grid), "sinc-dvr"
    if solver == "gaussian":
        if args.boundary == "periodic":
            raise ValueError("--solver gaussian is intended for nonperiodic paths")
        levels, vectors, expectations, info = solve_gaussian_basis(
            grid_au,
            potential_grid_cm,
            property_grid,
            n_levels=args.levels,
            basis_size=args.gaussian_basis_size,
            width_scale=args.gaussian_width_scale,
            quadrature_order=args.gaussian_quadrature,
            optimize_mode=args.gaussian_optimize,
            optimization_levels=args.gaussian_optimization_levels,
            optimization_maxiter=args.gaussian_optimization_maxiter,
            overlap_threshold=args.gaussian_overlap_threshold,
            condition_limit=args.gaussian_condition_limit,
        )
        args.solver_info = info
        return levels, vectors, expectations, "gaussian"
    raise ValueError(f"Unsupported solver: {solver}")


def resolve_axis_solver(boundary: str, solver: str) -> str:
    if solver == "auto":
        return "fourier" if boundary == "periodic" else "gaussian"
    if solver == "fourier" and boundary != "periodic":
        raise ValueError("Fourier solver requires a periodic coordinate")
    if solver == "gaussian" and boundary == "periodic":
        raise ValueError("Gaussian solver is intended for nonperiodic coordinates")
    return solver


def find_csv_column(fieldnames: list[str], requested: str | None, candidates: Iterable[str]) -> str:
    by_norm = {normalized_key(name): name for name in fieldnames}
    if requested:
        key = normalized_key(requested)
        if key not in by_norm:
            raise ValueError(f"Requested CSV column {requested!r} was not found")
        return by_norm[key]
    for candidate in candidates:
        key = normalized_key(candidate)
        if key in by_norm:
            return by_norm[key]
    raise ValueError(f"Could not infer CSV column; tried {', '.join(candidates)}")


def drop_periodic_duplicate_endpoint(
    q1: np.ndarray,
    q2: np.ndarray,
    potential: np.ndarray,
    properties: dict[str, np.ndarray],
    metric: dict[str, np.ndarray] | None = None,
    structure_grid: np.ndarray | None = None,
    *,
    boundary1: str,
    boundary2: str,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    np.ndarray | None,
]:
    metric = metric or {}
    scale = max(float(np.ptp(potential)), 1.0)
    if boundary1 == "periodic" and len(q1) > 2:
        if np.allclose(potential[0, :], potential[-1, :], rtol=1.0e-7, atol=1.0e-8 * scale):
            q1 = q1[:-1]
            potential = potential[:-1, :]
            properties = {key: values[:-1, :] for key, values in properties.items()}
            metric = {key: values[:-1, :] for key, values in metric.items()}
            if structure_grid is not None:
                structure_grid = structure_grid[:-1, :]
    if boundary2 == "periodic" and len(q2) > 2:
        if np.allclose(potential[:, 0], potential[:, -1], rtol=1.0e-7, atol=1.0e-8 * scale):
            q2 = q2[:-1]
            potential = potential[:, :-1]
            properties = {key: values[:, :-1] for key, values in properties.items()}
            metric = {key: values[:, :-1] for key, values in metric.items()}
            if structure_grid is not None:
                structure_grid = structure_grid[:, :-1]
    return q1, q2, potential, properties, metric, structure_grid


def assert_uniform_grid(values: np.ndarray, label: str) -> None:
    if len(values) < 2:
        raise ValueError(f"{label} needs at least two grid points")
    steps = np.diff(values)
    if not np.all(steps > 0.0):
        raise ValueError(f"{label} grid must be strictly increasing")
    if not np.allclose(steps, steps[0], rtol=1.0e-6, atol=1.0e-10):
        raise ValueError(f"{label} grid must be uniform for the current 2D solver")


def validate_metric_grid_2d(metric: dict[str, np.ndarray]) -> None:
    required = set(METRIC_KEYS)
    missing = required.difference(metric)
    if missing:
        raise ValueError(f"Missing 2D metric elements: {', '.join(sorted(missing))}")
    g11 = metric["g11"]
    g22 = metric["g22"]
    g12 = metric["g12"]
    if np.any(~np.isfinite(g11)) or np.any(~np.isfinite(g22)) or np.any(~np.isfinite(g12)):
        raise ValueError("The 2D metric contains non-finite values")
    determinant = g11 * g22 - g12 * g12
    if np.any(g11 <= 0.0) or np.any(g22 <= 0.0) or np.any(determinant <= 0.0):
        raise ValueError("The 2D metric must be positive definite at every grid point")


def read_metric_from_csv(
    rows: list[dict[str, str]],
    fieldnames: list[str],
    args: argparse.Namespace,
    q1_key: str,
    q2_key: str,
    index1: dict[float, int],
    index2: dict[float, int],
    shape: tuple[int, int],
) -> dict[str, np.ndarray]:
    if args.metric_mode != "csv":
        return {}
    g11_key = find_csv_column(fieldnames, args.g11_key, ("g11", "metric_g11", "G11"))
    g22_key = find_csv_column(fieldnames, args.g22_key, ("g22", "metric_g22", "G22"))
    g12_key = find_csv_column(fieldnames, args.g12_key, ("g12", "g21", "metric_g12", "G12"))
    metric = {key: np.full(shape, np.nan) for key in METRIC_KEYS}
    for row in rows:
        x = parse_float(row[q1_key]) * args.q1_scale_au
        y = parse_float(row[q2_key]) * args.q2_scale_au
        i = index1[x]
        j = index2[y]
        metric["g11"][i, j] = parse_float(row[g11_key])
        metric["g22"][i, j] = parse_float(row[g22_key])
        metric["g12"][i, j] = parse_float(row[g12_key])
    for key, values in metric.items():
        if np.isnan(values).any():
            raise ValueError(f"The 2D CSV does not contain a complete grid for metric column {key}")
    validate_metric_grid_2d(metric)
    return metric


def read_metric_structures(args: argparse.Namespace, expected_count: int) -> list[Structure]:
    sources = [args.grid2d_geom_xyz is not None, args.grid2d_geom_log is not None]
    if sum(sources) != 1:
        raise ValueError(
            "--metric-mode geometry needs exactly one of --grid2d-geom-xyz or --grid2d-geom-log"
        )
    if args.grid2d_geom_xyz is not None:
        structures = read_xyz(args.grid2d_geom_xyz.expanduser().resolve())
    else:
        structures = read_gaussian_log(
            args.grid2d_geom_log.expanduser().resolve(),
            selection=args.grid2d_geom_log_selection,
            energy_source=getattr(args, "gaussian_energy", "auto"),
        )
    if len(structures) != expected_count:
        raise ValueError(
            f"Metric geometry source contains {len(structures)} structures, "
            f"but the 2D CSV contains {expected_count} grid rows"
        )
    check_atom_order(structures)
    return structures


def structure_grid_from_rows(
    rows: list[dict[str, str]],
    structures: list[Structure],
    args: argparse.Namespace,
    q1_key: str,
    q2_key: str,
    index1: dict[float, int],
    index2: dict[float, int],
    shape: tuple[int, int],
) -> np.ndarray:
    grid = np.empty(shape, dtype=object)
    filled = np.zeros(shape, dtype=bool)
    for row, structure in zip(rows, structures):
        x = parse_float(row[q1_key]) * args.q1_scale_au
        y = parse_float(row[q2_key]) * args.q2_scale_au
        i = index1[x]
        j = index2[y]
        if filled[i, j]:
            raise ValueError(f"Duplicate geometry for 2D grid point ({x}, {y})")
        grid[i, j] = structure
        filled[i, j] = True
    if not np.all(filled):
        raise ValueError("The geometry source does not cover the complete 2D grid")
    return grid


def aligned_to_center(
    structure: Structure,
    center_coords: np.ndarray,
    masses: np.ndarray,
) -> np.ndarray:
    aligned, _ = mass_weighted_kabsch(structure.coords_angstrom, center_coords, masses)
    return aligned


def finite_difference_weights(offsets: np.ndarray, derivative_order: int = 1) -> np.ndarray:
    offsets = np.asarray(offsets, dtype=float)
    n_points = len(offsets)
    if n_points <= derivative_order:
        raise ValueError("Not enough stencil points for finite differences")
    matrix = np.vstack([offsets**power for power in range(n_points)])
    rhs = np.zeros(n_points, dtype=float)
    rhs[derivative_order] = math.factorial(derivative_order)
    return np.linalg.solve(matrix, rhs)


def finite_difference_stencil(
    grid: np.ndarray,
    index: int,
    *,
    boundary: str,
    requested_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_points = len(grid)
    if n_points < 2:
        raise ValueError("Metric finite differences need at least two points on each axis")
    requested = max(2, int(requested_points))
    if boundary == "periodic":
        if n_points < 3:
            raise ValueError("Periodic metric finite differences need at least three grid points")
        requested = min(requested, n_points if n_points % 2 == 1 else n_points - 1)
        requested = max(3, requested)
        if requested % 2 == 0:
            requested -= 1
        half = requested // 2
        offsets_index = np.arange(-half, half + 1, dtype=int)
        indices = (index + offsets_index) % n_points
        step = float(grid[1] - grid[0])
        offsets = offsets_index.astype(float) * step
        return indices, finite_difference_weights(offsets)

    requested = min(requested, n_points)
    start = index - requested // 2
    start = max(0, min(start, n_points - requested))
    indices = np.arange(start, start + requested, dtype=int)
    offsets = grid[indices] - grid[index]
    return indices, finite_difference_weights(offsets)


def derivative_from_structure_grid(
    structure_grid: np.ndarray,
    masses: np.ndarray,
    q1: np.ndarray,
    q2: np.ndarray,
    i: int,
    j: int,
    axis: int,
    boundary: str,
    stencil_points: int,
) -> np.ndarray:
    center = remove_translation(structure_grid[i, j].coords_angstrom, masses)
    coordinate = q1 if axis == 0 else q2
    center_index = i if axis == 0 else j
    indices, weights = finite_difference_stencil(
        coordinate,
        center_index,
        boundary=boundary,
        requested_points=stencil_points,
    )
    derivative = np.zeros_like(center)
    for index, weight in zip(indices, weights):
        if axis == 0:
            structure = structure_grid[index, j]
            is_center = index == i
        else:
            structure = structure_grid[i, index]
            is_center = index == j
        aligned = center if is_center else aligned_to_center(structure, center, masses)
        derivative += float(weight) * aligned
    return derivative


def metric_from_structure_grid(
    structure_grid: np.ndarray,
    q1: np.ndarray,
    q2: np.ndarray,
    *,
    boundary1: str,
    boundary2: str,
    stencil_points: int,
) -> dict[str, np.ndarray]:
    first = structure_grid[0, 0]
    masses = masses_for_atoms(first.atoms, first.symbols)
    shape = (len(q1), len(q2))
    cov11 = np.empty(shape, dtype=float)
    cov22 = np.empty(shape, dtype=float)
    cov12 = np.empty(shape, dtype=float)
    sqrt_masses = np.sqrt(masses)[:, None]

    for i in range(len(q1)):
        for j in range(len(q2)):
            d1 = derivative_from_structure_grid(
                structure_grid, masses, q1, q2, i, j, 0, boundary1, stencil_points
            )
            d2 = derivative_from_structure_grid(
                structure_grid, masses, q1, q2, i, j, 1, boundary2, stencil_points
            )
            mw1 = d1 * sqrt_masses * MW_ANGSTROM_TO_AU
            mw2 = d2 * sqrt_masses * MW_ANGSTROM_TO_AU
            cov11[i, j] = float(np.sum(mw1 * mw1))
            cov22[i, j] = float(np.sum(mw2 * mw2))
            cov12[i, j] = float(np.sum(mw1 * mw2))

    metric = metric_from_covariant(cov11, cov22, cov12)
    validate_metric_grid_2d(metric)
    return metric


def metric_from_covariant(
    cov11: np.ndarray,
    cov22: np.ndarray,
    cov12: np.ndarray,
) -> dict[str, np.ndarray]:
    determinant = cov11 * cov22 - cov12 * cov12
    if np.any(determinant <= 0.0):
        raise ValueError("The 2D covariant metric is singular")
    return {
        "g11": cov22 / determinant,
        "g22": cov11 / determinant,
        "g12": -cov12 / determinant,
        "cov11": cov11,
        "cov22": cov22,
        "cov12": cov12,
    }


def matrix_field_eigenvalues(a11: np.ndarray, a22: np.ndarray, a12: np.ndarray) -> np.ndarray:
    trace = a11 + a22
    delta = np.sqrt((a11 - a22) ** 2 + 4.0 * a12 * a12)
    return np.stack((0.5 * (trace - delta), 0.5 * (trace + delta)), axis=0)


def project_spd_field(
    a11: np.ndarray,
    a22: np.ndarray,
    a12: np.ndarray,
    *,
    floor: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float | int]]:
    projected11 = np.empty_like(a11, dtype=float)
    projected22 = np.empty_like(a22, dtype=float)
    projected12 = np.empty_like(a12, dtype=float)
    projection_count = 0
    max_relative_shift = 0.0

    for index in np.ndindex(a11.shape):
        matrix = np.array(
            [[float(a11[index]), float(a12[index])], [float(a12[index]), float(a22[index])]],
            dtype=float,
        )
        evals, evecs = np.linalg.eigh(matrix)
        scale = max(float(np.max(np.abs(evals))), 1.0)
        threshold = max(float(floor) * scale, 1.0e-14)
        clipped = np.maximum(evals, threshold)
        if np.any(clipped != evals):
            projection_count += 1
        projected = evecs @ np.diag(clipped) @ evecs.T
        denom = max(float(np.linalg.norm(matrix)), 1.0e-30)
        max_relative_shift = max(
            max_relative_shift,
            float(np.linalg.norm(projected - matrix) / denom),
        )
        projected11[index] = projected[0, 0]
        projected22[index] = projected[1, 1]
        projected12[index] = projected[0, 1]

    return projected11, projected22, projected12, {
        "projection_count": projection_count,
        "projection_max_relative_shift": max_relative_shift,
    }


def axis_fit_basis(
    values: np.ndarray,
    boundary: str,
    order: int,
) -> tuple[np.ndarray, np.ndarray]:
    if boundary == "periodic":
        step = float(values[1] - values[0])
        period = step * len(values)
        theta = 2.0 * math.pi * (values - values[0]) / period
        columns = [np.ones_like(values)]
        degrees = [0]
        for k in range(1, order + 1):
            columns.append(np.cos(k * theta))
            degrees.append(k)
            columns.append(np.sin(k * theta))
            degrees.append(k)
        return np.column_stack(columns), np.array(degrees, dtype=int)

    lower = float(values[0])
    upper = float(values[-1])
    if upper <= lower:
        raise ValueError("Metric smoothing needs increasing coordinate values")
    scaled = 2.0 * (values - lower) / (upper - lower) - 1.0
    basis = np.polynomial.chebyshev.chebvander(scaled, order)
    return basis, np.arange(order + 1, dtype=int)


def metric_fit_design(
    shape: tuple[int, int],
    q1: np.ndarray,
    q2: np.ndarray,
    *,
    boundary1: str,
    boundary2: str,
    order: int,
) -> tuple[np.ndarray, int, int]:
    max_order = max(0, int(order))
    point_count = int(np.prod(shape))

    for trial_order in range(max_order, -1, -1):
        basis1, degree1 = axis_fit_basis(q1, boundary1, trial_order)
        basis2, degree2 = axis_fit_basis(q2, boundary2, trial_order)
        columns = []
        for i, d1 in enumerate(degree1):
            for j, d2 in enumerate(degree2):
                if int(d1 + d2) <= trial_order:
                    columns.append(np.outer(basis1[:, i], basis2[:, j]).reshape(-1))
        if not columns:
            continue
        design = np.column_stack(columns)
        if design.shape[1] <= point_count:
            return design, trial_order, design.shape[1]

    return np.ones((point_count, 1), dtype=float), 0, 1


def fit_surface_with_design(
    values: np.ndarray,
    design: np.ndarray,
) -> np.ndarray:
    coeff, *_ = np.linalg.lstsq(design, values.reshape(-1), rcond=None)
    return (design @ coeff).reshape(values.shape)


def rms_relative_change(original: np.ndarray, fitted: np.ndarray) -> float:
    numerator = math.sqrt(float(np.mean((fitted - original) ** 2)))
    denominator = max(math.sqrt(float(np.mean(original**2))), 1.0e-30)
    return numerator / denominator


def smooth_metric_elements(
    metric: dict[str, np.ndarray],
    q1: np.ndarray,
    q2: np.ndarray,
    args: argparse.Namespace,
    *,
    source: str,
) -> tuple[dict[str, np.ndarray], dict[str, float | int | str]]:
    diagnostics: dict[str, float | int | str] = {
        "metric_mode": source,
        "metric_smoothing": args.metric_smoothing,
        "metric_smooth_order_requested": int(args.metric_smooth_order),
        "metric_spd_projection_floor": float(args.metric_eigenvalue_floor),
    }
    if source == "geometry":
        diagnostics["metric_stencil_points_requested"] = int(args.metric_stencil)

    smoothing = args.metric_smoothing
    if smoothing == "auto":
        smoothing = "fit" if source in {"geometry", "csv"} else "none"
    diagnostics["metric_smoothing_used"] = smoothing

    work = {key: np.array(value, dtype=float).copy() for key, value in metric.items()}
    if smoothing == "fit":
        if source == "geometry" and set(COVARIANT_METRIC_KEYS).issubset(work):
            target_keys = COVARIANT_METRIC_KEYS
            target_type = "covariant"
        else:
            target_keys = METRIC_KEYS
            target_type = "contravariant"
        diagnostics["metric_smoothing_target"] = target_type
        design, used_order, term_count = metric_fit_design(
            work[target_keys[0]].shape,
            q1,
            q2,
            boundary1=args.boundary1,
            boundary2=args.boundary2,
            order=args.metric_smooth_order,
        )
        for key in target_keys:
            original = work[key]
            fitted = fit_surface_with_design(original, design)
            work[key] = fitted
            diagnostics[f"metric_{key}_fit_rms_relative_change"] = rms_relative_change(
                original, fitted
            )
        diagnostics["metric_smooth_order_used"] = int(used_order)
        diagnostics["metric_smooth_terms"] = int(term_count)

    if source == "geometry" and set(COVARIANT_METRIC_KEYS).issubset(work):
        cov11, cov22, cov12, projection = project_spd_field(
            work["cov11"],
            work["cov22"],
            work["cov12"],
            floor=args.metric_eigenvalue_floor,
        )
        diagnostics.update({f"covariant_{key}": value for key, value in projection.items()})
        work = metric_from_covariant(cov11, cov22, cov12)
    else:
        g11, g22, g12, projection = project_spd_field(
            work["g11"],
            work["g22"],
            work["g12"],
            floor=args.metric_eigenvalue_floor,
        )
        diagnostics.update({f"contravariant_{key}": value for key, value in projection.items()})
        work["g11"] = g11
        work["g22"] = g22
        work["g12"] = g12

    validate_metric_grid_2d(work)
    eig = matrix_field_eigenvalues(work["g11"], work["g22"], work["g12"])
    diagnostics["metric_g_min_eigenvalue"] = float(np.min(eig))
    diagnostics["metric_g_max_eigenvalue"] = float(np.max(eig))
    diagnostics["metric_g_condition_max"] = float(np.max(eig[1] / eig[0]))
    diagnostics["metric_g11_min"] = float(np.min(work["g11"]))
    diagnostics["metric_g22_min"] = float(np.min(work["g22"]))
    diagnostics["metric_abs_g12_max"] = float(np.max(np.abs(work["g12"])))
    return work, diagnostics


def read_2d_csv(args: argparse.Namespace) -> Grid2DData:
    path = args.grid2d_csv.expanduser().resolve()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    if not rows:
        raise ValueError(f"No rows found in {path}")

    q1_key = find_csv_column(fieldnames, args.q1_key, ("q1", "coord1", "x", "phi"))
    q2_key = find_csv_column(fieldnames, args.q2_key, ("q2", "coord2", "y", "psi"))
    energy_key = find_csv_column(fieldnames, args.energy_key, ENERGY_KEYS)
    energy_unit = infer_energy_unit(energy_key, args.energy_unit)

    q1_values = sorted({parse_float(row[q1_key]) * args.q1_scale_au for row in rows})
    q2_values = sorted({parse_float(row[q2_key]) * args.q2_scale_au for row in rows})
    q1 = np.array(q1_values, dtype=float)
    q2 = np.array(q2_values, dtype=float)
    index1 = {value: i for i, value in enumerate(q1_values)}
    index2 = {value: i for i, value in enumerate(q2_values)}
    shape = (len(q1), len(q2))

    raw_energy = np.full(shape, np.nan)
    requested_props = args.property or []
    property_keys = [find_csv_column(fieldnames, item, (item,)) for item in requested_props]
    properties = {key: np.full(shape, np.nan) for key in property_keys}
    metric = read_metric_from_csv(rows, fieldnames, args, q1_key, q2_key, index1, index2, shape)
    structure_grid = None
    if args.metric_mode == "geometry":
        metric_structures = read_metric_structures(args, len(rows))
        structure_grid = structure_grid_from_rows(
            rows,
            metric_structures,
            args,
            q1_key,
            q2_key,
            index1,
            index2,
            shape,
        )

    for row in rows:
        x = parse_float(row[q1_key]) * args.q1_scale_au
        y = parse_float(row[q2_key]) * args.q2_scale_au
        i = index1[x]
        j = index2[y]
        raw_energy[i, j] = parse_float(row[energy_key])
        for key in property_keys:
            properties[key][i, j] = parse_float(row[key])

    if np.isnan(raw_energy).any():
        raise ValueError("The 2D CSV does not contain a complete rectangular energy grid")
    for key, values in properties.items():
        if np.isnan(values).any():
            raise ValueError(f"The 2D CSV does not contain a complete grid for property {key!r}")

    potential = energy_values_cm(raw_energy, energy_unit)
    potential = potential - np.nanmin(potential)
    q1, q2, potential, properties, metric, structure_grid = drop_periodic_duplicate_endpoint(
        q1,
        q2,
        potential,
        properties,
        metric,
        structure_grid,
        boundary1=args.boundary1,
        boundary2=args.boundary2,
    )
    assert_uniform_grid(q1, "q1")
    assert_uniform_grid(q2, "q2")
    metric_note = "constant"
    metric_diagnostics: dict[str, float | int | str] = {}
    if args.metric_mode == "csv":
        metric, metric_diagnostics = smooth_metric_elements(
            metric,
            q1,
            q2,
            args,
            source="csv",
        )
        metric_note = (
            "coordinate-dependent metric read from CSV columns; "
            f"smoothing={metric_diagnostics.get('metric_smoothing_used', 'none')}"
        )
    elif args.metric_mode == "geometry":
        if structure_grid is None:
            raise ValueError("Internal error: missing geometry grid for metric calculation")
        metric = metric_from_structure_grid(
            structure_grid,
            q1,
            q2,
            boundary1=args.boundary1,
            boundary2=args.boundary2,
            stencil_points=args.metric_stencil,
        )
        metric, metric_diagnostics = smooth_metric_elements(
            metric,
            q1,
            q2,
            args,
            source="geometry",
        )
        metric_note = (
            "coordinate-dependent metric from Eckart finite differences; "
            f"stencil={args.metric_stencil}; "
            f"smoothing={metric_diagnostics.get('metric_smoothing_used', 'none')}"
        )
    return Grid2DData(
        path=path,
        q1_key=q1_key,
        q2_key=q2_key,
        energy_key=energy_key,
        energy_unit=energy_unit,
        q1=q1,
        q2=q2,
        potential_cm=potential,
        properties=properties,
        metric=metric,
        metric_note=metric_note,
        metric_diagnostics=metric_diagnostics,
    )


def reference_potentials_2d(
    potential: np.ndarray,
    q1: np.ndarray,
    q2: np.ndarray,
    args: argparse.Namespace,
) -> tuple[np.ndarray, np.ndarray, str]:
    mode = args.prune_potential
    if mode == "min":
        return np.min(potential, axis=1), np.min(potential, axis=0), "minimum over the other coordinate"
    if mode == "mean":
        return np.mean(potential, axis=1), np.mean(potential, axis=0), "mean over the other coordinate"
    if mode == "cut":
        q1_ref = float(args.q1_ref) if args.q1_ref is not None else float(q1[np.argmin(np.min(potential, axis=1))])
        q2_ref = float(args.q2_ref) if args.q2_ref is not None else float(q2[np.argmin(np.min(potential, axis=0))])
        i = int(np.argmin(np.abs(q1 - q1_ref)))
        j = int(np.argmin(np.abs(q2 - q2_ref)))
        return potential[:, j], potential[i, :], f"cuts at q1={q1[i]:.8g}, q2={q2[j]:.8g}"
    raise ValueError(f"Unsupported --prune-potential value: {mode}")


def axis_args_for_2d(args: argparse.Namespace, boundary: str, solver: str, levels: int) -> argparse.Namespace:
    return argparse.Namespace(
        boundary=boundary,
        solver=solver,
        levels=levels,
        gaussian_basis_size=args.gaussian_basis_size,
        gaussian_width_scale=args.gaussian_width_scale,
        gaussian_quadrature=args.gaussian_quadrature,
        gaussian_optimize=args.gaussian_optimize,
        gaussian_optimization_levels=args.gaussian_optimization_levels,
        gaussian_optimization_maxiter=args.gaussian_optimization_maxiter,
        gaussian_overlap_threshold=args.gaussian_overlap_threshold,
        gaussian_condition_limit=args.gaussian_condition_limit,
    )


def solve_axis_basis_2d(
    grid: np.ndarray,
    reference_potential: np.ndarray,
    *,
    boundary: str,
    solver: str,
    basis_size: int,
    args: argparse.Namespace,
) -> AxisBasis2D:
    axis_args = axis_args_for_2d(args, boundary, solver, basis_size)
    levels, vectors, _, solver_used = solve_path_hamiltonian(
        grid,
        reference_potential,
        {},
        axis_args,
    )
    vref_matrix = vectors.T @ (reference_potential[:, None] * vectors)
    kinetic = symmetrize(np.diag(levels) - vref_matrix)
    derivative = derivative_matrix_in_selected_basis(vectors, grid, boundary)
    info = getattr(axis_args, "solver_info", {})
    return AxisBasis2D(
        grid=grid,
        boundary=boundary,
        solver=solver_used,
        reference_potential=reference_potential,
        energies=levels,
        vectors=vectors,
        kinetic=kinetic,
        derivative=derivative,
        info=info,
    )


def slice_axis_basis_2d(basis: AxisBasis2D, indices: np.ndarray) -> AxisBasis2D:
    return AxisBasis2D(
        grid=basis.grid,
        boundary=basis.boundary,
        solver=basis.solver,
        reference_potential=basis.reference_potential,
        energies=basis.energies[indices],
        vectors=basis.vectors[:, indices],
        kinetic=basis.kinetic[np.ix_(indices, indices)],
        derivative=basis.derivative[np.ix_(indices, indices)],
        info=basis.info,
    )


def apply_axis_energy_pruning(
    basis: AxisBasis2D,
    *,
    max_size: int,
    min_size: int,
    energy_window_cm: float,
) -> AxisBasis2D:
    max_size = max(1, min(int(max_size), len(basis.energies)))
    min_size = max(1, min(int(min_size), max_size))
    if energy_window_cm <= 0.0:
        keep_count = max_size
    else:
        threshold = float(basis.energies[0] + energy_window_cm)
        keep_count = int(np.count_nonzero(basis.energies <= threshold))
        keep_count = max(min_size, min(max_size, keep_count))
    return slice_axis_basis_2d(basis, np.arange(keep_count, dtype=int))


def derivative_matrix_in_selected_basis(
    vectors: np.ndarray,
    grid: np.ndarray,
    boundary: str,
) -> np.ndarray:
    dx = float(grid[1] - grid[0])
    if boundary == "periodic":
        derivative = (np.roll(vectors, -1, axis=0) - np.roll(vectors, 1, axis=0)) / (2.0 * dx)
    else:
        edge_order = 2 if len(grid) > 2 else 1
        derivative = np.gradient(vectors, grid, axis=0, edge_order=edge_order)
    matrix = vectors.T @ derivative
    return 0.5 * (matrix - matrix.T)


def product_potential_matrix_2d(
    potential: np.ndarray,
    vectors1: np.ndarray,
    vectors2: np.ndarray,
) -> np.ndarray:
    return product_grid_operator_matrix_2d(
        vectors1,
        vectors2,
        potential,
        vectors1,
        vectors2,
    )


def product_grid_operator_matrix_2d(
    left1: np.ndarray,
    left2: np.ndarray,
    grid_values: np.ndarray,
    right1: np.ndarray,
    right2: np.ndarray,
) -> np.ndarray:
    tensor = np.einsum(
        "ia,jb,ij,ic,jd->abcd",
        left1,
        left2,
        grid_values,
        right1,
        right2,
        optimize=True,
    )
    n1 = left1.shape[1]
    n2 = left2.shape[1]
    return symmetrize(tensor.reshape(n1 * n2, n1 * n2))


def basis_derivative_values(
    vectors: np.ndarray,
    grid: np.ndarray,
    boundary: str,
) -> np.ndarray:
    dx = float(grid[1] - grid[0])
    if boundary == "periodic":
        return (np.roll(vectors, -1, axis=0) - np.roll(vectors, 1, axis=0)) / (2.0 * dx)
    edge_order = 2 if len(grid) > 2 else 1
    return np.gradient(vectors, grid, axis=0, edge_order=edge_order)


def validate_constant_metric_2d(g11: float, g22: float, g12: float) -> None:
    if g11 <= 0.0 or g22 <= 0.0:
        raise ValueError("The 2D kinetic metric needs positive g11 and g22")
    determinant = g11 * g22 - g12 * g12
    if determinant <= 0.0:
        raise ValueError("The constant 2D kinetic metric must be positive definite")


def build_variable_metric_hamiltonian_2d(
    potential_cm: np.ndarray,
    basis1: AxisBasis2D,
    basis2: AxisBasis2D,
    metric: dict[str, np.ndarray],
) -> np.ndarray:
    validate_metric_grid_2d(metric)
    v1 = basis1.vectors
    v2 = basis2.vectors
    dv1 = basis_derivative_values(v1, basis1.grid, basis1.boundary)
    dv2 = basis_derivative_values(v2, basis2.grid, basis2.boundary)

    potential = product_potential_matrix_2d(potential_cm, v1, v2)
    kinetic = 0.5 * HARTREE_TO_CM * (
        product_grid_operator_matrix_2d(dv1, v2, metric["g11"], dv1, v2)
        + product_grid_operator_matrix_2d(v1, dv2, metric["g22"], v1, dv2)
        + product_grid_operator_matrix_2d(dv1, v2, metric["g12"], v1, dv2)
        + product_grid_operator_matrix_2d(v1, dv2, metric["g12"], dv1, v2)
    )
    return symmetrize(kinetic + potential)


def build_product_hamiltonian_2d(
    potential_cm: np.ndarray,
    basis1: AxisBasis2D,
    basis2: AxisBasis2D,
    *,
    g11: float,
    g22: float,
    g12: float,
    metric: dict[str, np.ndarray] | None = None,
) -> np.ndarray:
    if metric:
        return build_variable_metric_hamiltonian_2d(potential_cm, basis1, basis2, metric)
    validate_constant_metric_2d(g11, g22, g12)
    n1 = basis1.vectors.shape[1]
    n2 = basis2.vectors.shape[1]
    identity1 = np.eye(n1)
    identity2 = np.eye(n2)
    hamiltonian = (
        g11 * np.kron(basis1.kinetic, identity2)
        + g22 * np.kron(identity1, basis2.kinetic)
        + product_potential_matrix_2d(potential_cm, basis1.vectors, basis2.vectors)
    )
    if abs(g12) > 0.0:
        hamiltonian -= g12 * np.kron(basis1.derivative, basis2.derivative)
    return symmetrize(hamiltonian)


def prune_product_hamiltonian_by_coefficients(
    hamiltonian: np.ndarray,
    *,
    levels: int,
    threshold: float,
    states: int,
) -> tuple[np.ndarray, np.ndarray, str]:
    n_roots = min(int(levels), hamiltonian.shape[0])
    if threshold <= 0.0:
        return hamiltonian, np.arange(hamiltonian.shape[0], dtype=int), "none"
    values, coeff = eigh(hamiltonian, subset_by_index=[0, n_roots - 1])
    n_states = max(1, min(int(states), coeff.shape[1]))
    weights = np.sum(coeff[:, :n_states] ** 2, axis=1)
    keep = weights >= float(threshold)
    min_keep = min(hamiltonian.shape[0], n_roots)
    if int(np.count_nonzero(keep)) < min_keep:
        ranked = np.argsort(weights)[::-1][:min_keep]
        keep[ranked] = True
    indices = np.flatnonzero(keep)
    note = (
        f"coefficient threshold {threshold:g}; states={n_states}; "
        f"kept {len(indices)} of {hamiltonian.shape[0]}"
    )
    return symmetrize(hamiltonian[np.ix_(indices, indices)]), indices, note


def property_expectations_2d(
    coeff: np.ndarray,
    properties: dict[str, np.ndarray],
    basis1: AxisBasis2D,
    basis2: AxisBasis2D,
    product_indices: np.ndarray | None = None,
) -> list[dict[str, float]]:
    property_matrices = {
        key: product_potential_matrix_2d(values, basis1.vectors, basis2.vectors)
        for key, values in properties.items()
    }
    if product_indices is not None:
        property_matrices = {
            key: matrix[np.ix_(product_indices, product_indices)]
            for key, matrix in property_matrices.items()
        }
    rows: list[dict[str, float]] = []
    for state in range(coeff.shape[1]):
        c = coeff[:, state]
        rows.append({key: float(c @ matrix @ c) for key, matrix in property_matrices.items()})
    return rows


def calculate_2d_product_result(args: argparse.Namespace, data: Grid2DData | None = None) -> Result2D:
    if data is None:
        data = read_2d_csv(args)
    solver1 = resolve_axis_solver(args.boundary1, args.solver1)
    solver2 = resolve_axis_solver(args.boundary2, args.solver2)

    v1_ref, v2_ref, reference_note = reference_potentials_2d(
        data.potential_cm,
        data.q1,
        data.q2,
        args,
    )
    raw_basis1 = solve_axis_basis_2d(
        data.q1,
        v1_ref,
        boundary=args.boundary1,
        solver=solver1,
        basis_size=args.basis1,
        args=args,
    )
    raw_basis2 = solve_axis_basis_2d(
        data.q2,
        v2_ref,
        boundary=args.boundary2,
        solver=solver2,
        basis_size=args.basis2,
        args=args,
    )
    basis1 = apply_axis_energy_pruning(
        raw_basis1,
        max_size=args.basis1,
        min_size=args.prune_min_basis1,
        energy_window_cm=args.prune_energy_window1,
    )
    basis2 = apply_axis_energy_pruning(
        raw_basis2,
        max_size=args.basis2,
        min_size=args.prune_min_basis2,
        energy_window_cm=args.prune_energy_window2,
    )

    hamiltonian = build_product_hamiltonian_2d(
        data.potential_cm,
        basis1,
        basis2,
        g11=args.g11,
        g22=args.g22,
        g12=args.g12,
        metric=data.metric,
    )
    initial_size = hamiltonian.shape[0]
    hamiltonian, product_indices, pruning_note = prune_product_hamiltonian_by_coefficients(
        hamiltonian,
        levels=args.levels,
        threshold=args.prune_product_coeff_threshold,
        states=args.prune_product_states,
    )
    n_roots = min(args.levels, hamiltonian.shape[0])
    levels, coeff = eigh(hamiltonian, subset_by_index=[0, n_roots - 1])

    expectations = property_expectations_2d(
        coeff,
        data.properties,
        basis1,
        basis2,
        product_indices=product_indices,
    )
    return Result2D(
        data=data,
        basis1=basis1,
        basis2=basis2,
        reference_note=reference_note,
        levels=levels,
        coeff=coeff,
        expectations=expectations,
        product_basis_initial_size=initial_size,
        product_basis_final_size=hamiltonian.shape[0],
        product_pruning_note=pruning_note,
    )


def ensure_output_dirs(*paths: Path) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def solve_2d_hamiltonian_from_csv(args: argparse.Namespace) -> None:
    result = calculate_2d_product_result(args)
    data = result.data
    ensure_output_dirs(args.outdir, args.figdir)
    prefix = args.prefix
    write_2d_grid(args.outdir / f"{prefix}_2d_grid.csv", data)
    write_levels(args.outdir / f"{prefix}_2d_levels.csv", result.levels)
    if data.properties:
        write_expectations(args.outdir / f"{prefix}_2d_expectations.csv", result.levels, result.expectations)
    write_2d_summary(
        args.outdir / f"{prefix}_2d_summary.txt",
        args,
        data,
        result.basis1,
        result.basis2,
        result.reference_note,
        result.levels,
        product_basis_initial_size=result.product_basis_initial_size,
        product_basis_final_size=result.product_basis_final_size,
        product_pruning_note=result.product_pruning_note,
    )
    plot_2d_potential(
        args.figdir / f"{prefix}_2d_potential.pdf",
        args.figdir / f"{prefix}_2d_potential.png",
        data,
        result.levels,
    )

    print(args.outdir / f"{prefix}_2d_summary.txt")
    print(args.outdir / f"{prefix}_2d_levels.csv")
    if data.properties:
        print(args.outdir / f"{prefix}_2d_expectations.csv")
    print(args.figdir / f"{prefix}_2d_potential.pdf")


def write_2d_grid(path: Path, data: Grid2DData) -> None:
    keys = list(data.properties)
    metric_keys = list(data.metric)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([data.q1_key, data.q2_key, "relative_energy_cm-1", *keys, *metric_keys])
        for i, q1 in enumerate(data.q1):
            for j, q2 in enumerate(data.q2):
                writer.writerow(
                    [
                        q1,
                        q2,
                        data.potential_cm[i, j],
                        *[data.properties[key][i, j] for key in keys],
                        *[data.metric[key][i, j] for key in metric_keys],
                    ]
                )


def write_2d_summary(
    path: Path,
    args: argparse.Namespace,
    data: Grid2DData,
    basis1: AxisBasis2D,
    basis2: AxisBasis2D,
    reference_note: str,
    levels: np.ndarray,
    *,
    product_basis_initial_size: int,
    product_basis_final_size: int,
    product_pruning_note: str,
) -> None:
    lines = [
        "Two-dimensional product-basis Hamiltonian",
        f"Input: {data.path}",
        f"Coordinates: {data.q1_key}, {data.q2_key}",
        f"Grid: {len(data.q1)} x {len(data.q2)}",
        f"Boundaries: {args.boundary1}, {args.boundary2}",
        f"Solvers: {basis1.solver}, {basis2.solver}",
        f"Retained basis: {basis1.vectors.shape[1]} x {basis2.vectors.shape[1]}",
        f"Product basis: {product_basis_final_size} retained from {product_basis_initial_size}",
        f"Reference potentials for pruning: {reference_note}",
        f"Axis energy pruning windows: {args.prune_energy_window1:.10g}, {args.prune_energy_window2:.10g} cm^-1",
        f"Product coefficient pruning: {product_pruning_note}",
        (
            f"Kinetic metric: {data.metric_note}"
            if data.metric
            else f"Kinetic metric: g11={args.g11:.10g}, g22={args.g22:.10g}, g12={args.g12:.10g}"
        ),
        f"Energy key: {data.energy_key}",
        f"Energy unit: {data.energy_unit}",
        f"Properties: {', '.join(data.properties) if data.properties else 'none'}",
        "Planned extensions: 6-membered ring coordinate builders; 7-membered rings in 3D.",
        "",
        "Axis 1 diagnostics:",
        f"  q range: {data.q1[0]:.10g} to {data.q1[-1]:.10g}",
        f"  lowest pruning levels: {', '.join(f'{x:.6f}' for x in basis1.energies[:min(5, len(basis1.energies))])}",
        "Axis 2 diagnostics:",
        f"  q range: {data.q2[0]:.10g} to {data.q2[-1]:.10g}",
        f"  lowest pruning levels: {', '.join(f'{x:.6f}' for x in basis2.energies[:min(5, len(basis2.energies))])}",
        "",
        "Lowest 2D levels:",
    ]
    if data.metric_diagnostics:
        insert_at = lines.index("Axis 1 diagnostics:")
        metric_lines = ["Metric diagnostics:"]
        for key in sorted(data.metric_diagnostics):
            value = data.metric_diagnostics[key]
            if isinstance(value, float):
                metric_lines.append(f"  {key}: {value:.8g}")
            else:
                metric_lines.append(f"  {key}: {value}")
        metric_lines.append("")
        lines[insert_at:insert_at] = metric_lines
    ground = float(levels[0])
    for i, value in enumerate(levels[: min(12, len(levels))]):
        lines.append(f"  {i:3d} {value:14.6f} cm^-1  {value - ground:14.6f} from ground")
    path.write_text("\n".join(lines) + "\n")


def parse_int_csv(text: str, default: list[int]) -> list[int]:
    if not text.strip():
        return default
    values = sorted({int(item.strip()) for item in text.split(",") if item.strip()})
    if not values or min(values) <= 0:
        raise ValueError("Integer list entries must be positive")
    return values


def parse_string_csv(text: str, default: list[str], allowed: set[str]) -> list[str]:
    if not text.strip():
        return default
    values = [item.strip() for item in text.split(",") if item.strip()]
    invalid = [item for item in values if item not in allowed]
    if invalid:
        raise ValueError(f"Unsupported convergence value(s): {', '.join(invalid)}")
    return values


def downsample_indices(n_points: int, boundary: str, stride: int) -> np.ndarray:
    if stride <= 1:
        return np.arange(n_points, dtype=int)
    if boundary == "periodic":
        if n_points % stride != 0:
            raise ValueError(
                f"Periodic grid with {n_points} points cannot use stride {stride}"
            )
        return np.arange(0, n_points, stride, dtype=int)
    if (n_points - 1) % stride != 0:
        raise ValueError(
            f"Nonperiodic grid with {n_points} points cannot use stride {stride}"
        )
    return np.arange(0, n_points, stride, dtype=int)


def downsample_2d_data(data: Grid2DData, args: argparse.Namespace, stride: int) -> Grid2DData:
    i1 = downsample_indices(len(data.q1), args.boundary1, stride)
    i2 = downsample_indices(len(data.q2), args.boundary2, stride)
    properties = {
        key: values[np.ix_(i1, i2)]
        for key, values in data.properties.items()
    }
    metric = {
        key: values[np.ix_(i1, i2)]
        for key, values in data.metric.items()
    }
    return Grid2DData(
        path=data.path,
        q1_key=data.q1_key,
        q2_key=data.q2_key,
        energy_key=data.energy_key,
        energy_unit=data.energy_unit,
        q1=data.q1[i1],
        q2=data.q2[i2],
        potential_cm=data.potential_cm[np.ix_(i1, i2)],
        properties=properties,
        metric=metric,
        metric_note=data.metric_note,
        metric_diagnostics=dict(data.metric_diagnostics),
    )


def solve_2d_convergence_from_csv(args: argparse.Namespace) -> None:
    data = read_2d_csv(args)
    basis1_values = parse_int_csv(
        args.convergence_basis1,
        sorted({max(1, args.basis1 // 2), args.basis1}),
    )
    basis2_values = parse_int_csv(
        args.convergence_basis2,
        sorted({max(1, args.basis2 // 2), args.basis2}),
    )
    prune_values = parse_string_csv(
        args.convergence_prune_potentials,
        [args.prune_potential],
        {"min", "mean", "cut"},
    )
    strides = sorted(parse_int_csv(args.convergence_strides, [1]), reverse=True)
    n_levels = args.convergence_levels or args.levels
    rows: list[dict[str, object]] = []

    for stride in strides:
        for prune in prune_values:
            for basis1 in basis1_values:
                for basis2 in basis2_values:
                    run_args = argparse.Namespace(**vars(args))
                    run_args.basis1 = basis1
                    run_args.basis2 = basis2
                    run_args.prune_potential = prune
                    run_args.levels = n_levels
                    row: dict[str, object] = {
                        "stride": stride,
                        "prune_potential": prune,
                        "basis1_request": basis1,
                        "basis2_request": basis2,
                    }
                    try:
                        sub_data = downsample_2d_data(data, args, stride)
                        result = calculate_2d_product_result(run_args, data=sub_data)
                        row.update(
                            {
                                "status": "ok",
                                "grid1": len(sub_data.q1),
                                "grid2": len(sub_data.q2),
                                "basis1_retained": result.basis1.vectors.shape[1],
                                "basis2_retained": result.basis2.vectors.shape[1],
                                "product_initial": result.product_basis_initial_size,
                                "product_final": result.product_basis_final_size,
                                "message": "",
                                "levels": result.levels[:n_levels],
                            }
                        )
                    except Exception as exc:
                        row.update(
                            {
                                "status": "error",
                                "grid1": "",
                                "grid2": "",
                                "basis1_retained": "",
                                "basis2_retained": "",
                                "product_initial": "",
                                "product_final": "",
                                "message": str(exc),
                                "levels": np.array([], dtype=float),
                            }
                        )
                    rows.append(row)

    ok_rows = [row for row in rows if row["status"] == "ok"]
    reference_levels = ok_rows[-1]["levels"] if ok_rows else np.array([], dtype=float)

    ensure_output_dirs(args.outdir)
    path = args.outdir / f"{args.prefix}_2d_convergence.csv"
    with path.open("w", newline="") as handle:
        fieldnames = [
            "stride",
            "prune_potential",
            "basis1_request",
            "basis2_request",
            "status",
            "grid1",
            "grid2",
            "basis1_retained",
            "basis2_retained",
            "product_initial",
            "product_final",
            *[f"level_{i}_cm-1" for i in range(n_levels)],
            *[f"delta_{i}_cm-1" for i in range(n_levels)],
            "message",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            output = {key: row.get(key, "") for key in fieldnames}
            levels = np.asarray(row["levels"], dtype=float)
            for i in range(n_levels):
                if i < len(levels):
                    output[f"level_{i}_cm-1"] = levels[i]
                    if i < len(reference_levels):
                        output[f"delta_{i}_cm-1"] = levels[i] - reference_levels[i]
            writer.writerow(output)

    summary = args.outdir / f"{args.prefix}_2d_convergence_summary.txt"
    lines = [
        "Two-dimensional convergence scan",
        f"Input: {data.path}",
        f"Rows: {len(rows)}",
        f"Successful rows: {len(ok_rows)}",
        f"Reference row: last successful row in {path.name}",
        f"Basis1 values: {', '.join(map(str, basis1_values))}",
        f"Basis2 values: {', '.join(map(str, basis2_values))}",
        f"Grid strides: {', '.join(map(str, strides))}",
        f"Pruning potentials: {', '.join(prune_values)}",
        "Planned extensions: 6-membered ring coordinate builders; 7-membered rings in 3D.",
    ]
    summary.write_text("\n".join(lines) + "\n")
    print(path)
    print(summary)


def plot_2d_potential(
    path_pdf: Path,
    path_png: Path,
    data: Grid2DData,
    levels: np.ndarray,
) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.0), constrained_layout=True)
    x, y = np.meshgrid(data.q2, data.q1)
    contour = ax.contourf(x, y, data.potential_cm, levels=24, cmap="viridis")
    fig.colorbar(contour, ax=ax, label="Relative energy / cm$^{-1}$")
    vmin = float(np.min(data.potential_cm))
    vmax = float(np.max(data.potential_cm))
    for value in [level for level in levels[: min(6, len(levels))] if vmin <= level <= vmax]:
        ax.contour(x, y, data.potential_cm, levels=[value], colors="white", linewidths=0.6, alpha=0.75)
    ax.set_xlabel(data.q2_key)
    ax.set_ylabel(data.q1_key)
    ax.set_title("Two-dimensional potential")
    fig.savefig(path_pdf)
    fig.savefig(path_png, dpi=240)
    plt.close(fig)


def expectation_values(vectors: np.ndarray, property_grid: dict[str, np.ndarray]) -> list[dict[str, float]]:
    weights = vectors * vectors
    rows: list[dict[str, float]] = []
    for state in range(vectors.shape[1]):
        row = {}
        for key, values in property_grid.items():
            row[key] = float(weights[:, state] @ values)
        rows.append(row)
    return rows


def write_profile(
    path: Path,
    structures: list[Structure],
    s_sqrtamu_angstrom: np.ndarray,
    s_au: np.ndarray,
    step_lengths: np.ndarray,
    angular_residuals: np.ndarray,
    raw_energy: np.ndarray,
    rel_energy_cm: np.ndarray,
    energy_key: str,
    property_keys: list[str],
    cp_ring_indices: list[int] | None = None,
    gic_to_cp_bridges: dict[int, GicToCremerPopleBridge] | None = None,
) -> None:
    gic_keys = sorted(
        {
            key
            for structure in structures
            for key in structure.props
            if key.startswith("gic_")
        }
    )
    cp_keys: list[str] = []
    if cp_ring_indices is not None:
        cp_key_set = set()
        bridge_key_set = set()
        for structure in structures:
            cp_key_set.update(cremer_pople_mode_values(structure.coords_angstrom, cp_ring_indices))
            bridge_key_set.update(
                bridge_gic_to_cremer_pople_values(structure.props, gic_to_cp_bridges or {})
            )
        cp_keys = sorted(cp_key_set) + sorted(bridge_key_set)

    header = [
        "point",
        "s_sqrtamu_angstrom",
        "s_au",
        "step_sqrtamu_angstrom",
        "angular_residual_amu_angstrom2",
        energy_key,
        "relative_energy_cm-1",
    ]
    header.extend(gic_keys)
    header.extend(cp_keys)
    header.extend(property_keys)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for i, structure in enumerate(structures):
            cp_values: dict[str, float] = {}
            if cp_ring_indices is not None:
                cp_values.update(cremer_pople_mode_values(structure.coords_angstrom, cp_ring_indices))
                cp_values.update(
                    bridge_gic_to_cremer_pople_values(structure.props, gic_to_cp_bridges or {})
                )
            writer.writerow(
                [
                    i + 1,
                    s_sqrtamu_angstrom[i],
                    s_au[i],
                    step_lengths[i],
                    angular_residuals[i],
                    raw_energy[i],
                    rel_energy_cm[i],
                    *[structure.props.get(key, "") for key in gic_keys],
                    *[cp_values.get(key, "") for key in cp_keys],
                    *[structure.props.get(key, "") for key in property_keys],
                ]
            )


def write_oriented_xyz(
    path: Path,
    structures: list[Structure],
    oriented: list[np.ndarray],
    s_sqrtamu_angstrom: np.ndarray,
    s_au: np.ndarray,
    rel_energy_cm: np.ndarray,
    cp_ring_indices: list[int] | None = None,
) -> None:
    with path.open("w") as handle:
        for i, (structure, coords) in enumerate(zip(structures, oriented)):
            handle.write(f"{len(structure.atoms)}\n")
            comment = (
                f"point={i + 1} s_sqrtamu_angstrom={s_sqrtamu_angstrom[i]:.12f} "
                f"s_au={s_au[i]:.12f} relative_energy_cm-1={rel_energy_cm[i]:.12f}"
            )
            if cp_ring_indices is not None:
                cp_values = cremer_pople_mode_values(structure.coords_angstrom, cp_ring_indices)
                for key in sorted(cp_values):
                    comment += f" {key}={cp_values[key]:.12f}"
            handle.write(comment + "\n")
            for symbol, (x, y, z) in zip(structure.symbols, coords):
                handle.write(f"{symbol:2s} {x:18.10f} {y:18.10f} {z:18.10f}\n")


def write_grid(
    path: Path,
    grid_au: np.ndarray,
    potential_grid_cm: np.ndarray,
    property_grid: dict[str, np.ndarray],
    vectors: np.ndarray,
    save_states: int,
) -> None:
    dx = grid_au[1] - grid_au[0]
    state_count = min(save_states, vectors.shape[1])
    header = ["grid", "s_au", "s_sqrtamu_angstrom", "V_cm-1", *property_grid.keys()]
    for state in range(state_count):
        header.extend([f"psi_{state}", f"prob_density_{state}"])

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for i, s_value in enumerate(grid_au):
            row: list[float | int] = [
                i,
                s_value,
                s_value / MW_ANGSTROM_TO_AU,
                potential_grid_cm[i],
                *[values[i] for values in property_grid.values()],
            ]
            for state in range(state_count):
                psi = vectors[i, state]
                row.extend([psi, psi * psi / dx])
            writer.writerow(row)


def write_anharmonic_derivative_grid(
    path: Path,
    q: np.ndarray,
    potential_cm: np.ndarray,
    model_potential_cm: np.ndarray,
    taylor_potential_cm: np.ndarray,
    vectors: np.ndarray,
    save_states: int,
) -> None:
    state_count = min(save_states, vectors.shape[1])
    dx = float(q[1] - q[0])
    header = [
        "point",
        "Q_dimensionless",
        "V_cm-1",
        "V_model_cm-1",
        "V_taylor_cm-1",
    ]
    for state in range(state_count):
        header.extend([f"psi_{state}", f"prob_density_{state}"])
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for i, q_value in enumerate(q):
            row: list[float | int] = [
                i,
                q_value,
                potential_cm[i],
                model_potential_cm[i],
                taylor_potential_cm[i],
            ]
            for state in range(state_count):
                psi = vectors[i, state]
                row.extend([psi, psi * psi / dx])
            writer.writerow(row)


def write_1d_model_profile(path: Path, model: Grid1DModel) -> None:
    property_keys = list(model.model_properties)
    header = ["point", "s_au", "s_sqrtamu_angstrom", "V_cm-1", *property_keys]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for i, s_value in enumerate(model.model_points_au):
            writer.writerow(
                [
                    i + 1,
                    s_value,
                    s_value / MW_ANGSTROM_TO_AU,
                    model.model_potential_cm[i],
                    *[model.model_properties[key][i] for key in property_keys],
                ]
            )


def write_levels(path: Path, values: np.ndarray) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["state", "energy_cm-1", "energy_above_ground_cm-1"])
        ground = float(values[0])
        for i, value in enumerate(values):
            writer.writerow([i, value, value - ground])


def write_anharmonic_vpt2_comparison(
    path: Path,
    comparison: AnharmonicVPT2Comparison,
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "state",
                "harmonic_energy_cm-1",
                "vpt2_1d_energy_cm-1",
                "variational_energy_cm-1",
                "harmonic_transition_cm-1",
                "vpt2_1d_transition_cm-1",
                "variational_transition_cm-1",
                "variational_minus_vpt2_transition_cm-1",
            ]
        )
        harmonic_ground = float(comparison.harmonic_levels_cm[0])
        vpt2_ground = float(comparison.vpt2_levels_cm[0])
        variational_ground = float(comparison.variational_levels_cm[0])
        for state, (harmonic, vpt2, variational) in enumerate(
            zip(
                comparison.harmonic_levels_cm,
                comparison.vpt2_levels_cm,
                comparison.variational_levels_cm,
            )
        ):
            harmonic_transition = float(harmonic - harmonic_ground)
            vpt2_transition = float(vpt2 - vpt2_ground)
            variational_transition = float(variational - variational_ground)
            writer.writerow(
                [
                    state,
                    harmonic,
                    vpt2,
                    variational,
                    harmonic_transition,
                    vpt2_transition,
                    variational_transition,
                    variational_transition - vpt2_transition,
                ]
            )


def read_vpt2_property_csv(
    path: Path,
    property_keys: list[str],
    reference_values: dict[str, float],
) -> dict[str, dict[str, float | str]]:
    by_norm = {normalized_key(key): key for key in property_keys}
    rows: dict[str, dict[str, float]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row_number, row in enumerate(reader, start=2):
            name = (
                row.get("property")
                or row.get("name")
                or row.get("key")
                or row.get("label")
                or ""
            ).strip()
            if not name:
                raise ValueError(f"Missing property name in {path} row {row_number}")
            norm = normalized_key(name)
            if norm not in by_norm:
                raise ValueError(f"VPT2 property {name!r} in {path} row {row_number} was not requested")
            key = by_norm[norm]
            reference = reference_values.get(key)
            if reference is None:
                text = row.get("reference") or row.get("reference_value") or row.get("central_value")
                if text is None or not text.strip():
                    raise ValueError(f"Missing reference value for VPT2 property {name!r} in {path} row {row_number}")
                reference = parse_float(text)

            value_text = (
                row.get("vpt2_value")
                or row.get("perturbative_value")
                or row.get("value")
                or row.get("mean")
            )
            delta_text = row.get("vpt2_delta") or row.get("perturbative_delta") or row.get("delta")
            if value_text is not None and value_text.strip():
                vpt2_value = parse_float(value_text)
            elif delta_text is not None and delta_text.strip():
                vpt2_value = reference + parse_float(delta_text)
            else:
                raise ValueError(
                    f"VPT2 property {name!r} in {path} row {row_number} needs vpt2_value/value or vpt2_delta/delta"
                )

            item: dict[str, float | str] = {
                "reference_value": float(reference),
                "vpt2_1d_value": float(vpt2_value),
                "vpt2_1d_source": "csv",
            }
            total_value_text = row.get("total_perturbative_value") or row.get("total_value")
            total_delta_text = row.get("total_perturbative_delta") or row.get("total_delta")
            if total_value_text is not None and total_value_text.strip():
                item["total_perturbative_value"] = parse_float(total_value_text)
            if total_delta_text is not None and total_delta_text.strip():
                item["total_perturbative_delta"] = parse_float(total_delta_text)
            rows[key] = item
    return rows


def local_taylor_derivatives(
    coordinate_au: np.ndarray,
    values: np.ndarray,
    *,
    center_index: int,
    degree: int,
    fit_points: int,
) -> tuple[dict[int, float], dict[str, float | int | str]]:
    x = np.asarray(coordinate_au, dtype=float)
    y = np.asarray(values, dtype=float)
    if len(x) != len(y):
        raise ValueError("Taylor fit coordinate and values have different lengths")
    if len(x) < 3:
        raise ValueError("Taylor fit needs at least three points")
    center = float(x[center_index])
    order = np.argsort(np.abs(x - center), kind="stable")
    n_points = max(3, min(int(fit_points), len(x)))
    indices = np.sort(order[:n_points])
    q = x[indices] - center
    scale = max(float(np.max(np.abs(q))), 1.0e-8)
    fit_degree = max(1, min(int(degree), len(indices) - 1))
    t = q / scale
    coeff = np.polynomial.polynomial.polyfit(t, y[indices], fit_degree)
    derivatives: dict[int, float] = {}
    for n, value in enumerate(coeff):
        derivatives[n] = float(value) * math.factorial(n) / (scale**n)
    return derivatives, {
        "fit_center_au": center,
        "fit_degree": int(fit_degree),
        "fit_points": int(len(indices)),
        "fit_half_width_au": float(np.max(np.abs(q))),
    }


def harmonic_coordinate_matrix(size: int, omega_au: float) -> np.ndarray:
    sigma = math.sqrt(1.0 / (2.0 * omega_au))
    matrix = np.zeros((size, size), dtype=float)
    for n in range(size - 1):
        value = sigma * math.sqrt(n + 1.0)
        matrix[n, n + 1] = value
        matrix[n + 1, n] = value
    return matrix


def one_dimensional_perturbative_property_value(
    potential_derivatives: dict[int, float],
    property_derivatives: dict[int, float],
    *,
    basis_size: int,
) -> tuple[float, dict[str, float | int | str]]:
    f2 = float(potential_derivatives.get(2, 0.0))
    if f2 <= 0.0:
        raise ValueError("positive second derivative required for local property VPT2")
    omega_au = math.sqrt(f2 / HARTREE_TO_CM)
    omega_cm = HARTREE_TO_CM * omega_au
    size = max(8, int(basis_size))
    q = harmonic_coordinate_matrix(size, omega_au)
    powers = {0: np.eye(size), 1: q}
    for order in range(2, 5):
        powers[order] = powers[order - 1] @ q

    perturbation = np.zeros((size, size), dtype=float)
    for order in (3, 4):
        derivative = float(potential_derivatives.get(order, 0.0))
        if derivative:
            perturbation += derivative * powers[order] / math.factorial(order)

    prop_matrix = np.zeros((size, size), dtype=float)
    for order in range(0, 5):
        derivative = float(property_derivatives.get(order, 0.0))
        if derivative:
            prop_matrix += derivative * powers[order] / math.factorial(order)

    energies = omega_cm * (np.arange(size, dtype=float) + 0.5)
    coeff = np.zeros(size, dtype=float)
    for n in range(1, size):
        coeff[n] = perturbation[n, 0] / (energies[0] - energies[n])
    value = float(prop_matrix[0, 0] + 2.0 * np.dot(coeff[1:], prop_matrix[1:, 0]))
    return value, {
        "omega_cm-1": float(omega_cm),
        "omega_au": float(omega_au),
        "basis_size": int(size),
        "method": "harmonic_first_order_wavefunction_local_taylor",
    }


def automatic_property_vpt2(
    coordinate_au: np.ndarray,
    potential_cm: np.ndarray,
    property_values: dict[str, np.ndarray],
    reference_values: dict[str, float],
    *,
    fit_points: int,
    degree: int,
    basis_size: int,
) -> PropertyVPT2Result:
    comparisons: dict[str, dict[str, float | str]] = {}
    info: dict[str, float | int | str] = {}
    if not property_values:
        return PropertyVPT2Result(comparisons, info)
    x = np.asarray(coordinate_au, dtype=float)
    v = np.asarray(potential_cm, dtype=float)
    if len(x) < 3:
        return PropertyVPT2Result(comparisons, {"property_vpt2_status": "skipped: fewer than three points"})
    center_index = int(np.argmin(v))
    try:
        potential_derivatives, fit_info = local_taylor_derivatives(
            x,
            v,
            center_index=center_index,
            degree=max(2, min(int(degree), 4)),
            fit_points=fit_points,
        )
        f2 = float(potential_derivatives.get(2, 0.0))
        f3 = float(potential_derivatives.get(3, 0.0))
        f4 = float(potential_derivatives.get(4, 0.0))
        info.update({f"property_vpt2_{key}": value for key, value in fit_info.items()})
        info["property_vpt2_f2_cm_au2"] = f2
        info["property_vpt2_f3_cm_au3"] = f3
        info["property_vpt2_f4_cm_au4"] = f4
        if f2 <= 0.0:
            raise ValueError("non-positive local curvature")
    except Exception as exc:
        info["property_vpt2_status"] = f"unavailable: {exc}"
        for key in property_values:
            comparisons[key] = {
                "reference_value": float(reference_values.get(key, np.nan)),
                "vpt2_1d_source": "unavailable",
                "status": str(exc),
            }
        return PropertyVPT2Result(comparisons, info)

    for key, values in property_values.items():
        reference = float(reference_values.get(key, values[center_index]))
        try:
            property_derivatives, _ = local_taylor_derivatives(
                x,
                np.asarray(values, dtype=float),
                center_index=center_index,
                degree=max(0, min(int(degree), 4)),
                fit_points=fit_points,
            )
            property_derivatives[0] = reference
            vpt2_value, method_info = one_dimensional_perturbative_property_value(
                potential_derivatives,
                property_derivatives,
                basis_size=basis_size,
            )
            comparisons[key] = {
                "reference_value": reference,
                "vpt2_1d_value": float(vpt2_value),
                "vpt2_1d_source": "local_taylor",
                "status": "ok",
            }
            info[f"property_{normalized_key(key)}_vpt2_source"] = "local_taylor"
            for method_key, method_value in method_info.items():
                info[f"property_vpt2_{method_key}"] = method_value
        except Exception as exc:
            comparisons[key] = {
                "reference_value": reference,
                "vpt2_1d_source": "unavailable",
                "status": str(exc),
            }
            info[f"property_{normalized_key(key)}_vpt2_source"] = "unavailable"
    info.setdefault("property_vpt2_status", "ok")
    return PropertyVPT2Result(comparisons, info)


def merge_property_vpt2_comparisons(
    automatic: dict[str, dict[str, float | str]],
    external: dict[str, dict[str, float]],
    property_keys: list[str],
    reference_values: dict[str, float],
) -> dict[str, dict[str, float | str]]:
    merged: dict[str, dict[str, float | str]] = {}
    for key in property_keys:
        if key in automatic:
            merged[key] = dict(automatic[key])
        else:
            merged[key] = {
                "reference_value": float(reference_values.get(key, np.nan)),
                "vpt2_1d_source": "unavailable",
                "status": "no automatic perturbative value",
            }
    for key, item in external.items():
        merged[key] = {**merged.get(key, {}), **item, "vpt2_1d_source": "csv", "status": "ok"}
    return merged


def write_vpt2_property_correction(
    path: Path,
    comparisons: dict[str, dict[str, float | str]],
    variational_expectations: dict[str, float],
) -> None:
    header = [
        "property",
        "reference_value",
        "variational_ground_value",
        "vpt2_1d_value",
        "variational_delta",
        "vpt2_1d_delta",
        "variational_minus_vpt2_1d",
        "vpt2_1d_source",
        "status",
        "total_perturbative_value",
        "corrected_total_value",
        "total_perturbative_delta",
        "corrected_total_delta",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for key in comparisons:
            item = comparisons[key]
            reference = float(item.get("reference_value", np.nan))
            variational = variational_expectations[key]
            vpt2_value = item.get("vpt2_1d_value")
            vpt2_float = float(vpt2_value) if vpt2_value is not None else None
            correction = variational - vpt2_float if vpt2_float is not None else None
            total_value = item.get("total_perturbative_value")
            total_delta = item.get("total_perturbative_delta")
            total_value_float = float(total_value) if total_value is not None and total_value != "" else None
            total_delta_float = float(total_delta) if total_delta is not None and total_delta != "" else None
            writer.writerow(
                [
                    key,
                    reference,
                    variational,
                    vpt2_float if vpt2_float is not None else "",
                    variational - reference,
                    vpt2_float - reference if vpt2_float is not None else "",
                    correction if correction is not None else "",
                    item.get("vpt2_1d_source", ""),
                    item.get("status", ""),
                    total_value_float if total_value_float is not None else "",
                    total_value_float + correction if total_value_float is not None and correction is not None else "",
                    total_delta_float if total_delta_float is not None else "",
                    total_delta_float + correction if total_delta_float is not None and correction is not None else "",
                ]
            )


def write_expectations(
    path: Path,
    values: np.ndarray,
    expectations: list[dict[str, float]],
) -> None:
    keys = list(expectations[0].keys()) if expectations else []
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["state", "energy_cm-1", "energy_above_ground_cm-1", *keys])
        ground = float(values[0])
        for i, row in enumerate(expectations):
            writer.writerow([i, values[i], values[i] - ground, *[row[key] for key in keys]])


def normalized_temperatures(temperatures: list[float] | None) -> list[float]:
    values = [0.0]
    for temperature in temperatures or []:
        if not any(abs(float(temperature) - existing) <= 1.0e-12 for existing in values):
            values.append(float(temperature))
    return values


def thermal_expectation_rows(
    values: np.ndarray,
    expectations: list[dict[str, float]],
    temperatures: list[float] | None,
) -> tuple[list[str], list[list[float | int]]]:
    keys = list(expectations[0].keys()) if expectations else []
    if not keys:
        return keys, []
    relative = np.asarray(values, dtype=float) - float(values[0])
    rows: list[list[float | int]] = []
    for temperature in normalized_temperatures(temperatures):
        if temperature <= 0.0:
            weights = np.zeros_like(relative)
            weights[0] = 1.0
            partition = 1.0
        else:
            weights = np.exp(-relative / (KB_CM_PER_K * temperature))
            partition = float(np.sum(weights))
            weights /= partition
        row: list[float | int] = [
            float(temperature),
            int(len(relative)),
            partition,
            float(weights @ relative),
        ]
        for key in keys:
            state_values = np.array([item[key] for item in expectations], dtype=float)
            row.append(float(weights @ state_values))
        rows.append(row)
    return keys, rows


def write_thermal_expectations(
    path: Path,
    values: np.ndarray,
    expectations: list[dict[str, float]],
    temperatures: list[float] | None,
) -> None:
    keys, rows = thermal_expectation_rows(values, expectations, temperatures)
    if not keys:
        return
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "temperature_K",
                "states_used",
                "partition_function",
                "thermal_energy_above_ground_cm-1",
                *keys,
            ]
        )
        writer.writerows(rows)


def write_summary(
    path: Path,
    args: argparse.Namespace,
    n_grid: int,
    input_path: Path,
    structures: list[Structure],
    s_sqrtamu_angstrom: np.ndarray,
    s_au: np.ndarray,
    angular_residuals: np.ndarray,
    energy_key: str,
    energy_unit: str,
    property_keys: list[str],
    values: np.ndarray,
    model_info: dict[str, float | int | str] | None = None,
) -> None:
    lines = [
        "Mass-weighted path Hamiltonian",
        f"Input: {input_path}",
        f"Points: {len(structures)}",
        f"Boundary: {args.boundary}",
        f"Solver: {getattr(args, 'solver_used', args.solver)}",
        f"Grid points: {n_grid}",
        f"Energy key: {energy_key}",
        f"Energy unit: {energy_unit}",
        f"Periodic endpoint handling: {args.periodic_endpoints}",
        f"Path length: {s_sqrtamu_angstrom[-1]:.10f} sqrt(amu) Angstrom",
        f"Path length: {s_au[-1]:.10f} a.u.",
        "Orientation: local mass-weighted Eckart alignment of consecutive geometries",
        f"Maximum Eckart angular residual: {np.max(angular_residuals):.8e} amu Angstrom^2",
        "Reduced mass of path coordinate: 1",
        (
            f"Cremer-Pople labels: generalized ring coordinates from ring {args.ring}"
            if args.label_cremer_pople
            else "Cremer-Pople labels: not written"
        ),
        (
            "Gaussian GIC labels: written when present in the Gaussian log; "
            "QPck/PhiP are mapped to Cremer-Pople components when possible."
        ),
        f"Properties: {', '.join(property_keys) if property_keys else 'none'}",
        (
            "Thermal expectation temperatures: "
            + ", ".join(f"{value:g} K" for value in normalized_temperatures(args.temperature))
            if property_keys
            else "Thermal expectation temperatures: none"
        ),
    ]
    if model_info:
        lines.append("One-dimensional grid model:")
        for key in sorted(model_info):
            value = model_info[key]
            if isinstance(value, float):
                lines.append(f"  {key}: {value:.8g}")
            else:
                lines.append(f"  {key}: {value}")
    solver_info = getattr(args, "solver_info", None)
    if solver_info:
        lines.append("Solver diagnostics:")
        for key in sorted(solver_info):
            value = solver_info[key]
            if isinstance(value, float):
                lines.append(f"  {key}: {value:.8g}")
            else:
                lines.append(f"  {key}: {value}")
    lines.extend(["", "Lowest levels:"])
    ground = float(values[0])
    for i, value in enumerate(values[: min(10, len(values))]):
        lines.append(f"  {i:3d} {value:14.6f} cm^-1  {value - ground:14.6f} from ground")
    path.write_text("\n".join(lines) + "\n")


def write_anharmonic_derivative_summary(
    path: Path,
    args: argparse.Namespace,
    input_path: Path,
    derivatives: AnharmonicModeDerivatives,
    potential: AnharmonicDerivativePotential,
    kinetic_frequency_cm: float,
    values: np.ndarray,
    vpt2_comparison: AnharmonicVPT2Comparison | None = None,
) -> None:
    lines = [
        "Gaussian anharmonic one-mode variational calculation",
        f"Input: {input_path}",
        f"Requested mode: {derivatives.requested_mode}",
        f"Mode order: {derivatives.mode_order}",
        f"Force-table mode: {derivatives.force_mode}",
        f"Number of active modes: {derivatives.n_modes}",
        f"Derivative source: {derivatives.source}",
        f"Reference point: {potential.info.get('reference_point', 'unknown')}",
        f"Resolved well type: {potential.info.get('resolved_well_type', 'unknown')}",
        f"Derivative model: {potential.info.get('derivative_model', args.anharmonic_derivative_model)}",
        f"F2: {derivatives.f2_cm:.10f} cm^-1",
        f"F3: {derivatives.f3_cm:.10f} cm^-1",
        f"F4: {derivatives.f4_cm:.10f} cm^-1",
        f"Kinetic prefactor: {kinetic_frequency_cm:.10f} cm^-1",
        f"Grid points: {len(potential.q)}",
        f"Grid range: {potential.q[0]:.10f} to {potential.q[-1]:.10f} Q",
        "Coordinate: Gaussian dimensionless normal coordinate; reduced mass is one in this coordinate.",
    ]
    if derivatives.frequency_cm is not None:
        lines.append(f"Standard frequency for requested mode: {derivatives.frequency_cm:.10f} cm^-1")
    if derivatives.notes:
        lines.append(f"Parser notes: {', '.join(derivatives.notes)}")
    if potential.info:
        lines.append("Model diagnostics:")
        for key in sorted(potential.info):
            value = potential.info[key]
            if isinstance(value, float):
                lines.append(f"  {key}: {value:.8g}")
            else:
                lines.append(f"  {key}: {value}")
    if vpt2_comparison is None:
        if potential.info.get("resolved_well_type") == "double":
            lines.append("One-mode VPT2 comparison: skipped for a double-well/TS reference")
        else:
            lines.append("One-mode VPT2 comparison: skipped because F2 is not positive")
    else:
        lines.append("One-mode VPT2 comparison:")
        for key in sorted(vpt2_comparison.info):
            value = vpt2_comparison.info[key]
            if isinstance(value, float):
                lines.append(f"  {key}: {value:.8g}")
            else:
                lines.append(f"  {key}: {value}")
    lines.extend(["", "Lowest levels:"])
    ground = float(values[0])
    for i, value in enumerate(values[: min(10, len(values))]):
        lines.append(f"  {i:3d} {value:14.6f} cm^-1  {value - ground:14.6f} from ground")
    path.write_text("\n".join(lines) + "\n")


def plot_results(
    path_pdf: Path,
    path_png: Path,
    grid_au: np.ndarray,
    potential_grid_cm: np.ndarray,
    property_grid: dict[str, np.ndarray],
    property_expectations: list[dict[str, float]],
    values: np.ndarray,
    vectors: np.ndarray,
    *,
    states: int,
    plot_max_state: int,
    property_smooth_degree: int = 0,
    vpt2_properties: dict[str, dict[str, float]] | None = None,
) -> None:
    vpt2_properties = vpt2_properties or {}
    x = grid_au / MW_ANGSTROM_TO_AU
    property_priority = [
        "A_MHz",
        "B_MHz",
        "C_MHz",
        "A_GHz",
        "B_GHz",
        "C_GHz",
        "dipole_debye",
        "dipole_x_debye",
        "dipole_y_debye",
        "dipole_z_debye",
    ]
    property_names = [key for key in property_priority if key in property_grid]
    if not property_names:
        property_names = list(property_grid)[:3]
    property_names = property_names[:3]

    max_state = min(max(int(plot_max_state), 0), len(values) - 1)
    states_to_draw = min(max(states, max_state + 1), len(values))
    selected_energy = float(values[max_state])
    y_top = max(selected_energy * 1.15, selected_energy + 20.0, 25.0)
    visible = potential_grid_cm <= y_top
    if np.any(visible):
        indices = np.where(visible)[0]
        left = max(int(indices[0]) - 2, 0)
        right = min(int(indices[-1]) + 3, len(x))
    else:
        left, right = 0, len(x)

    if property_names:
        fig = plt.figure(figsize=(10.8, 5.9), constrained_layout=True)
        gs = fig.add_gridspec(len(property_names), 2, width_ratios=[1.45, 1.0])
        ax = fig.add_subplot(gs[:, 0])
        prop_axes = [fig.add_subplot(gs[i, 1]) for i in range(len(property_names))]
    else:
        fig, ax = plt.subplots(figsize=(6.8, 4.2), constrained_layout=True)
        prop_axes = []

    ax.plot(x, potential_grid_cm, color="#1f5a7a", lw=2.0)
    state_colors = {
        0: "#2563eb",
        1: "#d97706",
    }
    for state, value in enumerate(values[:states_to_draw]):
        if value > y_top:
            continue
        state_color = state_colors.get(state, "#8f2d2d")
        wave_color = state_colors.get(state, "#5f8a5a")
        ax.axhline(value, color=state_color, lw=1.0 if state < 2 else 0.8, alpha=0.85 if state < 2 else 0.65)
        scale = 0.08 * max(y_top, 1.0)
        psi = vectors[:, state]
        psi = psi / max(np.max(np.abs(psi)), 1.0e-12)
        ax.plot(
            x,
            value + scale * psi,
            color=wave_color,
            lw=1.2 if state < 2 else 0.8,
            alpha=0.95 if state < 2 else 0.8,
        )
        ax.text(
            x[right - 1],
            value,
            f"v={state}",
            ha="right",
            va="bottom",
            fontsize=8,
            color=state_color,
        )
    ax.set_xlim(float(x[left]), float(x[right - 1]))
    ax.set_ylim(min(-0.05 * y_top, -5.0), y_top)
    ax.set_xlabel(r"$s$ / $\sqrt{\mathrm{amu}}$ Angstrom")
    ax.set_ylabel(r"Energy / cm$^{-1}$")
    ax.set_title(f"Potential and levels up to state {max_state}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color="#d8d8d8", lw=0.5, alpha=0.7)

    colors = ["#0f766e", "#a16207", "#7c2d12"]
    prop_x_min = float(x[left])
    prop_x_max = float(x[right - 1])
    prop_x_span = max(prop_x_max - prop_x_min, 1.0e-12)
    prop_label_x = prop_x_max + 0.04 * prop_x_span
    prop_axis_xmax = prop_x_max + 0.48 * prop_x_span
    for i, (key, prop_ax) in enumerate(zip(property_names, prop_axes)):
        values_grid = property_grid[key]
        color = colors[i % len(colors)]
        if property_smooth_degree > 0 and len(x) > property_smooth_degree + 1:
            scaled_x = 2.0 * (x - float(np.min(x))) / max(float(np.ptp(x)), 1.0e-12) - 1.0
            degree = min(int(property_smooth_degree), len(x) - 1)
            polynomial = np.polynomial.Chebyshev.fit(scaled_x, values_grid, degree)
            plot_values = polynomial(scaled_x)
        else:
            plot_values = values_grid
        prop_ax.plot(x, plot_values, color=color, lw=1.8)
        if property_expectations:
            variational = property_expectations[0].get(key)
            if variational is not None:
                prop_ax.axhline(variational, color=color, ls="--", lw=1.1)
                prop_ax.annotate(
                    f"DVR {variational:.3f}",
                    xy=(prop_label_x, variational),
                    xycoords="data",
                    xytext=(0, 13),
                    textcoords="offset points",
                    ha="left",
                    va="bottom",
                    fontsize=8,
                    color=color,
                    bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.82},
                )
        if key in vpt2_properties and "vpt2_1d_value" in vpt2_properties[key]:
            vpt2_value = float(vpt2_properties[key]["vpt2_1d_value"])
            prop_ax.axhline(vpt2_value, color="#334155", ls=":", lw=1.1)
            if property_expectations and key in property_expectations[0]:
                correction = property_expectations[0][key] - vpt2_value
                prop_ax.annotate(
                    f"DVR-VPT2 {correction:+.3f}",
                    xy=(prop_label_x, vpt2_value),
                    xycoords="data",
                    xytext=(0, -13),
                    textcoords="offset points",
                    ha="left",
                    va="top",
                    fontsize=8,
                    color="#334155",
                    bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.82},
                )
        prop_ax.set_xlim(prop_x_min, prop_axis_xmax)
        prop_ax.set_ylabel(key)
        prop_ax.grid(color="#d8d8d8", lw=0.5, alpha=0.7)
        prop_ax.spines["top"].set_visible(False)
        prop_ax.spines["right"].set_visible(False)
        if i < len(prop_axes) - 1:
            prop_ax.tick_params(labelbottom=False)
        else:
            prop_ax.set_xlabel(r"$s$ / $\sqrt{\mathrm{amu}}$ Angstrom")

    fig.savefig(path_pdf)
    fig.savefig(path_png, dpi=240)
    plt.close(fig)


def plot_anharmonic_derivative_results(
    path_pdf: Path,
    path_png: Path,
    q: np.ndarray,
    potential_cm: np.ndarray,
    values: np.ndarray,
    vectors: np.ndarray,
    *,
    states: int,
) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 4.2), constrained_layout=True)
    ax.plot(q, potential_cm, color="#1f5a7a", lw=2.0)
    for state, value in enumerate(values[:states]):
        ax.axhline(value, color="#8f2d2d", lw=0.8, alpha=0.65)
        scale = 0.08 * max(float(np.ptp(potential_cm)), 1.0)
        psi = vectors[:, state]
        psi = psi / max(np.max(np.abs(psi)), 1.0e-12)
        ax.plot(q, value + scale * psi, color="#3f6f3b", lw=0.8, alpha=0.8)
    ax.set_xlabel(r"$Q$ / dimensionless normal coordinate")
    ax.set_ylabel(r"Energy / cm$^{-1}$")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color="#d8d8d8", lw=0.5, alpha=0.7)
    fig.savefig(path_pdf)
    fig.savefig(path_png, dpi=240)
    plt.close(fig)


def plot_potential_profile(
    path_pdf: Path,
    path_png: Path,
    rel_energy_cm: np.ndarray,
    *,
    period_deg: float,
) -> None:
    phase = np.linspace(0.0, period_deg, len(rel_energy_cm))
    energy = np.array(rel_energy_cm, dtype=float).copy()
    if len(energy) > 1:
        energy[-1] = energy[0]

    fig, ax = plt.subplots(figsize=(7.2, 4.3), constrained_layout=True)
    ax.plot(phase, energy, color="#1f5a7a", lw=2.4)
    ax.fill_between(phase, 0.0, energy, color="#1f5a7a", alpha=0.10)

    extrema: list[tuple[str, int]] = [("minimum", int(np.argmin(energy)))]
    local_maxima = [
        i
        for i in range(1, len(energy) - 1)
        if energy[i] >= energy[i - 1] and energy[i] >= energy[i + 1]
    ]
    for i in sorted(local_maxima, key=lambda item: energy[item], reverse=True)[:2]:
        extrema.append(("barrier", i))
    mid = int(np.argmin(np.abs(phase - 0.5 * period_deg)))
    extrema.append(("midpoint", mid))

    seen: set[int] = set()
    for kind, i in extrema:
        if i in seen:
            continue
        seen.add(i)
        marker = "o" if kind == "minimum" else "^" if kind == "barrier" else "s"
        color = "#1f5a7a" if kind == "minimum" else "#9b2f2f" if kind == "barrier" else "#555555"
        ax.scatter([phase[i]], [energy[i]], s=42, marker=marker, color=color, zorder=4)
        offset = -38 if kind == "barrier" else 18
        ax.annotate(
            f"{phase[i]:.0f} deg\n{energy[i]:.1f} cm$^{{-1}}$",
            xy=(phase[i], energy[i]),
            xytext=(0, offset),
            textcoords="offset points",
            ha="center",
            va="top" if kind == "barrier" else "bottom",
            fontsize=8,
            color=color,
            bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.82},
        )

    ax.annotate(
        f"{period_deg:.0f} deg equivalent to 0 deg",
        xy=(period_deg, energy[-1]),
        xytext=(-8, 28),
        textcoords="offset points",
        ha="right",
        va="bottom",
        fontsize=8,
        color="#1f5a7a",
        arrowprops={"arrowstyle": "->", "lw": 0.7, "color": "#1f5a7a"},
    )
    ax.set_xlim(0.0, period_deg)
    ax.set_ylim(min(-2.0, float(np.min(energy)) - 2.0), float(np.max(energy)) + 8.0)
    ax.set_xlabel(r"Relative pseudorotation phase / degrees")
    ax.set_ylabel(r"Relative energy / cm$^{-1}$")
    ax.set_title("Adiabatic pseudorotation potential", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color="#d8d8d8", lw=0.5, alpha=0.7)
    fig.savefig(path_pdf)
    fig.savefig(path_png, dpi=240)
    plt.close(fig)


def plot_circular_potential(
    path_pdf: Path,
    path_png: Path,
    grid_au: np.ndarray,
    potential_grid_cm: np.ndarray,
    values: np.ndarray,
    *,
    states: int,
) -> None:
    if len(grid_au) < 2:
        return
    span = float(grid_au[-1] - grid_au[0])
    if span <= 0.0:
        return

    theta = 2.0 * math.pi * (grid_au - grid_au[0]) / span
    theta = np.asarray(theta, dtype=float)
    potential = np.asarray(potential_grid_cm, dtype=float)
    vmin = float(np.min(potential))
    vmax = float(np.max(potential))
    vrange = max(vmax - vmin, 1.0)
    radius = 1.0 + 0.18 * (potential - vmin) / vrange

    phase_deg = 180.0 * (grid_au - grid_au[0]) / span

    fig = plt.figure(figsize=(5.8, 5.2), constrained_layout=True)
    ax = fig.add_subplot(111, projection="polar")
    ax.plot(theta, np.ones_like(theta), color="#b8b8b8", lw=0.8, zorder=1)
    points = ax.scatter(
        theta,
        radius,
        c=potential,
        cmap="viridis",
        s=18,
        linewidths=0.0,
        zorder=3,
    )

    extrema_indices = [
        int(np.argmin(potential)),
        int(np.argmax(potential)),
    ]
    if len(potential) >= 5:
        local_maxima = [
            i
            for i in range(1, len(potential) - 1)
            if potential[i] >= potential[i - 1] and potential[i] >= potential[i + 1]
        ]
        ranked = sorted(local_maxima, key=lambda i: potential[i], reverse=True)[:2]
        extrema_indices.extend(ranked)
    seen: set[int] = set()
    for i in extrema_indices:
        if i in seen:
            continue
        seen.add(i)
        if potential[i] - vmin < 1.0:
            label = f"0/180 deg\nminimum"
            text_theta = theta[i] - 0.12
            text_radius = radius[i] + 0.22
        else:
            label = f"{phase_deg[i]:.0f} deg\n{potential[i] - vmin:.1f} cm$^{{-1}}$"
            text_theta = theta[i]
            text_radius = radius[i] + 0.16
        ax.annotate(
            label,
            xy=(theta[i], radius[i]),
            xytext=(text_theta, text_radius),
            ha="center",
            va="center",
            fontsize=8,
            color="#222222",
            bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.82},
            arrowprops={"arrowstyle": "-", "lw": 0.5, "color": "#555555"},
        )

    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.set_xticks(np.deg2rad([90, 180, 270]))
    ax.set_xticklabels(["45", "90", "135"])
    ax.set_yticks([])
    ax.set_ylim(0.0, 1.42)
    ax.spines["polar"].set_visible(False)
    ax.grid(color="#d0d0d0", lw=0.5, alpha=0.6)
    cbar = fig.colorbar(points, ax=ax, shrink=0.72, pad=0.08)
    cbar.set_label(r"Energy / cm$^{-1}$")
    ax.set_title("Potential on the pseudorotation circle", pad=16)
    fig.savefig(path_pdf)
    fig.savefig(path_png, dpi=240)
    plt.close(fig)


def read_structures(args: argparse.Namespace) -> SourceData:
    if args.xyz:
        path = args.xyz.expanduser().resolve()
        return SourceData(read_xyz(path), path, charge=args.charge or 0, multiplicity=args.multiplicity or 1)
    if args.gaussian_log:
        path = args.gaussian_log.expanduser().resolve()
        charge, multiplicity = read_gaussian_charge_multiplicity(path)
        if args.charge is not None:
            charge = args.charge
        if args.multiplicity is not None:
            multiplicity = args.multiplicity
        return SourceData(
            read_gaussian_log(path, selection=args.log_selection, energy_source=args.gaussian_energy),
            path,
            charge=charge,
            multiplicity=multiplicity,
        )
    raise ValueError("Pass either --xyz or --gaussian-log")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="General one-dimensional Hamiltonian on a mass-weighted Cartesian path."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--xyz", type=Path, help="Multi-XYZ path file.")
    source.add_argument(
        "--gaussian-log",
        type=Path,
        help="Gaussian log with optimized scan/path points; the scanned coordinate does not have to be puckering.",
    )
    source.add_argument("--grid2d-csv", type=Path, help="Rectangular 2D grid CSV for product-basis analysis.")
    parser.add_argument(
        "--log-selection",
        choices=["all", "last", "last-per-link"],
        default="all",
        help="Which optimized structures to read from a Gaussian log.",
    )
    parser.add_argument(
        "--gaussian-energy",
        choices=["auto", "scf", "post-scf"],
        default="auto",
        help=(
            "Energy to attach to Gaussian geometries. auto uses post-SCF total energies "
            "such as E(method) or EUMP2 when present, otherwise SCF Done; use post-scf "
            "for double-hybrid/MP2 scans."
        ),
    )
    parser.add_argument(
        "--anharmonic-mode",
        type=int,
        help=(
            "Run a one-mode variational calculation from Gaussian anharmonic "
            "force constants for this mode instead of reading path geometries."
        ),
    )
    parser.add_argument(
        "--anharmonic-mode-order",
        choices=["frequency-ascending", "force-table"],
        default="frequency-ascending",
        help=(
            "Mode numbering used by --anharmonic-mode. frequency-ascending uses "
            "the standard Gaussian Frequencies block; force-table uses the "
            "quadratic/cubic/quartic force-constant tables directly."
        ),
    )
    parser.add_argument(
        "--anharmonic-derivative-model",
        choices=[
            "auto",
            "handy-morse",
            "handy-gaussian",
            "gaussian",
            "morse",
            "inverse-power",
            "taylor",
        ],
        default="auto",
        help=(
            "Analytic one-mode model fitted to Gaussian F2/F3/F4 derivatives. "
            "auto uses the Handy Gauss-like coordinate when F3 is zero, the "
            "Handy Morse-like coordinate for a single minimum with non-zero F3, "
            "and inverse-power for TS/double-well cases."
        ),
    )
    parser.add_argument(
        "--anharmonic-handy-beta",
        type=float,
        help=(
            "Optional beta for the Handy Gauss-like coordinate z=sqrt(1-exp(-beta Q^2)). "
            "If omitted, beta is estimated from F2/F4 or the grid width."
        ),
    )
    parser.add_argument(
        "--anharmonic-grid-half-width",
        type=float,
        default=6.0,
        help="Half-width of the dimensionless normal-coordinate DVR grid.",
    )
    parser.add_argument(
        "--anharmonic-cubic-threshold",
        type=float,
        default=1.0e-8,
        help="Absolute F3 threshold in cm^-1 used to choose the symmetric Gaussian model.",
    )
    parser.add_argument(
        "--anharmonic-wall-height-cm",
        type=float,
        default=10000.0,
        help=(
            "Confining high-order wall height at the grid edge in cm^-1. "
            "The wall has zero derivatives through fourth order at Q=0."
        ),
    )
    parser.add_argument(
        "--anharmonic-wall-degree",
        type=int,
        choices=[6, 8],
        default=8,
        help="Degree of the high-order confining wall used by --anharmonic-mode.",
    )
    parser.add_argument(
        "--anharmonic-kinetic-frequency-cm",
        type=float,
        help=(
            "Override the normal-coordinate kinetic prefactor in cm^-1. "
            "By default abs(F2) is used, which handles imaginary TS modes."
        ),
    )
    parser.add_argument("--properties-csv", type=Path, help="CSV with energies/properties along the path.")
    parser.add_argument(
        "--property-derivatives-csv",
        type=Path,
        help=(
            "CSV with one-dimensional property derivatives in the mass-scaled path coordinate. "
            "Columns: property,value,d1,d2,... with optional origin and parity."
        ),
    )
    parser.add_argument(
        "--vpt2-property-csv",
        type=Path,
        help=(
            "CSV with one-dimensional perturbative property values for comparison with "
            "variational averages. Columns: property plus vpt2_value/value or "
            "vpt2_delta/delta; optional total_perturbative_value/total_perturbative_delta."
        ),
    )
    parser.add_argument(
        "--property-vpt2",
        choices=["auto", "csv-only", "off"],
        default="auto",
        help=(
            "How to obtain one-dimensional perturbative property values. auto fits a local "
            "Taylor model near the 1D minimum and lets --vpt2-property-csv override selected "
            "properties; csv-only writes only values supplied by the CSV; off disables the comparison."
        ),
    )
    parser.add_argument(
        "--property-vpt2-fit-points",
        type=int,
        default=7,
        help="Number of nearest 1D model points used for automatic local property VPT2 fits.",
    )
    parser.add_argument(
        "--property-vpt2-degree",
        type=int,
        choices=[2, 3, 4],
        default=4,
        help="Taylor degree used for automatic local one-dimensional property VPT2 fits.",
    )
    parser.add_argument(
        "--property-vpt2-basis-size",
        type=int,
        default=24,
        help="Harmonic basis size used for automatic perturbative property expectation values.",
    )
    parser.add_argument("--energy-key", help="Column/property name containing the potential energy.")
    parser.add_argument(
        "--energy-unit",
        choices=["auto", "hartree", "cm-1", "kjmol", "kcalmol"],
        default="auto",
    )
    parser.add_argument(
        "--property",
        action="append",
        help="Scalar property to average over DVR states. Can be used multiple times.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        action="append",
        help=(
            "Temperature in K for Boltzmann-averaged one-dimensional properties. "
            "Can be used multiple times; 0 K ground-state averages are always written when properties exist."
        ),
    )
    parser.add_argument(
        "--plot-max-state",
        type=int,
        default=5,
        help="Highest vibrational state index used to set the energy range of 1D plots.",
    )
    parser.add_argument(
        "--plot-property-smooth-degree",
        type=int,
        default=0,
        help="Optional Chebyshev degree used only to smooth property curves in 1D plots; 0 disables.",
    )
    parser.add_argument(
        "--compute-rotconst",
        action="store_true",
        help="Compute A/B/C rotational constants from each Cartesian geometry.",
    )
    parser.add_argument(
        "--boundary",
        choices=["periodic", "nonperiodic"],
        default="nonperiodic",
        help="Boundary condition for the one-dimensional path Hamiltonian.",
    )
    parser.add_argument(
        "--solver",
        choices=["auto", "fourier", "gaussian", "sinc-dvr"],
        default="auto",
        help=(
            "Hamiltonian solver. auto uses Fourier for periodic paths and an optimized "
            "distributed Gaussian basis for nonperiodic paths; sinc-dvr is retained as a "
            "robust fallback."
        ),
    )
    parser.add_argument("--q1-key", help="First coordinate column for --grid2d-csv.")
    parser.add_argument("--q2-key", help="Second coordinate column for --grid2d-csv.")
    parser.add_argument(
        "--q1-scale-au",
        type=float,
        default=1.0,
        help="Scale applied to the first 2D coordinate before solving; default assumes a.u.",
    )
    parser.add_argument(
        "--q2-scale-au",
        type=float,
        default=1.0,
        help="Scale applied to the second 2D coordinate before solving; default assumes a.u.",
    )
    parser.add_argument(
        "--boundary1",
        choices=["periodic", "nonperiodic"],
        default="nonperiodic",
        help="Boundary condition for the first coordinate in --grid2d-csv mode.",
    )
    parser.add_argument(
        "--boundary2",
        choices=["periodic", "nonperiodic"],
        default="nonperiodic",
        help="Boundary condition for the second coordinate in --grid2d-csv mode.",
    )
    parser.add_argument(
        "--solver1",
        choices=["auto", "fourier", "gaussian", "sinc-dvr"],
        default="auto",
        help="One-dimensional pruning solver for the first 2D coordinate.",
    )
    parser.add_argument(
        "--solver2",
        choices=["auto", "fourier", "gaussian", "sinc-dvr"],
        default="auto",
        help="One-dimensional pruning solver for the second 2D coordinate.",
    )
    parser.add_argument("--basis1", type=int, default=12, help="Retained 1D functions for q1 in 2D mode.")
    parser.add_argument("--basis2", type=int, default=12, help="Retained 1D functions for q2 in 2D mode.")
    parser.add_argument(
        "--prune-energy-window1",
        type=float,
        default=0.0,
        help="Retain q1 pruning functions up to this energy above the 1D ground state; 0 keeps --basis1.",
    )
    parser.add_argument(
        "--prune-energy-window2",
        type=float,
        default=0.0,
        help="Retain q2 pruning functions up to this energy above the 1D ground state; 0 keeps --basis2.",
    )
    parser.add_argument("--prune-min-basis1", type=int, default=1, help="Minimum q1 functions retained by energy pruning.")
    parser.add_argument("--prune-min-basis2", type=int, default=1, help="Minimum q2 functions retained by energy pruning.")
    parser.add_argument(
        "--prune-product-coeff-threshold",
        type=float,
        default=0.0,
        help="Optional coefficient-weight threshold for pruning the 2D product basis after a first diagonalization.",
    )
    parser.add_argument(
        "--prune-product-states",
        type=int,
        default=3,
        help="Number of low-lying 2D states used for product coefficient pruning.",
    )
    parser.add_argument(
        "--prune-potential",
        choices=["min", "mean", "cut"],
        default="min",
        help="Reference 1D potentials used for 2D basis pruning.",
    )
    parser.add_argument("--q1-ref", type=float, help="q1 value for --prune-potential cut.")
    parser.add_argument("--q2-ref", type=float, help="q2 value for --prune-potential cut.")
    parser.add_argument("--g11", type=float, default=1.0, help="Constant kinetic metric element g11 for 2D mode.")
    parser.add_argument("--g22", type=float, default=1.0, help="Constant kinetic metric element g22 for 2D mode.")
    parser.add_argument("--g12", type=float, default=0.0, help="Constant kinetic metric element g12 for 2D mode.")
    parser.add_argument(
        "--metric-mode",
        choices=["constant", "csv", "geometry"],
        default="constant",
        help=(
            "Kinetic metric for --grid2d-csv. constant uses --g11/--g22/--g12; "
            "csv reads coordinate-dependent g11/g22/g12 columns; geometry computes "
            "the metric from a matching grid of Cartesian structures."
        ),
    )
    parser.add_argument(
        "--metric-stencil",
        type=int,
        default=5,
        help=(
            "Number of grid points used for each geometry-derived finite-difference "
            "derivative. The code reduces this automatically on small grids."
        ),
    )
    parser.add_argument(
        "--metric-smoothing",
        choices=["auto", "none", "fit"],
        default="auto",
        help=(
            "Smoothing applied to coordinate-dependent metrics. auto currently uses "
            "a low-order tensor polynomial/Fourier fit for csv and geometry metrics."
        ),
    )
    parser.add_argument(
        "--metric-smooth-order",
        type=int,
        default=3,
        help="Maximum total order/frequency used for metric smoothing fits.",
    )
    parser.add_argument(
        "--metric-eigenvalue-floor",
        type=float,
        default=1.0e-10,
        help="Relative eigenvalue floor used when projecting metric matrices to positive definite form.",
    )
    parser.add_argument("--g11-key", help="CSV column for g11 when --metric-mode csv is used.")
    parser.add_argument("--g22-key", help="CSV column for g22 when --metric-mode csv is used.")
    parser.add_argument("--g12-key", help="CSV column for g12 when --metric-mode csv is used.")
    parser.add_argument(
        "--grid2d-geom-xyz",
        type=Path,
        help="Multi-XYZ geometry grid, in the same row order as --grid2d-csv, for --metric-mode geometry.",
    )
    parser.add_argument(
        "--grid2d-geom-log",
        type=Path,
        help="Gaussian log geometry grid, in the same row order as --grid2d-csv, for --metric-mode geometry.",
    )
    parser.add_argument(
        "--grid2d-geom-log-selection",
        choices=["all", "last", "last-per-link"],
        default="last-per-link",
        help="Gaussian optimized structures read from --grid2d-geom-log for metric construction.",
    )
    parser.add_argument(
        "--convergence-scan",
        action="store_true",
        help="Run a 2D convergence sweep over grid strides, basis sizes, and pruning potentials.",
    )
    parser.add_argument(
        "--convergence-basis1",
        default="",
        help="Comma-separated q1 basis sizes for --convergence-scan. Default: half/current --basis1.",
    )
    parser.add_argument(
        "--convergence-basis2",
        default="",
        help="Comma-separated q2 basis sizes for --convergence-scan. Default: half/current --basis2.",
    )
    parser.add_argument(
        "--convergence-prune-potentials",
        default="",
        help="Comma-separated pruning potentials for convergence: min,mean,cut. Default: current --prune-potential.",
    )
    parser.add_argument(
        "--convergence-strides",
        default="1",
        help="Comma-separated grid downsampling strides for both axes in --convergence-scan.",
    )
    parser.add_argument(
        "--convergence-levels",
        type=int,
        default=0,
        help="Number of levels written to the convergence table. 0 uses --levels.",
    )
    parser.add_argument(
        "--gaussian-basis-size",
        type=int,
        default=40,
        help="Number of distributed Gaussians for --solver gaussian.",
    )
    parser.add_argument(
        "--gaussian-width-scale",
        type=float,
        default=1.0,
        help="Gaussian exponent scale alpha=scale/dx^2 for --solver gaussian.",
    )
    parser.add_argument(
        "--gaussian-quadrature",
        type=int,
        default=24,
        help="Gauss-Hermite quadrature order for distributed-Gaussian matrix elements.",
    )
    parser.add_argument(
        "--gaussian-optimize",
        choices=["none", "widths", "centers-widths"],
        default="centers-widths",
        help=(
            "Nonlinear optimization of the distributed-Gaussian basis. The default "
            "optimizes internal centers and widths variationally for the lowest states."
        ),
    )
    parser.add_argument(
        "--gaussian-optimization-levels",
        type=int,
        default=0,
        help=(
            "Number of low-lying levels included in the Gaussian basis optimization "
            "objective. 0 uses min(--levels, basis size)."
        ),
    )
    parser.add_argument(
        "--gaussian-optimization-maxiter",
        type=int,
        default=80,
        help="Maximum L-BFGS-B iterations for nonlinear Gaussian basis optimization.",
    )
    parser.add_argument(
        "--gaussian-overlap-threshold",
        type=float,
        default=1.0e-8,
        help="Relative overlap eigenvalue threshold used to remove linear dependencies.",
    )
    parser.add_argument(
        "--gaussian-condition-limit",
        type=float,
        default=1.0e8,
        help="Maximum allowed condition number of the retained Gaussian overlap subspace.",
    )
    parser.add_argument(
        "--periodic-endpoints",
        choices=["average", "first", "last", "none"],
        default="first",
        help="How to enforce equality of first/last scalar values for periodic DVR.",
    )
    parser.add_argument(
        "--path-symmetry",
        choices=["none", "half-even", "half-even-origin", "half-even-last", "center-even"],
        default="none",
        help=(
            "Optional one-dimensional symmetry preprocessing. half-even and "
            "half-even-origin mirror about the first point; half-even-last mirrors "
            "about the last point; center-even averages a complete path about its midpoint."
        ),
    )
    parser.add_argument(
        "--well-type",
        choices=["auto", "single", "double"],
        default="auto",
        help=(
            "Central potential topology for analytic diagnostics/extensions. "
            "auto classifies the symmetrized core; single enables Morse or "
            "inverse-power one-well tails; double enables the double-well criterion."
        ),
    )
    parser.add_argument(
        "--double-well-criterion",
        choices=["auto", "off"],
        default="auto",
        help=(
            "For symmetric nonperiodic 1D paths, estimate the Flanigan-de la Vega "
            "double-well criterion before analytic tail extension."
        ),
    )
    parser.add_argument(
        "--double-well-fit-points",
        type=int,
        default=7,
        help="Radial core points used for the local minimum curvature in the double-well criterion.",
    )
    parser.add_argument(
        "--core-model",
        choices=["auto", "sampled", "asymmetric-parabola-gaussian"],
        default="auto",
        help=(
            "Analytic model for the central one-dimensional potential before tail extension. "
            "auto keeps the sampled profile except for sparse unsymmetrized double wells, "
            "where it uses the Flanigan-de la Vega parabola+Gaussian form."
        ),
    )
    parser.add_argument(
        "--core-auto-max-points",
        type=int,
        default=7,
        help="Maximum number of core sample points for --core-model auto to replace a double well by an analytic core.",
    )
    parser.add_argument(
        "--potential-extension",
        choices=[
            "none",
            "repulsive-polynomial",
            "repulsive-exponential",
            "morse-polynomial",
            "single-morse",
            "single-inverse-power",
            "normal-mode-polynomial",
            "repulsive-quartic",
            "morse-quartic",
        ],
        default="none",
        help=(
            "Add analytic nonperiodic repulsive tails before interpolation. "
            "repulsive-polynomial matches endpoint slope/curvature and adds a "
            "convex degree-6/8 wall; repulsive-exponential uses a Born-Mayer-like "
            "exponential wall; morse-polynomial fits a symmetric double-Morse core; "
            "single-morse and single-inverse-power are for --well-type single; "
            "normal-mode-polynomial uses a Gaussian minimum frequency log to set "
            "the endpoint curvature from the normal mode most parallel to the scan. "
            "The old *-quartic names are accepted as legacy aliases."
        ),
    )
    parser.add_argument(
        "--extension-length-au",
        type=float,
        default=0.0,
        help="Length of each analytic repulsive tail in mass-weighted a.u.; 0 disables tails.",
    )
    parser.add_argument(
        "--extension-points",
        type=int,
        default=24,
        help="Number of analytic points added to each side when --potential-extension is enabled.",
    )
    parser.add_argument(
        "--extension-degree",
        type=int,
        choices=[6, 8],
        default=6,
        help="Degree of the repulsive polynomial wall added beyond the computed points.",
    )
    parser.add_argument(
        "--extension-target-cm",
        type=float,
        default=10000.0,
        help="Target potential height at the end of each analytic tail, in cm^-1.",
    )
    parser.add_argument(
        "--extension-fit-points",
        type=int,
        default=4,
        help="Endpoint points used to estimate slope/curvature for repulsive-polynomial tails.",
    )
    parser.add_argument(
        "--tail-normal-mode-log",
        type=Path,
        help=(
            "Gaussian harmonic/anharmonic log at the endpoint minimum used by "
            "--potential-extension normal-mode-polynomial."
        ),
    )
    parser.add_argument(
        "--tail-normal-mode-side",
        choices=["minimum", "left", "right", "both"],
        default="minimum",
        help="Endpoint side whose curvature is taken from the selected normal mode.",
    )
    parser.add_argument(
        "--tail-normal-mode-max-frequency",
        type=float,
        default=500.0,
        help="Maximum normal-mode frequency, in cm^-1, considered for endpoint-tail selection.",
    )
    parser.add_argument(
        "--tail-normal-mode-min-overlap",
        type=float,
        default=0.0,
        help="Minimum acceptable absolute mass-weighted overlap for the selected tail mode.",
    )
    parser.add_argument(
        "--tail-normal-mode-steps",
        type=int,
        default=2,
        help=(
            "Number of endpoint scan steps used to estimate the tail direction. "
            "The default 2 fits the endpoint tangent over the last three geometries."
        ),
    )
    parser.add_argument(
        "--tail-normal-mode-heavy-atoms",
        action="store_true",
        help="Compute endpoint/mode overlaps using only non-H atoms.",
    )
    parser.add_argument(
        "--tail-normal-mode-slope",
        choices=["zero", "endpoint"],
        default="zero",
        help=(
            "Slope used at a normal-mode-polynomial endpoint. zero is appropriate "
            "when the endpoint is an optimized minimum; endpoint uses the fitted scan slope."
        ),
    )
    parser.add_argument(
        "--potential-smoothing",
        choices=["none", "spline"],
        default="none",
        help="Interpolate the final one-dimensional potential with a smoothing spline.",
    )
    parser.add_argument(
        "--potential-spline-smoothing",
        type=float,
        default=0.0,
        help="UnivariateSpline smoothing factor used when --potential-smoothing spline.",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="Repeat a periodic path this many times before solving the DVR.",
    )
    parser.add_argument("--grid", type=int, default=401, help="Number of interpolation/output grid points.")
    parser.add_argument("--levels", type=int, default=12, help="Number of eigenstates to compute.")
    parser.add_argument("--save-states", type=int, default=6, help="Wavefunctions written to grid CSV.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--figdir", type=Path, default=DEFAULT_FIGS)
    parser.add_argument("--prefix", default="mw_path_dvr")
    parser.add_argument("--check-only", action="store_true", help="Read and orient the path only.")
    parser.add_argument("--charge", type=int, help="Override molecular charge metadata read from Gaussian.")
    parser.add_argument("--multiplicity", type=int, help="Override spin multiplicity metadata read from Gaussian.")
    parser.add_argument(
        "--ring",
        default="1,2,3,4,5",
        help=(
            "One-based atom indices of a ring, in cyclic order, used only for ring-puckering "
            "Cremer-Pople labeling of already computed Gaussian/XYZ geometries."
        ),
    )
    parser.add_argument(
        "--label-cremer-pople",
        action="store_true",
        help=(
            "Write generalized Cremer-Pople labels from Cartesian geometries. If Gaussian "
            "GIC QPck/PhiP values are present, also write the fitted GIC-to-Cremer-Pople map."
        ),
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    extension_mode = canonical_extension_mode(args.potential_extension)
    if args.property_derivatives_csv and (args.grid2d_csv or args.anharmonic_mode is not None):
        raise ValueError("--property-derivatives-csv is supported only for 1D path calculations")
    if args.vpt2_property_csv and (args.grid2d_csv or args.anharmonic_mode is not None):
        raise ValueError("--vpt2-property-csv is supported only for 1D path calculations")
    if args.plot_max_state < 0:
        raise ValueError("--plot-max-state must be non-negative")
    if args.plot_property_smooth_degree < 0:
        raise ValueError("--plot-property-smooth-degree must be non-negative")
    if args.property_vpt2_fit_points < 3:
        raise ValueError("--property-vpt2-fit-points must be at least 3")
    if args.property_vpt2_basis_size < 8:
        raise ValueError("--property-vpt2-basis-size must be at least 8")
    for temperature in args.temperature or []:
        if temperature < 0.0:
            raise ValueError("--temperature must be non-negative")
    if args.anharmonic_mode is not None:
        if not args.gaussian_log:
            raise ValueError("--anharmonic-mode requires --gaussian-log")
        if args.boundary != "nonperiodic":
            raise ValueError("--anharmonic-mode uses a nonperiodic normal-coordinate DVR")
        if args.solver not in {"auto", "sinc-dvr"}:
            raise ValueError("--anharmonic-mode supports --solver auto or --solver sinc-dvr")
        if args.anharmonic_mode < 1:
            raise ValueError("--anharmonic-mode must be positive")
        if args.anharmonic_grid_half_width <= 0.0:
            raise ValueError("--anharmonic-grid-half-width must be positive")
        if args.anharmonic_cubic_threshold < 0.0:
            raise ValueError("--anharmonic-cubic-threshold must be non-negative")
        if args.anharmonic_wall_height_cm < 0.0:
            raise ValueError("--anharmonic-wall-height-cm must be non-negative")
        if args.anharmonic_handy_beta is not None and args.anharmonic_handy_beta <= 0.0:
            raise ValueError("--anharmonic-handy-beta must be positive")
        if args.anharmonic_kinetic_frequency_cm is not None and args.anharmonic_kinetic_frequency_cm <= 0.0:
            raise ValueError("--anharmonic-kinetic-frequency-cm must be positive")
    if args.boundary == "periodic":
        if extension_mode != "none":
            raise ValueError("--potential-extension is supported only for nonperiodic 1D paths")
        if args.potential_smoothing != "none":
            raise ValueError("--potential-smoothing is supported only for nonperiodic 1D paths")
        if args.core_model == "asymmetric-parabola-gaussian":
            raise ValueError("--core-model asymmetric-parabola-gaussian is supported only for nonperiodic 1D paths")
    if extension_mode != "none":
        if args.extension_length_au <= 0.0:
            raise ValueError("--extension-length-au must be positive with --potential-extension")
        if args.extension_points < 1:
            raise ValueError("--extension-points must be positive")
        if args.extension_degree not in (6, 8):
            raise ValueError("--extension-degree must be 6 or 8")
        if args.extension_target_cm <= 0.0:
            raise ValueError("--extension-target-cm must be positive")
        if args.extension_fit_points < 3:
            raise ValueError("--extension-fit-points must be at least 3")
    if extension_mode == "morse-polynomial" and args.path_symmetry == "none":
        raise ValueError("--potential-extension morse-polynomial requires --path-symmetry")
    if extension_mode == "morse-polynomial" and args.well_type == "single":
        raise ValueError("--potential-extension morse-polynomial requires --well-type double or auto")
    if extension_mode == "normal-mode-polynomial":
        if args.path_symmetry != "none":
            raise ValueError("--potential-extension normal-mode-polynomial currently requires --path-symmetry none")
        if not args.tail_normal_mode_log:
            raise ValueError("--potential-extension normal-mode-polynomial requires --tail-normal-mode-log")
        if args.tail_normal_mode_max_frequency <= 0.0:
            raise ValueError("--tail-normal-mode-max-frequency must be positive")
        if args.tail_normal_mode_min_overlap < 0.0:
            raise ValueError("--tail-normal-mode-min-overlap must be non-negative")
        if args.tail_normal_mode_steps < 1:
            raise ValueError("--tail-normal-mode-steps must be positive")
    if extension_mode in {"single-morse", "single-inverse-power"} and args.well_type == "double":
        raise ValueError(f"--potential-extension {extension_mode} requires --well-type single or auto")
    if args.core_model == "asymmetric-parabola-gaussian" and args.path_symmetry != "none":
        raise ValueError("--core-model asymmetric-parabola-gaussian requires --path-symmetry none")
    if args.core_model == "asymmetric-parabola-gaussian" and args.well_type == "single":
        raise ValueError("--core-model asymmetric-parabola-gaussian requires --well-type double or auto")
    if args.core_auto_max_points < 5:
        raise ValueError("--core-auto-max-points must be at least 5")
    if args.double_well_fit_points < 5:
        raise ValueError("--double-well-fit-points must be at least 5")
    if args.potential_spline_smoothing < 0.0:
        raise ValueError("--potential-spline-smoothing must be non-negative")
    if not args.grid2d_csv and args.metric_mode != "constant":
        raise ValueError("--metric-mode other than constant requires --grid2d-csv")
    if args.metric_stencil < 2:
        raise ValueError("--metric-stencil must be at least 2")
    if args.metric_smooth_order < 0:
        raise ValueError("--metric-smooth-order must be non-negative")
    if args.metric_eigenvalue_floor <= 0.0:
        raise ValueError("--metric-eigenvalue-floor must be positive")
    if args.metric_mode != "geometry" and (args.grid2d_geom_xyz or args.grid2d_geom_log):
        raise ValueError("--grid2d-geom-xyz/--grid2d-geom-log require --metric-mode geometry")
    if args.metric_mode != "csv" and (args.g11_key or args.g22_key or args.g12_key):
        raise ValueError("--g11-key/--g22-key/--g12-key require --metric-mode csv")
    if args.metric_mode == "geometry":
        sources = [args.grid2d_geom_xyz is not None, args.grid2d_geom_log is not None]
        if sum(sources) != 1:
            raise ValueError(
                "--metric-mode geometry needs exactly one of --grid2d-geom-xyz or --grid2d-geom-log"
            )


def run_2d_mode(args: argparse.Namespace) -> None:
    if args.convergence_scan:
        solve_2d_convergence_from_csv(args)
    else:
        solve_2d_hamiltonian_from_csv(args)


def run_1d_analysis(source: SourceData, args: argparse.Namespace) -> None:
    structures = source.structures
    input_path = source.path
    check_atom_order(structures)

    if args.properties_csv:
        merge_property_csv(structures, args.properties_csv.expanduser().resolve())

    masses = masses_for_atoms(structures[0].atoms, structures[0].symbols)
    if args.compute_rotconst:
        add_rotational_constants(structures, masses)

    oriented, s_sqrtamu_angstrom, angular_residuals = orient_path(structures, masses)
    step_lengths = np.zeros(len(structures), dtype=float)
    step_lengths[1:] = np.diff(s_sqrtamu_angstrom)
    s_au = s_sqrtamu_angstrom * MW_ANGSTROM_TO_AU

    if args.check_only:
        print(f"points={len(structures)}")
        print(f"path_length_sqrtamu_angstrom={s_sqrtamu_angstrom[-1]:.10f}")
        print(f"path_length_au={s_au[-1]:.10f}")
        print(f"max_angular_residual={np.max(angular_residuals):.6e}")
        return

    energy_key = find_energy_key(structures, args.energy_key)
    energy_unit = infer_energy_unit(energy_key, args.energy_unit)
    raw_energy, rel_energy_cm = potential_cm(structures, energy_key, args.energy_unit)
    property_derivative_specs = (
        read_property_derivative_specs(args.property_derivatives_csv.expanduser().resolve())
        if args.property_derivatives_csv
        else {}
    )
    property_keys = scalar_property_keys(
        structures,
        energy_key,
        args.property,
        derivative_property_names=property_derivative_specs,
    )
    cp_ring_indices = None
    if args.label_cremer_pople:
        parsed_ring = parse_ring_indices(args.ring, len(structures[0].atoms))
        if len(parsed_ring) >= 4:
            cp_ring_indices = parsed_ring
    gic_to_cp_bridges = fit_gic_to_cremer_pople_bridges(structures, cp_ring_indices)
    repeat = max(args.repeat, 1)
    n_grid = args.grid * repeat if args.boundary == "periodic" else args.grid

    property_samples = {
        key: np.array([structure.props[key] for structure in structures], dtype=float)
        for key in property_keys
        if key not in property_derivative_specs
    }
    selected_property_derivative_specs = {
        key: spec for key, spec in property_derivative_specs.items() if key in property_keys
    }
    normal_mode_tail = None
    if canonical_extension_mode(args.potential_extension) == "normal-mode-polynomial":
        normal_mode_tail = build_normal_mode_tail_data(
            args.tail_normal_mode_log.expanduser().resolve(),
            oriented,
            s_au,
            rel_energy_cm,
            masses,
            args,
        )
    grid_model = build_1d_grid_model(
        s_au,
        rel_energy_cm,
        property_samples,
        selected_property_derivative_specs,
        args,
        n_grid=n_grid,
        repeat=repeat,
        normal_mode_tail=normal_mode_tail,
    )

    levels, vectors, expectations, solver_used = solve_path_hamiltonian(
        grid_model.grid_au,
        grid_model.potential_cm,
        grid_model.expectation_properties,
        args,
    )
    args.solver_used = solver_used

    min_index = int(np.argmin(rel_energy_cm))
    reference_values = {
        key: float(structures[min_index].props[key])
        for key in property_keys
        if key in structures[min_index].props
    }
    for key, spec in selected_property_derivative_specs.items():
        reference_values.setdefault(key, float(spec.value))
    automatic_vpt2 = PropertyVPT2Result({})
    if property_keys and args.property_vpt2 == "auto":
        automatic_vpt2 = automatic_property_vpt2(
            grid_model.model_points_au,
            grid_model.model_potential_cm,
            {
                key: grid_model.model_properties[key]
                for key in property_keys
                if key in grid_model.model_properties
            },
            reference_values,
            fit_points=args.property_vpt2_fit_points,
            degree=args.property_vpt2_degree,
            basis_size=args.property_vpt2_basis_size,
        )
        grid_model.info.update(automatic_vpt2.info)
    external_vpt2_property_comparison = (
        read_vpt2_property_csv(
            args.vpt2_property_csv.expanduser().resolve(),
            property_keys,
            reference_values,
        )
        if args.vpt2_property_csv and args.property_vpt2 != "off"
        else {}
    )
    vpt2_property_comparison = (
        merge_property_vpt2_comparisons(
            automatic_vpt2.comparisons if args.property_vpt2 == "auto" else {},
            external_vpt2_property_comparison,
            property_keys,
            reference_values,
        )
        if property_keys and args.property_vpt2 != "off"
        else {}
    )

    ensure_output_dirs(args.outdir, args.figdir)
    prefix = args.prefix
    write_profile(
        args.outdir / f"{prefix}_profile.csv",
        structures,
        s_sqrtamu_angstrom,
        s_au,
        step_lengths,
        angular_residuals,
        raw_energy,
        rel_energy_cm,
        energy_key,
        property_keys,
        cp_ring_indices,
        gic_to_cp_bridges,
    )
    write_oriented_xyz(
        args.outdir / f"{prefix}_oriented.xyz",
        structures,
        oriented,
        s_sqrtamu_angstrom,
        s_au,
        rel_energy_cm,
        cp_ring_indices,
    )
    write_grid(
        args.outdir / f"{prefix}_grid.csv",
        grid_model.grid_au,
        grid_model.potential_cm,
        grid_model.properties,
        vectors,
        args.save_states,
    )
    write_1d_model_profile(args.outdir / f"{prefix}_model_profile.csv", grid_model)
    write_levels(args.outdir / f"{prefix}_levels.csv", levels)
    write_expectations(args.outdir / f"{prefix}_expectations.csv", levels, expectations)
    if property_keys:
        write_thermal_expectations(
            args.outdir / f"{prefix}_thermal_expectations.csv",
            levels,
            expectations,
            args.temperature,
        )
    if vpt2_property_comparison:
        write_vpt2_property_correction(
            args.outdir / f"{prefix}_property_comparison.csv",
            vpt2_property_comparison,
            expectations[0],
        )
        write_vpt2_property_correction(
            args.outdir / f"{prefix}_vpt2_property_correction.csv",
            vpt2_property_comparison,
            expectations[0],
        )
    write_summary(
        args.outdir / f"{prefix}_summary.txt",
        args,
        n_grid,
        input_path,
        structures,
        s_sqrtamu_angstrom,
        s_au,
        angular_residuals,
        energy_key,
        energy_unit,
        property_keys,
        levels,
        model_info=grid_model.info,
    )
    plot_results(
        args.figdir / f"{prefix}_potential_levels.pdf",
        args.figdir / f"{prefix}_potential_levels.png",
        grid_model.grid_au,
        grid_model.potential_cm,
        grid_model.properties,
        expectations,
        levels,
        vectors,
        states=min(args.save_states, args.levels),
        plot_max_state=args.plot_max_state,
        property_smooth_degree=args.plot_property_smooth_degree,
        vpt2_properties=vpt2_property_comparison,
    )
    plot_potential_profile(
        args.figdir / f"{prefix}_potential_profile.pdf",
        args.figdir / f"{prefix}_potential_profile.png",
        rel_energy_cm,
        period_deg=180.0,
    )
    plot_circular_potential(
        args.figdir / f"{prefix}_circular_potential.pdf",
        args.figdir / f"{prefix}_circular_potential.png",
        grid_model.grid_au,
        grid_model.potential_cm,
        levels,
        states=min(args.save_states, args.levels),
    )

    print(args.outdir / f"{prefix}_summary.txt")
    print(args.outdir / f"{prefix}_levels.csv")
    print(args.outdir / f"{prefix}_model_profile.csv")
    if property_keys:
        print(args.outdir / f"{prefix}_expectations.csv")
        print(args.outdir / f"{prefix}_thermal_expectations.csv")
    if vpt2_property_comparison:
        print(args.outdir / f"{prefix}_property_comparison.csv")
        print(args.outdir / f"{prefix}_vpt2_property_correction.csv")
    print(args.figdir / f"{prefix}_potential_levels.pdf")
    print(args.figdir / f"{prefix}_potential_profile.pdf")
    print(args.figdir / f"{prefix}_circular_potential.pdf")


def run_anharmonic_mode_analysis(args: argparse.Namespace) -> None:
    input_path = args.gaussian_log.expanduser().resolve()
    derivatives = parse_gaussian_anharmonic_derivatives(
        input_path,
        args.anharmonic_mode,
        args.anharmonic_mode_order,
    )
    potential = build_anharmonic_derivative_potential(derivatives, args)
    kinetic_frequency_cm = (
        float(args.anharmonic_kinetic_frequency_cm)
        if args.anharmonic_kinetic_frequency_cm is not None
        else abs(float(derivatives.f2_cm))
    )
    if kinetic_frequency_cm <= 0.0:
        raise ValueError(
            "The kinetic prefactor is zero; pass --anharmonic-kinetic-frequency-cm"
        )
    levels, vectors = solve_anharmonic_derivative_dvr(
        potential.q,
        potential.potential_cm,
        kinetic_frequency_cm,
        args.levels,
    )
    vpt2_comparison = None
    if potential.info.get("resolved_well_type") != "double" and derivatives.f2_cm > 0.0:
        vpt2_comparison = build_one_mode_vpt2_comparison(derivatives, potential, levels)
    args.solver_used = "sinc-dvr-normal-coordinate"

    ensure_output_dirs(args.outdir, args.figdir)
    prefix = args.prefix
    write_anharmonic_derivative_grid(
        args.outdir / f"{prefix}_anharmonic_grid.csv",
        potential.q,
        potential.potential_cm,
        potential.model_potential_cm,
        potential.taylor_potential_cm,
        vectors,
        args.save_states,
    )
    write_levels(args.outdir / f"{prefix}_anharmonic_levels.csv", levels)
    write_anharmonic_derivative_summary(
        args.outdir / f"{prefix}_anharmonic_summary.txt",
        args,
        input_path,
        derivatives,
        potential,
        kinetic_frequency_cm,
        levels,
        vpt2_comparison=vpt2_comparison,
    )
    if vpt2_comparison is not None:
        write_anharmonic_vpt2_comparison(
            args.outdir / f"{prefix}_anharmonic_vpt2_comparison.csv",
            vpt2_comparison,
        )
    plot_anharmonic_derivative_results(
        args.figdir / f"{prefix}_anharmonic_potential_levels.pdf",
        args.figdir / f"{prefix}_anharmonic_potential_levels.png",
        potential.q,
        potential.potential_cm,
        levels,
        vectors,
        states=min(args.save_states, args.levels),
    )

    print(args.outdir / f"{prefix}_anharmonic_summary.txt")
    print(args.outdir / f"{prefix}_anharmonic_levels.csv")
    print(args.outdir / f"{prefix}_anharmonic_grid.csv")
    if vpt2_comparison is not None:
        print(args.outdir / f"{prefix}_anharmonic_vpt2_comparison.csv")
    print(args.figdir / f"{prefix}_anharmonic_potential_levels.pdf")


def run_1d_mode(args: argparse.Namespace) -> None:
    if args.anharmonic_mode is not None:
        run_anharmonic_mode_analysis(args)
        return
    source = read_structures(args)
    run_1d_analysis(source, args)


def main() -> None:
    args = parse_args()
    try:
        validate_args(args)
    except ValueError as exc:
        raise SystemExit(f"error: {exc}") from None
    if args.grid2d_csv:
        run_2d_mode(args)
    else:
        run_1d_mode(args)


if __name__ == "__main__":
    main()
