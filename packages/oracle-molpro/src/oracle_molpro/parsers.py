from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from oracle_chem import MolecularGeometry
from oracle_chem.geometry_io import GeometryParseError, normalize_atom_symbol


MOLPRO_OUTPUT_FORMAT = "molpro_output"


@dataclass(frozen=True)
class MolproOutputSummary:
    path: Path
    geometry: MolecularGeometry
    charge: int
    multiplicity: int
    atomic_coordinate_blocks: int


def summarize_molpro_output(path: Path | str) -> MolproOutputSummary:
    target = Path(path)
    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms, coords = _parse_atomic_coordinates(lines)
    charge, multiplicity = _parse_charge_multiplicity(lines)
    geometry = MolecularGeometry(
        atoms=tuple(atoms),
        coordinates_angstrom=np.asarray(coords, dtype=float),
        comment=target.name,
        source_format=MOLPRO_OUTPUT_FORMAT,
        source_path=target,
        charge=charge,
        multiplicity=multiplicity,
        metadata={"coordinate_block": "ATOMIC COORDINATES"},
    )
    return MolproOutputSummary(
        path=target,
        geometry=geometry,
        charge=charge,
        multiplicity=multiplicity,
        atomic_coordinate_blocks=sum(1 for line in lines if "ATOMIC COORDINATES" in line.upper()),
    )


def read_molpro_output_geometry(path: Path | str) -> MolecularGeometry:
    return summarize_molpro_output(path).geometry


def promote_molpro_output_to_xyzin(
    source: Path | str,
    xyzin: Path | str,
    *,
    symmetry_distance: float = 1.0e-3,
    symmetry_inertia: float = 1.0e-3,
    max_rotation_order: int = 6,
):
    from oracle_chem import SymmetryThresholds, preprocess_to_enriched_xyz

    return preprocess_to_enriched_xyz(
        Path(source),
        Path(xyzin),
        source_kind="molpro",
        symmetry_thresholds=SymmetryThresholds(
            distance_angstrom=symmetry_distance,
            inertia_relative=symmetry_inertia,
            max_rotation_order=max_rotation_order,
        ),
    )


def _parse_atomic_coordinates(lines: list[str]) -> tuple[list[str], list[list[float]]]:
    indices = [idx for idx, line in enumerate(lines) if "ATOMIC COORDINATES" in line.upper()]
    if not indices:
        raise GeometryParseError("Molpro output contains no ATOMIC COORDINATES section")

    idx = indices[-1] + 1
    while idx < len(lines) and not _coordinate_candidate(lines[idx]):
        idx += 1

    atoms: list[str] = []
    coords: list[list[float]] = []
    while idx < len(lines):
        parsed = _parse_coordinate_line(lines[idx])
        if parsed is None:
            if atoms:
                break
            idx += 1
            continue
        atom, xyz = parsed
        atoms.append(atom)
        coords.append(xyz)
        idx += 1

    if not atoms:
        raise GeometryParseError("Molpro ATOMIC COORDINATES section contains no atoms")
    return atoms, coords


def _parse_charge_multiplicity(lines: list[str]) -> tuple[int, int]:
    charge = 0
    multiplicity = 1
    charge_pattern = re.compile(r"\bCHARGE\b\s*(?:=|:)?\s*([-+]?\d+)", re.I)
    spin_pattern = re.compile(r"SPIN\s+QUANTUM\s+NUMBER\s*(?:=|:)?\s*([-+]?\d*\.?\d+)", re.I)
    mult_pattern = re.compile(r"\bMULTIPLICITY\b\s*(?:=|:)?\s*(\d+)", re.I)
    for line in lines:
        charge_match = charge_pattern.search(line)
        if charge_match is not None:
            charge = int(charge_match.group(1))
        mult_match = mult_pattern.search(line)
        if mult_match is not None:
            multiplicity = int(mult_match.group(1))
        spin_match = spin_pattern.search(line)
        if spin_match is not None:
            multiplicity = max(1, int(round(2.0 * float(spin_match.group(1)) + 1.0)))
    return charge, multiplicity


def _coordinate_candidate(line: str) -> bool:
    return _parse_coordinate_line(line) is not None


def _parse_coordinate_line(line: str) -> tuple[str, list[float]] | None:
    parts = line.split()
    if len(parts) < 4:
        return None
    atom_token = parts[1] if parts[0].isdigit() and len(parts) >= 5 else parts[0]
    try:
        atom = normalize_atom_symbol(atom_token)
    except GeometryParseError:
        return None
    numeric_values: list[float] = []
    for token in parts:
        value = _float_or_none(token)
        if value is not None:
            numeric_values.append(value)
    if len(numeric_values) < 3:
        return None
    return atom, [float(value) for value in numeric_values[-3:]]


def _float_or_none(token: str) -> float | None:
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None
