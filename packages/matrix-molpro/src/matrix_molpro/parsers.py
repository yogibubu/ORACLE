from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from matrix_chem import MolecularGeometry
from matrix_chem.geometry_io import GeometryParseError, normalize_atom_symbol


MOLPRO_OUTPUT_FORMAT = "molpro_output"


@dataclass(frozen=True)
class MolproOutputSummary:
    path: Path
    geometry: MolecularGeometry
    charge: int
    multiplicity: int
    atomic_coordinate_blocks: int


@dataclass(frozen=True)
class MolproQuadrupolePromotion:
    xyzin: Path
    output_path: Path
    wrote_properties: bool
    property_count: int


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
    from matrix_chem import SymmetryThresholds, preprocess_to_enriched_xyz

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


def parse_molpro_quadrupole_properties(
    path: Path | str,
    *,
    atom: int | None = None,
    isotope: str = "",
) -> tuple:
    from matrix_qm import infer_unique_quadrupolar_atom, quadrupole_property_records_from_efg

    target = Path(path)
    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
    method, relaxation, components = _parse_last_efg_components(lines)
    if not components:
        return ()
    input_atoms = _parse_input_geometry_atoms(lines)
    try:
        atoms, _ = _parse_atomic_coordinates(lines)
    except GeometryParseError:
        atoms = input_atoms
    selected_atom = atom
    selected_number = None
    if selected_atom is not None:
        atom_source = input_atoms if input_atoms and selected_atom <= len(input_atoms) else atoms
        if selected_atom < 1 or selected_atom > len(atom_source):
            raise ValueError(f"Molpro quadrupole atom index out of range: {selected_atom}")
        from matrix_chem.topology.elements import atomic_number

        selected_number = atomic_number(atom_source[selected_atom - 1])
    else:
        inferred = infer_unique_quadrupolar_atom(tuple(atoms))
        if inferred is None and input_atoms:
            inferred = infer_unique_quadrupolar_atom(tuple(input_atoms))
        if inferred is not None:
            selected_atom, selected_number = inferred
    if selected_number is None:
        if not isotope:
            raise ValueError(
                "Molpro EFG output does not identify the nucleus; provide --atom or --isotope"
            )
        from matrix_qm import parse_isotope_label

        selected_number, _, _ = parse_isotope_label(isotope)
    level = _parse_basis(lines)
    efg = (
        components.get("XX", 0.0),
        components.get("YY", 0.0),
        components.get("ZZ", 0.0),
        components.get("XY", 0.0),
        components.get("XZ", 0.0),
        components.get("YZ", 0.0),
    )
    return quadrupole_property_records_from_efg(
        atom=selected_atom,
        atomic_number_value=int(selected_number),
        efg_au=efg,
        program="Molpro",
        source=target,
        isotope=isotope,
        method=method,
        level=level,
        axes="MOLPRO_EFG:xx,yy,zz,xy,xz,yz",
        status="raw",
        comment=f"relaxation={relaxation}",
    )


def promote_molpro_quadrupole_properties_to_xyzin(
    output: Path | str,
    xyzin: Path | str,
    *,
    atom: int | None = None,
    isotope: str = "",
) -> MolproQuadrupolePromotion:
    from matrix_qm import merge_properties_section

    source = Path(output)
    target = Path(xyzin)
    records = parse_molpro_quadrupole_properties(source, atom=atom, isotope=isotope)
    if records:
        merge_properties_section(target, records)
    return MolproQuadrupolePromotion(
        xyzin=target,
        output_path=source,
        wrote_properties=bool(records),
        property_count=len(records),
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


def _parse_input_geometry_atoms(lines: list[str]) -> list[str]:
    starts = [idx for idx, line in enumerate(lines) if line.strip().lower().startswith("geometry={")]
    for start in reversed(starts):
        atoms: list[str] = []
        for raw in lines[start + 1 :]:
            text = raw.strip()
            if not text:
                continue
            if text.startswith("}"):
                break
            parts = text.split()
            if len(parts) == 1 and parts[0].isdigit():
                continue
            try:
                atoms.append(normalize_atom_symbol(parts[0]))
            except (GeometryParseError, IndexError):
                continue
        if atoms:
            return atoms
    return []


def _float_or_none(token: str) -> float | None:
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def _parse_last_efg_components(lines: list[str]) -> tuple[str, str, dict[str, float]]:
    groups: dict[str, tuple[str, dict[str, float]]] = {}
    for raw in lines:
        upper = raw.upper()
        if "FG" not in upper:
            continue
        component = _efg_component_from_line(upper)
        if component is None:
            continue
        value = _last_float(raw)
        if value is None:
            continue
        method = _method_from_efg_line(raw)
        relaxation = _relaxation_from_efg_line(upper)
        key = relaxation or "unknown"
        current_method, current_components = groups.get(key, (method, {}))
        if method:
            current_method = method
        current_components[component] = value
        groups[key] = (current_method, current_components)
    for key in ("relax", "norela", "unknown"):
        if key in groups:
            method, components = groups[key]
            if {"XX", "YY", "ZZ"}.issubset(components):
                return method, key, components
    return "", "", {}


def _efg_component_from_line(line: str) -> str | None:
    match = re.search(r"\bFG([XYZ][XYZ])\b", line)
    if match is None:
        return None
    component = match.group(1)
    if component in {"YX", "ZX", "ZY"}:
        component = component[::-1]
    return component


def _method_from_efg_line(line: str) -> str:
    match = re.search(r"!\s*([A-Za-z0-9()+\-]+)\s+\(", line)
    return "" if match is None else match.group(1)


def _relaxation_from_efg_line(line: str) -> str:
    if "ED,RELAX" in line:
        return "relax"
    if "ED,NORELA" in line:
        return "norela"
    return ""


def _last_float(line: str) -> float | None:
    for token in reversed(line.replace("=", " ").split()):
        value = _float_or_none(token)
        if value is not None:
            return value
    return None


def _parse_basis(lines: list[str]) -> str:
    basis_pattern = re.compile(r"\bBASIS\s*=\s*([A-Za-z0-9()+.,_\-/]+)", re.I)
    for raw in lines:
        match = basis_pattern.search(raw)
        if match is not None:
            return match.group(1)
    return ""
