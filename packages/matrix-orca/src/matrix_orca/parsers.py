from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

from matrix_chem import MolecularGeometry, atomic_mass
from matrix_chem.geometry_io import GeometryParseError, normalize_atom_symbol
from matrix_chem.topology.elements import atomic_number


ORCA_OUTPUT_FORMAT = "orca_output"
ANGSTROM_TO_BOHR = 1.0 / 0.52917721092
BOHR_TO_ANGSTROM = 0.52917721092
FINAL_ENERGY_RE = re.compile(
    r"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+([-+]?(?:\d+\.\d*|\.\d+)(?:[DEde][-+]?\d+)?)",
    re.I,
)
FREQUENCY_RE = re.compile(
    r"^\s*\d+\s*:\s*([-+]?(?:\d+\.\d*|\.\d+)(?:[DEde][-+]?\d+)?)\s*cm\*\*-1",
    re.I,
)
CHARGE_RE = re.compile(r"\bTotal\s+Charge\b.*?([-+]?\d+)\s*$", re.I)
MULT_RE = re.compile(r"\bMultiplicity\b.*?(\d+)\s*$", re.I)
XYZ_INPUT_RE = re.compile(r"^\s*\*\s+xyz(?:file)?\s+([-+]?\d+)\s+(\d+)\b", re.I)


@dataclass(frozen=True)
class OrcaOutputSummary:
    path: Path
    geometry: MolecularGeometry
    charge: int
    multiplicity: int
    final_energy_hartree: float | None
    frequencies_cm: tuple[float, ...]
    cartesian_hessian: np.ndarray | None
    cartesian_coordinate_blocks: int
    normal_termination: bool


@dataclass(frozen=True)
class OrcaOutputPromotion:
    xyzin: Path
    output_path: Path
    wrote_geometry: bool
    wrote_cartesian_hessian: bool


@dataclass(frozen=True)
class OrcaMoldenPromotion:
    xyzin: Path
    gbw_path: Path
    molden_path: Path
    command: tuple[str, ...]
    wrote_orbitals: bool


def summarize_orca_output(path: Path | str) -> OrcaOutputSummary:
    target = Path(path)
    text = target.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    geometry = _parse_last_cartesian_coordinates(lines, source_path=target)
    charge, multiplicity = _parse_charge_multiplicity(lines)
    frequencies = tuple(_parse_frequencies(lines))
    hessian = _parse_cartesian_hessian(lines, size=3 * geometry.natoms)
    if hessian is not None and not np.allclose(hessian, hessian.T):
        hessian = 0.5 * (hessian + hessian.T)
    return OrcaOutputSummary(
        path=target,
        geometry=geometry,
        charge=charge,
        multiplicity=multiplicity,
        final_energy_hartree=_parse_final_energy(text),
        frequencies_cm=frequencies,
        cartesian_hessian=hessian,
        cartesian_coordinate_blocks=_coordinate_block_count(lines),
        normal_termination="ORCA TERMINATED NORMALLY" in text.upper(),
    )


def read_orca_output_geometry(path: Path | str) -> MolecularGeometry:
    return summarize_orca_output(path).geometry


def hessian_input_from_orca_output(path: Path | str):
    from matrix_gf import HessianInput

    summary = summarize_orca_output(path)
    if summary.cartesian_hessian is None:
        raise GeometryParseError("ORCA output contains no readable Cartesian Hessian")
    geometry = summary.geometry
    numbers = np.asarray(_atomic_numbers_from_geometry(geometry), dtype=int)
    masses = np.asarray([atomic_mass(int(number)) for number in numbers], dtype=float)
    data = HessianInput(
        atomic_numbers=numbers,
        cartesian_coordinates_bohr=np.asarray(geometry.coordinates_angstrom, dtype=float)
        * ANGSTROM_TO_BOHR,
        masses_amu=masses,
        cartesian_hessian=np.asarray(summary.cartesian_hessian, dtype=float),
        harmonic_frequencies_cm=np.asarray(summary.frequencies_cm, dtype=float),
        source="orca-output",
    )
    data.validate()
    return data


def promote_orca_output_to_xyzin(
    source: Path | str,
    xyzin: Path | str,
    *,
    symmetry_distance: float = 1.0e-3,
    symmetry_inertia: float = 1.0e-3,
    max_rotation_order: int = 6,
) -> OrcaOutputPromotion:
    from matrix_chem import SymmetryThresholds, preprocess_to_enriched_xyz
    from matrix_qm import cartesian_hessian_section_from_hessian_input, write_cartesian_hessian_section

    output = Path(source)
    target = Path(xyzin)
    preprocess_to_enriched_xyz(
        output,
        target,
        source_kind="orca",
        symmetry_thresholds=SymmetryThresholds(
            distance_angstrom=symmetry_distance,
            inertia_relative=symmetry_inertia,
            max_rotation_order=max_rotation_order,
        ),
    )
    wrote_hessian = False
    try:
        hessian_input = hessian_input_from_orca_output(output)
    except GeometryParseError:
        pass
    else:
        write_cartesian_hessian_section(
            target,
            cartesian_hessian_section_from_hessian_input(hessian_input, source="orca-output"),
        )
        wrote_hessian = True
    return OrcaOutputPromotion(
        xyzin=target,
        output_path=output,
        wrote_geometry=True,
        wrote_cartesian_hessian=wrote_hessian,
    )


def convert_orca_gbw_to_molden(
    gbw: Path | str,
    *,
    output: Path | str | None = None,
    executable: str = "orca_2mkl",
    timeout: float | None = None,
) -> OrcaMoldenPromotion:
    """Run ORCA's orca_2mkl converter and return the generated Molden path.

    The returned object has an empty xyzin path because this helper only creates
    the external file. Use :func:`promote_orca_molden_to_xyzin` to register it.
    """
    source = Path(gbw)
    if source.suffix.lower() != ".gbw":
        raise ValueError("ORCA Molden conversion requires a .gbw file")
    if not source.is_file():
        raise FileNotFoundError(source)
    basename = source.with_suffix("")
    generated = basename.with_suffix(".molden.input")
    command = (executable, str(basename), "-molden")
    subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if not generated.is_file():
        raise FileNotFoundError(f"orca_2mkl did not create expected Molden file: {generated}")
    target = Path(output) if output is not None else generated
    if target != generated:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(generated, target)
    return OrcaMoldenPromotion(
        xyzin=Path(),
        gbw_path=source,
        molden_path=target,
        command=command,
        wrote_orbitals=False,
    )


def promote_orca_molden_to_xyzin(
    gbw: Path | str,
    xyzin: Path | str,
    *,
    output: Path | str | None = None,
    executable: str = "orca_2mkl",
    timeout: float | None = None,
) -> OrcaMoldenPromotion:
    """Convert an ORCA GBW file to Molden and register it in #ORBITALS."""
    from matrix_qm import merge_orbitals_section, orbital_file_record_from_path

    target = Path(xyzin)
    converted = convert_orca_gbw_to_molden(
        gbw,
        output=output,
        executable=executable,
        timeout=timeout,
    )
    merge_orbitals_section(
        target,
        (
            orbital_file_record_from_path(
                converted.molden_path,
                role="orbitals",
                label=converted.molden_path.stem,
                source="orca_2mkl",
            ),
        ),
    )
    return OrcaMoldenPromotion(
        xyzin=target,
        gbw_path=converted.gbw_path,
        molden_path=converted.molden_path,
        command=converted.command,
        wrote_orbitals=True,
    )


def _parse_last_cartesian_coordinates(
    lines: list[str],
    *,
    source_path: Path,
) -> MolecularGeometry:
    block = _last_coordinate_block(lines, unit="angstrom")
    unit = "angstrom"
    if block is None:
        block = _last_coordinate_block(lines, unit="bohr")
        unit = "bohr"
    if block is None:
        raise GeometryParseError("ORCA output contains no Cartesian coordinate block")
    atoms, coords = block
    coordinates = np.asarray(coords, dtype=float)
    if unit == "bohr":
        coordinates = coordinates * BOHR_TO_ANGSTROM
    charge, multiplicity = _parse_charge_multiplicity(lines)
    return MolecularGeometry(
        atoms=tuple(atoms),
        coordinates_angstrom=coordinates,
        comment=source_path.name,
        source_format=ORCA_OUTPUT_FORMAT,
        source_path=source_path,
        charge=charge,
        multiplicity=multiplicity,
        metadata={"coordinate_block": f"CARTESIAN COORDINATES ({unit.upper()})"},
    )


def _last_coordinate_block(
    lines: list[str],
    *,
    unit: str,
) -> tuple[list[str], list[list[float]]] | None:
    needle = (
        "CARTESIAN COORDINATES (ANGSTROEM)"
        if unit == "angstrom"
        else "CARTESIAN COORDINATES (A.U.)"
    )
    starts = [idx for idx, line in enumerate(lines) if needle in line.upper()]
    for start in reversed(starts):
        atoms: list[str] = []
        coords: list[list[float]] = []
        for raw in lines[start + 1 :]:
            text = raw.strip()
            if not text or set(text) <= {"-"}:
                if atoms:
                    break
                continue
            parsed = _parse_coordinate_line(raw)
            if parsed is None:
                if atoms:
                    break
                continue
            atom, xyz = parsed
            atoms.append(atom)
            coords.append(xyz)
        if atoms:
            return atoms, coords
    return None


def _parse_coordinate_line(line: str) -> tuple[str, list[float]] | None:
    parts = line.split()
    if len(parts) < 4:
        return None
    try:
        atom = normalize_atom_symbol(parts[0])
    except GeometryParseError:
        return None
    numeric = [_float_or_none(token) for token in parts[1:]]
    values = [value for value in numeric if value is not None]
    if len(values) < 3:
        return None
    return atom, [float(value) for value in values[-3:]]


def _parse_charge_multiplicity(lines: list[str]) -> tuple[int, int]:
    charge = 0
    multiplicity = 1
    for line in lines:
        charge_match = CHARGE_RE.search(line)
        if charge_match is not None:
            charge = int(charge_match.group(1))
        mult_match = MULT_RE.search(line)
        if mult_match is not None:
            multiplicity = int(mult_match.group(1))
        xyz_match = XYZ_INPUT_RE.search(line)
        if xyz_match is not None:
            charge = int(xyz_match.group(1))
            multiplicity = int(xyz_match.group(2))
    return charge, multiplicity


def _parse_final_energy(text: str) -> float | None:
    matches = tuple(FINAL_ENERGY_RE.finditer(text))
    if not matches:
        return None
    return float(matches[-1].group(1).replace("D", "E").replace("d", "e"))


def _parse_frequencies(lines: list[str]) -> list[float]:
    values: list[float] = []
    in_block = False
    for line in lines:
        upper = line.upper()
        if "VIBRATIONAL FREQUENCIES" in upper:
            values = []
            in_block = True
            continue
        if not in_block:
            continue
        match = FREQUENCY_RE.match(line)
        if match is not None:
            values.append(float(match.group(1).replace("D", "E").replace("d", "e")))
            continue
        if values and upper.strip() and upper == upper.upper() and ":" not in line:
            break
    return values


def _parse_cartesian_hessian(lines: list[str], *, size: int) -> np.ndarray | None:
    starts = [
        idx
        for idx, line in enumerate(lines)
        if "CARTESIAN FORCE CONSTANT MATRIX" in line.upper()
        or "CARTESIAN HESSIAN" in line.upper()
    ]
    for start in reversed(starts):
        matrix = _parse_block_matrix(lines[start + 1 :], size=size)
        if matrix is not None:
            return matrix
    return None


def _parse_block_matrix(lines: list[str], *, size: int) -> np.ndarray | None:
    matrix = np.zeros((size, size), dtype=float)
    filled: set[tuple[int, int]] = set()
    columns: list[int] = []
    for raw in lines:
        text = raw.strip()
        if not text or set(text) <= {"-"}:
            continue
        parts = text.split()
        int_tokens = [_int_or_none(token) for token in parts]
        if int_tokens and all(value is not None for value in int_tokens):
            columns = [int(value) for value in int_tokens if value is not None]
            continue
        row = _int_or_none(parts[0]) if parts else None
        if row is None:
            if filled and text.upper() == text and not _contains_float(text):
                break
            continue
        values = [_float_or_none(token) for token in parts[1:]]
        numeric = [value for value in values if value is not None]
        if not numeric:
            continue
        if columns and len(numeric) <= len(columns):
            selected_columns = columns[: len(numeric)]
        elif len(numeric) == size:
            selected_columns = list(range(size))
        else:
            continue
        if not 0 <= row < size:
            continue
        for column, value in zip(selected_columns, numeric):
            if 0 <= column < size:
                matrix[row, column] = float(value)
                filled.add((row, column))
        if len(filled) >= size * size:
            return matrix
    if not filled:
        return None
    if _looks_lower_triangular(filled, size):
        matrix = matrix + matrix.T - np.diag(np.diag(matrix))
    return matrix


def _looks_lower_triangular(filled: set[tuple[int, int]], size: int) -> bool:
    expected = {(row, col) for row in range(size) for col in range(row + 1)}
    return expected.issubset(filled) and not any(col > row for row, col in filled)


def _coordinate_block_count(lines: list[str]) -> int:
    return sum(
        1
        for line in lines
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line.upper()
        or "CARTESIAN COORDINATES (A.U.)" in line.upper()
    )


def _atomic_numbers_from_geometry(geometry: MolecularGeometry) -> tuple[int, ...]:
    numbers: list[int] = []
    for atom in geometry.atoms:
        number = atomic_number(atom)
        if number is None or number <= 0:
            raise GeometryParseError(f"unknown element symbol in ORCA output: {atom}")
        numbers.append(int(number))
    return tuple(numbers)


def _float_or_none(token: str) -> float | None:
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def _int_or_none(token: str) -> int | None:
    try:
        return int(token)
    except ValueError:
        return None


def _contains_float(text: str) -> bool:
    return any(_float_or_none(token) is not None for token in text.split())
