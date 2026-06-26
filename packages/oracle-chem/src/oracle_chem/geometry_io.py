from __future__ import annotations

from pathlib import Path
import re

import numpy as np

from .geometry import MolecularGeometry
from .topology.elements import atomic_number, atomic_symbol


class GeometryParseError(ValueError):
    """Raised when a geometry source cannot be parsed into ORACLE geometry."""


def normalize_atom_symbol(token: str) -> str:
    cleaned = str(token).split("(", 1)[0].split("-", 1)[0].strip()
    if not cleaned:
        raise GeometryParseError("empty atom token")
    if cleaned.isdigit():
        number = int(cleaned)
        symbol = atomic_symbol(number)
        if number <= 0 or symbol == "??":
            raise GeometryParseError(f"invalid atomic number: {token}")
        return symbol
    match = re.match(r"([A-Za-z]{1,3})", cleaned)
    if not match:
        raise GeometryParseError(f"invalid atom token: {token}")
    text = match.group(1)
    symbol = text[0].upper() + text[1:].lower()
    number = atomic_number(symbol)
    if number is None or number <= 0:
        raise GeometryParseError(f"invalid atom token: {token}")
    return symbol


def parse_xyz_lines(
    lines: list[str],
    *,
    source_path: Path | None = None,
    source_format: str = "xyz",
) -> MolecularGeometry:
    if len(lines) < 2:
        raise GeometryParseError("XYZ input has fewer than two lines")
    try:
        natoms = int(lines[0].strip())
    except ValueError as exc:
        raise GeometryParseError("first XYZ line must be an atom count") from exc
    if natoms < 0:
        raise GeometryParseError("XYZ atom count cannot be negative")

    comment = lines[1].rstrip()
    atoms: list[str] = []
    coords: list[list[float]] = []
    for raw in lines[2:]:
        if not raw.strip():
            continue
        parts = raw.split()
        if len(parts) < 4:
            continue
        atom = normalize_atom_symbol(parts[0])
        try:
            xyz = [float(parts[1]), float(parts[2]), float(parts[3])]
        except ValueError as exc:
            raise GeometryParseError(f"invalid XYZ coordinate line: {raw}") from exc
        atoms.append(atom)
        coords.append(xyz)
        if len(atoms) == natoms:
            break

    if len(atoms) != natoms:
        raise GeometryParseError(f"expected {natoms} atoms, found {len(atoms)}")

    return MolecularGeometry(
        atoms=tuple(atoms),
        coordinates_angstrom=np.asarray(coords, dtype=float),
        comment=comment,
        source_format=source_format,
        source_path=source_path,
    )


def read_xyz(path: Path) -> MolecularGeometry:
    target = Path(path)
    return parse_xyz_lines(
        target.read_text(encoding="utf-8", errors="replace").splitlines(),
        source_path=target,
        source_format="xyz",
    )


def read_enriched_xyz(path: Path) -> MolecularGeometry:
    target = Path(path)
    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
    geometry = parse_xyz_lines(lines, source_path=target, source_format="enriched_xyz")
    section_names = tuple(
        line.strip()[1:].strip().upper()
        for line in lines[2 + geometry.natoms :]
        if line.strip().startswith("#")
    )
    return MolecularGeometry(
        atoms=geometry.atoms,
        coordinates_angstrom=geometry.coordinates_angstrom,
        comment=geometry.comment,
        source_format="enriched_xyz",
        source_path=target,
        metadata={"sections": section_names},
    )


def read_geometry(path: Path) -> MolecularGeometry:
    target = Path(path)
    suffix = target.suffix.lower()
    if target.name.lower() == "xyzin":
        return read_enriched_xyz(target)
    if suffix == ".xyz" or suffix == "":
        return read_xyz(target)
    if suffix in {".zmat", ".zmt"}:
        from .zmatrix import read_zmatrix

        return read_zmatrix(target)
    if suffix in {".gjf", ".gau", ".com", ".inp"}:
        from oracle_babel import is_legacy_smiles_input, read_legacy_smiles_input
        from oracle_gaussian import read_gaussian_input

        if is_legacy_smiles_input(target):
            return read_legacy_smiles_input(target)
        return read_gaussian_input(target)
    raise GeometryParseError(f"unsupported geometry format: {target}")
