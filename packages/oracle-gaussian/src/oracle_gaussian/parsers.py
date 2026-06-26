from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from oracle_chem import MolecularGeometry
from oracle_chem.geometry_io import GeometryParseError, normalize_atom_symbol


SCF_RE = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s+([-+]?\d+\.\d+)")
NORMAL_TERMINATION = "Normal termination of Gaussian"
STANDARD_ORIENTATION = "Standard orientation:"
INPUT_ORIENTATION = "Input orientation:"


@dataclass(frozen=True)
class GaussianLogSummary:
    path: Path
    normal_termination: bool
    scf_energies_hartree: tuple[float, ...]
    standard_orientation_count: int
    input_orientation_count: int
    scan_marker_count: int
    puckering_marker_count: int
    frequencies_cm: tuple[float, ...] = ()
    last_orientation: MolecularGeometry | None = None


def read_gaussian_cartesian_input(path: Path) -> MolecularGeometry:
    target = Path(path)
    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
    route_end = _route_end(lines)
    if route_end is None:
        raise GeometryParseError("Gaussian input needs a route section starting with #")
    route_lines = _route_lines(lines)
    if _route_requests_zmatrix(route_lines):
        raise GeometryParseError("Gaussian Z-matrix input needs the unified Z-matrix adapter")

    idx = route_end + 1 if route_end < len(lines) else route_end
    title_start = idx
    while idx < len(lines) and lines[idx].strip():
        idx += 1
    title = " ".join(line.strip() for line in lines[title_start:idx] if line.strip())
    idx = _next_nonblank(lines, idx)
    if idx >= len(lines) or not _is_charge_multiplicity(lines[idx]):
        raise GeometryParseError("Gaussian input needs charge and multiplicity before coordinates")
    charge, multiplicity = (int(value) for value in lines[idx].split()[:2])
    idx += 1

    atoms: list[str] = []
    coords: list[list[float]] = []
    while idx < len(lines) and lines[idx].strip():
        atom, xyz = _parse_cartesian_line(lines[idx])
        atoms.append(atom)
        coords.append(xyz)
        idx += 1
    if not atoms:
        raise GeometryParseError("Gaussian input contains no Cartesian coordinate block")

    modredundant = lines[idx + 1 :] if idx < len(lines) else ()
    return MolecularGeometry(
        atoms=tuple(atoms),
        coordinates_angstrom=np.asarray(coords, dtype=float),
        comment=title or target.stem,
        source_format="gaussian_cartesian_input",
        source_path=target,
        charge=charge,
        multiplicity=multiplicity,
        fixed_parameters=_modredundant_fixed_patterns(modredundant),
        metadata={"route": tuple(route_lines)},
    )


def summarize_gaussian_log(path: Path) -> GaussianLogSummary:
    target = Path(path)
    text = target.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    scf = tuple(float(match.group(1)) for match in SCF_RE.finditer(text))
    orientation = _parse_last_orientation(lines, source_path=target)
    return GaussianLogSummary(
        path=target,
        normal_termination=NORMAL_TERMINATION in text,
        scf_energies_hartree=scf,
        standard_orientation_count=text.count(STANDARD_ORIENTATION),
        input_orientation_count=text.count(INPUT_ORIENTATION),
        scan_marker_count=text.lower().count("scan"),
        puckering_marker_count=text.count("QPck") + text.count("PhiP") + text.count("RPck"),
        frequencies_cm=tuple(_parse_frequencies(text)),
        last_orientation=orientation,
    )


def _route_end(lines: list[str]) -> int | None:
    idx = 0
    while idx < len(lines):
        if lines[idx].strip().startswith("#"):
            while idx < len(lines) and lines[idx].strip():
                idx += 1
            return idx
        idx += 1
    return None


def _route_lines(lines: list[str]) -> tuple[str, ...]:
    idx = 0
    while idx < len(lines):
        if lines[idx].strip().startswith("#"):
            route: list[str] = []
            while idx < len(lines) and lines[idx].strip():
                route.append(lines[idx].strip())
                idx += 1
            return tuple(route)
        idx += 1
    return ()


def _route_requests_zmatrix(route_lines: tuple[str, ...]) -> bool:
    route = " ".join(route_lines).lower()
    return "zmat" in route or "z-matrix" in route


def _next_nonblank(lines: list[str], idx: int) -> int:
    while idx < len(lines) and not lines[idx].strip():
        idx += 1
    return idx


def _is_charge_multiplicity(line: str) -> bool:
    parts = line.split()
    if len(parts) < 2:
        return False
    try:
        int(parts[0])
        int(parts[1])
    except ValueError:
        return False
    return True


def _parse_cartesian_line(line: str) -> tuple[str, list[float]]:
    parts = line.split()
    if len(parts) < 4:
        raise GeometryParseError(f"invalid Gaussian Cartesian line: {line}")
    try:
        coords = [float(parts[1]), float(parts[2]), float(parts[3])]
    except ValueError as exc:
        raise GeometryParseError(f"invalid Gaussian Cartesian coordinates: {line}") from exc
    return normalize_atom_symbol(parts[0]), coords


def _modredundant_fixed_patterns(lines: list[str] | tuple[str, ...]) -> tuple[str, ...]:
    patterns: list[str] = []
    seen: set[str] = set()
    for raw in lines:
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        parts = line.replace(",", " ").split()
        if not parts:
            continue
        kind = parts[0].upper()
        expected = {"B": 2, "A": 3, "D": 4, "O": 4, "L": 3}.get(kind)
        if expected is None:
            continue
        atoms: list[int] = []
        actions: list[str] = []
        for token in parts[1:]:
            try:
                atoms.append(int(token))
            except ValueError:
                actions.append(token.upper())
        if len(atoms) < expected or "F" not in actions:
            continue
        pattern = f"{kind}({','.join(str(atom) for atom in atoms[:expected])})"
        if pattern not in seen:
            patterns.append(pattern)
            seen.add(pattern)
    return tuple(patterns)


def _parse_frequencies(text: str) -> list[float]:
    values: list[float] = []
    for line in text.splitlines():
        if "Frequencies --" not in line:
            continue
        for token in line.split("--", 1)[1].split():
            try:
                values.append(float(token))
            except ValueError:
                continue
    return values


def _parse_last_orientation(lines: list[str], *, source_path: Path) -> MolecularGeometry | None:
    start = -1
    orientation_kind = ""
    for idx, line in enumerate(lines):
        if STANDARD_ORIENTATION in line or INPUT_ORIENTATION in line:
            start = idx
            orientation_kind = line.strip().rstrip(":")
    if start < 0:
        return None
    atoms: list[str] = []
    coords: list[list[float]] = []
    dash_count = 0
    for line in lines[start + 1 :]:
        if set(line.strip()) == {"-"}:
            dash_count += 1
            if dash_count >= 3:
                break
            continue
        if dash_count < 2:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            atom = normalize_atom_symbol(parts[1])
            xyz = [float(parts[3]), float(parts[4]), float(parts[5])]
        except ValueError:
            continue
        atoms.append(atom)
        coords.append(xyz)
    if not atoms:
        return None
    return MolecularGeometry(
        atoms=tuple(atoms),
        coordinates_angstrom=np.asarray(coords, dtype=float),
        comment=f"Gaussian {orientation_kind}",
        source_format="gaussian_log_orientation",
        source_path=source_path,
        metadata={"orientation": orientation_kind},
    )

