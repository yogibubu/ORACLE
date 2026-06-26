from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from oracle_gf.models import HessianInput

from .models import AnharmonicInput
from .vci import QuarticForceField


_HEADER_RE = re.compile(r"^(?P<label>.+?)\s+(?P<kind>[IRC])(?:\s+N=\s*(?P<count>\d+)|\s+(?P<scalar>[-+]?\d+))\s*$")


@dataclass(frozen=True)
class FCHKData:
    """Gaussian adapter output kept separate from canonical ORACLE models."""

    atomic_numbers: np.ndarray
    cartesian_coordinates_bohr: np.ndarray
    masses_amu: np.ndarray
    cartesian_hessian_lower: np.ndarray
    harmonic_frequencies_cm: np.ndarray
    anharmonic_frequencies_cm: np.ndarray
    anharmonic_e2: np.ndarray
    normal_modes: np.ndarray

    def to_hessian_input(self) -> HessianInput:
        hessian = lower_to_symmetric(self.cartesian_hessian_lower)
        data = HessianInput(
            atomic_numbers=self.atomic_numbers,
            cartesian_coordinates_bohr=self.cartesian_coordinates_bohr,
            masses_amu=self.masses_amu,
            cartesian_hessian=hessian,
            harmonic_frequencies_cm=self.harmonic_frequencies_cm,
            source="gaussian-fchk",
        )
        data.validate()
        return data

    def to_anharmonic_input(self) -> AnharmonicInput:
        data = AnharmonicInput(
            harmonic_frequencies_cm=self.harmonic_frequencies_cm,
            anharmonic_frequencies_cm=self.anharmonic_frequencies_cm,
            cubic_cm={},
            quartic_cm={},
            source="gaussian-fchk",
        )
        data.validate()
        return data


def _read_fchk_blocks(path: Path) -> dict[str, np.ndarray | int | float | str]:
    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    blocks: dict[str, np.ndarray | int | float | str] = {}
    i = 0
    while i < len(lines):
        match = _HEADER_RE.match(lines[i])
        if match is None:
            i += 1
            continue
        label = match.group("label").strip()
        kind = match.group("kind")
        count = match.group("count")
        scalar = match.group("scalar")
        if count is None:
            if scalar is not None:
                blocks[label] = int(scalar) if kind == "I" else float(scalar)
            i += 1
            continue

        nvalues = int(count)
        raw: list[str] = []
        i += 1
        while len(raw) < nvalues and i < len(lines):
            raw.extend(lines[i].split())
            i += 1
        if kind == "I":
            blocks[label] = np.array([int(x) for x in raw[:nvalues]], dtype=int)
        elif kind == "R":
            blocks[label] = np.array([float(x.replace("D", "E")) for x in raw[:nvalues]], dtype=float)
        else:
            blocks[label] = " ".join(raw[:nvalues])
    return blocks


def _first_array(blocks: dict[str, object], *labels: str) -> np.ndarray:
    for label in labels:
        value = blocks.get(label)
        if isinstance(value, np.ndarray):
            return value.astype(float, copy=False)
    raise ValueError(f"FCHK block not found: {' / '.join(labels)}")


def read_gaussian_fchk_qff(path: Path) -> FCHKData:
    """Read harmonic Hessian and Gaussian anharmonic arrays from an FCHK file."""
    blocks = _read_fchk_blocks(Path(path))
    atomic_numbers = _first_array(blocks, "Atomic numbers").astype(int)
    coords = _first_array(blocks, "Current cartesian coordinates").reshape((-1, 3))
    masses = _first_array(blocks, "Real atomic weights", "Atomic masses", "Vib-AtMass", "Anharmonic Vib-AtMass")
    hessian = _first_array(blocks, "Cartesian Force Constants")
    vib_e2 = _first_array(blocks, "Vib-E2") if "Vib-E2" in blocks else np.array((), dtype=float)
    anh_e2 = _first_array(blocks, "Anharmonic Vib-E2") if "Anharmonic Vib-E2" in blocks else np.array((), dtype=float)
    modes = _first_array(blocks, "Anharmonic Vib-Modes", "Vib-Modes") if (
        "Anharmonic Vib-Modes" in blocks or "Vib-Modes" in blocks
    ) else np.array((), dtype=float)

    harmonic = vib_e2[: int(blocks.get("Vib-NDim", len(masses) * 3))] if vib_e2.size else np.array((), dtype=float)
    n_anh = int(blocks.get("Anharmonic Vib-NDim", harmonic.size))
    anharmonic = anh_e2[:n_anh] if anh_e2.size else np.array((), dtype=float)
    return FCHKData(
        atomic_numbers=atomic_numbers,
        cartesian_coordinates_bohr=coords,
        masses_amu=masses,
        cartesian_hessian_lower=hessian,
        harmonic_frequencies_cm=harmonic,
        anharmonic_frequencies_cm=anharmonic,
        anharmonic_e2=anh_e2,
        normal_modes=modes,
    )


def hessian_input_from_gaussian_fchk(path: Path) -> HessianInput:
    """Gaussian FCHK adapter: return the canonical ORACLE Hessian input."""
    return read_gaussian_fchk_qff(path).to_hessian_input()


def anharmonic_input_from_gaussian_fchk(path: Path) -> AnharmonicInput:
    """Gaussian FCHK adapter: return the canonical ORACLE anharmonic input."""
    return read_gaussian_fchk_qff(path).to_anharmonic_input()


def lower_to_symmetric(lower: np.ndarray) -> np.ndarray:
    """Convert Gaussian lower-triangular packed storage to a full matrix."""
    n_float = (np.sqrt(8 * len(lower) + 1) - 1) / 2
    n = int(round(n_float))
    if n * (n + 1) // 2 != len(lower):
        raise ValueError("Packed lower-triangular array has an invalid length")
    mat = np.zeros((n, n), dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i + 1):
            mat[i, j] = lower[idx]
            mat[j, i] = lower[idx]
            idx += 1
    return mat


def read_indexed_qff_text(path: Path, frequencies_cm: np.ndarray | None = None) -> QuarticForceField:
    """Read an indexed cubic/quartic normal-coordinate force field.

    Accepted records are intentionally simple and solver-independent:

    - `FREQ i value_cm`
    - `CUBIC i j k value_cm`
    - `QUARTIC i j k l value_cm`

    Mode indices are one-based in the file and converted to zero-based tuples.
    Header lines from Gaussian-like sections are ignored, so extracted blocks can
    be pasted directly into a normalized `.qff` file.
    """
    freqs: dict[int, float] = {}
    cubic: dict[tuple[int, int, int], float] = {}
    quartic: dict[tuple[int, int, int, int], float] = {}
    for raw in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("!"):
            continue
        parts = line.replace(",", " ").split()
        tag = parts[0].upper()
        try:
            if tag in {"FREQ", "FREQUENCY"} and len(parts) >= 3:
                freqs[int(parts[1]) - 1] = float(parts[2].replace("D", "E"))
            elif tag in {"CUBIC", "C3"} and len(parts) >= 5:
                key = tuple(sorted(int(x) - 1 for x in parts[1:4]))
                cubic[key] = float(parts[4].replace("D", "E"))
            elif tag in {"QUARTIC", "C4"} and len(parts) >= 6:
                key = tuple(sorted(int(x) - 1 for x in parts[1:5]))
                quartic[key] = float(parts[5].replace("D", "E"))
            elif len(parts) in {4, 5} and all(p.lstrip("+-").isdigit() for p in parts[:-1]):
                key = tuple(sorted(int(x) - 1 for x in parts[:-1]))
                value = float(parts[-1].replace("D", "E"))
                if len(key) == 3:
                    cubic[key] = value
                else:
                    quartic[key] = value
        except ValueError:
            continue

    if frequencies_cm is None:
        if not freqs:
            raise ValueError("No frequencies found in indexed QFF text")
        n_modes = max(freqs) + 1
        frequencies = np.array([freqs[i] for i in range(n_modes)], dtype=float)
    else:
        frequencies = np.asarray(frequencies_cm, dtype=float)
        if freqs:
            frequencies = frequencies.copy()
            for index, value in freqs.items():
                if index >= len(frequencies):
                    raise ValueError("QFF frequency index exceeds supplied frequency array")
                frequencies[index] = value
    return QuarticForceField(frequencies, cubic, quartic)
