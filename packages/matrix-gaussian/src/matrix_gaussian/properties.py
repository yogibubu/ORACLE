from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(frozen=True)
class GaussianQuadrupolePromotion:
    xyzin: Path
    log_path: Path
    wrote_properties: bool
    property_count: int


def parse_gaussian_quadrupole_properties(path: Path | str) -> tuple:
    from matrix_qm import atomic_number_from_isotope_or_atom, quadrupole_property_records_from_nqcc

    target = Path(path)
    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
    records = []
    for row in _quadrupole_rows(lines):
        number = atomic_number_from_isotope_or_atom(
            isotope=row["isotope"],
            atom_symbol=row["atom_symbol"],
        )
        records.extend(
            quadrupole_property_records_from_nqcc(
                atom=row["atom"],
                atomic_number_value=number,
                nqcc_mhz=row["values"],
                program="Gaussian",
                source=target,
                isotope=row["isotope"],
                method=_parse_route_method(lines),
                level="",
                axes=row["axes"],
                status="raw",
                comment=row["comment"],
            )
        )
    return tuple(records)


def promote_gaussian_quadrupole_properties_to_xyzin(
    log_path: Path | str,
    xyzin: Path | str,
) -> GaussianQuadrupolePromotion:
    from matrix_qm import merge_properties_section

    source = Path(log_path)
    target = Path(xyzin)
    records = parse_gaussian_quadrupole_properties(source)
    if records:
        merge_properties_section(target, records)
    return GaussianQuadrupolePromotion(
        xyzin=target,
        log_path=source,
        wrote_properties=bool(records),
        property_count=len(records),
    )


def _quadrupole_rows(lines: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    in_section = False
    for raw in lines:
        text = raw.strip()
        upper = text.upper()
        if "QUADRUPOLE" in upper and ("COUPLING" in upper or "NQCC" in upper):
            in_section = True
            continue
        if in_section and not text:
            continue
        if in_section and _section_ended(upper):
            in_section = False
            continue
        if not in_section and not upper.startswith(("NQCC", "CHI")):
            continue
        row = _parse_direct_nqcc_row(text)
        if row is not None:
            rows.append(row)
    return rows


def _parse_direct_nqcc_row(text: str) -> dict[str, object] | None:
    clean = (
        text.replace("=", " ")
        .replace(",", " ")
        .replace(":", " ")
        .replace("(", " ")
        .replace(")", " ")
    )
    parts = clean.split()
    if not parts:
        return None
    if parts[0].upper() in {"NQCC", "CHI", "ATOM"}:
        parts = parts[1:]
    if len(parts) < 5 or not parts[0].isdigit():
        return None
    atom = int(parts[0])
    isotope = ""
    atom_symbol = ""
    idx = 1
    if idx < len(parts) and _looks_isotope(parts[idx]):
        isotope = _normalize_isotope(parts[idx])
        idx += 1
    elif idx < len(parts) and re.fullmatch(r"[A-Za-z]{1,2}", parts[idx]):
        atom_symbol = parts[idx].capitalize()
        idx += 1
        if idx < len(parts) and parts[idx].isdigit():
            isotope = f"{parts[idx]}{atom_symbol}"
            idx += 1
    values = _first_floats(parts[idx:], count=6)
    if len(values) < 3:
        return None
    axes = "GAUSSIAN:chi_xx,chi_yy,chi_zz"
    if len(values) >= 6:
        axes += ",chi_xy,chi_xz,chi_yz"
    return {
        "atom": atom,
        "isotope": isotope,
        "atom_symbol": atom_symbol,
        "values": tuple(values[:6] if len(values) >= 6 else values[:3]),
        "axes": axes,
        "comment": "Gaussian nuclear quadrupole coupling constants in MHz",
    }


def _first_floats(parts: list[str], *, count: int) -> list[float]:
    values: list[float] = []
    for part in parts:
        try:
            values.append(float(part.replace("D", "E").replace("d", "e")))
        except ValueError:
            continue
        if len(values) == count:
            break
    return values


def _looks_isotope(text: str) -> bool:
    return re.fullmatch(r"(?:\d+[A-Za-z]{1,2}|[A-Za-z]{1,2}\d+)", text) is not None


def _normalize_isotope(text: str) -> str:
    match = re.fullmatch(r"(\d+)([A-Za-z]{1,2})", text)
    if match is not None:
        return f"{int(match.group(1))}{match.group(2).capitalize()}"
    match = re.fullmatch(r"([A-Za-z]{1,2})(\d+)", text)
    if match is not None:
        return f"{int(match.group(2))}{match.group(1).capitalize()}"
    return text


def _section_ended(upper: str) -> bool:
    return upper.startswith(("----", "====", "DIPOLE", "THERMOCHEMISTRY", "NORMAL TERMINATION"))


def _parse_route_method(lines: list[str]) -> str:
    for raw in lines:
        text = raw.strip()
        if not text.startswith("#"):
            continue
        route = text.lstrip("#").strip()
        for token in route.split():
            if "/" in token:
                return token
    return ""
