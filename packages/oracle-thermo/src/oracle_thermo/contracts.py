from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from oracle_core import (
    normalize_key,
    parse_key_value_section,
    read_sectioned_lines,
    replace_section,
    section_content,
)

from .models import THERMO_KEYS, THERMO_LABELS, ThermoContribution, ThermoSection


ORACLE_XYZ_THERMO_SCHEMA = "oracle.xyz.thermo.v1"

_LABEL_ALIASES = {
    "trasl": ("TRASL", "TRANS", "TRANSL", "TRANSLATIONAL"),
    "rot": ("ROT", "ROTATIONAL"),
    "vib": ("VIB", "VIBRATIONAL"),
    "tot": ("TOT", "TOTAL"),
}
_FIELD_KEYS = {normalize_key(key): key for key in THERMO_KEYS}


def parse_thermo_section(lines: Iterable[str]) -> ThermoSection:
    filtered = [line for line in lines if not line.lstrip().startswith("#")]
    values = parse_key_value_section(filtered)
    return ThermoSection(
        translational=_parse_contribution(values, "trasl"),
        rotational=_parse_contribution(values, "rot"),
        vibrational=_parse_contribution(values, "vib"),
        total=_parse_contribution(values, "tot"),
        schema=values.get("SCHEMA", ORACLE_XYZ_THERMO_SCHEMA),
    )


def thermo_section_lines(section: ThermoSection) -> list[str]:
    lines = [f"SCHEMA {ORACLE_XYZ_THERMO_SCHEMA}"]
    for label in THERMO_LABELS:
        contribution = section.contribution(label)
        if contribution is None:
            continue
        for key in THERMO_KEYS:
            value = getattr(contribution, key)
            if value is not None:
                lines.append(f"{key}_{label} = {_format_float(value)}")
        if not contribution.available:
            lines.append(f"available_{label} = 0")
        if contribution.reason and contribution.reason != "ok":
            lines.append(f"reason_{label} = {contribution.reason}")
    return lines


def read_thermo_section(path: Path) -> ThermoSection:
    return parse_thermo_section(section_content(read_sectioned_lines(Path(path)), "THERMO"))


def write_thermo_section(path: Path, section: ThermoSection) -> None:
    replace_section(Path(path), "THERMO", thermo_section_lines(section))


def _parse_contribution(
    values: dict[str, str],
    label: str,
) -> ThermoContribution | None:
    aliases = _LABEL_ALIASES[label]
    parsed: dict[str, float | None] = {}
    present = False
    for normalized_field, attr in _FIELD_KEYS.items():
        value = _first_value(values, *(f"{normalized_field}_{alias}" for alias in aliases))
        parsed[attr] = _optional_float(value)
        present = present or value is not None

    available_text = _first_value(values, *(f"AVAILABLE_{alias}" for alias in aliases))
    reason = _first_value(values, *(f"REASON_{alias}" for alias in aliases)) or "ok"
    if not present and available_text is None and reason == "ok":
        return None
    return ThermoContribution(
        Q_dimless=parsed["Q_dimless"],
        U_kJmol=parsed["U_kJmol"],
        H_kJmol=parsed["H_kJmol"],
        S_JmolK=parsed["S_JmolK"],
        Cv_JmolK=parsed["Cv_JmolK"],
        Cp_JmolK=parsed["Cp_JmolK"],
        available=_optional_bool(available_text, default=True),
        reason=reason,
    )


def _first_value(values: dict[str, str], *keys: str) -> str | None:
    for key in keys:
        value = values.get(key)
        if value is not None:
            return value
    return None


def _optional_float(text: str | None) -> float | None:
    if text is None:
        return None
    return float(text.replace("D", "E").replace("d", "e"))


def _optional_bool(text: str | None, *, default: bool) -> bool:
    if text is None:
        return default
    key = normalize_key(text)
    if key in {"1", "TRUE", "YES", "Y"}:
        return True
    if key in {"0", "FALSE", "NO", "N"}:
        return False
    return default


def _format_float(value: float) -> str:
    return f"{float(value):.12g}"
