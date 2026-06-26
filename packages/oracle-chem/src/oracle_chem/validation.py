from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from oracle_core import read_sectioned_lines, replace_section, section_content

from .geometry_io import GeometryParseError, read_enriched_xyz


ORACLE_XYZ_VALIDATION_SCHEMA = "oracle.xyz.validation.v1"

REQUIRED_SYMMETRY_SCHEMA = "oracle.xyz.symmetry.v1"
REQUIRED_TOPOLOGY_SCHEMA = "oracle.xyz.topology.v1"
REQUIRED_SYNTHONS_SCHEMA = "oracle.xyz.synthons.v1"
OPTIONAL_FRAGMENTS_SCHEMA = "oracle.xyz.fragments.v1"


@dataclass(frozen=True)
class ValidationMessage:
    level: str
    code: str
    text: str


@dataclass(frozen=True)
class ValidationResult:
    status: str
    messages: tuple[ValidationMessage, ...]


def validate_enriched_molecule(path: Path, *, require_fragments: bool = False) -> ValidationResult:
    """Validate an enriched XYZ after topology/synthon preprocessing."""
    target = Path(path)
    lines = read_sectioned_lines(target)
    messages: list[ValidationMessage] = []

    _validate_geometry(target, messages)
    _require_schema(lines, "SYMMETRY", REQUIRED_SYMMETRY_SCHEMA, messages)
    _require_schema(lines, "TOPOLOGY", REQUIRED_TOPOLOGY_SCHEMA, messages)
    _require_schema(lines, "SYNTHONS", REQUIRED_SYNTHONS_SCHEMA, messages)
    if require_fragments:
        _require_schema(lines, "FRAGMENTS", OPTIONAL_FRAGMENTS_SCHEMA, messages)
    elif section_content(lines, "FRAGMENTS"):
        _require_schema(lines, "FRAGMENTS", OPTIONAL_FRAGMENTS_SCHEMA, messages)

    status = _status_from_messages(messages)
    if not messages:
        messages.append(
            ValidationMessage("INFO", "VALIDATION_PASS", "Molecule is ready for GICForge")
        )
    return ValidationResult(status=status, messages=tuple(messages))


def validation_section_lines(result: ValidationResult) -> list[str]:
    lines = [
        f"SCHEMA {ORACLE_XYZ_VALIDATION_SCHEMA}",
        f"STATUS {result.status}",
        "DEPENDENCIES SYMMETRY=oracle.xyz.symmetry.v1 "
        "TOPOLOGY=oracle.xyz.topology.v1 SYNTHONS=oracle.xyz.synthons.v1",
        "OPTIONAL FRAGMENTS=oracle.xyz.fragments.v1",
        "[MESSAGES]",
    ]
    lines.extend(f"{msg.level} {msg.code} {msg.text}" for msg in result.messages)
    return lines


def write_validation_section(path: Path, *, require_fragments: bool = False) -> ValidationResult:
    result = validate_enriched_molecule(path, require_fragments=require_fragments)
    replace_section(Path(path), "VALIDATION", validation_section_lines(result))
    return result


def _validate_geometry(path: Path, messages: list[ValidationMessage]) -> None:
    try:
        geometry = read_enriched_xyz(path)
    except GeometryParseError as exc:
        messages.append(ValidationMessage("ERROR", "INVALID_XYZ", str(exc)))
        return
    if geometry.natoms <= 0:
        messages.append(ValidationMessage("ERROR", "NO_ATOMS", "XYZ block contains no atoms"))
        return
    for idx, row in enumerate(geometry.coordinates_angstrom, start=1):
        if not all(math.isfinite(float(value)) for value in row):
            messages.append(
                ValidationMessage(
                    "ERROR",
                    "NONFINITE_COORDINATE",
                    f"Atom {idx} has invalid coordinates",
                )
            )
            return


def _require_schema(
    lines: list[str],
    section_name: str,
    schema: str,
    messages: list[ValidationMessage],
) -> None:
    content = section_content(lines, section_name)
    if not content:
        messages.append(
            ValidationMessage("ERROR", f"MISSING_{section_name}", f"Missing #{section_name}")
        )
        return
    expected = f"SCHEMA {schema}"
    if content[0].strip() != expected:
        messages.append(
            ValidationMessage(
                "ERROR",
                f"BAD_{section_name}_SCHEMA",
                f"#{section_name} must start with {expected}",
            )
        )


def _status_from_messages(messages: list[ValidationMessage]) -> str:
    levels = {message.level for message in messages}
    if "ERROR" in levels:
        return "FAIL"
    if "WARN" in levels:
        return "WARN"
    return "PASS"
