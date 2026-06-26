from __future__ import annotations

from pathlib import Path

from oracle_core import read_sectioned_lines, replace_section, section_content


ORACLE_XYZ_GIC_SCHEMA = "oracle.xyz.gic.v1"
ORACLE_XYZ_SYCART_SCHEMA = "oracle.xyz.sycart.v1"
REQUIRED_VALIDATION_SCHEMA = "oracle.xyz.validation.v1"


class GICForgeContractError(ValueError):
    """Raised when GICForge cannot consume the enriched XYZ state."""


def validate_gicforge_prerequisites(path: Path) -> None:
    lines = read_sectioned_lines(Path(path))
    validation = section_content(lines, "VALIDATION")
    if not validation:
        raise GICForgeContractError("missing #VALIDATION section")
    expected = f"SCHEMA {REQUIRED_VALIDATION_SCHEMA}"
    if validation[0].strip() != expected:
        raise GICForgeContractError(
            f"#VALIDATION must start with {expected!r}; found {validation[0]!r}"
        )
    status = _validation_status(validation)
    if status != "PASS":
        raise GICForgeContractError(f"#VALIDATION status must be PASS; found {status or 'UNKNOWN'}")


def gic_plan_section_lines(*, symmetrize: bool = False) -> list[str]:
    return [
        f"SCHEMA {ORACLE_XYZ_GIC_SCHEMA}",
        "STATUS PLANNED",
        "DEPENDENCIES VALIDATION=oracle.xyz.validation.v1 "
        "TOPOLOGY=oracle.xyz.topology.v1 SYNTHONS=oracle.xyz.synthons.v1 "
        "SYMMETRY=oracle.xyz.symmetry.v1",
        "INDEXING ATOMS=ONE_BASED",
        f"SYMMETRIZE {_bool_text(symmetrize)}",
        "BACKEND UNASSIGNED",
        "[FROZEN_GICS]",
        "PENDING GICFORGE_IMPLEMENTATION",
    ]


def sycart_plan_section_lines() -> list[str]:
    return [
        f"SCHEMA {ORACLE_XYZ_SYCART_SCHEMA}",
        "STATUS PLANNED",
        "DEPENDENCIES VALIDATION=oracle.xyz.validation.v1 GIC=oracle.xyz.gic.v1",
        "INDEXING ATOMS=ONE_BASED",
        "[SYCART]",
        "PENDING GICFORGE_IMPLEMENTATION",
    ]


def write_gicforge_plan_sections(
    path: Path,
    *,
    symmetrize: bool = False,
    sycart: bool = False,
) -> None:
    target = Path(path)
    validate_gicforge_prerequisites(target)
    replace_section(target, "GIC", gic_plan_section_lines(symmetrize=symmetrize))
    if sycart:
        replace_section(target, "SYCART", sycart_plan_section_lines())


def _validation_status(validation_lines: list[str]) -> str | None:
    for line in validation_lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0].upper() == "STATUS":
            return parts[1].upper()
    return None


def _bool_text(value: bool) -> str:
    return "TRUE" if value else "FALSE"
