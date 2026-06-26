from __future__ import annotations

from pathlib import Path

from oracle_chem import read_enriched_xyz
from oracle_core import read_sectioned_lines, section_content


ORACLE_GAUSSIAN_GIC_INPUT_SCHEMA = "oracle.gaussian.gic_input.v1"
REQUIRED_GIC_SCHEMA = "oracle.xyz.gic.v1"


class GaussianWriteError(ValueError):
    """Raised when ORACLE state cannot be exported to Gaussian input."""


def write_gicforge_gaussian_input(
    enriched_xyz: Path,
    output: Path,
    *,
    route: str = "#p hf/sto-3g",
    title: str | None = None,
    charge: int | None = None,
    multiplicity: int | None = None,
    link0: tuple[str, ...] = (),
) -> Path:
    """Write a Gaussian input file from an enriched XYZ carrying #GIC."""
    source = Path(enriched_xyz)
    _require_gic_section(source)
    geometry = read_enriched_xyz(source)
    job_charge = (
        charge
        if charge is not None
        else geometry.charge
        if geometry.charge is not None
        else 0
    )
    job_multiplicity = (
        multiplicity
        if multiplicity is not None
        else geometry.multiplicity
        if geometry.multiplicity is not None
        else 1
    )
    lines = [
        *[item.strip() for item in link0 if item.strip()],
        _normalize_route(route),
        "",
        title or geometry.comment or source.stem,
        "",
        f"{job_charge} {job_multiplicity}",
    ]
    for atom, (x, y, z) in zip(geometry.atoms, geometry.coordinates_angstrom):
        lines.append(f"{atom:2s} {x:15.8f} {y:15.8f} {z:15.8f}")
    lines.extend(
        [
            "",
            f"! ORACLE_SCHEMA {ORACLE_GAUSSIAN_GIC_INPUT_SCHEMA}",
            "! ORACLE_GIC_SOURCE enriched_xyz",
            "",
        ]
    )
    gic_lines = _gaussian_gic_lines(source)
    if gic_lines:
        lines.extend(gic_lines)
        lines.append("")
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines), encoding="utf-8")
    return target


def _require_gic_section(path: Path) -> None:
    lines = read_sectioned_lines(Path(path))
    gic = section_content(lines, "GIC")
    if not gic:
        raise GaussianWriteError("missing #GIC section")
    expected = f"SCHEMA {REQUIRED_GIC_SCHEMA}"
    if gic[0].strip() != expected:
        raise GaussianWriteError(f"#GIC must start with {expected!r}; found {gic[0]!r}")


def _gaussian_gic_lines(path: Path) -> list[str]:
    try:
        from oracle_gicforge import gaussian_gic_lines_from_xyzin
    except ImportError:
        return []
    return gaussian_gic_lines_from_xyzin(Path(path))


def _normalize_route(route: str) -> str:
    text = route.strip()
    if not text:
        raise GaussianWriteError("Gaussian route cannot be empty")
    return text if text.startswith("#") else f"# {text}"
