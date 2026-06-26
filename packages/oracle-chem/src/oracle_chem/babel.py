from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from oracle_core import (
    BasicSection,
    read_sectioned_lines,
    replace_section,
    replace_xyz_block,
    section_content,
    write_basic_section,
)

from .geometry import MolecularGeometry
from .geometry_io import GeometrySourceKind, read_geometry_with_kind
from .symmetry import analyze_molecular_symmetry, symmetry_section_lines
from .topology.pipeline import build_topology_objects


@dataclass(frozen=True)
class SymmetryThresholds:
    distance_angstrom: float = 1.0e-3
    inertia_relative: float = 1.0e-3
    max_rotation_order: int = 6


@dataclass(frozen=True)
class BabelPreprocessResult:
    path: Path
    geometry: MolecularGeometry
    point_group: str
    topology_bond_count: int
    ring_count: int


def preprocess_to_enriched_xyz(
    source: Path,
    target: Path,
    *,
    source_kind: GeometrySourceKind = "auto",
    symmetry_thresholds: SymmetryThresholds = SymmetryThresholds(),
) -> BabelPreprocessResult:
    """Import a geometry source and materialize initial ORACLE sections.

    This is the initial ORACLE-Babel scaffold. Program-specific and SMILES
    adapters should call into this once they have produced `MolecularGeometry`.
    """
    geometry = read_geometry_with_kind(Path(source), source_kind)
    write_enriched_geometry(target, geometry)
    write_source_section(target, source=Path(source), source_kind=source_kind, geometry=geometry)
    write_gaussian_topology_section(target, source=Path(source))
    symmetry = determine_initial_symmetry(geometry, symmetry_thresholds)
    point_group = symmetry.point_group
    write_basic_section_from_geometry(target, geometry=geometry, point_group=point_group)
    write_symmetry_section(target, symmetry=symmetry, thresholds=symmetry_thresholds)
    bond_count, ring_count = write_topology_and_synthons_sections(target, geometry)
    return BabelPreprocessResult(
        path=Path(target),
        geometry=geometry,
        point_group=point_group,
        topology_bond_count=bond_count,
        ring_count=ring_count,
    )


def write_enriched_geometry(path: Path, geometry: MolecularGeometry) -> None:
    replace_xyz_block(Path(path), geometry.xyz_lines())


def write_source_section(
    path: Path,
    *,
    source: Path,
    source_kind: str,
    geometry: MolecularGeometry,
) -> None:
    replace_section(
        Path(path),
        "SOURCE",
        [
            "SCHEMA oracle.xyz.source.v1",
            f"KIND {source_kind}",
            f"FORMAT {geometry.source_format}",
            f"PATH {source}",
        ],
    )


def write_basic_section_from_geometry(
    path: Path,
    *,
    geometry: MolecularGeometry,
    point_group: str,
) -> None:
    write_basic_section(
        Path(path),
        BasicSection(
            charge=0 if geometry.charge is None else int(geometry.charge),
            multiplicity=1 if geometry.multiplicity is None else int(geometry.multiplicity),
            point_group=point_group,
        ),
    )


def write_gaussian_topology_section(path: Path, *, source: Path) -> int:
    """Write Gaussian CM5/Mayer annotations when the import source provides them."""
    suffix = Path(source).suffix.lower()
    if suffix not in {".log", ".out"}:
        return 0
    from oracle_gaussian import gaussian_topology_section_lines

    lines = gaussian_topology_section_lines(Path(source))
    if not lines:
        return 0
    replace_section(Path(path), "GAUSSIAN_TOPOLOGY", lines)
    return len(lines)


def determine_initial_symmetry(
    geometry: MolecularGeometry,
    thresholds: SymmetryThresholds,
) -> object:
    return analyze_molecular_symmetry(
        geometry,
        distance_tolerance=thresholds.distance_angstrom,
        inertia_tolerance=thresholds.inertia_relative,
        max_rotation_order=thresholds.max_rotation_order,
    )


def determine_initial_point_group(
    geometry: MolecularGeometry,
    thresholds: SymmetryThresholds,
) -> str:
    return str(determine_initial_symmetry(geometry, thresholds).point_group)


def write_symmetry_section(
    path: Path,
    *,
    symmetry,
    thresholds: SymmetryThresholds,
) -> None:
    replace_section(
        Path(path),
        "SYMMETRY",
        symmetry_section_lines(symmetry, thresholds=thresholds),
    )


def write_topology_and_synthons_sections(
    path: Path,
    geometry: MolecularGeometry,
) -> tuple[int, int]:
    atomic_numbers = [_atomic_number(atom) for atom in geometry.atoms]
    gaussian = gaussian_topology_overrides_from_xyzin(Path(path))
    continuous, discrete, ringset, synthons, aromaticity = build_topology_objects(
        geometry.coordinates_angstrom,
        atomic_numbers,
        bond_order_overrides=gaussian["bond_orders"],
        external_charges=gaussian["charges"],
        charge_source=gaussian["charge_source"],
        bond_order_source=gaussian["bond_order_source"],
    )
    topology_lines = [
        "SCHEMA oracle.xyz.topology.v1",
        "INDEXING ATOMS=ONE_BASED",
        f"BOND_ORDER_SOURCE {gaussian['bond_order_source']}",
        "[BONDS]",
    ]
    if discrete.bonds:
        topology_lines.extend(f"{i + 1} {j + 1}" for i, j in discrete.bonds)
    else:
        topology_lines.append("NONE")
    topology_lines.append("[RINGS]")
    if ringset.rings:
        for idx, ring in enumerate(ringset.rings, start=1):
            atoms = " ".join(str(atom + 1) for atom in ring.atoms)
            topology_lines.append(f"{idx} SIZE={len(ring)} ATOMS={atoms}")
    else:
        topology_lines.append("NONE")
    topology_lines.append("[AROMATICITY]")
    aromatic_atoms = sorted(getattr(aromaticity, "aromatic_atoms", set()))
    topology_lines.append(
        "ATOMS "
        + (" ".join(str(atom + 1) for atom in aromatic_atoms) if aromatic_atoms else "NONE")
    )
    replace_section(Path(path), "TOPOLOGY", topology_lines)

    synthon_lines = [
        "SCHEMA oracle.xyz.synthons.v1",
        "INDEXING ATOMS=ONE_BASED",
        f"CHARGE_SOURCE {gaussian['charge_source']}",
        f"BOND_ORDER_SOURCE {gaussian['bond_order_source']}",
        "COLUMNS ATOM Z ZEFF CHARGE COVALENCY DELOCALIZATION STRAIN SIGNATURE",
    ]
    for idx, atom in enumerate(geometry.atoms):
        signature = synthons.canonical_signature(idx)
        signature_text = ",".join(str(item) for item in signature)
        synthon_lines.append(
            f"{idx + 1} {atom} "
            f"{float(synthons.Zeff(idx)):.8g} "
            f"{float(synthons.charge(idx)):.8g} "
            f"{float(synthons.covalency(idx)):.8g} "
            f"{float(synthons.delocalization(idx)):.8g} "
            f"{float(synthons.strain(idx)):.8g} "
            f"{signature_text}"
        )
    replace_section(Path(path), "SYNTHONS", synthon_lines)
    return len(discrete.bonds), len(ringset.rings)


def gaussian_topology_overrides_from_xyzin(path: Path) -> dict[str, object]:
    """Read #GAUSSIAN_TOPOLOGY as ORACLE topology overrides."""
    content = section_content(read_sectioned_lines(Path(path)), "GAUSSIAN_TOPOLOGY")
    charges: dict[int, float] = {}
    bond_orders: dict[tuple[int, int], float] = {}
    bo_source: str | None = None
    for raw in content:
        text = raw.strip()
        if not text or text.upper().startswith(("SCHEMA ", "INDEXING ", "CM5_COUNT", "BO_COUNT")):
            continue
        parts = text.replace("=", " = ").split()
        key = parts[0].upper() if parts else ""
        if key == "CM5" and len(parts) >= 3:
            idx = int(parts[1]) - 1
            if idx >= 0:
                charges[idx] = float(parts[2])
            continue
        if key == "BO_SOURCE" and len(parts) >= 2:
            bo_source = parts[2] if len(parts) >= 3 and parts[1] == "=" else parts[1]
            continue
        if key == "BO" and len(parts) >= 4:
            i = int(parts[1]) - 1
            j = int(parts[2]) - 1
            if i >= 0 and j >= 0 and i != j:
                pair = (i, j) if i < j else (j, i)
                bond_orders[pair] = float(parts[3])
    return {
        "charges": charges,
        "bond_orders": bond_orders,
        "charge_source": "Gaussian CM5" if charges else "Synthons electronegativity model",
        "bond_order_source": (
            f"Gaussian {bo_source}" if bond_orders and bo_source else "Topology Pauling continuous model"
        ),
    }


def _atomic_number(symbol: str) -> int:
    from .topology.elements import atomic_number

    number = atomic_number(symbol)
    if number is None or number <= 0:
        raise ValueError(f"unknown element symbol: {symbol}")
    return int(number)
