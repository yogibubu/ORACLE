from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from matrix_core import (
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
from .topology.contracts import (
    MATRIX_XYZ_SYNTHONS_SCHEMA,
    MATRIX_XYZ_TOPOLOGY_SCHEMA,
    ORACLE_XYZ_SYNTHONS_SCHEMA,
    ORACLE_XYZ_TOPOLOGY_SCHEMA,
)
from .topology.pipeline import build_topology_objects


@dataclass(frozen=True)
class SymmetryThresholds:
    distance_angstrom: float = 1.0e-3
    inertia_relative: float = 1.0e-3
    max_rotation_order: int = 6


@dataclass(frozen=True)
class LinkPreprocessResult:
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
) -> LinkPreprocessResult:
    """Import a geometry source and materialize initial MATRIX sections.

    Program-specific and SMILES adapters call into LINK once they have produced
    a `MolecularGeometry`.
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
    return LinkPreprocessResult(
        path=Path(target),
        geometry=geometry,
        point_group=point_group,
        topology_bond_count=bond_count,
        ring_count=ring_count,
    )


BabelPreprocessResult = LinkPreprocessResult


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
    from matrix_gaussian import gaussian_topology_section_lines

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
    explicit_bond_orders = dict(geometry.metadata.get("explicit_bond_orders", {}))
    bond_orders = {**explicit_bond_orders, **gaussian["bond_orders"]}
    bond_order_source = (
        gaussian["bond_order_source"]
        if gaussian["bond_orders"]
        else geometry.metadata.get("bond_order_source", gaussian["bond_order_source"])
    )
    continuous, discrete, ringset, synthons, aromaticity = build_topology_objects(
        geometry.coordinates_angstrom,
        atomic_numbers,
        bond_order_overrides=bond_orders,
        external_charges=gaussian["charges"],
        charge_source=gaussian["charge_source"],
        bond_order_source=str(bond_order_source),
    )
    topology_lines = [
        f"SCHEMA {MATRIX_XYZ_TOPOLOGY_SCHEMA}",
        f"ALIAS_SCHEMA {ORACLE_XYZ_TOPOLOGY_SCHEMA}",
        "INDEXING ATOMS=ONE_BASED",
        f"BOND_ORDER_SOURCE {bond_order_source}",
        "RING_BASIS_POLICY CHORDLESS_NONMETAL_MINIMUM_CYCLE_BASIS",
        *_ring_basis_diagnostic_lines(ringset),
        "[BONDS]",
    ]
    if discrete.bonds:
        topology_lines.extend(f"{i + 1} {j + 1}" for i, j in discrete.bonds)
    else:
        topology_lines.append("NONE")
    topology_lines.append("[BOND_ORDERS]")
    bond_order_rows = []
    for i, j in discrete.bonds:
        try:
            value = float(synthons.bond_order(i, j))
        except Exception:
            continue
        bond_order_rows.append(f"{i + 1} {j + 1} {value:.10g}")
    if bond_order_rows:
        topology_lines.extend(bond_order_rows)
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
        f"SCHEMA {MATRIX_XYZ_SYNTHONS_SCHEMA}",
        f"ALIAS_SCHEMA {ORACLE_XYZ_SYNTHONS_SCHEMA}",
        "INDEXING ATOMS=ONE_BASED",
        f"CHARGE_SOURCE {gaussian['charge_source']}",
        f"BOND_ORDER_SOURCE {bond_order_source}",
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


def _ring_basis_diagnostic_lines(ringset) -> list[str]:
    diagnostics = getattr(ringset, "cycle_basis_diagnostics", None)
    if diagnostics is None:
        return []
    excluded = (
        ",".join(str(atom + 1) for atom in diagnostics.excluded_atoms)
        if diagnostics.excluded_atoms
        else "NONE"
    )
    return [
        f"RING_CANDIDATE_COUNT {diagnostics.candidate_cycle_count}",
        f"RING_BASIS_RANK {diagnostics.cycle_rank}",
        f"RING_BASIS_COUNT {diagnostics.selected_cycle_count}",
        f"RING_BASIS_ALLOWED_ATOMS {diagnostics.allowed_atom_count}",
        f"RING_BASIS_ALLOWED_EDGES {diagnostics.allowed_edge_count}",
        f"RING_BASIS_EXCLUDED_ATOMS {excluded}",
    ]


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
            f"Gaussian {bo_source}"
            if bond_orders and bo_source
            else "Topology Pauling continuous model"
        ),
    }


def _atomic_number(symbol: str) -> int:
    from .topology.elements import atomic_number

    number = atomic_number(symbol)
    if number is None or number <= 0:
        raise ValueError(f"unknown element symbol: {symbol}")
    return int(number)
