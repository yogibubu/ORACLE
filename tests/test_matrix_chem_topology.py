from __future__ import annotations

import json
from io import StringIO
from pathlib import Path

from matrix_chem import (
    AtomicSynthons,
    MATRIX_XYZ_TOPOLOGY_SCHEMA,
    TOPOLOGY_SNAPSHOT_SCHEMA,
    build_topology_objects,
    compare_topology_snapshot_entry,
    preprocess_to_enriched_xyz,
    topology_snapshot_from_xyzin,
)
from matrix_chem.topology.topology_writer import write_topology_section
from matrix_neo.definition import _topology_bonds, _topology_rings, topology_bond_orders_from_lines

CORPUS = Path(__file__).resolve().parent / "fixtures" / "test_molecules" / "molecules"
GOLDEN_TOPOLOGY = (
    Path(__file__).resolve().parent / "fixtures" / "golden_corpus" / "matrix_topology_snapshots.json"
)


def test_topology_builds_graph_and_atomic_synthons_for_water():
    coords = [
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 1.0, 0.0),
    ]
    atomic_numbers = [8, 1, 1]

    continuous, discrete, ringset, synthons, aromaticity = build_topology_objects(
        coords,
        atomic_numbers,
    )

    assert continuous.natoms == 3
    assert sorted(discrete.bonds) == [(0, 1), (0, 2)]
    assert ringset.rings == []
    assert isinstance(synthons, AtomicSynthons)
    assert synthons.Zeff(0) > 0.0
    assert synthons.canonical_signature(0)[0] == 8
    assert aromaticity.aromatic_atoms == set()


def test_topology_writer_uses_modern_xyzin_contract_for_neo_consumers():
    coords = [
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 1.0, 0.0),
    ]
    atomic_numbers = [8, 1, 1]
    continuous, discrete, ringset, synthons, aromaticity = build_topology_objects(
        coords,
        atomic_numbers,
    )

    stream = StringIO()
    write_topology_section(
        stream,
        cg=continuous,
        dg=discrete,
        ringset=ringset,
        synthons=synthons,
        aromaticity=aromaticity,
    )

    lines = [line.rstrip("\n") for line in stream.getvalue().splitlines()]
    topology_start = lines.index("#TOPOLOGY") + 1
    topology = lines[topology_start:]
    assert topology[0] == f"SCHEMA {MATRIX_XYZ_TOPOLOGY_SCHEMA}"
    assert "ALIAS_SCHEMA oracle.xyz.topology.v1" in topology
    assert "[BONDS]" in topology
    assert "[BOND_ORDERS]" in topology
    assert "[RINGS]" in topology
    assert _topology_bonds(lines, natoms=3) == ((1, 2), (1, 3))
    assert topology_bond_orders_from_lines(lines, natoms=3).keys() == {(1, 2), (1, 3)}
    assert _topology_rings(lines, natoms=3) == ()


def test_topology_golden_snapshots_are_stable(tmp_path):
    snapshot = json.loads(GOLDEN_TOPOLOGY.read_text(encoding="utf-8"))

    assert snapshot["schema"] == TOPOLOGY_SNAPSHOT_SCHEMA
    assert snapshot["rounding_decimals"] == 8
    assert len(snapshot["entries"]) >= 10

    for entry in snapshot["entries"]:
        xyzin = tmp_path / f"{entry['id']}.xyzin"
        preprocess_to_enriched_xyz(Path(entry["source"]), xyzin)
        comparison = compare_topology_snapshot_entry(
            entry,
            xyzin,
            rounding_decimals=snapshot["rounding_decimals"],
        )

        assert comparison.ok, "\n".join(comparison.messages)


def test_topology_uses_elementary_ring_basis_for_fused_cage_and_metallocene(tmp_path):
    expected = {
        "azulene.inp": 2,
        "pyrene.inp": 4,
        "cubane.inp": 5,
        "ferrocene.inp": 2,
        "ferrocene_staggered.inp": 2,
    }

    for molecule, ring_count in expected.items():
        xyzin = tmp_path / f"{Path(molecule).stem}.xyzin"
        preprocess_to_enriched_xyz(CORPUS / molecule, xyzin)
        snapshot = topology_snapshot_from_xyzin(xyzin)

        assert snapshot["ring_count"] == ring_count, molecule
        if molecule.startswith("ferrocene"):
            assert all(2 not in ring["atoms"] for ring in snapshot["rings"]), molecule
