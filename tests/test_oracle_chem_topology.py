from __future__ import annotations

from oracle_chem import AtomicSynthons, build_topology_objects


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
