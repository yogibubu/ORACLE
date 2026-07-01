from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from matrix_core import section_content
from matrix_chem import preprocess_to_enriched_xyz
from matrix_fragments import (
    FragmentContractError,
    read_interaction_center_definition,
    read_fragment_records,
    write_fragment_build_section,
    write_interaction_center_section,
    write_fragment_plan_section,
)


def _test_molecule_path(name: str) -> Path:
    return Path(__file__).resolve().parent / "fixtures" / "test_molecules" / "molecules" / name


def test_fragment_plan_requires_topology_and_synthons(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(FragmentContractError, match="missing #SYNTHONS"):
        write_fragment_plan_section(path)


def test_fragment_plan_section_preserves_existing_sections(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "2",
                "h2",
                "H 0 0 0",
                "H 0 0 1",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
                "INDEXING ATOMS=ONE_BASED",
                "[BONDS]",
                "1 2",
                "",
                "#SYNTHONS",
                "SCHEMA oracle.xyz.synthons.v1",
                "INDEXING ATOMS=ONE_BASED",
                "",
                "#GIC",
                "SCHEMA oracle.xyz.gic.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    write_fragment_plan_section(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    assert section_content(lines, "GIC")[0] == "SCHEMA oracle.xyz.gic.v1"
    fragments = section_content(lines, "FRAGMENTS")
    assert fragments[0] == "SCHEMA matrix.xyz.fragments.v1"
    assert "ALIAS_SCHEMA oracle.xyz.fragments.v1" in fragments
    assert (
        "DEPENDENCIES TOPOLOGY=matrix.xyz.topology.v1 SYNTHONS=matrix.xyz.synthons.v1"
        in fragments
    )
    assert "PENDING ROBUST_TOPOLOGY_CONTRACT" in fragments


def test_fragment_build_materializes_components_centers_and_frames(tmp_path):
    source = tmp_path / "dimer.xyz"
    source.write_text(
        "\n".join(
            [
                "6",
                "water dimer fragments",
                "O 0.00 0.00 0.00",
                "H 0.96 0.00 0.00",
                "H -0.24 0.93 0.00",
                "O 0.00 0.00 3.20",
                "H 0.96 0.00 3.20",
                "H -0.24 0.93 3.20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "dimer.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    definition = write_fragment_build_section(xyzin)
    lines = xyzin.read_text(encoding="utf-8").splitlines()
    fragments = section_content(lines, "FRAGMENTS")
    records = read_fragment_records(xyzin)

    assert definition.reference_fragment == "F001"
    assert len(definition.fragments) == 2
    assert "STATUS BUILT" in fragments
    assert "FRAGMENT_COUNT 2" in fragments
    assert "F001 LABEL=component_1 SIZE=3 ATOMS=1,2,3" in fragments
    assert "F002 LABEL=component_2 SIZE=3 ATOMS=4,5,6" in fragments
    assert any(line.startswith("F001 X=") for line in fragments)
    assert any(line.startswith("F001 X=") and " Y=" in line for line in fragments)
    assert len(records) == 2
    assert records[0].atoms == (1, 2, 3)


def test_fragment_workflow_handles_non_covalent_formic_acid_water_probe(tmp_path):
    source = tmp_path / "formic_acid_water.xyz"
    source.write_text(
        "\n".join(
            [
                "8",
                "formic acid-water non-covalent probe",
                "6    -1.171727   -0.018999   -0.001370",
                "1    -2.256869    0.130369   -0.015107",
                "8    -0.651376   -1.113045    0.004964",
                "8    -0.533544    1.143228    0.007603",
                "1     0.435203    0.960633    0.019343",
                "8     1.930145   -0.025976   -0.052269",
                "1     2.594643   -0.128917    0.633345",
                "1     1.351317   -0.802634    0.008831",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "formic_acid_water.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    fragments = write_fragment_build_section(xyzin)
    centers = write_interaction_center_section(xyzin)
    read_back = read_interaction_center_definition(xyzin)

    assert len(fragments.fragments) == 2
    assert fragments.fragments[0].atoms == (1, 2, 3, 4, 5)
    assert fragments.fragments[1].atoms == (6, 7, 8)
    assert len(centers.centers) == 6
    assert centers.interactions == ()
    assert [center.identifier for center in read_back.centers] == [
        center.identifier for center in centers.centers
    ]
    assert [center.kind for center in read_back.centers] == [
        center.kind for center in centers.centers
    ]
    assert [center.atoms for center in read_back.centers] == [
        center.atoms for center in centers.centers
    ]
    assert np.allclose(
        [center.center for center in read_back.centers],
        [center.center for center in centers.centers],
    )
    assert read_back.interactions == ()


def test_interaction_centers_materialize_ring_center_atom_candidate(tmp_path):
    xyzin = tmp_path / "ring_metal.xyzin"
    xyzin.write_text(
        "\n".join(
            [
                "7",
                "ring center interaction",
                "C 1.000000 0.000000 0.000000",
                "C 0.500000 0.866025 0.000000",
                "C -0.500000 0.866025 0.000000",
                "C -1.000000 0.000000 0.000000",
                "C -0.500000 -0.866025 0.000000",
                "C 0.500000 -0.866025 0.000000",
                "Na 0.000000 0.000000 2.200000",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
                "INDEXING ATOMS=ONE_BASED",
                "[BONDS]",
                "1 2",
                "2 3",
                "3 4",
                "4 5",
                "5 6",
                "1 6",
                "[RINGS]",
                "1 SIZE=6 ATOMS=1 2 3 4 5 6",
                "",
                "#SYNTHONS",
                "SCHEMA oracle.xyz.synthons.v1",
                "INDEXING ATOMS=ONE_BASED",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    definition = write_interaction_center_section(xyzin)
    lines = xyzin.read_text(encoding="utf-8").splitlines()
    section = section_content(lines, "INTERACTION_CENTERS")
    read_back = read_interaction_center_definition(xyzin)

    assert definition.centers[-1].kind == "RING_CENTER"
    assert definition.centers[-1].atoms == (1, 2, 3, 4, 5, 6)
    assert any(
        interaction.atom == 7 and interaction.center_id == definition.centers[-1].identifier
        for interaction in definition.interactions
    )
    assert section[0] == "SCHEMA oracle.xyz.interaction_centers.v1"
    assert any("KIND=RING_CENTER" in line for line in section)
    assert any("ATOM=7" in line and "CENTER=" in line for line in section)
    assert read_back.centers[-1].kind == "RING_CENTER"


def test_interaction_centers_for_ferrocene_keep_metal_out_of_ring_centers(tmp_path):
    source = _test_molecule_path("ferrocene.inp")
    xyzin = tmp_path / "ferrocene.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    definition = write_interaction_center_section(xyzin)

    ring_centers = tuple(center for center in definition.centers if center.kind == "RING_CENTER")
    ring_atom_sets = {frozenset(center.atoms) for center in ring_centers}
    ring_center_ids = {center.identifier for center in ring_centers}

    assert ring_atom_sets == {
        frozenset((1, 3, 4, 5, 6)),
        frozenset((7, 8, 9, 10, 11)),
    }
    assert all(2 not in center.atoms for center in ring_centers)
    assert len(definition.interactions) == 2
    assert all(interaction.kind == "ATOM_RING_CENTER" for interaction in definition.interactions)
    assert all(interaction.atom == 2 for interaction in definition.interactions)
    assert {interaction.center_id for interaction in definition.interactions} == ring_center_ids
