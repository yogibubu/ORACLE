from __future__ import annotations

import pytest

from oracle_core import section_content
from oracle_chem import preprocess_to_enriched_xyz
from oracle_fragments import (
    FragmentContractError,
    read_fragment_records,
    write_fragment_build_section,
    write_fragment_plan_section,
)


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
    assert fragments[0] == "SCHEMA oracle.xyz.fragments.v1"
    assert (
        "DEPENDENCIES TOPOLOGY=oracle.xyz.topology.v1 SYNTHONS=oracle.xyz.synthons.v1"
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
