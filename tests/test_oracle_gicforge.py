from __future__ import annotations

import pytest

from oracle_core import section_content
from oracle_chem import preprocess_to_enriched_xyz, write_validation_section
from oracle_gicforge import GICForgeContractError, write_gicforge_plan_sections


def test_gicforge_plan_requires_validation_pass(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#VALIDATION",
                "SCHEMA oracle.xyz.validation.v1",
                "STATUS FAIL",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(GICForgeContractError, match="status must be PASS"):
        write_gicforge_plan_sections(path)


def test_gicforge_plan_writes_gic_and_sycart_sections(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#VALIDATION",
                "SCHEMA oracle.xyz.validation.v1",
                "STATUS PASS",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    write_gicforge_plan_sections(path, symmetrize=True, sycart=True)
    lines = path.read_text(encoding="utf-8").splitlines()

    gic = section_content(lines, "GIC")
    sycart = section_content(lines, "SYCART")
    assert gic[0] == "SCHEMA oracle.xyz.gic.v1"
    assert "STATUS PLANNED" in gic
    assert "SYMMETRIZE TRUE" in gic
    assert "PENDING GICFORGE_IMPLEMENTATION" in gic
    assert sycart[0] == "SCHEMA oracle.xyz.sycart.v1"


def test_preprocess_validate_then_gicforge_plan_pipeline(tmp_path):
    source = tmp_path / "water.xyz"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    target = tmp_path / "water.xyzin"

    preprocess_to_enriched_xyz(source, target)
    validation = write_validation_section(target)
    write_gicforge_plan_sections(target, symmetrize=True)
    lines = target.read_text(encoding="utf-8").splitlines()

    assert validation.status == "PASS"
    assert "STATUS PASS" in section_content(lines, "VALIDATION")
    assert "SYMMETRIZE TRUE" in section_content(lines, "GIC")
