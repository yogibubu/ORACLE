from __future__ import annotations

from pathlib import Path

import pytest

from oracle_core import section_content
from oracle_chem import preprocess_to_enriched_xyz, write_validation_section
from oracle_gaussian import read_gaussian_cartesian_input
from oracle_gicforge import (
    GICForgeContractError,
    build_gic_definition_from_xyzin,
    gaussian_gic_lines_from_xyzin,
    write_gicforge_build_sections,
    write_gicforge_gaussian_input,
    write_gicforge_plan_sections,
)


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


def test_gicforge_writes_gaussian_input_after_gic_plan(tmp_path):
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
    xyzin = tmp_path / "water.xyzin"
    gjf = tmp_path / "water.gjf"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_gicforge_plan_sections(xyzin, symmetrize=True)
    write_gicforge_gaussian_input(
        xyzin,
        gjf,
        route="#p b3lyp/sto-3g opt",
        title="water from ORACLE",
    )

    text = gjf.read_text(encoding="utf-8")
    geometry = read_gaussian_cartesian_input(gjf)
    assert "#p b3lyp/sto-3g opt" in text
    assert "! ORACLE_SCHEMA oracle.gaussian.gic_input.v1" in text
    assert geometry.atoms == ("O", "H", "H")
    assert geometry.comment == "water from ORACLE"


def test_gicforge_build_writes_frozen_gics_and_sycart(tmp_path):
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
    xyzin = tmp_path / "water.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True, sycart=True)
    lines = xyzin.read_text(encoding="utf-8").splitlines()
    gic = section_content(lines, "GIC")
    sycart = section_content(lines, "SYCART")

    assert definition.target_rank == 3
    assert definition.rank == 3
    assert gic[0] == "SCHEMA oracle.xyz.gic.v1"
    assert "STATUS BUILT" in gic
    assert "BACKEND oracle-native-primitive.v1" in gic
    assert "SYMMETRY_MODE IDENTITY_C1" in gic
    assert "GIC_COUNT 3" in gic
    assert not any("PENDING" in line for line in gic)
    assert "GIC001 = R(1,2)" in gic
    assert "GIC003 = A(2,1,3)" in gic
    assert sycart[0] == "SCHEMA oracle.xyz.sycart.v1"
    assert "STATUS BUILT" in sycart
    assert "COORD_COUNT 3" in sycart
    assert any(line.startswith("SYC001 IRREP=A COMPONENTS=") for line in sycart)


def test_gicforge_build_definition_uses_saved_topology(tmp_path):
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
    xyzin = tmp_path / "water.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = build_gic_definition_from_xyzin(xyzin)

    assert [primitive.gaussian_expression() for primitive in definition.primitives] == [
        "R(1,2)",
        "R(1,3)",
        "A(2,1,3)",
    ]


def test_gicforge_build_gaussian_input_includes_readgic_block(tmp_path):
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
    xyzin = tmp_path / "water.xyzin"
    gjf = tmp_path / "water.gjf"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_gicforge_build_sections(xyzin)
    write_gicforge_gaussian_input(xyzin, gjf, route="#p hf/sto-3g opt=GIC")

    text = gjf.read_text(encoding="utf-8")
    assert gaussian_gic_lines_from_xyzin(xyzin) == [
        "GIC001 = R(1,2)",
        "GIC002 = R(1,3)",
        "GIC003 = A(2,1,3)",
    ]
    assert "GIC001 = R(1,2)" in text
    assert "GIC003 = A(2,1,3)" in text


def test_gicforge_build_handles_corpus_zmatrix_case(tmp_path):
    source = (
        Path(__file__).resolve().parents[1]
        / "tests"
        / "fixtures"
        / "test_molecules"
        / "molecules"
        / "h2cozmat.inp"
    )
    xyzin = tmp_path / "h2co.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)

    assert definition.target_rank == 6
    assert definition.rank == 6
    assert len(definition.gics) == 6
