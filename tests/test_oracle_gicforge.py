from __future__ import annotations

from pathlib import Path

import pytest

from oracle_core import section_content
from oracle_chem import preprocess_to_enriched_xyz, write_validation_section
from oracle_fragments import write_fragment_build_section
from oracle_gaussian import read_gaussian_cartesian_input
from oracle_gicforge import (
    GICForgeContractError,
    build_gic_b_matrix_from_xyzin,
    build_gic_definition_from_xyzin,
    gaussian_gic_lines_from_xyzin,
    gic_b_matrix_lines,
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


def test_gicforge_build_uses_built_fragments_for_relative_coordinates(tmp_path):
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
    write_validation_section(xyzin)
    write_fragment_build_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    gic = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "GIC")
    primitive_lines = [line for line in gic if line.startswith("P")]
    gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin)

    assert definition.target_rank == 12
    assert definition.rank == 12
    assert len(definition.gics) == 12
    assert any("FAMILY=FRAG_DISTANCE FUNCTION=FC_DIST" in line for line in primitive_lines)
    assert any(
        "FAMILY=FRAG_CENTER_ATOM_DISTANCE FUNCTION=FCA_DIST" in line
        for line in primitive_lines
    )
    assert any("FAMILY=FRAG_ORIENTATION FUNCTION=FROT" in line for line in primitive_lines)
    assert any("REFS=F002,F001" in line for line in primitive_lines)
    assert "F001=Fragment(1-3)" in gaussian_lines
    assert "F002=Fragment(4-6)" in gaussian_lines
    assert "CxF001(Inactive)=XCntr(F001)" in gaussian_lines
    assert any(line.startswith("GIC005 = SQRT((CxF001-CxF002)**2") for line in gaussian_lines)
    assert any(line.startswith("GIC006 = SQRT((CxF002-X(1))**2") for line in gaussian_lines)
    assert any("KxF002F001" in line for line in gaussian_lines)
    assert any(line.startswith("GIC011 = KxF002F001") for line in gaussian_lines)


def test_gicforge_b_matrix_includes_fragment_coordinate_rows(tmp_path):
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
    write_validation_section(xyzin)
    write_fragment_build_section(xyzin)
    write_gicforge_build_sections(xyzin)
    matrix = build_gic_b_matrix_from_xyzin(xyzin)

    assert matrix.coordinate_labels == tuple(f"GIC{idx:03d}" for idx in range(1, 13))
    assert matrix.cartesian_columns[:6] == ("1:X", "1:Y", "1:Z", "2:X", "2:Y", "2:Z")
    assert len(matrix.rows) == 12
    assert len(matrix.rows[0]) == 18
    assert all(value == value for row in matrix.rows for value in row)

    center_distance_row = matrix.rows[matrix.coordinate_labels.index("GIC005")]
    assert center_distance_row[2] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[5] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[8] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[11] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[14] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[17] == pytest.approx(1.0 / 3.0, abs=1.0e-8)

    translation_row = matrix.rows[matrix.coordinate_labels.index("GIC009")]
    assert translation_row[0] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert translation_row[3] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert translation_row[6] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert translation_row[9] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert translation_row[12] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert translation_row[15] == pytest.approx(1.0 / 3.0, abs=1.0e-8)

    lines = gic_b_matrix_lines(matrix)
    assert lines[0] == "SCHEMA oracle.gic.bmatrix.v1"
    assert "ROW_COUNT 12" in lines
    assert "COLUMN_COUNT 18" in lines
    assert any(line.startswith("GIC005 NAME=FCDi0001 IRREP=A VALUES=") for line in lines)
