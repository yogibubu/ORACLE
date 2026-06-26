from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from oracle_core import section_content
from oracle_chem import preprocess_to_enriched_xyz, write_validation_section
from oracle_fragments import write_fragment_build_section
from oracle_gaussian import read_gaussian_cartesian_input
from oracle_gicforge import (
    FrozenGIC,
    GICForgeContractError,
    GICDefinition,
    GICPointGroupOperation,
    GICPrimitive,
    build_gic_b_matrix_from_xyzin,
    build_gic_definition_from_xyzin,
    gaussian_gic_lines_from_xyzin,
    gic_report_from_xyzin,
    gic_b_matrix_lines,
    gic_definition_section_lines,
    special_symmetry_source_blocks,
    symmetrize_gic_definition,
    total_symmetric_gic_names,
    write_gicforge_build_sections,
    write_gicforge_gaussian_input,
    write_gicforge_plan_sections,
)
from oracle_gicforge.definition import (
    _analytic_b_row,
    _finite_difference_b_row,
    _select_ranked_primitives,
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
    assert "POINT_GROUP C2v" in gic
    assert "SYMMETRY_MODE POINT_GROUP_PROJECTOR" in gic
    assert "SYMMETRY_GROUP C2v" in gic
    assert "TOTAL_SYMMETRIC_IRREP A1" in gic
    assert "TOTAL_SYMMETRIC_GIC_COUNT 2" in gic
    assert "TOTAL_SYMMETRIC_GICS A1Str001,A1Bend001" in gic
    assert "GIC_COUNT 3" in gic
    assert not any("PENDING" in line for line in gic)
    assert "[SYMMETRY_DIAGNOSTICS]" in gic
    assert "METHOD POINT_GROUP_PROJECTOR" in gic
    assert "STATUS APPLIED" in gic
    assert any(
        line.startswith(
            "GROUP 1 BLOCK=STRETCH FAMILY=STRETCH "
            "SIGNATURE=OPS=E,sigma_yz,sigma_xy,C2y^1 "
            "SOURCES=Str0001,Str0002 OUTPUTS=A1Str001,B2Str001"
        )
        for line in gic
    )
    assert any(
        line.startswith(
            "GIC001 NAME=A1Str001 FAMILY=STRETCH IRREP=A1 "
            "COEFFS=P001:0.707106781187,P002:0.707106781187"
        )
        for line in gic
    )
    assert "A1Str001 = 0.707106781187*(R(1,2))+0.707106781187*(R(1,3))" in gic
    assert "B2Str001 = 0.707106781187*(R(1,2))-0.707106781187*(R(1,3))" in gic
    assert "A1Bend001 = A(2,1,3)" in gic
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.status == "APPLIED"
    assert definition.symmetry_diagnostics.symmetry_group == "C2v"
    assert definition.symmetry_diagnostics.total_symmetric_irrep == "A1"
    assert definition.symmetry_diagnostics.total_symmetric_gics == (
        "A1Str001",
        "A1Bend001",
    )
    assert total_symmetric_gic_names(definition) == (
        "A1Str001",
        "A1Bend001",
    )
    assert definition.gics[0].name == "A1Str001"
    assert definition.gics[0].coefficients[0][1] == pytest.approx(1.0 / np.sqrt(2.0))
    matrix = build_gic_b_matrix_from_xyzin(xyzin)
    assert matrix.coordinate_names[:2] == ("A1Str001", "B2Str001")
    assert matrix.rows[0][5] == pytest.approx(1.0 / np.sqrt(2.0))
    assert matrix.rows[0][7] == pytest.approx(1.0 / np.sqrt(2.0))
    assert matrix.rows[1][5] == pytest.approx(1.0 / np.sqrt(2.0))
    assert matrix.rows[1][7] == pytest.approx(-1.0 / np.sqrt(2.0))
    report_lines = gic_report_from_xyzin(xyzin)
    assert "Method: POINT_GROUP_PROJECTOR" in report_lines
    assert "Total irrep: A1" in report_lines
    assert (
        "STRETCH OPS=E,sigma_yz,sigma_xy,C2y^1: "
        "Str0001,Str0002 -> A1Str001,B2Str001"
    ) in report_lines
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


def test_gicforge_native_analytic_b_rows_match_diagnostic_finite_difference():
    coords = np.asarray(
        [
            [0.00, 0.00, 0.00],
            [1.20, 0.00, 0.00],
            [1.80, 1.05, 0.20],
            [2.55, 1.20, 1.05],
        ],
        dtype=float,
    )
    linear_coords = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.001, 0.0],
        ],
        dtype=float,
    )
    primitives = [
        (GICPrimitive("P001", "Str", "STRETCH", "R", (1, 2)), coords),
        (GICPrimitive("P002", "Bend", "BEND", "A", (1, 2, 3)), coords),
        (GICPrimitive("P003", "Tors", "TORSION", "D", (1, 2, 3, 4)), coords),
        (GICPrimitive("P004", "OuPl", "OUT_OF_PLANE", "U", (2, 1, 3, 4)), coords),
        (
            GICPrimitive("P005", "LinB", "LINEAR_BEND", "L", (1, 2, 3), mode=-1),
            linear_coords,
        ),
    ]

    for primitive, primitive_coords in primitives:
        analytic = _analytic_b_row(primitive, primitive_coords)
        diagnostic_fd = _finite_difference_b_row(
            primitive,
            primitive_coords,
            step_angstrom=1.0e-6,
        )
        assert np.max(np.abs(analytic - diagnostic_fd)) < 1.0e-5


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
    assert "RANK_METHOD analytic_b_matrix_mgs_greedy" in gic
    assert (
        "REDUCTION_POLICY SPECIAL_PROTECTED_FIRST_THEN_ORDINARY_ANALYTIC_RANK"
        in gic
    )
    assert "[REDUCTION_DIAGNOSTICS]" in gic
    assert any(line.startswith("SELECTED P") for line in gic)
    assert any(
        "FAMILY=FRAG_DISTANCE CLASS=SPECIAL_PROTECTED FUNCTION=FC_DIST" in line
        for line in primitive_lines
    )
    assert any(
        "FAMILY=FRAG_CENTER_ATOM_DISTANCE CLASS=SPECIAL_PROTECTED FUNCTION=FCA_DIST"
        in line
        for line in primitive_lines
    )
    assert any(
        "FAMILY=FRAG_ORIENTATION CLASS=SPECIAL_PROTECTED FUNCTION=FROT" in line
        for line in primitive_lines
    )
    assert any("REFS=F002,F001" in line for line in primitive_lines)
    assert "F001=Fragment(1-3)" in gaussian_lines
    assert "F002=Fragment(4-6)" in gaussian_lines
    assert "CxF001(Inactive)=XCntr(F001)" in gaussian_lines
    assert any("= SQRT((CxF001-CxF002)**2" in line for line in gaussian_lines)
    assert any("= SQRT((CxF002-X(1))**2" in line for line in gaussian_lines)
    assert any("KxF002F001" in line for line in gaussian_lines)
    assert any("ExF002F001" in line for line in gaussian_lines)
    assert any("1.0D-24" in line and line.startswith("KnF002F001") for line in gaussian_lines)
    assert any(" = ExF002F001" in line for line in gaussian_lines)


def test_gicforge_symmetrizes_special_fragment_coordinates(tmp_path):
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
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    matrix = build_gic_b_matrix_from_xyzin(xyzin)
    gic = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "GIC")

    assert definition.symmetry_diagnostics is not None
    special_groups = [
        group
        for group in definition.symmetry_diagnostics.groups
        if group.block.startswith("SPECIAL_")
    ]
    assert special_groups
    assert special_groups[0].block == "SPECIAL_FRAGMENT_CENTER_ATOM"
    assert special_groups[0].output_gics == ("A1FCAtS001", "A2FCAtD001")
    assert "TOTAL_SYMMETRIC_IRREP A1" in gic
    assert "A1FCAtS001" in total_symmetric_gic_names(definition)
    assert "A2FCAtD001" not in total_symmetric_gic_names(definition)
    assert any(
        line.startswith(
            "GIC003 NAME=A1FCAtS001 FAMILY=FRAG_CENTER_ATOM_DISTANCE IRREP=A1 "
            "COEFFS=P007:0.707106781187,P008:0.707106781187"
        )
        for line in gic
    )
    assert any(line.startswith("A1FCAtS001 = 0.707106781187*(") for line in gic)
    assert matrix.coordinate_names[2:4] == ("A1FCAtS001", "A2FCAtD001")
    assert len(matrix.rows) == definition.rank


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
    definition = write_gicforge_build_sections(xyzin)
    matrix = build_gic_b_matrix_from_xyzin(xyzin)

    assert matrix.coordinate_labels == tuple(f"GIC{idx:03d}" for idx in range(1, 13))
    assert matrix.cartesian_columns[:6] == ("1:X", "1:Y", "1:Z", "2:X", "2:Y", "2:Z")
    assert len(matrix.rows) == 12
    assert len(matrix.rows[0]) == 18
    assert all(value == value for row in matrix.rows for value in row)

    center_distance_row = matrix.rows[matrix.coordinate_names.index("FCDi0001")]
    assert center_distance_row[2] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[5] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[8] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[11] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[14] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[17] == pytest.approx(1.0 / 3.0, abs=1.0e-8)

    translation_row = matrix.rows[matrix.coordinate_names.index("FTrn0001")]
    assert translation_row[0] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert translation_row[3] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert translation_row[6] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert translation_row[9] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert translation_row[12] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert translation_row[15] == pytest.approx(1.0 / 3.0, abs=1.0e-8)

    lines = gic_b_matrix_lines(matrix)
    assert lines[0] == "SCHEMA oracle.gic.bmatrix.v1"
    assert "BACKEND oracle-native-analytic-bmatrix.v1" in lines
    assert "DERIVATIVE_MODE ANALYTIC" in lines
    assert "ROW_COUNT 12" in lines
    assert "COLUMN_COUNT 18" in lines
    center_distance_label = matrix.coordinate_labels[
        matrix.coordinate_names.index("FCDi0001")
    ]
    assert any(
        line.startswith(f"{center_distance_label} NAME=FCDi0001 IRREP=UNASSIGNED VALUES=")
        for line in lines
    )

    coords = np.asarray(definition.reference_coordinates_angstrom, dtype=float)
    frot = next(primitive for primitive in definition.primitives if primitive.function == "FROT")
    analytic = _analytic_b_row(frot, coords)
    diagnostic_fd = _finite_difference_b_row(frot, coords, step_angstrom=1.0e-6)
    assert np.max(np.abs(analytic - diagnostic_fd)) < 1.0e-5

    blocks = special_symmetry_source_blocks(definition)
    assert {block.block for block in blocks} >= {
        "SPECIAL_FRAGMENT_DISTANCE",
        "SPECIAL_FRAGMENT_CENTER_ATOM",
        "SPECIAL_FRAGMENT_TRANSLATION",
        "SPECIAL_FRAGMENT_ORIENTATION",
    }

    report_lines = gic_report_from_xyzin(xyzin)
    assert "ORACLE GICForge Report" in report_lines
    assert "Rank method: analytic_b_matrix_mgs_greedy" in report_lines
    assert any("SPECIAL_FRAGMENT_DISTANCE" in line for line in report_lines)


def test_gicforge_center_atom_distance_uses_analytic_b_and_gaussian_center():
    primitive = GICPrimitive(
        identifier="P001",
        name="CnAt0001",
        family="CENTER_ATOM_DISTANCE",
        function="CENTER_ATOM_DIST",
        atoms=(1, 2),
        ref_atoms=(3,),
        refs=("C001", "A3"),
    )
    coords = np.asarray(
        [
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 2.0],
        ],
        dtype=float,
    )

    analytic = _analytic_b_row(primitive, coords)
    diagnostic_fd = _finite_difference_b_row(primitive, coords, step_angstrom=1.0e-6)
    assert np.max(np.abs(analytic - diagnostic_fd)) < 1.0e-8
    assert analytic[2] == pytest.approx(-0.5)
    assert analytic[5] == pytest.approx(-0.5)
    assert analytic[8] == pytest.approx(1.0)

    definition = GICDefinition(
        backend="test",
        point_group="C1",
        symmetrize=False,
        target_rank=1,
        rank=1,
        candidate_count=1,
        reference_coordinates_angstrom=tuple(tuple(row) for row in coords),
        primitives=(primitive,),
        gics=(
            FrozenGIC(
                identifier="GIC001",
                name="CnAt0001",
                family="CENTER_ATOM_DISTANCE",
                irrep="A",
                primitive_id="P001",
                gaussian_expression="NONE",
                coefficients=(("P001", 1.0),),
            ),
        ),
    )
    gaussian_lines = gic_definition_section_lines(definition)
    assert "CxC001(Inactive)=(X(1)+X(2))/2" in gaussian_lines
    assert any(line.startswith("GIC001 = SQRT((CxC001-X(3))**2") for line in gaussian_lines)


def test_gicforge_rank_reduction_protects_special_coordinates():
    ordinary = GICPrimitive(
        identifier="P001",
        name="Str0001",
        family="STRETCH",
        function="R",
        atoms=(1, 2),
    )
    special = GICPrimitive(
        identifier="P002",
        name="CnAt0001",
        family="CENTER_ATOM_DISTANCE",
        function="CENTER_ATOM_DIST",
        atoms=(1,),
        ref_atoms=(2,),
        refs=("C001", "A2"),
    )
    coords = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )

    selected, rank = _select_ranked_primitives(
        (ordinary, special),
        coords,
        target_rank=1,
        rank_tolerance=1.0e-7,
    )

    assert rank == 1
    assert selected == (special,)
    assert selected[0].reduction_class == "SPECIAL_PROTECTED"


def test_gicforge_symmetry_labels_use_point_group_irreps():
    primitives = tuple(
        GICPrimitive(
            identifier=f"P{idx:03d}",
            name=f"Str{idx:04d}",
            family="STRETCH",
            function="R",
            atoms=(1, idx + 1),
        )
        for idx in range(1, 5)
    )
    definition = GICDefinition(
        backend="test",
        point_group="C2v",
        symmetrize=False,
        target_rank=4,
        rank=4,
        candidate_count=4,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 5,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("O", "H", "H", "H", "H"),
    )

    assert [gic.irrep for gic in symmetrized.gics] == ["A1", "A2", "B1", "B2"]
    assert [gic.name for gic in symmetrized.gics] == [
        "A1StrS001",
        "A2StrD001",
        "B1StrD001",
        "B2StrD001",
    ]
    assert total_symmetric_gic_names(symmetrized) == ("A1StrS001",)


def test_gicforge_point_group_projector_uses_operations_without_type_mixing():
    primitives = (
        GICPrimitive("P001", "Str0001", "STRETCH", "R", (1, 2)),
        GICPrimitive("P002", "Str0002", "STRETCH", "R", (1, 3)),
        GICPrimitive("P003", "Bend0001", "BEND", "A", (2, 1, 3)),
    )
    definition = GICDefinition(
        backend="test",
        point_group="C2v",
        symmetrize=False,
        target_rank=3,
        rank=3,
        candidate_count=3,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 3,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )
    operations = (
        GICPointGroupOperation(
            "E",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (1, 2, 3),
        ),
        GICPointGroupOperation(
            "sigma_yz",
            ((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (1, 3, 2),
        ),
        GICPointGroupOperation(
            "sigma_xy",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)),
            (1, 2, 3),
        ),
        GICPointGroupOperation(
            "C2y^1",
            ((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)),
            (1, 3, 2),
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("O", "H", "H"),
        symmetry_operations=operations,
    )

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert [gic.name for gic in symmetrized.gics] == [
        "A1Str001",
        "B2Str001",
        "A1Bend001",
    ]
    assert [gic.irrep for gic in symmetrized.gics] == ["A1", "B2", "A1"]
    assert [
        (group.block, group.family, group.output_gics)
        for group in symmetrized.symmetry_diagnostics.groups
    ] == [
        ("STRETCH", "STRETCH", ("A1Str001", "B2Str001")),
        ("BEND", "BEND", ("A1Bend001",)),
    ]
    assert total_symmetric_gic_names(symmetrized) == ("A1Str001", "A1Bend001")


def test_gicforge_projector_symmetrizes_special_fragment_center_atom_coordinates():
    primitives = (
        GICPrimitive(
            "P001",
            "FCAt0001",
            "FRAG_CENTER_ATOM_DISTANCE",
            "FCA_DIST",
            (1, 2),
            ref_atoms=(3,),
        ),
        GICPrimitive(
            "P002",
            "FCAt0002",
            "FRAG_CENTER_ATOM_DISTANCE",
            "FCA_DIST",
            (4, 5),
            ref_atoms=(6,),
        ),
    )
    definition = GICDefinition(
        backend="test",
        point_group="C2",
        symmetrize=False,
        target_rank=2,
        rank=2,
        candidate_count=2,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 6,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )
    operations = (
        GICPointGroupOperation(
            "E",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (1, 2, 3, 4, 5, 6),
        ),
        GICPointGroupOperation(
            "C2z^1",
            ((-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
            (4, 5, 6, 1, 2, 3),
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("O", "H", "H", "O", "H", "H"),
        symmetry_operations=operations,
    )

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert [gic.name for gic in symmetrized.gics] == ["AFCAt001", "BFCAt001"]
    assert [gic.irrep for gic in symmetrized.gics] == ["A", "B"]
    assert symmetrized.symmetry_diagnostics.groups[0].block == (
        "SPECIAL_FRAGMENT_CENTER_ATOM"
    )
    assert total_symmetric_gic_names(symmetrized) == ("AFCAt001",)
