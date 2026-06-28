from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
from matrix_gaussian import hessian_input_from_gaussian_fchk, lower_to_symmetric, read_gaussian_fchk
from matrix_gf import (
    format_gf_scaling_preview,
    gf_csv_tables,
    gf_scaling_preview_from_xyzin,
    gf_from_cartesian_hessian_and_gic_b_matrix,
    large_amplitude_analysis_from_gf_matrices,
    nonbonded_cartesian_hessian_correction,
    pulay_scaling_factors,
    pulay_scaling_preview,
    read_gf_ped_section,
    run_xyzin_gf_report_from_fchk,
    solve_wilson_gf,
    write_csv_tables,
    write_gf_ped_section_from_report,
)
from matrix_neo import write_gicforge_build_sections


ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ROOT / "tests" / "fixtures" / "test_molecules" / "molecules"
GF_FIXTURES = ROOT / "tests" / "fixtures" / "gf"


def test_gaussian_fchk_hessian_adapter_reads_merlino_blocks():
    data = read_gaussian_fchk(GF_FIXTURES / "h2o.fchk")
    hessian = lower_to_symmetric(data.cartesian_hessian_lower)
    canonical = hessian_input_from_gaussian_fchk(GF_FIXTURES / "h2o.fchk")

    assert data.atomic_numbers.tolist() == [1, 8, 1]
    assert data.cartesian_coordinates_bohr.shape == (3, 3)
    assert hessian.shape == (9, 9)
    assert np.allclose(hessian, hessian.T)
    assert canonical.cartesian_hessian.shape == (9, 9)


def test_wilson_gf_solver_matches_diagonal_reference():
    result = solve_wilson_gf(
        np.diag([4.0, 9.0]),
        np.eye(2),
        scale_to_cm=False,
    )

    assert result.eigenvalues.tolist() == [4.0, 9.0]
    assert result.frequencies_cm.tolist() == [2.0, 3.0]


def test_wilson_gf_solver_routes_dense_diagonalization_through_matrix_core(monkeypatch):
    from matrix_gf import harmonic

    calls = []

    def fake_eigh(matrix):
        calls.append(matrix.shape)
        return np.linalg.eigh(matrix)

    monkeypatch.setattr(harmonic, "eigh_arrays", fake_eigh)

    solve_wilson_gf(np.diag([4.0, 9.0]), np.eye(2), scale_to_cm=False)

    assert calls == [(2, 2), (2, 2)]


def test_gf_can_solve_separated_symmetry_blocks():
    result = gf_from_cartesian_hessian_and_gic_b_matrix(
        np.asarray(
            [
                [4.0, 0.0, 0.0],
                [0.0, 9.0, 0.0],
                [0.0, 0.0, 0.0],
            ],
            dtype=float,
        ),
        np.asarray([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float),
        np.asarray([1.0], dtype=float),
        gic_labels=("GIC001 R(1,2)", "GIC002 A(1,2,3)"),
        gic_names=("A1Str001", "B2Bend001"),
        gic_irreps=("A1", "B2"),
        block_by_irrep=True,
    )

    assert result.matrix_model == "FULL+IRREP_BLOCKS"
    assert result.block_labels == ("A1", "B2")
    assert result.force_constants[0, 1] == 0.0
    assert result.frequencies_cm[0] < result.frequencies_cm[1]
    assert np.allclose(result.ped.values, np.diag([100.0, 100.0]))
    assert np.allclose(np.sum(result.ped.values, axis=0), 100.0)


def test_gf_large_amplitude_subspaces_use_existing_gic_blocks_without_projection():
    result = gf_from_cartesian_hessian_and_gic_b_matrix(
        np.diag([4.0, 9.0, 16.0]),
        np.eye(3),
        np.asarray([1.0], dtype=float),
        gic_labels=("GIC001 R(1,2)", "GIC002 D(1,2,3,4)", "GIC003 U(1,2,3,4)"),
        gic_names=("A1Str001", "A2Tors001", "B1OuPl001"),
        gic_irreps=("A1", "A2", "B1"),
    )
    large = result.large_amplitude

    assert large is not None
    assert [(item.index, item.family) for item in large.coordinates] == [
        (2, "torsion"),
        (3, "oop"),
    ]
    assert [(block.label, block.indices) for block in large.blocks] == [
        ("torsion", (2,)),
        ("oop", (3,)),
        ("all_large_amplitude", (2, 3)),
    ]
    assert np.allclose(large.blocks[0].frequencies_cm, [3.0 * 5140.487143715055])
    assert large.blocks[0].max_f_coupling_to_rest == 0.0
    assert large.mode_contributions[1].ped_percent == pytest.approx(100.0)
    assert large.mode_contributions[2].ped_percent == pytest.approx(100.0)


def test_large_amplitude_coupling_diagnostic_uses_both_f_and_g():
    base_kwargs = {
        "frequencies_cm": np.asarray([100.0, 200.0, 300.0], dtype=float),
        "ped": np.eye(3) * 100.0,
        "gic_labels": ("GIC001 R(1,2)", "GIC002 D(1,2,3,4)", "GIC003 U(1,2,3,4)"),
        "gic_names": ("A1Str001", "A2Tors001", "B1OuPl001"),
        "gic_irreps": ("A1", "A2", "B1"),
    }

    f_dominated = large_amplitude_analysis_from_gf_matrices(
        force_constants=np.asarray(
            [[9.0, 3.0, 0.0], [3.0, 4.0, 0.0], [0.0, 0.0, 1.0]],
            dtype=float,
        ),
        g_matrix=np.eye(3),
        **base_kwargs,
    )
    f_block = next(block for block in f_dominated.blocks if block.label == "torsion")
    assert f_block.relative_f_coupling_to_rest == pytest.approx(3.0 / 9.0)
    assert f_block.relative_g_coupling_to_rest == pytest.approx(0.0)
    assert f_block.relative_fg_coupling_to_rest == pytest.approx(3.0 / 9.0)

    g_dominated = large_amplitude_analysis_from_gf_matrices(
        force_constants=np.diag([9.0, 4.0, 1.0]),
        g_matrix=np.asarray(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.25], [0.0, 0.25, 1.0]],
            dtype=float,
        ),
        **base_kwargs,
    )
    g_block = next(block for block in g_dominated.blocks if block.label == "torsion")
    assert g_block.relative_f_coupling_to_rest == pytest.approx(0.0)
    assert g_block.relative_g_coupling_to_rest == pytest.approx(0.25)
    assert g_block.relative_fg_coupling_to_rest == pytest.approx(0.25)


def test_gf_rejects_cross_irrep_couplings_when_symmetry_blocks_requested():
    with pytest.raises(ValueError, match="not block diagonal"):
        gf_from_cartesian_hessian_and_gic_b_matrix(
            np.asarray(
                [
                    [4.0, 1.0, 0.0],
                    [1.0, 9.0, 0.0],
                    [0.0, 0.0, 0.0],
                ],
                dtype=float,
            ),
            np.asarray([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float),
            np.asarray([1.0], dtype=float),
            gic_labels=("GIC001 R(1,2)", "GIC002 A(1,2,3)"),
            gic_names=("A1Str001", "B2Bend001"),
            gic_irreps=("A1", "B2"),
            block_by_irrep=True,
        )


def test_pulay_scaling_classes_match_multiple_gics_and_reject_mixed_types(tmp_path):
    labels = (
        "GIC001 R(1,2)",
        "GIC002 R(1,3)",
        "GIC003 A(2,1,3)",
        "GIC004 D(2,1,3,4)",
    )
    names = ("A1Str0001", "A1Str0002", "A1Bend0001", "B2Tors0001")
    scale_file = tmp_path / "scale.txt"
    scale_file.write_text(
        "default 1.0\nclass CH_stretches 0.95 R(1,2)|R(1,3)\nclass torsions:0.80:D(\n",
        encoding="utf-8",
    )

    factors = pulay_scaling_factors(
        4,
        labels=labels,
        names=names,
        scale_path=scale_file,
        scale_class_records=("bends:0.90:Bend",),
    )

    assert factors is not None
    assert factors.tolist() == pytest.approx([0.95, 0.95, 0.90, 0.80])
    with pytest.raises(ValueError, match="mixes coordinate types"):
        pulay_scaling_factors(
            4,
            labels=labels,
            names=names,
            scale_class_records=("bad:0.9:R(|A(",),
        )


def test_pulay_scaling_preview_reports_rules_and_final_assignments(tmp_path):
    labels = ("GIC001 R(1,2)", "GIC002 R(1,3)", "GIC003 A(2,1,3)")
    names = ("A1Str0001", "A1Str0002", "A1Bend0001")
    scale_file = tmp_path / "scale.txt"
    scale_file.write_text(
        "default 1.0\nclass CH_stretches 0.95 R(1,2)|R(1,3)\n",
        encoding="utf-8",
    )

    preview = pulay_scaling_preview(
        3,
        labels=labels,
        names=names,
        scale_path=scale_file,
        scale_class_records=("bends:0.90:Bend",),
    )
    text = format_gf_scaling_preview(preview)

    assert preview.changed_count == 3
    assert [rule.name for rule in preview.rules] == ["default", "CH_stretches", "bends"]
    assert preview.rules[1].matches == (1, 2)
    assert preview.rules[1].family == "stretch"
    assert [(item.identifier, item.factor, item.source) for item in preview.assignments] == [
        ("GIC001", pytest.approx(0.95), "class CH_stretches"),
        ("GIC002", pytest.approx(0.95), "class CH_stretches"),
        ("GIC003", pytest.approx(0.90), "class bends"),
    ]
    assert "GF/PED Pulay scaling preview" in text
    assert "class    CH_stretches" in text
    assert "GIC003" in text


def test_gf_scaling_preview_reads_frozen_xyzin_gics(tmp_path):
    source = MOLECULES / "h2ocart.inp"
    xyzin = tmp_path / "h2o.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_gicforge_build_sections(xyzin)

    preview = gf_scaling_preview_from_xyzin(
        xyzin,
        scale_class_records=("stretches:0.97:Str",),
    )
    text = format_gf_scaling_preview(preview)

    assert preview.xyzin_path == xyzin
    assert preview.assignments
    assert any(item.source == "class stretches" for item in preview.assignments)
    assert any(item.factor == pytest.approx(0.97) for item in preview.assignments)
    assert "Frozen xyzin:" in text


def test_nonbonded_correction_excludes_12_13_and_scales_14_terms():
    coords = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    atomic_numbers = np.asarray([1, 1, 1, 1], dtype=int)
    charges = np.asarray([1.0, -1.0, 1.0, -1.0], dtype=float)
    bonds = ((1, 2), (2, 3), (3, 4))

    zero_14 = nonbonded_cartesian_hessian_correction(
        coords,
        atomic_numbers,
        bonds,
        charges=charges,
        electrostatic=True,
        uff_vdw=True,
        one_four_scale=0.0,
    )
    half_14 = nonbonded_cartesian_hessian_correction(
        coords,
        atomic_numbers,
        bonds,
        charges=charges,
        electrostatic=True,
        uff_vdw=True,
        one_four_scale=0.5,
    )
    full_14 = nonbonded_cartesian_hessian_correction(
        coords,
        atomic_numbers,
        bonds,
        charges=charges,
        electrostatic=True,
        uff_vdw=True,
        one_four_scale=1.0,
    )

    assert np.allclose(zero_14, 0.0)
    assert np.allclose(half_14, 0.5 * full_14)
    assert np.allclose(half_14, half_14.T)


def test_xyzin_gf_report_runs_from_fchk_and_frozen_gics(tmp_path):
    source = MOLECULES / "h2ocart.inp"
    xyzin = tmp_path / "h2o.xyzin"
    scale = tmp_path / "scale.txt"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    scale.write_text("GIC003 0.90\n", encoding="utf-8")

    report = run_xyzin_gf_report_from_fchk(
        GF_FIXTURES / "h2o.fchk",
        xyzin,
        scale_path=scale,
        local=True,
        force_threshold=1.0e-12,
    )
    written = write_csv_tables(report, tmp_path / "csv", prefix="gic_gf")
    section = write_gf_ped_section_from_report(
        xyzin,
        report,
        source_kind="fchk",
        source_path=GF_FIXTURES / "h2o.fchk",
        report_path=tmp_path / "gf.report",
        csv_dir=tmp_path / "csv",
    )
    reread = read_gf_ped_section(xyzin)

    assert "GF/PED from ORACLE non-redundant GICs" in report.text
    assert "Frozen xyzin:" in report.text
    assert "Pulay Hessian scaling: applied" in report.text
    assert "Matrix model: LOCAL" in report.text
    assert "Force-constant threshold:" in report.text
    assert len(report.result.frequencies_cm) == definition.rank
    assert np.all(np.isfinite(report.result.frequencies_cm))
    assert report.result.ped.values.shape == (definition.rank, definition.rank)
    assert (tmp_path / "csv" / "gic_gf_frequencies.csv").is_file()
    assert "ped.csv" in written
    assert "large_amplitude_coordinates.csv" in written
    assert "large_amplitude_blocks.csv" in written
    assert section.source_kind == "fchk"
    assert reread.source_path == GF_FIXTURES / "h2o.fchk"
    assert reread.report_path == tmp_path / "gf.report"
    assert len(reread.modes) == definition.rank
    assert len(reread.gics) == definition.rank
    assert len(reread.gics[0].ped) == definition.rank
    assert reread.large_amplitude_blocks == section.large_amplitude_blocks


def test_gf_report_and_csv_include_large_amplitude_ring_coordinates(tmp_path):
    source = MOLECULES / "pyrrole.inp"
    xyzin = tmp_path / "pyrrole.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_gicforge_build_sections(xyzin, symmetrize=True)

    report = run_xyzin_gf_report_from_fchk(MOLECULES / "pyrrole.fchk", xyzin)
    tables = gf_csv_tables(report)
    section = write_gf_ped_section_from_report(xyzin, report, source_kind="fchk")
    reread = read_gf_ped_section(xyzin)
    large = report.result.large_amplitude

    assert large is not None
    assert large.coordinate_count > 0
    assert any(coordinate.family == "ring_puckering" for coordinate in large.coordinates)
    assert any(coordinate.family == "oop" for coordinate in large.coordinates)
    assert "Large-amplitude block GF frequencies" in report.text
    assert "no projection" in report.text
    assert "large_amplitude_mode_ped.csv" in tables
    assert section.large_amplitude_blocks
    assert reread.large_amplitude_coordinates == section.large_amplitude_coordinates
    assert any(block.family == "ring_puckering" for block in reread.large_amplitude_blocks)


def test_xyzin_gf_report_can_use_symmetry_blocks(tmp_path):
    source = MOLECULES / "h2ocart.inp"
    xyzin = tmp_path / "h2o.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_gicforge_build_sections(xyzin, symmetrize=True)

    report = run_xyzin_gf_report_from_fchk(
        GF_FIXTURES / "h2o.fchk",
        xyzin,
        block_by_irrep=True,
    )

    assert "IRREP_BLOCKS" in report.result.matrix_model
    assert report.result.block_labels


def test_xyzin_gf_report_subtracts_nonbonded_terms_from_gaussian_hessian(tmp_path):
    source = MOLECULES / "pyrrole.inp"
    xyzin = tmp_path / "pyrrole.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)

    report = run_xyzin_gf_report_from_fchk(
        MOLECULES / "pyrrole.fchk",
        xyzin,
        subtract_electrostatic=True,
        subtract_uff_vdw=True,
    )

    assert "Cartesian Hessian correction:" in report.text
    assert "1-4 scale=0.5" in report.text
    assert "Synthons electronegativity model" in report.result.hessian_correction
    assert len(report.result.frequencies_cm) == definition.rank
    assert np.all(np.isfinite(report.result.frequencies_cm))
