from __future__ import annotations

from pathlib import Path

import numpy as np

from oracle_chem import preprocess_to_enriched_xyz, write_validation_section
from oracle_gaussian import hessian_input_from_gaussian_fchk, lower_to_symmetric, read_gaussian_fchk
from oracle_gf import (
    gf_from_cartesian_hessian_and_gic_b_matrix,
    nonbonded_cartesian_hessian_correction,
    read_gf_ped_section,
    run_xyzin_gf_report_from_fchk,
    solve_wilson_gf,
    write_csv_tables,
    write_gf_ped_section_from_report,
)
from oracle_gicforge import write_gicforge_build_sections


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


def test_wilson_gf_solver_routes_dense_diagonalization_through_oracle_core(monkeypatch):
    from oracle_gf import harmonic

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

    assert result.matrix_model == "FULL+IRREP_BLOCKS"
    assert result.block_labels == ("A1", "B2")
    assert result.force_constants[0, 1] == 0.0
    assert result.frequencies_cm[0] < result.frequencies_cm[1]


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
    assert section.source_kind == "fchk"
    assert reread.source_path == GF_FIXTURES / "h2o.fchk"
    assert reread.report_path == tmp_path / "gf.report"
    assert len(reread.modes) == definition.rank
    assert len(reread.gics) == definition.rank
    assert len(reread.gics[0].ped) == definition.rank


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
