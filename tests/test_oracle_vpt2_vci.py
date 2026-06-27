from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from oracle_gf import gf_from_gaussian_fchk_with_oracle_gics, gf_from_hessian_input_with_oracle_gics
from oracle_vpt2_vci import (
    AnharmonicInput,
    QuarticForceField,
    ScientificValidationError,
    VCIOptions,
    collect_vpt2_vci_outputs,
    compare_vpt2_vci,
    davidson_lowest,
    generate_vibrational_basis,
    inventory_vpt2_vci_backends,
    load_force_field,
    lower_to_symmetric,
    read_vpt2_vci_section,
    refresh_vpt2_vci_section,
    read_gaussian_fchk_qff,
    read_indexed_qff_text,
    run_python_vci_from_gaussian_fchk,
    run_vpt2_vci_report,
    solve_vci,
    solve_vci_from_anharmonic_input,
    solve_vpt2_from_anharmonic_input,
    validate_force_field,
    vpt2_vci_section_from_run,
    write_vpt2_vci_section,
    zero_anharmonic_force_field,
)
from oracle_vpt2_vci.gaussian_qff import hessian_input_from_gaussian_fchk


ROOT = Path(__file__).resolve().parents[1]
FCHK = ROOT / "tests" / "fixtures" / "gf" / "h2o.fchk"


def test_vci_harmonic_basis_and_energies_are_deterministic():
    qff = zero_anharmonic_force_field(np.array([100.0, 200.0]))
    result = solve_vci(qff, max_quanta=1)

    assert result.basis == ((0, 0), (1, 0), (0, 1))
    assert np.allclose(result.energies_cm, [150.0, 250.0, 350.0])


def test_vci_dense_solver_routes_diagonalization_through_oracle_core(monkeypatch):
    from oracle_vpt2_vci import vci

    calls = []

    def fake_eigh(matrix):
        calls.append(matrix.shape)
        return np.linalg.eigh(matrix)

    monkeypatch.setattr(vci, "eigh_arrays", fake_eigh)

    qff = zero_anharmonic_force_field(np.array([100.0, 200.0]))
    solve_vci(qff, max_quanta=1)

    assert calls == [(3, 3)]


def test_vci_quartic_term_shifts_oscillator_ground_state():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([100.0]),
        cubic_cm={},
        quartic_cm={(0, 0, 0, 0): 4.0},
    )
    result = solve_vci(qff, max_quanta=0)

    assert np.allclose(result.energies_cm, [53.0])


def test_gaussian_fchk_qff_reader_uses_real_anharmonic_blocks():
    data = read_gaussian_fchk_qff(FCHK)

    assert data.atomic_numbers.tolist() == [1, 8, 1]
    assert data.cartesian_coordinates_bohr.shape == (3, 3)
    assert data.masses_amu.shape == (3,)
    assert data.cartesian_hessian_lower.shape == (45,)
    assert np.allclose(data.anharmonic_frequencies_cm[:3], [2123.50470, 4016.61987, 4266.73074])
    assert lower_to_symmetric(data.cartesian_hessian_lower).shape == (9, 9)
    hessian_input = data.to_hessian_input()
    assert hessian_input.source == "gaussian-fchk"
    assert hessian_input.cartesian_hessian.shape == (9, 9)


def test_python_fchk_workflow_produces_vci_levels():
    run = run_python_vci_from_gaussian_fchk(FCHK, max_quanta=1)

    assert run.gf.frequencies_cm.size == 9
    assert len(run.vci.basis) == len(generate_vibrational_basis(3, 1))
    assert run.vci.energies_cm[0] > 0.0


def test_indexed_qff_text_feeds_quartic_vci(tmp_path):
    qff_file = tmp_path / "field.qff"
    qff_file.write_text("FREQ 1 100.0\nQUARTIC 1 1 1 1 4.0\n", encoding="utf-8")

    qff = read_indexed_qff_text(qff_file)
    result = solve_vci(qff, max_quanta=0)

    assert np.allclose(result.energies_cm, [53.0])


def test_independent_davidson_matches_dense_symmetric_diagonalization():
    matrix = np.array(
        [
            [2.0, 0.1, 0.0, 0.0],
            [0.1, 3.0, 0.2, 0.0],
            [0.0, 0.2, 5.0, 0.3],
            [0.0, 0.0, 0.3, 8.0],
        ]
    )
    expected, _ = np.linalg.eigh(matrix)

    result = davidson_lowest(lambda vector: matrix @ vector, np.diag(matrix), n_roots=2)

    assert result.converged
    assert np.allclose(result.eigenvalues, expected[:2], atol=1.0e-8)


def test_vci_davidson_path_matches_dense_path():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([100.0, 140.0]),
        cubic_cm={(0, 0, 1): 1.5},
        quartic_cm={(0, 0, 0, 0): 0.2, (0, 1, 1, 1): -0.1},
    )

    dense = solve_vci(qff, max_quanta=3, n_roots=3)
    iterative = solve_vci(qff, max_quanta=3, n_roots=3, method="davidson")

    assert iterative.davidson is not None
    assert iterative.davidson.converged
    assert np.allclose(iterative.energies_cm, dense.energies_cm, atol=1.0e-7)


def test_canonical_anharmonic_input_runs_vci_and_shifts_levels():
    harmonic = solve_vci_from_anharmonic_input(
        AnharmonicInput(
            harmonic_frequencies_cm=np.array([1000.0, 1500.0]),
            anharmonic_frequencies_cm=np.array([]),
            cubic_cm={},
            quartic_cm={},
        ),
        max_quanta=2,
        n_roots=4,
    )
    anharmonic = solve_vci_from_anharmonic_input(
        AnharmonicInput(
            harmonic_frequencies_cm=np.array([1000.0, 1500.0]),
            anharmonic_frequencies_cm=np.array([]),
            cubic_cm={(0, 0, 1): -25.0},
            quartic_cm={(0, 0, 0, 0): 8.0, (0, 0, 1, 1): -3.0},
        ),
        max_quanta=2,
        n_roots=4,
    )

    assert np.allclose(harmonic.excitation_energies_cm[:4], [0.0, 1000.0, 1500.0, 2000.0])
    assert not np.allclose(anharmonic.excitation_energies_cm, harmonic.excitation_energies_cm)
    assert np.all(anharmonic.excitation_energies_cm[1:] > 0.0)


def test_vci_basis_cutoff_frequency_window_pruning_blocks_and_contributions():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([500.0, 1000.0, 1800.0]),
        cubic_cm={(0, 0, 1): 0.001, (1, 1, 1): 20.0},
        quartic_cm={(1, 1, 1, 1): 5.0, (0, 0, 1, 1): 0.0001},
    )
    options = VCIOptions(
        frequency_min_cm=800.0,
        frequency_max_cm=1600.0,
        basis_energy_cutoff_cm=2200.0,
        force_constant_threshold_cm=0.01,
        mode_symmetries=("a", "b", "c"),
        separate_symmetry_blocks=True,
        coefficient_threshold=0.05,
    )

    result = solve_vci(qff, max_quanta=4, n_roots=4, options=options)

    assert result.basis == ((0,), (1,), (2,))
    assert {block.label for block in result.blocks} == {"A", "b"}
    assert len(result.state_contributions) == 3
    assert np.allclose(result.state_contributions[1].mode_quanta, [1.0], atol=1.0e-8)
    assert result.state_contributions[1].dominant_basis_states[0][0] == (1,)


def test_vci_per_mode_and_excitation_class_limits():
    basis = generate_vibrational_basis(
        3,
        max_quanta=4,
        frequencies_cm=np.array([100.0, 200.0, 300.0]),
        mode_max_quanta=(3, 2, 1),
        excitation_class_limits={1: (1, 2), 2: (2, 3), 3: (3, 3)},
    )

    assert (0, 0, 0) in basis
    assert (3, 0, 0) not in basis
    assert (2, 2, 0) not in basis
    assert (2, 1, 0) in basis
    assert (1, 1, 1) in basis
    assert (2, 1, 1) not in basis
    assert all(state[1] <= 2 and state[2] <= 1 for state in basis)


def test_vpt2_quartic_first_order_matches_ground_state_shift():
    result = solve_vpt2_from_anharmonic_input(
        AnharmonicInput(
            harmonic_frequencies_cm=np.array([100.0]),
            anharmonic_frequencies_cm=np.array([]),
            cubic_cm={},
            quartic_cm={(0, 0, 0, 0): 4.0},
        ),
        max_quanta=0,
    )

    assert np.allclose(result.energies_cm, [53.0])
    assert np.allclose(result.states[0].first_order_cm, 3.0)


def test_vpt2_vci_comparison_uses_same_reduced_mode_selection():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([500.0, 1000.0, 1500.0]),
        cubic_cm={(1, 1, 2): -2.0, (0, 1, 1): 100.0},
        quartic_cm={(1, 1, 1, 1): 0.8, (2, 2, 2, 2): 0.5, (0, 0, 1, 1): 200.0},
    )
    options = VCIOptions(
        active_modes=(1, 2),
        mode_max_quanta=(0, 3, 2),
        excitation_class_limits={1: (1, 2), 2: (2, 3)},
        force_constant_threshold_cm=1.0,
    )

    comparison = compare_vpt2_vci(qff, max_quanta=3, n_roots=4, options=options)

    assert comparison.vpt2.basis == comparison.vci.basis[: len(comparison.vpt2.basis)]
    assert all(len(state) == 2 for state in comparison.vpt2.basis)
    assert all(state[0] <= 3 and state[1] <= 2 for state in comparison.vpt2.basis)
    assert comparison.energy_differences_cm.shape == (4,)
    assert np.max(np.abs(comparison.excitation_differences_cm)) < 20.0


def test_gf_from_gaussian_cartesian_hessian_uses_oracle_nonredundant_gics():
    canonical = hessian_input_from_gaussian_fchk(FCHK)
    result = gf_from_hessian_input_with_oracle_gics(canonical)
    adapter_result = gf_from_gaussian_fchk_with_oracle_gics(FCHK)

    assert result.b_matrix.shape == (3, 9)
    assert result.force_constants.shape == (3, 3)
    assert result.g_matrix.shape == (3, 3)
    assert result.ped.values.shape == (3, 3)
    assert result.gic_labels == ("GIC001", "GIC002", "GIC003")
    assert result.primitive_labels == ("R(1,2)", "R(2,3)", "A(1,2,3)")
    assert np.all(result.frequencies_cm > 0.0)
    assert np.allclose(result.ped.values.sum(axis=0), np.full(3, 100.0))
    assert np.allclose(result.frequencies_cm, [2169.878, 4141.256, 4392.363], atol=1.0e-3)
    assert np.allclose(adapter_result.frequencies_cm, result.frequencies_cm)


def test_gui_service_reports_are_independent_from_qt(tmp_path):
    qff_file = tmp_path / "field.qff"
    qff_file.write_text("FREQ 1 1000.0\nFREQ 2 1500.0\nQUARTIC 1 1 1 1 0.8\n", encoding="utf-8")

    qff = load_force_field(qff_path=qff_file)
    report = run_vpt2_vci_report(qff, max_quanta=2, roots=3, options=VCIOptions(active_modes=(0, 1)))

    assert "VPT2/VCI comparison on canonical ORACLE QFF" in report.text
    assert len(report.comparison.vci.basis) >= 3


def test_report_service_exposes_davidson_vci_path():
    qff = zero_anharmonic_force_field(np.array([100.0, 140.0]))

    report = run_vpt2_vci_report(qff, max_quanta=2, roots=3, vci_method="davidson")

    assert report.comparison.vci.davidson is not None
    assert "VPT2/VCI comparison on canonical ORACLE QFF" in report.text


def test_vpt2_vci_section_roundtrip_and_collect_outputs(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    run_dir = tmp_path / "vpt2_vci"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")
    _write_vpt2_vci_outputs(run_dir)
    section = vpt2_vci_section_from_run(
        source_kind="xyzin",
        source_path=xyzin,
        run_dir=run_dir,
        report_path=run_dir / "vpt2_vci.report",
        csv_dir=run_dir,
        manifest_path=run_dir / "vpt2_vci_manifest.json",
        max_quanta=3,
        roots=4,
        vci_method="davidson",
        outputs={},
        status="prepared",
    )

    write_vpt2_vci_section(xyzin, section)
    snapshot = refresh_vpt2_vci_section(xyzin)
    parsed = read_vpt2_vci_section(xyzin)

    assert snapshot.status == "complete"
    assert snapshot.comparison[1].delta_exc_cm == 1.5
    assert snapshot.frequencies[0].harmonic_frequency_cm == 100.0
    assert snapshot.mode_contributions[0].expected_quanta == 0.25
    assert parsed.status == "complete"
    assert parsed.outputs["comparison"] == run_dir / "vpt2_vci_comparison.csv"


def test_vpt2_vci_collect_marks_partial_when_csvs_are_missing(tmp_path):
    run_dir = tmp_path / "vpt2_vci"
    run_dir.mkdir()
    (run_dir / "vpt2_vci.report").write_text("report only\n", encoding="utf-8")
    section = vpt2_vci_section_from_run(
        source_kind="xyzin",
        source_path=tmp_path / "molecule.xyzin",
        run_dir=run_dir,
        report_path=run_dir / "vpt2_vci.report",
        csv_dir=run_dir,
        manifest_path=None,
        max_quanta=2,
        roots=2,
        vci_method="dense",
    )

    snapshot = collect_vpt2_vci_outputs(section)

    assert snapshot.status == "partial"
    assert snapshot.missing_outputs == ("comparison", "frequencies", "mode_contributions")


def test_force_field_validation_rejects_bad_modes():
    qff = QuarticForceField(
        harmonic_frequencies_cm=np.array([1000.0]),
        cubic_cm={(0, 0, 2): 1.0},
        quartic_cm={},
    )

    with pytest.raises(ScientificValidationError):
        validate_force_field(qff)


def test_vpt2_vci_fortran_inventory_uses_oracle_engines_layout():
    inventory = inventory_vpt2_vci_backends(ROOT)

    assert inventory.harmonic_source == ROOT / "engines" / "fortran" / "vpt2_vci" / "gf_core.f"
    assert {path.name for path in inventory.active_fortran_sources} == {
        "gf_core.f",
        "vpt2_core.f",
        "vci_core.f",
        "davidson_core.f",
    }
    assert inventory.davidson_backend == ROOT / "engines" / "fortran" / "vpt2_vci" / "davidson_core.f"


def _write_vpt2_vci_outputs(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "vpt2_vci.report").write_text(
        "VPT2/VCI comparison on canonical ORACLE QFF\n",
        encoding="utf-8",
    )
    (run_dir / "vpt2_vci_frequencies.csv").write_text(
        "mode,harmonic_frequency_cm-1\n1,100.0\n2,150.0\n",
        encoding="utf-8",
    )
    (run_dir / "vpt2_vci_comparison.csv").write_text(
        "root,vpt2_abs_cm-1,vci_abs_cm-1,delta_abs_cm-1,"
        "vpt2_exc_cm-1,vci_exc_cm-1,delta_exc_cm-1\n"
        "1,125.0,125.0,0.0,0.0,0.0,0.0\n"
        "2,225.0,226.5,1.5,100.0,101.5,1.5\n",
        encoding="utf-8",
    )
    (run_dir / "vpt2_vci_mode_contributions.csv").write_text(
        "root,mode,expected_quanta\n1,1,0.25\n1,2,0.75\n",
        encoding="utf-8",
    )
