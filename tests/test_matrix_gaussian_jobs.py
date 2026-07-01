from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from matrix_gaussian import (
    GaussianInputError,
    anharmonic_input_from_gaussian_fchk,
    ensure_gjf_input,
    gaussian_completion_message,
    gaussian_job_status,
    read_gaussian_fchk_qff,
    read_gaussian_fchk_geometry,
    read_indexed_qff_text,
    run_gaussian_job,
    select_latest_log,
)


ROOT = Path(__file__).resolve().parents[1]
FCHK = ROOT / "tests" / "fixtures" / "gf" / "h2o.fchk"


def test_ensure_gjf_input_prefers_existing_gjf(tmp_path):
    (tmp_path / "gauin").write_text("raw\n", encoding="utf-8")
    gjf = tmp_path / "gauin.gjf"
    gjf.write_text("gjf\n", encoding="utf-8")

    assert ensure_gjf_input(tmp_path) == gjf
    assert gjf.read_text(encoding="utf-8") == "gjf\n"


def test_ensure_gjf_input_copies_raw_gauin(tmp_path):
    (tmp_path / "gauin").write_text("raw\n", encoding="utf-8")

    assert ensure_gjf_input(tmp_path) == tmp_path / "gauin.gjf"
    assert (tmp_path / "gauin.gjf").read_text(encoding="utf-8") == "raw\n"


def test_ensure_gjf_input_requires_input(tmp_path):
    with pytest.raises(GaussianInputError):
        ensure_gjf_input(tmp_path)


def test_select_latest_log_uses_mtime_then_size(tmp_path):
    old_log = tmp_path / "gauin.log"
    new_log = tmp_path / "gauout.log"
    old_log.write_text("old\n", encoding="utf-8")
    new_log.write_text("newer\n", encoding="utf-8")
    os.utime(old_log, (1, 1))
    os.utime(new_log, (2, 2))

    assert select_latest_log(tmp_path) == new_log


def test_gaussian_status_and_completion_detect_normal_termination(tmp_path):
    (tmp_path / "gauin.log").write_text(
        "Normal termination of Gaussian 16\n",
        encoding="utf-8",
    )

    success, message = gaussian_completion_message(tmp_path, exit_code=1)
    status = gaussian_job_status(tmp_path)

    assert success
    assert message == "Gaussian completed successfully"
    assert status.status == "completed"
    assert status.normal_termination is True


def test_run_gaussian_job_uses_non_gui_backend(tmp_path):
    (tmp_path / "gauin").write_text("#p hf/sto-3g\n\njob\n\n0 1\nH 0 0 0\n\n", encoding="utf-8")
    executable = tmp_path / "fake_gaussian.sh"
    executable.write_text(
        "#!/bin/sh\n"
        'test -f "$1" || exit 9\n'
        "printf 'Normal termination of Gaussian 16\\n' > gauin.log\n"
        "exit 0\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)

    result = run_gaussian_job(tmp_path, executable=str(executable))

    assert result.success is True
    assert result.exit_code == 0
    assert result.input_path == tmp_path / "gauin.gjf"
    assert (tmp_path / "gauin.log").exists()


def test_matrix_gaussian_owns_fchk_qff_adapter():
    data = read_gaussian_fchk_qff(FCHK)
    anh = anharmonic_input_from_gaussian_fchk(FCHK)

    assert data.atomic_numbers.tolist() == [1, 8, 1]
    assert data.cartesian_hessian_lower.shape == (45,)
    assert np.allclose(anh.anharmonic_frequencies_cm[:3], [2123.50470, 4016.61987, 4266.73074])


def test_gaussian_fchk_geometry_reader_does_not_require_hessian_consumers():
    geometry = read_gaussian_fchk_geometry(FCHK)

    assert geometry.source_format == "gaussian_fchk"
    assert geometry.atoms == ("H", "O", "H")
    assert geometry.coordinates_angstrom.shape == (3, 3)


def test_indexed_qff_text_reader_lives_in_matrix_gaussian(tmp_path):
    qff_file = tmp_path / "field.qff"
    qff_file.write_text(
        "FREQ 1 100.0\nCUBIC 1 1 1 -2.0\nQUARTIC 1 1 1 1 4.0\n",
        encoding="utf-8",
    )

    qff = read_indexed_qff_text(qff_file)

    assert np.allclose(qff.harmonic_frequencies_cm, [100.0])
    assert qff.cubic_cm == {(0, 0, 0): -2.0}
    assert qff.quartic_cm == {(0, 0, 0, 0): 4.0}
