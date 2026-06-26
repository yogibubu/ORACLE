from __future__ import annotations

import numpy as np

from oracle_rovib import (
    VibinData,
    compute_coriolis_from_vibin,
    compute_qcent,
    direct_vibrational_dos_from_xyzin,
    qcent_report_lines,
    read_vibin,
    rovib_pipeline,
    write_vibin,
)


def _sample_vibin_data() -> VibinData:
    return VibinData(
        symbols=("O", "H", "H"),
        coords_A=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.75, 0.0, 0.5],
                [-0.75, 0.0, 0.5],
            ],
            dtype=float,
        ),
        masses_amu=np.array([15.9949, 1.0078, 1.0078], dtype=float),
        freq_cm1=np.array([1000.0, 1500.0, 3000.0], dtype=float),
        modes_mw=np.array(
            [
                [[0.01, 0.02, 0.00], [0.03, 0.00, 0.01], [0.00, 0.04, 0.02]],
                [[0.02, 0.00, 0.01], [0.00, 0.05, 0.02], [0.01, 0.00, 0.03]],
                [[0.00, 0.03, 0.02], [0.02, 0.01, 0.00], [0.03, 0.02, 0.01]],
            ],
            dtype=float,
        ),
        representation="Ir",
        linear=False,
        didq_sym6=np.ones((6, 3), dtype=float) * 0.01,
    )


def test_vibin_round_trip_and_coriolis(tmp_path):
    path = tmp_path / "vibin"
    written = write_vibin(path, _sample_vibin_data())
    parsed = read_vibin(written)
    coriolis = compute_coriolis_from_vibin(parsed, 1000.0, 800.0, 600.0, Geff_thr_cm1=0.0)

    assert parsed.natoms == 3
    assert parsed.nvib == 3
    assert parsed.didq_sym6 is not None
    assert parsed.didq_sym6.shape == (6, 3)
    assert len(coriolis.entries) > 0


def test_qcent_report_from_vibin_data():
    result = compute_qcent(_sample_vibin_data())
    lines = qcent_report_lines(result)

    assert result.nvib == 3
    assert len(result.QA_MHz) == 5
    assert "A-reduction_MHz DelJ DelJK DelK delJ delK" in lines


def test_vibrational_and_rovibrational_dos_from_xyzin(tmp_path):
    xyzin = tmp_path / "mol.xyzin"
    xyzin.write_text(
        "\n".join(
            [
                "2",
                "linear",
                "H 0 0 0",
                "H 0 0 0.74",
                "",
                "#BASIC",
                "T_K = 298.15",
                "",
                "#ROTATIONAL",
                "rotor_type = linear",
                "A_MHz = 0.0",
                "B_MHz = 60000.0",
                "C_MHz = 60000.0",
                "Symm. Number = 2",
                "",
                "#VIBRATIONAL",
                "freq_cm1 = 1000.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    vib_dos = direct_vibrational_dos_from_xyzin(
        xyzin,
        vmax=2,
        emax_cm1=2500.0,
        bin_cm1=100.0,
        out=tmp_path / "dos_vib.dat",
    )
    rovib = rovib_pipeline(
        xyzin,
        vib_dos=vib_dos.path,
        out=tmp_path / "dos_rovib.dat",
        rot_out=tmp_path / "dos_rot.dat",
        q_out=tmp_path / "rovib_qt.dat",
        jmax=3,
    )

    assert vib_dos.path is not None
    assert vib_dos.path.exists()
    assert len(vib_dos.bins_logg) > 0
    assert rovib.dos_rovib.exists()
    assert rovib.q_path is not None
    assert rovib.q_path.exists()
    assert rovib.Q_rovib > 0.0
