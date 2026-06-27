from __future__ import annotations

from pathlib import Path

import numpy as np

from oracle_gaussian import compute_deltavib_from_alpha, parse_gaussian_rovib_log, promote_gaussian_rovib_to_xyzin
from oracle_rovib import (
    DeltaBVibAlphaRow,
    read_deltabvib_section,
    read_rotational_section,
    read_vibrational_section,
)


def _write_xyzin(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "2",
                "h2",
                "H 0 0 0",
                "H 0 0 0.74",
                "",
                "#BASIC",
                "T_K = 298.15",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _gaussian_rovib_log() -> str:
    return (
        "Full point group                 C2V     NOp   4\n"
        " T =  298.15 K; Pressure =  1.00000 Atm.\n"
        "Harmonic frequencies (cm**-1)\n"
        " -------------------\n"
        " Frequencies --  -50.0000  1000.0000\n"
        " IR Inten    --    1.5000     2.5000\n"
        "Total Anharmonic X Matrix (in cm^-1)\n"
        "      1       2\n"
        "  1  -10.0000   0.5000\n"
        "  2    0.5000 -20.0000\n"
        "\n"
        "Vibro-Rot alpha Matrix (in MHz)\n"
        "       Q(  1)      2.0000      4.0000      6.0000\n"
        "       Q(  2)     10.0000     20.0000     30.0000\n"
        "\n"
        "Rotational Constants (in MHz)\n"
        "Ae=  100.000     A00=  111.000     A0=  110.000\n"
        "Be=  200.000     B00=  222.000     B0=  220.000\n"
        "Ce=  300.000     C00=  333.000     C0=  330.000\n"
        "Dipole moment (Debye):\n"
        "      1.1000000      2.2000000      3.3000000  Tot=      4.1158231\n"
        "Normal termination of Gaussian 16\n"
    )


def test_parse_gaussian_rovib_log_builds_oracle_sections(tmp_path):
    log = tmp_path / "job.log"
    log.write_text(_gaussian_rovib_log(), encoding="utf-8")

    data = parse_gaussian_rovib_log(log)

    assert data.vibrational is not None
    assert data.vibrational.frequencies_cm1 == (-50.0, 1000.0)
    assert data.vibrational.ir_intensities_km_mol == (1.5, 2.5)
    assert data.vibrational.chi_cm1 == ((1, 1, -10.0), (2, 1, 0.5), (2, 2, -20.0))
    assert data.deltabvib is not None
    assert data.deltabvib.source == "gaussian-rotational-constants"
    assert data.deltabvib.delta_A_MHz == 11.0
    assert data.deltabvib.delta_B_MHz == 22.0
    assert data.deltabvib.delta_C_MHz == 33.0
    assert len(data.deltabvib.alpha_rows_MHz) == 2
    assert data.rotational is not None
    assert data.rotational.point_group == "C2V"
    assert data.rotational.temperature_K == 298.15
    assert data.rotational.A_MHz == 110.0
    assert data.rotational.dipole_debye == (1.1, 2.2, 3.3)
    assert data.rotational.delta_vib_MHz == (11.0, 22.0, 33.0)


def test_promote_gaussian_rovib_to_xyzin_writes_shared_sections(tmp_path):
    log = tmp_path / "job.log"
    xyzin = tmp_path / "molecule.xyzin"
    log.write_text(_gaussian_rovib_log(), encoding="utf-8")
    _write_xyzin(xyzin)

    result = promote_gaussian_rovib_to_xyzin(log, xyzin)
    vib = read_vibrational_section(xyzin)
    delta = read_deltabvib_section(xyzin)
    rot = read_rotational_section(xyzin)

    assert result.wrote_vibrational is True
    assert result.wrote_deltabvib is True
    assert result.wrote_rotational is True
    assert vib.nvib == 2
    assert delta.alpha_rows_MHz[0].mode == 1
    assert rot.A_MHz == 110.0
    assert rot.dipole_debye == (1.1, 2.2, 3.3)
    assert rot.delta_vib_MHz == (11.0, 22.0, 33.0)


def test_alpha_delta_fallback_respects_imaginary_and_excluded_modes():
    rows = (
        DeltaBVibAlphaRow(1, 2.0, 4.0, 6.0),
        DeltaBVibAlphaRow(2, 10.0, 20.0, 30.0),
        DeltaBVibAlphaRow(3, 100.0, 200.0, 300.0),
    )

    delta = compute_deltavib_from_alpha(
        rows,
        (-50.0, 1000.0, 1500.0),
        exclude_modes=(3,),
    )

    assert delta is not None
    assert np.allclose(delta, ((-2.0 + 10.0) / 2.0, (-4.0 + 20.0) / 2.0, (-6.0 + 30.0) / 2.0))

    assert compute_deltavib_from_alpha(rows, exclude_modes=(1, 2, 3)) == (0.0, 0.0, 0.0)


def test_anharmonic_x_matrix_handles_gaussian_column_blocks(tmp_path):
    log = tmp_path / "job.log"
    log.write_text(
        "\n".join(
            [
                "Harmonic frequencies (cm**-1)",
                " Frequencies --  100.0000  200.0000  300.0000",
                "Total Anharmonic X Matrix (in cm^-1)",
                "      1       2",
                "  1   -1.0000",
                "  2    0.2000  -2.0000",
                "  3    0.3000   0.4000",
                "      3",
                "  3   -3.0000",
                "Fundamental Bands",
            ]
        ),
        encoding="utf-8",
    )

    data = parse_gaussian_rovib_log(log)

    assert data.vibrational is not None
    assert data.vibrational.chi_cm1 == (
        (1, 1, -1.0),
        (2, 1, 0.2),
        (2, 2, -2.0),
        (3, 1, 0.3),
        (3, 2, 0.4),
        (3, 3, -3.0),
    )
