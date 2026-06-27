from __future__ import annotations

from oracle_core import section_content
from oracle_rovib import (
    ORACLE_XYZ_ROTATIONAL_SCHEMA,
    ORACLE_XYZ_VIBRATIONAL_SCHEMA,
    RotationalSection,
    WMSRotInputOptions,
    WMSRotSimulationOptions,
    VibrationalSection,
    parse_rotational_section,
    parse_vibrational_section,
    read_rotational_section,
    read_vibrational_section,
    rovib_summary_lines,
    simulate_wmsrot_spectrum,
    summarize_xyzin,
    wmsrot_input_text_from_xyzin,
    write_rotational_section,
    write_vibrational_section,
)
from oracle_thermo import (
    ORACLE_XYZ_THERMO_SCHEMA,
    ThermoContribution,
    ThermoSection,
    parse_thermo_section,
    read_thermo_section,
    run_thermo_on_xyzin,
    thermo_section_lines,
    write_thermo_section,
)


def test_rovib_and_thermo_schema_constants_are_oracle_xyz_sections():
    assert ORACLE_XYZ_ROTATIONAL_SCHEMA == "oracle.xyz.rotational.v1"
    assert ORACLE_XYZ_VIBRATIONAL_SCHEMA == "oracle.xyz.vibrational.v1"
    assert ORACLE_XYZ_THERMO_SCHEMA == "oracle.xyz.thermo.v1"


def test_rotational_section_accepts_merlino_keys():
    section = parse_rotational_section(
        [
            "rotor_type = asymmetric_top_quasi_prolate",
            "Point Group = C2v",
            "Watson Reduction = A",
            "Symm. Number = 2",
            "A_MHz = 1000.0",
            "B_MHz = 800.0",
            "C_MHz = 600.0",
            "Dipole_a_D = 1.1",
            "Dipole_b_D = 2.2",
            "Dipole_c_D = 3.3",
            "DVibA_MHz=1.0",
            "DVibB_MHz=2.0",
            "DVibC_MHz=3.0",
            "Q_rot = 10.5",
        ]
    )

    assert section.point_group == "C2v"
    assert section.watson_reduction == "A"
    assert section.symmetry_number == 2
    assert section.A_MHz == 1000.0
    assert section.dipole_debye == (1.1, 2.2, 3.3)
    assert section.delta_vib_MHz == (1.0, 2.0, 3.0)
    assert section.q_rot == 10.5


def test_rotational_section_writer_preserves_other_sections(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text("1\nh\nH 0 0 0\n\n#GIC\nSCHEMA oracle.xyz.gic.v1\n", encoding="utf-8")

    write_rotational_section(
        path,
        RotationalSection(
            rotor_type="linear_top",
            B_MHz=42.0,
            symmetry_number=1,
            temperature_K=298.15,
            pressure_atm=1.0,
        ),
    )
    parsed = read_rotational_section(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    assert parsed.B_MHz == 42.0
    assert section_content(lines, "GIC")[0] == "SCHEMA oracle.xyz.gic.v1"


def test_wmsrot_input_export_uses_normalized_rotational_sections(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#BASIC",
                "POINT_GROUP = C2v",
                "WATSON_REDUCTION = A",
                "T_K = 100.0",
                "",
                "#ROTATIONAL",
                "rotor_type = asymmetric_top_quasi_prolate",
                "representation = Ir",
                "A_MHz = 1000.0",
                "B_MHz = 800.0",
                "C_MHz = 600.0",
                "Dipole_a_D = 1.0",
                "Dipole_b_D = 0.0",
                "Dipole_c_D = 2.0",
                "",
                "#DELTABVIB",
                "DVibA_MHz = 1.0",
                "DVibB_MHz = 2.0",
                "DVibC_MHz = 3.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    text = wmsrot_input_text_from_xyzin(
        path,
        options=WMSRotInputOptions(j_min=1, j_max=8),
    )

    assert "rotor_type = asymmetric" in text
    assert "Watson Reduction = A" in text
    assert "Point Group = C2v" in text
    assert "T_K = 100" in text
    assert "DVibC_MHz= 3" in text
    assert "Dipole_c_D = 2" in text
    assert "J_MIN = 1" in text
    assert "J_MAX = 8" in text
    assert "DELTA J_MHz =  0" in text


def test_wmsrot_simulation_wrapper_calls_vendored_engine_with_deltabvib():
    class FakeEngine:
        REDUCTION = "S"
        REPRESENTATION = "Ir"
        FREQ_UNIT = "MHz"

        def __init__(self):
            self.calls = []

        def simulate_rigid_spectrum(self, *args):
            self.calls.append(args)
            return [{"Frequency (MHz)": 1.0, "Relative intensity": 2.0}]

    engine = FakeEngine()

    rows = simulate_wmsrot_spectrum(
        RotationalSection(
            rotor_type="asymmetric",
            representation="IIIl",
            point_group="C2v",
            watson_reduction="A",
            temperature_K=100.0,
            A_MHz=1000.0,
            B_MHz=800.0,
            C_MHz=600.0,
            dipole_debye=(1.0, 0.0, 2.0),
            delta_vib_MHz=(1.0, 2.0, 3.0),
        ),
        options=WMSRotSimulationOptions(j_min=1, j_max=4),
        engine=engine,
    )

    assert rows == [{"Frequency (MHz)": 1.0, "Relative intensity": 2.0}]
    call = engine.calls[0]
    assert call[:4] == (100.0, 1001.0, 802.0, 603.0)
    assert call[16:19] == (1.0, 0.0, 2.0)
    assert call[19:23] == (4, 1.0e-20, "C2v", "asymmetric")
    assert call[-1] == 1


def test_vibrational_section_reads_frequencies_and_chi_block():
    section = parse_vibrational_section(
        [
            "linear = 0",
            "nvib = 3",
            "n_imag_like = 1",
            "symmetry_group = Cs",
            "freq_cm1 = 100.0 200.0 300.0",
            "ir_inten_km_mol = 1.0 2.0 3.0",
            "chi_cm1 = [",
            "1 1 -0.5",
            "2 1 0.1",
            "]",
        ]
    )

    assert section.linear is False
    assert section.nvib == 3
    assert section.frequencies_cm1 == (100.0, 200.0, 300.0)
    assert section.ir_intensities_km_mol == (1.0, 2.0, 3.0)
    assert section.chi_cm1 == ((1, 1, -0.5), (2, 1, 0.1))


def test_vibrational_section_writer_round_trips(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text("1\nh\nH 0 0 0\n", encoding="utf-8")

    write_vibrational_section(
        path,
        VibrationalSection(
            linear=False,
            nvib=2,
            frequencies_cm1=(100.0, 200.0),
            chi_cm1=((1, 1, -0.5),),
        ),
    )
    parsed = read_vibrational_section(path)

    assert parsed.nvib == 2
    assert parsed.frequencies_cm1 == (100.0, 200.0)
    assert parsed.chi_cm1 == ((1, 1, -0.5),)


def test_rovib_summary_reads_standalone_xyzin(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#BASIC",
                "CHARGE 0",
                "SPIN_MULTIPLICITY 1",
                "POINT_GROUP C1",
                "",
                "#ROTATIONAL",
                "A_MHz = 1000.0",
                "B_MHz = 900.0",
                "C_MHz = 800.0",
                "",
                "#VIBRATIONAL",
                "freq_cm1 = 100.0 200.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = summarize_xyzin(path)
    lines = rovib_summary_lines(summary)

    assert summary.basic.point_group == "C1"
    assert summary.rotational is not None
    assert summary.vibrational is not None
    assert "rotational: A=1000MHz B=900MHz C=800MHz" in lines
    assert "vibrational: 2 frequencies" in lines


def test_thermo_section_round_trips_merlino_labels(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text("1\nh\nH 0 0 0\n", encoding="utf-8")

    write_thermo_section(
        path,
        ThermoSection(
            translational=ThermoContribution(Q_dimless=2.0, U_kJmol=1.0),
            rotational=ThermoContribution(Q_dimless=3.0, S_JmolK=4.0),
            vibrational=ThermoContribution(
                Q_dimless=1.0,
                available=False,
                reason="no frequencies above cutoff",
            ),
            total=ThermoContribution(Q_dimless=6.0, U_kJmol=1.0, S_JmolK=4.0),
        ),
    )
    parsed = read_thermo_section(path)
    lines = thermo_section_lines(parsed)

    assert parsed.translational is not None
    assert parsed.translational.Q_dimless == 2.0
    assert parsed.rotational is not None
    assert parsed.rotational.S_JmolK == 4.0
    assert parsed.vibrational is not None
    assert parsed.vibrational.available is False
    assert parsed.total is not None
    assert parsed.total.Q_dimless == 6.0
    assert "Q_dimless_trasl = 2" in lines


def test_thermo_section_accepts_legacy_merlino_comment_lines():
    parsed = parse_thermo_section(
        [
            "# keys: Q_dimless U_kJmol H_kJmol S_JmolK Cv_JmolK Cp_JmolK",
            "Q_dimless_trasl = 2.0",
            "U_kJmol_tot = 5.0",
        ]
    )

    assert parsed.translational is not None
    assert parsed.translational.Q_dimless == 2.0
    assert parsed.total is not None
    assert parsed.total.U_kJmol == 5.0


def test_section_content_keeps_thermo_comment_until_real_next_section(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#THERMO",
                "# keys: legacy Merlino comment",
                "Q_dimless_trasl = 2.0",
                "",
                "#GIC",
                "SCHEMA oracle.xyz.gic.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    thermo = read_thermo_section(path)
    content = section_content(path.read_text(encoding="utf-8").splitlines(), "THERMO")

    assert "# keys: legacy Merlino comment" in content
    assert thermo.translational is not None
    assert thermo.translational.Q_dimless == 2.0


def test_thermo_pipeline_writes_thermo_section_and_report(tmp_path):
    path = tmp_path / "water.xyzin"
    path.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.758602 0.000000 0.504284",
                "H -0.758602 0.000000 0.504284",
                "",
                "#BASIC",
                "SCHEMA oracle.xyz.basic.v1",
                "T_K = 298.15",
                "P_ATM = 1.0",
                "",
                "#ROTATIONAL",
                "A_MHz = 835000.0",
                "B_MHz = 435000.0",
                "C_MHz = 287000.0",
                "Symm. Number = 2",
                "rotor_type = asymmetric_prolate",
                "",
                "#VIBRATIONAL",
                "freq_cm1 = 1595.0 3657.0 3756.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = run_thermo_on_xyzin(path, report=True)
    parsed = read_thermo_section(path)
    report = tmp_path / "thermo.report"

    assert result.total is not None
    assert result.total.Q_dimless is not None
    assert result.total.Q_dimless > 0.0
    assert parsed.total is not None
    assert parsed.total.H_kJmol is not None
    assert report.exists()
    assert "THERMO PIPELINE REPORT" in report.read_text(encoding="utf-8")
