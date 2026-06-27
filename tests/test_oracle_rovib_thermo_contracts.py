from __future__ import annotations

import numpy as np

from oracle_core import section_content
from oracle_qm import normal_modes_section_from_arrays, write_normal_modes_section
from oracle_rovib import (
    ORACLE_XYZ_ROTATIONAL_SCHEMA,
    ORACLE_XYZ_VIBRATIONAL_SCHEMA,
    RotationalSection,
    WMSRotInputOptions,
    WMSRotSimulationOptions,
    VibrationalSection,
    VibrationalSpectrumOptions,
    build_hybrid_vibrational_spectrum_from_xyzin,
    build_vibrational_spectrum,
    compare_vibrational_spectra,
    fetch_nist_ir_gas_phase_csv,
    nist_ir_points_to_spectrum,
    parse_nist_jcamp_ir_points,
    parse_rotational_section,
    parse_vibrational_section,
    read_rotational_section,
    read_vibrational_section,
    rovib_summary_lines,
    simulate_wmsrot_spectrum,
    summarize_xyzin,
    wmsrot_input_text_from_xyzin,
    write_vibrational_spectrum_comparison_outputs,
    write_normal_mode_match_csv,
    write_rotational_section,
    write_vibrational_spectrum_outputs,
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


def test_wmsrot_vendor_engine_routes_eigh_through_oracle_diagonalizer(monkeypatch):
    import numpy as np

    from oracle_rovib.vendor import wmsrot_engine

    calls = []

    def fake_eigh(matrix):
        calls.append(matrix.shape)
        return np.linalg.eigh(matrix)

    monkeypatch.setattr(wmsrot_engine, "_oracle_eigh_arrays", fake_eigh)

    values, vectors = wmsrot_engine._oracle_eigh(np.diag([2.0, 1.0]))

    assert calls == [(2, 2)]
    assert np.allclose(values, [1.0, 2.0])
    assert vectors.shape == (2, 2)


def test_vibrational_section_reads_frequencies_and_chi_block():
    section = parse_vibrational_section(
        [
            "linear = 0",
            "nvib = 3",
            "n_imag_like = 1",
            "symmetry_group = Cs",
            "freq_cm1 = 100.0 200.0 300.0",
            "anharmonic_freq_cm1 = 98.0 197.0 296.0",
            "ir_inten_km_mol = 1.0 2.0 3.0",
            "raman_act_a4_amu = 4.0 5.0 6.0",
            "vcd_rot_strength = -0.1 0.0 0.2",
            "roa_inten = 0.01 -0.02 0.03",
            "chi_cm1 = [",
            "1 1 -0.5",
            "2 1 0.1",
            "]",
        ]
    )

    assert section.linear is False
    assert section.nvib == 3
    assert section.frequencies_cm1 == (100.0, 200.0, 300.0)
    assert section.anharmonic_frequencies_cm1 == (98.0, 197.0, 296.0)
    assert section.ir_intensities_km_mol == (1.0, 2.0, 3.0)
    assert section.raman_activities_A4_amu == (4.0, 5.0, 6.0)
    assert section.vcd_rot_strengths == (-0.1, 0.0, 0.2)
    assert section.roa_intensities == (0.01, -0.02, 0.03)
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
            anharmonic_frequencies_cm1=(98.0, 197.0),
            ir_intensities_km_mol=(1.0, 2.0),
            raman_activities_A4_amu=(3.0, 4.0),
            chi_cm1=((1, 1, -0.5),),
        ),
    )
    parsed = read_vibrational_section(path)

    assert parsed.nvib == 2
    assert parsed.frequencies_cm1 == (100.0, 200.0)
    assert parsed.anharmonic_frequencies_cm1 == (98.0, 197.0)
    assert parsed.raman_activities_A4_amu == (3.0, 4.0)
    assert parsed.chi_cm1 == ((1, 1, -0.5),)


def test_vibrational_spectrum_service_draws_and_exports_ir_raman_vcd_roa(tmp_path):
    section = VibrationalSection(
        frequencies_cm1=(100.0, 200.0),
        anharmonic_frequencies_cm1=(98.0, 197.0),
        ir_intensities_km_mol=(1.0, 2.0),
        raman_activities_A4_amu=(3.0, 4.0),
        vcd_rot_strengths=(-1.0, 2.0),
        roa_intensities=(0.5, -0.25),
    )

    ir = build_vibrational_spectrum(
        section,
        observable="IR",
        options=VibrationalSpectrumOptions(fwhm_cm1=8.0, step_cm1=2.0),
    )
    vcd = build_vibrational_spectrum(
        section,
        observable="VCD",
        options=VibrationalSpectrumOptions(fwhm_cm1=8.0, step_cm1=2.0),
    )
    assert ir.observable == "IR"
    assert len(ir.peaks) == 2
    assert max(ir.y) == 1.0
    assert min(vcd.y) < 0.0

    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text("1\nh\nH 0 0 0\n", encoding="utf-8")
    write_vibrational_section(xyzin, section)
    spectrum = write_vibrational_spectrum_outputs(
        xyzin,
        csv_path=tmp_path / "ir.csv",
        plot_path=tmp_path / "ir.svg",
        peaks_path=tmp_path / "ir_peaks.csv",
        observable="RAMAN",
        source="harmonic",
        options=VibrationalSpectrumOptions(fwhm_cm1=8.0, step_cm1=2.0),
    )

    assert spectrum.observable == "RAMAN"
    assert (tmp_path / "ir.csv").is_file()
    assert (tmp_path / "ir.svg").is_file()
    assert (tmp_path / "ir_peaks.csv").is_file()


def test_vibrational_spectrum_comparison_mirrors_ir_but_not_vcd(tmp_path):
    section = VibrationalSection(
        frequencies_cm1=(100.0, 200.0),
        anharmonic_frequencies_cm1=(101.0, 198.0),
        ir_intensities_km_mol=(1.0, 2.0),
        anharmonic_ir_intensities_km_mol=(0.5, 1.5),
        vcd_rot_strengths=(-1.0, 2.0),
        anharmonic_vcd_rot_strengths=(-0.8, 1.8),
    )
    options = VibrationalSpectrumOptions(fwhm_cm1=8.0, step_cm1=2.0)
    harmonic_ir = build_vibrational_spectrum(section, observable="IR", options=options)
    anharmonic_ir = build_vibrational_spectrum(
        section,
        observable="IR",
        source="anharmonic",
        options=options,
    )
    harmonic_vcd = build_vibrational_spectrum(section, observable="VCD", options=options)
    anharmonic_vcd = build_vibrational_spectrum(
        section,
        observable="VCD",
        source="anharmonic",
        options=options,
    )

    ir_comparison = compare_vibrational_spectra(harmonic_ir, anharmonic_ir)
    vcd_comparison = compare_vibrational_spectra(harmonic_vcd, anharmonic_vcd)

    assert ir_comparison.mirror_second
    assert np.allclose(ir_comparison.plotted_second_y, -ir_comparison.second_y)
    assert not vcd_comparison.mirror_second
    assert np.allclose(vcd_comparison.plotted_second_y, vcd_comparison.second_y)

    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text("1\nh\nH 0 0 0\n", encoding="utf-8")
    write_vibrational_section(xyzin, section)
    written = write_vibrational_spectrum_comparison_outputs(
        xyzin,
        csv_path=tmp_path / "compare.csv",
        plot_path=tmp_path / "compare.svg",
        observable="IR",
        options=options,
    )

    assert written.mirror_second
    assert (tmp_path / "compare.csv").is_file()
    assert (tmp_path / "compare.svg").is_file()


def test_vibrational_spectrum_comparison_accepts_two_xyzin_files(tmp_path):
    first_section = VibrationalSection(
        frequencies_cm1=(100.0, 200.0),
        ir_intensities_km_mol=(1.0, 2.0),
    )
    second_section = VibrationalSection(
        frequencies_cm1=(110.0, 210.0),
        ir_intensities_km_mol=(0.8, 1.8),
    )
    first_xyzin = tmp_path / "b3lyp.xyzin"
    second_xyzin = tmp_path / "ccsd_t.xyzin"
    first_xyzin.write_text("1\nfirst\nH 0 0 0\n", encoding="utf-8")
    second_xyzin.write_text("1\nsecond\nH 0 0 0\n", encoding="utf-8")
    write_vibrational_section(first_xyzin, first_section)
    write_vibrational_section(second_xyzin, second_section)

    comparison = write_vibrational_spectrum_comparison_outputs(
        first_xyzin,
        second_xyzin=second_xyzin,
        csv_path=tmp_path / "two_files.csv",
        observable="IR",
        first_source="harmonic",
        second_source="harmonic",
        options=VibrationalSpectrumOptions(fwhm_cm1=8.0, step_cm1=2.0),
    )

    assert comparison.mirror_second
    assert comparison.first.peaks[0].frequency_cm1 == 100.0
    assert comparison.second.peaks[0].frequency_cm1 == 110.0
    assert (tmp_path / "two_files.csv").is_file()


def test_hybrid_anharmonic_spectrum_matches_modes_before_applying_correction(tmp_path):
    level1 = tmp_path / "level1.xyzin"
    level2 = tmp_path / "level2.xyzin"
    level1.write_text("1\nlevel1\nH 0 0 0\n", encoding="utf-8")
    level2.write_text("1\nlevel2\nH 0 0 0\n", encoding="utf-8")
    write_vibrational_section(
        level1,
        VibrationalSection(
            frequencies_cm1=(100.0, 200.0, 300.0),
            ir_intensities_km_mol=(1.0, 2.0, 3.0),
        ),
    )
    write_vibrational_section(
        level2,
        VibrationalSection(
            frequencies_cm1=(195.0, 305.0, 95.0),
            anharmonic_frequencies_cm1=(190.0, 315.0, 93.0),
            ir_intensities_km_mol=(9.0, 9.0, 9.0),
        ),
    )
    write_normal_modes_section(
        level1,
        normal_modes_section_from_arrays(
            (100.0, 200.0, 300.0),
            np.eye(3),
            coordinate_count=3,
        ),
    )
    write_normal_modes_section(
        level2,
        normal_modes_section_from_arrays(
            (195.0, 305.0, 95.0),
            np.asarray(
                [
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, -1.0],
                    [-1.0, 0.0, 0.0],
                ]
            ),
            coordinate_count=3,
        ),
    )

    result = build_hybrid_vibrational_spectrum_from_xyzin(
        level1,
        level2,
        observable="IR",
        options=VibrationalSpectrumOptions(fwhm_cm1=8.0, step_cm1=2.0),
        min_mode_overlap=0.95,
    )
    match_csv = write_normal_mode_match_csv(tmp_path / "matches.csv", result.matches)

    assert [match.level2_mode for match in result.matches] == [3, 1, 2]
    assert np.allclose([peak.frequency_cm1 for peak in result.spectrum.peaks], [98.0, 195.0, 310.0])
    assert [peak.source for peak in result.spectrum.peaks] == ["hybrid", "hybrid", "hybrid"]
    assert np.allclose([peak.intensity for peak in result.spectrum.peaks], [1.0, 2.0, 3.0])
    assert match_csv.is_file()


def test_nist_jcamp_parser_and_gas_phase_download_policy(tmp_path, monkeypatch):
    from oracle_rovib import vibspec

    jcamp_gas = "\n".join(
        [
            "##TITLE=METHANE",
            "##STATE=GAS (150 mmHg)",
            "##XFACTOR=1",
            "##YFACTOR=1",
            "##DELTAX=1",
            "##XYDATA=(X++(Y..Y))",
            "100.0 0.9 0.8",
            "102.0 0.7",
            "##END=",
        ]
    )
    jcamp_liquid = jcamp_gas.replace("##STATE=GAS (150 mmHg)", "##STATE=LIQUID")
    page = '<a href="/cgi/cbook.cgi?JCAMP=C74828&Index=1&Type=IR">JCAMP</a>'

    points = parse_nist_jcamp_ir_points(jcamp_gas)
    assert [point.wavenumber_cm1 for point in points] == [100.0, 101.0, 102.0]
    assert [point.value for point in points] == [0.9, 0.8, 0.7]
    experimental = nist_ir_points_to_spectrum(points)
    assert experimental.observable == "IR"
    assert experimental.source == "nist-gas-experiment"
    assert max(experimental.y) == 1.0

    def fake_fetch(url, *, timeout, encoding="utf-8"):
        if "JCAMP" in url:
            return jcamp_gas
        return page

    monkeypatch.setattr(vibspec, "_fetch_text", fake_fetch)
    result = fetch_nist_ir_gas_phase_csv("74-82-8", tmp_path / "nist.csv")

    assert result.status == "downloaded"
    assert result.csv_path == tmp_path / "nist.csv"
    assert "wavenumber_cm-1,transmittance" in (tmp_path / "nist.csv").read_text(encoding="utf-8")

    def fake_fetch_liquid(url, *, timeout, encoding="utf-8"):
        if "JCAMP" in url:
            return jcamp_liquid
        return page

    monkeypatch.setattr(vibspec, "_fetch_text", fake_fetch_liquid)
    liquid = fetch_nist_ir_gas_phase_csv("74-82-8", tmp_path / "liquid.csv")

    assert liquid.status == "not_gas_phase"
    assert liquid.needs_user_instruction
    assert not (tmp_path / "liquid.csv").exists()


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
