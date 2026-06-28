from __future__ import annotations

from pathlib import Path

import numpy as np

from matrix_chem import preprocess_to_enriched_xyz, write_validation_section
from matrix_gaussian import (
    hessian_input_from_gaussian_fchk,
    promote_gaussian_electronic_log_to_xyzin,
    promote_gaussian_fchk_to_xyzin,
)
from matrix_gf import run_xyzin_gf_report_from_fchk, run_xyzin_gf_report_from_xyzin
from matrix_neo import write_gicforge_build_sections
from matrix_qm import (
    QFFSection,
    ElectronicSection,
    ElectronicStateRecord,
    ElectronicTransitionRecord,
    OrbitalFileRecord,
    OrbitalsSection,
    PropertiesSection,
    PropertyRecord,
    TransitionsSection,
    electronic_section_lines,
    merge_properties_section,
    hessian_input_from_xyzin,
    orbital_file_record_from_path,
    parse_electronic_section,
    parse_orbitals_section,
    parse_properties_section,
    parse_transitions_section,
    properties_section_lines,
    property_records_by_name,
    property_records_for_atom,
    qff_section_from_quartic_force_field,
    quartic_force_field_from_qff_section,
    read_cartesian_hessian_section,
    read_electronic_section,
    read_normal_modes_section,
    read_orbitals_section,
    read_properties_section,
    read_qff_section,
    read_transitions_section,
    write_properties_section,
    write_qff_section,
)
from matrix_vpt2_vci import QuarticForceField, load_force_field
from tools import matrix_run


ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ROOT / "tests" / "fixtures" / "test_molecules" / "molecules"
FCHK = ROOT / "tests" / "fixtures" / "gf" / "h2o.fchk"


def _prepared_h2o_xyzin(tmp_path: Path) -> Path:
    xyzin = tmp_path / "h2o.xyzin"
    preprocess_to_enriched_xyz(MOLECULES / "h2ocart.inp", xyzin)
    write_validation_section(xyzin)
    write_gicforge_build_sections(xyzin)
    return xyzin


def test_gaussian_fchk_promotion_writes_qm_xyzin_sections(tmp_path):
    xyzin = _prepared_h2o_xyzin(tmp_path)

    result = promote_gaussian_fchk_to_xyzin(FCHK, xyzin)
    hessian_section = read_cartesian_hessian_section(xyzin)
    hessian_input = hessian_input_from_xyzin(xyzin)
    qff = read_qff_section(xyzin)
    electronic = read_electronic_section(xyzin)
    orbitals = read_orbitals_section(xyzin)

    assert result.wrote_cartesian_hessian is True
    assert result.wrote_normal_modes is True
    assert result.wrote_qff is True
    assert result.wrote_electronic is True
    assert result.wrote_orbitals is True
    assert hessian_section.source == "gaussian-fchk"
    assert hessian_input.source == "gaussian-fchk"
    assert np.allclose(
        hessian_input.cartesian_hessian,
        hessian_input_from_gaussian_fchk(FCHK).cartesian_hessian,
    )
    assert qff.source == "gaussian-fchk"
    assert np.allclose(qff.anharmonic_frequencies_cm[:3], [2123.50470, 4016.61987, 4266.73074])
    assert read_normal_modes_section(xyzin).modes.shape[1] == 9
    assert electronic.states[0].label == "S0"
    assert np.isclose(electronic.states[0].energy_hartree, -74.96590119079237)
    assert any(record.format == "FCHK" and record.role == "orbitals" for record in orbitals.files)
    assert any(record.format == "FCHK" and record.role == "density" for record in orbitals.files)


def test_electronic_xyzin_sections_roundtrip(tmp_path):
    molden = tmp_path / "mol.molden"
    cube = tmp_path / "density.cube"
    state_lines = electronic_section_lines(
        ElectronicSection(
            (
                ElectronicStateRecord(
                    "S0",
                    energy_hartree=-75.0,
                    energy_ev=0.0,
                    multiplicity="1",
                    symmetry="A1",
                    source="unit",
                ),
            )
        )
    )
    transition_lines = [
        "SCHEMA oracle.xyz.transitions.v1",
        "COLUMNS FROM TO ENERGY_EV WAVELENGTH_NM OSC STRENGTH SOURCE",
        "S0 S1 4.2 295.2 0.1 electric-dipole unit",
    ]
    orbital_lines = [
        "SCHEMA oracle.xyz.orbitals.v1",
        "COLUMNS KIND FORMAT ROLE PATH LABEL SOURCE",
        f"FILE MOLDEN orbitals {molden} mol unit",
        f"FILE CUBE density {cube} density unit",
    ]

    electronic = parse_electronic_section(state_lines)
    transitions = parse_transitions_section(transition_lines)
    orbitals = parse_orbitals_section(orbital_lines)

    assert isinstance(electronic, ElectronicSection)
    assert isinstance(transitions, TransitionsSection)
    assert isinstance(orbitals, OrbitalsSection)
    assert electronic.states[0].symmetry == "A1"
    assert isinstance(transitions.transitions[0], ElectronicTransitionRecord)
    assert transitions.transitions[0].oscillator_strength == 0.1
    assert isinstance(orbitals.files[0], OrbitalFileRecord)
    assert orbitals.files[0] == orbital_file_record_from_path(molden, source="unit")
    assert orbitals.files[1].role == "density"


def test_properties_section_roundtrip_preserves_conversion_metadata(tmp_path):
    xyzin = tmp_path / "properties.xyzin"
    xyzin.write_text("2\nnqcc\nN 0 0 0\nH 0 0 1\n", encoding="utf-8")
    record = PropertyRecord(
        name="NUCLEAR_QUADRUPOLE_COUPLING",
        target="atom",
        target_id="N1",
        atom=1,
        isotope="14N",
        value=(-4.123456, 2.061728, 2.061728),
        unit="MHz",
        axes="PAS:chi_aa,chi_bb,chi_cc",
        program="Molpro",
        method="CCSD(T)",
        level="cc-pVTZ",
        source="molpro.out",
        status="converted",
        conversion="EFG_AU_TO_MHZ_WITH_14N_Q",
        uncertainty=0.002,
        comment="converted from Molpro EFG",
    )

    lines = properties_section_lines(PropertiesSection((record,)))
    parsed = parse_properties_section(lines)
    write_properties_section(xyzin, parsed)
    restored = read_properties_section(xyzin)

    assert restored.schema == "oracle.xyz.properties.v1"
    assert restored.records == (record,)
    assert restored.records[0].target == "ATOM"
    assert property_records_by_name(restored, "nuclear_quadrupole_coupling") == (record,)
    assert property_records_for_atom(restored, 1) == (record,)


def test_properties_merge_and_cli_summary(tmp_path, capsys):
    xyzin = tmp_path / "properties.xyzin"
    xyzin.write_text("1\nnqcc\nN 0 0 0\n", encoding="utf-8")
    raw = PropertyRecord(
        name="EFG_TENSOR",
        target="ATOM",
        atom=1,
        isotope="14N",
        value=(-1.0, 0.5, 0.5),
        unit="a.u.",
        axes="PAS",
        program="Molpro",
        method="HF",
        level="cc-pVDZ",
        source="molpro.out",
    )
    converted = PropertyRecord(
        name="EFG_TENSOR",
        target="ATOM",
        atom=1,
        isotope="14N",
        value=(-2.0, 1.0, 1.0),
        unit="a.u.",
        axes="PAS",
        program="Molpro",
        method="HF",
        level="cc-pVDZ",
        source="molpro-rerun.out",
    )

    merge_properties_section(xyzin, (raw,))
    merge_properties_section(xyzin, (converted,))
    section = read_properties_section(xyzin)

    assert isinstance(section, PropertiesSection)
    assert len(section.records) == 1
    assert section.records[0].value == (-2.0, 1.0, 1.0)

    rc = matrix_run.main(["properties", "summary", str(xyzin), "--atom", "1"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "properties: 1" in out
    assert "EFG_TENSOR" in out


def test_gaussian_electronic_log_promotion_writes_sections(tmp_path):
    xyzin = tmp_path / "mol.xyzin"
    log = tmp_path / "td.log"
    cube = tmp_path / "density.cube"
    xyzin.write_text("1\nmol\nH 0 0 0\n", encoding="utf-8")
    cube.write_text("cube\n", encoding="utf-8")
    log.write_text(
        "\n".join(
            [
                " SCF Done:  E(RHF) =  -7.500000000000E+01     A.U. after 1 cycles",
                " Excited State   1:      Singlet-A1     4.2000 eV  295.20 nm  f=0.1000",
                " Excited State   2:      Triplet-B2     5.3000 eV  233.93 nm  f=0.0000",
            ]
        ),
        encoding="utf-8",
    )

    result = promote_gaussian_electronic_log_to_xyzin(log, xyzin, orbital_files=(cube,))
    electronic = read_electronic_section(xyzin)
    transitions = read_transitions_section(xyzin)
    orbitals = read_orbitals_section(xyzin)

    assert result.wrote_electronic
    assert result.wrote_transitions
    assert result.wrote_orbitals
    assert [state.label for state in electronic.states] == ["S0", "S1", "S2"]
    assert transitions.transitions[0].to_state == "S1"
    assert transitions.transitions[0].wavelength_nm == 295.2
    assert orbitals.files[0].format == "CUBE"


def test_gf_can_run_from_frozen_xyzin_hessian_section(tmp_path):
    xyzin = _prepared_h2o_xyzin(tmp_path)
    promote_gaussian_fchk_to_xyzin(FCHK, xyzin)

    from_xyzin = run_xyzin_gf_report_from_xyzin(xyzin)
    from_fchk = run_xyzin_gf_report_from_fchk(FCHK, xyzin)

    assert "#CARTESIAN_HESSIAN" in from_xyzin.text
    assert np.allclose(from_xyzin.result.frequencies_cm, from_fchk.result.frequencies_cm)
    assert from_xyzin.result.gic_labels == from_fchk.result.gic_labels


def test_qff_section_roundtrip_preserves_force_terms(tmp_path):
    xyzin = tmp_path / "field.xyzin"
    xyzin.write_text("1\nqff\nH 0 0 0\n", encoding="utf-8")
    source = QuarticForceField(
        harmonic_frequencies_cm=np.array([100.0, 200.0]),
        cubic_cm={(0, 0, 1): -2.0},
        quartic_cm={(0, 1, 1, 1): 0.5},
    )

    write_qff_section(xyzin, qff_section_from_quartic_force_field(source, source="unit-test"))
    section = read_qff_section(xyzin)
    restored = quartic_force_field_from_qff_section(section)

    assert isinstance(section, QFFSection)
    assert section.source == "unit-test"
    assert restored.cubic_cm == source.cubic_cm
    assert restored.quartic_cm == source.quartic_cm
    assert np.allclose(restored.harmonic_frequencies_cm, source.harmonic_frequencies_cm)
    assert np.allclose(load_force_field(xyzin_path=xyzin).harmonic_frequencies_cm, [100.0, 200.0])
