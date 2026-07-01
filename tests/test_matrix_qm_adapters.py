from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from matrix_chem import (
    detect_qm_output_format,
    preprocess_to_enriched_xyz,
    read_geometry,
    write_validation_section,
)
from matrix_core import read_basic_section, read_sectioned_lines, section_content
from matrix_gaussian import (
    parse_gaussian_quadrupole_properties,
    promote_gaussian_quadrupole_properties_to_xyzin,
)
from matrix_molpro import (
    parse_molpro_quadrupole_properties,
    promote_molpro_molden_to_xyzin,
    promote_molpro_quadrupole_properties_to_xyzin,
    read_molpro_output_geometry,
    summarize_molpro_output,
)
from matrix_mrcc import read_mrcc_output_geometry, summarize_mrcc_output
from matrix_orca import (
    convert_orca_gbw_to_molden,
    hessian_input_from_orca_output,
    parse_orca_quadrupole_properties,
    promote_orca_molden_to_xyzin,
    promote_orca_output_to_xyzin,
    promote_orca_quadrupole_properties_to_xyzin,
    read_orca_output_geometry,
    summarize_orca_output,
)
from matrix_gf import read_gf_ped_section, run_xyzin_gf_report_from_xyzin, write_gf_ped_section_from_report
from matrix_neo import write_gicforge_build_sections
from matrix_qm import (
    EFG_AU_TO_NQCC_MHZ_PER_BARN,
    read_cartesian_hessian_section,
    read_orbitals_section,
    read_properties_section,
)
from tools import matrix_run


ROOT = Path(__file__).resolve().parents[1]
NH3_QUADRUPOLE_FIXTURE = (
    ROOT / "tests" / "fixtures" / "qm_properties" / "nh3_hf_ccpvtz_quadrupole"
)


def _single_nqcc_record(records):
    matches = [record for record in records if record.name == "NUCLEAR_QUADRUPOLE_COUPLING"]
    assert len(matches) == 1
    return matches[0]


def _molpro_output() -> str:
    return "\n".join(
        [
            " PROGRAM SYSTEM MOLPRO",
            " Charge = -1",
            " Spin quantum number = 0.5",
            " ATOMIC COORDINATES",
            " NR  ATOM    CHARGE       X              Y              Z",
            "  1  O       8.000000    0.00000000     0.00000000     0.00000000",
            "  2  H       1.000000    0.00000000     0.00000000     0.96000000",
            "  3  H       1.000000    0.00000000     0.76000000    -0.58000000",
            "",
        ]
    )


def _molpro_quadrupole_output() -> str:
    return "\n".join(
        [
            " PROGRAM SYSTEM MOLPRO",
            " SETTING BASIS          =    CC-PVTZ",
            " ATOMIC COORDINATES",
            " NR  ATOM    CHARGE       X              Y              Z",
            "  1  N       7.000000   -1.37109100     0.00000000     0.00000000",
            "  2  H       1.000000   -1.74414300    -0.47145500     0.81658400",
            "  3  H       1.000000   -1.74414300    -0.47145500    -0.81658400",
            "  4  H       1.000000   -1.74414300     0.94291100     0.00000000",
            "  5  F       9.000000    1.27143200     0.00000000     0.00000000",
            "",
            " !CCSD(T) (ED,norela <1.1|FGXX|1.1>    -9.091670249869",
            " !CCSD(T) (ED,norela <1.1|FGYY|1.1>     4.545835124971",
            " !CCSD(T) (ED,norela <1.1|FGZZ|1.1>     4.545835124898",
            "",
        ]
    )


def _gaussian_quadrupole_output() -> str:
    return "\n".join(
        [
            "#p b3lyp/6-31g nmr",
            "",
            " Nuclear quadrupole coupling constants (MHz)",
            " Atom Isotope ChiXX ChiYY ChiZZ ChiXY ChiXZ ChiYZ",
            " 1 14N -4.100000 2.050000 2.050000 0.010000 0.020000 0.030000",
            "",
            " Normal termination of Gaussian",
        ]
    )


def _gaussian_pickett_quadrupole_output() -> str:
    return "\n".join(
        [
            "#p hf/cc-pVTZ scf=tight prop=efg output=pickett nosymm",
            "",
            "Nuclear quadrupole coupling constants [Chi] (MHz):",
            " 1  N(14) ",
            "  aa=     2.3869   ba=    -0.0000   ca=    -0.0000",
            "  ab=    -0.0000   bb=     2.3828   cb=     0.0008",
            "  ac=    -0.0000   bc=     0.0008   cc=    -4.7697",
            "Dipole moment (Debye):",
            "",
            "Normal termination of Gaussian DV",
        ]
    )


def _gaussian_efg_quadrupole_output() -> str:
    return "\n".join(
        [
            "#p hf/cc-pVTZ scf=tight prop=efg nosymm",
            "",
            "Input orientation:",
            " ---------------------------------------------------------------------",
            " Center     Atomic      Atomic             Coordinates (Angstroms)",
            " Number     Number       Type             X           Y           Z",
            " ---------------------------------------------------------------------",
            "      1          7           0        0.000000    0.000000    0.000000",
            " ---------------------------------------------------------------------",
            " -----------------------------------------------------",
            "    Center         ---- Electric Field Gradient ----",
            "                     XY            XZ            YZ",
            " -----------------------------------------------------",
            "    1 Atom        0.000000      0.000000     -0.000174",
            " -----------------------------------------------------",
            " -----------------------------------------------------",
            "    Center         ---- Electric Field Gradient ----",
            "                       ( tensor representation )",
            "                   3XX-RR        3YY-RR        3ZZ-RR",
            " -----------------------------------------------------",
            "    1 Atom       -0.496965     -0.496105      0.993070",
            " -----------------------------------------------------",
            "",
            "Normal termination of Gaussian",
        ]
    )


def _orca_quadrupole_output(*, direct: bool) -> str:
    rows = (
        [
            " NUCLEAR QUADRUPOLE COUPLING CONSTANTS (MHz)",
            " Atom Isotope ChiXX ChiYY ChiZZ ChiXY ChiXZ ChiYZ",
            " 1 14N -4.100000 2.050000 2.050000 0.010000 0.020000 0.030000",
        ]
        if direct
        else [
            " ELECTRIC FIELD GRADIENTS (a.u.)",
            " Atom Isotope VXX VYY VZZ VXY VXZ VYZ",
            " 1 14N -9.091670249869 4.545835124971 4.545835124898 0.0 0.0 0.0",
        ]
    )
    return "\n".join(["Program ORCA", *rows, "****ORCA TERMINATED NORMALLY****"])


def _orca_eprnmr_quadrupole_output() -> str:
    return "\n".join(
        [
            "Program ORCA",
            "Nucleus   0N : A  : Isotope=   14 I=  1.0 P= 38.5677 MHz/au**3",
            "               Q  : Isotope=   14 I=  1.0 Q=  0.0204 barn",
            "Raw EFG matrix (all values in a.u.**-3):",
            " ------------------------------------------------------",
            "             0.4970207      -0.0000000       0.0000000",
            "            -0.0000000       0.4961606       0.0001744",
            "             0.0000000       0.0001744      -0.9931813",
            "",
            "V(Tot)       0.4961606       0.4970207      -0.9931813",
            "",
            "Quadrupole tensor eigenvalues (in MHz;Q= 0.0204 I=  1.0)",
            " e**2qQ            =    -4.769937 MHz",
            " eta               =     0.000866",
            "****ORCA TERMINATED NORMALLY****",
        ]
    )


def _mrcc_output() -> str:
    return "\n".join(
        [
            " MRCC program system",
            " Charge of the system: 1",
            " Spin multiplicity: 2",
            " Cartesian coordinates",
            " ---------------------",
            " O       0.00000000     0.00000000     0.00000000",
            " H       0.00000000     0.00000000     0.96000000",
            " H       0.00000000     0.76000000    -0.58000000",
            "",
        ]
    )


def _orca_output() -> str:
    hessian = np.eye(9) * 0.01
    hessian_rows = ["          " + " ".join(str(idx) for idx in range(9))]
    for row in range(9):
        values = " ".join(f"{hessian[row, col]: .8f}" for col in range(9))
        hessian_rows.append(f"{row:5d} {values}")
    return "\n".join(
        [
            "Program ORCA",
            "* xyz 0 1",
            "Total Charge           Charge    ....    0",
            "Multiplicity           Mult      ....    1",
            "FINAL SINGLE POINT ENERGY     -75.000000000000",
            "CARTESIAN COORDINATES (ANGSTROEM)",
            "---------------------------------",
            "O       0.00000000     0.00000000     0.00000000",
            "H       0.00000000     0.00000000     0.96000000",
            "H       0.00000000     0.76000000    -0.58000000",
            "",
            "VIBRATIONAL FREQUENCIES",
            "-----------------------",
            "   0:          0.00 cm**-1",
            "   1:          0.00 cm**-1",
            "   2:          0.00 cm**-1",
            "   3:          0.00 cm**-1",
            "   4:          0.00 cm**-1",
            "   5:          0.00 cm**-1",
            "   6:       1600.00 cm**-1",
            "   7:       3650.00 cm**-1",
            "   8:       3750.00 cm**-1",
            "",
            "CARTESIAN HESSIAN",
            "-----------------",
            *hessian_rows,
            "",
            "****ORCA TERMINATED NORMALLY****",
            "",
        ]
    )


def test_molpro_output_adapter_returns_shared_geometry(tmp_path):
    output = tmp_path / "molpro.out"
    output.write_text(_molpro_output(), encoding="utf-8")

    summary = summarize_molpro_output(output)
    geometry = read_molpro_output_geometry(output)

    assert summary.atomic_coordinate_blocks == 1
    assert geometry.source_format == "molpro_output"
    assert geometry.atoms == ("O", "H", "H")
    assert geometry.charge == -1
    assert geometry.multiplicity == 2
    assert np.allclose(geometry.coordinates_angstrom[1], [0.0, 0.0, 0.96])


def test_molpro_quadrupole_adapter_converts_efg_to_properties(tmp_path):
    output = tmp_path / "molpro.out"
    xyzin = tmp_path / "mol.xyzin"
    output.write_text(_molpro_quadrupole_output(), encoding="utf-8")
    xyzin.write_text("1\nn\nN 0 0 0\n", encoding="utf-8")

    records = parse_molpro_quadrupole_properties(output)
    result = promote_molpro_quadrupole_properties_to_xyzin(output, xyzin)
    section = read_properties_section(xyzin)

    nqcc = [record for record in records if record.name == "NUCLEAR_QUADRUPOLE_COUPLING"][0]
    expected = -EFG_AU_TO_NQCC_MHZ_PER_BARN * 0.02044 * -9.091670249869
    assert result.wrote_properties is True
    assert result.property_count == 3
    assert len(section.records) == len(records)
    assert [record.name for record in section.records] == [record.name for record in records]
    assert nqcc.atom == 1
    assert nqcc.isotope == "14N"
    assert nqcc.status == "converted"
    assert np.isclose(nqcc.value[0], expected)
    stored_nqcc = [
        record for record in section.records if record.name == "NUCLEAR_QUADRUPOLE_COUPLING"
    ][0]
    assert np.allclose(stored_nqcc.value, nqcc.value)


def test_molpro_quadrupole_adapter_uses_input_geometry_fallback(tmp_path):
    output = tmp_path / "molpro.out"
    output.write_text(
        "\n".join(
            [
                " geometry={",
                " 6",
                " N    -1.371091    0.000000   -0.000000",
                " H    -1.744143   -0.471455    0.816584",
                " F     1.271432    0.000000    0.000000",
                " }",
                " ATOMIC COORDINATES:q",
                "   2     0.00092229     0.46181898     8.14e-05      2  2    0.42",
                "   3     0.00092684     0.46294355     4.81e-06      3  3    0.66",
                " !CCSD(T) (ED,norela <1.1|FGXX|1.1>    -9.091670249869",
                " !CCSD(T) (ED,norela <1.1|FGYY|1.1>     4.545835124971",
                " !CCSD(T) (ED,norela <1.1|FGZZ|1.1>     4.545835124898",
            ]
        ),
        encoding="utf-8",
    )

    records = parse_molpro_quadrupole_properties(output)

    assert records[0].atom == 1
    assert records[0].isotope == "14N"


def test_gaussian_quadrupole_adapter_promotes_direct_constants(tmp_path):
    output = tmp_path / "gaussian.log"
    xyzin = tmp_path / "mol.xyzin"
    output.write_text(_gaussian_quadrupole_output(), encoding="utf-8")
    xyzin.write_text("1\nn\nN 0 0 0\n", encoding="utf-8")

    records = parse_gaussian_quadrupole_properties(output)
    result = promote_gaussian_quadrupole_properties_to_xyzin(output, xyzin)
    section = read_properties_section(xyzin)

    assert result.wrote_properties is True
    assert result.property_count == 1
    assert section.records == records
    assert records[0].program == "Gaussian"
    assert records[0].method == "b3lyp/6-31g"
    assert records[0].value == (-4.1, 2.05, 2.05, 0.01, 0.02, 0.03)


def test_gaussian_quadrupole_adapter_prefers_pickett_constants(tmp_path):
    output = tmp_path / "gaussian_pickett.log"
    output.write_text(_gaussian_pickett_quadrupole_output(), encoding="utf-8")

    records = parse_gaussian_quadrupole_properties(output)

    assert len(records) == 1
    assert records[0].program == "Gaussian"
    assert records[0].axes == "PICKETT:chi_aa,chi_bb,chi_cc,chi_ab,chi_ac,chi_bc"
    assert records[0].isotope == "14N"
    assert np.allclose(records[0].value, (2.3869, 2.3828, -4.7697, -0.0, -0.0, 0.0008))


def test_gaussian_quadrupole_adapter_converts_standard_efg_to_pickett(tmp_path):
    output = tmp_path / "gaussian_efg.log"
    output.write_text(_gaussian_efg_quadrupole_output(), encoding="utf-8")

    records = parse_gaussian_quadrupole_properties(output)
    converted = [record for record in records if record.name == "NUCLEAR_QUADRUPOLE_COUPLING"]

    assert len(records) == 3
    assert len(converted) == 1
    assert np.isclose(converted[0].value[2], -EFG_AU_TO_NQCC_MHZ_PER_BARN * 0.02044 * 0.993070)
    assert "convention=Pickett/Gaussian-EFG" in converted[0].conversion


def test_orca_quadrupole_adapter_promotes_direct_or_converted_constants(tmp_path):
    direct = tmp_path / "orca_direct.out"
    efg = tmp_path / "orca_efg.out"
    direct.write_text(_orca_quadrupole_output(direct=True), encoding="utf-8")
    efg.write_text(_orca_quadrupole_output(direct=False), encoding="utf-8")

    direct_records = parse_orca_quadrupole_properties(direct)
    efg_records = parse_orca_quadrupole_properties(efg)

    assert len(direct_records) == 1
    assert direct_records[0].status == "raw"
    assert direct_records[0].value[:3] == (-4.1, 2.05, 2.05)
    converted = [record for record in efg_records if record.name == "NUCLEAR_QUADRUPOLE_COUPLING"]
    assert len(converted) == 1
    assert converted[0].status == "converted"
    assert np.isclose(converted[0].value[0], EFG_AU_TO_NQCC_MHZ_PER_BARN * 0.02044 * -9.091670249869)


def test_orca_quadrupole_adapter_reads_eprnmr_vtot_block(tmp_path):
    output = tmp_path / "orca_eprnmr.out"
    output.write_text(_orca_eprnmr_quadrupole_output(), encoding="utf-8")

    records = parse_orca_quadrupole_properties(output)
    converted = [record for record in records if record.name == "NUCLEAR_QUADRUPOLE_COUPLING"]

    assert len(records) == 3
    assert len(converted) == 1
    assert converted[0].axes == "ORCA_EFG_PAS:Vxx,Vyy,Vzz"
    assert np.isclose(converted[0].value[2], -4.769937, atol=5.0e-5)
    assert "convention=Pickett/ORCA-EFG" in converted[0].conversion


def test_nh3_hf_ccpvtz_quadrupole_golden_matches_across_qm_codes():
    gaussian = _single_nqcc_record(
        parse_gaussian_quadrupole_properties(NH3_QUADRUPOLE_FIXTURE / "gaussian_pickett.log")
    )
    molpro = _single_nqcc_record(
        parse_molpro_quadrupole_properties(NH3_QUADRUPOLE_FIXTURE / "molpro.out")
    )
    orca = _single_nqcc_record(
        parse_orca_quadrupole_properties(NH3_QUADRUPOLE_FIXTURE / "orca.out")
    )

    assert np.allclose(gaussian.value[:3], (2.3869, 2.3828, -4.7697), atol=5.0e-5)
    assert np.allclose(molpro.value[:3], (2.3868979, 2.3827673, -4.7696652), atol=5.0e-5)
    assert np.allclose(orca.value[:3], (2.3828998, 2.3870306, -4.7699304), atol=5.0e-5)
    assert np.allclose(gaussian.value[:3], molpro.value[:3], atol=5.0e-5)
    assert np.allclose(
        (*sorted(gaussian.value[:2]), gaussian.value[2]),
        (*sorted(orca.value[:2]), orca.value[2]),
        atol=3.0e-4,
    )


def test_mrcc_output_adapter_returns_shared_geometry(tmp_path):
    output = tmp_path / "mrcc.out"
    output.write_text(_mrcc_output(), encoding="utf-8")

    summary = summarize_mrcc_output(output)
    geometry = read_mrcc_output_geometry(output)

    assert summary.cartesian_coordinate_blocks == 1
    assert geometry.source_format == "mrcc_output"
    assert geometry.atoms == ("O", "H", "H")
    assert geometry.charge == 1
    assert geometry.multiplicity == 2
    assert np.allclose(geometry.coordinates_angstrom[2], [0.0, 0.76, -0.58])


def test_orca_output_adapter_returns_shared_geometry_and_hessian(tmp_path):
    output = tmp_path / "orca.out"
    output.write_text(_orca_output(), encoding="utf-8")

    summary = summarize_orca_output(output)
    geometry = read_orca_output_geometry(output)
    hessian = hessian_input_from_orca_output(output)

    assert summary.cartesian_coordinate_blocks == 1
    assert summary.normal_termination is True
    assert summary.final_energy_hartree == -75.0
    assert geometry.source_format == "orca_output"
    assert geometry.atoms == ("O", "H", "H")
    assert geometry.charge == 0
    assert geometry.multiplicity == 1
    assert np.allclose(geometry.coordinates_angstrom[1], [0.0, 0.0, 0.96])
    assert np.allclose(summary.frequencies_cm[-3:], [1600.0, 3650.0, 3750.0])
    assert hessian.source == "orca-output"
    assert hessian.cartesian_hessian.shape == (9, 9)


def test_read_geometry_auto_dispatches_qm_outputs(tmp_path):
    molpro = tmp_path / "molpro.out"
    mrcc = tmp_path / "mrcc.out"
    orca = tmp_path / "orca.out"
    molpro.write_text(_molpro_output(), encoding="utf-8")
    mrcc.write_text(_mrcc_output(), encoding="utf-8")
    orca.write_text(_orca_output(), encoding="utf-8")

    assert detect_qm_output_format(molpro) == "molpro"
    assert detect_qm_output_format(mrcc) == "mrcc"
    assert detect_qm_output_format(orca) == "orca"
    assert read_geometry(molpro).source_format == "molpro_output"
    assert read_geometry(mrcc).source_format == "mrcc_output"
    assert read_geometry(orca).source_format == "orca_output"


def test_link_preprocess_imports_molpro_with_basic_section(tmp_path):
    output = tmp_path / "molpro.out"
    xyzin = tmp_path / "molecule.xyzin"
    output.write_text(_molpro_output(), encoding="utf-8")

    result = preprocess_to_enriched_xyz(output, xyzin, source_kind="molpro")
    basic = read_basic_section(xyzin)
    lines = read_sectioned_lines(xyzin)

    assert result.geometry.source_format == "molpro_output"
    assert basic.charge == -1
    assert basic.multiplicity == 2
    assert section_content(lines, "SOURCE")[2] == "FORMAT molpro_output"
    assert section_content(lines, "SYMMETRY")[0] == "SCHEMA oracle.xyz.symmetry.v1"
    assert section_content(lines, "TOPOLOGY")[0] == "SCHEMA matrix.xyz.topology.v1"


def test_orca_promote_writes_geometry_and_cartesian_hessian_sections(tmp_path):
    output = tmp_path / "orca.out"
    xyzin = tmp_path / "molecule.xyzin"
    output.write_text(_orca_output(), encoding="utf-8")

    result = promote_orca_output_to_xyzin(output, xyzin)
    basic = read_basic_section(xyzin)
    hessian = read_cartesian_hessian_section(xyzin)
    lines = read_sectioned_lines(xyzin)

    assert result.wrote_geometry is True
    assert result.wrote_cartesian_hessian is True
    assert basic.charge == 0
    assert basic.multiplicity == 1
    assert section_content(lines, "SOURCE")[2] == "FORMAT orca_output"
    assert hessian.source == "orca-output"
    assert hessian.cartesian_hessian.shape == (9, 9)
    assert np.allclose(hessian.harmonic_frequencies_cm[-3:], [1600.0, 3650.0, 3750.0])


def test_molpro_promote_cli_calls_adapter(tmp_path, monkeypatch, capsys):
    calls = {}
    output = tmp_path / "molpro.out"
    xyzin = tmp_path / "molecule.xyzin"

    def fake_promote(source, target, *, symmetry_distance, symmetry_inertia, max_rotation_order):
        calls["source"] = source
        calls["target"] = target
        calls["symmetry_distance"] = symmetry_distance
        calls["symmetry_inertia"] = symmetry_inertia
        calls["max_rotation_order"] = max_rotation_order
        return SimpleNamespace(
            path=target,
            geometry=SimpleNamespace(natoms=3),
            point_group="C1",
            topology_bond_count=2,
            ring_count=0,
        )

    monkeypatch.setattr("matrix_molpro.promote_molpro_output_to_xyzin", fake_promote)

    rc = matrix_run.main(
        [
            "molpro",
            "promote",
            str(output),
            str(xyzin),
            "--symmetry-distance",
            "0.002",
            "--symmetry-inertia",
            "0.003",
            "--max-rotation-order",
            "8",
        ]
    )

    assert rc == 0
    assert calls == {
        "source": output,
        "target": xyzin,
        "symmetry_distance": 0.002,
        "symmetry_inertia": 0.003,
        "max_rotation_order": 8,
    }
    assert "Promoted Molpro output" in capsys.readouterr().out


def test_molpro_molden_registers_orbitals_without_properties(tmp_path):
    output = tmp_path / "molpro.out"
    molden = tmp_path / "molpro.molden"
    xyzin = tmp_path / "molecule.xyzin"
    output.write_text(_molpro_output(), encoding="utf-8")
    molden.write_text("[Molden Format]\n", encoding="utf-8")
    xyzin.write_text("1\nn\nN 0 0 0\n", encoding="utf-8")

    result = promote_molpro_molden_to_xyzin(output, xyzin, molden=molden)
    orbitals = read_orbitals_section(xyzin)

    assert result.wrote_orbitals
    assert orbitals.files[0].format == "MOLDEN"
    assert orbitals.files[0].source == "molpro-molden"
    assert read_properties_section(xyzin).records == ()


def test_molpro_molden_cli_registers_orbitals(tmp_path, capsys):
    output = tmp_path / "molpro.out"
    molden = tmp_path / "molpro.molden.input"
    xyzin = tmp_path / "molecule.xyzin"
    output.write_text(_molpro_output(), encoding="utf-8")
    molden.write_text("[Molden Format]\n", encoding="utf-8")
    xyzin.write_text("1\nn\nN 0 0 0\n", encoding="utf-8")

    assert matrix_run.main(["molpro", "molden", str(output), str(xyzin), "--molden", str(molden)]) == 0
    out = capsys.readouterr().out

    assert "Registered Molpro Molden file" in out
    assert read_orbitals_section(xyzin).files[0].path == molden


def test_mrcc_summary_cli_calls_adapter(tmp_path, monkeypatch, capsys):
    calls = {}
    output = tmp_path / "mrcc.out"

    def fake_summary(source):
        calls["source"] = source
        return SimpleNamespace(
            path=source,
            geometry=SimpleNamespace(natoms=3),
            charge=1,
            multiplicity=2,
            cartesian_coordinate_blocks=1,
        )

    monkeypatch.setattr("matrix_mrcc.summarize_mrcc_output", fake_summary)

    rc = matrix_run.main(["mrcc", "summary", str(output)])
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {"source": output}
    assert "multiplicity: 2" in out


def test_orca_summary_and_promote_cli_use_adapter(tmp_path, capsys):
    output = tmp_path / "orca.out"
    xyzin = tmp_path / "molecule.xyzin"
    output.write_text(_orca_output(), encoding="utf-8")

    assert matrix_run.main(["orca", "summary", str(output)]) == 0
    summary_out = capsys.readouterr().out
    assert "final_energy_hartree: -75" in summary_out
    assert "cartesian_hessian: 1" in summary_out

    assert matrix_run.main(["orca", "promote", str(output), str(xyzin)]) == 0
    promote_out = capsys.readouterr().out
    assert "Promoted ORCA output" in promote_out
    assert "wrote_cartesian_hessian: 1" in promote_out
    assert read_cartesian_hessian_section(xyzin).cartesian_hessian.shape == (9, 9)


def test_orca_molden_converter_runs_orca_2mkl_and_registers_orbitals(tmp_path, monkeypatch):
    gbw = tmp_path / "orca.gbw"
    xyzin = tmp_path / "molecule.xyzin"
    gbw.write_text("gbw placeholder", encoding="utf-8")
    xyzin.write_text("1\nn\nN 0 0 0\n", encoding="utf-8")
    calls = {}

    def fake_run(command, *, check, text, capture_output, timeout):
        calls["command"] = command
        calls["check"] = check
        calls["text"] = text
        calls["capture_output"] = capture_output
        calls["timeout"] = timeout
        gbw.with_suffix(".molden.input").write_text("[Molden Format]\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("matrix_orca.parsers.subprocess.run", fake_run)

    conversion = convert_orca_gbw_to_molden(gbw, executable="orca_2mkl", timeout=12.0)
    result = promote_orca_molden_to_xyzin(gbw, xyzin, executable="orca_2mkl", timeout=12.0)
    orbitals = read_orbitals_section(xyzin)

    assert calls["command"] == ("orca_2mkl", str(gbw.with_suffix("")), "-molden")
    assert calls["timeout"] == 12.0
    assert conversion.molden_path == gbw.with_suffix(".molden.input")
    assert result.wrote_orbitals
    assert orbitals.files[0].format == "MOLDEN"
    assert orbitals.files[0].source == "orca_2mkl"
    assert read_properties_section(xyzin).records == ()


def test_orca_molden_cli_registers_orbitals(tmp_path, monkeypatch, capsys):
    gbw = tmp_path / "orca.gbw"
    xyzin = tmp_path / "molecule.xyzin"
    explicit_molden = tmp_path / "orbitals.molden"
    gbw.write_text("gbw placeholder", encoding="utf-8")
    xyzin.write_text("1\nn\nN 0 0 0\n", encoding="utf-8")

    def fake_run(command, *, check, text, capture_output, timeout):
        gbw.with_suffix(".molden.input").write_text("[Molden Format]\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("matrix_orca.parsers.subprocess.run", fake_run)

    assert (
        matrix_run.main(
            [
                "orca",
                "molden",
                str(gbw),
                str(xyzin),
                "--output",
                str(explicit_molden),
                "--timeout",
                "10",
            ]
        )
        == 0
    )
    out = capsys.readouterr().out

    assert "Converted ORCA GBW to Molden" in out
    assert explicit_molden.is_file()
    assert read_orbitals_section(xyzin).files[0].path == explicit_molden


def test_quadrupole_promote_cli_updates_properties(tmp_path, capsys):
    gaussian = tmp_path / "gaussian.log"
    molpro = tmp_path / "molpro.out"
    orca = tmp_path / "orca.out"
    gaussian_xyzin = tmp_path / "gaussian.xyzin"
    molpro_xyzin = tmp_path / "molpro.xyzin"
    orca_xyzin = tmp_path / "orca.xyzin"
    gaussian.write_text(_gaussian_quadrupole_output(), encoding="utf-8")
    molpro.write_text(_molpro_quadrupole_output(), encoding="utf-8")
    orca.write_text(_orca_quadrupole_output(direct=False), encoding="utf-8")
    for xyzin in (gaussian_xyzin, molpro_xyzin, orca_xyzin):
        xyzin.write_text("1\nn\nN 0 0 0\n", encoding="utf-8")

    assert matrix_run.main(["gaussian", "promote-quadrupole", str(gaussian), str(gaussian_xyzin)]) == 0
    assert matrix_run.main(["molpro", "promote-quadrupole", str(molpro), str(molpro_xyzin)]) == 0
    assert matrix_run.main(["orca", "promote-quadrupole", str(orca), str(orca_xyzin)]) == 0
    out = capsys.readouterr().out

    assert "Promoted Gaussian quadrupole properties" in out
    assert "Promoted Molpro quadrupole properties" in out
    assert "Promoted ORCA quadrupole properties" in out
    assert len(read_properties_section(gaussian_xyzin).records) == 1
    assert len(read_properties_section(molpro_xyzin).records) == 3
    assert len(read_properties_section(orca_xyzin).records) == 3


def test_orca_promoted_hessian_feeds_gf_without_orca_parsing_in_gf(tmp_path):
    output = tmp_path / "orca.out"
    xyzin = tmp_path / "molecule.xyzin"
    output.write_text(_orca_output(), encoding="utf-8")

    promote_orca_output_to_xyzin(output, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    report = run_xyzin_gf_report_from_xyzin(xyzin)
    write_gf_ped_section_from_report(
        xyzin,
        report,
        source_kind="xyzin",
        source_path=xyzin,
        report_path=tmp_path / "gf.report",
    )
    ped = read_gf_ped_section(xyzin)

    assert report.hessian_source == f"#CARTESIAN_HESSIAN in {xyzin}"
    assert report.result.coordinate_source.startswith("xyzin-frozen-gic:")
    assert len(report.result.frequencies_cm) == definition.rank
    assert ped.source_kind == "xyzin"
    assert ped.source_path == xyzin
    assert len(ped.modes) == definition.rank
