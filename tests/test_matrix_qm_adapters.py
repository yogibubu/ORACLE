from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from matrix_chem import (
    detect_qm_output_format,
    preprocess_to_enriched_xyz,
    read_geometry,
    write_validation_section,
)
from matrix_core import read_basic_section, read_sectioned_lines, section_content
from matrix_molpro import read_molpro_output_geometry, summarize_molpro_output
from matrix_mrcc import read_mrcc_output_geometry, summarize_mrcc_output
from matrix_orca import (
    hessian_input_from_orca_output,
    promote_orca_output_to_xyzin,
    read_orca_output_geometry,
    summarize_orca_output,
)
from matrix_gf import read_gf_ped_section, run_xyzin_gf_report_from_xyzin, write_gf_ped_section_from_report
from matrix_neo import write_gicforge_build_sections
from matrix_qm import read_cartesian_hessian_section
from tools import matrix_run


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
    assert section_content(lines, "TOPOLOGY")[0] == "SCHEMA oracle.xyz.topology.v1"


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
