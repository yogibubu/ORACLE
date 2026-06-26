from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from oracle_chem import detect_qm_output_format, preprocess_to_enriched_xyz, read_geometry
from oracle_core import read_basic_section, read_sectioned_lines, section_content
from oracle_molpro import read_molpro_output_geometry, summarize_molpro_output
from oracle_mrcc import read_mrcc_output_geometry, summarize_mrcc_output
from tools import oracle_run


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


def test_read_geometry_auto_dispatches_qm_outputs(tmp_path):
    molpro = tmp_path / "molpro.out"
    mrcc = tmp_path / "mrcc.out"
    molpro.write_text(_molpro_output(), encoding="utf-8")
    mrcc.write_text(_mrcc_output(), encoding="utf-8")

    assert detect_qm_output_format(molpro) == "molpro"
    assert detect_qm_output_format(mrcc) == "mrcc"
    assert read_geometry(molpro).source_format == "molpro_output"
    assert read_geometry(mrcc).source_format == "mrcc_output"


def test_babel_preprocess_imports_molpro_with_basic_section(tmp_path):
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

    monkeypatch.setattr("oracle_molpro.promote_molpro_output_to_xyzin", fake_promote)

    rc = oracle_run.main(
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

    monkeypatch.setattr("oracle_mrcc.summarize_mrcc_output", fake_summary)

    rc = oracle_run.main(["mrcc", "summary", str(output)])
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {"source": output}
    assert "multiplicity: 2" in out
