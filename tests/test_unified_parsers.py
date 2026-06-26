from __future__ import annotations

import pytest

from oracle_chem import MolecularGeometry, read_enriched_xyz, read_xyz
from oracle_gaussian import read_gaussian_cartesian_input, summarize_gaussian_log


def test_xyz_parser_normalizes_symbols_and_numbers(tmp_path):
    path = tmp_path / "mixed.xyz"
    path.write_text(
        "\n".join(
            [
                "3",
                "mixed tokens",
                "8 0.0 0.0 0.0",
                "h 0.0 0.0 1.0",
                "CL 1.0 0.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    geometry = read_xyz(path)

    assert isinstance(geometry, MolecularGeometry)
    assert geometry.atoms == ("O", "H", "Cl")
    assert geometry.coordinates_angstrom.shape == (3, 3)
    assert geometry.source_format == "xyz"


def test_enriched_xyz_parser_reads_first_block_and_reports_sections(tmp_path):
    path = tmp_path / "xyzin"
    path.write_text(
        "\n".join(
            [
                "2",
                "h2",
                "H 0 0 0",
                "H 0 0 1",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
                "",
                "#GIC",
                "SCHEMA oracle.xyz.gic.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    geometry = read_enriched_xyz(path)

    assert geometry.atoms == ("H", "H")
    assert geometry.metadata["sections"] == ("TOPOLOGY", "GIC")


def test_gaussian_cartesian_input_parser_returns_shared_geometry(tmp_path):
    path = tmp_path / "water.gjf"
    path.write_text(
        "\n".join(
            [
                "%nproc=4",
                "#p hf/sto-3g opt",
                "",
                "water title",
                "",
                "0 1",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
                "",
                "B 1 2 F",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    geometry = read_gaussian_cartesian_input(path)

    assert isinstance(geometry, MolecularGeometry)
    assert geometry.atoms == ("O", "H", "H")
    assert geometry.charge == 0
    assert geometry.multiplicity == 1
    assert geometry.fixed_parameters == ("B(1,2)",)
    assert geometry.metadata["route"] == ("#p hf/sto-3g opt",)


def test_gaussian_zmatrix_input_is_rejected_until_unified_adapter_exists(tmp_path):
    path = tmp_path / "zmat.gjf"
    path.write_text(
        "\n".join(
            [
                "#p hf/sto-3g geom=zmat",
                "",
                "zmatrix",
                "",
                "0 1",
                "O",
                "H 1 r1",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Z-matrix"):
        read_gaussian_cartesian_input(path)


def test_gaussian_log_summary_uses_shared_geometry_for_last_orientation(tmp_path):
    path = tmp_path / "test.log"
    path.write_text(
        "\n".join(
            [
                " SCF Done:  E(RHF) = -1.000000 A.U.",
                " Standard orientation:",
                " ---------------------------------------------------------------------",
                " Center     Atomic      Atomic             Coordinates (Angstroms)",
                " Number     Number       Type             X           Y           Z",
                " ---------------------------------------------------------------------",
                "      1          8           0        0.000000    0.000000    0.000000",
                "      2          1           0        0.000000    0.000000    1.000000",
                " ---------------------------------------------------------------------",
                " Frequencies -- 100.0 200.0 300.0",
                " Normal termination of Gaussian 16",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = summarize_gaussian_log(path)

    assert summary.normal_termination is True
    assert summary.scf_energies_hartree == (-1.0,)
    assert summary.frequencies_cm == (100.0, 200.0, 300.0)
    assert summary.last_orientation is not None
    assert summary.last_orientation.atoms == ("O", "H")
