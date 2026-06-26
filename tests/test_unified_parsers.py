from __future__ import annotations

from pathlib import Path

from oracle_chem import MolecularGeometry, read_enriched_xyz, read_geometry, read_xyz, read_zmatrix
from oracle_gaussian import (
    parse_gaussian_topology,
    read_gaussian_cartesian_input,
    read_gaussian_input,
    read_gaussian_log_geometry,
    read_gaussian_zmatrix_input,
    summarize_gaussian_log,
)


CORPUS = Path(__file__).resolve().parent / "fixtures" / "test_molecules" / "molecules"


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


def test_zmatrix_parser_handles_variables_and_dummy_atoms(tmp_path):
    path = tmp_path / "water.zmat"
    path.write_text(
        "\n".join(
            [
                "O",
                "X 1 1.0",
                "H 1 rOH 2 aHOX",
                "H 1 rOH 2 aHOX 3 dih",
                "",
                "Variables:",
                "rOH = 9.6D-01",
                "aHOX = 52.25",
                "dih = 180.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    geometry = read_zmatrix(path)

    assert geometry.atoms == ("O", "H", "H")
    assert geometry.metadata["dummy_atoms"] == 1
    assert geometry.coordinates_angstrom.shape == (3, 3)


def test_gaussian_zmatrix_input_uses_unified_zmatrix_adapter(tmp_path):
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
                "H 1 rOH",
                "H 1 rOH 2 aHOH",
                "",
                "rOH=0.96",
                "aHOH=104.5",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    geometry = read_gaussian_zmatrix_input(path)

    assert geometry.atoms == ("O", "H", "H")
    assert geometry.charge == 0
    assert geometry.multiplicity == 1
    assert geometry.source_format == "gaussian_zmatrix_input"


def test_gaussian_input_auto_detects_cartesian_with_zmat_route():
    geometry = read_gaussian_input(CORPUS / "h2o2zmat.inp")

    assert geometry.atoms == ("H", "O", "O", "H")
    assert geometry.charge == 0
    assert geometry.multiplicity == 1
    assert geometry.source_format == "gaussian_cartesian_input"


def test_gaussian_input_ignores_following_gic_block_after_zmatrix():
    geometry = read_gaussian_input(CORPUS / "c6h6_zmat_gic.inp")

    assert geometry.natoms == 12
    assert geometry.atoms.count("C") == 6
    assert geometry.atoms.count("H") == 6
    assert geometry.source_format == "gaussian_zmatrix_input"


def test_read_geometry_dispatches_gaussian_inp_to_unified_adapter():
    geometry = read_geometry(CORPUS / "pyridine_zmat.inp")

    assert geometry.natoms == 11
    assert geometry.atoms.count("N") == 1
    assert geometry.source_format == "gaussian_zmatrix_input"


def test_gaussian_zmatrix_variables_allow_blank_separated_blocks():
    geometry = read_gaussian_input(CORPUS / "h2co_zmat.inp")

    assert geometry.atoms == ("O", "C", "H", "H")
    assert geometry.source_format == "gaussian_zmatrix_input"


def test_gaussian_zmatrix_variables_are_case_insensitive():
    geometry = read_gaussian_input(CORPUS / "alaIIN.inp")

    assert geometry.natoms == 13
    assert geometry.atoms.count("N") == 1
    assert geometry.source_format == "gaussian_zmatrix_input"


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


def test_gaussian_log_geometry_dispatch_reads_last_orientation(tmp_path):
    path = tmp_path / "job.out"
    path.write_text(
        "\n".join(
            [
                " Input orientation:",
                " ---------------------------------------------------------------------",
                " Center     Atomic      Atomic             Coordinates (Angstroms)",
                " Number     Number       Type             X           Y           Z",
                " ---------------------------------------------------------------------",
                "      1          6           0        0.100000    0.200000    0.300000",
                "      2          1           0        0.100000    0.200000    1.300000",
                " ---------------------------------------------------------------------",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    geometry = read_gaussian_log_geometry(path)

    assert geometry.atoms == ("C", "H")
    assert geometry.source_format == "gaussian_log_orientation"


def test_gaussian_topology_parser_uses_cm5_and_mayer_only():
    data = parse_gaussian_topology(
        [
            " Hirshfeld charges, spin densities, dipoles, and CM5 charges",
            "   1  C   -0.0100   0.0000  -0.1234",
            "   2  H    0.0100   0.0000   0.0456",
            "",
            " Mayer bond orders and valences:",
            "  B( 1-C, 2-H)=0.9123",
            "",
            " Total bond order between atoms:",
            "  B( 1-C, 2-H)=1.5000",
        ]
    )

    assert data.cm5_charges == {1: -0.1234, 2: 0.0456}
    assert data.mayer_bond_orders == {(1, 2): 0.9123}
