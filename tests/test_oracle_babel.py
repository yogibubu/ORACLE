from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from oracle_babel import (
    RDKitUnavailableError,
    extract_legacy_smiles_input,
    rdkit_available,
    smiles_to_geometry,
)
from oracle_chem import (
    MolecularGeometry,
    preprocess_to_enriched_xyz,
    read_enriched_xyz,
    read_geometry,
)
from oracle_core import read_sectioned_lines, section_content
from oracle_chem.symmetry import candidate_ops


CORPUS = Path(__file__).resolve().parent / "fixtures" / "test_molecules" / "molecules"


def test_babel_preprocess_writes_avogadro_compatible_enriched_xyz(tmp_path):
    source = tmp_path / "water.xyz"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    target = tmp_path / "oracle.xyz"

    result = preprocess_to_enriched_xyz(source, target)
    lines = read_sectioned_lines(target)

    assert result.path == target
    assert lines[:5] == [
        "3",
        "water",
        "O       0.00000000      0.00000000      0.00000000",
        "H       0.00000000      0.00000000      1.00000000",
        "H       0.00000000      1.00000000      0.00000000",
    ]
    assert section_content(lines, "SOURCE")[0] == "SCHEMA oracle.xyz.source.v1"
    symmetry = section_content(lines, "SYMMETRY")
    assert "POINT_GROUP C2v" in symmetry
    assert "OPERATION_COUNT 4" in symmetry
    assert "[OPERATIONS]" in symmetry
    assert any("LABEL=sigma_yz" in line for line in symmetry)
    assert section_content(lines, "TOPOLOGY")[0] == "SCHEMA oracle.xyz.topology.v1"
    assert section_content(lines, "SYNTHONS")[0] == "SCHEMA oracle.xyz.synthons.v1"
    assert result.topology_bond_count == 2

    geometry = read_enriched_xyz(target)
    assert geometry.atoms == ("O", "H", "H")
    assert set(geometry.metadata["sections"]) >= {"SOURCE", "SYMMETRY", "TOPOLOGY", "SYNTHONS"}


def test_babel_preprocess_imports_gaussian_cm5_and_mayer_topology(tmp_path):
    source = tmp_path / "water.log"
    source.write_text(
        "\n".join(
            [
                " Standard orientation:",
                " ---------------------------------------------------------------------",
                " Center     Atomic      Atomic             Coordinates (Angstroms)",
                " Number     Number       Type             X           Y           Z",
                " ---------------------------------------------------------------------",
                "      1          8           0        0.000000    0.000000    0.000000",
                "      2          1           0        0.000000    0.000000    0.960000",
                "      3          1           0        0.000000    0.760000   -0.580000",
                " ---------------------------------------------------------------------",
                " Hirshfeld charges, spin densities, dipoles, and CM5 charges",
                "   1  O   -0.0100   0.0000  -0.4000",
                "   2  H    0.0100   0.0000   0.2000",
                "   3  H    0.0100   0.0000   0.2000",
                "",
                " Mayer bond orders and valences:",
                "  B( 1-O, 2-H)=0.9000 B( 1-O, 3-H)=0.9100",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    target = tmp_path / "water.xyzin"

    result = preprocess_to_enriched_xyz(source, target)
    lines = read_sectioned_lines(target)
    gaussian = section_content(lines, "GAUSSIAN_TOPOLOGY")
    synthons = section_content(lines, "SYNTHONS")
    topology = section_content(lines, "TOPOLOGY")

    assert result.topology_bond_count == 2
    assert "CM5 1 -0.4000000000" in gaussian
    assert "BO_SOURCE = Mayer" in gaussian
    assert "CHARGE_SOURCE Gaussian CM5" in synthons
    assert "BOND_ORDER_SOURCE Gaussian Mayer" in topology
    oxygen_line = next(line for line in synthons if line.startswith("1 O "))
    assert float(oxygen_line.split()[3]) == -0.4


def test_symmetry_candidate_operations_include_nonprimitive_rotation_powers():
    labels = {label for label, _operation in candidate_ops(max_n=4)}

    assert "C4z^1" in labels
    assert "C4z^2" in labels
    assert "C4z^3" in labels
    assert "sigma_h*C4z^1" in labels


def test_legacy_smiles_input_extraction_reads_title_charge_and_smiles():
    smiles_input = extract_legacy_smiles_input(CORPUS / "pyrrole_smile1.inp")

    assert smiles_input.title == "Pyrrole"
    assert smiles_input.charge == 0
    assert smiles_input.multiplicity == 1
    assert smiles_input.smiles == "[nH]1cccc1"
    assert any("SMILES" in line for line in smiles_input.route_lines)


def test_read_geometry_dispatches_marked_smiles_to_oracle_babel(monkeypatch):
    calls = {}

    def fake_smiles_to_geometry(
        smiles,
        *,
        title="",
        charge=None,
        multiplicity=None,
        source_path=None,
        route_lines=(),
        random_seed=0,
    ):
        calls["smiles"] = smiles
        calls["title"] = title
        calls["charge"] = charge
        calls["multiplicity"] = multiplicity
        calls["source_path"] = source_path
        calls["route_lines"] = route_lines
        return MolecularGeometry(
            atoms=("N", "H"),
            coordinates_angstrom=np.zeros((2, 3)),
            comment=title,
            source_format="smiles_rdkit",
            source_path=source_path,
            charge=charge,
            multiplicity=multiplicity,
        )

    monkeypatch.setattr("oracle_babel.smiles.smiles_to_geometry", fake_smiles_to_geometry)

    geometry = read_geometry(CORPUS / "pyrrole_smile1.inp")

    assert calls["smiles"] == "[nH]1cccc1"
    assert calls["title"] == "Pyrrole"
    assert calls["charge"] == 0
    assert geometry.source_format == "smiles_rdkit"


def test_smiles_to_geometry_reports_missing_rdkit_when_unavailable():
    if rdkit_available():
        pytest.skip("RDKit is installed in this environment")

    with pytest.raises(RDKitUnavailableError, match="RDKit is required"):
        smiles_to_geometry("NCC(=O)O")
