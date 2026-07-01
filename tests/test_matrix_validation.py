from __future__ import annotations

from matrix_chem import MATRIX_XYZ_VALIDATION_SCHEMA, write_validation_section
from matrix_core import section_content


def test_validation_section_passes_after_required_sections(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "2",
                "h2",
                "H 0 0 0",
                "H 0 0 1",
                "",
                "#SYMMETRY",
                "SCHEMA oracle.xyz.symmetry.v1",
                "POINT_GROUP C1",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
                "INDEXING ATOMS=ONE_BASED",
                "[BONDS]",
                "1 2",
                "[BOND_ORDERS]",
                "1 2 1.0",
                "[RINGS]",
                "NONE",
                "",
                "#SYNTHONS",
                "SCHEMA oracle.xyz.synthons.v1",
                "INDEXING ATOMS=ONE_BASED",
                "",
                "#GIC",
                "SCHEMA oracle.xyz.gic.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = write_validation_section(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    assert result.status == "PASS"
    validation = section_content(lines, "VALIDATION")
    assert validation[0] == f"SCHEMA {MATRIX_XYZ_VALIDATION_SCHEMA}"
    assert "STATUS PASS" in validation
    assert section_content(lines, "GIC")[0] == "SCHEMA oracle.xyz.gic.v1"


def test_validation_section_fails_when_topology_is_missing(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#SYMMETRY",
                "SCHEMA oracle.xyz.symmetry.v1",
                "",
                "#SYNTHONS",
                "SCHEMA oracle.xyz.synthons.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = write_validation_section(path)
    validation = section_content(path.read_text(encoding="utf-8").splitlines(), "VALIDATION")

    assert result.status == "FAIL"
    assert "STATUS FAIL" in validation
    assert any("MISSING_TOPOLOGY" in line for line in validation)


def test_validation_section_fails_when_topology_lacks_bond_order(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "2",
                "h2",
                "H 0 0 0",
                "H 0 0 1",
                "",
                "#SYMMETRY",
                "SCHEMA oracle.xyz.symmetry.v1",
                "",
                "#TOPOLOGY",
                "SCHEMA matrix.xyz.topology.v1",
                "INDEXING ATOMS=ONE_BASED",
                "[BONDS]",
                "1 2",
                "[BOND_ORDERS]",
                "NONE",
                "",
                "#SYNTHONS",
                "SCHEMA matrix.xyz.synthons.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = write_validation_section(path)
    validation = section_content(path.read_text(encoding="utf-8").splitlines(), "VALIDATION")

    assert result.status == "FAIL"
    assert any("MISSING_BOND_ORDERS" in line for line in validation)


def test_validation_section_checks_fragment_atom_coverage(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "3",
                "h3 fragment coverage",
                "H 0 0 0",
                "H 0 0 1",
                "H 0 0 2",
                "",
                "#SYMMETRY",
                "SCHEMA oracle.xyz.symmetry.v1",
                "",
                "#TOPOLOGY",
                "SCHEMA matrix.xyz.topology.v1",
                "INDEXING ATOMS=ONE_BASED",
                "[BONDS]",
                "1 2",
                "2 3",
                "[BOND_ORDERS]",
                "1 2 1.0",
                "2 3 1.0",
                "[RINGS]",
                "NONE",
                "",
                "#SYNTHONS",
                "SCHEMA matrix.xyz.synthons.v1",
                "",
                "#FRAGMENTS",
                "SCHEMA matrix.xyz.fragments.v1",
                "INDEXING ATOMS=ONE_BASED",
                "[FRAGMENTS]",
                "F001 LABEL=frag ATOMS=1,2 CENTER=0,0,0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = write_validation_section(path)
    validation = section_content(path.read_text(encoding="utf-8").splitlines(), "VALIDATION")

    assert result.status == "FAIL"
    assert any("INCOMPLETE_FRAGMENT_COVERAGE" in line for line in validation)
