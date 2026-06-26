from __future__ import annotations

from oracle_chem import write_validation_section
from oracle_core import section_content


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
    assert validation[0] == "SCHEMA oracle.xyz.validation.v1"
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
