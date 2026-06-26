from __future__ import annotations

from oracle_chem import preprocess_to_enriched_xyz, read_enriched_xyz
from oracle_core import read_sectioned_lines, section_content


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
    assert "POINT_GROUP C1" in section_content(lines, "SYMMETRY")
    assert section_content(lines, "TOPOLOGY")[0] == "SCHEMA oracle.xyz.topology.v1"
    assert section_content(lines, "SYNTHONS")[0] == "SCHEMA oracle.xyz.synthons.v1"
    assert result.topology_bond_count == 2

    geometry = read_enriched_xyz(target)
    assert geometry.atoms == ("O", "H", "H")
    assert set(geometry.metadata["sections"]) >= {"SOURCE", "SYMMETRY", "TOPOLOGY", "SYNTHONS"}
