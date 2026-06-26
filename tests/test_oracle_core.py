from __future__ import annotations

from oracle_core import (
    ORACLE_MANIFEST_SCHEMA,
    ORACLE_XYZ_ISOTOPOLOGUES_SCHEMA,
    MERLINO_XYZIN_ISOTOPOLOGUES_SCHEMA,
    XyzinIsotopologueRecord,
    build_run_manifest,
    ensure_workspace,
    parse_xyzin_isotopologue_records,
    replace_section,
    section_content,
    validate_xyzin_isotopologue_records,
    xyzin_isotopologue_section_lines,
)


def test_workspace_layout(tmp_path):
    layout = ensure_workspace(tmp_path / "project")
    for name in ("inputs", "runs", "outputs", "reports", "cache", "logs"):
        assert (layout.root / name).is_dir()


def test_manifest_schema(tmp_path):
    manifest = build_run_manifest(workflow="smoke", status="completed", run_dir=tmp_path)
    assert manifest.schema_version == ORACLE_MANIFEST_SCHEMA
    assert manifest.to_dict()["workflow"] == "smoke"


def test_section_replacement_preserves_other_sections(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0 0 0",
                "H 0 0 1",
                "H 0 1 0",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
                "BOND 1 2",
                "",
                "#GIC",
                "SCHEMA oracle.xyz.gic.v1",
                "OLD true",
                "",
                "#MORPHEUS",
                "SCHEMA oracle.xyz.morpheus.v1",
                "STATUS draft",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    replace_section(path, "GIC", ["SCHEMA oracle.xyz.gic.v1", "UPDATED true"])
    lines = path.read_text(encoding="utf-8").splitlines()

    assert section_content(lines, "TOPOLOGY") == [
        "SCHEMA oracle.xyz.topology.v1",
        "BOND 1 2",
        "",
    ]
    assert section_content(lines, "GIC") == [
        "SCHEMA oracle.xyz.gic.v1",
        "UPDATED true",
    ]
    assert section_content(lines, "MORPHEUS") == [
        "SCHEMA oracle.xyz.morpheus.v1",
        "STATUS draft",
    ]


def test_isotopologue_records_write_oracle_schema_and_validate():
    records = (
        XyzinIsotopologueRecord(
            label="parent",
            rotational_MHz=(1000.0, 800.0, 600.0),
        ),
        XyzinIsotopologueRecord(
            label="D2",
            substitutions={2: 2},
            rotational_MHz=(990.0, 790.0, 590.0),
            deltavib_MHz=(0.1, 0.2, 0.3),
        ),
    )

    lines = xyzin_isotopologue_section_lines(records)
    assert lines[0] == f"SCHEMA {ORACLE_XYZ_ISOTOPOLOGUES_SCHEMA}"
    parsed = parse_xyzin_isotopologue_records(lines)

    assert parsed == records
    assert validate_xyzin_isotopologue_records(parsed, atom_count=3, require_rotational=True) == ()


def test_isotopologue_parser_accepts_merlino_schema():
    lines = [
        f"SCHEMA {MERLINO_XYZIN_ISOTOPOLOGUES_SCHEMA}",
        "UNITS ROTATIONAL=MHz DELTAVIB=MHz DELTAEL=MHz SIGMA=MHz",
        "INDEXING ATOMS=ONE_BASED",
        "BEGIN parent",
        "DEFINITION parent",
        "ROTATIONAL_MHZ A=1000 B=800 C=600",
        "END",
    ]

    parsed = parse_xyzin_isotopologue_records(lines)

    assert parsed[0].label == "parent"
    assert parsed[0].rotational_MHz == (1000.0, 800.0, 600.0)
