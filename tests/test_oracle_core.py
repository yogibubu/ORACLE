from __future__ import annotations

from oracle_core import (
    ORACLE_MANIFEST_SCHEMA,
    ORACLE_XYZ_ISOTOPOLOGUES_SCHEMA,
    MERLINO_XYZIN_ISOTOPOLOGUES_SCHEMA,
    PLANNED_FRAMEWORK_EXPANSION,
    PLANNED_FRAMEWORK_NAME,
    XyzinIsotopologueRecord,
    build_run_manifest,
    BasicSection,
    ensure_workspace,
    parse_basic_section,
    parse_xyzin_isotopologue_records,
    read_basic_section,
    replace_section,
    replace_xyz_block,
    section_content,
    tool_contract,
    tool_contract_markdown_table,
    tool_contract_readiness,
    tool_contracts,
    validate_xyzin_isotopologue_records,
    write_basic_section,
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


def test_tool_contract_registry_records_standalone_sections_and_future_names():
    contracts = {contract.key: contract for contract in tool_contracts()}

    assert PLANNED_FRAMEWORK_NAME == "MATRIX"
    assert PLANNED_FRAMEWORK_EXPANSION == (
        "Molecular Analysis Toolkit for Reusable Integrated eXperiments"
    )
    assert contracts["gicforge"].planned_name == "NEO"
    assert contracts["gicforge"].expanded_name == "Nonredundant Equivariant Orthogonalizer"
    assert contracts["gicforge"].produced_sections == ("GIC", "SYCART")
    assert contracts["gui"].planned_name == "ORACLE"
    assert contracts["gui"].expanded_name == (
        "Operator for Routing, Analysis, Control, Launch and Exploration"
    )
    assert contracts["trinity"].status == "prepare-only"
    assert contracts["trinity"].produced_sections == ("TRINITY",)
    assert tool_contract("NEO").key == "gicforge"
    assert "TRINITY" in tool_contract_markdown_table()


def test_tool_contract_readiness_checks_required_xyzin_sections(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0.0 0.0 0.0",
                "",
                "#BASIC",
                "SCHEMA oracle.xyz.basic.v1",
                "",
                "#ROTATIONAL",
                "SCHEMA oracle.xyz.rotational.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    thermo = tool_contract_readiness(path, "thermo")
    gf = tool_contract_readiness(path, "gf")

    assert thermo.ready
    assert thermo.missing_required_sections == ()
    assert not gf.ready
    assert gf.missing_required_sections == ("GIC", "CARTESIAN_HESSIAN")


def test_basic_section_accepts_merlino_aligned_key_values():
    section = parse_basic_section(
        [
            "CHARGE              -1",
            "SPIN_MULTIPLICITY   2",
            "POINT_GROUP         Cs",
            "Watson Reduction A",
            "T_K = 150.0",
            "P_atm = 0.5",
        ]
    )

    assert section.charge == -1
    assert section.multiplicity == 2
    assert section.point_group == "Cs"
    assert section.watson_reduction == "A"
    assert section.temperature_K == 150.0
    assert section.pressure_atm == 0.5


def test_basic_section_writer_preserves_other_sections(tmp_path):
    path = tmp_path / "molecule.xyzin"
    path.write_text("1\nh\nH 0 0 0\n\n#GIC\nSCHEMA oracle.xyz.gic.v1\n", encoding="utf-8")

    write_basic_section(path, BasicSection(charge=1, multiplicity=2, point_group="C1"))
    lines = path.read_text(encoding="utf-8").splitlines()
    parsed = read_basic_section(path)

    assert parsed.charge == 1
    assert parsed.multiplicity == 2
    assert section_content(lines, "GIC")[0] == "SCHEMA oracle.xyz.gic.v1"


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


def test_replace_xyz_block_preserves_sections_after_avogadro_edit(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "2",
                "old",
                "H 0 0 0",
                "H 0 0 1",
                "",
                "#SMILES",
                "[H][H]",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    replace_xyz_block(path, ["2", "edited in Avogadro", "H 0 0 0", "H 0 0 2"])
    lines = path.read_text(encoding="utf-8").splitlines()

    assert lines[:4] == ["2", "edited in Avogadro", "H 0 0 0", "H 0 0 2"]
    assert section_content(lines, "SMILES") == ["[H][H]", ""]
    assert section_content(lines, "TOPOLOGY") == ["SCHEMA oracle.xyz.topology.v1"]
