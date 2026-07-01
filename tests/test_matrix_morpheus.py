from __future__ import annotations

import json
from pathlib import Path


def test_matrix_morpheus_imports_and_data_library():
    import matrix_morpheus as morpheus

    assert morpheus.SEMIEXP_JOB_SCHEMA == "oracle.semiexp.job.v1"
    assert morpheus.XYZIN_ISOTOPOLOGUES_SCHEMA == "oracle.xyz.isotopologues.v1"
    assert morpheus.DEFAULT_SE_GEOMETRY_LIBRARY.is_dir()
    assert (morpheus.DEFAULT_SE_GEOMETRY_LIBRARY / "manifest.csv").is_file()
    assert morpheus.DEFAULT_SEMIEXP_BENCHMARK_SNAPSHOT.parts[:2] == (
        "benchmarks",
        "semiexp_msr",
    )


def test_semiexp_paper_benchmark_paths_resolve_from_repo_root():
    from matrix_morpheus.paper_benchmarks import repository_root

    root = Path(__file__).resolve().parents[1]

    assert repository_root() == root


def test_semiexp_paper_benchmark_snapshot_generates_tables(tmp_path):
    from matrix_morpheus import generate_paper_benchmark_artifacts

    root = Path(__file__).resolve().parents[1]
    snapshot, artifacts = generate_paper_benchmark_artifacts(
        outdir=tmp_path / "paper",
        refresh_from_outputs=False,
    )

    assert snapshot["schema"] == "oracle.semiexp.paper_regression.v1"
    assert tuple(snapshot["cases"]) == (
        "glycolaldehyde",
        "cyclopentadiene",
        "nitrobenzene",
        "p-EBN",
        "azulene",
        "norcamphor",
        "glycine_I",
        "glycine_II",
    )
    assert set(snapshot["planar_pair_diagnostics"]) == {"nitrobenzene", "azulene"}
    assert artifacts["summary_tex"].is_file()
    assert artifacts["planar_tex"].is_file()
    assert "Nitrobenzene" in artifacts["summary_tex"].read_text(encoding="utf-8")
    paper_generated = root / "docs/papers/morpheus_semiexp/generated"
    for name, artifact in (
        ("paper_benchmark_summary.csv", artifacts["summary_csv"]),
        ("paper_planar_pair_diagnostics.csv", artifacts["planar_csv"]),
        ("benchmark_summary.tex", artifacts["summary_tex"]),
        ("planar_pair_diagnostics.tex", artifacts["planar_tex"]),
    ):
        assert artifact.read_text(encoding="utf-8") == (
            paper_generated / name
        ).read_text(encoding="utf-8")


def test_legacy_msr_zmatrix_constraints_import_for_p_ebn():
    from matrix_morpheus.msr_legacy import read_msr_legacy_input

    root = Path(__file__).resolve().parents[1]
    legacy = read_msr_legacy_input(
        root / "benchmarks/semiexp_msr/inputs/legacy_msr_import/p-EBN.msr.inp"
    )

    assert legacy.geometry.source_format == "msr_legacy_zmatrix"
    assert len(legacy.geometry.atoms) == 15
    assert len(legacy.observations) == 11
    assert "CC2=R(1,2)" in legacy.geometry.fixed_parameters
    assert "cc2 - cc4 Frozen" in legacy.geometry.fixed_parameters


def test_semiexp_job_accepts_oracle_and_legacy_schemas(tmp_path):
    from matrix_morpheus import read_semiexperimental_job

    body = """
title = "water inline"

[geometry]
units = "angstrom"
atoms = [
  ["O", 0.0, 0.0, 0.0],
  ["H", 0.0, 0.0, 0.9572],
  ["H", 0.9266, 0.0, -0.2396],
]
"""
    oracle_job = tmp_path / "water_oracle.mse.toml"
    legacy_job = tmp_path / "water_legacy.mse.toml"
    oracle_job.write_text('schema = "oracle.semiexp.job.v1"\n' + body, encoding="utf-8")
    legacy_job.write_text('schema = "merlino.semiexp.job.v1"\n' + body, encoding="utf-8")

    assert read_semiexperimental_job(oracle_job).geometry.source_format == "matrix_morpheus_job"
    assert read_semiexperimental_job(legacy_job).geometry.source_format == "matrix_morpheus_job"


def test_prepare_semiexp_xyzin_supports_custom_path(tmp_path):
    from matrix_morpheus import prepare_semiexperimental_xyzin, read_geometry_input

    root = Path(__file__).resolve().parents[1]
    xyzin = tmp_path / "water_container"
    result = prepare_semiexperimental_xyzin(
        root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz",
        observations_source=root
        / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml",
        xyzin_path=xyzin,
    )

    assert result.xyzin == xyzin
    assert "SCHEMA oracle.xyz.isotopologues.v1" in xyzin.read_text(encoding="utf-8")
    assert read_geometry_input(xyzin).source_format == "xyzin"


def test_primitive_class_advisor_maps_user_classes_to_gics():
    from matrix_morpheus import derive_primitive_class_plan, parse_primitive_class_spec

    labels = (
        "GIC001 AStr [ 0.75*R(1,2)+0.10*R(7,13)]",
        "GIC002 AStr [ 0.25*R(1,2)+0.74*R(7,13)]",
        "GIC003 AStr [ 0.74*R(1,2)+0.74*R(7,13)]",
        "GIC004 AStr [ 0.60*R(1,2)]",
    )
    classes = (
        parse_primitive_class_spec("CC_skeleton:R(2,1)"),
        parse_primitive_class_spec("CO_stretches:R(7,13)"),
    )

    plan = derive_primitive_class_plan(
        labels,
        classes,
        min_fraction=0.70,
        cross_fraction_max=0.20,
    )

    assert [(item.name, item.patterns) for item in plan.parameter_classes] == [
        ("CO_stretches", ("GIC002", "GIC003")),
        ("CC_skeleton", ("GIC001",)),
    ]
    assert plan.fixed_patterns == ("GIC004",)
    assert plan.rejected_labels == ("GIC004",)


def test_primitive_class_advisor_respects_data_budget():
    from matrix_morpheus import derive_primitive_class_plan, parse_primitive_class_spec

    labels = (
        "GIC001 AStr [ 0.90*R(1,2)]",
        "GIC002 AStr [ 0.80*R(2,3)]",
        "GIC003 AStr [ 0.75*R(3,4)]",
        "GIC004 AStr [ 0.78*R(4,5)]",
    )
    classes = (
        parse_primitive_class_spec("class_a:R(1,2)|R(2,3)"),
        parse_primitive_class_spec("class_b:R(3,4)"),
        parse_primitive_class_spec("class_c:R(4,5)"),
    )

    plan = derive_primitive_class_plan(labels, classes, max_classes=1)

    assert [(item.name, item.patterns) for item in plan.parameter_classes] == [
        ("class_a", ("GIC001", "GIC002")),
    ]
    assert plan.fixed_patterns == ("GIC003", "GIC004")


def test_primitive_class_advisor_supports_nonbond_primitives():
    from matrix_morpheus import derive_primitive_class_plan, parse_primitive_class_spec

    labels = (
        "GIC001 ABend [ 0.80*A(1,2,3)]",
        "GIC002 ATor [ 0.72*D(1,2,3,4)]",
        "GIC003 ABend [ 0.69*A(5,6,7)]",
    )
    classes = (
        parse_primitive_class_spec("bends:A(3,2,1)"),
        parse_primitive_class_spec("torsions:D(4,3,2,1)"),
    )

    plan = derive_primitive_class_plan(labels, classes)

    assert [(item.name, item.patterns) for item in plan.parameter_classes] == [
        ("bends", ("GIC001",)),
        ("torsions", ("GIC002",)),
    ]
    assert plan.fixed_patterns == ("GIC003",)


def test_oracle_semiexp_cli_runs_water_gic(tmp_path):
    from matrix_core.cli import main
    from matrix_morpheus import read_morpheus_section

    root = Path(__file__).resolve().parents[1]
    outdir = tmp_path / "run"
    xyzin = tmp_path / "xyzin"
    status = main(
        [
            "semiexp",
            "--xyz",
            str(root / "packages/matrix-morpheus/examples/semiexp/water/parent.xyz"),
            "--observations",
            str(root / "packages/matrix-morpheus/examples/semiexp/water/isotopologues.toml"),
            "--xyzin",
            str(xyzin),
            "--outdir",
            str(outdir),
            "--coordinate-model",
            "gic",
            "--max-iter",
            "2",
        ],
        repo_root=root,
    )

    manifest = json.loads((outdir / "semiexp_manifest.json").read_text(encoding="utf-8"))
    checkpoint = json.loads((outdir / "semiexp_checkpoint.json").read_text(encoding="utf-8"))
    section = read_morpheus_section(xyzin)

    assert status == 0
    assert manifest["schema_version"] == "oracle.run.v1"
    assert checkpoint["schema"] == "oracle.semiexp.checkpoint.v1"
    assert section.status == "complete"
    assert section.run_dir == outdir
    assert section.manifest_path == outdir / "semiexp_manifest.json"
    assert section.html_report_path == outdir / "semiexp_report.html"
    assert section.latex_tables_path == outdir / "semiexp_tables.tex"
    assert section.coordinate_model == "gic"
    assert section.iterations <= 2
    assert section.parameter_count > 0
    assert "ORACLE semiexperimental equilibrium-geometry fit" in (
        outdir / "semiexp_report.txt"
    ).read_text(encoding="utf-8")
    assert (outdir / "semiexp_report.html").is_file()
