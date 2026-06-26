from __future__ import annotations

import json
from pathlib import Path


def test_oracle_morpheus_imports_and_data_library():
    import oracle_morpheus as morpheus

    assert morpheus.SEMIEXP_JOB_SCHEMA == "oracle.semiexp.job.v1"
    assert morpheus.XYZIN_ISOTOPOLOGUES_SCHEMA == "oracle.xyz.isotopologues.v1"
    assert morpheus.DEFAULT_SE_GEOMETRY_LIBRARY.is_dir()
    assert (morpheus.DEFAULT_SE_GEOMETRY_LIBRARY / "manifest.csv").is_file()


def test_semiexp_job_accepts_oracle_and_legacy_schemas(tmp_path):
    from oracle_morpheus import read_semiexperimental_job

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

    assert read_semiexperimental_job(oracle_job).geometry.source_format == "oracle_morpheus_job"
    assert read_semiexperimental_job(legacy_job).geometry.source_format == "oracle_morpheus_job"


def test_prepare_semiexp_xyzin_supports_custom_path(tmp_path):
    from oracle_morpheus import prepare_semiexperimental_xyzin, read_geometry_input

    root = Path(__file__).resolve().parents[1]
    xyzin = tmp_path / "water_container"
    result = prepare_semiexperimental_xyzin(
        root / "packages/oracle-morpheus/examples/semiexp/water/parent.xyz",
        observations_source=root / "packages/oracle-morpheus/examples/semiexp/water/isotopologues.toml",
        xyzin_path=xyzin,
    )

    assert result.xyzin == xyzin
    assert "SCHEMA oracle.xyz.isotopologues.v1" in xyzin.read_text(encoding="utf-8")
    assert read_geometry_input(xyzin).source_format == "xyzin"


def test_oracle_semiexp_cli_runs_water_gic(tmp_path):
    from oracle_core.cli import main

    root = Path(__file__).resolve().parents[1]
    outdir = tmp_path / "run"
    xyzin = tmp_path / "xyzin"
    status = main(
        [
            "semiexp",
            "--xyz",
            str(root / "packages/oracle-morpheus/examples/semiexp/water/parent.xyz"),
            "--observations",
            str(root / "packages/oracle-morpheus/examples/semiexp/water/isotopologues.toml"),
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

    assert status == 0
    assert manifest["schema_version"] == "oracle.run.v1"
    assert checkpoint["schema"] == "oracle.semiexp.checkpoint.v1"
    assert "ORACLE semiexperimental equilibrium-geometry fit" in (
        outdir / "semiexp_report.txt"
    ).read_text(encoding="utf-8")
    assert (outdir / "semiexp_report.html").is_file()
