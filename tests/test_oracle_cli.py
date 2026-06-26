from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

from oracle_babel import rdkit_available
from tools import oracle_run


ROOT = Path(__file__).resolve().parents[1]
GIC_CORPUS = ROOT / "tests" / "fixtures" / "test_molecules" / "molecules"


def test_python_module_entrypoint_prints_help():
    result = subprocess.run(
        [sys.executable, "-m", "oracle", "--help"],
        cwd=ROOT,
        env=os.environ.copy(),
        check=True,
        capture_output=True,
        text=True,
    )

    assert "ORACLE workflow CLI" in result.stdout
    assert "babel" in result.stdout
    assert "gicforge" in result.stdout


def test_validate_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyz"

    def fake_validate(target, *, require_fragments=False):
        calls["target"] = target
        calls["require_fragments"] = require_fragments
        return SimpleNamespace(status="PASS")

    monkeypatch.setattr("oracle_chem.write_validation_section", fake_validate)

    rc = oracle_run.main(["validate", str(path), "--require-fragments"])

    assert rc == 0
    assert calls == {"target": path, "require_fragments": True}
    assert "(PASS)" in capsys.readouterr().out


def test_lcb25_fetch_cli_calls_sync(tmp_path, monkeypatch, capsys):
    calls = {}

    def fake_sync(root, *, datasets=None, force=False):
        calls["root"] = root
        calls["datasets"] = datasets
        calls["force"] = force
        return Path(root) / "manifest.json"

    monkeypatch.setattr("oracle_babel.sync_lcb25_library", fake_sync)

    rc = oracle_run.main(
        ["lcb25", "fetch", "--root", str(tmp_path / "cache"), "--dataset", "se", "--force"]
    )

    assert rc == 0
    assert calls == {"root": tmp_path / "cache", "datasets": ["se"], "force": True}
    assert "Synced LCB25 library" in capsys.readouterr().out


def test_babel_preprocess_cli_calls_shared_pipeline(tmp_path, monkeypatch, capsys):
    calls = {}
    source = tmp_path / "source.inp"
    output = tmp_path / "molecule.xyzin"

    def fake_preprocess(target_source, target_output, *, source_kind="auto", symmetry_thresholds):
        calls["source"] = target_source
        calls["output"] = target_output
        calls["source_kind"] = source_kind
        calls["thresholds"] = symmetry_thresholds
        return SimpleNamespace(
            path=target_output,
            geometry=SimpleNamespace(natoms=3),
            point_group="C1",
            topology_bond_count=2,
            ring_count=0,
        )

    monkeypatch.setattr("oracle_chem.preprocess_to_enriched_xyz", fake_preprocess)

    rc = oracle_run.main(
        [
            "babel",
            "preprocess",
            str(source),
            str(output),
            "--source-kind",
            "auto",
            "--symmetry-distance",
            "0.002",
            "--symmetry-inertia",
            "0.003",
            "--max-rotation-order",
            "8",
        ]
    )

    assert rc == 0
    assert calls["source"] == source
    assert calls["output"] == output
    assert calls["source_kind"] == "auto"
    assert calls["thresholds"].distance_angstrom == 0.002
    assert calls["thresholds"].inertia_relative == 0.003
    assert calls["thresholds"].max_rotation_order == 8
    assert "Preprocessed ORACLE-Babel molecule" in capsys.readouterr().out


def test_fragments_plan_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyz"

    def fake_write(target):
        calls["target"] = target

    monkeypatch.setattr("oracle_fragments.write_fragment_plan_section", fake_write)

    rc = oracle_run.main(["fragments", "plan", str(path)])

    assert rc == 0
    assert calls == {"target": path}
    assert "Planned ORACLE fragment workflow" in capsys.readouterr().out


def test_fragments_build_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyzin"

    class FakeDefinition:
        reference_fragment = "F001"
        fragments = (object(), object())

    def fake_write(target):
        calls["target"] = target
        return FakeDefinition()

    monkeypatch.setattr("oracle_fragments.write_fragment_build_section", fake_write)

    rc = oracle_run.main(["fragments", "build", str(path)])

    assert rc == 0
    assert calls == {"target": path}
    assert "Built ORACLE fragments" in capsys.readouterr().out


def test_fragments_centers_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyzin"

    class FakeDefinition:
        centers = (object(), object(), object())
        interactions = (object(),)

    def fake_write(target):
        calls["target"] = target
        return FakeDefinition()

    monkeypatch.setattr("oracle_fragments.write_interaction_center_section", fake_write)

    rc = oracle_run.main(["fragments", "centers", str(path)])

    assert rc == 0
    assert calls == {"target": path}
    assert "Built ORACLE interaction centers" in capsys.readouterr().out


def test_rovib_summarize_cli_prints_summary(tmp_path, capsys):
    path = tmp_path / "molecule.xyzin"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#ROTATIONAL",
                "A_MHz = 1000.0",
                "B_MHz = 900.0",
                "C_MHz = 800.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rc = oracle_run.main(["rovib", "summarize", str(path)])

    assert rc == 0
    assert "rotational: A=1000MHz B=900MHz C=800MHz" in capsys.readouterr().out


def test_gicforge_plan_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyz"

    def fake_write(target, *, symmetrize=False, sycart=False):
        calls["target"] = target
        calls["symmetrize"] = symmetrize
        calls["sycart"] = sycart

    monkeypatch.setattr("oracle_gicforge.write_gicforge_plan_sections", fake_write)

    rc = oracle_run.main(["gicforge", "plan", str(path), "--symmetrize", "--sycart"])

    assert rc == 0
    assert calls == {"target": path, "symmetrize": True, "sycart": True}
    assert "Planned GICForge workflow" in capsys.readouterr().out


def test_gicforge_build_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyzin"

    class FakeDefinition:
        gics = (object(), object())
        rank = 2

    def fake_write(target, *, symmetrize=False, sycart=False):
        calls["target"] = target
        calls["symmetrize"] = symmetrize
        calls["sycart"] = sycart
        return FakeDefinition()

    monkeypatch.setattr("oracle_gicforge.write_gicforge_build_sections", fake_write)

    rc = oracle_run.main(["gicforge", "build", str(path), "--symmetrize", "--sycart"])

    assert rc == 0
    assert calls == {"target": path, "symmetrize": True, "sycart": True}
    assert "Built GICForge definition" in capsys.readouterr().out


def test_gicforge_cli_passes_improper_dihedrals_flag(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyzin"

    class FakeDefinition:
        gics = (object(),)
        rank = 1

    def fake_plan(target, *, symmetrize=False, sycart=False, improper_dihedrals=False):
        calls["plan"] = (target, symmetrize, sycart, improper_dihedrals)

    def fake_build(target, *, symmetrize=False, sycart=False, improper_dihedrals=False):
        calls["build"] = (target, symmetrize, sycart, improper_dihedrals)
        return FakeDefinition()

    monkeypatch.setattr("oracle_gicforge.write_gicforge_plan_sections", fake_plan)
    monkeypatch.setattr("oracle_gicforge.write_gicforge_build_sections", fake_build)

    rc_plan = oracle_run.main(["gicforge", "plan", str(path), "--improper-dihedrals"])
    rc_build = oracle_run.main(["gicforge", "build", str(path), "--improper-dihedrals"])

    assert rc_plan == 0
    assert rc_build == 0
    assert calls["plan"] == (path, False, False, True)
    assert calls["build"] == (path, False, False, True)
    assert "Built GICForge definition" in capsys.readouterr().out


def test_gf_cli_runs_xyzin_report_and_csv_export(tmp_path, monkeypatch, capsys):
    calls = {}
    fchk = tmp_path / "job.fchk"
    xyzin = tmp_path / "molecule.xyzin"
    out = tmp_path / "gf.txt"
    csv_dir = tmp_path / "csv"

    def fake_run(
        fchk_path,
        xyzin_path,
        *,
        scale_path=None,
        scale_records=(),
        local=False,
        force_threshold=None,
        block_by_irrep=False,
        subtract_electrostatic=False,
        subtract_uff_vdw=False,
        nonbonded_14_scale=0.5,
    ):
        calls["run"] = (
            fchk_path,
            xyzin_path,
            scale_path,
            scale_records,
            local,
            force_threshold,
            block_by_irrep,
            subtract_electrostatic,
            subtract_uff_vdw,
            nonbonded_14_scale,
        )
        return SimpleNamespace(text="gf report")

    def fake_csv(report, outdir, *, prefix="gf"):
        calls["csv"] = (report.text, outdir, prefix)
        return {"frequencies.csv": outdir / "gic_gf_frequencies.csv"}

    monkeypatch.setattr("oracle_gf.run_xyzin_gf_report_from_fchk", fake_run)
    monkeypatch.setattr("oracle_gf.write_csv_tables", fake_csv)

    rc = oracle_run.main(
        [
            "gf",
            "--fchk",
            str(fchk),
            "--xyzin",
            str(xyzin),
            "--out",
            str(out),
            "--csv-dir",
            str(csv_dir),
            "--scale",
            "GIC003=0.9",
            "--local",
            "--symmetry-blocks",
            "--force-threshold",
            "1e-8",
            "--subtract-electrostatic",
            "--subtract-uff-vdw",
            "--nonbonded-14-scale",
            "0.25",
        ]
    )

    assert rc == 0
    assert calls["run"] == (
        fchk,
        xyzin,
        None,
        ("GIC003=0.9",),
        True,
        1.0e-8,
        True,
        True,
        True,
        0.25,
    )
    assert calls["csv"] == ("gf report", csv_dir, "gic_gf")
    assert out.read_text(encoding="utf-8") == "gf report\n"
    assert "Wrote GF/PED report" in capsys.readouterr().out


def test_gicforge_bmatrix_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyzin"
    output = tmp_path / "bmat.out"

    class FakeMatrix:
        rows = ((1.0, 2.0, 3.0),)
        cartesian_columns = ("1:X", "1:Y", "1:Z")

    def fake_write(target, out):
        calls["target"] = target
        calls["output"] = out
        return FakeMatrix()

    monkeypatch.setattr("oracle_gicforge.write_gic_b_matrix", fake_write)

    rc = oracle_run.main(["gicforge", "bmatrix", str(path), str(output)])

    assert rc == 0
    assert calls == {"target": path, "output": output}
    assert "Wrote GIC B matrix" in capsys.readouterr().out


def test_gicforge_report_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyzin"
    output = tmp_path / "gicforge_report.txt"

    def fake_write(target, out):
        calls["target"] = target
        calls["output"] = out
        return out

    monkeypatch.setattr("oracle_gicforge.write_gic_report", fake_write)

    rc = oracle_run.main(["gicforge", "report", str(path), str(output)])

    assert rc == 0
    assert calls == {"target": path, "output": output}
    assert "Wrote GICForge report" in capsys.readouterr().out


def test_gicforge_corpus_cli_prints_inventory(capsys):
    rc = oracle_run.main(["gicforge", "corpus", "--root", str(GIC_CORPUS)])

    out = capsys.readouterr().out
    assert rc == 0
    assert "TOTAL_FILES 153" in out
    assert "SUFFIX .inp 126" in out
    assert "ROLE legacy_gic_input 126" in out


def test_gicforge_corpus_cli_filters_paths(capsys):
    rc = oracle_run.main(
        [
            "gicforge",
            "corpus",
            "--root",
            str(GIC_CORPUS),
            "--format",
            "paths",
            "--suffix",
            ".inp",
            "--limit",
            "2",
        ]
    )

    lines = capsys.readouterr().out.splitlines()
    assert rc == 0
    assert len(lines) == 2
    assert all(line.endswith(".inp") for line in lines)


def test_gicforge_corpus_audit_cli_prints_parser_budget(capsys):
    rc = oracle_run.main(["gicforge", "corpus-audit", "--root", str(GIC_CORPUS)])

    out = capsys.readouterr().out
    assert rc == 0
    assert "TOTAL_FILES 129" in out
    if rdkit_available():
        assert "PASS " in out
        assert "FAIL " in out
    else:
        assert "PASS 114" in out
        assert "FAIL 15" in out
    assert "SOURCE_FORMAT gaussian_zmatrix_input " in out


def test_gicforge_corpus_audit_cli_prints_failures(capsys):
    rc = oracle_run.main(
        [
            "gicforge",
            "corpus-audit",
            "--root",
            str(GIC_CORPUS),
            "--format",
            "failures",
            "--limit",
            "1",
        ]
    )

    lines = capsys.readouterr().out.splitlines()
    assert rc == 0
    assert len(lines) == 1
    assert lines[0].startswith("FAIL ")
    assert "Error" in lines[0]


def test_gicforge_gaussian_input_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    xyzin = tmp_path / "molecule.xyzin"
    output = tmp_path / "molecule.gjf"

    def fake_write(target, out, *, route, title=None, charge=None, multiplicity=None):
        calls["target"] = target
        calls["output"] = out
        calls["route"] = route
        calls["title"] = title
        calls["charge"] = charge
        calls["multiplicity"] = multiplicity
        return out

    monkeypatch.setattr("oracle_gicforge.write_gicforge_gaussian_input", fake_write)

    rc = oracle_run.main(
        [
            "gicforge",
            "gaussian-input",
            str(xyzin),
            str(output),
            "--route",
            "#p hf/3-21g opt",
            "--title",
            "job",
            "--charge",
            "1",
            "--multiplicity",
            "2",
        ]
    )

    assert rc == 0
    assert calls == {
        "target": xyzin,
        "output": output,
        "route": "#p hf/3-21g opt",
        "title": "job",
        "charge": 1,
        "multiplicity": 2,
    }
    assert "Wrote Gaussian input" in capsys.readouterr().out
