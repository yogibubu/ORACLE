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


def test_gaussian_summary_cli_prints_log_markers(tmp_path, capsys):
    log = tmp_path / "job.log"
    log.write_text(
        "\n".join(
            [
                "SCF Done:  E(RHF) =  -75.0",
                "Frequencies -- 100.0 200.0 300.0",
                "Normal termination of Gaussian 16",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rc = oracle_run.main(["gaussian", "summary", str(log)])
    out = capsys.readouterr().out

    assert rc == 0
    assert "normal_termination: 1" in out
    assert "scf_count: 1" in out
    assert "frequencies: 3" in out


def test_gaussian_status_cli_calls_adapter(tmp_path, monkeypatch, capsys):
    calls = {}

    def fake_status(workdir):
        calls["workdir"] = workdir
        return SimpleNamespace(
            status="running",
            workdir=workdir,
            log_path=workdir / "gauin.log",
            input_path=workdir / "gauin.gjf",
            pid=123,
            normal_termination=False,
            error_termination=False,
            message="running",
        )

    monkeypatch.setattr("oracle_gaussian.gaussian_job_status", fake_status)

    rc = oracle_run.main(["gaussian", "status", str(tmp_path)])
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {"workdir": tmp_path}
    assert "status: running" in out
    assert "pid: 123" in out


def test_gaussian_run_cli_calls_adapter(tmp_path, monkeypatch, capsys):
    calls = {}

    def fake_run(workdir, *, executable=None, input_path=None, background=False, timeout=None):
        calls["workdir"] = workdir
        calls["executable"] = executable
        calls["input_path"] = input_path
        calls["background"] = background
        calls["timeout"] = timeout
        return SimpleNamespace(
            input_path=workdir / "gauin.gjf",
            log_path=workdir / "gauin.log",
            pid=None,
            exit_code=0,
            success=True,
            message="ok",
        )

    monkeypatch.setattr("oracle_gaussian.run_gaussian_job", fake_run)

    rc = oracle_run.main(
        [
            "gaussian",
            "run",
            str(tmp_path),
            "--executable",
            "g16",
            "--input",
            str(tmp_path / "job.gjf"),
            "--timeout",
            "1.5",
        ]
    )
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {
        "workdir": tmp_path,
        "executable": "g16",
        "input_path": tmp_path / "job.gjf",
        "background": False,
        "timeout": 1.5,
    }
    assert "success: 1" in out
    assert "message: ok" in out


def test_gaussian_formchk_cli_calls_adapter(tmp_path, monkeypatch, capsys):
    calls = {}
    chk = tmp_path / "job.chk"
    fchk = tmp_path / "job.fchk"

    def fake_formchk(source, target=None, *, executable, timeout=None):
        calls["source"] = source
        calls["target"] = target
        calls["executable"] = executable
        calls["timeout"] = timeout
        return target

    monkeypatch.setattr("oracle_gaussian.formchk_checkpoint", fake_formchk)

    rc = oracle_run.main(
        [
            "gaussian",
            "formchk",
            str(chk),
            str(fchk),
            "--executable",
            "formchk-test",
            "--timeout",
            "2",
        ]
    )

    assert rc == 0
    assert calls == {
        "source": chk,
        "target": fchk,
        "executable": "formchk-test",
        "timeout": 2.0,
    }
    assert f"Wrote formatted checkpoint: {fchk}" in capsys.readouterr().out


def test_gaussian_fchk_summary_cli_calls_qff_reader(tmp_path, monkeypatch, capsys):
    calls = {}
    fchk = tmp_path / "job.fchk"

    def fake_read(path):
        calls["path"] = path
        return SimpleNamespace(
            atomic_numbers=[1, 8, 1],
            cartesian_hessian_lower=[0.0] * 45,
            harmonic_frequencies_cm=[1.0, 2.0],
            anharmonic_frequencies_cm=[1.1],
            anharmonic_e2=[0.0] * 6,
            normal_modes=[0.0] * 9,
        )

    monkeypatch.setattr("oracle_gaussian.read_gaussian_fchk_qff", fake_read)

    rc = oracle_run.main(["gaussian", "fchk-summary", str(fchk)])
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {"path": fchk}
    assert "atoms: 3" in out
    assert "anharmonic_frequencies: 1" in out


def test_gaussian_promote_fchk_cli_calls_adapter(tmp_path, monkeypatch, capsys):
    calls = {}
    fchk = tmp_path / "job.fchk"
    xyzin = tmp_path / "mol.xyzin"

    def fake_promote(
        fchk_path,
        xyzin_path,
        *,
        write_cartesian_hessian=True,
        write_normal_modes=True,
        write_qff=True,
    ):
        calls["fchk"] = fchk_path
        calls["xyzin"] = xyzin_path
        calls["write_cartesian_hessian"] = write_cartesian_hessian
        calls["write_normal_modes"] = write_normal_modes
        calls["write_qff"] = write_qff
        return SimpleNamespace(
            fchk_path=fchk_path,
            xyzin=xyzin_path,
            wrote_cartesian_hessian=write_cartesian_hessian,
            wrote_normal_modes=write_normal_modes,
            wrote_qff=write_qff,
        )

    monkeypatch.setattr("oracle_gaussian.promote_gaussian_fchk_to_xyzin", fake_promote)

    rc = oracle_run.main(
        [
            "gaussian",
            "promote-fchk",
            str(fchk),
            str(xyzin),
            "--no-normal-modes",
        ]
    )
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {
        "fchk": fchk,
        "xyzin": xyzin,
        "write_cartesian_hessian": True,
        "write_normal_modes": False,
        "write_qff": True,
    }
    assert "wrote_normal_modes: 0" in out


def test_gaussian_promote_rovib_cli_calls_adapter(tmp_path, monkeypatch, capsys):
    calls = {}
    log = tmp_path / "job.log"
    xyzin = tmp_path / "mol.xyzin"

    def fake_promote(
        log_path,
        xyzin_path,
        *,
        write_vibrational=True,
        write_rotational=True,
        write_deltabvib=True,
        invert_imaginary_modes=True,
        exclude_modes=(),
    ):
        calls["log"] = log_path
        calls["xyzin"] = xyzin_path
        calls["write_vibrational"] = write_vibrational
        calls["write_rotational"] = write_rotational
        calls["write_deltabvib"] = write_deltabvib
        calls["invert_imaginary_modes"] = invert_imaginary_modes
        calls["exclude_modes"] = exclude_modes
        return SimpleNamespace(
            log_path=log_path,
            xyzin=xyzin_path,
            wrote_vibrational=write_vibrational,
            wrote_rotational=write_rotational,
            wrote_deltabvib=write_deltabvib,
        )

    monkeypatch.setattr("oracle_gaussian.promote_gaussian_rovib_to_xyzin", fake_promote)

    rc = oracle_run.main(
        [
            "gaussian",
            "promote-rovib",
            str(log),
            str(xyzin),
            "--no-rotational",
            "--no-invert-imaginary",
            "--exclude-mode",
            "3",
        ]
    )
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {
        "log": log,
        "xyzin": xyzin,
        "write_vibrational": True,
        "write_rotational": False,
        "write_deltabvib": True,
        "invert_imaginary_modes": False,
        "exclude_modes": (3,),
    }
    assert "wrote_rotational: 0" in out


def test_gf_cli_can_use_hessian_from_xyzin_without_fchk(tmp_path, monkeypatch, capsys):
    calls = {}
    xyzin = tmp_path / "molecule.xyzin"

    def fake_report(
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
        calls["xyzin"] = xyzin_path
        calls["local"] = local
        calls["block_by_irrep"] = block_by_irrep
        return SimpleNamespace(text="GF from xyzin", result=SimpleNamespace())

    monkeypatch.setattr("oracle_gf.run_xyzin_gf_report_from_xyzin", fake_report)

    rc = oracle_run.main(["gf", "--xyzin", str(xyzin), "--local", "--symmetry-blocks"])
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {"xyzin": xyzin, "local": True, "block_by_irrep": True}
    assert "GF from xyzin" in out


def test_vpt2_vci_cli_uses_qff_section_from_xyzin(tmp_path, monkeypatch, capsys):
    calls = {}
    xyzin = tmp_path / "molecule.xyzin"
    report_path = tmp_path / "vpt2_vci.report"
    csv_dir = tmp_path / "csv"

    def fake_load_force_field(*, fchk_path=None, qff_path=None, xyzin_path=None):
        calls["load"] = (fchk_path, qff_path, xyzin_path)
        return SimpleNamespace(harmonic_frequencies_cm=(100.0,))

    def fake_run(force_field, *, max_quanta, roots, options=None, vci_method="dense"):
        calls["run"] = (force_field.harmonic_frequencies_cm, max_quanta, roots, vci_method)
        return SimpleNamespace(text="vpt2 vci report")

    def fake_csv(report, outdir, *, prefix="vpt2_vci"):
        calls["csv"] = (report.text, outdir, prefix)
        return {"comparison.csv": outdir / "vpt2_vci_comparison.csv"}

    monkeypatch.setattr("oracle_vpt2_vci.load_force_field", fake_load_force_field)
    monkeypatch.setattr("oracle_vpt2_vci.run_vpt2_vci_report", fake_run)
    monkeypatch.setattr("oracle_vpt2_vci.write_csv_tables", fake_csv)

    rc = oracle_run.main(
        [
            "vpt2-vci",
            "--xyzin",
            str(xyzin),
            "--max-quanta",
            "3",
            "--roots",
            "4",
            "--vci-method",
            "davidson",
            "--out",
            str(report_path),
            "--csv-dir",
            str(csv_dir),
        ]
    )

    assert rc == 0
    assert calls["load"] == (None, None, xyzin)
    assert calls["run"] == ((100.0,), 3, 4, "davidson")
    assert calls["csv"] == ("vpt2 vci report", csv_dir, "vpt2_vci")
    assert report_path.read_text(encoding="utf-8") == "vpt2 vci report\n"
    assert "Wrote VPT2/VCI report" in capsys.readouterr().out


def test_vpt2_vci_cli_writes_manifest_and_xyzin_section(tmp_path, monkeypatch, capsys):
    from oracle_vpt2_vci import read_vpt2_vci_section

    calls = {}
    xyzin = tmp_path / "molecule.xyzin"
    run_dir = tmp_path / "run"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")

    def fake_load_force_field(*, fchk_path=None, qff_path=None, xyzin_path=None):
        return SimpleNamespace(harmonic_frequencies_cm=(100.0,))

    def fake_run(force_field, *, max_quanta, roots, options=None, vci_method="dense"):
        calls["run"] = (max_quanta, roots, vci_method)
        return SimpleNamespace(text="vpt2 vci report")

    def fake_csv(report, outdir, *, prefix="vpt2_vci"):
        paths = {
            "frequencies.csv": outdir / f"{prefix}_frequencies.csv",
            "comparison.csv": outdir / f"{prefix}_comparison.csv",
            "mode_contributions.csv": outdir / f"{prefix}_mode_contributions.csv",
        }
        paths["frequencies.csv"].write_text("mode,harmonic_frequency_cm-1\n1,100.0\n", encoding="utf-8")
        paths["comparison.csv"].write_text(
            "root,vpt2_abs_cm-1,vci_abs_cm-1,delta_abs_cm-1,"
            "vpt2_exc_cm-1,vci_exc_cm-1,delta_exc_cm-1\n"
            "1,50.0,50.0,0.0,0.0,0.0,0.0\n",
            encoding="utf-8",
        )
        paths["mode_contributions.csv"].write_text("root,mode,expected_quanta\n1,1,0.5\n", encoding="utf-8")
        return paths

    monkeypatch.setattr("oracle_vpt2_vci.load_force_field", fake_load_force_field)
    monkeypatch.setattr("oracle_vpt2_vci.run_vpt2_vci_report", fake_run)
    monkeypatch.setattr("oracle_vpt2_vci.write_csv_tables", fake_csv)

    rc = oracle_run.main(
        [
            "vpt2-vci",
            "--xyzin",
            str(xyzin),
            "--run-dir",
            str(run_dir),
            "--max-quanta",
            "3",
            "--roots",
            "4",
            "--vci-method",
            "davidson",
        ]
    )
    out = capsys.readouterr().out
    section = read_vpt2_vci_section(xyzin)

    assert rc == 0
    assert calls["run"] == (3, 4, "davidson")
    assert (run_dir / "vpt2_vci_manifest.json").is_file()
    assert section.status == "complete"
    assert section.outputs["comparison"] == run_dir / "vpt2_vci_comparison.csv"
    assert "manifest:" in out
    assert "Updated #VPT2_VCI" in out


def test_vpt2_vci_collect_cli_refreshes_section(tmp_path, capsys):
    from oracle_vpt2_vci import read_vpt2_vci_section, vpt2_vci_section_from_run, write_vpt2_vci_section

    xyzin = tmp_path / "molecule.xyzin"
    run_dir = tmp_path / "vpt2_vci"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")
    run_dir.mkdir()
    (run_dir / "vpt2_vci.report").write_text("report\n", encoding="utf-8")
    (run_dir / "vpt2_vci_frequencies.csv").write_text(
        "mode,harmonic_frequency_cm-1\n1,100.0\n",
        encoding="utf-8",
    )
    (run_dir / "vpt2_vci_comparison.csv").write_text(
        "root,vpt2_abs_cm-1,vci_abs_cm-1,delta_abs_cm-1,"
        "vpt2_exc_cm-1,vci_exc_cm-1,delta_exc_cm-1\n"
        "1,50.0,50.0,0.0,0.0,0.0,0.0\n",
        encoding="utf-8",
    )
    (run_dir / "vpt2_vci_mode_contributions.csv").write_text(
        "root,mode,expected_quanta\n1,1,0.5\n",
        encoding="utf-8",
    )
    write_vpt2_vci_section(
        xyzin,
        vpt2_vci_section_from_run(
            source_kind="xyzin",
            source_path=xyzin,
            run_dir=run_dir,
            report_path=run_dir / "vpt2_vci.report",
            csv_dir=run_dir,
            manifest_path=None,
            max_quanta=2,
            roots=1,
            vci_method="dense",
            status="prepared",
        ),
    )

    rc = oracle_run.main(["vpt2-vci", "--collect", str(xyzin)])
    out = capsys.readouterr().out
    section = read_vpt2_vci_section(xyzin)

    assert rc == 0
    assert section.status == "complete"
    assert section.outputs["frequencies"] == run_dir / "vpt2_vci_frequencies.csv"
    assert "status: complete" in out
    assert "Updated #VPT2_VCI" in out


def test_dvr_prepare_cli_writes_manifest_and_xyzin_section(tmp_path, monkeypatch, capsys):
    calls = {}
    log = tmp_path / "scan.log"
    xyzin = tmp_path / "molecule.xyzin"
    outdir = tmp_path / "dvr"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")

    def fake_manifest(request, args):
        calls["manifest_request"] = request
        calls["manifest_args"] = tuple(args)
        manifest = request.outdir / f"{request.normalized_prefix}_manifest.json"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text("{}\n", encoding="utf-8")
        return manifest

    monkeypatch.setattr("oracle_dvr.write_dvr_manifest", fake_manifest)

    rc = oracle_run.main(
        [
            "dvr",
            "prepare",
            str(log),
            "--outdir",
            str(outdir),
            "--prefix",
            "demo",
            "--boundary",
            "nonperiodic",
            "--no-rotconst",
            "--check-only",
            "--xyzin",
            str(xyzin),
        ]
    )
    out = capsys.readouterr().out
    request = calls["manifest_request"]

    assert rc == 0
    assert request.log_path == log
    assert request.outdir == outdir
    assert request.normalized_prefix == "demo"
    assert request.boundary == "nonperiodic"
    assert not request.compute_rotconst
    assert request.check_only
    assert "--gaussian-log" in calls["manifest_args"]
    assert "Updated #DVR" in out
    assert "#DVR" in xyzin.read_text(encoding="utf-8")


def test_dvr_collect_cli_refreshes_post_run_outputs(tmp_path, capsys):
    from oracle_dvr import DVRRequest, dvr_section_from_request, read_dvr_section, write_dvr_section

    xyzin = tmp_path / "molecule.xyzin"
    outdir = tmp_path / "dvr"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=outdir,
        figdir=tmp_path / "fig",
        prefix="demo",
    )
    write_dvr_section(xyzin, dvr_section_from_request(request))
    outdir.mkdir()
    (outdir / "demo_summary.txt").write_text("Mass-weighted path Hamiltonian\n", encoding="utf-8")
    (outdir / "demo_levels.csv").write_text(
        "state,energy_cm-1,energy_above_ground_cm-1\n0,10.0,0.0\n",
        encoding="utf-8",
    )
    (outdir / "demo_grid.csv").write_text(
        "grid,s_au,s_sqrtamu_angstrom,V_cm-1\n0,0.0,0.0,0.0\n",
        encoding="utf-8",
    )

    rc = oracle_run.main(["dvr", "collect", str(xyzin)])
    out = capsys.readouterr().out
    section = read_dvr_section(xyzin)

    assert rc == 0
    assert section.status == "complete"
    assert section.outputs["levels"] == outdir / "demo_levels.csv"
    assert "status: complete" in out
    assert "Updated #DVR" in out


def test_dvr_run_cli_uses_xyzin_request_and_collects_outputs(tmp_path, monkeypatch, capsys):
    from oracle_dvr import DVRRequest, dvr_section_from_request, read_dvr_section, write_dvr_section

    xyzin = tmp_path / "molecule.xyzin"
    outdir = tmp_path / "dvr"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=outdir,
        figdir=tmp_path / "fig",
        prefix="demo",
    )
    write_dvr_section(xyzin, dvr_section_from_request(request))
    calls = {}

    def fake_run(run_request, *, timeout=None):
        calls["request"] = run_request
        calls["timeout"] = timeout
        _write_cli_dvr_outputs(run_request.outdir, run_request.normalized_prefix)
        manifest = run_request.outdir / f"{run_request.normalized_prefix}_manifest.json"
        manifest.write_text("{}\n", encoding="utf-8")
        return SimpleNamespace(manifest_path=manifest, status="complete")

    monkeypatch.setattr("oracle_dvr.run_dvr_request", fake_run)

    rc = oracle_run.main(["dvr", "run", "--xyzin", str(xyzin), "--timeout", "3"])
    out = capsys.readouterr().out
    section = read_dvr_section(xyzin)

    assert rc == 0
    assert calls["request"].log_path == tmp_path / "scan.log"
    assert calls["request"].outdir == outdir
    assert calls["timeout"] == 3.0
    assert section.status == "complete"
    assert section.outputs["levels"] == outdir / "demo_levels.csv"
    assert "manifest:" in out
    assert "status: complete" in out
    assert "Updated #DVR" in out


def test_dvr_run_cli_builds_log_request_and_updates_xyzin(tmp_path, monkeypatch, capsys):
    from oracle_dvr import read_dvr_section

    log = tmp_path / "scan.log"
    xyzin = tmp_path / "molecule.xyzin"
    outdir = tmp_path / "dvr"
    log.write_text("scan\n", encoding="utf-8")
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")
    calls = {}

    def fake_run(run_request, *, timeout=None):
        calls["request"] = run_request
        _write_cli_dvr_outputs(run_request.outdir, run_request.normalized_prefix)
        manifest = run_request.outdir / f"{run_request.normalized_prefix}_manifest.json"
        manifest.write_text("{}\n", encoding="utf-8")
        return SimpleNamespace(manifest_path=manifest, status="complete")

    monkeypatch.setattr("oracle_dvr.run_dvr_request", fake_run)

    rc = oracle_run.main(
        [
            "dvr",
            "run",
            str(log),
            "--outdir",
            str(outdir),
            "--prefix",
            "demo",
            "--solver",
            "sinc-dvr",
            "--xyzin",
            str(xyzin),
        ]
    )
    section = read_dvr_section(xyzin)

    assert rc == 0
    assert calls["request"].log_path == log
    assert calls["request"].outdir == outdir
    assert calls["request"].solver == "sinc-dvr"
    assert section.status == "complete"
    assert section.manifest_path == outdir / "demo_manifest.json"
    assert "Updated #DVR" in capsys.readouterr().out


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


def test_thermo_cli_calls_shared_pipeline(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyzin"
    report = tmp_path / "custom.report"

    def fake_run(target, *, report=True, report_path=None, write_section=True, **kwargs):
        calls["target"] = target
        calls["report"] = report
        calls["report_path"] = report_path
        calls["write_section"] = write_section
        calls.update(kwargs)
        return SimpleNamespace(total=SimpleNamespace(Q_dimless=12.5))

    monkeypatch.setattr("oracle_thermo.run_thermo_on_xyzin", fake_run)

    rc = oracle_run.main(
        [
            "thermo",
            str(path),
            "--out",
            str(report),
            "--no-write-section",
            "--cutoff-cm1",
            "20",
            "--keep-low-positive",
        ]
    )

    assert rc == 0
    assert calls == {
        "target": path,
        "report": True,
        "report_path": report,
        "write_section": False,
        "cutoff_cm1": 20.0,
        "keep_low_positive": True,
    }
    out = capsys.readouterr().out
    assert "Ran ORACLE Thermo" in out
    assert "Q=12.5" in out


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


def test_gicforge_fortran_audit_cli_prints_summary(tmp_path, monkeypatch, capsys):
    from oracle_gicforge import GICForgeFortranAudit, GICForgeFortranAuditResult

    calls = {}

    def fake_audit(**kwargs):
        calls.update(kwargs)
        return GICForgeFortranAudit(
            root=kwargs["root"],
            workdir=kwargs["workdir"],
            tolerance=kwargs["tolerance"],
            results=(
                GICForgeFortranAuditResult(
                    molecule="naphtalene.inp",
                    source=kwargs["root"] / "naphtalene.inp",
                    status="PASS",
                    oracle_rank=48,
                    fortran_rank=48,
                    oracle_row_rank=48,
                    fortran_row_rank=48,
                    row_space_residual=1.0e-10,
                ),
            ),
        )

    monkeypatch.setattr("oracle_gicforge.audit_gicforge_fortran_corpus", fake_audit)

    rc = oracle_run.main(
        [
            "gicforge",
            "fortran-audit",
            "--root",
            str(GIC_CORPUS),
            "--workdir",
            str(tmp_path / "audit"),
            "--molecule",
            "naphtalene.inp",
            "--tolerance",
            "1e-7",
        ]
    )
    out = capsys.readouterr().out

    assert rc == 0
    assert calls["root"] == GIC_CORPUS
    assert calls["workdir"] == tmp_path / "audit"
    assert calls["molecules"] == ["naphtalene.inp"]
    assert calls["tolerance"] == 1.0e-7
    assert "PASS 1" in out
    assert "CASE PASS naphtalene.inp" in out


def test_semiexp_benchmark_cli_calls_paper_artifact_generator(tmp_path, monkeypatch, capsys):
    calls = {}
    snapshot = tmp_path / "snapshot.json"
    outdir = tmp_path / "paper"

    def fake_generate(**kwargs):
        calls.update(kwargs)
        return (
            {"cases": {"water": {}}, "planar_pair_diagnostics": {"azulene": {}}},
            {"summary_tex": outdir / "benchmark_summary.tex"},
        )

    monkeypatch.setattr("oracle_morpheus.generate_paper_benchmark_artifacts", fake_generate)

    rc = oracle_run.main(
        [
            "semiexp-benchmark",
            "--snapshot",
            str(snapshot),
            "--outdir",
            str(outdir),
            "--no-refresh",
            "--update-snapshot",
        ]
    )
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {
        "snapshot_path": snapshot,
        "outdir": outdir,
        "refresh_from_outputs": False,
        "update_snapshot": True,
    }
    assert "cases: 1" in out
    assert f"summary_tex: {outdir / 'benchmark_summary.tex'}" in out


def test_semiexp_ensemble_paper_cli_calls_artifact_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    job = tmp_path / "ensemble.toml"
    paper_dir = tmp_path / "paper"
    outdir = tmp_path / "analysis"

    def fake_write(job_path, target_paper_dir, **kwargs):
        calls["job"] = job_path
        calls["paper_dir"] = target_paper_dir
        calls.update(kwargs)
        return {"prior_comparison": paper_dir / "generated" / "prior.tex"}

    monkeypatch.setattr("oracle_morpheus.write_ensemble_jpcl_artifacts", fake_write)

    rc = oracle_run.main(
        [
            "semiexp-ensemble-paper",
            "--job",
            str(job),
            "--paper-dir",
            str(paper_dir),
            "--outdir",
            str(outdir),
            "--soft-prior-sigma",
            "0.002",
        ]
    )
    out = capsys.readouterr().out

    assert rc == 0
    assert calls == {
        "job": job,
        "paper_dir": paper_dir,
        "outdir": outdir,
        "soft_prior_sigma": 0.002,
    }
    assert "prior_comparison:" in out


def test_semiexp_ensemble_scan_cli_calls_scan_helpers(tmp_path, monkeypatch, capsys):
    calls = {}
    job = tmp_path / "ensemble.toml"
    outdir = tmp_path / "scan"

    def fake_prior(job_path, target, **kwargs):
        calls["prior"] = (job_path, target, kwargs)
        return [{"prior_sigma": 1.0e-4}, {"prior_sigma": 1.0e-3}]

    def fake_synthon(job_path, target, **kwargs):
        calls["synthon"] = (job_path, target, kwargs)
        return [{"synthon_threshold": 0.01}]

    monkeypatch.setattr("oracle_morpheus.run_ensemble_prior_scan", fake_prior)
    monkeypatch.setattr("oracle_morpheus.run_ensemble_synthon_threshold_scan", fake_synthon)

    rc_prior = oracle_run.main(
        [
            "semiexp-ensemble-prior-scan",
            "--job",
            str(job),
            "--outdir",
            str(outdir / "prior"),
            "--sigma",
            "0.0001",
            "--sigma",
            "0.001",
        ]
    )
    out_prior = capsys.readouterr().out
    rc_synthon = oracle_run.main(
        [
            "semiexp-ensemble-synthon-scan",
            "--job",
            str(job),
            "--outdir",
            str(outdir / "synthon"),
            "--threshold",
            "0.01",
        ]
    )
    out_synthon = capsys.readouterr().out

    assert rc_prior == 0
    assert rc_synthon == 0
    assert calls["prior"] == (job, outdir / "prior", {"sigmas": (1.0e-4, 1.0e-3)})
    assert calls["synthon"] == (job, outdir / "synthon", {"thresholds": (0.01,)})
    assert "rows: 2" in out_prior
    assert "rows: 1" in out_synthon


def test_multistructure_reference_search_cli_calls_library(tmp_path, monkeypatch, capsys):
    calls = {}
    query = tmp_path / "query.xyz"
    outdir = tmp_path / "matches"

    def fake_search(query_xyz, **kwargs):
        calls["query"] = query_xyz
        calls.update(kwargs)
        return SimpleNamespace(
            library_root=tmp_path / "se_geometries",
            matches=(SimpleNamespace(slug="water", similarity_combined=0.99),),
            skipped=(),
        )

    monkeypatch.setattr("oracle_morpheus.search_reference_library", fake_search)

    rc = oracle_run.main(
        [
            "multistructure-reference-search",
            "--query-xyz",
            str(query),
            "--outdir",
            str(outdir),
            "--top-k",
            "3",
            "--ring-weight",
            "0.4",
        ]
    )
    out = capsys.readouterr().out

    assert rc == 0
    assert calls["query"] == query
    assert calls["outdir"] == outdir
    assert calls["top_k"] == 3
    assert calls["include_ring_comparison"] is True
    assert calls["ring_weight"] == 0.4
    assert "top_match: water similarity=0.99" in out


def test_multistructure_build_reference_geometry_cli_calls_library(tmp_path, monkeypatch, capsys):
    calls = {}
    query = tmp_path / "query.xyz"
    outdir = tmp_path / "assisted"

    def fake_build(query_xyz, **kwargs):
        calls["query"] = query_xyz
        calls.update(kwargs)
        return SimpleNamespace(
            library_root=tmp_path / "se_geometries",
            targets=(object(), object()),
            unmatched=(object(),),
            iterations=4,
            rms_target_residual_final=1.2e-4,
        )

    monkeypatch.setattr("oracle_morpheus.build_reference_assisted_geometry", fake_build)

    rc = oracle_run.main(
        [
            "multistructure-build-reference-geometry",
            "--query-xyz",
            str(query),
            "--outdir",
            str(outdir),
            "--top-library-matches",
            "5",
            "--max-fragment-matches",
            "2",
            "--apply-kind",
            "bond",
            "--apply-kind",
            "angle",
            "--no-ring-comparison",
        ]
    )
    out = capsys.readouterr().out

    assert rc == 0
    assert calls["query"] == query
    assert calls["outdir"] == outdir
    assert calls["top_library_matches"] == 5
    assert calls["max_fragment_matches"] == 2
    assert calls["apply_kinds"] == ("bond", "angle")
    assert calls["include_ring_comparison"] is False
    assert "targets: 2" in out
    assert "unmatched: 1" in out


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


def _write_cli_dvr_outputs(outdir: Path, prefix: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{prefix}_summary.txt").write_text(
        "Mass-weighted path Hamiltonian\n",
        encoding="utf-8",
    )
    (outdir / f"{prefix}_levels.csv").write_text(
        "state,energy_cm-1,energy_above_ground_cm-1\n0,10.0,0.0\n",
        encoding="utf-8",
    )
    (outdir / f"{prefix}_grid.csv").write_text(
        "grid,s_au,s_sqrtamu_angstrom,V_cm-1\n0,0.0,0.0,0.0\n",
        encoding="utf-8",
    )
