from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from oracle_dvr import (
    DVRRequest,
    build_fortran_bridge_args,
    build_fortran_shell_command,
    build_path_analysis_args,
    collect_dvr_outputs,
    dvr_output_summary_lines,
    dvr_request_from_section,
    dvr_section_from_request,
    read_dvr_section,
    refresh_dvr_section,
    run_dvr_request,
    write_dvr_section,
)
from oracle_engines import (
    DVR_FORTRAN_FILES,
    dvr_fortran_layout,
    puckering_dvr_layout,
    validate_dvr_sources,
    validate_puckering_dvr_backend,
)


ROOT = Path(__file__).resolve().parents[1]


def test_dvr_request_builds_python_args(tmp_path):
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=tmp_path / "out",
        figdir=tmp_path / "fig",
        prefix="demo",
        boundary="nonperiodic",
        solver="fourier",
        compute_rotconst=False,
        label_cremer_pople=False,
    )

    args = build_path_analysis_args(request)

    assert args[0] == str(
        tmp_path / "engines" / "puckering_dvr" / "scripts" / "mw_path_dvr.py"
    )
    assert "--gaussian-log" in args
    assert str(tmp_path / "scan.log") in args
    assert args[args.index("--solver") + 1] == "fourier"
    assert "--compute-rotconst" not in args
    assert "--label-cremer-pople" not in args


def test_fortran_solver_uses_sinc_python_grid_then_bridge(tmp_path):
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=tmp_path / "out",
        figdir=tmp_path / "fig",
        prefix="demo",
        boundary="nonperiodic",
        solver="fortran-gaussian",
        python_executable="/usr/bin/python3",
    )

    python_args = build_path_analysis_args(request)
    bridge_args = build_fortran_bridge_args(request, tmp_path / "bin" / "path_dvr.x")
    shell = build_fortran_shell_command(request, python_args, bridge_args)

    assert python_args[python_args.index("--solver") + 1] == "sinc-dvr"
    assert bridge_args[bridge_args.index("--grid-csv") + 1] == str(
        tmp_path / "out" / "demo_grid.csv"
    )
    assert bridge_args[-2:] == ["--mode", "gaussian"]
    assert "'/usr/bin/python3'" in shell
    assert "run_fortran_dvr.py" in shell


def test_dvr_section_roundtrip_preserves_workflow_request(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=tmp_path / "out",
        figdir=tmp_path / "fig",
        prefix="demo",
        boundary="nonperiodic",
        solver="sinc-dvr",
        compute_rotconst=False,
        label_cremer_pople=True,
        check_only=True,
    )

    write_dvr_section(
        xyzin,
        dvr_section_from_request(request, manifest_path=tmp_path / "out" / "demo_manifest.json"),
    )
    section = read_dvr_section(xyzin)

    assert section.log_path == tmp_path / "scan.log"
    assert section.grid_csv == tmp_path / "out" / "demo_grid.csv"
    assert section.manifest_path == tmp_path / "out" / "demo_manifest.json"
    assert section.boundary == "nonperiodic"
    assert section.solver == "sinc-dvr"
    assert not section.compute_rotconst
    assert section.label_cremer_pople
    assert section.check_only
    restored = dvr_request_from_section(section, repo_root=tmp_path, python_executable="py")
    assert restored.repo_root == tmp_path
    assert restored.log_path == request.log_path
    assert restored.outdir == request.outdir
    assert restored.figdir == request.figdir
    assert restored.solver == request.solver
    assert restored.python_executable == "py"


def test_run_dvr_request_executes_python_backend_and_writes_manifest(tmp_path, monkeypatch):
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=tmp_path / "out",
        figdir=tmp_path / "fig",
        prefix="demo",
        python_executable="/usr/bin/python3",
    )
    calls = []

    def fake_run(command, **kwargs):
        calls.append((tuple(command), kwargs))
        _write_dvr_outputs(request.outdir, "demo")
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr("oracle_dvr.workflow.subprocess.run", fake_run)

    result = run_dvr_request(request, timeout=2.5)
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))

    assert result.status == "complete"
    assert result.python_args[0].endswith("mw_path_dvr.py")
    assert result.commands[0].stdout == "ok"
    assert calls[0][0][0] == "/usr/bin/python3"
    assert calls[0][1]["timeout"] == 2.5
    assert manifest["workflow"] == "dvr"
    assert manifest["status"] == "complete"
    assert manifest["outputs"]["levels"] == str(request.outdir / "demo_levels.csv")


def test_run_dvr_request_executes_fortran_bridge_after_grid(tmp_path, monkeypatch):
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=tmp_path / "out",
        figdir=tmp_path / "fig",
        prefix="demo",
        solver="fortran-gaussian",
        python_executable="/usr/bin/python3",
    )
    calls = []

    def fake_run(command, **kwargs):
        calls.append(tuple(command))
        _write_dvr_outputs(request.outdir, "demo")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("oracle_dvr.workflow.subprocess.run", fake_run)
    monkeypatch.setattr(
        "oracle_dvr.workflow.resolve_dvr_executable",
        lambda repo_root: tmp_path / "bin" / "path_dvr.x",
    )

    result = run_dvr_request(request)

    assert len(calls) == 2
    assert calls[0][1].endswith("mw_path_dvr.py")
    assert calls[1][1].endswith("run_fortran_dvr.py")
    assert "--mode" in result.bridge_args
    assert "gaussian" in result.bridge_args


def test_dvr_output_reader_normalizes_post_run_files(tmp_path):
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=tmp_path / "out",
        figdir=tmp_path / "fig",
        prefix="demo",
    )
    _write_dvr_outputs(request.outdir, "demo")
    section = dvr_section_from_request(request, manifest_path=tmp_path / "out" / "demo_manifest.json")

    snapshot = collect_dvr_outputs(section)

    assert snapshot.status == "complete"
    assert snapshot.levels[1].energy_above_ground_cm == 12.5
    assert snapshot.grid[0].potential_cm == 0.0
    assert snapshot.grid[1].values["rotA"] == 9.9
    assert snapshot.expectations[0].values["rotA"] == 10.0
    assert "output:levels" in "\n".join(dvr_output_summary_lines(snapshot))


def test_refresh_dvr_section_adds_detected_outputs(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text("1\ncomment\nH 0.0 0.0 0.0\n", encoding="utf-8")
    request = DVRRequest(
        repo_root=tmp_path,
        log_path=tmp_path / "scan.log",
        outdir=tmp_path / "out",
        figdir=tmp_path / "fig",
        prefix="demo",
    )
    write_dvr_section(xyzin, dvr_section_from_request(request))
    _write_dvr_outputs(request.outdir, "demo")

    snapshot = refresh_dvr_section(xyzin)
    section = read_dvr_section(xyzin)

    assert snapshot.status == "complete"
    assert section.status == "complete"
    assert section.outputs["levels"] == request.outdir / "demo_levels.csv"
    assert section.outputs["expectations"] == request.outdir / "demo_expectations.csv"


def test_dvr_fortran_sources_are_vendored_under_oracle_engines():
    layout = dvr_fortran_layout(ROOT)

    assert validate_dvr_sources(ROOT) == ()
    assert set(DVR_FORTRAN_FILES) == {"path_dvr.f", "dvr_hqrii.f"}
    assert layout.source_dir == ROOT / "engines" / "fortran" / "dvr"
    assert layout.compile_script.is_file()


def test_puckering_dvr_backend_is_vendored_under_oracle_engines():
    layout = puckering_dvr_layout(ROOT)

    assert validate_puckering_dvr_backend(ROOT) == ()
    assert layout.engine_dir == ROOT / "engines" / "puckering_dvr"
    assert layout.path_analysis_script.is_file()
    assert layout.fortran_bridge_script.is_file()


def _write_dvr_outputs(outdir: Path, prefix: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{prefix}_summary.txt").write_text(
        "Mass-weighted path Hamiltonian\nLowest levels:\n",
        encoding="utf-8",
    )
    (outdir / f"{prefix}_levels.csv").write_text(
        "state,energy_cm-1,energy_above_ground_cm-1\n"
        "0,100.0,0.0\n"
        "1,112.5,12.5\n",
        encoding="utf-8",
    )
    (outdir / f"{prefix}_grid.csv").write_text(
        "grid,s_au,s_sqrtamu_angstrom,V_cm-1,rotA,psi_0,prob_density_0\n"
        "0,0.0,0.0,0.0,10.0,0.1,0.01\n"
        "1,0.5,0.1,3.5,9.9,0.2,0.04\n",
        encoding="utf-8",
    )
    (outdir / f"{prefix}_expectations.csv").write_text(
        "state,energy_cm-1,energy_above_ground_cm-1,rotA\n"
        "0,100.0,0.0,10.0\n"
        "1,112.5,12.5,9.8\n",
        encoding="utf-8",
    )
