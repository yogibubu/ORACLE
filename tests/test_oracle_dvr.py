from __future__ import annotations

from pathlib import Path

from oracle_dvr import (
    DVRRequest,
    build_fortran_bridge_args,
    build_fortran_shell_command,
    build_path_analysis_args,
    dvr_section_from_request,
    read_dvr_section,
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
