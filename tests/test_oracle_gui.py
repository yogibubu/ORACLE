from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

from oracle_dvr import DVRRequest, dvr_section_from_request, write_dvr_section
from oracle_gui import (
    ORACLE_GUI_WINDOWS,
    WorkflowStatus,
    all_known_sections,
    avogadro_command,
    dvr_gui_state_lines,
    gaussian_promote_fchk_command,
    gaussian_promote_rovib_command,
    gicforge_build_command,
    load_dvr_gui_state,
    load_oracle_project_state,
    load_vpt2_vci_gui_state,
    molden_command,
    preprocess_command,
    project_state_lines,
    rovib_density_command,
    rovib_summary_command,
    rovib_vibrational_dos_command,
    vpt2_vci_gui_state_lines,
    window_spec,
)
from oracle_vpt2_vci import vpt2_vci_section_from_run, write_vpt2_vci_section


ROOT = Path(__file__).resolve().parents[1]


def test_gui_python_module_entrypoint_prints_help():
    result = subprocess.run(
        [sys.executable, "-m", "oracle_gui", "--help"],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
        check=True,
        capture_output=True,
        text=True,
    )

    assert "ORACLE desktop GUI" in result.stdout
    assert "xyzin" in result.stdout


def test_gui_dvr_state_reads_refreshed_dvr_section(tmp_path):
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
        "state,energy_cm-1,energy_above_ground_cm-1\n0,25.0,0.0\n",
        encoding="utf-8",
    )
    (outdir / "demo_grid.csv").write_text(
        "grid,s_au,s_sqrtamu_angstrom,V_cm-1\n0,0.0,0.0,0.0\n",
        encoding="utf-8",
    )

    state = load_dvr_gui_state(xyzin, refresh=True)
    lines = dvr_gui_state_lines(state)

    assert state.ready
    assert state.level_count == 1
    assert "status: complete" in lines
    assert "ready: 1" in lines


def test_gui_vpt2_vci_state_reads_refreshed_section(tmp_path):
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

    state = load_vpt2_vci_gui_state(xyzin, refresh=True)
    lines = vpt2_vci_gui_state_lines(state)

    assert state.ready
    assert state.root_count == 1
    assert state.mode_count == 1
    assert "status: complete" in lines
    assert "ready: 1" in lines


def test_gui_project_state_reports_xyzin_sections_and_workflows(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 0.9",
                "H 0.0 0.8 -0.2",
                "",
                "#BASIC",
                "SCHEMA oracle.xyz.basic.v1",
                "POINT_GROUP = C2V",
                "",
                "#SYMMETRY",
                "SCHEMA oracle.xyz.symmetry.v1",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
                "",
                "#SYNTHONS",
                "SCHEMA oracle.xyz.synthons.v1",
                "",
                "#GIC",
                "SCHEMA oracle.xyz.gic.v1",
                "STATUS BUILT",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    state = load_oracle_project_state(xyzin)
    lines = project_state_lines(state)

    assert state.exists
    assert state.atom_count == 3
    assert state.point_group == "C2V"
    assert state.section("GIC").schema == "oracle.xyz.gic.v1"
    assert state.workflow("babel").status == WorkflowStatus.COMPLETE
    assert state.workflow("gicforge").status == WorkflowStatus.COMPLETE
    assert any("workflow gicforge: complete" in line for line in lines)


def test_gui_project_state_missing_file_is_reported_without_crashing(tmp_path):
    state = load_oracle_project_state(tmp_path / "missing.xyzin")

    assert not state.exists
    assert state.validation_status == "MISSING"
    assert state.workflow("dashboard").status == WorkflowStatus.MISSING


def test_gui_window_specs_cover_primary_oracle_workflows():
    keys = {spec.key for spec in ORACLE_GUI_WINDOWS}

    assert {
        "dashboard",
        "babel",
        "topology",
        "gicforge",
        "gf",
        "sefit",
        "rovib_thermo",
        "anharmonic",
        "rotational_spectroscopy",
        "vibrational_spectroscopy",
        "electronic_spectroscopy",
        "thermochemistry_kinetics",
        "diagnostics",
    }.issubset(keys)
    assert window_spec("gicforge").produced_sections == ("GIC", "SYCART")
    assert window_spec("rotational_spectroscopy").category == "spectroscopy"
    assert "mode heat-map" in window_spec("vibrational_spectroscopy").publication_exports
    assert "Molden" in window_spec("electronic_spectroscopy").external_viewers
    assert "KINETICS" in window_spec("thermochemistry_kinetics").produced_sections
    assert "GIC" in all_known_sections()
    assert "KINETICS" in all_known_sections()


def test_gui_command_specs_use_oracle_cli_without_shell(tmp_path):
    source = tmp_path / "mol.xyz"
    xyzin = tmp_path / "mol.xyzin"
    fchk = tmp_path / "mol.fchk"
    log = tmp_path / "mol.log"

    preprocess = preprocess_command(source, xyzin, source_kind="xyz")
    avogadro = avogadro_command(xyzin, executable="avogadro")
    molden = molden_command(fchk)
    promote_fchk = gaussian_promote_fchk_command(fchk, xyzin, qff=False)
    promote_rovib = gaussian_promote_rovib_command(log, xyzin, exclude_modes=(1, 3))
    rovib_summary = rovib_summary_command(xyzin)
    vib_dos = rovib_vibrational_dos_command(xyzin, out=tmp_path / "dos_vib.dat")
    rovib_dos = rovib_density_command(
        xyzin,
        vib_dos=tmp_path / "dos_vib.dat",
        out=tmp_path / "dos_rovib.dat",
        jmax=20,
    )
    gic = gicforge_build_command(xyzin)

    assert preprocess.argv[1:4] == ("-m", "oracle", "babel")
    assert preprocess.argv[-2:] == ("--source-kind", "xyz")
    assert avogadro.argv == ("avogadro", str(xyzin))
    assert molden.argv == ("molden", str(fchk))
    assert promote_fchk.argv[-1] == "--no-qff"
    assert promote_fchk.produced_sections == ("CARTESIAN_HESSIAN", "NORMAL_MODES")
    assert promote_rovib.argv[-4:] == ("--exclude-mode", "1", "--exclude-mode", "3")
    assert promote_rovib.produced_sections == ("VIBRATIONAL", "ROTATIONAL", "DELTABVIB")
    assert rovib_summary.argv[-3:] == ("rovib", "summarize", str(xyzin))
    assert "dos" in vib_dos.argv
    assert rovib_dos.argv[-2:] == ("--jmax", "20")
    assert " ".join(gic.argv).endswith(
        "gicforge build "
        f"{xyzin} --symmetrize --sycart --improper-dihedrals"
    )
    assert gic.required_sections == ("SYMMETRY", "TOPOLOGY", "SYNTHONS")
