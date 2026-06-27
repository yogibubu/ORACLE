from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

from oracle_chem import preprocess_to_enriched_xyz, write_validation_section
from oracle_core import XyzinIsotopologueRecord, write_xyzin_isotopologue_records
from oracle_dvr import DVRRequest, dvr_section_from_request, write_dvr_section
from oracle_gf import GFGICRow, GFModeRow, GFPEDSection, write_gf_ped_section
from oracle_gicforge import write_gicforge_build_sections
from oracle_morpheus import MorpheusSection, write_morpheus_section
from oracle_trinity import TrinitySection, write_trinity_section
from oracle_gui import (
    DashboardActionTemplate,
    ORACLE_GUI_WINDOWS,
    OracleDashboardController,
    OracleGICForgeController,
    OracleGFController,
    OracleGuiCommand,
    OracleSEFitController,
    OracleStructureController,
    OracleTrinityController,
    WorkflowStatus,
    all_known_sections,
    avogadro_command,
    default_gf_csv_dir,
    default_gf_report_output,
    default_gicforge_bmatrix_output,
    default_gicforge_gaussian_output,
    default_gicforge_report_output,
    dvr_gui_state_lines,
    default_sefit_outdir,
    default_trinity_run_dir,
    gaussian_promote_fchk_command,
    gaussian_promote_rovib_command,
    gf_command,
    gf_gui_state_lines,
    gicforge_gui_state_lines,
    gicforge_build_command,
    load_dvr_gui_state,
    load_gf_gui_state,
    load_gicforge_gui_state,
    load_oracle_project_state,
    load_structure_gui_state,
    load_tool_contract_gui_state,
    load_vpt2_vci_gui_state,
    load_workbench_gui_state,
    molden_command,
    preprocess_command,
    project_state_lines,
    rovib_density_command,
    rovib_summary_command,
    rovib_vibrational_dos_command,
    semiexp_command,
    sefit_gui_state_lines,
    load_sefit_gui_state,
    vpt2_vci_gui_state_lines,
    workbench_gui_state_lines,
    window_spec,
    trinity_prepare_command,
    trinity_gui_state_lines,
    load_trinity_gui_state,
    tool_contract_gui_state_lines,
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


def test_gui_tool_contract_state_reports_standalone_readiness(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text(
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

    state = load_tool_contract_gui_state(xyzin)
    lines = tool_contract_gui_state_lines(state)
    rows = {row[0]: row for row in state.table.rows}

    assert state.exists
    assert "ready tools:" in "\n".join(lines)
    assert rows["thermo"][3] == "yes"
    assert rows["thermo"][5] == "none"
    assert rows["gf"][3] == "no"
    assert rows["gf"][5] == "GIC, CARTESIAN_HESSIAN"
    assert rows["gicforge"][2] == "NEO"


def test_gui_workbench_state_reads_window_specs_and_project_sections(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0.0 0.0 0.0",
                "",
                "#VIBRATIONAL",
                "SCHEMA oracle.xyz.vibrational.v1",
                "",
                "#QFF",
                "SCHEMA oracle.xyz.qff.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    vibrational = load_workbench_gui_state(xyzin, "vibrational_spectroscopy")
    anharmonic = load_workbench_gui_state(xyzin, "anharmonic")

    assert vibrational.title == "Vibrational Spectroscopy"
    assert ("required", "VIBRATIONAL", "yes") in vibrational.sections.rows
    assert ("produced", "GF_PED", "no") in vibrational.sections.rows
    assert any(row[0] == "vpt2_vci_run" and row[3] == "yes" for row in vibrational.actions.rows)
    assert any(row[1] == "mode heat-map" for row in vibrational.exports.rows)
    assert anharmonic.status == "ready"
    assert "window: Anharmonic: VPT2 / VCI / DVR" in workbench_gui_state_lines(anharmonic)


def test_dashboard_controller_builds_ready_actions_from_xyzin_sections(tmp_path):
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
    controller = OracleDashboardController(xyzin, repo_root=ROOT)
    state = controller.state()

    actions = {action.key: action for action in controller.actions(state)}

    assert actions["validate"].enabled
    assert actions["gicforge_build"].enabled
    assert actions["gicforge_report"].enabled
    assert actions["rovib_summary"].enabled is False
    assert actions["rovib_summary"].reason == "missing ROTATIONAL"
    assert actions["gicforge_build"].command.cwd == ROOT
    assert actions["avogadro"].command.cwd is None


def test_dashboard_controller_runs_action_and_records_log(tmp_path, monkeypatch):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text("1\nh\nH 0.0 0.0 0.0\n", encoding="utf-8")
    command = OracleGuiCommand("Fake action", (sys.executable, "-c", "print('ok')"))
    controller = OracleDashboardController(
        xyzin,
        actions=(
            DashboardActionTemplate(
                "fake",
                "dashboard",
                "Fake action",
                lambda _xyzin: command,
            ),
        ),
        repo_root=ROOT,
    )

    def fake_run(argv, **kwargs):
        assert argv == command.argv
        assert kwargs["cwd"] is None
        return subprocess.CompletedProcess(argv, 0, stdout="stdout\n", stderr="stderr\n")

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = controller.run_action("fake")

    assert result.ok
    assert result.stdout == "stdout\n"
    assert controller.log_lines[0].startswith("$ ")
    assert "stdout" in controller.log_lines
    assert "stderr" in controller.log_lines
    assert controller.log_lines[-1] == "[exit 0] Fake action"


def test_structure_gui_state_reads_topology_synthons_and_fragments(tmp_path):
    xyzin = tmp_path / "dimer.xyzin"
    xyzin.write_text(
        "\n".join(
            [
                "3",
                "fragmented",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 0.9",
                "Na 3.0 0.0 0.0",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
                "INDEXING ATOMS=ONE_BASED",
                "[BONDS]",
                "1 2",
                "[RINGS]",
                "1 SIZE=3 ATOMS=1 2 3",
                "",
                "#SYNTHONS",
                "SCHEMA oracle.xyz.synthons.v1",
                "INDEXING ATOMS=ONE_BASED",
                "COLUMNS ATOM Z ZEFF CHARGE COVALENCY DELOCALIZATION STRAIN SIGNATURE",
                "1 O 7.5 -0.4 2.0 0.1 0.0 O,water",
                "2 H 1.1 0.2 1.0 0.0 0.0 H,water",
                "",
                "#FRAGMENTS",
                "SCHEMA oracle.xyz.fragments.v1",
                "STATUS BUILT",
                "INDEXING ATOMS=ONE_BASED",
                "[FRAGMENTS]",
                "F001 LABEL=component_1 SIZE=2 ATOMS=1,2",
                "[CENTERS]",
                "F001 X=0.0 Y=0.0 Z=0.45",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    state = load_structure_gui_state(xyzin)

    assert state.exists
    assert state.topology_bonds.rows == (("1", "2"),)
    assert state.topology_rings.rows == (("1", "3", "1 2 3"),)
    assert state.synthons.rows[0] == ("1", "O", "7.5", "-0.4", "2.0", "0.1", "0.0", "O,water")
    assert state.fragments.rows == (("F001", "component_1", "2", "1,2", "0.0, 0.0, 0.45"),)


def test_structure_controller_builds_babel_and_fragment_commands(tmp_path):
    source = tmp_path / "water.xyz"
    xyzin = tmp_path / "water.xyzin"
    controller = OracleStructureController(xyzin)

    preprocess = controller.preprocess_command(source, source_kind="xyz")
    fragments = controller.fragments_command("build")

    assert preprocess.argv[-2:] == ("--source-kind", "xyz")
    assert preprocess.argv[3:6] == ("babel", "preprocess", str(source))
    assert preprocess.argv[6] == str(source.with_suffix(".xyzin"))
    assert fragments.argv[-3:] == ("fragments", "build", str(xyzin))
    assert fragments.required_sections == ("TOPOLOGY", "SYNTHONS")


def test_gicforge_gui_state_reports_missing_gic_without_crashing(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text("1\nh\nH 0.0 0.0 0.0\n", encoding="utf-8")

    state = load_gicforge_gui_state(xyzin)
    lines = gicforge_gui_state_lines(state)

    assert state.exists
    assert not state.ready
    assert "missing #GIC section" in state.messages
    assert "ready: 0" in lines


def test_gicforge_gui_state_reads_frozen_definition_tables(tmp_path):
    source = tmp_path / "water.xyz"
    xyzin = tmp_path / "water.xyzin"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.000000 0.000000 0.000000",
                "H 0.000000 0.757000 0.586000",
                "H 0.000000 -0.757000 0.586000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True, sycart=True)

    state = load_gicforge_gui_state(xyzin)
    lines = gicforge_gui_state_lines(state)

    assert state.ready
    assert state.summary.point_group == definition.point_group
    assert state.summary.rank == definition.rank
    assert state.summary.primitive_count == len(definition.primitives)
    assert state.summary.gic_count == len(definition.gics)
    assert state.primitives.rows
    assert state.frozen_gics.rows
    assert state.diagnostics.rows
    assert "ready: 1" in lines
    assert any(line.startswith("total symmetric irrep:") for line in lines)


def test_gicforge_controller_builds_dedicated_commands(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    controller = OracleGICForgeController(xyzin)

    build = controller.build_command(
        symmetrize=False,
        sycart=False,
        improper_dihedrals=False,
    )
    report = controller.report_command(default_gicforge_report_output(xyzin))
    bmatrix = controller.bmatrix_command(default_gicforge_bmatrix_output(xyzin))
    gaussian = controller.gaussian_input_command(
        default_gicforge_gaussian_output(xyzin),
        route="#p b3lyp/def2svp",
    )

    assert build.argv[-1] == str(xyzin)
    assert "--symmetrize" not in build.argv
    assert "--sycart" not in build.argv
    assert "--improper-dihedrals" not in build.argv
    assert report.argv[-1].endswith(".gicforge_report.txt")
    assert bmatrix.argv[-1].endswith(".gic_bmatrix.txt")
    assert gaussian.argv[-3:] == (
        str(default_gicforge_gaussian_output(xyzin)),
        "--route",
        "#p b3lyp/def2svp",
    )


def test_gf_gui_state_reports_missing_section_without_crashing(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text("1\nh\nH 0.0 0.0 0.0\n", encoding="utf-8")

    state = load_gf_gui_state(xyzin)
    lines = gf_gui_state_lines(state)

    assert state.exists
    assert not state.ready
    assert "missing #GF_PED section" in state.messages
    assert "ready: 0" in lines


def test_gf_gui_state_reads_normalized_gf_ped_section(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    report = tmp_path / "gf.report"
    csv_dir = tmp_path / "gf_csv"
    xyzin.write_text("1\nh\nH 0.0 0.0 0.0\n", encoding="utf-8")
    write_gf_ped_section(
        xyzin,
        GFPEDSection(
            source_kind="fchk",
            source_path=tmp_path / "job.fchk",
            hessian_source="FCHK job.fchk",
            coordinate_source="xyzin-frozen-gic:C2V",
            report_path=report,
            csv_dir=csv_dir,
            status="complete",
            point_group="C2V",
            symmetrized_gics=True,
            matrix_model="LOCAL+IRREP_BLOCKS",
            hessian_correction="ELECTROSTATIC(Synthons)+UFF_VDW; 1-4 scale=0.5",
            force_threshold=1.0e-8,
            modes=(
                GFModeRow(1, 1000.0, "A1"),
                GFModeRow(2, 1500.0, "B2"),
            ),
            gics=(
                GFGICRow("GIC001", "A1Str001", "A1", "R(1,2)", (90.0, 10.0), 0.9),
                GFGICRow("GIC002", "B2Bend001", "B2", "A(1,2,3)", (5.0, 95.0)),
            ),
        ),
    )

    state = load_gf_gui_state(xyzin)
    lines = gf_gui_state_lines(state)

    assert state.ready
    assert state.summary.point_group == "C2V"
    assert state.summary.mode_count == 2
    assert state.summary.gic_count == 2
    assert state.frequencies.rows[0] == ("1", "1000", "A1")
    assert state.gics.rows[0][:6] == ("GIC001", "A1Str001", "A1", "0.9", "1", "90")
    assert state.ped.columns == ("GIC", "Name", "Irrep", "M01", "M02")
    assert state.ped.rows[1] == ("GIC002", "B2Bend001", "B2", "5", "95")
    assert "matrix model: LOCAL+IRREP_BLOCKS" in lines
    assert f"report: {report}" in lines


def test_gf_controller_builds_run_command_with_options(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    fchk = tmp_path / "job.fchk"
    scale = tmp_path / "scale.txt"
    controller = OracleGFController(xyzin)

    command = controller.run_command(
        fchk=fchk,
        out=default_gf_report_output(xyzin),
        csv_dir=default_gf_csv_dir(xyzin),
        scale_file=scale,
        scale_records=("GIC003=0.9",),
        local=True,
        symmetry_blocks=False,
        force_threshold=1.0e-8,
        subtract_electrostatic=True,
        subtract_uff_vdw=True,
        nonbonded_14_scale=0.25,
        write_section=False,
    )

    assert command.argv[:4] == (sys.executable, "-m", "oracle", "gf")
    assert ("--xyzin", str(xyzin)) == command.argv[4:6]
    assert "--fchk" in command.argv
    assert "--local" in command.argv
    assert "--symmetry-blocks" not in command.argv
    assert command.argv[command.argv.index("--scale") + 1] == "GIC003=0.9"
    assert command.argv[command.argv.index("--force-threshold") + 1] == "1e-08"
    assert "--subtract-electrostatic" in command.argv
    assert "--subtract-uff-vdw" in command.argv
    assert command.argv[command.argv.index("--nonbonded-14-scale") + 1] == "0.25"
    assert "--no-write-section" in command.argv


def test_sefit_gui_state_reports_missing_morpheus_without_crashing(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text("1\nh\nH 0.0 0.0 0.0\n", encoding="utf-8")

    state = load_sefit_gui_state(xyzin)
    lines = sefit_gui_state_lines(state)

    assert state.exists
    assert not state.ready
    assert "missing #MORPHEUS section" in state.messages
    assert "ready: 0" in lines


def test_sefit_gui_state_reads_morpheus_and_isotopologues(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    outdir = tmp_path / "semiexp"
    xyzin.write_text("1\nh\nH 0.0 0.0 0.0\n", encoding="utf-8")
    write_xyzin_isotopologue_records(
        xyzin,
        (
            XyzinIsotopologueRecord(
                label="parent",
                rotational_MHz=(1000.0, 900.0, 800.0),
                sigma_MHz=(0.01, 0.01, 0.01),
            ),
            XyzinIsotopologueRecord(
                label="d1",
                substitutions={1: 2},
                rotational_MHz=(950.0, 850.0, 750.0),
            ),
        ),
    )
    write_morpheus_section(
        xyzin,
        MorpheusSection(
            status="complete",
            xyzin_path=xyzin,
            run_dir=outdir,
            manifest_path=outdir / "semiexp_manifest.json",
            text_report_path=outdir / "semiexp_report.txt",
            html_report_path=outdir / "semiexp_report.html",
            latex_tables_path=outdir / "semiexp_tables.tex",
            geometry_path=outdir / "semiexp_geometry.xyz",
            parameters_path=outdir / "semiexp_parameters.csv",
            residuals_path=outdir / "semiexp_residuals.csv",
            rotational_constants_path=outdir / "semiexp_rotational_constants.csv",
            diagnostics_path=outdir / "semiexp_diagnostics.csv",
            backend="python",
            coordinate_model="gic",
            observable="moments",
            components=("A", "B", "C"),
            rms_MHz=0.02,
            rotational_rms_MHz=0.03,
            iterations=2,
            stationary_point="minimum",
            convergence="max_iter",
            rank=3,
            condition_number=123.0,
            isotopologue_count=2,
            parameter_count=3,
            active_parameter_count=2,
            warning_count=1,
        ),
    )

    state = load_sefit_gui_state(xyzin)
    lines = sefit_gui_state_lines(state)

    assert state.ready
    assert state.summary.run_dir == outdir
    assert state.summary.coordinate_model == "gic"
    assert state.summary.components == ("A", "B", "C")
    assert state.isotopologues.rows[0] == (
        "parent",
        "parent",
        "1000,900,800",
        "",
        "",
        "0.01,0.01,0.01",
    )
    assert state.isotopologues.rows[1][1] == "1:2"
    assert ("Run dir", str(outdir)) in state.outputs.rows
    assert ("Condition number", "123") in state.diagnostics.rows
    assert "coordinate model: gic" in lines
    assert "warnings: 1" in lines


def test_sefit_controller_builds_semiexp_command(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    job = tmp_path / "job.toml"
    controller = OracleSEFitController(xyzin)

    command = controller.run_command(
        job=job,
        outdir=default_sefit_outdir(xyzin),
        backend="python",
        write_section=False,
        extra_args=("--max-iter", "2"),
    )

    assert command.argv[:4] == (sys.executable, "-m", "oracle", "semiexp")
    assert command.argv[command.argv.index("--job") + 1] == str(job)
    assert command.argv[command.argv.index("--xyzin") + 1] == str(xyzin)
    assert command.argv[command.argv.index("--outdir") + 1] == str(default_sefit_outdir(xyzin))
    assert "--no-write-section" in command.argv
    assert command.argv[-2:] == ("--max-iter", "2")


def test_trinity_gui_state_reports_missing_section_without_crashing(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    xyzin.write_text("1\nh\nH 0.0 0.0 0.0\n", encoding="utf-8")

    state = load_trinity_gui_state(xyzin)
    lines = trinity_gui_state_lines(state)

    assert state.exists
    assert not state.ready
    assert "missing #TRINITY section" in state.messages
    assert "ready: 0" in lines


def test_trinity_gui_state_reads_prepared_section(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    run_dir = tmp_path / "trinity"
    xyzin.write_text("1\nh\nH 0.0 0.0 0.0\n", encoding="utf-8")
    write_trinity_section(
        xyzin,
        TrinitySection(
            source_path=xyzin,
            run_dir=run_dir,
            manifest_path=run_dir / "trinity_manifest.json",
            engine_command="engine --gradient",
            coordinate_model="gic",
            active_space="total_symmetric",
            max_steps=20,
            trust_radius=0.25,
            trajectory_path=run_dir / "trinity_trajectory.xyz",
            final_geometry_path=run_dir / "trinity_final.xyz",
            energy_gradient_log_path=run_dir / "trinity_energy_gradient.jsonl",
        ),
    )

    state = load_trinity_gui_state(xyzin)
    lines = trinity_gui_state_lines(state)

    assert state.ready
    assert state.summary.run_dir == run_dir
    assert state.summary.engine_command == "engine --gradient"
    assert ("Coordinate model", "gic") in state.settings.rows
    assert ("Final geometry", str(run_dir / "trinity_final.xyz")) in state.outputs.rows
    assert "active space: total_symmetric" in lines


def test_trinity_controller_builds_prepare_command(tmp_path):
    xyzin = tmp_path / "molecule.xyzin"
    controller = OracleTrinityController(xyzin)

    command = controller.prepare_command(
        run_dir=default_trinity_run_dir(xyzin),
        engine_command="engine --gradient",
        coordinate_model="cartesian",
        active_space="cartesian",
        max_steps=5,
        trust_radius=0.1,
    )

    assert command.argv[:5] == (sys.executable, "-m", "oracle", "trinity", "prepare")
    assert command.argv[5] == str(xyzin)
    assert command.argv[command.argv.index("--run-dir") + 1] == str(default_trinity_run_dir(xyzin))
    assert command.argv[command.argv.index("--engine-command") + 1] == "engine --gradient"
    assert command.argv[command.argv.index("--coordinate-model") + 1] == "cartesian"
    assert command.argv[command.argv.index("--active-space") + 1] == "cartesian"
    assert command.produced_sections == ("TRINITY",)


def test_gui_window_specs_cover_primary_oracle_workflows():
    keys = {spec.key for spec in ORACLE_GUI_WINDOWS}

    assert {
        "dashboard",
        "babel",
        "topology",
        "gicforge",
        "gf",
        "sefit",
        "trinity",
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
    assert window_spec("trinity").produced_sections == ("TRINITY",)
    assert "GIC" in all_known_sections()
    assert "KINETICS" in all_known_sections()
    assert "TRINITY" in all_known_sections()


def test_gui_command_specs_use_oracle_cli_without_shell(tmp_path):
    source = tmp_path / "mol.xyz"
    xyzin = tmp_path / "mol.xyzin"
    fchk = tmp_path / "mol.fchk"
    log = tmp_path / "mol.log"
    job = tmp_path / "job.toml"

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
    gf = gf_command(
        xyzin,
        fchk=fchk,
        out=tmp_path / "gf.report",
        csv_dir=tmp_path / "gf_csv",
        scale_records=("GIC001=0.95",),
        local=True,
        subtract_electrostatic=True,
        subtract_uff_vdw=True,
        nonbonded_14_scale=0.25,
    )
    sefit = semiexp_command(job, tmp_path / "semiexp", xyzin=xyzin, extra_args=("--max-iter", "2"))
    trinity = trinity_prepare_command(
        xyzin,
        run_dir=tmp_path / "trinity",
        engine_command="engine --gradient",
        max_steps=3,
    )

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
    assert gf.argv[3:6] == ("gf", "--xyzin", str(xyzin))
    assert "--local" in gf.argv
    assert "--subtract-electrostatic" in gf.argv
    assert "--subtract-uff-vdw" in gf.argv
    assert gf.argv[gf.argv.index("--nonbonded-14-scale") + 1] == "0.25"
    assert gf.produced_sections == ("GF_PED",)
    assert sefit.argv[3:5] == ("semiexp", "--job")
    assert sefit.argv[sefit.argv.index("--xyzin") + 1] == str(xyzin)
    assert sefit.argv[-2:] == ("--max-iter", "2")
    assert sefit.produced_sections == ("MORPHEUS",)
    assert trinity.argv[3:6] == ("trinity", "prepare", str(xyzin))
    assert trinity.argv[trinity.argv.index("--engine-command") + 1] == "engine --gradient"
    assert trinity.argv[trinity.argv.index("--max-steps") + 1] == "3"
    assert trinity.produced_sections == ("TRINITY",)
