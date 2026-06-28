from __future__ import annotations

import os
from pathlib import Path
import json
import subprocess

import pytest

from matrix_engines import (
    RemoteQMMetadata,
    remote_qm_fetch,
    remote_qm_submit,
)
from matrix_core import build_run_manifest, read_basic_section, read_sectioned_lines, section_content
from matrix_core.cli import main
from matrix_molpro import molpro_job_status, run_molpro_job
from matrix_oracle.qm_jobs import OracleQMJobsController
from matrix_orca import orca_job_status, run_orca_job


ROOT = Path(__file__).resolve().parents[1]
QM_LAUNCHERS = ROOT / "tests" / "fixtures" / "qm_launchers"


def test_qm_launcher_golden_manifest_covers_all_engines(tmp_path):
    manifest = json.loads((QM_LAUNCHERS / "golden_manifest.json").read_text(encoding="utf-8"))
    engines = {case["engine"] for case in manifest["cases"]}

    assert manifest["schema"] == "matrix.qm_launcher_golden.v1"
    assert engines == {"gdv32", "g16", "molpro", "orca"}

    for case in manifest["cases"]:
        input_path = QM_LAUNCHERS / case["input"]
        output_path = tmp_path / f"{case['engine']}{case['native_output_suffix']}"
        output_path.write_text("golden output placeholder\n", encoding="utf-8")
        run_manifest = build_run_manifest(
            workflow=f"qm.launcher.{case['engine']}",
            status="golden",
            run_dir=tmp_path / case["engine"],
            inputs={"input": input_path},
            outputs={"native_output": output_path},
            parameters={"promote_modes": case["promote_modes"]},
            backend={"name": case["engine"], "executable": case["engine"]},
            command=case["local_command"],
        )

        assert input_path.is_file()
        assert run_manifest.framework == "MATRIX"
        assert run_manifest.inputs["input"].endswith(case["input"])
        assert run_manifest.outputs["native_output"].endswith(case["native_output_suffix"])


@pytest.mark.parametrize("case_name", ("gdv32_h2.gjf", "g16_h2.gjf", "molpro_h2.com", "orca_h2.inp"))
def test_qm_launcher_golden_inputs_can_be_remote_submitted(case_name):
    source = QM_LAUNCHERS / case_name
    runner = _FakeRemoteRunner(QM_LAUNCHERS)
    engine = case_name.split("_", maxsplit=1)[0]

    result = remote_qm_submit(source, engine=engine, host="enzo@oracle", runner=runner)

    assert result.engine == engine
    assert result.job == f"job-{engine}"
    assert result.native_output.endswith(".log" if engine in {"gdv32", "g16"} else ".out")
    assert any(call[0] == "scp" and call[1] == str(source.resolve()) for call in runner.calls)


def test_molpro_launcher_uses_shared_external_qm_runner(tmp_path):
    (tmp_path / "molpro.com").write_text("***,water\ngeometry={H}\n", encoding="utf-8")
    executable = _fake_executable(
        tmp_path / "fake_molpro.sh",
        "printf 'Molpro calculation terminated\\nVariable memory released\\n'\n",
    )

    result = run_molpro_job(tmp_path, executable=str(executable))
    status = molpro_job_status(tmp_path)

    assert result.success is True
    assert result.exit_code == 0
    assert result.output_path == tmp_path / "molpro.out"
    assert status.status == "completed"
    assert status.normal_termination is True
    assert "Variable memory released" in result.output_path.read_text(encoding="utf-8")


def test_orca_launcher_uses_shared_external_qm_runner(tmp_path):
    (tmp_path / "orca.inp").write_text("! HF STO-3G\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
    executable = _fake_executable(
        tmp_path / "fake_orca.sh",
        "printf '****ORCA TERMINATED NORMALLY****\\n'\n",
    )

    result = run_orca_job(tmp_path, executable=str(executable))
    status = orca_job_status(tmp_path)

    assert result.success is True
    assert result.exit_code == 0
    assert result.output_path == tmp_path / "orca.out"
    assert status.status == "completed"
    assert status.normal_termination is True


def test_molpro_and_orca_cli_run_commands(tmp_path, capsys):
    molpro_dir = tmp_path / "molpro"
    molpro_dir.mkdir()
    (molpro_dir / "molpro.com").write_text("***,h\ngeometry={H}\n", encoding="utf-8")
    molpro_exe = _fake_executable(
        tmp_path / "fake_molpro.sh",
        "printf 'Molpro calculation terminated\\n'\n",
    )

    assert main(["molpro", "run", str(molpro_dir), "--executable", str(molpro_exe)]) == 0
    molpro_output = capsys.readouterr().out
    assert "program: Molpro" in molpro_output
    assert "success: 1" in molpro_output

    orca_dir = tmp_path / "orca"
    orca_dir.mkdir()
    (orca_dir / "orca.inp").write_text("! SP\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
    orca_exe = _fake_executable(
        tmp_path / "fake_orca.sh",
        "printf 'ORCA TERMINATED NORMALLY\\n'\n",
    )

    assert main(["orca", "run", str(orca_dir), "--executable", str(orca_exe)]) == 0
    orca_output = capsys.readouterr().out
    assert "program: ORCA" in orca_output
    assert "success: 1" in orca_output


def test_qm_jobs_controller_exposes_molpro_and_orca_launchers(tmp_path):
    controller = OracleQMJobsController(tmp_path / "molecule.xyzin")

    molpro_run = controller.molpro_run_command(
        tmp_path / "molpro",
        executable="molpro",
        input_path="input.com",
        background=True,
    )
    orca_status = controller.orca_status_command(tmp_path / "orca", output_path="orca.out")

    assert molpro_run.argv[1:4] == ("-m", "matrix", "molpro")
    assert "run" in molpro_run.argv
    assert "--background" in molpro_run.argv
    assert orca_status.argv[1:4] == ("-m", "matrix", "orca")
    assert "status" in orca_status.argv
    assert "--output" in orca_status.argv


def test_remote_metadata_parser_uses_native_output_when_present():
    metadata = RemoteQMMetadata.from_text(
        "job=job1\n"
        "engine=molpro\n"
        "pid=123\n"
        "workdir=/home/enzo/matrix/jobs/job1\n"
        "log=/home/enzo/matrix/logs/job1.log\n"
        "native_output=/home/enzo/matrix/jobs/job1/molecule.out\n"
        "input=/home/enzo/matrix/jobs/job1/molecule.com\n"
    )

    assert metadata.job == "job1"
    assert metadata.engine == "molpro"
    assert metadata.native_output.endswith("molecule.out")


def test_remote_submit_uses_ssh_scp_and_parses_submit_output(tmp_path):
    source = tmp_path / "input.gjf"
    source.write_text("#p hf/sto-3g sp\n", encoding="utf-8")
    runner = _FakeRemoteRunner(tmp_path)

    result = remote_qm_submit(
        source,
        engine="gdv32",
        host="enzo@oracle",
        extra_args=("--flag",),
        runner=runner,
    )

    assert result.job == "job-gdv32"
    assert result.remote_input == "~/matrix/inputs/input.gjf"
    assert any("matrix-submit" in " ".join(call) for call in runner.calls)
    assert any(call[0] == "scp" for call in runner.calls)


def test_remote_fetch_writes_output_metadata_and_manifest(tmp_path):
    runner = _FakeRemoteRunner(tmp_path)

    result = remote_qm_fetch(
        "job-molpro",
        host="enzo@oracle",
        destination=tmp_path / "fetched",
        promote="none",
        runner=runner,
    )

    assert result.output_path.read_text(encoding="utf-8") == "Molpro calculation terminated\n"
    assert result.metadata_path.exists()
    assert result.manifest_path.exists()
    assert "matrix.remote_qm.fetch.v1" in result.manifest_path.read_text(encoding="utf-8")


def test_remote_fetch_can_promote_orca_output(tmp_path):
    runner = _FakeRemoteRunner(tmp_path)
    xyzin = tmp_path / "molecule.xyzin"

    result = remote_qm_fetch(
        "job-orca",
        host="enzo@oracle",
        destination=tmp_path / "fetched",
        promote="auto",
        xyzin=xyzin,
        runner=runner,
    )
    basic = read_basic_section(xyzin)
    lines = read_sectioned_lines(xyzin)

    assert result.engine == "orca"
    assert result.promotion is not None
    assert result.promotion.mode == "orca"
    assert result.promotion.status == "promoted"
    assert basic.charge == 0
    assert section_content(lines, "SOURCE")[2] == "FORMAT orca_output"


def test_qm_jobs_controller_exposes_remote_qm_commands(tmp_path):
    controller = OracleQMJobsController(tmp_path / "molecule.xyzin")

    submit = controller.remote_submit_command(tmp_path / "job.gjf", engine="gdv32")
    status = controller.remote_status_command(host="enzo@oracle")
    fetch = controller.remote_fetch_command("job-gdv32", promote="gaussian-log-hessian")

    assert submit.argv[1:5] == ("-m", "matrix", "qm", "remote-submit")
    assert "--engine" in submit.argv
    assert status.argv[1:5] == ("-m", "matrix", "qm", "remote-status")
    assert fetch.argv[1:5] == ("-m", "matrix", "qm", "remote-fetch")
    assert "--xyzin" in fetch.argv


def test_qm_jobs_controller_exposes_orca_promote(tmp_path):
    controller = OracleQMJobsController(tmp_path / "molecule.xyzin")

    summary = controller.orca_summary_command(tmp_path / "orca.out")
    promote = controller.orca_promote_command(tmp_path / "orca.out")
    fetch = controller.remote_fetch_command("job-orca", promote="orca")

    assert summary.argv[1:4] == ("-m", "matrix", "orca")
    assert "summary" in summary.argv
    assert promote.argv[1:4] == ("-m", "matrix", "orca")
    assert "promote" in promote.argv
    assert "CARTESIAN_HESSIAN" in promote.produced_sections
    assert "--promote" in fetch.argv
    assert "orca" in fetch.argv


@pytest.mark.skipif(
    os.environ.get("MATRIX_ORACLE_REMOTE_SMOKE") != "1",
    reason="remote oracle smoke test is opt-in",
)
def test_remote_oracle_status_smoke_opt_in():
    ssh = os.environ.get("MATRIX_ORACLE_SSH", "ssh")
    host = os.environ.get("MATRIX_ORACLE_HOST", "enzo@oracle")
    assert main(["qm", "remote-status", "--host", host, "--ssh", ssh]) == 0


def _fake_executable(path: Path, body: str) -> Path:
    path.write_text(
        "#!/bin/sh\n"
        'test -f "$1" || exit 9\n'
        f"{body}"
        "exit 0\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


class _FakeRemoteRunner:
    def __init__(self, tmp_path: Path) -> None:
        self.tmp_path = tmp_path
        self.calls: list[list[str]] = []

    def __call__(
        self,
        argv,
        *,
        text: bool,
        capture_output: bool,
        check: bool,
    ) -> subprocess.CompletedProcess[str]:
        del text, capture_output, check
        args = [str(item) for item in argv]
        self.calls.append(args)
        if args[0] == "ssh" and "matrix-submit" in args[-1]:
            engine = next(
                (item for item in ("gdv32", "g16", "molpro", "orca") if f" {item} " in args[-1]),
                "gdv32",
            )
            suffix = ".log" if engine in {"gdv32", "g16"} else ".out"
            return subprocess.CompletedProcess(
                args,
                0,
                stdout=(
                    f"Submitted job-{engine}\n"
                    "PID: 123\n"
                    f"Workdir: /home/enzo/matrix/jobs/job-{engine}\n"
                    f"Log: /home/enzo/matrix/logs/job-{engine}.log\n"
                    f"Native output: /home/enzo/matrix/jobs/job-{engine}/input{suffix}\n"
                ),
                stderr="",
            )
        if args[0] == "ssh" and "metadata.txt" in args[-1]:
            if "job-orca" in args[-1]:
                return subprocess.CompletedProcess(
                    args,
                    0,
                    stdout=(
                        "job=job-orca\n"
                        "engine=orca\n"
                        "pid=123\n"
                        "workdir=/home/enzo/matrix/jobs/job-orca\n"
                        "log=/home/enzo/matrix/logs/job-orca.log\n"
                        "native_output=/home/enzo/matrix/jobs/job-orca/orca.out\n"
                        "input=/home/enzo/matrix/jobs/job-orca/orca.inp\n"
                    ),
                    stderr="",
                )
            return subprocess.CompletedProcess(
                args,
                0,
                stdout=(
                    "job=job-molpro\n"
                    "engine=molpro\n"
                    "pid=123\n"
                    "workdir=/home/enzo/matrix/jobs/job-molpro\n"
                    "log=/home/enzo/matrix/logs/job-molpro.log\n"
                    "native_output=/home/enzo/matrix/jobs/job-molpro/molecule.out\n"
                    "input=/home/enzo/matrix/jobs/job-molpro/molecule.com\n"
                ),
                stderr="",
            )
        if args[0] == "ssh":
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")
        if args[0] == "scp" and args[1].endswith("molecule.out"):
            Path(args[2]).write_text("Molpro calculation terminated\n", encoding="utf-8")
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")
        if args[0] == "scp" and args[1].endswith("orca.out"):
            Path(args[2]).write_text(
                "\n".join(
                    [
                        "Program ORCA",
                        "* xyz 0 1",
                        "Total Charge           Charge    ....    0",
                        "Multiplicity           Mult      ....    1",
                        "FINAL SINGLE POINT ENERGY     -1.000000000000",
                        "CARTESIAN COORDINATES (ANGSTROEM)",
                        "---------------------------------",
                        "H       0.00000000     0.00000000     0.00000000",
                        "H       0.00000000     0.00000000     0.74000000",
                        "",
                        "****ORCA TERMINATED NORMALLY****",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")
        if args[0] == "scp":
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")
        return subprocess.CompletedProcess(args, 1, stdout="", stderr="unexpected command")
