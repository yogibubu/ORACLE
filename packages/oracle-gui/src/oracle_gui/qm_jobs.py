from __future__ import annotations

from pathlib import Path

from .commands import (
    OracleGuiCommand,
    gaussian_fchk_summary_command,
    gaussian_formchk_command,
    gaussian_promote_fchk_command,
    gaussian_promote_rovib_command,
    gaussian_run_command,
    gaussian_status_command,
    gaussian_summary_command,
    gicforge_gaussian_input_command,
    molpro_promote_command,
    molpro_summary_command,
    mrcc_promote_command,
    mrcc_summary_command,
)


class OracleQMJobsController:
    def __init__(self, xyzin: Path | str | None = None) -> None:
        self.xyzin = None if xyzin is None else Path(xyzin)

    def set_xyzin(self, xyzin: Path | str | None) -> None:
        self.xyzin = None if xyzin is None else Path(xyzin)

    def gaussian_input_command(
        self,
        output: Path | str,
        *,
        route: str,
        title: str | None = None,
    ) -> OracleGuiCommand:
        if self.xyzin is None:
            raise ValueError("no ORACLE xyzin project is loaded")
        return gicforge_gaussian_input_command(self.xyzin, output, route=route, title=title)

    def gaussian_status_command(self, workdir: Path | str) -> OracleGuiCommand:
        return gaussian_status_command(workdir)

    def gaussian_run_command(
        self,
        workdir: Path | str,
        *,
        executable: str | None = None,
        input_path: Path | str | None = None,
        background: bool = False,
        timeout: float | None = None,
    ) -> OracleGuiCommand:
        return gaussian_run_command(
            workdir,
            executable=executable,
            input_path=input_path,
            background=background,
            timeout=timeout,
        )

    def gaussian_formchk_command(
        self,
        chk: Path | str,
        fchk: Path | str | None = None,
        *,
        executable: str | None = None,
        timeout: float | None = None,
    ) -> OracleGuiCommand:
        return gaussian_formchk_command(chk, fchk, executable=executable, timeout=timeout)

    def gaussian_fchk_summary_command(self, fchk: Path | str) -> OracleGuiCommand:
        return gaussian_fchk_summary_command(fchk)

    def gaussian_log_summary_command(self, log: Path | str) -> OracleGuiCommand:
        return gaussian_summary_command(log)

    def gaussian_promote_fchk_command(
        self,
        fchk: Path | str,
        *,
        cartesian_hessian: bool = True,
        normal_modes: bool = True,
        qff: bool = True,
    ) -> OracleGuiCommand:
        if self.xyzin is None:
            raise ValueError("no ORACLE xyzin project is loaded")
        return gaussian_promote_fchk_command(
            fchk,
            self.xyzin,
            cartesian_hessian=cartesian_hessian,
            normal_modes=normal_modes,
            qff=qff,
        )

    def gaussian_promote_rovib_command(
        self,
        log: Path | str,
        *,
        vibrational: bool = True,
        rotational: bool = True,
        deltabvib: bool = True,
        invert_imaginary: bool = True,
        exclude_modes: tuple[int, ...] = (),
    ) -> OracleGuiCommand:
        if self.xyzin is None:
            raise ValueError("no ORACLE xyzin project is loaded")
        return gaussian_promote_rovib_command(
            log,
            self.xyzin,
            vibrational=vibrational,
            rotational=rotational,
            deltabvib=deltabvib,
            invert_imaginary=invert_imaginary,
            exclude_modes=exclude_modes,
        )

    def molpro_summary_command(self, output: Path | str) -> OracleGuiCommand:
        return molpro_summary_command(output)

    def molpro_promote_command(self, output: Path | str) -> OracleGuiCommand:
        if self.xyzin is None:
            raise ValueError("no ORACLE xyzin project is loaded")
        return molpro_promote_command(output, self.xyzin)

    def mrcc_summary_command(self, output: Path | str) -> OracleGuiCommand:
        return mrcc_summary_command(output)

    def mrcc_promote_command(self, output: Path | str) -> OracleGuiCommand:
        if self.xyzin is None:
            raise ValueError("no ORACLE xyzin project is loaded")
        return mrcc_promote_command(output, self.xyzin)


def default_qm_gaussian_input_output(xyzin: Path | str) -> Path:
    target = Path(xyzin)
    return target.with_name(f"{target.stem}.gic.gjf")


def default_qm_gaussian_workdir(xyzin: Path | str) -> Path:
    target = Path(xyzin)
    return target.with_name(f"{target.stem}.gaussian")


def default_qm_formchk_output(chk: Path | str) -> Path:
    target = Path(chk)
    return target.with_suffix(".fchk")
