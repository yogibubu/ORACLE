#!/usr/bin/env python3
"""Small Tkinter front-end for the mass-weighted path DVR command line tool."""

from __future__ import annotations

import shlex
import subprocess
import sys
import threading
from pathlib import Path
from tkinter import (
    BooleanVar,
    Button,
    Checkbutton,
    Entry,
    Frame,
    Label,
    LabelFrame,
    OptionMenu,
    StringVar,
    Tk,
    Text,
    filedialog,
    messagebox,
)
from tkinter.scrolledtext import ScrolledText


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "mw_path_dvr.py"


class PathDVRGui:
    def __init__(self) -> None:
        self.root = Tk()
        self.root.title("MW Path DVR")
        self.root.geometry("1040x760")

        self.workflow = StringVar(value="Anharmonic Gaussian mode")
        self.source_path = StringVar()
        self.properties_csv = StringVar()
        self.property_derivatives_csv = StringVar()
        self.vpt2_property_csv = StringVar()
        self.outdir = StringVar(value=str(ROOT / "gaussian_outputs"))
        self.figdir = StringVar(value=str(ROOT / "figs"))
        self.prefix = StringVar(value="mw_path_dvr")
        self.grid = StringVar(value="401")
        self.levels = StringVar(value="12")
        self.save_states = StringVar(value="6")
        self.plot_max_state = StringVar(value="5")
        self.plot_property_smooth_degree = StringVar(value="0")
        self.property_vpt2 = StringVar(value="auto")
        self.property_vpt2_fit_points = StringVar(value="7")
        self.property_vpt2_degree = StringVar(value="4")
        self.property_vpt2_basis_size = StringVar(value="24")
        self.solver = StringVar(value="auto")
        self.boundary = StringVar(value="nonperiodic")
        self.log_selection = StringVar(value="all")
        self.gaussian_energy = StringVar(value="auto")
        self.energy_key = StringVar()
        self.energy_unit = StringVar(value="auto")
        self.property_name = StringVar()
        self.temperature = StringVar()
        self.compute_rotconst = BooleanVar(value=False)
        self.path_symmetry = StringVar(value="none")
        self.well_type = StringVar(value="auto")
        self.potential_extension = StringVar(value="none")
        self.extension_length = StringVar(value="0")
        self.extension_points = StringVar(value="24")
        self.extension_degree = StringVar(value="6")
        self.extension_target = StringVar(value="10000")
        self.potential_smoothing = StringVar(value="none")
        self.potential_spline_smoothing = StringVar(value="0.0")
        self.anharmonic_mode = StringVar(value="1")
        self.anharmonic_mode_order = StringVar(value="frequency-ascending")
        self.anharmonic_model = StringVar(value="auto")
        self.anharmonic_width = StringVar(value="6")
        self.anharmonic_cubic_threshold = StringVar(value="1e-8")
        self.anharmonic_handy_beta = StringVar()
        self.anharmonic_wall = StringVar(value="10000")
        self.anharmonic_wall_degree = StringVar(value="8")
        self.anharmonic_kinetic_frequency = StringVar()
        self.q1_key = StringVar()
        self.q2_key = StringVar()
        self.boundary1 = StringVar(value="nonperiodic")
        self.boundary2 = StringVar(value="nonperiodic")
        self.basis1 = StringVar(value="12")
        self.basis2 = StringVar(value="12")
        self.metric_mode = StringVar(value="constant")
        self.g11 = StringVar(value="1.0")
        self.g22 = StringVar(value="1.0")
        self.g12 = StringVar(value="0.0")
        self.metric_xyz = StringVar()

        self.command_box: Text
        self.output_box: ScrolledText
        self.run_button: Button

        self.build()
        self.refresh_command()

    def build(self) -> None:
        root = self.root
        top = Frame(root)
        top.pack(fill="x", padx=10, pady=8)
        Label(top, text="Workflow").grid(row=0, column=0, sticky="w")
        OptionMenu(
            top,
            self.workflow,
            "Anharmonic Gaussian mode",
            "1D path",
            "2D grid",
            command=lambda _value: self.refresh_command(),
        ).grid(row=0, column=1, sticky="ew", padx=6)
        Button(top, text="Refresh command", command=self.refresh_command).grid(row=0, column=2, padx=6)
        self.run_button = Button(top, text="Run", command=self.run_command)
        self.run_button.grid(row=0, column=3, padx=6)
        top.columnconfigure(1, weight=1)

        source = LabelFrame(root, text="Input")
        source.pack(fill="x", padx=10, pady=4)
        self.path_row(source, 0, "Source", self.source_path)
        self.path_row(source, 1, "Properties CSV", self.properties_csv, optional=True)
        self.path_row(source, 2, "Property derivatives CSV", self.property_derivatives_csv, optional=True)
        self.path_row(source, 3, "VPT2 property CSV", self.vpt2_property_csv, optional=True)
        self.option_row(source, 4, "Log selection", self.log_selection, ["all", "last", "last-per-link"])
        self.option_row(source, 5, "Gaussian energy", self.gaussian_energy, ["auto", "scf", "post-scf"])
        source.columnconfigure(1, weight=1)

        common = LabelFrame(root, text="Common")
        common.pack(fill="x", padx=10, pady=4)
        self.entry_row(common, 0, "Prefix", self.prefix)
        self.path_row(common, 1, "Outdir", self.outdir, directory=True)
        self.path_row(common, 2, "Figdir", self.figdir, directory=True)
        self.entry_row(common, 3, "Grid", self.grid)
        self.entry_row(common, 4, "Levels", self.levels)
        self.entry_row(common, 5, "Save states", self.save_states)
        self.entry_row(common, 6, "Plot max state", self.plot_max_state)
        self.entry_row(common, 7, "Property plot smooth degree", self.plot_property_smooth_degree)
        self.option_row(common, 8, "Property VPT2", self.property_vpt2, ["auto", "csv-only", "off"])
        self.entry_row(common, 9, "Property VPT2 fit points", self.property_vpt2_fit_points)
        self.option_row(common, 10, "Property VPT2 degree", self.property_vpt2_degree, ["2", "3", "4"])
        self.entry_row(common, 11, "Property VPT2 basis size", self.property_vpt2_basis_size)
        self.option_row(common, 12, "Solver", self.solver, ["auto", "fourier", "gaussian", "sinc-dvr"])
        common.columnconfigure(1, weight=1)

        panels = Frame(root)
        panels.pack(fill="both", expand=True, padx=10, pady=4)
        panels.columnconfigure(0, weight=1)
        panels.columnconfigure(1, weight=1)

        one_d = LabelFrame(panels, text="1D Path")
        one_d.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        self.option_row(one_d, 0, "Boundary", self.boundary, ["nonperiodic", "periodic"])
        self.entry_row(one_d, 1, "Energy key", self.energy_key)
        self.option_row(one_d, 2, "Energy unit", self.energy_unit, ["auto", "hartree", "cm-1", "kjmol", "kcalmol"])
        self.entry_row(one_d, 3, "Property", self.property_name)
        self.entry_row(one_d, 4, "Temperature K", self.temperature)
        Checkbutton(one_d, text="Compute rotational constants", variable=self.compute_rotconst, command=self.refresh_command).grid(
            row=5, column=1, sticky="w"
        )
        self.option_row(one_d, 6, "Path symmetry", self.path_symmetry, ["none", "half-even", "half-even-origin", "half-even-last", "center-even"])
        self.option_row(one_d, 7, "Well type", self.well_type, ["auto", "single", "double"])
        self.option_row(one_d, 8, "Extension", self.potential_extension, ["none", "repulsive-polynomial", "repulsive-exponential", "morse-polynomial", "single-morse", "single-inverse-power"])
        self.entry_row(one_d, 9, "Extension length au", self.extension_length)
        self.entry_row(one_d, 10, "Extension points", self.extension_points)
        self.option_row(one_d, 11, "Extension degree", self.extension_degree, ["6", "8"])
        self.entry_row(one_d, 12, "Extension target cm-1", self.extension_target)
        self.option_row(one_d, 13, "Potential smoothing", self.potential_smoothing, ["none", "spline"])
        self.entry_row(one_d, 14, "Spline smoothing", self.potential_spline_smoothing)

        anh = LabelFrame(panels, text="Anharmonic Mode")
        anh.grid(row=0, column=1, sticky="nsew", padx=(5, 0))
        self.entry_row(anh, 0, "Mode number", self.anharmonic_mode)
        self.option_row(anh, 1, "Mode order", self.anharmonic_mode_order, ["frequency-ascending", "force-table"])
        self.option_row(anh, 2, "Well type", self.well_type, ["auto", "single", "double"])
        self.option_row(
            anh,
            3,
            "Derivative model",
            self.anharmonic_model,
            ["auto", "handy-morse", "handy-gaussian", "gaussian", "morse", "inverse-power", "taylor"],
        )
        self.entry_row(anh, 4, "Grid half-width", self.anharmonic_width)
        self.entry_row(anh, 5, "F3 zero threshold", self.anharmonic_cubic_threshold)
        self.entry_row(anh, 6, "Handy beta", self.anharmonic_handy_beta)
        self.entry_row(anh, 7, "Kinetic freq cm-1", self.anharmonic_kinetic_frequency)
        self.entry_row(anh, 8, "Wall height cm-1", self.anharmonic_wall)
        self.option_row(anh, 9, "Wall degree", self.anharmonic_wall_degree, ["6", "8"])

        grid2d = LabelFrame(root, text="2D Grid")
        grid2d.pack(fill="x", padx=10, pady=4)
        self.entry_row(grid2d, 0, "q1 key", self.q1_key)
        self.entry_row(grid2d, 1, "q2 key", self.q2_key)
        self.option_row(grid2d, 2, "Boundary 1", self.boundary1, ["nonperiodic", "periodic"])
        self.option_row(grid2d, 3, "Boundary 2", self.boundary2, ["nonperiodic", "periodic"])
        self.entry_row(grid2d, 4, "Basis 1", self.basis1)
        self.entry_row(grid2d, 5, "Basis 2", self.basis2)
        self.option_row(grid2d, 6, "Metric", self.metric_mode, ["constant", "csv", "geometry"])
        self.entry_row(grid2d, 7, "g11/g22/g12", self.g11, self.g22, self.g12)
        self.path_row(grid2d, 8, "Metric XYZ", self.metric_xyz, optional=True)

        command_frame = LabelFrame(root, text="Command")
        command_frame.pack(fill="x", padx=10, pady=4)
        self.command_box = Text(command_frame, height=4, wrap="word")
        self.command_box.pack(fill="x", padx=4, pady=4)

        output_frame = LabelFrame(root, text="Output")
        output_frame.pack(fill="both", expand=True, padx=10, pady=(4, 10))
        self.output_box = ScrolledText(output_frame, height=10, wrap="word")
        self.output_box.pack(fill="both", expand=True, padx=4, pady=4)

    def entry_row(self, parent: Frame, row: int, label: str, *variables: StringVar) -> None:
        Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=2)
        for column, variable in enumerate(variables, start=1):
            entry = Entry(parent, textvariable=variable)
            entry.grid(row=row, column=column, sticky="ew", padx=4, pady=2)
            entry.bind("<KeyRelease>", lambda _event: self.refresh_command())
            parent.columnconfigure(column, weight=1)

    def option_row(self, parent: Frame, row: int, label: str, variable: StringVar, options: list[str]) -> None:
        Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=2)
        OptionMenu(parent, variable, *options, command=lambda _value: self.refresh_command()).grid(
            row=row, column=1, sticky="ew", padx=4, pady=2
        )
        parent.columnconfigure(1, weight=1)

    def path_row(
        self,
        parent: Frame,
        row: int,
        label: str,
        variable: StringVar,
        *,
        directory: bool = False,
        optional: bool = False,
    ) -> None:
        Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=2)
        entry = Entry(parent, textvariable=variable)
        entry.grid(row=row, column=1, sticky="ew", padx=4, pady=2)
        entry.bind("<KeyRelease>", lambda _event: self.refresh_command())
        Button(
            parent,
            text="Browse",
            command=lambda: self.browse(variable, directory=directory, optional=optional),
        ).grid(row=row, column=2, padx=4, pady=2)

    def browse(self, variable: StringVar, *, directory: bool, optional: bool) -> None:
        if directory:
            value = filedialog.askdirectory(initialdir=str(ROOT))
        else:
            value = filedialog.askopenfilename(initialdir=str(ROOT))
        if value or not optional:
            variable.set(value)
            self.refresh_command()

    def add_if_value(self, command: list[str], option: str, value: StringVar) -> None:
        text = value.get().strip()
        if text:
            command.extend([option, text])

    def build_command(self) -> list[str]:
        workflow = self.workflow.get()
        command = [sys.executable, str(DRIVER)]
        source = self.source_path.get().strip()
        if workflow == "2D grid":
            command.extend(["--grid2d-csv", source or "<grid.csv>"])
        elif source.lower().endswith(".xyz"):
            command.extend(["--xyz", source])
        else:
            command.extend(["--gaussian-log", source or "<gaussian.log>"])
        if "--gaussian-log" in command and self.gaussian_energy.get() != "auto":
            command.extend(["--gaussian-energy", self.gaussian_energy.get()])

        self.add_if_value(command, "--prefix", self.prefix)
        self.add_if_value(command, "--outdir", self.outdir)
        self.add_if_value(command, "--figdir", self.figdir)
        self.add_if_value(command, "--grid", self.grid)
        self.add_if_value(command, "--levels", self.levels)
        self.add_if_value(command, "--save-states", self.save_states)
        self.add_if_value(command, "--plot-max-state", self.plot_max_state)
        self.add_if_value(command, "--plot-property-smooth-degree", self.plot_property_smooth_degree)
        if self.property_vpt2.get() != "auto":
            command.extend(["--property-vpt2", self.property_vpt2.get()])
        self.add_if_value(command, "--property-vpt2-fit-points", self.property_vpt2_fit_points)
        self.add_if_value(command, "--property-vpt2-degree", self.property_vpt2_degree)
        self.add_if_value(command, "--property-vpt2-basis-size", self.property_vpt2_basis_size)
        if self.solver.get() != "auto":
            command.extend(["--solver", self.solver.get()])

        if workflow == "Anharmonic Gaussian mode":
            command.extend(["--anharmonic-mode", self.anharmonic_mode.get().strip() or "1"])
            command.extend(["--anharmonic-mode-order", self.anharmonic_mode_order.get()])
            command.extend(["--well-type", self.well_type.get()])
            command.extend(["--anharmonic-derivative-model", self.anharmonic_model.get()])
            self.add_if_value(command, "--anharmonic-grid-half-width", self.anharmonic_width)
            self.add_if_value(command, "--anharmonic-cubic-threshold", self.anharmonic_cubic_threshold)
            self.add_if_value(command, "--anharmonic-handy-beta", self.anharmonic_handy_beta)
            self.add_if_value(command, "--anharmonic-kinetic-frequency-cm", self.anharmonic_kinetic_frequency)
            self.add_if_value(command, "--anharmonic-wall-height-cm", self.anharmonic_wall)
            command.extend(["--anharmonic-wall-degree", self.anharmonic_wall_degree.get()])
        elif workflow == "1D path":
            command.extend(["--log-selection", self.log_selection.get()])
            command.extend(["--boundary", self.boundary.get()])
            if self.properties_csv.get().strip():
                command.extend(["--properties-csv", self.properties_csv.get().strip()])
            if self.property_derivatives_csv.get().strip():
                command.extend(["--property-derivatives-csv", self.property_derivatives_csv.get().strip()])
            if self.vpt2_property_csv.get().strip():
                command.extend(["--vpt2-property-csv", self.vpt2_property_csv.get().strip()])
            self.add_if_value(command, "--energy-key", self.energy_key)
            if self.energy_unit.get() != "auto":
                command.extend(["--energy-unit", self.energy_unit.get()])
            self.add_if_value(command, "--property", self.property_name)
            self.add_if_value(command, "--temperature", self.temperature)
            if self.compute_rotconst.get():
                command.append("--compute-rotconst")
            if self.path_symmetry.get() != "none":
                command.extend(["--path-symmetry", self.path_symmetry.get()])
            command.extend(["--well-type", self.well_type.get()])
            if self.potential_extension.get() != "none":
                command.extend(["--potential-extension", self.potential_extension.get()])
                self.add_if_value(command, "--extension-length-au", self.extension_length)
                self.add_if_value(command, "--extension-points", self.extension_points)
                command.extend(["--extension-degree", self.extension_degree.get()])
                self.add_if_value(command, "--extension-target-cm", self.extension_target)
            if self.potential_smoothing.get() != "none":
                command.extend(["--potential-smoothing", self.potential_smoothing.get()])
                self.add_if_value(command, "--potential-spline-smoothing", self.potential_spline_smoothing)
        else:
            self.add_if_value(command, "--q1-key", self.q1_key)
            self.add_if_value(command, "--q2-key", self.q2_key)
            self.add_if_value(command, "--energy-key", self.energy_key)
            if self.energy_unit.get() != "auto":
                command.extend(["--energy-unit", self.energy_unit.get()])
            command.extend(["--boundary1", self.boundary1.get(), "--boundary2", self.boundary2.get()])
            self.add_if_value(command, "--basis1", self.basis1)
            self.add_if_value(command, "--basis2", self.basis2)
            command.extend(["--metric-mode", self.metric_mode.get()])
            if self.metric_mode.get() == "constant":
                self.add_if_value(command, "--g11", self.g11)
                self.add_if_value(command, "--g22", self.g22)
                self.add_if_value(command, "--g12", self.g12)
            elif self.metric_mode.get() == "geometry" and self.metric_xyz.get().strip():
                command.extend(["--grid2d-geom-xyz", self.metric_xyz.get().strip()])
        return command

    def refresh_command(self) -> None:
        command = self.build_command()
        text = " ".join(shlex.quote(part) for part in command)
        self.command_box.delete("1.0", "end")
        self.command_box.insert("1.0", text)

    def run_command(self) -> None:
        command = self.build_command()
        if any(part.startswith("<") and part.endswith(">") for part in command):
            messagebox.showerror("Missing input", "Choose an input file before running.")
            return
        self.run_button.configure(state="disabled")
        self.output_box.insert("end", "$ " + " ".join(shlex.quote(part) for part in command) + "\n")
        self.output_box.see("end")
        threading.Thread(target=self.run_worker, args=(command,), daemon=True).start()

    def run_worker(self, command: list[str]) -> None:
        try:
            process = subprocess.Popen(
                command,
                cwd=str(ROOT),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            assert process.stdout is not None
            for line in process.stdout:
                self.root.after(0, self.append_output, line)
            code = process.wait()
            self.root.after(0, self.append_output, f"\nExit code: {code}\n")
        except Exception as exc:
            self.root.after(0, self.append_output, f"\nError: {exc}\n")
        finally:
            self.root.after(0, lambda: self.run_button.configure(state="normal"))

    def append_output(self, text: str) -> None:
        self.output_box.insert("end", text)
        self.output_box.see("end")

    def run(self) -> None:
        self.root.mainloop()


def main() -> None:
    PathDVRGui().run()


if __name__ == "__main__":
    main()
