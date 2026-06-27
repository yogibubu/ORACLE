from __future__ import annotations

import argparse
from pathlib import Path

from .commands import OracleGuiCommand
from .dashboard import DashboardAction, OracleDashboardController
from .gicforge import (
    GICForgeTable,
    OracleGICForgeController,
    default_gicforge_bmatrix_output,
    default_gicforge_gaussian_output,
    default_gicforge_report_output,
    gicforge_gui_state_lines,
    load_gicforge_gui_state,
)
from .project import load_oracle_project_state
from .structure import (
    OracleStructureController,
    StructureTable,
    default_preprocess_output,
    load_structure_gui_state,
)
from .workflows import ORACLE_GUI_WINDOWS


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="ORACLE desktop GUI")
    parser.add_argument("xyzin", nargs="?", type=Path, help="ORACLE enriched XYZ project")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return _run_qt(args.xyzin)
    except ModuleNotFoundError as exc:
        if exc.name and exc.name.startswith("PySide6"):
            raise SystemExit(
                "PySide6 is not installed. Install the GUI extra with "
                "`pip install -e packages/oracle-gui[qt]` or run the headless "
                "oracle_gui project controllers from Python."
            ) from exc
        raise


def _run_qt(initial_xyzin: Path | None) -> int:
    from PySide6.QtCore import QProcess, Qt
    from PySide6.QtWidgets import (
        QApplication,
        QFileDialog,
        QHBoxLayout,
        QCheckBox,
        QLabel,
        QComboBox,
        QLineEdit,
        QListWidget,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QTableWidget,
        QTableWidgetItem,
        QTabWidget,
        QTextEdit,
        QVBoxLayout,
        QWidget,
    )

    class OracleDashboardWindow(QMainWindow):
        def __init__(self, xyzin: Path | None = None) -> None:
            super().__init__()
            self.controller = OracleDashboardController(xyzin)
            self.structure_controller = OracleStructureController(xyzin)
            self.gicforge_controller = OracleGICForgeController(xyzin)
            self.process: QProcess | None = None
            self.current_actions: tuple[DashboardAction, ...] = ()
            self.pending_xyzin_after_run: Path | None = None
            self.setWindowTitle("ORACLE Project Dashboard")
            self.resize(1280, 780)

            root = QWidget()
            self.setCentralWidget(root)
            layout = QVBoxLayout(root)

            toolbar = QHBoxLayout()
            self.path_label = QLabel("No ORACLE project loaded")
            self.path_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
            open_button = QPushButton("Open")
            open_button.clicked.connect(self.open_project)
            refresh_button = QPushButton("Refresh")
            refresh_button.clicked.connect(self.refresh)
            self.run_button = QPushButton("Run Selected")
            self.run_button.clicked.connect(self.run_selected_action)
            toolbar.addWidget(self.path_label, stretch=1)
            toolbar.addWidget(open_button)
            toolbar.addWidget(refresh_button)
            toolbar.addWidget(self.run_button)
            layout.addLayout(toolbar)

            tabs = QTabWidget()
            state_tab = QWidget()
            body = QHBoxLayout(state_tab)
            self.workflow_list = QListWidget()
            self.workflow_list.itemDoubleClicked.connect(self.show_workflow_details)
            self.section_table = QTableWidget(0, 3)
            self.section_table.setHorizontalHeaderLabels(("Section", "Lines", "Schema"))
            body.addWidget(self.workflow_list, stretch=1)
            body.addWidget(self.section_table, stretch=2)
            tabs.addTab(state_tab, "Project State")

            actions_tab = QWidget()
            actions_layout = QVBoxLayout(actions_tab)
            self.action_table = QTableWidget(0, 4)
            self.action_table.setHorizontalHeaderLabels(("Action", "Status", "Window", "Command"))
            self.action_table.itemDoubleClicked.connect(self.run_selected_action)
            actions_layout.addWidget(self.action_table)
            tabs.addTab(actions_tab, "Actions")
            structure_tab = self._build_structure_tab()
            tabs.addTab(structure_tab, "Structure")
            gicforge_tab = self._build_gicforge_tab()
            tabs.addTab(gicforge_tab, "GICForge")
            layout.addWidget(tabs, stretch=2)

            self.details = QTextEdit()
            self.details.setReadOnly(True)
            self.log_output = QTextEdit()
            self.log_output.setReadOnly(True)

            bottom_tabs = QTabWidget()
            bottom_tabs.addTab(self.details, "Details")
            bottom_tabs.addTab(self.log_output, "Log")
            layout.addWidget(bottom_tabs, stretch=1)

            self.refresh()

        def open_project(self) -> None:
            path, _selected = QFileDialog.getOpenFileName(
                self,
                "Open ORACLE xyzin",
                str(Path.cwd()),
                "ORACLE xyzin (*.xyzin xyzin *.xyz);;All files (*)",
            )
            if path:
                self.controller.set_xyzin(Path(path))
                self.structure_controller.set_xyzin(Path(path))
                self.gicforge_controller.set_xyzin(Path(path))
                self.refresh()

        def refresh(self) -> None:
            self.workflow_list.clear()
            self.section_table.setRowCount(0)
            self.action_table.setRowCount(0)
            self._clear_structure_tables()
            self._clear_gicforge_tables()
            if self.controller.xyzin is None:
                self.path_label.setText("No ORACLE project loaded")
                self.gicforge_summary.setPlainText("No ORACLE project loaded")
                self.details.setPlainText(
                    "\n".join(
                        f"{spec.title}: {spec.description}" for spec in ORACLE_GUI_WINDOWS
                    )
                )
                self.log_output.setPlainText("\n".join(self.controller.log_lines))
                self.run_button.setEnabled(False)
                return

            state = load_oracle_project_state(self.controller.xyzin)
            self.structure_controller.set_xyzin(state.xyzin)
            self.gicforge_controller.set_xyzin(state.xyzin)
            self.structure_output_path.setText(str(state.xyzin))
            self._set_default_gicforge_outputs(state.xyzin)
            self.path_label.setText(str(state.xyzin))
            for workflow in state.workflows:
                self.workflow_list.addItem(
                    f"{workflow.status.value.upper():8s}  {workflow.title}"
                )

            self.section_table.setRowCount(len(state.sections))
            for row, section in enumerate(state.sections):
                self.section_table.setItem(row, 0, QTableWidgetItem(section.name))
                self.section_table.setItem(row, 1, QTableWidgetItem(str(section.line_count)))
                self.section_table.setItem(row, 2, QTableWidgetItem(section.schema or ""))
            self.section_table.resizeColumnsToContents()
            self._populate_actions(state)

            self.details.setPlainText(
                "\n".join(
                    [
                        f"atoms: {state.atom_count}",
                        f"point group: {state.point_group or 'unknown'}",
                        f"validation: {state.validation_status}",
                        "",
                        *state.validation_messages,
                    ]
                )
            )
            self.log_output.setPlainText("\n".join(self.controller.log_lines))
            self.run_button.setEnabled(
                self.process is None and any(action.enabled for action in self.current_actions)
            )
            self._populate_structure_tables(state.xyzin)
            self._populate_gicforge_tables(state.xyzin)

        def show_workflow_details(self, _item=None) -> None:
            if self.controller.xyzin is None:
                return
            state = load_oracle_project_state(self.controller.xyzin)
            row = self.workflow_list.currentRow()
            if row < 0 or row >= len(state.workflows):
                return
            workflow = state.workflows[row]
            QMessageBox.information(
                self,
                workflow.title,
                "\n".join(
                    [
                        f"Status: {workflow.status.value}",
                        workflow.message,
                        "",
                        "Required: "
                        + (", ".join(workflow.required_sections) or "none"),
                        "Produced: "
                        + (", ".join(workflow.produced_sections) or "none"),
                        "",
                        *self._spec_extra_lines(workflow.key),
                    ]
                ),
            )

        def _spec_extra_lines(self, key: str) -> list[str]:
            for spec in ORACLE_GUI_WINDOWS:
                if spec.key != key:
                    continue
                lines = [f"Category: {spec.category}"]
                if spec.capabilities:
                    lines.append("Capabilities:")
                    lines.extend(f"- {item}" for item in spec.capabilities)
                if spec.publication_exports:
                    lines.append("Publication export:")
                    lines.extend(f"- {item}" for item in spec.publication_exports)
                if spec.external_viewers:
                    lines.append("External viewers:")
                    lines.extend(f"- {item}" for item in spec.external_viewers)
                return lines
            return []

        def _populate_actions(self, state) -> None:
            self.current_actions = self.controller.actions(state)
            self.action_table.setRowCount(len(self.current_actions))
            for row, action in enumerate(self.current_actions):
                values = (
                    action.label,
                    "ready" if action.enabled else action.reason,
                    action.window_key,
                    action.shell_line,
                )
                for column, value in enumerate(values):
                    item = QTableWidgetItem(value)
                    if not action.enabled:
                        item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
                    self.action_table.setItem(row, column, item)
            self.action_table.resizeColumnsToContents()

        def run_selected_action(self, _item=None) -> None:
            row = self.action_table.currentRow()
            if row < 0 or row >= len(self.current_actions):
                QMessageBox.information(self, "ORACLE", "Select an action first.")
                return
            action = self.current_actions[row]
            if not action.enabled:
                QMessageBox.warning(self, action.label, action.reason)
                return
            if self.process is not None:
                QMessageBox.information(self, "ORACLE", "A command is already running.")
                return
            self._start_process(action)

        def _start_process(self, action: DashboardAction) -> None:
            self._start_command(action.command, action.label)

        def _start_command(self, command: OracleGuiCommand, label: str) -> None:
            prepared = self.controller.prepare_command(command)
            self.controller.log(f"$ {prepared.shell_line()}")
            self.log_output.setPlainText("\n".join(self.controller.log_lines))
            process = QProcess(self)
            if prepared.cwd is not None:
                process.setWorkingDirectory(str(prepared.cwd))
            process.setProgram(prepared.argv[0])
            process.setArguments(list(prepared.argv[1:]))
            process.readyReadStandardOutput.connect(self._read_process_stdout)
            process.readyReadStandardError.connect(self._read_process_stderr)
            process.errorOccurred.connect(lambda error: self._process_error(label, error))
            process.finished.connect(lambda code, _status: self._process_finished(label, code))
            self.process = process
            self.run_button.setEnabled(False)
            process.start()

        def _read_process_stdout(self) -> None:
            if self.process is None:
                return
            text = bytes(self.process.readAllStandardOutput()).decode(errors="replace")
            self.controller.log(text.rstrip())
            self.log_output.setPlainText("\n".join(self.controller.log_lines))

        def _read_process_stderr(self) -> None:
            if self.process is None:
                return
            text = bytes(self.process.readAllStandardError()).decode(errors="replace")
            self.controller.log(text.rstrip())
            self.log_output.setPlainText("\n".join(self.controller.log_lines))

        def _process_finished(self, label: str, code: int) -> None:
            self.controller.log(f"[exit {int(code)}] {label}")
            self.process = None
            if self.pending_xyzin_after_run is not None and self.pending_xyzin_after_run.exists():
                self.controller.set_xyzin(self.pending_xyzin_after_run)
                self.structure_controller.set_xyzin(self.pending_xyzin_after_run)
                self.gicforge_controller.set_xyzin(self.pending_xyzin_after_run)
            self.pending_xyzin_after_run = None
            self.refresh()

        def _process_error(self, label: str, error) -> None:
            self.controller.log(f"[process error {int(error)}] {label}")
            if self.process is not None and self.process.state() == QProcess.NotRunning:
                self.process = None
                self.refresh()

        def _build_structure_tab(self) -> QWidget:
            tab = QWidget()
            layout = QVBoxLayout(tab)

            source_row = QHBoxLayout()
            source_row.addWidget(QLabel("Source"))
            self.structure_source_path = QLineEdit()
            browse_source = QPushButton("Browse")
            browse_source.clicked.connect(self.browse_structure_source)
            source_row.addWidget(self.structure_source_path, stretch=1)
            source_row.addWidget(browse_source)
            layout.addLayout(source_row)

            output_row = QHBoxLayout()
            output_row.addWidget(QLabel("Output xyzin"))
            self.structure_output_path = QLineEdit()
            browse_output = QPushButton("Save As")
            browse_output.clicked.connect(self.browse_structure_output)
            output_row.addWidget(self.structure_output_path, stretch=1)
            output_row.addWidget(browse_output)
            layout.addLayout(output_row)

            controls = QHBoxLayout()
            self.structure_source_kind = QComboBox()
            self.structure_source_kind.addItems(("auto", "xyz", "enriched_xyz", "gaussian", "molpro", "mrcc"))
            preprocess_button = QPushButton("Preprocess")
            preprocess_button.clicked.connect(self.run_structure_preprocess)
            avogadro_button = QPushButton("Avogadro")
            avogadro_button.clicked.connect(self.run_structure_avogadro)
            fragments_button = QPushButton("Build Fragments")
            fragments_button.clicked.connect(self.run_structure_fragments)
            controls.addWidget(QLabel("Kind"))
            controls.addWidget(self.structure_source_kind)
            controls.addWidget(preprocess_button)
            controls.addWidget(avogadro_button)
            controls.addWidget(fragments_button)
            controls.addStretch(1)
            layout.addLayout(controls)

            self.structure_table_tabs = QTabWidget()
            self.structure_bond_table = QTableWidget(0, 2)
            self.structure_ring_table = QTableWidget(0, 3)
            self.structure_synthon_table = QTableWidget(0, 8)
            self.structure_fragment_table = QTableWidget(0, 5)
            self.structure_table_tabs.addTab(self.structure_bond_table, "Bonds")
            self.structure_table_tabs.addTab(self.structure_ring_table, "Rings")
            self.structure_table_tabs.addTab(self.structure_synthon_table, "Synthons")
            self.structure_table_tabs.addTab(self.structure_fragment_table, "Fragments")
            layout.addWidget(self.structure_table_tabs, stretch=1)
            return tab

        def browse_structure_source(self) -> None:
            path, _selected = QFileDialog.getOpenFileName(
                self,
                "Import molecule source",
                str(Path.cwd()),
                "Molecule sources (*.xyz *.xyzin *.gjf *.gau *.com *.inp *.log *.out *.zmat *.zmt);;All files (*)",
            )
            if not path:
                return
            source = Path(path)
            self.structure_source_path.setText(str(source))
            if not self.structure_output_path.text().strip():
                self.structure_output_path.setText(str(default_preprocess_output(source)))

        def browse_structure_output(self) -> None:
            path, _selected = QFileDialog.getSaveFileName(
                self,
                "Write ORACLE xyzin",
                self.structure_output_path.text().strip() or str(Path.cwd() / "molecule.xyzin"),
                "ORACLE xyzin (*.xyzin *.xyz);;All files (*)",
            )
            if path:
                self.structure_output_path.setText(path)

        def run_structure_preprocess(self) -> None:
            source_text = self.structure_source_path.text().strip()
            output_text = self.structure_output_path.text().strip()
            if not source_text or not output_text:
                QMessageBox.warning(self, "ORACLE-Babel", "Select source and output files first.")
                return
            if self.process is not None:
                QMessageBox.information(self, "ORACLE", "A command is already running.")
                return
            output = Path(output_text)
            command = self.structure_controller.preprocess_command(
                Path(source_text),
                output,
                source_kind=self.structure_source_kind.currentText(),
            )
            self.pending_xyzin_after_run = output
            self._start_command(command, command.label)

        def run_structure_avogadro(self) -> None:
            if self.controller.xyzin is None:
                QMessageBox.warning(self, "Avogadro", "Open or preprocess an ORACLE xyzin first.")
                return
            if self.process is not None:
                QMessageBox.information(self, "ORACLE", "A command is already running.")
                return
            self.structure_controller.set_xyzin(self.controller.xyzin)
            command = self.structure_controller.avogadro_command()
            self._start_command(command, command.label)

        def run_structure_fragments(self) -> None:
            if self.controller.xyzin is None:
                QMessageBox.warning(self, "Fragments", "Open or preprocess an ORACLE xyzin first.")
                return
            if self.process is not None:
                QMessageBox.information(self, "ORACLE", "A command is already running.")
                return
            self.structure_controller.set_xyzin(self.controller.xyzin)
            command = self.structure_controller.fragments_command("build")
            self._start_command(command, command.label)

        def _populate_structure_tables(self, xyzin: Path) -> None:
            state = load_structure_gui_state(xyzin)
            self._fill_table(self.structure_bond_table, state.topology_bonds)
            self._fill_table(self.structure_ring_table, state.topology_rings)
            self._fill_table(self.structure_synthon_table, state.synthons)
            self._fill_table(self.structure_fragment_table, state.fragments)

        def _clear_structure_tables(self) -> None:
            for table in (
                getattr(self, "structure_bond_table", None),
                getattr(self, "structure_ring_table", None),
                getattr(self, "structure_synthon_table", None),
                getattr(self, "structure_fragment_table", None),
            ):
                if table is not None:
                    table.setRowCount(0)

        def _build_gicforge_tab(self) -> QWidget:
            tab = QWidget()
            layout = QVBoxLayout(tab)

            option_row = QHBoxLayout()
            self.gicforge_symmetrize = QCheckBox("Symmetrize")
            self.gicforge_symmetrize.setChecked(True)
            self.gicforge_sycart = QCheckBox("SYCART")
            self.gicforge_sycart.setChecked(True)
            self.gicforge_improper = QCheckBox("Improper dihedrals")
            self.gicforge_improper.setChecked(True)
            build_button = QPushButton("Build GICs")
            build_button.clicked.connect(self.run_gicforge_build)
            option_row.addWidget(self.gicforge_symmetrize)
            option_row.addWidget(self.gicforge_sycart)
            option_row.addWidget(self.gicforge_improper)
            option_row.addWidget(build_button)
            option_row.addStretch(1)
            layout.addLayout(option_row)

            report_row = QHBoxLayout()
            report_row.addWidget(QLabel("Report"))
            self.gicforge_report_output = QLineEdit()
            report_browse = QPushButton("Save As")
            report_browse.clicked.connect(self.browse_gicforge_report_output)
            report_button = QPushButton("Write Report")
            report_button.clicked.connect(self.run_gicforge_report)
            report_row.addWidget(self.gicforge_report_output, stretch=1)
            report_row.addWidget(report_browse)
            report_row.addWidget(report_button)
            layout.addLayout(report_row)

            bmatrix_row = QHBoxLayout()
            bmatrix_row.addWidget(QLabel("B Matrix"))
            self.gicforge_bmatrix_output = QLineEdit()
            bmatrix_browse = QPushButton("Save As")
            bmatrix_browse.clicked.connect(self.browse_gicforge_bmatrix_output)
            bmatrix_button = QPushButton("Evaluate B")
            bmatrix_button.clicked.connect(self.run_gicforge_bmatrix)
            bmatrix_row.addWidget(self.gicforge_bmatrix_output, stretch=1)
            bmatrix_row.addWidget(bmatrix_browse)
            bmatrix_row.addWidget(bmatrix_button)
            layout.addLayout(bmatrix_row)

            gaussian_row = QHBoxLayout()
            gaussian_row.addWidget(QLabel("Gaussian"))
            self.gicforge_gaussian_output = QLineEdit()
            gaussian_browse = QPushButton("Save As")
            gaussian_browse.clicked.connect(self.browse_gicforge_gaussian_output)
            self.gicforge_gaussian_route = QLineEdit("#p hf/sto-3g")
            gaussian_button = QPushButton("Write Input")
            gaussian_button.clicked.connect(self.run_gicforge_gaussian_input)
            gaussian_row.addWidget(self.gicforge_gaussian_output, stretch=1)
            gaussian_row.addWidget(gaussian_browse)
            gaussian_row.addWidget(QLabel("Route"))
            gaussian_row.addWidget(self.gicforge_gaussian_route)
            gaussian_row.addWidget(gaussian_button)
            layout.addLayout(gaussian_row)

            self.gicforge_summary = QTextEdit()
            self.gicforge_summary.setReadOnly(True)
            self.gicforge_summary.setMaximumHeight(150)
            layout.addWidget(self.gicforge_summary)

            self.gicforge_table_tabs = QTabWidget()
            self.gicforge_primitive_table = QTableWidget(0, 8)
            self.gicforge_frozen_table = QTableWidget(0, 7)
            self.gicforge_symmetry_table = QTableWidget(0, 5)
            self.gicforge_diagnostics_table = QTableWidget(0, 2)
            self.gicforge_table_tabs.addTab(self.gicforge_primitive_table, "Primitives")
            self.gicforge_table_tabs.addTab(self.gicforge_frozen_table, "Frozen GICs")
            self.gicforge_table_tabs.addTab(self.gicforge_symmetry_table, "Symmetry")
            self.gicforge_table_tabs.addTab(self.gicforge_diagnostics_table, "Diagnostics")
            layout.addWidget(self.gicforge_table_tabs, stretch=1)
            return tab

        def browse_gicforge_report_output(self) -> None:
            self._browse_gicforge_output(
                self.gicforge_report_output,
                "Write GICForge report",
                "Text reports (*.txt);;All files (*)",
            )

        def browse_gicforge_bmatrix_output(self) -> None:
            self._browse_gicforge_output(
                self.gicforge_bmatrix_output,
                "Write GIC B matrix",
                "Text files (*.txt);;All files (*)",
            )

        def browse_gicforge_gaussian_output(self) -> None:
            self._browse_gicforge_output(
                self.gicforge_gaussian_output,
                "Write Gaussian GIC input",
                "Gaussian input (*.gjf *.com);;All files (*)",
            )

        def _browse_gicforge_output(
            self,
            field: QLineEdit,
            title: str,
            file_filter: str,
        ) -> None:
            path, _selected = QFileDialog.getSaveFileName(
                self,
                title,
                field.text().strip() or str(Path.cwd()),
                file_filter,
            )
            if path:
                field.setText(path)

        def run_gicforge_build(self) -> None:
            if not self._ensure_gicforge_ready_for_command("GICForge"):
                return
            self.gicforge_controller.set_xyzin(self.controller.xyzin)
            command = self.gicforge_controller.build_command(
                symmetrize=self.gicforge_symmetrize.isChecked(),
                sycart=self.gicforge_sycart.isChecked(),
                improper_dihedrals=self.gicforge_improper.isChecked(),
            )
            self._start_command(command, command.label)

        def run_gicforge_report(self) -> None:
            if not self._ensure_gicforge_ready_for_command("GICForge report"):
                return
            if not self._ensure_gicforge_section_ready("GICForge report"):
                return
            output = self._gicforge_output_path(
                self.gicforge_report_output,
                default_gicforge_report_output,
            )
            command = self.gicforge_controller.report_command(output)
            self._start_command(command, command.label)

        def run_gicforge_bmatrix(self) -> None:
            if not self._ensure_gicforge_ready_for_command("GICForge B matrix"):
                return
            if not self._ensure_gicforge_section_ready("GICForge B matrix"):
                return
            output = self._gicforge_output_path(
                self.gicforge_bmatrix_output,
                default_gicforge_bmatrix_output,
            )
            command = self.gicforge_controller.bmatrix_command(output)
            self._start_command(command, command.label)

        def run_gicforge_gaussian_input(self) -> None:
            if not self._ensure_gicforge_ready_for_command("Gaussian GIC input"):
                return
            if not self._ensure_gicforge_section_ready("Gaussian GIC input"):
                return
            output = self._gicforge_output_path(
                self.gicforge_gaussian_output,
                default_gicforge_gaussian_output,
            )
            route = self.gicforge_gaussian_route.text().strip() or "#p hf/sto-3g"
            command = self.gicforge_controller.gaussian_input_command(output, route=route)
            self._start_command(command, command.label)

        def _ensure_gicforge_ready_for_command(self, title: str) -> bool:
            if self.controller.xyzin is None:
                QMessageBox.warning(self, title, "Open or preprocess an ORACLE xyzin first.")
                return False
            if self.process is not None:
                QMessageBox.information(self, "ORACLE", "A command is already running.")
                return False
            self.gicforge_controller.set_xyzin(self.controller.xyzin)
            return True

        def _ensure_gicforge_section_ready(self, title: str) -> bool:
            if self.controller.xyzin is None:
                return False
            state = load_gicforge_gui_state(self.controller.xyzin)
            if state.ready:
                return True
            message = "\n".join(state.messages) or "Build GICs first."
            QMessageBox.warning(self, title, message)
            return False

        def _gicforge_output_path(self, field: QLineEdit, default_factory) -> Path:
            if self.controller.xyzin is None:
                raise ValueError("no ORACLE xyzin project is loaded")
            text = field.text().strip()
            output = Path(text) if text else default_factory(self.controller.xyzin)
            field.setText(str(output))
            return output

        def _set_default_gicforge_outputs(self, xyzin: Path) -> None:
            self.gicforge_report_output.setText(str(default_gicforge_report_output(xyzin)))
            self.gicforge_bmatrix_output.setText(str(default_gicforge_bmatrix_output(xyzin)))
            self.gicforge_gaussian_output.setText(str(default_gicforge_gaussian_output(xyzin)))

        def _populate_gicforge_tables(self, xyzin: Path) -> None:
            state = load_gicforge_gui_state(xyzin)
            self.gicforge_summary.setPlainText("\n".join(gicforge_gui_state_lines(state)))
            self._fill_table(self.gicforge_primitive_table, state.primitives)
            self._fill_table(self.gicforge_frozen_table, state.frozen_gics)
            self._fill_table(self.gicforge_symmetry_table, state.symmetry_groups)
            self._fill_table(self.gicforge_diagnostics_table, state.diagnostics)

        def _clear_gicforge_tables(self) -> None:
            summary = getattr(self, "gicforge_summary", None)
            if summary is not None:
                summary.clear()
            for table in (
                getattr(self, "gicforge_primitive_table", None),
                getattr(self, "gicforge_frozen_table", None),
                getattr(self, "gicforge_symmetry_table", None),
                getattr(self, "gicforge_diagnostics_table", None),
            ):
                if table is not None:
                    table.setRowCount(0)

        def _fill_table(self, widget: QTableWidget, table: StructureTable | GICForgeTable) -> None:
            widget.setColumnCount(len(table.columns))
            widget.setHorizontalHeaderLabels(table.columns)
            widget.setRowCount(len(table.rows))
            for row_index, row in enumerate(table.rows):
                for column_index, value in enumerate(row):
                    widget.setItem(row_index, column_index, QTableWidgetItem(str(value)))
            widget.resizeColumnsToContents()

    app = QApplication([])
    window = OracleDashboardWindow(initial_xyzin)
    window.show()
    return int(app.exec())


if __name__ == "__main__":
    raise SystemExit(main())
