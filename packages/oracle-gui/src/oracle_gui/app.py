from __future__ import annotations

import argparse
from pathlib import Path

from .commands import OracleGuiCommand
from .dashboard import DashboardAction, OracleDashboardController
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
                self.refresh()

        def refresh(self) -> None:
            self.workflow_list.clear()
            self.section_table.setRowCount(0)
            self.action_table.setRowCount(0)
            self._clear_structure_tables()
            if self.controller.xyzin is None:
                self.path_label.setText("No ORACLE project loaded")
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
            self.structure_output_path.setText(str(state.xyzin))
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

        def _fill_table(self, widget: QTableWidget, table: StructureTable) -> None:
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
