from __future__ import annotations

import argparse
from pathlib import Path

from .commands import OracleGuiCommand
from .contracts import (
    ToolContractTable,
    load_tool_contract_gui_state,
    tool_contract_gui_state_lines,
)
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
from .gf import (
    GFTable,
    OracleGFController,
    default_gf_csv_dir,
    default_gf_report_output,
    gf_gui_state_lines,
    load_gf_gui_state,
)
from .project import load_oracle_project_state
from .sefit import (
    OracleSEFitController,
    SEFitTable,
    default_sefit_outdir,
    load_sefit_gui_state,
    sefit_gui_state_lines,
)
from .structure import (
    OracleStructureController,
    StructureTable,
    default_preprocess_output,
    load_structure_gui_state,
)
from .trinity import (
    OracleTrinityController,
    TrinityTable,
    default_trinity_run_dir,
    load_trinity_gui_state,
    trinity_gui_state_lines,
)
from .workbench import (
    WorkbenchTable,
    load_workbench_gui_state,
    workbench_gui_state_lines,
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
            self.gf_controller = OracleGFController(xyzin)
            self.sefit_controller = OracleSEFitController(xyzin)
            self.trinity_controller = OracleTrinityController(xyzin)
            self.process: QProcess | None = None
            self.current_actions: tuple[DashboardAction, ...] = ()
            self.pending_xyzin_after_run: Path | None = None
            self.workbench_tabs = (
                ("rovib_thermo", "Rovib/Thermo"),
                ("anharmonic", "Anharmonic"),
                ("qm_jobs", "QM Jobs"),
                ("diagnostics", "Diagnostics"),
                ("rotational_spectroscopy", "Rotational"),
                ("vibrational_spectroscopy", "Vibrational"),
                ("electronic_spectroscopy", "Electronic"),
                ("thermochemistry_kinetics", "Thermo/Kinetics"),
            )
            self.workbench_summaries: dict[str, QTextEdit] = {}
            self.workbench_section_tables: dict[str, QTableWidget] = {}
            self.workbench_action_tables: dict[str, QTableWidget] = {}
            self.workbench_capability_tables: dict[str, QTableWidget] = {}
            self.workbench_export_tables: dict[str, QTableWidget] = {}
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
            contracts_tab = self._build_contracts_tab()
            tabs.addTab(contracts_tab, "Tool Contracts")
            structure_tab = self._build_structure_tab()
            tabs.addTab(structure_tab, "Structure")
            gicforge_tab = self._build_gicforge_tab()
            tabs.addTab(gicforge_tab, "GICForge")
            gf_tab = self._build_gf_tab()
            tabs.addTab(gf_tab, "GF/PED")
            sefit_tab = self._build_sefit_tab()
            tabs.addTab(sefit_tab, "SEFit")
            trinity_tab = self._build_trinity_tab()
            tabs.addTab(trinity_tab, "TRINITY")
            for key, label in self.workbench_tabs:
                tabs.addTab(self._build_workbench_tab(key), label)
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
                self.gf_controller.set_xyzin(Path(path))
                self.sefit_controller.set_xyzin(Path(path))
                self.trinity_controller.set_xyzin(Path(path))
                self.refresh()

        def refresh(self) -> None:
            self.workflow_list.clear()
            self.section_table.setRowCount(0)
            self.action_table.setRowCount(0)
            self._clear_contracts_table()
            self._clear_structure_tables()
            self._clear_gicforge_tables()
            self._clear_gf_tables()
            self._clear_sefit_tables()
            self._clear_trinity_tables()
            self._clear_workbench_tables()
            if self.controller.xyzin is None:
                self.path_label.setText("No ORACLE project loaded")
                self.gicforge_summary.setPlainText("No ORACLE project loaded")
                self.gf_summary.setPlainText("No ORACLE project loaded")
                self.sefit_summary.setPlainText("No ORACLE project loaded")
                self.trinity_summary.setPlainText("No ORACLE project loaded")
                self.contracts_summary.setPlainText("No ORACLE project loaded")
                for summary in self.workbench_summaries.values():
                    summary.setPlainText("No ORACLE project loaded")
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
            self.gf_controller.set_xyzin(state.xyzin)
            self.sefit_controller.set_xyzin(state.xyzin)
            self.trinity_controller.set_xyzin(state.xyzin)
            self.structure_output_path.setText(str(state.xyzin))
            self._set_default_gicforge_outputs(state.xyzin)
            self._set_default_gf_outputs(state.xyzin)
            self._set_default_sefit_outputs(state.xyzin)
            self._set_default_trinity_outputs(state.xyzin)
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
            self._populate_contracts_table(state.xyzin)

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
            self._populate_gf_tables(state.xyzin)
            self._populate_sefit_tables(state.xyzin)
            self._populate_trinity_tables(state.xyzin)
            self._populate_workbench_tables(state.xyzin)

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

        def _build_contracts_tab(self) -> QWidget:
            tab = QWidget()
            layout = QVBoxLayout(tab)
            self.contracts_summary = QTextEdit()
            self.contracts_summary.setReadOnly(True)
            self.contracts_summary.setMaximumHeight(90)
            layout.addWidget(self.contracts_summary)
            self.contracts_table = QTableWidget(0, 8)
            layout.addWidget(self.contracts_table, stretch=1)
            return tab

        def _populate_contracts_table(self, xyzin: Path) -> None:
            state = load_tool_contract_gui_state(xyzin)
            self.contracts_summary.setPlainText(
                "\n".join(tool_contract_gui_state_lines(state))
            )
            self._fill_table(self.contracts_table, state.table)

        def _clear_contracts_table(self) -> None:
            summary = getattr(self, "contracts_summary", None)
            if summary is not None:
                summary.clear()
            table = getattr(self, "contracts_table", None)
            if table is not None:
                table.setRowCount(0)

        def _build_workbench_tab(self, key: str) -> QWidget:
            tab = QWidget()
            layout = QVBoxLayout(tab)

            summary = QTextEdit()
            summary.setReadOnly(True)
            summary.setMaximumHeight(120)
            self.workbench_summaries[key] = summary
            layout.addWidget(summary)

            tables = QTabWidget()
            section_table = QTableWidget(0, 3)
            action_table = QTableWidget(0, 5)
            capability_table = QTableWidget(0, 1)
            export_table = QTableWidget(0, 2)
            self.workbench_section_tables[key] = section_table
            self.workbench_action_tables[key] = action_table
            self.workbench_capability_tables[key] = capability_table
            self.workbench_export_tables[key] = export_table
            tables.addTab(section_table, "Sections")
            tables.addTab(action_table, "Actions")
            tables.addTab(capability_table, "Capabilities")
            tables.addTab(export_table, "Exports")
            layout.addWidget(tables, stretch=1)
            return tab

        def _populate_workbench_tables(self, xyzin: Path) -> None:
            for key, _label in self.workbench_tabs:
                state = load_workbench_gui_state(xyzin, key)
                self.workbench_summaries[key].setPlainText(
                    "\n".join(workbench_gui_state_lines(state))
                )
                self._fill_table(self.workbench_section_tables[key], state.sections)
                self._fill_table(self.workbench_action_tables[key], state.actions)
                self._fill_table(self.workbench_capability_tables[key], state.capabilities)
                self._fill_table(self.workbench_export_tables[key], state.exports)

        def _clear_workbench_tables(self) -> None:
            for summary in getattr(self, "workbench_summaries", {}).values():
                summary.clear()
            for tables in (
                getattr(self, "workbench_section_tables", {}),
                getattr(self, "workbench_action_tables", {}),
                getattr(self, "workbench_capability_tables", {}),
                getattr(self, "workbench_export_tables", {}),
            ):
                for table in tables.values():
                    table.setRowCount(0)

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
                self.gf_controller.set_xyzin(self.pending_xyzin_after_run)
                self.sefit_controller.set_xyzin(self.pending_xyzin_after_run)
                self.trinity_controller.set_xyzin(self.pending_xyzin_after_run)
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

        def _build_gf_tab(self) -> QWidget:
            tab = QWidget()
            layout = QVBoxLayout(tab)

            source_row = QHBoxLayout()
            source_row.addWidget(QLabel("FCHK"))
            self.gf_fchk_path = QLineEdit()
            fchk_browse = QPushButton("Browse")
            fchk_browse.clicked.connect(self.browse_gf_fchk)
            source_row.addWidget(self.gf_fchk_path, stretch=1)
            source_row.addWidget(fchk_browse)
            layout.addLayout(source_row)

            output_row = QHBoxLayout()
            output_row.addWidget(QLabel("Report"))
            self.gf_report_output = QLineEdit()
            report_browse = QPushButton("Save As")
            report_browse.clicked.connect(self.browse_gf_report_output)
            output_row.addWidget(self.gf_report_output, stretch=1)
            output_row.addWidget(report_browse)
            output_row.addWidget(QLabel("CSV dir"))
            self.gf_csv_dir = QLineEdit()
            csv_browse = QPushButton("Browse")
            csv_browse.clicked.connect(self.browse_gf_csv_dir)
            output_row.addWidget(self.gf_csv_dir, stretch=1)
            output_row.addWidget(csv_browse)
            layout.addLayout(output_row)

            scale_row = QHBoxLayout()
            scale_row.addWidget(QLabel("Scale file"))
            self.gf_scale_file = QLineEdit()
            scale_browse = QPushButton("Browse")
            scale_browse.clicked.connect(self.browse_gf_scale_file)
            self.gf_scale_records = QLineEdit()
            scale_row.addWidget(self.gf_scale_file, stretch=1)
            scale_row.addWidget(scale_browse)
            scale_row.addWidget(QLabel("Scale"))
            scale_row.addWidget(self.gf_scale_records, stretch=1)
            layout.addLayout(scale_row)

            options_row = QHBoxLayout()
            self.gf_symmetry_blocks = QCheckBox("Symmetry blocks")
            self.gf_symmetry_blocks.setChecked(True)
            self.gf_local = QCheckBox("Local")
            self.gf_subtract_electrostatic = QCheckBox("Subtract electrostatic")
            self.gf_subtract_uff_vdw = QCheckBox("Subtract UFF vdW")
            self.gf_write_section = QCheckBox("Update #GF_PED")
            self.gf_write_section.setChecked(True)
            self.gf_force_threshold = QLineEdit()
            self.gf_force_threshold.setPlaceholderText("force threshold")
            self.gf_nonbonded_14_scale = QLineEdit("0.5")
            run_button = QPushButton("Run GF/PED")
            run_button.clicked.connect(self.run_gf)
            options_row.addWidget(self.gf_symmetry_blocks)
            options_row.addWidget(self.gf_local)
            options_row.addWidget(self.gf_subtract_electrostatic)
            options_row.addWidget(self.gf_subtract_uff_vdw)
            options_row.addWidget(self.gf_write_section)
            options_row.addWidget(QLabel("Threshold"))
            options_row.addWidget(self.gf_force_threshold)
            options_row.addWidget(QLabel("1-4"))
            options_row.addWidget(self.gf_nonbonded_14_scale)
            options_row.addWidget(run_button)
            layout.addLayout(options_row)

            self.gf_summary = QTextEdit()
            self.gf_summary.setReadOnly(True)
            self.gf_summary.setMaximumHeight(150)
            layout.addWidget(self.gf_summary)

            self.gf_table_tabs = QTabWidget()
            self.gf_frequency_table = QTableWidget(0, 3)
            self.gf_gic_table = QTableWidget(0, 7)
            self.gf_ped_table = QTableWidget(0, 3)
            self.gf_diagnostics_table = QTableWidget(0, 2)
            self.gf_table_tabs.addTab(self.gf_frequency_table, "Frequencies")
            self.gf_table_tabs.addTab(self.gf_gic_table, "GICs")
            self.gf_table_tabs.addTab(self.gf_ped_table, "PED")
            self.gf_table_tabs.addTab(self.gf_diagnostics_table, "Diagnostics")
            layout.addWidget(self.gf_table_tabs, stretch=1)
            return tab

        def browse_gf_fchk(self) -> None:
            path, _selected = QFileDialog.getOpenFileName(
                self,
                "Select Gaussian FCHK",
                str(Path.cwd()),
                "Gaussian FCHK (*.fchk *.fch);;All files (*)",
            )
            if path:
                self.gf_fchk_path.setText(path)

        def browse_gf_report_output(self) -> None:
            path, _selected = QFileDialog.getSaveFileName(
                self,
                "Write GF/PED report",
                self.gf_report_output.text().strip() or str(Path.cwd() / "gf_ped_report.txt"),
                "Text reports (*.txt *.report);;All files (*)",
            )
            if path:
                self.gf_report_output.setText(path)

        def browse_gf_csv_dir(self) -> None:
            path = QFileDialog.getExistingDirectory(
                self,
                "Select GF/PED CSV directory",
                self.gf_csv_dir.text().strip() or str(Path.cwd()),
            )
            if path:
                self.gf_csv_dir.setText(path)

        def browse_gf_scale_file(self) -> None:
            path, _selected = QFileDialog.getOpenFileName(
                self,
                "Select GF scaling file",
                str(Path.cwd()),
                "Scaling files (*.txt *.csv *.dat);;All files (*)",
            )
            if path:
                self.gf_scale_file.setText(path)

        def run_gf(self) -> None:
            if not self._ensure_gf_ready_for_command("GF/PED"):
                return
            if not self._ensure_gf_inputs_ready("GF/PED"):
                return
            force_threshold = self._optional_float_field(
                self.gf_force_threshold,
                "GF/PED",
                "force threshold",
            )
            if force_threshold is False:
                return
            nonbonded_14_scale = self._optional_float_field(
                self.gf_nonbonded_14_scale,
                "GF/PED",
                "1-4 scale",
                default=0.5,
            )
            if nonbonded_14_scale is False:
                return
            self.gf_controller.set_xyzin(self.controller.xyzin)
            command = self.gf_controller.run_command(
                fchk=self._optional_path_text(self.gf_fchk_path),
                out=self._gf_output_path(self.gf_report_output, default_gf_report_output),
                csv_dir=self._gf_output_path(self.gf_csv_dir, default_gf_csv_dir),
                scale_file=self._optional_path_text(self.gf_scale_file),
                scale_records=self._scale_records(),
                local=self.gf_local.isChecked(),
                symmetry_blocks=self.gf_symmetry_blocks.isChecked(),
                force_threshold=force_threshold,
                subtract_electrostatic=self.gf_subtract_electrostatic.isChecked(),
                subtract_uff_vdw=self.gf_subtract_uff_vdw.isChecked(),
                nonbonded_14_scale=nonbonded_14_scale,
                write_section=self.gf_write_section.isChecked(),
            )
            self._start_command(command, command.label)

        def _ensure_gf_ready_for_command(self, title: str) -> bool:
            if self.controller.xyzin is None:
                QMessageBox.warning(self, title, "Open or preprocess an ORACLE xyzin first.")
                return False
            if self.process is not None:
                QMessageBox.information(self, "ORACLE", "A command is already running.")
                return False
            self.gf_controller.set_xyzin(self.controller.xyzin)
            return True

        def _ensure_gf_inputs_ready(self, title: str) -> bool:
            if self.controller.xyzin is None:
                return False
            state = load_oracle_project_state(self.controller.xyzin)
            sections = set(state.section_names)
            missing = []
            if "GIC" not in sections:
                missing.append("#GIC")
            if not self.gf_fchk_path.text().strip() and "CARTESIAN_HESSIAN" not in sections:
                missing.append("#CARTESIAN_HESSIAN or FCHK")
            if not missing:
                return True
            QMessageBox.warning(self, title, "Missing " + ", ".join(missing))
            return False

        def _optional_float_field(
            self,
            field: QLineEdit,
            title: str,
            label: str,
            *,
            default=None,
        ):
            text = field.text().strip()
            if not text:
                return default
            try:
                return float(text)
            except ValueError:
                QMessageBox.warning(self, title, f"Invalid {label}: {text}")
                return False

        def _optional_path_text(self, field: QLineEdit) -> Path | None:
            text = field.text().strip()
            return Path(text) if text else None

        def _scale_records(self) -> tuple[str, ...]:
            text = self.gf_scale_records.text().strip()
            if not text:
                return ()
            return tuple(item.strip() for item in text.split(";") if item.strip())

        def _gf_output_path(self, field: QLineEdit, default_factory) -> Path:
            if self.controller.xyzin is None:
                raise ValueError("no ORACLE xyzin project is loaded")
            text = field.text().strip()
            output = Path(text) if text else default_factory(self.controller.xyzin)
            field.setText(str(output))
            return output

        def _set_default_gf_outputs(self, xyzin: Path) -> None:
            self.gf_report_output.setText(str(default_gf_report_output(xyzin)))
            self.gf_csv_dir.setText(str(default_gf_csv_dir(xyzin)))

        def _populate_gf_tables(self, xyzin: Path) -> None:
            state = load_gf_gui_state(xyzin)
            self.gf_summary.setPlainText("\n".join(gf_gui_state_lines(state)))
            self._fill_table(self.gf_frequency_table, state.frequencies)
            self._fill_table(self.gf_gic_table, state.gics)
            self._fill_table(self.gf_ped_table, state.ped)
            self._fill_table(self.gf_diagnostics_table, state.diagnostics)

        def _clear_gf_tables(self) -> None:
            summary = getattr(self, "gf_summary", None)
            if summary is not None:
                summary.clear()
            for table in (
                getattr(self, "gf_frequency_table", None),
                getattr(self, "gf_gic_table", None),
                getattr(self, "gf_ped_table", None),
                getattr(self, "gf_diagnostics_table", None),
            ):
                if table is not None:
                    table.setRowCount(0)

        def _build_sefit_tab(self) -> QWidget:
            tab = QWidget()
            layout = QVBoxLayout(tab)

            job_row = QHBoxLayout()
            job_row.addWidget(QLabel("Job"))
            self.sefit_job_path = QLineEdit()
            job_browse = QPushButton("Browse")
            job_browse.clicked.connect(self.browse_sefit_job)
            job_row.addWidget(self.sefit_job_path, stretch=1)
            job_row.addWidget(job_browse)
            layout.addLayout(job_row)

            out_row = QHBoxLayout()
            out_row.addWidget(QLabel("Outdir"))
            self.sefit_outdir = QLineEdit()
            out_browse = QPushButton("Browse")
            out_browse.clicked.connect(self.browse_sefit_outdir)
            self.sefit_backend = QComboBox()
            self.sefit_backend.addItems(("python", "fortran77"))
            self.sefit_write_section = QCheckBox("Update #MORPHEUS")
            self.sefit_write_section.setChecked(True)
            out_row.addWidget(self.sefit_outdir, stretch=1)
            out_row.addWidget(out_browse)
            out_row.addWidget(QLabel("Backend"))
            out_row.addWidget(self.sefit_backend)
            out_row.addWidget(self.sefit_write_section)
            layout.addLayout(out_row)

            option_row = QHBoxLayout()
            self.sefit_extra_args = QLineEdit()
            self.sefit_extra_args.setPlaceholderText("--max-iter 20 --coordinate-model gic")
            run_button = QPushButton("Run SEFit")
            run_button.clicked.connect(self.run_sefit)
            option_row.addWidget(QLabel("Options"))
            option_row.addWidget(self.sefit_extra_args, stretch=1)
            option_row.addWidget(run_button)
            layout.addLayout(option_row)

            self.sefit_summary = QTextEdit()
            self.sefit_summary.setReadOnly(True)
            self.sefit_summary.setMaximumHeight(150)
            layout.addWidget(self.sefit_summary)

            self.sefit_table_tabs = QTabWidget()
            self.sefit_isotopologue_table = QTableWidget(0, 6)
            self.sefit_output_table = QTableWidget(0, 2)
            self.sefit_diagnostics_table = QTableWidget(0, 2)
            self.sefit_table_tabs.addTab(self.sefit_isotopologue_table, "Isotopologues")
            self.sefit_table_tabs.addTab(self.sefit_output_table, "Outputs")
            self.sefit_table_tabs.addTab(self.sefit_diagnostics_table, "Diagnostics")
            layout.addWidget(self.sefit_table_tabs, stretch=1)
            return tab

        def browse_sefit_job(self) -> None:
            path, _selected = QFileDialog.getOpenFileName(
                self,
                "Select SEFit / MORPHEUS job",
                str(Path.cwd()),
                "SEFit jobs (*.toml *.mfit *.msr *.inp);;All files (*)",
            )
            if path:
                self.sefit_job_path.setText(path)

        def browse_sefit_outdir(self) -> None:
            path = QFileDialog.getExistingDirectory(
                self,
                "Select SEFit output directory",
                self.sefit_outdir.text().strip() or str(Path.cwd()),
            )
            if path:
                self.sefit_outdir.setText(path)

        def run_sefit(self) -> None:
            if not self._ensure_sefit_ready_for_command("SEFit / MORPHEUS"):
                return
            job_text = self.sefit_job_path.text().strip()
            if not job_text:
                QMessageBox.warning(self, "SEFit / MORPHEUS", "Select a SEFit job first.")
                return
            self.sefit_controller.set_xyzin(self.controller.xyzin)
            command = self.sefit_controller.run_command(
                job=Path(job_text),
                outdir=self._sefit_outdir(),
                backend=self.sefit_backend.currentText(),
                write_section=self.sefit_write_section.isChecked(),
                extra_args=self._sefit_extra_args(),
            )
            self._start_command(command, command.label)

        def _ensure_sefit_ready_for_command(self, title: str) -> bool:
            if self.controller.xyzin is None:
                QMessageBox.warning(self, title, "Open or preprocess an ORACLE xyzin first.")
                return False
            if self.process is not None:
                QMessageBox.information(self, "ORACLE", "A command is already running.")
                return False
            self.sefit_controller.set_xyzin(self.controller.xyzin)
            return True

        def _sefit_outdir(self) -> Path:
            if self.controller.xyzin is None:
                raise ValueError("no ORACLE xyzin project is loaded")
            text = self.sefit_outdir.text().strip()
            outdir = Path(text) if text else default_sefit_outdir(self.controller.xyzin)
            self.sefit_outdir.setText(str(outdir))
            return outdir

        def _sefit_extra_args(self) -> tuple[str, ...]:
            import shlex

            return tuple(shlex.split(self.sefit_extra_args.text().strip()))

        def _set_default_sefit_outputs(self, xyzin: Path) -> None:
            self.sefit_outdir.setText(str(default_sefit_outdir(xyzin)))

        def _populate_sefit_tables(self, xyzin: Path) -> None:
            state = load_sefit_gui_state(xyzin)
            self.sefit_summary.setPlainText("\n".join(sefit_gui_state_lines(state)))
            self._fill_table(self.sefit_isotopologue_table, state.isotopologues)
            self._fill_table(self.sefit_output_table, state.outputs)
            self._fill_table(self.sefit_diagnostics_table, state.diagnostics)

        def _clear_sefit_tables(self) -> None:
            summary = getattr(self, "sefit_summary", None)
            if summary is not None:
                summary.clear()
            for table in (
                getattr(self, "sefit_isotopologue_table", None),
                getattr(self, "sefit_output_table", None),
                getattr(self, "sefit_diagnostics_table", None),
            ):
                if table is not None:
                    table.setRowCount(0)

        def _build_trinity_tab(self) -> QWidget:
            tab = QWidget()
            layout = QVBoxLayout(tab)

            engine_row = QHBoxLayout()
            engine_row.addWidget(QLabel("Engine"))
            self.trinity_engine_command = QLineEdit()
            self.trinity_engine_command.setPlaceholderText("external-qm --input step.xyz --gradient gradient.json")
            engine_row.addWidget(self.trinity_engine_command, stretch=1)
            layout.addLayout(engine_row)

            out_row = QHBoxLayout()
            out_row.addWidget(QLabel("Run dir"))
            self.trinity_run_dir = QLineEdit()
            out_browse = QPushButton("Browse")
            out_browse.clicked.connect(self.browse_trinity_run_dir)
            self.trinity_coordinate_model = QComboBox()
            self.trinity_coordinate_model.addItems(("gic", "cartesian"))
            self.trinity_active_space = QComboBox()
            self.trinity_active_space.addItems(("total_symmetric", "all", "cartesian"))
            out_row.addWidget(self.trinity_run_dir, stretch=1)
            out_row.addWidget(out_browse)
            out_row.addWidget(QLabel("Coordinates"))
            out_row.addWidget(self.trinity_coordinate_model)
            out_row.addWidget(QLabel("Active"))
            out_row.addWidget(self.trinity_active_space)
            layout.addLayout(out_row)

            numeric_row = QHBoxLayout()
            self.trinity_max_steps = QLineEdit("50")
            self.trinity_trust_radius = QLineEdit("0.2")
            self.trinity_gradient_tolerance = QLineEdit("1e-5")
            self.trinity_step_tolerance = QLineEdit("1e-5")
            self.trinity_energy_tolerance = QLineEdit("1e-8")
            prepare_button = QPushButton("Prepare TRINITY")
            prepare_button.clicked.connect(self.run_trinity_prepare)
            numeric_row.addWidget(QLabel("Steps"))
            numeric_row.addWidget(self.trinity_max_steps)
            numeric_row.addWidget(QLabel("Trust"))
            numeric_row.addWidget(self.trinity_trust_radius)
            numeric_row.addWidget(QLabel("Grad tol"))
            numeric_row.addWidget(self.trinity_gradient_tolerance)
            numeric_row.addWidget(QLabel("Step tol"))
            numeric_row.addWidget(self.trinity_step_tolerance)
            numeric_row.addWidget(QLabel("Energy tol"))
            numeric_row.addWidget(self.trinity_energy_tolerance)
            numeric_row.addWidget(prepare_button)
            layout.addLayout(numeric_row)

            self.trinity_summary = QTextEdit()
            self.trinity_summary.setReadOnly(True)
            self.trinity_summary.setMaximumHeight(150)
            layout.addWidget(self.trinity_summary)

            self.trinity_table_tabs = QTabWidget()
            self.trinity_settings_table = QTableWidget(0, 2)
            self.trinity_output_table = QTableWidget(0, 2)
            self.trinity_table_tabs.addTab(self.trinity_settings_table, "Settings")
            self.trinity_table_tabs.addTab(self.trinity_output_table, "Outputs")
            layout.addWidget(self.trinity_table_tabs, stretch=1)
            return tab

        def browse_trinity_run_dir(self) -> None:
            path = QFileDialog.getExistingDirectory(
                self,
                "Select TRINITY run directory",
                self.trinity_run_dir.text().strip() or str(Path.cwd()),
            )
            if path:
                self.trinity_run_dir.setText(path)

        def run_trinity_prepare(self) -> None:
            if not self._ensure_trinity_ready_for_command("TRINITY"):
                return
            engine_command = self.trinity_engine_command.text().strip()
            if not engine_command:
                QMessageBox.warning(self, "TRINITY", "Enter an external engine command first.")
                return
            max_steps = self._required_positive_int_field(
                self.trinity_max_steps,
                "TRINITY",
                "max steps",
            )
            if max_steps is None:
                return
            trust_radius = self._optional_float_field(
                self.trinity_trust_radius,
                "TRINITY",
                "trust radius",
                default=0.2,
            )
            if trust_radius is False:
                return
            gradient_tolerance = self._optional_float_field(
                self.trinity_gradient_tolerance,
                "TRINITY",
                "gradient tolerance",
                default=1.0e-5,
            )
            if gradient_tolerance is False:
                return
            step_tolerance = self._optional_float_field(
                self.trinity_step_tolerance,
                "TRINITY",
                "step tolerance",
                default=1.0e-5,
            )
            if step_tolerance is False:
                return
            energy_tolerance = self._optional_float_field(
                self.trinity_energy_tolerance,
                "TRINITY",
                "energy tolerance",
                default=1.0e-8,
            )
            if energy_tolerance is False:
                return
            self.trinity_controller.set_xyzin(self.controller.xyzin)
            command = self.trinity_controller.prepare_command(
                run_dir=self._trinity_run_dir(),
                engine_command=engine_command,
                coordinate_model=self.trinity_coordinate_model.currentText(),
                active_space=self.trinity_active_space.currentText(),
                max_steps=max_steps,
                trust_radius=trust_radius,
                gradient_tolerance=gradient_tolerance,
                step_tolerance=step_tolerance,
                energy_tolerance=energy_tolerance,
            )
            self._start_command(command, command.label)

        def _ensure_trinity_ready_for_command(self, title: str) -> bool:
            if self.controller.xyzin is None:
                QMessageBox.warning(self, title, "Open or preprocess an ORACLE xyzin first.")
                return False
            if self.process is not None:
                QMessageBox.information(self, "ORACLE", "A command is already running.")
                return False
            self.trinity_controller.set_xyzin(self.controller.xyzin)
            return True

        def _trinity_run_dir(self) -> Path:
            if self.controller.xyzin is None:
                raise ValueError("no ORACLE xyzin project is loaded")
            text = self.trinity_run_dir.text().strip()
            run_dir = Path(text) if text else default_trinity_run_dir(self.controller.xyzin)
            self.trinity_run_dir.setText(str(run_dir))
            return run_dir

        def _required_positive_int_field(
            self,
            field: QLineEdit,
            title: str,
            label: str,
        ) -> int | None:
            text = field.text().strip()
            try:
                value = int(text)
            except ValueError:
                QMessageBox.warning(self, title, f"Invalid {label}: {text}")
                return None
            if value <= 0:
                QMessageBox.warning(self, title, f"{label} must be positive.")
                return None
            return value

        def _set_default_trinity_outputs(self, xyzin: Path) -> None:
            self.trinity_run_dir.setText(str(default_trinity_run_dir(xyzin)))

        def _populate_trinity_tables(self, xyzin: Path) -> None:
            state = load_trinity_gui_state(xyzin)
            self.trinity_summary.setPlainText("\n".join(trinity_gui_state_lines(state)))
            self._fill_table(self.trinity_settings_table, state.settings)
            self._fill_table(self.trinity_output_table, state.outputs)

        def _clear_trinity_tables(self) -> None:
            summary = getattr(self, "trinity_summary", None)
            if summary is not None:
                summary.clear()
            for table in (
                getattr(self, "trinity_settings_table", None),
                getattr(self, "trinity_output_table", None),
            ):
                if table is not None:
                    table.setRowCount(0)

        def _fill_table(
            self,
            widget: QTableWidget,
            table: (
                StructureTable
                | GICForgeTable
                | GFTable
                | SEFitTable
                | TrinityTable
                | ToolContractTable
                | WorkbenchTable
            ),
        ) -> None:
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
