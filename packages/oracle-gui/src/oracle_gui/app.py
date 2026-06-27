from __future__ import annotations

import argparse
from pathlib import Path

from .dashboard import DashboardAction, OracleDashboardController
from .project import load_oracle_project_state
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
            self.process: QProcess | None = None
            self.current_actions: tuple[DashboardAction, ...] = ()
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
                self.refresh()

        def refresh(self) -> None:
            self.workflow_list.clear()
            self.section_table.setRowCount(0)
            self.action_table.setRowCount(0)
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
            self.controller.log(f"$ {action.shell_line}")
            self.log_output.setPlainText("\n".join(self.controller.log_lines))
            process = QProcess(self)
            if action.command.cwd is not None:
                process.setWorkingDirectory(str(action.command.cwd))
            process.setProgram(action.command.argv[0])
            process.setArguments(list(action.command.argv[1:]))
            process.readyReadStandardOutput.connect(self._read_process_stdout)
            process.readyReadStandardError.connect(self._read_process_stderr)
            process.errorOccurred.connect(lambda error: self._process_error(action, error))
            process.finished.connect(lambda code, _status: self._process_finished(action, code))
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

        def _process_finished(self, action: DashboardAction, code: int) -> None:
            self.controller.log(f"[exit {int(code)}] {action.label}")
            self.process = None
            self.refresh()

        def _process_error(self, action: DashboardAction, error) -> None:
            self.controller.log(f"[process error {int(error)}] {action.label}")
            if self.process is not None and self.process.state() == QProcess.NotRunning:
                self.process = None
                self.refresh()

    app = QApplication([])
    window = OracleDashboardWindow(initial_xyzin)
    window.show()
    return int(app.exec())


if __name__ == "__main__":
    raise SystemExit(main())
