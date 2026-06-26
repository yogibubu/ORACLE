# ORACLE Environment Helpers

ORACLE keeps Merlino-style shell helpers in the repository instead of editing a
personal shell startup file implicitly. Source the helper once per shell:

```bash
source /Users/vincenzobarone/ORACLE/scripts/oracle_env.sh
oracle-set
```

The helper defines:

- `oracle-set`: activate the first available ORACLE environment, export
  `ORACLE_HOME`, add package `src` directories to `PYTHONPATH`, add local tools
  to `PATH`, and `cd` to the repo.
- `oracle-run`: launch `oracle_gui.app` when it exists, otherwise dispatch to
  `tools/oracle_run.py`.
- `oracle-cli`: run `tools/oracle_run.py` directly.
- `oracle-run-bg`: launch the same target in background with a log file.
- `oracle-run-check`: verify the current Python can import the core scientific
  runtime stack.
- `oracle-test`, `oracle-test-all`: run focused or full tests in the activated
  environment.
- `oracle-clean`: remove Python cache by default; with `outputs` or `all`, also
  remove ordinary runtime logs while preserving `tests/fixtures`.
- `oracle-gic-corpus-list`: list imported demanding GIC input files.

Default variables mirror the Merlino helpers but use ORACLE names:

```bash
ORACLE_HOME=/Users/vincenzobarone/ORACLE
ORACLE_VENV=$HOME/.venvs/oracle
ORACLE_CONDA_ENV=oracle_26
ORACLE_AUTO_INSTALL_GUI_DEPS=0
```

Set `ORACLE_AUTO_INSTALL_GUI_DEPS=1` before `oracle-set` only when automatic
GUI dependency installation is desired. For routine development, explicit
dependency installation is safer:

```bash
oracle-install-gui-deps
```

To make the commands permanent in bash, add only this line to `~/.bashrc`:

```bash
source /Users/vincenzobarone/ORACLE/scripts/oracle_env.sh
```

For zsh, put the same line in `~/.zshrc`.
