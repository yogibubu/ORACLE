# ORACLE Environment Helpers

ORACLE keeps sourceable shell helpers in the repository instead of editing a
personal shell startup file implicitly. Source the helper once per shell:

```bash
source /Users/vincenzobarone/ORACLE/scripts/oracle_env.sh
oracle-set
```

The helper defines:

- `oracle-set`: create `ORACLE_VENV` when missing, activate the first available
  ORACLE environment, install runtime dependencies including RDKit when needed,
  export `ORACLE_HOME`, add package `src` directories to `PYTHONPATH`, add local
  tools to `PATH`, and `cd` to the repo.
- `oracle-unset`: deactivate the ORACLE virtualenv/conda env when active,
  restore the previous `PATH`, `PYTHONPATH` and working directory saved by
  `oracle-set`.
- `oracle-run`: launch `oracle_gui.app` when it exists, otherwise dispatch to
  `python -m oracle`.
- `oracle-cli`: run `python -m oracle` directly.
- `oracle-run-bg`: launch the same target in background with a log file.
- `oracle-run-check`: verify the current Python can import the core scientific
  runtime stack and report optional external viewers such as Molden.
- `oracle-install-runtime-deps`: explicitly install/upgrade the core runtime
  stack in the active ORACLE environment.
- `oracle-test`, `oracle-test-all`: run focused or full tests in the activated
  environment.
- `oracle-clean`: remove Python cache by default; with `outputs` or `all`, also
  remove ordinary runtime logs while preserving `tests/fixtures`.
- `oracle-gic-corpus-list`: list imported demanding GIC input files.

Default variables use ORACLE names:

```bash
ORACLE_HOME=/Users/vincenzobarone/ORACLE
ORACLE_VENV=$HOME/.venvs/oracle
ORACLE_CONDA_ENV=oracle_26
ORACLE_PYTHON=python3
ORACLE_AUTO_CREATE_VENV=1
ORACLE_AUTO_INSTALL_RUNTIME_DEPS=1
ORACLE_AUTO_INSTALL_GUI_DEPS=0
ORACLE_RUNTIME_DEPS="numpy scipy matplotlib pandas sympy pytest rdkit"
ORACLE_GUI_DEPS="PySide6 pytest-qt"
```

Dense Hermitian diagonalization is routed through
`oracle_core.diagonalizer`. The default policy uses SciPy/NumPy on CPU and can
use CuPy or PyTorch GPU backends when they are already installed in the active
environment. GPU packages are not installed by `oracle-set` because the correct
wheel depends on the local hardware and driver stack.

Useful controls:

```bash
ORACLE_DIAGONALIZER_BACKEND=auto
ORACLE_DIAGONALIZER_GPU_MIN_SIZE=128
ORACLE_DIAGONALIZER_STRICT_GPU=0
```

Molden is an optional external viewer for the Electronic Spectroscopy GUI. On
macOS it needs a Molden executable in `PATH` and XQuartz for the X11 display.
`oracle-run-check` reports both conditions, but it does not install XQuartz
because the macOS package installer requires an interactive administrator
password.

Set `ORACLE_AUTO_CREATE_VENV=0` to prevent `oracle-set` from creating a
virtualenv. Set `ORACLE_AUTO_INSTALL_RUNTIME_DEPS=0` to prevent automatic
runtime dependency installation. RDKit is part of the runtime dependency set
because ORACLE-Babel uses it for SMILES imports.

Set `ORACLE_AUTO_INSTALL_GUI_DEPS=1` before `oracle-set` only when automatic
GUI dependency installation is desired. For routine development, explicit
dependency installation is safer:

```bash
oracle-install-runtime-deps
oracle-install-gui-deps
```

To make the commands permanent in bash, add only this line to `~/.bashrc`:

```bash
source /Users/vincenzobarone/ORACLE/scripts/oracle_env.sh
```

For zsh, put the same line in `~/.zshrc`.
