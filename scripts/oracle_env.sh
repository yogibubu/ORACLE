# Source this file from bash or zsh to install ORACLE shell helpers.
#
# Example:
#   source /Users/vincenzobarone/ORACLE/scripts/oracle_env.sh
#   oracle-set

if [ -n "${BASH_SOURCE:-}" ]; then
    _oracle_env_file="${BASH_SOURCE[0]}"
elif [ -n "${ZSH_VERSION:-}" ]; then
    _oracle_env_file="${(%):-%x}"
else
    _oracle_env_file="$0"
fi

_oracle_env_dir="$(cd "$(dirname "$_oracle_env_file")/.." 2>/dev/null && pwd)"

export ORACLE_HOME="${ORACLE_HOME:-${_oracle_env_dir:-$HOME/ORACLE}}"
export ORACLE_VENV="${ORACLE_VENV:-$HOME/.venvs/oracle}"
export ORACLE_CONDA_ENV="${ORACLE_CONDA_ENV:-oracle_26}"
export ORACLE_PYTHON="${ORACLE_PYTHON:-python3}"
export ORACLE_AUTO_CREATE_VENV="${ORACLE_AUTO_CREATE_VENV:-1}"
export ORACLE_AUTO_INSTALL_RUNTIME_DEPS="${ORACLE_AUTO_INSTALL_RUNTIME_DEPS:-1}"
export ORACLE_AUTO_INSTALL_GUI_DEPS="${ORACLE_AUTO_INSTALL_GUI_DEPS:-0}"
export ORACLE_RUNTIME_DEPS="${ORACLE_RUNTIME_DEPS:-numpy scipy matplotlib pytest rdkit}"
export ORACLE_GUI_DEPS="${ORACLE_GUI_DEPS:-PySide6 pytest-qt}"

oracle-package-path() {
    find "$ORACLE_HOME/packages" -mindepth 2 -maxdepth 2 -type d -name src 2>/dev/null \
        | sort \
        | awk '
            BEGIN { first = 1 }
            { printf "%s%s", first ? "" : ":", $0; first = 0 }
            END { print "" }
        '
}

oracle-ensure-gui-deps() {
    if [ "$ORACLE_AUTO_INSTALL_GUI_DEPS" = "0" ]; then
        return 0
    fi
    python - <<'PY' >/dev/null 2>&1
import matplotlib
import scipy
import PySide6
PY
    if [ $? -eq 0 ]; then
        return 0
    fi
    echo "Dipendenze ORACLE GUI/DVR mancanti: installo $ORACLE_GUI_DEPS..."
    python -m pip install $ORACLE_GUI_DEPS
}

oracle-create-venv() {
    if [ "$ORACLE_AUTO_CREATE_VENV" = "0" ]; then
        return 1
    fi
    if [ -z "$ORACLE_VENV" ]; then
        echo "ORACLE_VENV non impostato: impossibile creare il virtualenv."
        return 1
    fi
    if ! command -v "$ORACLE_PYTHON" >/dev/null 2>&1; then
        echo "Python per ORACLE non trovato: $ORACLE_PYTHON"
        return 1
    fi
    echo "Creo virtualenv ORACLE: $ORACLE_VENV"
    mkdir -p "$(dirname "$ORACLE_VENV")" || return
    "$ORACLE_PYTHON" -m venv "$ORACLE_VENV"
}

oracle-activate-venv() {
    if [ -n "$ORACLE_VENV" ] && [ -f "$ORACLE_VENV/bin/activate" ]; then
        # shellcheck disable=SC1091
        source "$ORACLE_VENV/bin/activate"
        return 0
    fi
    if [ -f "$ORACLE_HOME/.venv/bin/activate" ]; then
        # shellcheck disable=SC1091
        source "$ORACLE_HOME/.venv/bin/activate"
        return 0
    fi
    if oracle-create-venv && [ -f "$ORACLE_VENV/bin/activate" ]; then
        # shellcheck disable=SC1091
        source "$ORACLE_VENV/bin/activate"
        return 0
    fi
    return 1
}

oracle-ensure-runtime-deps() {
    if [ "$ORACLE_AUTO_INSTALL_RUNTIME_DEPS" = "0" ]; then
        return 0
    fi
    python - <<'PY' >/dev/null 2>&1
import matplotlib
import numpy
import pytest
import rdkit
import scipy
PY
    if [ $? -eq 0 ]; then
        return 0
    fi
    echo "Dipendenze ORACLE runtime mancanti: installo $ORACLE_RUNTIME_DEPS..."
    python -m pip install --upgrade pip setuptools wheel || return
    python -m pip install $ORACLE_RUNTIME_DEPS
}

oracle-conda-hook() {
    if ! command -v conda >/dev/null 2>&1; then
        return 1
    fi
    local conda_base
    conda_base="$(conda info --base 2>/dev/null)" || return 1
    if [ -f "$conda_base/etc/profile.d/conda.sh" ]; then
        # shellcheck disable=SC1090
        source "$conda_base/etc/profile.d/conda.sh"
        return 0
    fi
    conda activate --help >/dev/null 2>&1
}

oracle-set() {
    if oracle-activate-venv; then
        :
    elif oracle-conda-hook && conda env list | awk '{print $1}' | grep -qx "$ORACLE_CONDA_ENV"; then
        conda activate "$ORACLE_CONDA_ENV"
    elif oracle-conda-hook && conda env list | awk '{print $1}' | grep -qx "oracle"; then
        conda activate oracle
    else
        echo "Warning: ORACLE environment not found, continuing with current Python."
    fi

    local package_path
    package_path="$(oracle-package-path)"
    if [ -n "$package_path" ]; then
        case ":${PYTHONPATH:-}:" in
            *":$package_path:"*) ;;
            *) export PYTHONPATH="$package_path${PYTHONPATH:+:$PYTHONPATH}" ;;
        esac
    fi

    export PATH="$ORACLE_HOME/tools:$ORACLE_HOME/scripts:$ORACLE_HOME/bin:$PATH"
    oracle-ensure-runtime-deps || return
    oracle-ensure-gui-deps || return
    cd "$ORACLE_HOME" || return
}

oracle-cli() {
    oracle-set || return
    python -m oracle "$@"
}

oracle-run() {
    oracle-set || return
    if python - <<'PY' >/dev/null 2>&1
import importlib.util
raise SystemExit(0 if importlib.util.find_spec("oracle_gui.app") else 1)
PY
    then
        python -m oracle_gui.app "$@"
    else
        python -m oracle "$@"
    fi
}

oracle-run-bg() {
    oracle-set || return
    local log="${1:-$ORACLE_HOME/logs/oracle_run.log}"
    mkdir -p "$(dirname "$log")"
    if python - <<'PY' >/dev/null 2>&1
import importlib.util
raise SystemExit(0 if importlib.util.find_spec("oracle_gui.app") else 1)
PY
    then
        nohup python -m oracle_gui.app > "$log" 2>&1 &
    else
        nohup python -m oracle > "$log" 2>&1 &
    fi
    echo "ORACLE avviato in background. PID: $! Log: $log"
}

oracle-test() {
    oracle-set || return
    if [ -z "$1" ]; then
        echo "Usage: oracle-test <test_script.py>"
        return 1
    fi
    python "$1"
}

oracle-test-all() {
    oracle-set || return
    python -m pytest "$ORACLE_HOME/tests" "$ORACLE_HOME/packages"
}

oracle-run-check() {
    oracle-set || return
    python - <<'PY'
import importlib

mods = ["oracle_core", "numpy", "scipy", "matplotlib", "PySide6", "rdkit"]
ok = True
for name in mods:
    try:
        mod = importlib.import_module(name)
        version = getattr(mod, "__version__", "n/a")
        print(f"[OK] {name} {version}")
    except Exception as exc:
        ok = False
        print(f"[FAIL] {name}: {exc}")
if ok:
    print("ORACLE runtime check: OK")
else:
    raise SystemExit(1)
PY
}

oracle-install-gui-deps() {
    ORACLE_AUTO_INSTALL_GUI_DEPS=0 oracle-set || return
    python -m pip install $ORACLE_GUI_DEPS
}

oracle-install-runtime-deps() {
    ORACLE_AUTO_INSTALL_RUNTIME_DEPS=0 oracle-set || return
    python -m pip install --upgrade pip setuptools wheel || return
    python -m pip install $ORACLE_RUNTIME_DEPS
}

oracle-clean() {
    if [ "$1" = "outputs" ] || [ "$1" = "all" ]; then
        echo "Cleaning ORACLE output files outside tests/fixtures (.log .out .err)..."
        find . -path "./tests/fixtures" -prune -o -type f \
            \( -name "*.log" -o -name "*.out" -o -name "*.err" \) -delete
    else
        echo "Cleaning Python cache only..."
    fi

    find . -type d -name "__pycache__" -prune -exec rm -rf {} +
    find . -type f \( -name "*.pyc" -o -name "*.pyo" \) -delete

    echo "Done."
}

oracle-gic-corpus-list() {
    oracle-cli gicforge corpus --format paths --suffix .inp "$@"
}

oracle-gic-corpus-summary() {
    oracle-cli gicforge corpus "$@"
}

oracle-gic-corpus-audit() {
    oracle-cli gicforge corpus-audit "$@"
}

unset _oracle_env_file
unset _oracle_env_dir
