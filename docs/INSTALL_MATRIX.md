# MATRIX Installation Notes

MATRIX is the framework name for the refactored repository. During the
transition the checkout path and Python packages still use `oracle-*` names,
while the installation contract already covers the whole MATRIX workflow and ORACLE GUI.

These notes describe a macOS-first scientific workstation setup. Linux systems
use the same Python environment and package layout, but install external
viewers through the local distribution package manager.

## Required Base Tools

Install the command-line and build utilities first:

```bash
xcode-select --install
brew install git python gcc make cmake pkg-config
```

`gcc` provides `gfortran`, which is required for strict Fortran77 backends and
for building Molden from source when no binary package is available.

## MATRIX / ORACLE Environment

From a fresh checkout:

```bash
cd /Users/vincenzobarone/MATRIX
source scripts/matrix_env.sh
matrix-set
matrix-run-check
```

`matrix-set` creates or activates the MATRIX virtual environment, adds package
`src` directories to `PYTHONPATH`, and installs the core runtime dependencies
when missing:

- `numpy`
- `scipy`
- `matplotlib`
- `pandas`
- `sympy`
- `pytest`
- `rdkit`

RDKit is part of the runtime stack because LINK uses it for SMILES
imports. Pandas and SymPy are part of the runtime stack because the vendored
WMS-Rot Hamiltonian engine uses them for line-list tables and Wigner algebra.
GUI dependencies are optional and can be installed explicitly:

```bash
matrix-install-gui-deps
```

This installs the declared GUI stack, including `PySide6`.

## Development Tools

`ruff` is the preferred local linter/formatter gate for Python code. Install it
inside the active ORACLE environment:

```bash
cd /Users/vincenzobarone/MATRIX
source scripts/matrix_env.sh
matrix-set
python -m pip install ruff
```

Then run:

```bash
python -m ruff check .
python -m ruff format --check .
```

## Shared Diagonalizer And GPU Acceleration

GF, DVR, VCI/Davidson, vibro-rotational normal-mode extraction and the local
WMS-Rot rotational engine use the shared MATRIX diagonalizer in
`matrix_core.diagonalizer`. By default it chooses GPU acceleration only when a
supported GPU backend is already available and the matrix is large enough to
justify transfer overhead. Otherwise it uses SciPy/NumPy LAPACK on CPU.

The main controls are:

```bash
export MATRIX_DIAGONALIZER_BACKEND=auto
export MATRIX_DIAGONALIZER_GPU_MIN_SIZE=128
export MATRIX_DIAGONALIZER_STRICT_GPU=0
```

Use `MATRIX_DIAGONALIZER_BACKEND=gpu` only when a GPU package is installed.
Install CuPy or PyTorch manually in the active ORACLE environment because the
right package depends on the workstation:

- NVIDIA/CUDA: install the CuPy or PyTorch build matching the CUDA driver.
- Apple Silicon: install PyTorch with MPS support.
- CPU workstations: no extra package is needed.

## Visualization Programs

MATRIX delegates visualization to external programs. The GUI must not
parse or render private chemistry data when a shared section or external viewer
already owns that task.

| Program | Role | Required by |
| --- | --- | --- |
| Avogadro 2 | 3D structure editing, visual inspection and coordinate cleanup | Structure tab, molecular editing |
| Molden | X11 molecular/orbital/vibration viewer for Molden, Gaussian and related files | Electronic viewer, legacy-style orbital/vibrational inspection |
| MOrbVis browser | Browser-based WebGPU/CPU viewer for Molden and Cube files | Electronic viewer fallback or publication-quality orbital images |
| WMS-Rot local engine | Vendored first-party rotational Hamiltonian and line-list engine | Rotational spectroscopy simulation |
| WMS-Rot browser | Browser/Pyodide rotational spectrum reference workflow | Rotational spectroscopy comparison and WMS-Rot input bridge |
| XQuartz | X11 display server needed by Molden on macOS | Molden GUI |
| WebGPU-capable browser | MOrbVis acceleration; CPU fallback is available when WebGPU is absent | MOrbVis |

### Avogadro 2

Install Avogadro from Homebrew:

```bash
brew install --cask avogadro
```

Homebrew installs `Avogadro2.app`. If the `avogadro2` command is not placed in
`PATH`, either set the GUI executable field to:

```text
/Applications/Avogadro2.app/Contents/MacOS/Avogadro2
```

or create a local command shim:

```bash
ln -sf /Applications/Avogadro2.app/Contents/MacOS/Avogadro2 /opt/homebrew/bin/avogadro2
```

### MOrbVis Browser Fallback

If Molden is not installed, or if XQuartz is not available, use MOrbVis in a
browser for orbital and density visualization from Molden or Cube files.

The ORACLE Electronic GUI exposes a `MOrbVis` button that opens:

```text
https://yasuaki-ito.github.io/morbvis/
```

Because browsers protect local files, ORACLE opens the application URL and the
user loads the `.molden`, `.molden.input`, `.cube` or `.cub` file from inside
the browser UI. MOrbVis uses WebGPU when available and falls back to CPU/Web
Workers when WebGPU is unavailable.

The `Open Selected` button in the Electronic tab chooses MOrbVis automatically
for selected Molden/Cube records when Molden or XQuartz is not available.

The reference paper is stored in:

```text
bibliography/morbvis-browser-based-molecular-orbital-visualization-with-webgpu-accelerated-on-the-fly-evaluation.pdf
```

### WMS-Rot Local Engine And Browser Reference

MATRIX includes a first-party WMS-Rot snapshot in:

```text
external/wmsrot-site/
```

The callable Python Hamiltonian engine is vendored as:

```text
packages/matrix-rovib/src/matrix_rovib/vendor/wmsrot_engine.py
```

It requires the same scientific Python stack as the browser Pyodide engine:

```bash
python -m pip install numpy pandas sympy matplotlib
```

Run the local engine from normalized `xyzin` data with:

```bash
matrix rovib wmsrot-run molecule.xyzin --out molecule.rotational.csv
```

The Rotational Spectroscopy workbench can also open the WMS-Rot browser
application:

```text
https://www.skies-village.it/webtools/wmsrot/
```

Generate a compatible input file from the MATRIX container with:

```bash
matrix rovib wmsrot-input molecule.xyzin --out molecule.wmsrot.txt
```

WMS-Rot browser code is treated as a reference and compatibility target.
Production MATRIX rotational spectroscopy calls the vendored WMS-Rot Python
engine through `matrix-rovib` over shared `xyzin` sections.

### Molden On macOS

Molden is not currently available from the default Homebrew or conda channels
on this Apple Silicon workstation. Build it from the official source tarball:

```bash
brew install gcc libx11 libxmu
mkdir -p /Users/vincenzobarone/MATRIX/.external/molden-build
cd /Users/vincenzobarone/MATRIX/.external/molden-build
curl -L -o molden7.3.tar.gz https://ftp.science.ru.nl/Molden/molden7.3.tar.gz
tar -xzf molden7.3.tar.gz
cd molden7.3
make clean
rm -f src/atomdens.o
make molden \
  FFLAGS='-O3 -funroll-loops -DDARWIN -fallow-argument-mismatch' \
  CFLAGS='-g -std=gnu90 -Wno-return-type -Wno-error=implicit-function-declaration -Wno-implicit-function-declaration -DDARWIN -I/opt/homebrew/include -DDOBACK -DHASTIMER -DCOLOR_OFFSET=0.0 -Wno-logical-op-parentheses -Wno-tautological-pointer-compare -Wno-tautological-constant-out-of-range-compare' \
  LIBS='-L/opt/homebrew/lib -lX11 -lm' \
  LDR=gfortran
cp bin/molden /opt/homebrew/bin/molden
chmod 755 /opt/homebrew/bin/molden
molden -h
```

The `rm -f src/atomdens.o` step is intentional: the upstream tarball can carry
a precompiled Linux/x86-64 object that cannot be linked into a macOS arm64
binary.

Molden needs XQuartz to open windows on macOS.

## XQuartz On macOS

Install XQuartz from Terminal:

```bash
brew install --cask xquartz
```

This runs the official macOS package installer and may ask for the administrator
password. A non-interactive shell cannot provide that password; run the command
yourself in Terminal.app, not from an automated Codex command.

Alternative manual installation:

1. Open `https://www.xquartz.org/releases/`.
2. Download the current stable `.dmg`.
3. Open the `.dmg`.
4. Double-click `XQuartz.pkg`.
5. Follow the installer prompts and enter the administrator password.
6. Log out and log back in, or reboot.
7. Start XQuartz once from `/Applications/Utilities/XQuartz.app`.

Verify after logging back in:

```bash
test -d /Applications/Utilities/XQuartz.app && echo "XQuartz app installed"
test -x /opt/X11/bin/Xquartz && echo "XQuartz server installed"
open -a XQuartz
molden -h
```

If `molden molecule.molden` reports that it cannot open the display, log out
and back in again so XQuartz startup files and `DISPLAY` handling are loaded.

## QM And External Calculation Utilities

These tools are optional but used by specific MATRIX workflow and ORACLE GUIs:

| Utility | Role | Notes |
| --- | --- | --- |
| Gaussian executable (`g16`, `gdv`, or site wrapper) | QM jobs, Hessians, FCHK, rovibrational logs | Configure in the QM Jobs GUI or CLI `--executable` |
| `formchk` | Convert Gaussian checkpoint files to FCHK | Used by `matrix gaussian formchk` |
| Molpro | External QM job launcher and output source | Launched by `matrix molpro run`; parsed only by `matrix-molpro` adapters |
| ORCA | External QM job launcher and output source | Launched by `matrix orca run`; `matrix orca promote` writes final geometry and `#CARTESIAN_HESSIAN` when ORCA prints a Cartesian Hessian |
| MRCC | External QM output source | Parsed only by `matrix-mrcc` adapters |
| Browser | Opens MOrbVis and future HTML reports | Prefer Chrome/Edge/Safari versions with WebGPU support |

All QM output parsing must go through the single adapter for that external
format. Scientific tools consume normalized `xyzin` sections and must not parse
Gaussian, Molpro, ORCA or MRCC output privately.

For the private `oracle` Linux host, MATRIX provides SSH/SCP wrappers around
the remote `~/matrix/bin/matrix-submit` helper:

```bash
matrix qm remote-submit job.gjf --engine gdv32 --host enzo@oracle
matrix qm remote-status --host enzo@oracle
matrix qm remote-fetch JOB_NAME --host enzo@oracle --dest runs
matrix qm remote-fetch JOB_NAME --host enzo@oracle --dest runs --promote orca --xyzin molecule.xyzin
```

`remote-fetch` copies the `native_output` recorded in the remote `metadata.txt`
and writes a local `remote_qm_fetch_manifest.json`.  Promotion into `xyzin`
remains adapter-based: `molpro`, `orca` and the explicit Gaussian modes
(`gaussian-log-hessian`, `gaussian-rovib`, `gaussian-electronic`,
`gaussian-fchk`) are the supported promotion choices.

## Post-install Checks

Run:

```bash
source /Users/vincenzobarone/MATRIX/scripts/matrix_env.sh
matrix-set
matrix-run-check
PYTHONPATH=. pytest -q
```

For viewer checks:

```bash
command -v avogadro2 || echo "Use /Applications/Avogadro2.app/Contents/MacOS/Avogadro2 in the GUI"
command -v molden || echo "Molden missing; use MOrbVis browser for Molden/Cube orbital files"
open -a XQuartz
python -m webbrowser -t https://yasuaki-ito.github.io/morbvis/
```
