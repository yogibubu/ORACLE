# Development State

Last updated: 2026-06-20.

This repository is the standalone development area for the puckering/path
Hamiltonian code. It is intentionally separated from the tetrose manuscript
folder so that the code, examples, and developer documentation can be versioned
and shared independently.

## Entry Point

The main program is:

```bash
python3 scripts/mw_path_dvr.py --help
```

The code currently lives in one script to keep distribution simple. The most
likely future split is:

- Gaussian log parsing;
- ring-puckering diagnostic label builders;
- one-dimensional solvers;
- two-dimensional product-basis solvers;
- command-line and graphical user interfaces.

Earlier tetrahydrofuran benchmark helpers are kept in the frozen ORACLE
tree. New development should use `scripts/mw_path_dvr.py` as the reference
implementation.

## Implemented Features

- Gaussian log reading for scan geometries, energies, and optional properties.
- Gaussian log reading for any optimized scan/path; the scanned coordinate does
  not have to be puckering.
- Four-membered-ring puckering scans with one out-of-plane atom.
- Five-membered-ring puckering scans using endocyclic-dihedral combinations.
- Five-membered-ring labels such as `E1`, `E2`, `T12`, and related variants.
- Cremer-Pople quantities as diagnostic labels, not as the working coordinate.
- One-dimensional periodic calculations with a Fourier basis.
- One-dimensional nonperiodic calculations with either an optimized distributed
  Gaussian basis or a sinc discrete-variable-representation fallback.
- One-dimensional mass-weighted path analysis from Cartesian scan points.
- Mandatory one-dimensional convention: production Hamiltonians use the
  cumulative mass-weighted Cartesian distance as the dynamical coordinate, with
  reduced mass equal to one. Angles, puckering phases, and Cremer-Pople
  quantities are labels or scan generators only. Consecutive geometries are
  first translated and locally Eckart-oriented so that the distance is not
  contaminated by overall rotation.
- Expectation values of scalar properties, including rotational constants when
  they are available from the Gaussian output.
- Symmetry preprocessing for one-dimensional paths:
  `half-even-origin`/`half-even-last` for one-branch scans and `center-even`
  for complete even paths.
- Analytic endpoint continuation for nonperiodic 1D potentials with degree-6
  or degree-8 repulsive polynomial walls, plus an optional double-Morse
  endpoint model.
- Born-Mayer-like exponential endpoint continuation for steep repulsive tails.
- Explicit `--well-type` handling for single- versus double-minimum 1D
  potentials, including direct single-Morse and shifted inverse-power
  continuations for one-well coordinates.
- Explicit symmetry reference choices for half scans: first point,
  last point, or center of a complete path.
- Symmetry preprocessing can be combined with periodic one-dimensional
  Hamiltonians; analytic tails remain restricted to nonperiodic paths.
- Sparse asymmetric double-well core fitting with the Flanigan-de la Vega
  parabola+Gaussian form, kept separate from endpoint tail continuation.
- Flanigan-de la Vega double-well criterion reported for symmetric 1D scans to
  guide parabola+Gaussian versus double-Morse central representations.
- Optional spline interpolation of the symmetrized and extended 1D potential,
  with the full model profile and tail diagnostics written to output files.
- One-dimensional scalar properties from either point values or local
  derivatives in the mass-scaled path coordinate. For symmetrized potentials,
  even-degree derivative properties are treated as symmetric and odd-degree
  properties as antisymmetric, giving the central value for diagonal
  expectations of the odd part.
- State-by-state, 0 K, and finite-temperature Boltzmann averages for
  one-dimensional scalar properties.
- One-mode variational calculations from Gaussian anharmonic outputs via
  `--anharmonic-mode`. The parser reads F2/F3/F4 from the Gaussian
  quadratic/cubic/quartic normal-mode force-constant tables, maps standard
  frequency numbering to force-table numbering, and treats double-well
  references as transition states at `Q=0`.
- Analytic derivative-only one-mode models: Handy Gauss-like coordinates for
  vanishing F3, Handy Morse-like coordinates for one-well asymmetric minima,
  and shifted `1/r + 1/r^2 + 1/r^3` when the F2/F3/F4 real-pole condition is
  satisfied. The Handy transformations preserve the local derivatives through
  fourth order while improving the long-range behavior of normal-coordinate
  cuts. A degree-6/8 confining wall is available without changing derivatives
  through fourth order at `Q=0`.
- One-mode VPT2 comparison for non-double-well references, using the standard
  normal-coordinate perturbative levels from F2/F3/F4 beside the variational
  levels of the selected Handy-Morse or Handy-Gaussian analytic potential.
- Two-dimensional rectangular-grid Hamiltonians with independently periodic or
  nonperiodic axes.
- Constant two-dimensional kinetic metric with optional constant cross term.
- Coordinate-dependent two-dimensional kinetic metrics read from CSV columns.
- Automatic two-dimensional kinetic metric construction from a matching grid of
  molecular structures, using mass-weighted Eckart finite differences and the
  Hermitian `p_i g_ij(q1,q2) p_j` operator.
- Robust geometry-derived metric construction with configurable finite-
  difference stencil, automatic low-order polynomial/Fourier smoothing,
  positive-definite projection, and summary diagnostics.
- Two-dimensional product-basis pruning from one-dimensional reference
  Hamiltonians.
- Energy-window pruning on each one-dimensional axis.
- Product-basis coefficient pruning after an initial two-dimensional
  diagonalization.
- Automated two-dimensional convergence scans over basis sizes, pruning
  potentials, and grid downsampling.
- Lightweight Tkinter GUI launcher in `scripts/mw_path_dvr_gui.py`. It builds
  command lines for Gaussian anharmonic one-mode analyses, 1D paths, and 2D
  grids, then runs the same CLI backend.

## Default Production Procedure

- Use the sampled mass-weighted path as the default potential core.
- Specify `--well-type single` or `--well-type double` explicitly in production
  commands.
- Use `--path-symmetry half-even-origin`, `half-even-last`, or `center-even`
  only when the scan ordering and physics justify that symmetry.
- Keep `--core-model auto`; it remains sampled except for sparse unsymmetrized
  double wells.
- Use `--potential-extension repulsive-exponential` as the first nonperiodic
  repulsive-tail option when endpoint continuation is required, and compare
  degree-8 polynomial tails when the result is sensitive.
- For periodic 1D potentials, symmetry preprocessing is allowed, but analytic
  tails and spline smoothing remain disabled.
- For derivative-only VPT2 comparisons, use `--anharmonic-mode` on the Gaussian
  anharmonic output. Set `--well-type double` only when the reference structure
  is the TS; otherwise use `--well-type single` for a minimum.

## Known Limitations

- Coordinate-dependent two-dimensional metrics neglect the small Podolsky
  pseudopotential term. It can be added later as an optional correction if a
  target application needs it.
- Stationary-point-only 1D fits require curvatures in the mass-weighted
  coordinate or an explicitly supplied reduced mass. For curvilinear paths, a
  sparse set of stationary structures does not by itself determine the path.
- Gaussian input generation is intentionally handled by ORACLE/`oracle`,
  not by this DVR backend.
- Six- and seven-membered-ring puckering coordinates are not yet implemented.
- The GUI is a first practical launcher, not yet a full project manager. It
  does not replace the command line for advanced or scripted production runs.
- Regression tests are currently a compact shell smoke test plus manual
  application checks rather than a formal unit-test suite.

## Minimal Restart Procedure

From a clean checkout:

```bash
python3 -m pip install -r requirements.txt
python3 -m py_compile scripts/mw_path_dvr.py
python3 scripts/mw_path_dvr.py --help
./scripts/smoke_test.sh
```

Generate Gaussian inputs with the ORACLE GUI Advanced puckering workflow or
with `oracle pucker-gaussian`, then analyze the completed Gaussian log with
this DVR backend.

Run a two-dimensional convergence scan once a rectangular CSV grid is available:

```bash
python3 scripts/mw_path_dvr.py \
  --grid2d-csv examples/gaussian_outputs/phi_q_grid.csv \
  --q1-key phi --q2-key q \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --boundary1 periodic --boundary2 nonperiodic \
  --metric-mode constant \
  --solver1 fourier --solver2 gaussian \
  --basis1 16 --basis2 12 \
  --convergence-scan \
  --convergence-basis1 8,12,16 \
  --convergence-basis2 8,12 \
  --convergence-prune-potentials min,mean \
  --convergence-strides 2,1 \
  --levels 10 \
  --prefix phi_q_2d
```

## Files to Read First

1. `README.md` for installation and quick commands.
2. `docs/mw_path_dvr_manual.md` for user-level details and equations.
3. `docs/ROADMAP.md` for the ordered development plan.
4. `scripts/mw_path_dvr.py` for implementation details.

## Current Validation Commands

```bash
python3 -m py_compile scripts/mw_path_dvr.py
python3 scripts/mw_path_dvr.py --help
./scripts/smoke_test.sh
```

The smoke test covers command-line parsing, Gaussian-log path analysis,
one-dimensional nonperiodic symmetry/tail/spline handling, two-dimensional
geometry-derived metrics, two-dimensional CSV metrics, and constant metrics. For
any future code change, add a small reproducible command here or replace this
section with a formal test runner.
