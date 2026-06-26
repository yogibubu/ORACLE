# Roadmap

This file records planned development items so the project can be resumed
without reconstructing decisions from the tetrose manuscript work.

## Immediate Next Steps

1. Refactor the current single script into modules before extending the
   graphical interface further. This avoids duplicating command-line logic
   inside the interface.
2. Promote the current shell smoke test into a formal regression-test suite with
   small analytic one- and two-dimensional potentials.
3. Extend the first Tkinter launcher into a fuller project interface:
   - save and reload calculation presets;
   - inspect generated CSV summaries and figures;
   - expose convergence scans;
   - validate incompatible options before launch.

## Methodological Extensions

1. Assess whether to add the small Podolsky pseudopotential term as an optional
   correction for coordinate-dependent 2D metrics. The current implementation
   reads or constructs `g_ij(q1,q2)`, smooths it, monitors positive
   definiteness, and uses the Hermitian `p_i g_ij p_j` kinetic operator.
2. Extend the 1D property workflow conventions to 2D: derivative-defined
   properties, symmetry/parity handling, and finite-temperature averages with a
   coordinate-dependent metric.
3. Extend the derivative-only 1D input beyond Gaussian anharmonic outputs:
   direct stationary-point cards for minima, transition states, energies,
   reduced-coordinate derivatives, and optional reduced masses for converting
   unreduced curvatures.
4. Generalize the product-basis machinery from two dimensions to an arbitrary
   number of coordinates.
5. Keep Gaussian input/path generation in ORACLE/`oracle`; the DVR backend
   should consume completed Gaussian outputs for levels, profiles, and
   properties.
6. Add native puckering-coordinate builders for six- and seven-membered rings.
   Six-membered rings motivate two-dimensional puckering Hamiltonians; seven-
   membered rings motivate the later three-dimensional version.
7. Extend the pruning strategy to higher dimensions by solving the one-
   dimensional problems first, retaining compact axis bases, and then pruning
   the product basis from low-state coefficient weights.

## Graphical Interface Notes

The current `scripts/mw_path_dvr_gui.py` launcher is intentionally a thin layer
over the CLI. A fuller interface should still be a thin layer over stable
library functions, not a separate implementation. The required screens are:

- project setup and file selection;
- Gaussian input preparation;
- Gaussian log analysis;
- one-dimensional solver setup;
- two-dimensional solver setup;
- convergence/pruning setup;
- output viewer for CSV summaries and generated figures.

The command-line interface should remain available as the reproducible backend.
Every graphical action should be exportable as the equivalent command-line
command.

## Documentation Requirements

Before each release or shared zip, check that the repository contains:

- a current `README.md`;
- the full user manual in `docs/mw_path_dvr_manual.md`;
- this roadmap;
- `docs/DEVELOPMENT_STATE.md` with current limitations and validation commands;
- example inputs that can be run without access to the tetrose manuscript
  folder.
