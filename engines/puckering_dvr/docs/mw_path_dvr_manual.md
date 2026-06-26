# Short Manual: `mw_path_dvr.py`

This program treats one-dimensional large-amplitude motion along a Cartesian
geometry path read from Gaussian output. Gaussian input/path generation is not
part of this DVR backend in ORACLE; it is handled by ORACLE/`oracle`.

Program file:

```bash
scripts/mw_path_dvr.py
```

## 1. Scope

The program analyzes completed Gaussian outputs. After the scan has finished,
it reads optimized geometries and energies from the Gaussian log, translates
each structure to its center of mass, reorients consecutive structures with a
local mass-weighted Eckart alignment, constructs the mass-weighted path
coordinate, solves a one-dimensional path Hamiltonian, and computes expectation
values of scalar properties such as rotational constants.

The dynamical coordinate is not the puckering phase itself. The coordinate used
in the Hamiltonian is the cumulative mass-weighted Cartesian distance,

```text
s_i = s_{i-1} + [ sum_a m_a |R_i r_{a,i} - r_{a,i-1}|^2 ]^(1/2)
```

where `R_i` is the local mass-weighted rotation between adjacent geometries.
The rotation is chosen so that the rotational component of the displacement is
removed, equivalently minimizing the angular momentum of the step under the
discrete Eckart condition

```text
sum_a m_a rbar_{a,i} x (R_i r_{a,i} - r_{a,i-1}) = 0
```

with `rbar_{a,i}` the midpoint geometry of the aligned step. Only after this
reorientation is the path increment accumulated. In this coordinate the reduced
mass is one by construction.

This is mandatory for all one-dimensional calculations. The Gaussian scan
coordinate, an angle, a puckering phase, or a Cremer-Pople label may identify
the relaxed structures generated upstream by ORACLE, but it must not replace the
mass-weighted Cartesian path coordinate in the Hamiltonian. If the relaxed
points are not equally spaced in Eckart-reoriented mass-weighted distance, the
program must keep the nonuniform abscissa, interpolate `V(s)`, and solve the
one-dimensional problem with reduced mass equal to one. A constant effective
inertia or a linear relation such as `s = k phi` is acceptable only as a
diagnostic, not as the production one-dimensional vibrational model.

## 2. Ring Coordinate Conventions

The ring atoms must be supplied in cyclic Prelog order through `--ring`.

For five-membered rings:

- `phi = 0` is the envelope in which ring atom 1 is out of the mean ring plane
  (`E1`);
- `phi = 180` is the opposite envelope of the same atom;
- the planar structure has zero puckering amplitude, `q = 0`, and the phase is
  undefined;
- the minimum can occur at any value of `phi`.

For a ring ordered as atoms 1-5, the endocyclic dihedrals are:

```text
tau0 = D(1,2,3,4)
tau1 = D(2,3,4,5)
tau2 = D(3,4,5,1)
tau3 = D(4,5,1,2)
tau4 = D(5,1,2,3)
```

The operational scan coordinate is built from linear combinations of these
dihedrals. Cremer-Pople coordinates are written only as diagnostic labels when
`--label-cremer-pople` is used; they are not used as the scan coordinate and are
not averaged as spectroscopic properties.

For four-membered rings there is only one puckering coordinate in this simple
description: atom 1 moves out of the plane defined by the other three ring
atoms. The Gaussian input uses the improper/endocyclic torsion

```text
Puck = D(2,3,4,1)
```

as the generalized internal coordinate. If a different atom should be the
out-of-plane atom, rotate the `--ring` list so that atom is first. For example,
to make the original ring atom 3 the reference atom, pass a cyclic ordering with
that atom first.

## 3. Preparing Gaussian Input

Prepare the scan in ORACLE, not in this DVR backend. In the GUI use the
Advanced puckering Gaussian workflow; from the command line use
`oracle pucker-gaussian`. ORACLE writes the Gaussian `.gjf` and any scan
manifest. After Gaussian has run, pass the resulting `.log` to
`scripts/mw_path_dvr.py`.

## 6. Running Gaussian

Run a generated input with the local Gaussian command, for example:

```bash
g16 < examples/gaussian_inputs/thf_phi_function_scan.gjf > examples/gaussian_outputs/thf_phi_function_scan.log
```

For erythrose test inputs:

```bash
g16 < examples/gaussian_inputs/erythrose_alpha_E2_dpcs3_phi_scan_10deg.gjf > examples/gaussian_outputs/erythrose_alpha_E2_dpcs3_phi_scan_10deg.log
```

The generated files contain one Gaussian link for each target phase.

## 7. Analyzing a Gaussian Log

The solver is selected with `--solver`.

- `--solver auto`: default. Uses Fourier for periodic paths and an optimized
  distributed Gaussian basis for nonperiodic paths.
- `--solver fourier`: periodic Fourier representation. This is the recommended
  solver for pseudorotations and other cyclic paths.
- `--solver gaussian`: optimized distributed Gaussian basis. This is the
  recommended solver for localized nonperiodic paths.
- `--solver sinc-dvr`: sinc discrete variable representation. This is retained
  as a simple, robust fallback for nonperiodic paths.

For a periodic potential:

```bash
python3 scripts/mw_path_dvr.py \
  --gaussian-log examples/gaussian_outputs/thf_phi_function_scan.log \
  --log-selection last-per-link \
  --boundary periodic \
  --solver fourier \
  --compute-rotconst \
  --label-cremer-pople \
  --prefix thf_phi_mass_weighted
```

For a nonperiodic open path with the distributed Gaussian basis:

```bash
python3 scripts/mw_path_dvr.py \
  --gaussian-log examples/gaussian_outputs/path.log \
  --boundary nonperiodic \
  --solver gaussian \
  --compute-rotconst \
  --prefix path_mass_weighted
```

Important options:

- `--log-selection all`: use all optimized geometries found in a continuous
  Gaussian scan. If the log contains Gaussian relaxed-scan records of the form
  `Step number ... on scan point ...`, the program keeps the last step record
  for each scan point, separates scan blocks when point numbering restarts, and
  uses the corresponding `Input orientation` or `Z-Matrix orientation` when
  available to avoid discontinuities from changing standard-orientation axes.
- `--log-selection last-per-link`: use the last optimized geometry of each
  Gaussian link. Use this for inputs generated by `functional-targets` or the
  multi-link fallback.
- `--gaussian-energy auto`: use post-SCF total energies such as `E(method)` or
  `EUMP2` when they are present, otherwise use `SCF Done`.
- `--gaussian-energy post-scf`: force MP2/double-hybrid total energies and fail
  if any selected structure lacks them.
- `--gaussian-energy scf`: force the `SCF Done` energies.
- `--boundary periodic`: periodic path. Use `--solver fourier`.
- `--boundary nonperiodic`: open path. Use `--solver gaussian` or
  `--solver sinc-dvr`.
- `--periodic-endpoints first`: identify the final periodic point with the
  first point. This is the default.
- `--compute-rotconst`: compute `A`, `B`, and `C` from each Cartesian geometry
  and average them over the vibrational states.
- Without explicit `--property` options, the 1D workflow averages only the
  standard rotational constants `A_MHz`, `B_MHz`, `C_MHz` when available and the
  invariant Gaussian dipole magnitude `dipole_debye` when present. Other scalar
  properties, including dipole components, are included only on request.
- `--property A_MHz --property B_MHz --property C_MHz`: average the three
  rotational constants read from a Gaussian log or produced by
  `--compute-rotconst`. Use `--property dipole_x_debye` etc. to request dipole
  components, or `--property all` to restore the older all-numeric-column
  behavior.
- `--vpt2-property-csv file.csv`: compare DVR ground-state property averages
  with one-dimensional VPT2 values. The CSV must contain `property` and either
  `vpt2_value`/`value` or `vpt2_delta`/`delta`, with optional
  `total_perturbative_value` or `total_perturbative_delta`. The output
  `*_property_comparison.csv` and compatibility copy
  `*_vpt2_property_correction.csv` report both the variational 1D and
  perturbative 1D values plus the variational-minus-VPT2 correction.
- `--property-vpt2 auto`: default. For every selected property not supplied in
  `--vpt2-property-csv`, fit a local Taylor expansion around the 1D minimum and
  estimate the perturbative expectation value with a first-order anharmonic
  wavefunction correction in a harmonic basis. Use `csv-only` to write only CSV
  supplied values or `off` to disable property VPT2 comparisons.
- `--property-vpt2-fit-points`, `--property-vpt2-degree`, and
  `--property-vpt2-basis-size`: control the local Taylor fit and harmonic basis
  used by the automatic perturbative property estimate.
- `--plot-max-state 5`: highest state index used to set the energy window of
  the potential plot. The default is 5.
- `--plot-property-smooth-degree N`: optional Chebyshev degree used only to
  smooth the plotted property curves. The numerical averages still use the
  unsmoothed property grid.
- `--label-cremer-pople`: add `CP_q2_angstrom` and `CP_phi2_deg` diagnostic
  labels to profile outputs.
- `--grid`: number of interpolation/output grid points; default 401.
- `--gaussian-basis-size`: number of distributed Gaussians for
  `--solver gaussian`; default 40.
- `--gaussian-width-scale`: Gaussian exponent scale, defined as
  `alpha = scale/dx^2`; default 1.0.
- `--gaussian-quadrature`: Gauss-Hermite quadrature order for Gaussian matrix
  elements; default 24.
- `--gaussian-optimize`: nonlinear Gaussian-basis optimization. The default is
  `centers-widths`, which optimizes the internal centers and all widths.
  `widths` keeps the centers fixed and optimizes only the exponents. `none`
  uses a fixed Hamilton-Light distributed Gaussian basis.
- `--gaussian-optimization-levels`: number of low-lying states minimized during
  basis optimization. The default value `0` means `min(--levels, basis size)`.
- `--gaussian-optimization-maxiter`: maximum number of L-BFGS-B iterations for
  the Gaussian basis optimization; default 80.
- `--gaussian-overlap-threshold`: relative overlap eigenvalue threshold used to
  remove numerical linear dependencies; default `1e-8`.
- `--gaussian-condition-limit`: maximum condition number retained in the
  Gaussian overlap subspace; default `1e8`.
- `--levels`: number of eigenvalues/eigenvectors to compute; default 12.

### Symmetric Paths and Analytic Tails

Recommended production procedure:

1. Use the sampled mass-weighted Cartesian path as the core potential whenever
   enough relaxed structures are available.
2. Declare the topology explicitly with `--well-type single` or
   `--well-type double`.
3. Apply `--path-symmetry half-even-origin`, `half-even-last`, or `center-even`
   only when that symmetry is physically present.
4. Keep `--core-model auto`; it resolves to `sampled` except for very sparse
   unsymmetrized double wells, where the Flanigan-de la Vega
   parabola+Gaussian core is used.
5. Add endpoint tails only for nonperiodic paths whose computed endpoints are
   not high enough. The robust first choice is
   `--potential-extension repulsive-exponential`; compare polynomial degree 8
   if the tail affects low levels or expectation values.
6. For periodic paths, use symmetry preprocessing if needed, but do not add
   tails or spline smoothing.

For inversion or puckering coordinates sampled on only one branch, the
1D workflow can construct the missing branch before interpolation. This
symmetry preprocessing can be used with either nonperiodic or periodic boundary
conditions:

- `--path-symmetry none`: use the path as supplied.
- `--path-symmetry half-even` or `half-even-origin`: treat the first point as
  the symmetry origin and mirror the positive branch. This is appropriate when
  the computed scan starts at a planar structure or transition state and follows
  one equivalent branch.
- `--path-symmetry half-even-last`: treat the last point as the symmetry origin
  and mirror the branch supplied before it. This is appropriate when the scan is
  ordered from one minimum toward the planar structure or transition state.
- `--path-symmetry center-even`: treat the midpoint of a complete path as the
  symmetry origin and average the two branches.

For periodic paths the symmetrized sample is used as one periodic cell before
Fourier interpolation. Analytic tail extensions and spline smoothing are still
nonperiodic operations and are rejected with `--boundary periodic`.

Use `--well-type auto`, `single`, or `double` to state the topology of the
central potential. `auto` classifies the symmetrized core. Explicit `single` and
`double` are safer for production scans because they prevent using a one-well
analytic model on a double-minimum potential, or a double-Morse central model on
a single-minimum potential.

The central potential and the endpoint tails are separate choices. The default
`--core-model auto` keeps the sampled core in most cases. It switches to the
Flanigan-de la Vega asymmetric parabola+Gaussian form only for sparse,
unsymmetrized double-minimum scans with at most `--core-auto-max-points` points.
Use `--core-model sampled` to force the numerical core, or
`--core-model asymmetric-parabola-gaussian` to force the analytic asymmetric
double-well core. The fitted form is
`c + a2*(s-sp)^2 + Vg*exp[-ag*(s-sg)^2]`, i.e. the 1974
parabola+Gaussian double-well model with an offset and a movable parabolic
origin so the result does not depend on the arbitrary path origin.

If only stationary structures are supplied, the second derivatives must be in
the same mass-weighted coordinate used by the Hamiltonian, or a reduced mass
must be supplied before converting them. With curvilinear coordinates, stationary
structures alone are normally insufficient; provide additional intermediate
structures so the mass-weighted path length and curvature are determined from
the actual relaxed path.

Scalar properties can be averaged in two ways. Pointwise values are read from
the path comments or from `--properties-csv` and selected with repeated
`--property` options. Alternatively, `--property-derivatives-csv` reads local
derivatives in the final mass-scaled path coordinate, using rows such as

```text
property,value,d1,d2,d3,d4,origin,parity
B_MHz,4859.0,,0.25,,0.001,zero,auto
mu_z,1.20,0.05,,0.0002,,zero,auto
```

The derivatives are ordinary derivatives with respect to `s_au`; the expansion
uses `value + d1*q + d2*q^2/2! + ...`, so the path reduced mass remains one and
no separate kinetic-energy parameter is introduced. `origin` can be `zero`,
`first`, `last`, `minimum`, or a numeric coordinate in a.u. For symmetrized
potentials, `parity=auto` expects the non-zero property derivatives to be either
all even or all odd. The even case keeps only even powers; the odd case keeps
only odd powers around the central value. In the odd case the diagonal
expectation value is exactly the central value for every eigenstate of the
symmetric potential, and therefore also for every Boltzmann average. Use
explicit `parity=even`, `odd`, or `full` only when the automatic rule is not
appropriate.

`*_expectations.csv` reports state-by-state diagonal expectations. When scalar
properties are present, `*_thermal_expectations.csv` is also written: it always
contains the 0 K ground-state row and adds any temperatures requested with
`--temperature T` using Boltzmann weights over the computed one-dimensional
levels.

Terminal points can be continued analytically before the interpolation used by
the DVR:

- `--potential-extension repulsive-polynomial`: estimate endpoint slope and
  curvature from the last `--extension-fit-points` points and add a convex
  degree-6 or degree-8 wall. The old `*-quartic` option names are accepted only
  as legacy aliases; quartic walls are not recommended for production.
- `--potential-extension repulsive-exponential`: match endpoint slope and
  curvature, then add a Born-Mayer-like exponential repulsive wall adjusted to
  reach `--extension-target-cm` at the tail endpoint.
- `--potential-extension morse-polynomial`: fit an even double-Morse profile to
  the symmetrized potential, then use the larger of Morse and local endpoint
  slope/curvature for the final polynomial wall.
- `--potential-extension single-morse`: for `--well-type single`, fit a direct
  one-well Morse form and use it to continue the scan.
- `--potential-extension single-inverse-power`: for `--well-type single`, fit a
  shifted Kratzer-like `1/r + 1/r^2` form and continue the scan with that
  inverse-power model.
- For polynomial tails, `--extension-degree 6` is the default production
  choice; `--extension-degree 8` is recommended as a sensitivity check.
- `--extension-length-au`, `--extension-points`, and `--extension-target-cm`
  control the wall length, sampling density, and terminal height.
- `--potential-smoothing spline` applies a cubic `UnivariateSpline` to the final
  extended sample; `--potential-spline-smoothing 0.0` gives an interpolating
  spline.

The endpoint derivatives are obtained from a quadratic least-squares fit
constrained to pass through the terminal point and scaled by the local interval.
This avoids a small endpoint offset in the fit from changing the C2 analytic
tail. The summary file reports fitted slopes, curvatures, polynomial
coefficients, and the minimum first and second outward finite differences of
each tail. The model points used for interpolation are written to
`*_model_profile.csv`. For repulsive tails these finite differences should be
positive. For single-well Morse or inverse-power long-range continuations, a
soft or slightly decreasing asymptotic side can be physical and should be judged
against the intended coordinate.

For symmetric double-minimum potentials, the summary also reports the
Flanigan-de la Vega criterion when `--double-well-criterion auto` is active. The
diagnostic is evaluated on the symmetrized ab initio core before analytic tails
or spline smoothing are applied. With barrier height `E`, distance `D` between
the two equivalent minima, and minimum curvature `k_m`, the reported value is
`D/4 * sqrt(k_m/(2E))`. Values below one favor a parabola+Gaussian
representation, while values above one favor a double-Morse representation.
This is a diagnostic for the central double-well representation, not a criterion
for choosing the endpoint tails. It is useful for fitting sparse scans and
reducing the number of additional central points to compute; the DVR still uses
the selected numerical profile, extension, and smoothing options. Use
`--double-well-fit-points` to change the local radial fit used for `k_m`, or
`--double-well-criterion off` to disable the diagnostic.

### One-Mode Gaussian Anharmonic Derivatives

For comparison with VPT2, the code can read the one-mode derivatives directly
from a Gaussian anharmonic output and solve a variational one-dimensional
problem in the dimensionless normal coordinate:

```bash
python3 scripts/mw_path_dvr.py \
  --gaussian-log examples/gaussian_outputs/alpha_VPT2.log \
  --anharmonic-mode 1 \
  --well-type double \
  --solver sinc-dvr \
  --grid 401 \
  --levels 8 \
  --prefix alpha_mode1_derivative_dvr
```

The parser reads the `QUADRATIC`, `CUBIC`, and `QUARTIC FORCE CONSTANTS IN
NORMAL MODES` tables. The derivative values are the Gaussian reduced constants
in `cm^-1`: `F2`, diagonal `F3(i,i,i)`, and diagonal `F4(i,i,i,i)`. If Gaussian
does not print a diagonal cubic term, it is treated as zero. If it does not
print a quartic term, the summary records that zero was assumed.

By default, `--anharmonic-mode` is the ordinary mode number that should be
treated variationally, using the standard Gaussian `Frequencies --` numbering.
Gaussian's anharmonic force-constant tables can be printed in the opposite
order; the code maps the requested frequency index to the table row by matching
`F2`. Use `--anharmonic-mode-order force-table` only when the requested number
is already the force-constant table index.

For `--well-type double`, the reference point is the transition state: `Q=0` is
the barrier geometry and an imaginary mode is handled with kinetic prefactor
`abs(F2)`. For `--well-type single`, `Q=0` is the minimum. If `--well-type auto`
is used, a negative `F2` is treated as a double-well/TS reference and a positive
`F2` as a single-minimum reference.

The default derivative model follows the Burcl-Carter-Handy normal-coordinate
transformations for non-double-well references:

- `handy-gaussian` when `abs(F3)` is below
  `--anharmonic-cubic-threshold`. This uses
  `z = sqrt(1 - exp(-beta Q^2))` and rewrites the local F2/F4 terms in powers
  of `z^2` and `z^4`, preserving the derivatives through fourth order.
- `handy-morse` for a single minimum with non-zero `F3`. This uses
  `y = 1 - exp(-a Q)`, with `a = -F3/(3 F2)`, and rewrites the local F2/F3/F4
  expansion in powers of `y` through `y^4`, again preserving the derivatives.
- `inverse-power` for a double-well/TS case with non-zero `F3`; this fits
  `1/r + 1/r^2 + 1/r^3` analytically when the real-pole condition
  `F3^2 - F2*F4 >= 0` is satisfied. Otherwise `auto` falls back to the
  Handy-Morse form when possible, or Taylor for diagnostics.

All derivative models can be selected explicitly with
`--anharmonic-derivative-model`; the older `gaussian` and `morse` choices are
retained as simpler non-Handy forms. `--anharmonic-handy-beta` can be used to
override beta for the Handy Gauss-like coordinate. A high-order confining wall
is added by default using `--anharmonic-wall-degree 8` and
`--anharmonic-wall-height-cm 10000`. This wall starts at sixth or eighth order
and therefore does not change the derivatives through fourth order at `Q=0`.

For single-well references, the summary and
`*_anharmonic_vpt2_comparison.csv` compare the standard one-mode
normal-coordinate VPT2 levels with the variational levels computed on the
selected Handy-Morse or Handy-Gaussian analytic potential. The local F2/F3/F4
constants define the one-mode perturbative anharmonicity

```text
xe = 5 F3^2/(48 w) - F4/16
```

and the VPT2 levels are written as

```text
E_n = w(n+1/2) - xe(n+1/2)^2 + constant.
```

The variational column is obtained by diagonalizing the actual selected
analytic potential. This is the intended isomorphism check: the normal-mode
perturbative result is compared directly with a variational calculation on the
Handy Morse-like potential, and the same output can be generated for the
Gauss-like potential when F3 vanishes. The comparison is skipped for
double-well/TS references because there is no stable harmonic reference
oscillator at the barrier.

Outputs are written separately from Cartesian path analyses:
`*_anharmonic_summary.txt`, `*_anharmonic_levels.csv`,
`*_anharmonic_grid.csv`, optional `*_anharmonic_vpt2_comparison.csv`, and
`*_anharmonic_potential_levels.pdf/png`.

### Graphical Interface

The repository also includes a lightweight Tkinter launcher:

```bash
python3 scripts/mw_path_dvr_gui.py
```

It exposes the main groups of options for Gaussian anharmonic one-mode
calculations, 1D path calculations, and 2D grid Hamiltonians. The GUI always
shows the exact command line before running it; this keeps graphical runs
reproducible and makes it easy to transfer a successful setup back to scripts.

Example for a half-path inversion scan:

```bash
python3 scripts/mw_path_dvr.py \
  --gaussian-log path_scan.log \
  --log-selection all \
  --gaussian-energy post-scf \
  --boundary nonperiodic \
  --solver sinc-dvr \
  --compute-rotconst \
  --path-symmetry half-even \
  --well-type double \
  --potential-extension repulsive-exponential \
  --extension-length-au 50 \
  --extension-target-cm 10000 \
  --extension-points 32 \
  --potential-smoothing spline \
  --potential-spline-smoothing 0.0 \
  --prefix path_extended
```

### Gaussian Solver Details

The nonperiodic Gaussian solver follows the distributed-Gaussian formulation of
Hamilton and Light and the optimized-basis philosophy of Kelemen and Luber. The
primitive functions are normalized real Gaussians,

```text
g_i(s) = (2 alpha_i/pi)^(1/4) exp[-alpha_i (s-s_i)^2].
```

Overlap and kinetic-energy matrix elements are analytic. Potential and scalar
property matrix elements are evaluated with Gauss-Hermite quadrature over the
product Gaussian. The generalized eigenvalue problem is solved after
orthogonalizing the overlap matrix and discarding eigenvectors below
`--gaussian-overlap-threshold` or above the requested condition limit. This is
the same numerical safeguard needed for higher-dimensional product-basis
calculations.

By default the code performs a compact variational optimization of the Gaussian
basis before the final diagonalization. The current one-dimensional
implementation minimizes the sum of a selected set of low-lying eigenvalues
with respect to the internal Gaussian centers and widths. The end-point centers
are kept as anchors for the finite coordinate interval. If the nonlinear
optimization does not give a stable variational improvement, the code falls back
to the initial distributed Gaussian basis and reports this in `*_summary.txt`.

## 8. Two-Dimensional Product-Basis Calculations

The program can also solve a two-dimensional Hamiltonian from a rectangular CSV
grid. This mode is intended for cases such as two torsions, a torsion plus a
nonperiodic coordinate, or two nonperiodic coordinates. The CSV must contain one
row per grid point, with two coordinate columns and one energy column:

```text
q1,q2,energy_cm-1
0.000000,-1.000000,120.0
0.000000,-0.750000,80.0
...
```

The coordinates are assumed to be in the units used by the kinetic operator. If
the CSV contains angles or another coordinate convention, use `--q1-scale-au`
and `--q2-scale-au` to convert the tabulated values before solving.

Example with one periodic coordinate and one nonperiodic coordinate:

```bash
python3 scripts/mw_path_dvr.py \
  --grid2d-csv examples/gaussian_outputs/phi_q_grid.csv \
  --q1-key phi --q2-key q \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --boundary1 periodic --boundary2 nonperiodic \
  --solver1 fourier --solver2 gaussian \
  --basis1 16 --basis2 12 \
  --levels 20 \
  --prefix phi_q_2d
```

`--solver1 auto` and `--solver2 auto` select Fourier for periodic axes and the
optimized Gaussian solver for nonperiodic axes. `--solver1 sinc-dvr` or
`--solver2 sinc-dvr` can be used as a nonperiodic fallback.

The 2D calculation uses pruning. First, two one-dimensional reference problems
are built from the 2D potential. The default reference potential is the minimum
over the other coordinate:

```text
V1(q1) = min_q2 V(q1,q2)
V2(q2) = min_q1 V(q1,q2)
```

Alternative choices are `--prune-potential mean` and `--prune-potential cut`.
The selected one-dimensional Hamiltonians are diagonalized, the first
`--basis1` and `--basis2` eigenfunctions are retained, and the full 2D
Hamiltonian is then diagonalized in the reduced product basis. The final
potential matrix always uses the complete tabulated `V(q1,q2)`, not the sum of
the pruning potentials.

Additional pruning options are available. `--prune-energy-window1` and
`--prune-energy-window2` keep only one-dimensional functions within the selected
energy window above the corresponding one-dimensional ground state. The
requested `--basis1` and `--basis2` values remain upper limits, and
`--prune-min-basis1` and `--prune-min-basis2` enforce lower limits. After an
initial two-dimensional diagonalization, `--prune-product-coeff-threshold` can
discard product functions whose summed coefficient weight in the first
`--prune-product-states` states is below the threshold. This gives a compact
second diagonalization while retaining at least enough product functions for the
requested roots.

For convergence checks, use `--convergence-scan`. The scan can vary retained
basis sizes, reference pruning potentials, and downsampled grid strides:

```bash
python3 scripts/mw_path_dvr.py \
  --grid2d-csv examples/gaussian_outputs/phi_q_grid.csv \
  --q1-key phi --q2-key q \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --boundary1 periodic --boundary2 nonperiodic \
  --solver1 fourier --solver2 gaussian \
  --convergence-scan \
  --convergence-basis1 8,12,16 \
  --convergence-basis2 8,12 \
  --convergence-prune-potentials min,mean \
  --convergence-strides 2,1 \
  --levels 12 \
  --prefix phi_q_2d
```

The output `*_2d_convergence.csv` reports the first levels for each run and
their differences from the last successful reference row in the table. When a
single pruning strategy is used and the strides/bases are ordered from coarse to
fine, this last row is the most refined calculation.

### Coordinate-Dependent Kinetic Metric

The 2D Hamiltonian can use either a constant metric or a coordinate-dependent
metric. The coordinate-dependent form follows the Podolsky/Laane treatment of
large-amplitude vibrations, neglecting the small pseudopotential term. For two
coordinates `q1` and `q2`, the kinetic operator is written in Hermitian form as

```text
2 T = p1 g11(q1,q2) p1
    + p2 g22(q1,q2) p2
    + p1 g12(q1,q2) p2
    + p2 g12(q1,q2) p1
```

with `pi = -i d/dqi` in atomic units. This is the numerical analogue of the
`g44`, `g45`, and `g55` kinetic-energy functions used in Laane-type
two-dimensional large-amplitude Hamiltonians. If all `gij` are constants, this
reduces to

```text
T = -1/2 [ g11 d2/dq1^2 + 2 g12 d2/(dq1 dq2) + g22 d2/dq2^2 ].
```

The default is still the constant metric,

```text
--metric-mode constant --g11 1 --g22 1 --g12 0
```

This corresponds to two mass-weighted coordinates with no kinetic cross term.
A constant cross term can be supplied with `--g12`.

If the kinetic elements have already been generated elsewhere, provide them as
CSV columns:

```bash
python3 scripts/mw_path_dvr.py \
  --grid2d-csv examples/gaussian_outputs/phi_q_grid_with_metric.csv \
  --q1-key phi --q2-key q \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --boundary1 periodic --boundary2 nonperiodic \
  --metric-mode csv \
  --g11-key g11 --g12-key g12 --g22-key g22 \
  --basis1 16 --basis2 12 \
  --levels 20 \
  --prefix phi_q_2d_variable_metric
```

The recommended production route is automatic metric construction from the
optimized Cartesian structures on the same rectangular grid:

```bash
python3 scripts/mw_path_dvr.py \
  --grid2d-csv examples/gaussian_outputs/phi_q_grid.csv \
  --grid2d-geom-log examples/gaussian_outputs/phi_q_grid.log \
  --grid2d-geom-log-selection last-per-link \
  --q1-key phi --q2-key q \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --boundary1 periodic --boundary2 nonperiodic \
  --metric-mode geometry \
  --metric-stencil 5 \
  --metric-smoothing auto \
  --basis1 16 --basis2 12 \
  --levels 20 \
  --prefix phi_q_2d_geometry_metric
```

The CSV defines the ordering and coordinate values. The geometry source must
contain exactly one optimized structure for each CSV row, in the same row order.
It may be either a Gaussian log (`--grid2d-geom-log`) or a multi-XYZ file
(`--grid2d-geom-xyz`). For each grid point the program forms finite-difference
Cartesian derivatives with respect to the two selected coordinates. The default
finite-difference stencil uses five points when available; it is reduced
automatically on smaller grids or near nonperiodic boundaries. Neighboring
structures are translated to their centers of mass and mass-weighted
Eckart-aligned to the central structure before the derivative is formed. The raw
covariant metric is then

```text
B_ij(q1,q2) = sum_a dX_a/dqi . dX_a/dqj
```

where `X_a = sqrt(m_a) r_a` is expressed in mass-weighted atomic units.

By default, coordinate-dependent metrics are postprocessed before the
Hamiltonian is assembled:

- `--metric-smoothing auto` fits a smooth tensor-product surface to the metric
  elements. Periodic axes use Fourier terms and nonperiodic axes use Chebyshev
  polynomial terms. The maximum total order is set by
  `--metric-smooth-order` (default 3), and is reduced automatically if the grid
  is too small.
- For geometry-derived metrics the fit is applied to the covariant elements
  `cov11`, `cov12`, and `cov22`, after which the matrix is inverted. This keeps
  the smoothing closer to the directly differentiated Cartesian data.
- `--metric-smoothing none` skips the fit and uses the raw finite-difference
  metric.
- `--metric-eigenvalue-floor` sets a relative floor used to project any nearly
  singular pointwise metric back to positive-definite form. The summary reports
  how many points were projected and the maximum relative shift.

The contravariant kinetic metric used in the Hamiltonian is

```text
g(q1,q2) = B(q1,q2)^-1.
```

Thus the reduced masses and kinetic coupling are not fitted parameters: they
come directly from the molecular structures and the chosen coordinate grid. The
output `*_2d_grid.csv` includes `g11`, `g12`, `g22`, and, for geometry-derived
metrics, also the covariant elements `cov11`, `cov12`, and `cov22`. The summary
file records the smoothing mode, fitted order, relative root-mean-square change
of each fitted element, metric condition number, and positive-definite
projection diagnostics.

For five-membered rings, a 2D model can use the two endocyclic-dihedral
components already used to define the puckering phase,

```text
QC = sum_i c_i tau_i
QS = sum_i s_i tau_i
```

or any equivalent pair such as radial/phase coordinates. For six-membered
rings, the same machinery should be used with two explicitly chosen puckering
coordinates, for example two Cremer-Pople-like Cartesian components or a
pair such as `(q2, q3)` on a defined phase cut. The code does not assume a
unique six-membered-ring coordinate system: the selected `q1,q2` grid and the
Cartesian derivatives define the metric.

2D output files use the suffixes `*_2d_grid.csv`, `*_2d_levels.csv`,
`*_2d_expectations.csv`, `*_2d_summary.txt`, and `*_2d_potential.pdf`.

## 9. TODO and Planned Extensions

The 2D code now supports coordinate-dependent metrics from either CSV columns
or molecular structures, including smoothing and positive-definite projection.
The remaining methodological TODO is to decide whether the small Podolsky
pseudopotential term is negligible for the target application or should be added
as an optional correction.

Six- and seven-membered rings are a main motivation for the product-basis
design. Their puckering cannot in general be represented by the single
five-membered-ring phase used for tetrahydrofuran-like systems. Six-membered
rings naturally require at least a two-dimensional puckering description, and
seven-membered rings motivate a three-dimensional extension. The current 2D CSV
workflow is therefore intended as the intermediate layer between the present
four-/five-membered-ring scans and future higher-dimensional ring-puckering
Hamiltonians.

Other planned items:

- generalize the 2D product-basis machinery to an arbitrary number of
  coordinates;
- generate Gaussian inputs for 2D and 3D grids directly from generalized
  internal coordinates;
- add native puckering-coordinate builders for six- and seven-membered rings;
- add regression tests based on analytic one- and two-dimensional potentials.

## 10. Output Files

With `--prefix thf_phi_mass_weighted`, the main output files are:

```text
examples/gaussian_outputs/thf_phi_mass_weighted_profile.csv
examples/gaussian_outputs/thf_phi_mass_weighted_oriented.xyz
examples/gaussian_outputs/thf_phi_mass_weighted_grid.csv
examples/gaussian_outputs/thf_phi_mass_weighted_levels.csv
examples/gaussian_outputs/thf_phi_mass_weighted_expectations.csv
examples/gaussian_outputs/thf_phi_mass_weighted_thermal_expectations.csv
examples/gaussian_outputs/thf_phi_mass_weighted_summary.txt
examples/figs/thf_phi_mass_weighted_potential_levels.pdf
examples/figs/thf_phi_mass_weighted_potential_profile.pdf
examples/figs/thf_phi_mass_weighted_circular_potential.pdf
```

Contents:

- `*_profile.csv`: original path points, mass-weighted coordinate, energies,
  optional Cremer-Pople labels, and pointwise scalar properties.
- `*_oriented.xyz`: geometries reoriented along the path; with
  `--label-cremer-pople`, each frame comment also contains `CP_q2` and
  `CP_phi2`.
- `*_grid.csv`: discrete-variable-representation grid, interpolated potential,
  interpolated properties, and wavefunctions.
- `*_levels.csv`: eigenvalues.
- `*_expectations.csv`: expectation values of scalar properties, for example
  rotational constants.
- `*_thermal_expectations.csv`: 0 K and requested finite-temperature Boltzmann
  averages of scalar properties over the computed one-dimensional states.
- `*_summary.txt`: readable calculation summary.
- `*_potential_levels.pdf`: potential plot with levels and wavefunctions.
- `*_potential_profile.pdf`: potential profile versus phase, useful for barrier
  discussions.
- `*_circular_potential.pdf`: circular pseudorotation plot.

## 11. Method References

The implementation is intentionally lightweight, but the numerical choices are
based on the following references.

- J. C. Light, I. P. Hamilton, and J. V. Lill, "Generalized discrete
  variable approximation in quantum mechanics," Journal of Chemical Physics
  82, 1400-1409 (1985).
- I. P. Hamilton and J. C. Light, "On distributed Gaussian bases for simple
  model multidimensional vibrational problems," Journal of Chemical Physics
  84, 306-317 (1986).
- B. Poirier and J. C. Light, "Efficient distributed Gaussian basis for
  rovibrational spectroscopy calculations," Journal of Chemical Physics 113,
  211-217 (2000).
- E. Matyus, G. Czako, B. T. Sutcliffe, and A. G. Csaszar, "Vibrational
  energy levels with arbitrary potentials using the Eckart-Watson Hamiltonians
  and the discrete variable representation," Journal of Chemical Physics 127,
  084102 (2007).
- R. Burcl, S. Carter, and N. C. Handy, "The vibrational spectroscopy of
  polyatomic molecules in normal coordinates: II. Parameterisation of the force
  field," Chemical Physics Letters 373, 357-365 (2003).
- A. K. Kelemen and S. Luber, "Compact vibrational wave functions via linear
  optimization," International Journal of Quantum Chemistry 126, e70177
  (2026), DOI: 10.1002/qua.70177.
- M. A. Harthcock and J. Laane, "Calculation of two-dimensional vibrational
  potential energy surfaces utilizing prediagonalized basis sets and Van Vleck
  perturbation methods," Journal of Physical Chemistry 89, 4231-4240 (1985).
- N. Meinander and J. Laane, "Computation of the energy levels of
  large-amplitude low-frequency vibrations. Comparison of the prediagonalized
  harmonic basis and the prediagonalized distributed Gaussian basis," Journal
  of Molecular Structure 569, 1-24 (2001).

## 12. Quick Checks

To check only that a Gaussian log is read and aligned correctly:

```bash
python3 scripts/mw_path_dvr.py \
  --gaussian-log examples/gaussian_outputs/thf_phi_scan.log \
  --log-selection last-per-link \
  --check-only
```

The program prints:

```text
points=...
path_length_sqrtamu_angstrom=...
path_length_au=...
max_angular_residual=...
```

`max_angular_residual` should be small. It measures how well the local Eckart
reorientation has removed the rotational component of the path before the
mass-weighted distance was accumulated.

## 13. Common Problems

If the number of path points is wrong, check that:

- Gaussian finished all links;
- `--log-selection last-per-link` was used for multi-link inputs;
- `--log-selection all` was used for a single Gaussian relaxed-scan log with
  `Step number ... on scan point ...` records;
- the analyzed log is the complete scan log.

If the potential barrier is inconsistent with a double-hybrid or MP2 scan,
check the energy source. Use `--gaussian-energy post-scf` to force post-SCF
total energies; otherwise `auto` should use them when present, and `scf` will
explicitly reproduce the `SCF Done` curve.

If expected properties are missing from `*_expectations.csv`, check that:

- the properties are present in the Gaussian log;
- or `--compute-rotconst` was used to compute rotational constants from the
  geometries.

If the ring is not numbered `1,2,3,4,5`, change `--ring`, for example:

```bash
--ring 5,1,2,3,4
```

The atoms must still be supplied in cyclic order.
