# ORACLE Enriched XYZ Container

The enriched XYZ file is the canonical communication file between ORACLE
modules.

It starts as a normal XYZ file and is enriched step by step by tools that append
or replace named sections. This is the pipeline spine: each tool can be run
independently, but it consumes the same container and preserves sections owned
by other tools.

## Non-Negotiable Rules

- Use shared ORACLE APIs to read, write and replace sections.
- Preserve unrelated sections exactly when updating a file.
- Do not duplicate parsers in GUI, scripts or individual scientific modules.
- Treat external file formats as adapters. Convert them into enriched XYZ
  sections before downstream consumption.
- Give every section a `SCHEMA` line.
- Keep section names uppercase.

## Minimal Flow

```text
plain XYZ
  -> oracle-chem adds #BASIC, #TOPOLOGY and #SYMMETRY
  -> oracle-fragments adds #FRAGMENTS, #FRAGMENT_LIBRARY or #ASSEMBLY
  -> oracle-chem adds #VALIDATION
  -> oracle-gicforge adds #GIC and optionally #SYCART
  -> branches:
       oracle-gaussian writes Gaussian inputs or imports QM outputs through oracle-qm
       oracle-qm adds #CARTESIAN_HESSIAN, #NORMAL_MODES or #QFF
       oracle-gf adds #GF_PED from #CARTESIAN_HESSIAN plus frozen #GIC
       oracle-rovib adds #ROTATIONAL, #VIBRATIONAL, #DELTABVIB, #CORIOLIS or #QCENT
       oracle-thermo adds #THERMO from BASIC/ROTATIONAL/VIBRATIONAL state
       oracle-morpheus adds #ISOTOPOLOGUES and #MORPHEUS
       oracle-trinity adds #TRINITY for external-gradient geometry optimization
       oracle-vpt2-vci adds #VPT2_VCI from normalized anharmonic data
       oracle-dvr adds #DVR
```

Molpro and MRCC outputs enter the same flow through `oracle-molpro` and
`oracle-mrcc`: their adapters normalize geometry, charge and multiplicity into
`MolecularGeometry`, then ORACLE-Babel writes the shared sections.

The GUI should display and orchestrate this state, not own a parallel data
model.

The same file is also the standalone input contract. A tool may be run directly
from an already prepared `xyzin` file when the sections it needs are present.

## Gaussian Topology Overrides

When ORACLE-Babel imports a Gaussian `.log` or `.out`, `oracle-gaussian` may
write `#GAUSSIAN_TOPOLOGY` before the topology/synthon pass. Only two Gaussian
quantities are accepted as topology overrides:

- CM5 atomic charges;
- Mayer bond orders.

If CM5 charges are present, `#SYNTHONS` uses them as the atomic charge column
and records `CHARGE_SOURCE Gaussian CM5`. If CM5 charges are absent, synthons
use the ORACLE electronegativity charge model and record
`CHARGE_SOURCE Synthons electronegativity model`.

If Mayer bond orders are present, topology and synthons use them and record
`BOND_ORDER_SOURCE Gaussian Mayer`. If Mayer bond orders are absent, topology
uses the ORACLE continuous Pauling bond-order model and records
`BOND_ORDER_SOURCE Topology Pauling continuous model`. Gaussian total bond
orders are not a fallback source.

## Gaussian Rovibrational Promotion

Gaussian log/out text is converted once by `oracle-gaussian` before downstream
tools run. `oracle gaussian promote-rovib` promotes harmonic frequencies, IR
intensities, anharmonic chi matrices, rotational constants and vibrational
rotational corrections into shared ORACLE sections:

- `#VIBRATIONAL` stores frequencies, intensities and optional chi values.
- `#ROTATIONAL` stores rotational constants, point group, Watson reduction,
  temperature, dipole components and the DeltaBvib bridge values consumed by
  rovibrational utilities and rotational spectroscopy.
- `#DELTABVIB` stores the normalized DeltaBvib values and, when available, the
  Gaussian alpha rows used to compute them.

GF, Thermo, SEfit/MORPHEUS and anharmonic workflows consume these sections
rather than reparsing Gaussian output.

`oracle rovib wmsrot-input molecule.xyzin` exports the normalized rotational
state to the browser WMS-Rot input format. `oracle rovib wmsrot-run
molecule.xyzin --out molecule.rotational.csv` calls the copied first-party
WMS-Rot Hamiltonian engine locally and writes the generated line list. Future
GUI state should record generated line-list/broadened-spectrum artifacts under
`#ROTATIONAL_SPECTRUM`.

## SEFit / MORPHEUS State

`oracle semiexp --xyzin molecule.xyzin` updates `#MORPHEUS` by default after a
successful semiexperimental fit. The section uses schema
`oracle.xyz.morpheus.v1` and records the run directory, manifest, report paths,
fitted geometry, parameter/residual/rotational-constant CSV files and summary
diagnostics such as backend, coordinate model, observable, components, RMS,
rank, condition number, iteration count and warning count.

Use `--no-write-section` for standalone benchmark or scratch runs that should
leave the active container unchanged. GUI controllers consume this section
rather than scanning the MORPHEUS output directory.

## TRINITY Geometry Optimization State

`#TRINITY` stores the request/state for geometry optimizations that call an
external program at each step to obtain energy and gradient. TRINITY stands for
Trust-Region Interface for Numerical Iterative Trajectories with external
energY/gradients. The first implementation is intentionally a skeleton:
`oracle trinity prepare molecule.xyzin --run-dir run --engine-command "..."`
writes schema `oracle.xyz.trinity.v1`, a prepared run manifest and the
persistent request.

The section records:

- external engine command and protocol;
- coordinate model (`gic` or `cartesian`) and active space, with
  `total_symmetric` as the default GIC active-space policy;
- energy/gradient units;
- max-step and trust-region tolerances;
- expected trajectory, final-geometry and energy/gradient log paths.

The future TRINITY runner must be autonomous: if `#TRINITY` is present, it reads
that section and the shared sections it names, then writes outputs back into the
same container. GUI controllers consume `#TRINITY` through `oracle-trinity` and
must not duplicate optimizer input parsing.

## QM Tensor Sections

`oracle-qm` owns the shared tensor sections used after external QM adapters have
finished parsing:

- `#CARTESIAN_HESSIAN` stores atomic numbers, Cartesian coordinates in bohr,
  atomic masses, the packed lower Cartesian Hessian and harmonic frequencies.
- `#NORMAL_MODES` stores normal-mode vectors as a mode by Cartesian-coordinate
  matrix, plus their frequencies.
- `#QFF` stores harmonic and anharmonic frequencies plus indexed cubic and
  quartic normal-coordinate force constants.

`oracle gaussian promote-fchk` writes these sections from Gaussian FCHK data.
GF/PED can then run from `oracle gf --xyzin molecule.xyzin` without reparsing the
FCHK. VPT2/VCI loaders can read `#QFF` directly from the same container.

`oracle-qm` also owns the normalized electronic sections:

- `#ELECTRONIC` stores electronic-state records with canonical columns
  `LABEL`, `ENERGY_HARTREE`, `ENERGY_EV`, `MULTIPLICITY`, `SYMMETRY` and
  `SOURCE`.
- `#TRANSITIONS` stores transition records with canonical columns `FROM`, `TO`,
  `ENERGY_EV`, `WAVELENGTH_NM`, `OSC`, `STRENGTH` and `SOURCE`.
- `#ORBITALS` stores external orbital/density/geometry file records with
  canonical columns `KIND`, `FORMAT`, `ROLE`, `PATH`, `LABEL` and `SOURCE`.
  Supported external file formats include FCHK/FCH, Molden, Cube/Cub and XYZ.

`oracle gaussian promote-fchk` registers the FCHK file in `#ORBITALS` and writes
the ground electronic state when the FCHK contains a total energy. `oracle
gaussian promote-electronic` promotes excited-state energies and oscillator
strengths from Gaussian logs and can register associated Molden/Cube/FCHK files.
GUI and scientific tools consume these sections; they must not parse Gaussian
logs or FCHK files privately.

`oracle gf --xyzin molecule.xyzin` updates `#GF_PED` by default. The section
uses schema `oracle.xyz.gf_ped.v1` and stores the Hessian source, coordinate
source, point group, matrix model, nonbonded correction label, frequencies,
GIC labels/irreps and the PED matrix. If a readable report or CSV directory is
requested, their paths are recorded in the same section. Use
`--no-write-section` only for standalone report/benchmark runs that should not
modify the project container.

`oracle vpt2-vci --xyzin molecule.xyzin` consumes `#QFF` as the canonical
standalone input. `--fchk` and `--qff-file` remain adapter/compatibility entry
points; scientific VPT2/VCI code should not reparse Gaussian output directly.

`oracle vpt2-vci --xyzin molecule.xyzin --run-dir runs/vpt2_vci` writes the
readable report, CSV comparison tables, an `oracle.run.v1` manifest and a
`#VPT2_VCI` section. `oracle vpt2-vci --collect molecule.xyzin` rereads those
outputs, records `OUTPUT_*` pointers and sets the section status to `complete`,
`partial` or `prepared`. GUI controllers consume this normalized state through
`oracle-vpt2-vci`, not by scanning CSV files directly.

## DVR State

`#DVR` stores the normalized DVR workflow request and output pointers:

- source Gaussian scan/path log handled by the adapter workflow;
- run directory, figure directory, prefix, boundary and selected solver;
- rotconst/Cremer-Pople/check-only switches;
- manifest path and expected grid, summary and level-table paths.

The raw Gaussian log remains an external adapter input. `oracle dvr prepare`
normalizes the request, writes an `oracle.run.v1` manifest and, when an `xyzin`
is supplied, updates `#DVR` so the GUI and later workflow steps can discover
the same run state without duplicating parser logic. `oracle dvr run` executes
the same normalized request directly: either from `LOG --outdir` or from an
existing `#DVR` section with `--xyzin`.

After the backend has run, `oracle dvr collect molecule.xyzin` reads the
produced summary/levels/grid/expectation files, detects optional 2D,
anharmonic and Fortran outputs, updates `#DVR` with `OUTPUT_*` pointers and
sets the section status to `complete`, `partial` or `prepared`. The command
`oracle dvr run --xyzin molecule.xyzin` performs this collection immediately
after a successful backend execution. GUI controllers must consume this
normalized state through `oracle-dvr` rather than scanning the run directory
themselves.

## GIC State

`#GIC` stores the frozen coordinate contract, not only a list of Gaussian input
lines. In built files it includes:

- selected primitive coordinates;
- final frozen GICs, including linear-combination coefficients;
- point-group/symmetry-group metadata;
- total-symmetric irrep and active total-symmetric GIC list;
- reduction diagnostics;
- symmetrization diagnostics;
- Gaussian ReadGIC text generated from the frozen state.

Downstream modules should consume the frozen GICs and diagnostics directly and
use `[GAUSSIAN_GIC]` only when writing Gaussian inputs. Optimizers and
least-squares refinements should use `TOTAL_SYMMETRIC_GICS` for
symmetry-preserving active variables and reevaluate the B matrix from the frozen
definition at each geometry step.

## Compatibility

The ORACLE `xyzin` format is the historical source:

- ordinary XYZ block first;
- uppercase appended sections;
- section replacement must preserve all unrelated sections.

ORACLE keeps this behavior and upgrades legacy section schemas from
`oracle.xyzin.*` to `oracle.xyz.*` as modules migrate.
