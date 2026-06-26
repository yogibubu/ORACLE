# ORACLE Package Architecture

ORACLE is organized as narrow packages with explicit data contracts. The GUI
orchestrates these packages; it must not duplicate scientific logic.

## Package Boundaries

- `oracle_core`: CLI, configuration, manifests, workspace handling and common
  validation/errors.
- `oracle_gicforge`: GIC definition and B-matrix evaluation. It owns topology-driven
  GIC construction, optional symmetry adaptation, frozen GIC schemas and
  Gaussian-readable GIC blocks.
- `oracle_gf`: harmonic Cartesian-Hessian to internal-coordinate GF/PED. It
  owns `HessianInput`, Wilson GF linear algebra, frozen-GIC Hessian
  transformation, Pulay scaling and GF/PED report/CSV services.
- `oracle_vpt2_vci`: anharmonic normal-mode VPT2/VCI. It owns QFF data,
  normal-mode basis selection, VPT2, VCI and Davidson diagonalization. It does
  not build or evaluate GICs.
- `oracle_morpheus`: semiexperimental geometry refinement. Its default model
  calls `oracle_gicforge` for frozen GIC definitions and B matrices, then performs
  least-squares refinement and uncertainty propagation. Its optional
  symmetry-Cartesian model builds a Hessian-free translation/rotation-free
  Cartesian displacement basis, filters totally symmetric directions and
  propagates final errors to primitive internal coordinates without using a GIC
  B matrix.
- `oracle_dvr`: scan/grid to DVR levels and wavefunctions.
- `oracle_gaussian`: Gaussian input/output adapters. Gaussian is a file-format
  source, not a solver dependency.
- `oracle_engines`: active executable/source discovery and build checks for
  Fortran kernels.
- `oracle_gui` and `advanced`: PySide6 windows/controllers only. They call the
  package services and CLI workflows.

## GIC Data Flow

```text
Cartesian reference geometry
  -> oracle_gicforge.define_gics_from_cartesian(symmetrize=True|False)
  -> oracle.gic.definition.v1
       primitives, U matrix, labels, names, irreps, point group, symmetrized flag

Current Cartesian geometry
  -> oracle_gicforge.evaluate_gic_definition(schema, geometry)
  -> values, primitive B, GIC B, labels, names, irreps, point group
```

Symmetrization is an option of the definition stage only. Downstream programs
must consume the frozen labels and irreps; they must not re-symmetrize.

## GF Flow

```text
Hessian adapter
  -> oracle_gf.HessianInput
Frozen GIC definition + current geometry
  -> oracle_gicforge.evaluate_gic_definition
HessianInput + B(GIC)
  -> oracle_gf internal Hessian/G matrix
  -> optional Pulay scaling
  -> Wilson GF, frequencies, normal modes, PED
```

`oracle_gf` is physically separate from `oracle_vpt2_vci` so additional
harmonic analyses can be added without touching VPT2/VCI.

## VPT2/VCI Flow

```text
QFF adapter or normalized QFF text
  -> oracle_vpt2_vci.QuarticForceField
  -> active-mode selection/pruning/symmetry blocks
  -> VPT2 and/or VCI
  -> Davidson for large VCI spaces
```

VPT2/VCI works in Cartesian normal modes. It is intentionally independent from
GIC definition and B-matrix evaluation.

## SEfit Flow

Default GIC model:

```text
Cartesian parent geometry
  -> oracle_gicforge.define_gics_from_cartesian(symmetrize=True)
  -> GICForge post-pruning GICSYM schema/manifest
  -> totally symmetric frozen GIC subspace
  -> analytic B / Cartesian projector as needed
  -> semiexperimental least-squares refinement
  -> Cartesian rSE_e + primitive internal coordinates/errors
```

Optional symmetry-Cartesian model:

```text
Cartesian parent geometry
  -> GICForge SYCART symmetrized Cartesian coordinates
  -> translations/rotations projected out
  -> Cartesian displacement basis with assigned irreps
  -> totally symmetric symmetry-Cartesian subspace
  -> semiexperimental least-squares refinement in Cartesian amplitudes
  -> Cartesian rSE_e + primitive internal coordinates/errors
```

The two SEfit coordinate models share observations, correction conventions,
parameter classes, primitive constraints, covariance analysis and reporting.
Only the working-coordinate basis and Cartesian update map differ.
Neither model regenerates coordinates inside SEfit: both call the public
GICForge API and consume the frozen post-pruning schema or SYCART coordinates.

Reusable SEfit submodules:

- `oracle_morpheus.constraints`: public Gaussian-style constraint parsing,
  expression values/targets and analytic-vs-finite-difference B-matrix checks.
- `oracle_morpheus.diagnostics`: public SVD, uncertainty and iteration-trace
  CSV helpers.
- `oracle_morpheus.performance`: cached isotope-aware mass vectors used by
  observable and Jacobian builders.
