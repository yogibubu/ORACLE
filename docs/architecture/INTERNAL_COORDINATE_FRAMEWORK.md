# MATRIX Internal Coordinate Framework

This document is the MATRIX adaptation of the internal-coordinate architecture
notes from the `newcoord` repository.  It is a refactoring contract, not a
request to redesign the mathematics.  Merlino/GICForge numerical behaviour is
preserved unless a change is explicitly justified, tested, and documented.

## Guiding Rules

- Topology is built before geometry-dependent analysis.
- Geometry is evaluated before coordinate values and B rows.
- Chemically meaningful coordinate definitions are generated before numerical
  reduction.
- Numerical linear algebra is confined to reduction, symmetrization and GF-like
  stages.
- Every stage appends information to the shared molecular/xyzin model and must
  not silently invalidate earlier stages.
- Deterministic output is mandatory whenever chemically possible.
- New coordinate types are integrated by adding a generator, registering it and
  adding regression tests.

## Pipeline

The canonical coordinate workflow is:

```text
Input
  -> topology perception
  -> topology analysis and validation
  -> geometry evaluation
  -> coordinate generation
  -> coordinate classification
  -> redundancy elimination
  -> symmetrization
  -> validation and serialization
```

In MATRIX this maps to:

| Stage | MATRIX owner | Notes |
| --- | --- | --- |
| Input/import | LINK and QM adapters | Reads XYZ, QM output or SMILES; writes shared xyzin sections. |
| Topology | LINK / matrix-chem | Owns bonds, rings, fragments, synthons and validation. |
| Geometry | LINK / matrix-chem | Owns coordinates, inertia, symmetry thresholds and auxiliary-node positions. |
| Coordinate generation | NEO/GICForge | Builds primitive and chemically adapted coordinates. |
| Classification | NEO policy layer | Assigns coordinate family, reduction class and symmetry block. |
| Redundancy elimination | NEO reducer | Uses analytic B rows and type-local rank logic. |
| Symmetrization | NEO symmetrizer | Projects only after a globally non-redundant basis exists. |
| Downstream use | GF, MORPHEUS, DVR, VPT2/VCI | Consume frozen GIC/SYCART definitions; they do not rebuild coordinates. |

## Data Ownership

Topology contains graph-only data:

- atoms and atom properties;
- bonds and bond orders;
- rings, fused rings, spiro and bridged systems;
- fragments and pseudo-bonds;
- auxiliary-node definitions;
- synthon/equivalence classes.

Geometry contains coordinate-dependent data:

- Cartesian coordinates;
- centers of mass and geometric centers;
- inertia tensor and principal axes;
- fragment orientations;
- instantiated auxiliary-node positions.

Coordinate definitions contain:

- coordinate type and family;
- participating atoms or auxiliary nodes;
- parent primitives for composite coordinates;
- generating algorithm;
- reduction class;
- symmetry block;
- provenance.

The frozen `#GIC` and `#SYCART` xyzin sections are the contract consumed by
later tools.  GF, MORPHEUS and anharmonic modules must not reparse Gaussian
ReadAllGIC blocks or regenerate topology to reinterpret those coordinates.

## Coordinate Generator Boundary

Each coordinate type has one logical generator.  The current production path is
still the Merlino-compatible NEO/GICForge implementation; the generator registry
is a migration boundary for gradual refactoring.

A generator may:

- read topology and geometry;
- create coordinate definitions;
- attach metadata and provenance.

A generator must not:

- mutate topology;
- perform redundancy elimination;
- perform symmetrization;
- diagonalize matrices or solve eigenproblems;
- manipulate Hessians.

Current extraction status:

| Generator | Status | Behavioural contract |
| --- | --- | --- |
| `StretchGenerator` | Extracted in `matrix_neo.generators.generate_stretch_coordinates` | Must remain identical to the stretch part of `definition._primitive_candidates`. |
| `LocalXHStretchGenerator` | Extracted in `matrix_neo.generators.generate_stretch_coordinates` | Uses the same opt-in X-H policy as NEO; GF treats only these rows as local X-H. |
| `LocalSymmetryAngleSALCGenerator` | Refactor boundary around the Merlino local-angle path | Must build SALCs for local angle environments through coordination 9 before it can replace the legacy implementation. |
| Other generators | Still embedded in the Merlino-compatible NEO path | Must be extracted one family at a time with golden tests before behaviour changes. |

## Reduction and Symmetry Boundary

Reduction receives generated coordinates and analytic B rows.  It may use
modified Gram-Schmidt, rank-revealing logic, BBt diagnostics or transformed
Hessian information when available.  Tie-breaking must prefer chemically
meaningful coordinates when numerically stable.

Symmetrization receives the already non-redundant coordinate set.  It must
preserve coordinate families and avoid mixing unrelated physical motions.  A
symmetry-adapted coordinate basis that cannot reproduce the expected
vibrational representation is a contract failure, not a downstream GF problem.

## Regression Corpus

Changes to coordinate construction, reduction, symmetry labels, B rows or
downstream GF use require tests on representative systems:

- acyclic molecules;
- isolated rings;
- fused and polycyclic rings;
- spiro and bridged systems;
- linear molecules;
- weak complexes and pseudo-bonds;
- metal/auxiliary-center cases when relevant;
- Python/Fortran comparison where the Fortran path exists.

Golden tests should compare at least coordinate counts, labels/families,
selected primitive/GIC rows, symmetry diagnostics and B-row numerical values
where available.
