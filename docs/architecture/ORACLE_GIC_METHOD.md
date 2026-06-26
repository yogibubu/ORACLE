# ORACLE GIC Method

Date: 2026-06-26

## Scope

This note documents the ORACLE-native generalized internal coordinate method
implemented by GICForge. It covers coordinate construction, non-redundant
reduction, B-matrix ownership and symmetry policy. The same frozen definition is
used by Gaussian export, GF, SEfit, thermochemistry and anharmonic modules.

GICForge is not allowed to rediscover molecular state. It consumes the enriched
`xyzin` container produced upstream:

- geometry from ORACLE-Babel;
- validation status;
- topology, bonds, rings and synthons;
- point-group assignment;
- fragment records;
- interaction-center records for bond midpoints, ring centers and future
  eta-n/contact centers.

The output is a frozen `#GIC` section. Downstream tools must consume that
section and must not regenerate topology, fragments, centers or GIC rows.

The Python source of truth for the contract is `oracle_gicforge.policy`. It
defines backend names, rank method, reduction policy, primitive families,
reduction classes and symmetry source blocks. Strict Fortran77 kernels must
mirror this module rather than introduce independent family names or pruning
rules.

## Primitive Families

Ordinary primitives describe local molecular deformations:

- `STRETCH` / `R`;
- `BEND` / `A`;
- `LINEAR_BEND` / `L`;
- `TORSION` / `D`;
- `OUT_OF_PLANE` / `U`.

Special primitives describe objects that are not ordinary valence graph edges:

- `FRAG_DISTANCE` / `FC_DIST`, the distance between two fragment centers;
- `FRAG_CENTER_ATOM_DISTANCE` / `FCA_DIST`, the distance between a fragment
  center and an atom;
- `FRAG_TRANSLATION` / `FTRANS`, one Cartesian component of a fragment-center
  displacement;
- `FRAG_ORIENTATION` / `FROT`, one exponential-map component of a relative
  fragment orientation;
- `CENTER_ATOM_DISTANCE` / `CENTER_ATOM_DIST`, the distance between an atom and
  a virtual center such as a ring centroid, bond midpoint or eta/contact center.

Every primitive stored in `#GIC` carries a `CLASS` field. Ordinary primitives
are written as `CLASS=ORDINARY`. Special primitives are written as
`CLASS=SPECIAL_PROTECTED`.

## B Matrix

The Wilson B matrix is evaluated analytically from the frozen `#GIC` section.
Distance, angle and center-distance rows use closed-form derivatives. Fragment
translations use direct centroid derivatives. Fragment orientations use the
relative frame quaternion and its exponential-map components, with the same
small-rotation limit used by the Gaussian symbolic export.

This makes the B matrix a shared service. Python and Fortran backends must use
the same mathematical definitions, and downstream programs must call the
GICForge B-matrix API rather than carrying private derivative code.

## Non-Redundant Reduction

The reduction target is the vibrational rank:

- `3N-6` for nonlinear molecules;
- `3N-5` for linear molecules.

All rank decisions use normalized analytic B rows, a modified Gram-Schmidt
incremental rank basis and a single numerical tolerance. The important
ORACLE-specific rule is the reduction order:

1. Special protected primitives are tested first.
2. Ordinary primitives are used only to complete the rank after the protected
   set has been handled.
3. A protected primitive can be skipped only if its B row is singular or
   linearly dependent on previously accepted protected rows.
4. A protected primitive is never eliminated because an ordinary stretch, bend
   or torsion was accepted earlier.
5. If independent protected primitives alone exceed the vibrational rank,
   GICForge stops with a contract error instead of silently discarding one.

This policy is different from a global SVD or a pure candidate-order greedy
selection. It preserves chemically meaningful inter-fragment, atom-center and
fragment-orientation coordinates even when an ordinary primitive spans a similar
first-order displacement. Python uses the same MGS selector as the Fortran
`ORCGSEL` kernel in `engines/fortran/gicforge/frag_tric_bmat.f`.

The `#GIC` header records the policy as:

```text
RANK_METHOD analytic_b_matrix_mgs_greedy
REDUCTION_POLICY SPECIAL_PROTECTED_FIRST_THEN_ORDINARY_ANALYTIC_RANK
```

The frozen `#GIC` section also stores `[REDUCTION_DIAGNOSTICS]` with selected
primitive ids and primitives skipped because their analytic rows were singular
or linearly dependent. Human reports are generated from that frozen diagnostic
state, not by recomputing the reduction.

## Symmetrization

Symmetrization is applied after a valid non-redundant basis exists. It is not a
repair step for a redundant coordinate set. The symmetry layer may combine only
homogeneous source blocks: stretches with stretches, bends with bends,
fragment-orientation coordinates with fragment-orientation coordinates, and so
on. It must not mix unrelated coordinate families just because they share an
irreducible representation.

The first imported ORACLE implementation mirrors Merlino's type-local
GICForge symmetrizer. When `--symmetrize` is requested, GICForge groups
one-term frozen GICs by:

- symmetry source block;
- primitive family;
- atom-type signature.

For each eligible group, the first output coordinate is the normalized totally
symmetric sum of the source primitives. The remaining output coordinates are
normalized adjacent differences. The transformation is written directly in
`[FROZEN_GICS]` as `COEFFS=primitive_id:coefficient` and the B matrix is built
analytically as the same linear combination of primitive B rows. The final GIC
count therefore remains equal to the vibrational rank.

Symmetrized GIC names start with their assigned irreducible representation:
for example `A1StrS001`, `B2BendD001`, `AgFCDiS001` or `AStrS001` in `C1`.
The exact irrep label is stored separately in `IRREP=...`, while the coordinate
name uses a filesystem- and Gaussian-safe prefix.

For special coordinates, the intended symmetry source blocks are the protected
classes themselves:

- fragment-center distances;
- fragment-center to atom distances;
- fragment translations;
- fragment orientations;
- atom to virtual-center distances.

The symmetrized block must preserve the family counts and the protected
semantic role. If symmetry projection would destroy these counts or produce
ambiguous labels, the correct behavior is a clean stop.

The frozen `#GIC` section records `[SYMMETRY_DIAGNOSTICS]` with the method,
policy, status, transformed source groups and output labels. Human reports read
this stored diagnostic state; they do not recompute symmetry.

For optimization and least-squares refinement, the built `#GIC` header also
stores:

```text
SYMMETRY_GROUP ...
TOTAL_SYMMETRIC_IRREP ...
TOTAL_SYMMETRIC_GIC_COUNT ...
TOTAL_SYMMETRIC_GICS ...
```

Downstream modules must use this total-symmetric subset for
symmetry-preserving equilibrium refinement. In `C1`, all GICs transform as
`A`; in higher symmetry, only the GICs whose `IRREP` matches
`TOTAL_SYMMETRIC_IRREP` are active.

The helper `oracle_gicforge.symmetry.gic_symmetry_source_blocks` exposes these
homogeneous blocks to future GICSYM implementations and reports. It is a
preparatory layer: non-C1 projection remains constrained by backend
availability, but the grouping contract and local sum/difference
symmetrization are now explicit.

## Python And Fortran Backends

ORACLE keeps two Fortran layers:

- `engines/fortran/gicforge/frag_tric_bmat.f`, the ORACLE fragment/TRIC kernel
  that mirrors Python definitions for special coordinates and protected-first
  reduction;
- `engines/fortran/gicforge/legacy_merlino`, the imported Merlino3.0 strict
  Fortran77 GICForge source tree.

The legacy tree is compiled with `engines/fortran/gicforge/compile_legacy` and
is treated as a reference identity target. Python and Fortran paths are allowed
to coexist only when they share the frozen `#GIC` contract, family names,
protected-class policy, B-matrix definitions and regression fixtures.

The native Python API deliberately separates the stages that run at different
frequencies:

- `construct_gic_definition_from_xyzin`: topology-driven candidate generation,
  protected-first non-redundant reduction and frozen unsymmetrized GICs;
- `symmetrize_gic_definition`: post-reduction symmetry adaptation and active
  total-symmetric manifest;
- `build_gic_b_matrix`: analytic evaluation of the frozen definition for the
  current Cartesian geometry.

Optimizers and least-squares solvers call the first two stages once at setup
time, then call the B-matrix evaluator at every geometry iteration.

## Reports

`oracle gicforge report molecule.xyzin [report.txt]` emits a readable report
from the frozen `#GIC` state. It includes backend, point group, rank, candidate
count, selected family counts, protected counts, reduction diagnostics, symmetry
source blocks and final frozen GIC labels. This report is intended for method
debugging and scientific review before GF, SEfit or anharmonic steps consume the
coordinates.

## Gaussian Export

Gaussian output is a representation of the frozen ORACLE GICs, not the owner of
the method. Native Gaussian primitives are emitted directly where possible.
Fragment and virtual-center coordinates are emitted through Gaussian ReadGIC
symbolic expressions using `Fragment(...)`, `XCntr/YCntr/ZCntr(...)`,
coordinate references and regularized exponential-map expressions.

When a frozen GIC is a linear combination, Gaussian export emits the same
linear combination of the primitive symbolic expressions. It does not silently
fall back to the first primitive in the combination.

The Gaussian block is therefore reproducible from the frozen `#GIC` section and
the saved geometry. It is not used as the internal source of truth for rank,
symmetry or B-matrix construction.

## Guarantees

For a successful build, GICForge guarantees:

- the final GIC count equals the vibrational rank;
- all selected rows have finite analytic B rows;
- special protected coordinates are selected before ordinary coordinates;
- the frozen primitive list, Gaussian export and B-matrix rows refer to the same
  coordinate labels;
- downstream modules can restart from `xyzin` without rerunning topology,
  fragmentation or coordinate selection.

Any violation of these guarantees is a contract failure and should stop the
pipeline before GF, SEfit, thermochemistry or anharmonic analysis starts.
