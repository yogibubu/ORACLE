# ORACLE GICForge Fortran Backend

This directory contains the strict Fortran77-side building blocks used while
porting the legacy Merlino GICForge backend into ORACLE.

`legacy_merlino/` vendors the complete Merlino3.0 GICForge Fortran source set.
It is built with:

```sh
./compile_legacy
```

The wrapper writes `build/gicforge_legacy` and keeps build logs under
`build/`. This backend is treated as a reference/identity target until the
ORACLE-native and legacy implementations agree on the frozen `#GIC` contract.

`frag_tric_bmat.f` is the first shared fragment-coordinate library. It mirrors
the native Python GICForge definitions for:

- `FC_DIST`: fragment-center / fragment-center distance;
- `FCA_DIST`: fragment-center / atom distance;
- `FTRANS`: Cartesian component of a fragment-center displacement;
- `FROT`: TRIC/geomeTRIC-style exponential-map rotation component between two
  fragment local frames.

It also mirrors the Python non-redundant reduction contract:

- `ORCGSPC` classifies the same special protected primitive families used by
  Python (`FRAG_DISTANCE`, `FRAG_CENTER_ATOM_DISTANCE`, `FRAG_TRANSLATION`,
  `FRAG_ORIENTATION`, `CENTER_ATOM_DISTANCE`);
- `ORCGSEL` performs the same protected-first modified Gram-Schmidt rank
  selection used by `oracle_gicforge.definition`;
- `ORCGLCB` builds a frozen GIC B row as the linear combination of primitive B
  rows stored in `COEFFS=`.

`legacy_merlino/gic_type_symmetry.f` is the strict Fortran77 reference for the
local GIC symmetrizer imported into Python. Both paths use the same policy:
group one-term coordinates only inside homogeneous type blocks, write the first
coordinate as a normalized symmetric sum, and write the remaining coordinates as
normalized adjacent differences. ORACLE serializes that transformation in
`#GIC/[FROZEN_GICS]` as primitive coefficients so Python and Fortran downstream
tools can build the same B rows. Symmetry labels are part of the frozen
contract: symmetrized coordinate names begin with their assigned irrep and the
header records the total-symmetric active subset used by optimizers and
least-squares fits.

`legacy_merlino/mkcyc.f`, `legacy_merlino/mksalc.f` and
`legacy_merlino/gicprune.f` are the strict Fortran77 reference for ring
coordinate construction, cyclic SALCs, fused-ring butterfly coordinates and
block-local pruning. ORACLE-native Python now exposes the corresponding family
names `CYCLIC_BEND`, `CYCLIC_TORSION`, `CONDENSED_RING_TORSION` and
`BUTTERFLY`; these are ordinary rank coordinates, but they must remain separate
homogeneous blocks and must not be collapsed into generic bend/torsion families.

`legacy_merlino/symm.f` is the reference source for Merlino-style point-group
families and operation builders. ORACLE keeps its classifier families aligned
with the Python character/projector layer for `Cn`, `Cnv`, `Cnh`, `Dn` and
`Dnh`; `Dnd` is covered through the same Merlino-style `D_n + sigma_d D_n`
operation construction. Python now has matrix-classified projector character
tests for `Td`, `Oh` and `Ih`, and the strict Fortran77 source now exposes
`OPS_I`/`OPS_IH` builders for the icosahedral groups. Python projector tests
now include selected `RPck` ring-puckering sources and verify that derived
`QPck`/`PhiP` Gaussian labels remain attached after symmetrization.

`oracle_engines.run_legacy_gicforge` is the executable parity harness for this
vendored backend. It writes normalized `provin` and lowercase `xyzin` files,
runs `build/gicforge_legacy`, parses `provout` and `bmat.out`, and lets tests
compare final GIC ranks, ring-puckering labels and Wilson-B row spaces against
ORACLE-native GICForge. The harness writes eight-decimal XYZ coordinates because
the legacy tokenizer is unstable with longer decimal tokens.

Construction, symmetrization and B evaluation are intentionally separate.
Fortran optimizers should read the frozen primitive/GIC coefficient table once,
then call `ORCGLCB` at each geometry iteration to combine the current primitive
B rows into the current frozen GIC B matrix.

The `FROT` B row is analytic: derivatives are propagated through the centroid,
local frame, relative rotation matrix, quaternion and exponential-map
small-rotation limit. It is intended to be called from the imported legacy
`MkBNew` path instead of adding another local derivative implementation.
