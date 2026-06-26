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
  selection used by `oracle_gicforge.definition`.

The `FROT` B row is analytic: derivatives are propagated through the centroid,
local frame, relative rotation matrix, quaternion and exponential-map
small-rotation limit. It is intended to be called from the imported legacy
`MkBNew` path instead of adding another local derivative implementation.
