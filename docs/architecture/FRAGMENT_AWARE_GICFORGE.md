# Fragment-Aware GICForge

Date: 2026-06-26

## Contract

`oracle-fragments` owns fragment discovery and writes concrete `#FRAGMENTS`
records. GICForge consumes those records; it does not rediscover fragments.

The first implemented strategy is `CONNECTED_COMPONENTS`, which materializes:

- fragment atom lists;
- geometric centers;
- deterministic local frames compatible with Gaussian rotor GICs;
- a reference fragment, chosen as the largest component.

## GICForge Primitive Families

When `#FRAGMENTS` has `STATUS BUILT`, GICForge adds fragment primitive
candidates before ordinary angle/torsion pruning:

- `FRAG_DISTANCE` / `FC_DIST`: distance between two fragment centers;
- `FRAG_CENTER_ATOM_DISTANCE` / `FCA_DIST`: distance between a fragment center
  and an anchor atom;
- `FRAG_TRANSLATION` / `FTRANS`: Cartesian component of a center-center
  displacement relative to the reference fragment;
- `FRAG_ORIENTATION` / `FROT`: quaternion-vector component of a fragment local
  frame relative to the reference fragment frame.

These primitives are frozen in `#GIC` and available to ORACLE downstream tools.
`[GAUSSIAN_GIC]` exports them as Gaussian symbolic GIC expressions using
`Fragment(...)`, `XCntr/YCntr/ZCntr(...)`, Cartesian `X/Y/Z`, and the
fragment-frame `P/Q/S` plus quaternion construction documented in Gaussian's
GIC guide.

## Legacy Origin

The translation and rotation definitions follow the Merlino/Gaussian rotor
model: fragment centers are geometric centroids, frames are built from fragment
anchor atoms, and relative orientations are expressed through quaternion
components.
