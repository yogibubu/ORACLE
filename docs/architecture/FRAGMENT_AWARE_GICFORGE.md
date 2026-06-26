# Fragment-Aware GICForge

Date: 2026-06-26

## Contract

`oracle-fragments` owns fragment discovery and writes concrete `#FRAGMENTS`
records. GICForge consumes those records; it does not rediscover fragments.
The same package owns `#INTERACTION_CENTERS`, where bond midpoints, ring
centroids and future eta-n/H-bond centers are stored as virtual centers. These
centers are metadata over real atoms, not dummy atoms in the molecular graph.

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
- `FRAG_ORIENTATION` / `FROT`: TRIC/geomeTRIC-style exponential-map rotation
  component of a fragment local frame relative to the reference fragment frame.

When `#INTERACTION_CENTERS` has `STATUS BUILT`, GICForge also adds:

- `CENTER_ATOM_DISTANCE` / `CENTER_ATOM_DIST`: distance between a real atom and
  a virtual center defined by a bond midpoint, ring centroid, eta-n ligand
  centroid or future H-bond/contact center.

This covers atom-to-ring-center and atom-to-bond-center interactions used for
eta coordination, pi/H-bond contacts and transition-metal skeletal coordinates.
Additional tilt, torsion and frame-orientation primitives should reuse the same
center records rather than introducing separate center parsers.

These primitives are frozen in `#GIC` and available to ORACLE downstream tools.
They are serialized with `CLASS=SPECIAL_PROTECTED` and are reduced before
ordinary primitives. This prevents a stretch, bend or torsion from eliminating
an inter-fragment or atom-center coordinate that was explicitly requested by
fragment/topology state. Redundancy among special coordinates is still tested by
analytic B-matrix rank.

`[GAUSSIAN_GIC]` exports them as Gaussian symbolic GIC expressions using
`Fragment(...)`, `XCntr/YCntr/ZCntr(...)`, Cartesian `X/Y/Z`, and the
fragment-frame `P/Q/S` plus quaternion/exponential-map construction documented
in Gaussian's GIC guide. The Gaussian text export regularizes `||Kvec||` by a
tiny constant so exactly aligned fragment frames use the same small-rotation
limit instead of producing a `0/0` expression.

## B Matrix

The Wilson B matrix is evaluated from the frozen `#GIC` section, not by
regenerating the GIC basis. For fragment-aware coordinates, ORACLE differentiates
the same scalar definitions used for the Gaussian symbolic export:

- center-center and center-atom distances;
- atom to virtual bond/ring/contact center distances;
- Cartesian center translations;
- exponential-map components of relative fragment orientation.

The backend is analytic and intentionally centralized in `oracle-gicforge`.
Downstream GF, MORPHEUS and refinement tools must call this service instead of
carrying their own fragment-coordinate derivatives.

## Legacy Origin

The translation and rotation definitions follow the Merlino TRIC roadmap and
geomeTRIC reference model: fragment centers are geometric centroids, frames are
built from fragment anchor atoms, and relative orientations are expressed as the
three-component exponential map derived from the relative-orientation
quaternion. The native Python and Fortran B-matrix paths use the same analytic
definition, including the regular small-rotation limit `2*Kvec`.
