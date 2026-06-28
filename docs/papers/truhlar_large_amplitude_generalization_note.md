# MATRIX Technical Note: From Truhlar Torsional Projection To GIC/GF/DVR Large-Amplitude Coordinates

## Scope

This note records the MATRIX strategy for generalizing the torsional-mode
identification approach of Chen, Zhang, Truhlar, Zheng and Xu
(`bibliography/torsion_truhlar.pdf`) to a broader large-amplitude framework.
The goal is not to reproduce their projection machinery inside MATRIX.  MATRIX
already owns a frozen, non-redundant, chemically classified and optionally
symmetry-adapted GIC basis.  GF therefore works directly in that basis and
passes the correct kinetic metric to DVR or thermochemistry without a second
torsional-identification layer.

## What The Truhlar Work Solves

The JCTC 2022 MS-T(CD) paper addresses a practical bottleneck in
multistructural torsional anharmonicity: in complex molecules, torsions are
often mixed with other internal motions, and choosing a good non-redundant
torsional coordinate set by hand is fragile.  Their solution starts from
redundant internal coordinates, constructs delocalized constrained coordinates,
and separates torsional and non-torsional subspaces so that coupled torsional
frequencies, barriers and partition-function corrections can be evaluated more
robustly.

This is an important reference point for MATRIX, but our starting point is
different.  In MATRIX the coordinate problem is solved upstream by NEO: the
same topology, symmetry, primitive families, protected coordinates and analytic
B rows are frozen once and reused by GF, SEFit/MORPHEUS, DVR, VPT2/VCI and the
GUI.

## MATRIX Generalization

MATRIX generalizes the torsional-only problem in three directions.

First, the large-amplitude set is not limited to ordinary dihedrals.  It may
contain torsions, ring puckering coordinates, butterfly coordinates,
out-of-plane or improper coordinates, linear-bend coordinates and protected
fragment/center coordinates.  The selected set is denoted \(S\).

Second, MATRIX does not identify large-amplitude modes by projecting Cartesian
normal modes after the GF calculation.  The relevant coordinates are explicit
rows of the frozen GIC basis:

\[
  q = (q_1, \ldots, q_N), \qquad S \subset \{1,\ldots,N\}.
\]

GF transforms the Cartesian Hessian to this basis and obtains

\[
  F = J^T H_x J, \qquad G = B M^{-1} B^T,
\]

where \(B = \partial q / \partial x\), \(M\) is the Cartesian mass matrix and
\(J\) is the internal-to-Cartesian backtransform used by MATRIX GF.

Third, the kinetic information for DVR is obtained from the inverse metric of
the final non-redundant basis:

\[
  G^{-1} = (B M^{-1} B^T)^{-1}.
\]

For a one-dimensional torsion, \((G^{-1})_{ii}\) is the reduced moment in the
current GIC metric.  For a multidimensional large-amplitude DVR the object to
carry forward is not the list of diagonal terms but the complete selected
sub-block

\[
  G^{-1}_{SS}.
\]

The off-diagonal elements are the kinetic couplings between large-amplitude
coordinates and are required whenever the DVR is multidimensional.

## Hamiltonian Contract

For a selected large-amplitude coordinate vector \(Q_S\), the local harmonic
information delivered by GF is

\[
  F_{SS}, \qquad G_{SS}, \qquad G^{-1}_{SS}.
\]

The downstream DVR Hamiltonian should be built from the potential
\(V(Q_S)\) and the kinetic metric.  In the simplest constant-metric form near
the reference structure,

\[
  \hat{T}_S =
  -\frac{1}{2}\sum_{i,j \in S}
  (G^{-1}_{SS})_{ij}
  \frac{\partial^2}{\partial Q_i \partial Q_j}.
\]

This is already more general than a scalar hindered-rotor reduced moment.  It
supports coupled torsions, coupled puckering coordinates, torsion-puckering
blocks and fragment-orientation coordinates, provided that the coordinates are
kept in a homogeneous physical block when required by the coordinate model.

For a 1D periodic torsion with periodicity \(n\), equilibrium position \(q_0\)
and local curvature \(F_{ii}\), the minimal one-term potential estimate is

\[
  V(q) = A[1-\cos(n(q-q_0))], \qquad A = F_{ii}/n^2.
\]

This estimate is only a DVR/bootstrap model.  Production DVR should prefer
explicit scan data when available.  The same idea can be applied to non-torsion
large-amplitude angles only after the coordinate periodicity and physical range
are explicitly defined.

## Relation To The Truhlar Projection Strategy

The Truhlar MS-T(CD) strategy solves a subspace-identification problem: locate
torsional degrees of freedom in a redundant internal-coordinate description and
separate them from the rest.  MATRIX moves that decision earlier:

- NEO builds protected coordinate families before reduction.
- Reduction is performed within chemically meaningful families.
- Symmetry projection does not mix coordinate types.
- GF validates the resulting basis through \(F\), \(G\), symmetry blocks and
  PED diagnostics.
- DVR receives the actual \(F_{SS}\) and \(G^{-1}_{SS}\) blocks.

Thus MATRIX avoids a second torsional projection step.  The scientific
generalization is that the same machinery applies to any explicitly protected
large-amplitude coordinate family, not only internal rotations around single
bonds.

## Periodicity And Coordinate Ranges

For ordinary torsions, MATRIX estimates periodicity from topology and synthon
equivalence around the central bond.  Double or high-order bonds can be
excluded from hindered-rotor treatment by bond-order thresholds.  For rings and
other special coordinates, periodicity is not assumed unless the coordinate
definition supplies it.  Ring puckering can therefore be handled in two ways:

- use the protected puckering coordinates directly in a multidimensional DVR;
- compute derived \(Q,\phi\) descriptors in post-processing when their symmetry
  or periodicity is not suitable for constrained optimization.

For weak complexes and fragments, the coordinate policy is selected upstream:
fragment-center coordinates retain rigid-body relative translations and
orientations, whereas pseudo-bond/H-bond mode connects fragments and then uses
ordinary internal coordinates without creating artificial rings.

## Data Contract In MATRIX

GF now exports the following large-amplitude kinetic data:

- `g_inverse.csv`: the full global \(G^{-1}\) matrix in frozen GIC order;
- `large_amplitude_blocks.csv`: each selected large-amplitude block with its
  compact \(G^{-1}_{SS}\) block;
- `large_amplitude_g_inverse_blocks.csv`: the same \(G^{-1}_{SS}\) blocks in
  row/column long form;
- `#GF_PED [LARGE_AMPLITUDE_BLOCKS]`: restartable block definitions including
  `G_INV_BLOCK`;
- `#GF_PED [LARGE_AMPLITUDE_DVR_PLAN]`: per-coordinate status, periodicity,
  one-dimensional Fourier estimate when justified, and the diagonal
  \((G^{-1})_{ii}\) for quick inspection.

The diagonal field is diagnostic.  A multidimensional DVR must read the block.

## Validation Rules

A MATRIX large-amplitude treatment should fail or warn when any of these
conditions is violated:

- the frozen GIC rank is inconsistent with \(3N-6\) or \(3N-5\);
- \(F\) or \(G\) has non-negligible cross-irrep couplings in a symmetry-adapted
  basis;
- a large-amplitude block mixes coordinate families that the selected physical
  model says must remain separate;
- a 1D torsional estimate is requested despite strong \(F\) or \(G\) coupling
  to other large-amplitude coordinates;
- periodicity is undefined for a coordinate that is being placed on a periodic
  DVR grid;
- scan data use a coordinate definition different from the frozen MATRIX GIC.

## Development Consequences

The practical path is:

1. Let NEO define and freeze the coordinate families.
2. Let GF compute \(F\), \(G\) and global \(G^{-1}\) once.
3. Select large-amplitude blocks using coordinate family, frequency threshold,
   bond order, topology and \(F/G\) coupling diagnostics.
4. Pass \(F_{SS}\), \(G^{-1}_{SS}\), periodicity/range metadata and scan
   instructions to DVR.
5. Use DVR levels or partition functions in Thermo/Kinetics and spectroscopy.

This gives a direct extension of the Truhlar torsional idea to localized,
symmetry-aware and chemically protected MATRIX coordinates.  The key point is
that the kinetic metric is already available in GF; it should be reused, not
reconstructed by a separate rotor algorithm.
