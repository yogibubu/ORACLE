# Pyrrole Gaussian ReadAllGIC Example

This directory contains a successful Gaussian DV J32 run generated from MATRIX
NEO/GICForge symmetrized GICs.

Files:

- `pyrrole.gjf`: Gaussian input using `opt=readallgic freq`.
- `pyrrole.log`: Gaussian output showing successful optimization and frequency
  calculation.

The example demonstrates that:

- Gaussian reads the MATRIX-generated ReadAllGIC block.
- The full 24-coordinate nonredundant basis is present (`NTRed=24`,
  `NRank=24`).
- Only the totally symmetric `A1` coordinates remain active.
- Non-totally-symmetric coordinates, including `B2` stretches/bends,
  `A2/B1` ring puckering and `A2/B1` out-of-plane coordinates, are marked
  `Frozen`.
- Gaussian preserves `C2V` symmetry, finds a stationary point and reports
  `NImag=0`.
