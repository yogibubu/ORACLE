Fortran77 VPT2/VCI Core
=======================

This directory contains standalone Fortran77 kernels for the ORACLE
VPT2/VCI backend. They are not copied from GDV.

- `gf_core.f`: small symmetric Jacobi kernel and a Wilson-GF helper for already
  independent coordinates.
- `vci_core.f`: product-basis generation and dense small-space VCI helpers.
  - `M4VCIBasis` keeps the simple total-quanta basis contract.
  - `M4VCIBasisCtl` adds per-mode min/max quanta and min/max total quanta for
    one-, two-, three- and four-mode excitation classes.
- `vpt2_core.f`: small-space VPT2 helper using quartic first-order and cubic
  second-order corrections on a supplied basis.
- `davidson_core.f`: independent Davidson support routines based only on
  residual vectors, diagonals and orthogonalization.

Gaussian parsing, coordinate construction and tensor normalization are handled
by Python. The Fortran code receives only numerical arrays and shares the same
solver contracts as the Python implementation.

Use:

```bash
cd engines/fortran/vpt2_vci
./compile_check
```
