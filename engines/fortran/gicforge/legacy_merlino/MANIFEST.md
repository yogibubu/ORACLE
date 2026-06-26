# Merlino GICForge Legacy Source Import

Date: 2026-06-26

Source root:
`/Users/vincenzobarone/merlino3.0/fortran/gicforge`

This directory vendors the strict Fortran77 GICForge backend from Merlino3.0 so
ORACLE can compare Python and Fortran implementations against one frozen GIC
contract. Files are kept as source material and should not be edited casually.
ORACLE-specific Fortran kernels live one directory above this one.

## Imported Runtime Files

- `dina25.f`: main legacy GICForge driver.
- `coord.f`: Cartesian/Z-matrix reader and connectivity setup.
- `mkprim.f`: primitive generation and B-matrix routines.
- `mkcyc.f`: ring and puckering coordinate generation.
- `mksalc.f`: local SALC/GNIC construction helpers.
- `gicprune.f`: legacy pruning and block packing.
- `gic_type_symmetry.f`: type-local GIC symmetry helpers.
- `locsvd.f`: local SVD/Jacobi helper.
- `pcsgeo.f`: internal-to-Cartesian reconstruction helpers.
- `entors.f`, `symang.f`, `tools1.f`, `tools2.f`, `symm.f`: support kernels.

## Imported Data And Reference Files

- `gicforge_keywords`
- `masses.txt`
- `bondout.f`
- `cycdat.f`
- `rd_xyz.f`
- `symdih.f`

The active ORACLE build wrapper is `../compile_legacy`; it mirrors the Merlino
`compile_MAC` source list but uses ORACLE-local include paths.
