# GICForge Legacy Inventory

Date: 2026-06-26

## Scope

This note records the Merlino GICForge material that ORACLE must preserve while
the implementation is split into shared libraries and independent tools.

## Legacy Sources

- `/Users/vincenzobarone/merlino3.0/fortran/gicforge`
  - `mkprim.f`: stretch, bend, linear-bend, torsion and out-of-plane primitive
    generation.
  - `mkcyc.f`, `cycdat.f`: ring data and ring-coordinate construction.
  - `gicprune.f`, `locsvd.f`: rank pruning and local SVD selection.
  - `mksalc.f`, `symm.f`, `symang.f`, `symdih.f`,
    `gic_type_symmetry.f`: symmetry-adapted GICs and SYCART support.
  - `rd_xyz.f`, `coord.f`, `bondout.f`, `pcsgeo.f`, `dina25.f`,
    `tools1.f`, `tools2.f`: input, geometry and driver support.
- `/Users/vincenzobarone/merlino3.0/merlino_gic`
  - `model.py`: frozen Python `GICDefinition`, Gaussian `gauin` parser, and
    service wrapper.
  - `gicforge_service.py`: Fortran executable runner and generated-file
    collection.
  - `gicforge_python.py`: Python GIC generator tied to Merlino topology,
    primitive and B-matrix modules.
  - `gic_symmetry.py`: Python post-check and promotion of symmetrized `gauin`.

## ORACLE Decisions

- ORACLE-GICForge must consume the shared `xyzin` state: `#VALIDATION`,
  `#TOPOLOGY`, `#SYNTHONS` and `#SYMMETRY`. It must not rediscover topology in a
  private parser.
- The Fortran77 backend remains intentional and should later be wrapped behind
  the same frozen `#GIC` and `#SYCART` contracts.
- The initial ORACLE-native backend generates primitive candidates from the
  saved topology, selects a non-redundant set by finite-difference B-matrix rank,
  and freezes the result in the enriched XYZ file.
- Full non-C1 symmetry adaptation and ring-specialized coordinates remain
  porting targets from the legacy Fortran/Python stack. Until those projectors
  are migrated, non-C1 `--symmetrize` requests must stop cleanly.

## Porting Order

1. Freeze a real `#GIC` schema in ORACLE with primitive definitions, GIC rows and
   Gaussian ReadGIC lines.
2. Freeze `#SYCART` as an external-mode-free Cartesian basis for the same
   validated molecule.
3. Route `oracle gicforge build` and Gaussian input generation through those
   frozen sections.
4. Replace the native primitive backend with the faithful Fortran77/Python
   GICForge backends without changing downstream contracts.
5. Add corpus regressions for Z-matrix, SMILES/RDKit, ring and symmetry cases.
