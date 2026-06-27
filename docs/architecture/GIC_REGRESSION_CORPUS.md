# GIC Regression Corpus

`tests/fixtures/test_molecules` is the imported MATRIX `test_molecules`
corpus. It contains demanding molecular inputs for GICForge, topology, ring,
out-of-plane, Z-matrix and Gaussian-adapter regressions.

The corpus is intentionally versioned because it is small enough for the MATRIX
repository and because GIC regressions need stable input files. Runtime caches
and generated project outputs still belong outside git.

Current import source:

```text
/Users/vincenzobarone/test_molecules
```

The copied corpus excludes local OS artifacts such as `.DS_Store`, but keeps
scientific `.log` and `.out` files when they are part of the regression input
or adapter baseline.

## Official Golden Registry

The official first-release NEO/GICForge golden set is machine-readable:

```text
tests/fixtures/golden_corpus/neo_gic_golden.json
```

This registry is the stable gate for molecules that must not drift without an
explicit review. It currently anchors benzene, pyrrole, azulene, the fused PAH
set, norbornane/norbornene/norbornadiene/norcamphor, thujone, ribose, cubane,
ferrocene in D5h and D5d forms, spiro rings, cyclottane and the
formic-acid--water H-bond complex. Each entry records the source file and the
scientific role it protects: ring coordinates, fused/polycyclic behavior,
bridged rings, special metal/ring-center coordinates, H-bond pseudo-bonds,
fragment coordinates, symmetry projectors and Python/Fortran parity.

Inputs in the registry are deliberately versioned. Larger generated audit
workdirs, numerical reports and run outputs are not committed; their provenance
belongs in run manifests.

Use this corpus as an input library, not as a second implementation:

1. LINK parses each legacy input into the canonical enriched XYZ.
2. Topology and symmetry are computed once and saved as sections.
3. GICForge generates and freezes the GIC definition.
4. Python and Fortran77 backends are compared against the same frozen contract.
5. GF/PED, SEFit/MORPHEUS, Thermo and anharmonic tools consume the same file
   without rebuilding GICs unless an explicit restart requests it.

Quick inventory:

```bash
source /Users/vincenzobarone/MATRIX/scripts/matrix_env.sh
matrix-gic-corpus-list
```

Numerical Python/Fortran77 parity audit:

```bash
python -m matrix gicforge fortran-audit \
  --root tests/fixtures/test_molecules/molecules \
  --workdir runs/gicforge_fortran_audit
```

The audit preprocesses each selected corpus molecule through LINK,
builds the frozen MATRIX GIC/B matrix, runs the vendored Merlino Fortran77
GICForge harness, and compares final rank, row-space rank and Wilson-B row
span. Add `--molecule name.inp` repeatedly to pin a smaller periodic audit set.

Current required gate:

- The default audit set is pyrrole, benzene, pyridine, pyrimidine, naphthalene,
  phenanthrene, anthracene, pyrene, fluorene, azulene, norbornane, norbornene,
  norbornadiene, norcamphor, spiro, acetylene, linear C4S, thujone, ribose,
  cubane and cyclottane.
- The default audit is tied to the golden registry roles. Every entry tagged
  `fused_ring`, `bridged_ring`, `spiro_ring` or `python_fortran_parity` must
  remain in `DEFAULT_FORTRAN_AUDIT_MOLECULES` unless the registry is explicitly
  reviewed.
- This gate must pass with matching final rank, matching row-space rank and
  Wilson-B row-space residual below the audit tolerance. The finite
  point-group projector must be active for every non-`C1`, non-linear molecule
  in the gate.
- The audit summary reports projector status, symmetry block counts, mixed
  symmetry-family counts, total-symmetric GIC counts, nontrivial SALC
  coefficient counts and the largest SALC coefficient normalization residual.
  Exact coefficient-vector comparison with executable Merlino is still a TODO
  until the strict Fortran backend emits projector coefficients in a stable
  machine-readable form.
- Pyrrole is an explicit regression: point group `C2v`, rank 24, retained ring
  coordinates, and symmetrized `A1RPck001`/`B2RPck001`.
- Fused PAHs are explicit regressions: naphthalene, phenanthrene, anthracene,
  pyrene and fluorene must retain protected ring `RDef`/`RPck` sources, keep
  `BtFl` bridge coordinates when present, and export fused ring puckerings as
  `RPck` rather than polar `QPck`/`PhiP` functionals.
- Non-benzenoid and bridged rings are explicit regressions: azulene must keep a
  `BtFl` bridge and point-group projector symmetry, while norbornane,
  norbornene, norbornadiene and norcamphor must retain the shared-cycle `RDef`
  and `RPck` source spaces without falling back in the Fortran harness.
- Spiro rings are an explicit regression: `spiro.inp` must detect point group
  `D2`, build the Merlino `Spir` inter-ring angle block around the shared atom,
  retain all six `RPck` ring components, and close the point-group projector
  without falling back to local SALCs.
- Linear molecules are explicit regressions: `c2h2.inp` and `c4s.inp` must run
  the legacy Merlino harness through the `ECKART` path so that Fortran sets its
  linear-top flag and retains the physical `3N-5` rank with both components of
  each linear bend.
- Thujone is an explicit regression for a larger nonaromatic bridged aliphatic
  molecule: `thujone.inp` must keep the Merlino-equivalent final rank 75 and
  Wilson-B row space.
- Ribose is an explicit regression for a flexible oxygenated sugar ring:
  `ribose.inp` must keep the Merlino-equivalent final rank 51, preserve the
  ring `RDef`/`RPck`/`QPck`/`PhiP` coordinate families, and match the Wilson-B
  row space.
- Cubane is an explicit high-symmetry regression: `cubane.inp` must detect
  point group `Oh`, keep the Merlino-equivalent final rank 42, use the full
  point-group projector, keep `BUTTERFLY` and `CYCLIC_BEND` as separate
  homogeneous blocks, and preserve Wilson-B rank 42 after symmetrization.
- Ferrocene is an explicit metal-center projector regression in the Python
  corpus: eclipsed `ferrocene.inp` must detect `D5h`, staggered
  `ferrocene_staggered.inp` must detect `D5d`, both must serialize 20 closed
  finite operations, protect the two Fe-to-Cp-ring-center coordinates, keep the
  full point-group projector active, and preserve Wilson-B rank 57.
- Cyclottane is an explicit large noncondensed ring regression: `cyclottane.inp`
  must detect point group `D2`, keep the Merlino-equivalent final rank 66, use
  the full point-group projector, preserve five `RDef`/`RPck` components, and
  export the expected polar `QPck`/`PhiP` ring-puckering functionals.

Full-corpus status from the 127 imported `.inp` files is not yet a release
gate. The current broad audit separates:

- passing cases;
- input/backend errors that must be fixed in LINK, topology or legacy harness
  normalization before they can be used as GIC parity tests;
- real GIC parity mismatches.

Open GIC parity triage:

- `biphenylene.inp`: nonredundant rank and B span match Merlino, but the D2h
  point-group projector falls back because the generated RPck candidate space
  is not closed under all D2h operations.
- `9cyanophenantrene.inp` and `cyanopyridine.inp`: rank matches, but the
  Python nonredundant space differs from Merlino because the ring-puckering /
  polar `QPck/PhiP` selection is not yet Merlino-equivalent.
- The remaining residual-only cases near the tolerance boundary must be
  reviewed after the structural mismatches above are fixed.
