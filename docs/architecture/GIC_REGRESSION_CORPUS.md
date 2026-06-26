# GIC Regression Corpus

`tests/fixtures/test_molecules` is the imported Merlino `test_molecules`
corpus. It contains demanding molecular inputs for GICForge, topology, ring,
out-of-plane, Z-matrix and Gaussian-adapter regressions.

The corpus is intentionally versioned because it is small enough for the ORACLE
repository and because GIC regressions need stable input files. Runtime caches
and generated project outputs still belong outside git.

Current import source:

```text
/Users/vincenzobarone/test_molecules
```

The copied corpus excludes local OS artifacts such as `.DS_Store`, but keeps
scientific `.log` and `.out` files when they are part of the regression input
or adapter baseline.

Use this corpus as an input library, not as a second implementation:

1. ORACLE-Babel parses each legacy input into the canonical enriched XYZ.
2. Topology and symmetry are computed once and saved as sections.
3. GICForge generates and freezes the GIC definition.
4. Python and Fortran77 backends are compared against the same frozen contract.
5. GF/PED, SEFit/MORPHEUS, Thermo and anharmonic tools consume the same file
   without rebuilding GICs unless an explicit restart requests it.

Quick inventory:

```bash
source /Users/vincenzobarone/ORACLE/scripts/oracle_env.sh
oracle-gic-corpus-list
```
