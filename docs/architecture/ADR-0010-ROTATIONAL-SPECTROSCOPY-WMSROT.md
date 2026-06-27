# ADR-0010: Rotational Spectroscopy And WMS-Rot

## Status

Accepted.

## Context

WMS-Rot is available at <https://www.skies-village.it/webtools/wmsrot/> and
already implements rotational spectrum simulation, WMS-Prep, WMS-FitRot and
wavefunction post-processing workflows. The page is a browser application: it
loads Pyodide, Plotly and the WMS-Rot Python engine into the user's browser and
states that calculations are executed locally, without a remote calculation
server.

The Hamiltonian and fitting code were developed in the group. Rewriting that
work inside MATRIX would be wasteful and would create a second implementation
to validate. The architectural issue is not the WMS-Rot science code; it is the
boundary. MATRIX is organized around one enriched `xyzin` container, one
adapter for each external format, and scientific modules that can run
standalone when their required `xyzin` sections are present.

## Decision

MATRIX will copy the first-party WMS-Rot code locally and reuse its Python
Hamiltonian engine as the initial native rotational spectroscopy backend.
WMS-Rot is therefore kept in two forms:

- a vendored Python engine used by `oracle-rovib`;
- a first-party browser-suite snapshot for WMS-Rot/WMS-Prep/WMS-FitRot/
  WMS-Wavefunctions reference;
- an external browser reference;
- a compatibility target for generated input files;
- a numerical comparison source when future extensions are added.

MATRIX will not depend on the public WMS-Rot URL for production calculations.
The GUI may still open the WMS-Rot page, and `oracle rovib wmsrot-input`
exports a WMS-Rot input file from normalized `xyzin` sections.

The browser snapshot lives under `external/wmsrot-site/`. The callable Python
engine is copied to `oracle_rovib.vendor.wmsrot_engine` and is invoked through
`oracle_rovib.wmsrot`. Vendored browser JavaScript must not become a private
parser or data model for MATRIX; shared adapters still own Gaussian/QM parsing.

## Contract

The internal source of truth remains:

- `#ROTATIONAL` for rotor type, representation, point group, Watson reduction,
  rotational constants, temperature, dipole components and optional `Q_rot`;
- `#DELTABVIB` for rovibrational corrections;
- `#ROTATIONAL_SPECTRUM` or equivalent workflow outputs for generated line
  lists, assignment tables and broadened/publication spectra.

The WMS-Rot bridge is an export adapter:

```bash
oracle rovib wmsrot-input molecule.xyzin --out molecule.wmsrot.txt
```

The local WMS-Rot engine is a run backend:

```bash
oracle rovib wmsrot-run molecule.xyzin --out molecule.rotational.csv
```

The adapter/backend must not parse Gaussian logs, infer topology, or duplicate
QM reader logic. Missing hyperfine constants or centrifugal distortion terms
are exported or passed as zeros until those data are represented in shared
MATRIX sections.

## Consequences

- The Rotational Spectroscopy window can immediately open WMS-Rot, export a
  compatible input file and run the local WMS-Rot Hamiltonian.
- Reproducible MATRIX calculations use copied first-party code inside
  `oracle-rovib`, not an external website call.
- Future extensions, including internal hindered rotors, should extend this
  local backend rather than starting a second rotational Hamiltonian stack.
