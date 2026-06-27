# WMS-Rot Site Snapshot

This directory contains a first-party snapshot of the WMS-Rot browser suite
from:

```text
https://www.skies-village.it/webtools/wmsrot/
```

Snapshot date: 2026-06-27.

The scientific Hamiltonian engine is also copied into:

```text
packages/oracle-rovib/src/oracle_rovib/vendor/wmsrot_engine.py
```

MATRIX uses that vendored Python engine through `oracle_rovib.wmsrot`; the
HTML/JavaScript files here are retained as a reference for WMS-Rot input,
WMS-Prep, WMS-FitRot and wavefunction workflows.

Excluded from this snapshot:

- CDN libraries such as Plotly, MathJax and Pyodide;
- large browser runtime/vendor bundles;
- image/icon assets not needed for scientific porting.
