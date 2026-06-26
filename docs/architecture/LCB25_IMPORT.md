# LCB25 Import

The LCB25 website provides a database of PCS2, SE and HPCS2 geometries. The
public page exposes three ZIP downloads:

- `PCS2.zip`
- `SE.zip`
- `HPCS2.zip`

ORACLE treats LCB25 as a remote geometry library source. The import flow is:

```text
LCB25 ZIP
  -> local library cache
  -> ORACLE-Babel import for each XYZ
  -> enriched XYZ with #SOURCE/#TOPOLOGY/#SYNTHONS
  -> fragment extraction or full-molecule reference search
```

The managed local cache is created with:

```bash
python -m oracle lcb25 fetch
```

By default it writes `data/lcb25/manifest.json`, `data/lcb25/archives/*.zip`
and `data/lcb25/xyz/<dataset>/*.xyz`. The cache directory is intentionally
ignored by git: code and manifests are reproducible, while redistribution of the
downloaded database remains an explicit project decision.

LCB25 molecules can be used in two directions:

- as whole-molecule references for MORPHEUS/reference-assisted workflows;
- as fragment libraries through `oracle-fragments` after ORACLE
  topology/synthon preprocessing.

Conversely, an arbitrary query molecule can be fragmented by the same
`#TOPOLOGY/#SYNTHONS` state and compared against LCB25-derived fragments.

The adapter intentionally starts with URL planning and local archive extraction.
Search/index metadata should be added after the downloaded XYZ naming and any
sidecar metadata files are inspected.
