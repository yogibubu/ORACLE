# MATRIX Topology Contract

LINK owns topology perception. Downstream tools must consume the frozen `xyzin`
sections and must not rebuild a private graph unless an explicit restart asks
LINK to regenerate topology from a new geometry.

## Sections

- `#TOPOLOGY` uses `SCHEMA matrix.xyz.topology.v1`.
- `#SYNTHONS` uses `SCHEMA matrix.xyz.synthons.v1`.
- `#FRAGMENTS` uses `SCHEMA matrix.xyz.fragments.v1` when present.
- `oracle.xyz.*.v1` schemas remain accepted as compatibility aliases for
  imported first-release files.

`#TOPOLOGY` stores one-based atom indexes, the frozen bond graph, one bond
order for every stored bond, canonical elementary ring records and aromaticity
annotations. Bond orders are sourced from Gaussian Mayer data when available.
When the source is MOL, SDF or MOL2, explicit file bonds are preserved and used
as topology overrides. Otherwise LINK uses the shared Pauling continuous
topology model. Charges are sourced from Gaussian CM5 when available; otherwise
they come from the synthon electronegativity model.

Ring perception is deliberately not "all simple cycles". LINK first enumerates
chordless cycles in the non-metal covalent graph, then selects a deterministic
minimum cycle basis over GF(2). This removes perimeter cycles in fused PAHs,
diagonal cycles in cages and metal-containing artificial cycles in metallocenes.
Metal-ring interactions are represented later by interaction centers and
special coordinates, not by adding metal atoms to `#TOPOLOGY [RINGS]`.

The ring-basis policy is written into every generated topology section:

```text
RING_BASIS_POLICY CHORDLESS_NONMETAL_MINIMUM_CYCLE_BASIS
RING_CANDIDATE_COUNT ...
RING_BASIS_RANK ...
RING_BASIS_COUNT ...
RING_BASIS_ALLOWED_ATOMS ...
RING_BASIS_ALLOWED_EDGES ...
RING_BASIS_EXCLUDED_ATOMS ...
```

These diagnostics are part of the topology golden hash, so a change in ring
candidate enumeration, rank, selected basis size or metal exclusion is visible
before NEO/GICForge sees the graph.

## Geometry Inputs

LINK can currently create the initial geometry from:

- plain XYZ;
- Gaussian input/log-style geometries;
- Gaussian formatted checkpoints (`.fchk`, `.fch`);
- MOL and SDF V2000 records;
- Tripos MOL2 files;
- SMILES through the RDKit adapter.

SDF/MOL/MOL2 are structure-file adapters, not independent topology engines.
They provide atoms, Cartesian coordinates and optional explicit bond orders;
LINK still owns the normalized topology, synthons, fragments and validation
sections written to `xyzin`.

## Validation

`matrix validate` now checks topology content, not only section presence:

- bond and bond-order indexes must be valid;
- every frozen bond must have a bond order;
- nontrivial molecules must not contain isolated topology atoms;
- ring edges must be present in the frozen bond graph;
- declared fragments must use valid atom indexes and cover the molecule.

Weak-complex policy remains explicit. Pseudo-bonds may connect fragments for
internal-coordinate construction, but they must not silently create artificial
ring topology. Fragment special coordinates and pseudo-bonds are therefore
separate policy choices above the frozen graph.

## Golden Snapshots

Topology has its own golden corpus:

```text
tests/fixtures/golden_corpus/matrix_topology_snapshots.json
```

The snapshot hash covers bonds, bond orders, elementary rings, ring relations,
aromaticity and connected fragments. A small readable subset of ring relations
is kept in the file for diagnosis; the full relation set is still protected by
the hash.

Create a snapshot or report with:

```bash
python -m matrix topology snapshot molecule.xyzin topology_snapshot.json
python -m matrix topology report molecule.xyzin topology_report.txt
```

NEO/GICForge, GF, MORPHEUS and future fragment/nano-LEGO tools should all use
this same topology contract rather than duplicating graph perception.
