# ORACLE GUI Architecture

The ORACLE GUI is a project dashboard over the shared enriched XYZ container.
It must not implement private chemistry, parser, fitting or Fortran logic.

The allowed call path is:

```text
GUI view -> oracle_gui controller -> oracle-* service/CLI -> xyzin sections
```

The forbidden call paths are:

```text
GUI view -> private XYZ/Gaussian/Z-matrix parser
GUI view -> private topology or synthon perception
GUI view -> direct Fortran executable without an oracle-* service boundary
```

## Project Windows

The first window is `ORACLE Project Dashboard`. It opens or creates the active
project, displays the molecule state and lists all known `xyzin` sections.
Every other window is reachable from this state.

The project/workflow windows are:

| Window | Responsibility | Main sections |
| --- | --- | --- |
| ORACLE Project Dashboard | Project state, validation, section/workflow status | all sections |
| Tool Contracts | Standalone tool readiness from the current `xyzin` sections | all tool-owned sections |
| ORACLE-Babel / Preprocessing | Import XYZ, QM formats, Z-matrix, SMILES/RDKit and LCB25 | `#SOURCE`, `#BASIC`, `#SYMMETRY`, `#TOPOLOGY`, `#SYNTHONS` |
| Molecule Editor / Avogadro Bridge | Open/edit the first XYZ block and reimport edited coordinates | XYZ block |
| Molecular Structure / Synthons | Inspect bonds, rings, charges, synthons and fragments | `#TOPOLOGY`, `#SYNTHONS`, `#FRAGMENTS` |
| GICForge | Build, symmetrize and diagnose GICs and B matrices | `#GIC`, `#SYCART` |
| GF / PED | Harmonic force-field analysis from Hessian plus GICs | `#CARTESIAN_HESSIAN`, `#GF_PED` |
| SEFit / MORPHEUS | Single-molecule and ensemble semiexperimental refinement | `#ISOTOPOLOGUES`, `#MORPHEUS` |
| TRINITY Geometry Optimization | Prepare external energy/gradient geometry optimization requests | `#TRINITY` |
| Rovib / Thermo Utilities | Rotational/vibrational summaries and thermochemistry utilities | `#ROTATIONAL`, `#VIBRATIONAL`, `#THERMO` |
| Anharmonic: VPT2 / VCI / DVR | Run and collect anharmonic workflow state | `#QFF`, `#VPT2_VCI`, `#DVR` |
| QM Jobs | Generate Gaussian inputs and normalize QM output sections | `#CARTESIAN_HESSIAN`, `#NORMAL_MODES`, `#QFF` |
| Diagnostics / Regression | Corpus, Python/Fortran and benchmark audits | reports/artifacts |

## Scientific Workbenches

The final user-facing GUI must also expose domain workbenches. These collect
the project workflow outputs, perform domain-specific analysis and draw spectra
or figures. They still use the same `xyzin` sections and service boundaries.

| Workbench | Responsibility | Examples |
| --- | --- | --- |
| Rotational Spectroscopy | Collect rotational constants, corrections and assignments; draw rotational spectra. | isotopologue comparison, SEFit residuals, line-list/stick/envelope plots, publication CSV/SVG/PDF/LaTeX export |
| Vibrational Spectroscopy | Collect harmonic, GF/PED, VPT2/VCI and DVR data; draw vibrational spectra. | normal-mode overlap heat maps, force-constant scaling, IR/Raman-style peak tables, publication plots |
| Electronic Spectroscopy | Collect electronic-state and transition data; visualize orbitals/densities. | UV/visible stick/broadened spectra, orbital viewing through Avogadro/Avogadro2/Molden, transition tables |
| Molecular Structure / Synthons | Inspect and publish structure, topology, fragments and synthon classifications. | structure tables, synthon maps, fragment maps, Avogadro editing |
| Thermochemistry / Kinetics | Collect thermochemical functions and kinetic-model outputs. | thermo tables, rovibrational DOS, kinetic comparisons, publication-ready plots/tables |

Publication export is a GUI contract: workbenches must export the final
accepted spectrum or table as machine-readable data plus vector formats suited
for papers. The first implementation records these export targets in
`oracle_gui.workflows`; plotting backends can then implement them without
changing the scientific services.

Orbital visualization is delegated to external viewers. ORACLE should prepare
or pass through supported files, then launch Avogadro/Avogadro2, Molden or
another configured viewer through `oracle_gui.commands`, rather than embedding
a private orbital renderer in the first GUI pass.

## Implementation Boundary

`oracle-gui` owns only:

- project and section view-models;
- GUI command specifications;
- optional Qt widgets;
- file selection and user interaction state.

Scientific behavior remains in the owning packages:

- `oracle-chem` for ORACLE-Babel, symmetry, topology and synthons;
- `oracle-gicforge` for GIC construction, symmetrization and B matrices;
- `oracle-gf` for GF/PED;
- `oracle-morpheus` for SEFit and MORPHEUS;
- `oracle-trinity` for external energy/gradient geometry optimization state;
- `oracle-rovib`, `oracle-thermo`, `oracle-vpt2-vci` and `oracle-dvr` for
  their corresponding sections.
- spectrum drawing and publication export code should consume normalized
  ORACLE reports/CSV/sections and must not parse QM outputs privately.

The source-tree entry point is:

```bash
python -m oracle_gui [molecule.xyzin]
```

The installable console entry point is:

```bash
oracle-gui [molecule.xyzin]
```

Qt is optional and loaded only when the GUI process starts. Headless tests cover
the controllers without requiring a display server.

## Dashboard Runtime

The dashboard runtime is split in two layers:

- `oracle_gui.dashboard.OracleDashboardController` is headless and testable. It
  loads the current project state, builds the available action list, checks
  required `xyzin` sections, assigns the repository working directory for
  `python -m oracle` commands and records process logs.
- `oracle_gui.app` is the optional Qt layer. It displays project state,
  sections, actions and logs, then launches selected commands through
  `QProcess` so long-running calculations do not block the desktop window.

The first operational dashboard actions are intentionally conservative:
validation, Avogadro launch, fragment build, GICForge build/report/B-matrix,
rovibrational summary, Thermo, VPT2/VCI collection and DVR collection. Actions
that need additional file choices, such as ORACLE-Babel import or FCHK
promotion, remain command specifications until their dedicated windows provide
the required file selectors.

## Tool Contracts Tab

The Tool Contracts tab uses `oracle_gui.contracts` and the central
`oracle_core.tool_contracts` registry. It lists each standalone tool, its
current name, future MATRIX-era name where defined, required `xyzin` sections,
missing sections and canonical CLI entry point.

The tab is read-only. It must not maintain a private readiness model or infer
scientific prerequisites from GUI state. The same readiness logic is exposed by
`python -m oracle contracts --check-xyzin molecule.xyzin`.

## Structure And ORACLE-Babel Tab

The first dedicated workflow tab is the Structure tab inside the dashboard. It
uses `oracle_gui.structure` and provides:

- source and output file selectors for ORACLE-Babel preprocessing;
- source-kind selection passed directly to `oracle babel preprocess`;
- Avogadro launch for the active `xyzin`;
- fragment build through the shared `oracle fragments build` CLI;
- read-only tables for saved `#TOPOLOGY` bonds/rings, `#SYNTHONS` rows and
  built `#FRAGMENTS`.

The Structure tab only displays saved ORACLE sections. It must not rediscover
bonds, rings, synthons or fragments in GUI code. Any refresh of molecular
state must go through `oracle-chem`, `oracle-babel` or `oracle-fragments`.

## GICForge Tab

The GICForge tab uses `oracle_gui.gicforge` and exposes the coordinate workflow
without duplicating GIC construction in Qt:

- build GICs through `oracle gicforge build`, with explicit controls for
  symmetrization, `#SYCART` generation and improper-dihedral out-of-plane mode;
- write the GICForge report, evaluate the analytic B matrix and generate a
  Gaussian input with the saved Gaussian GIC block;
- display the frozen `#GIC` definition as read-only tables for primitives,
  frozen GICs, symmetry projector groups and reduction/symmetry diagnostics.

The tab reads `#GIC` with `oracle_gicforge.read_gic_definition_from_xyzin`.
It must not construct primitives, symmetrize coordinates or evaluate B rows in
GUI code. Those operations remain in `oracle-gicforge`, because the same
utilities are required by optimizers and least-squares fitting at each geometry
iteration.

## GF / PED Tab

The GF/PED tab uses `oracle_gui.gf` and the `#GF_PED` section written by
`oracle gf`. It provides:

- optional Gaussian FCHK selection when the project does not already contain a
  `#CARTESIAN_HESSIAN` section;
- report and CSV output locations;
- controls for symmetry blocks, local force-field filtering, Pulay scaling,
  force-constant thresholding, electrostatic/UFF-vdW subtraction and 1-4
  scaling;
- read-only tables for frequencies, GIC labels/irreps, dominant PED
  contributions, the full PED matrix and run diagnostics.

The tab must not parse Gaussian/FCHK data, reconstruct Hessians or solve GF
itself. It launches `oracle gf`, then reloads `#GF_PED` through
`oracle_gf.read_gf_ped_section`.

## SEFit / MORPHEUS Tab

The SEFit tab uses `oracle_gui.sefit` and the `#MORPHEUS` section written by
`oracle semiexp`. It provides:

- SEFit/MORPHEUS job selection, output directory selection and backend choice;
- explicit control over whether the active `xyzin` is updated with
  `#MORPHEUS`;
- optional CLI arguments for advanced MORPHEUS controls such as iteration
  limits, coordinate model, robust loss and constraints;
- read-only tables for `#ISOTOPOLOGUES`, output artifacts and fit diagnostics.

The tab must not run least-squares logic or parse MORPHEUS CSV reports itself.
It launches `oracle semiexp`, then reloads `#MORPHEUS` through
`oracle_morpheus.read_morpheus_section`.

## TRINITY Tab

The TRINITY tab uses `oracle_gui.trinity` and the `#TRINITY` section written by
`oracle trinity prepare`. It provides:

- external engine command entry;
- run-directory selection;
- coordinate model and active-space selection;
- trust-region and convergence tolerance fields;
- read-only tables for saved settings and expected output paths.

The tab is currently a skeleton for the future optimizer branch. It must not
implement geometry optimization or external energy/gradient parsing in Qt. It
only prepares the autonomous `#TRINITY` request and then reloads it through
`oracle_trinity.read_trinity_section`.

## Generic Workbench Tabs

The Rovib/Thermo, Anharmonic, QM Jobs, Diagnostics, Rotational, Vibrational,
Electronic and Thermo/Kinetics tabs use `oracle_gui.workbench`. They expose the
central `oracle_gui.workflows.WindowSpec` contract as four read-only tables:
required/produced sections, available actions, capabilities and publication or
viewer exports.

These tabs are intentionally thin until their plotting or launcher controllers
are implemented. They must not duplicate parser, projector, GF, SEFit,
VPT2/VCI, DVR, thermo or spectroscopy logic. A workbench becomes operational by
adding or improving the owning service/CLI and its `xyzin` sections, then
teaching the GUI to call that service.
