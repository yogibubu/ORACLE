from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class WorkflowStatus(str, Enum):
    MISSING = "missing"
    READY = "ready"
    WARNING = "warning"
    COMPLETE = "complete"


@dataclass(frozen=True)
class WorkflowActionSpec:
    key: str
    label: str
    command: str
    required_sections: tuple[str, ...] = ()
    produced_sections: tuple[str, ...] = ()


@dataclass(frozen=True)
class WindowSpec:
    key: str
    title: str
    description: str
    category: str = "workflow"
    required_sections: tuple[str, ...] = ()
    produced_sections: tuple[str, ...] = ()
    capabilities: tuple[str, ...] = ()
    publication_exports: tuple[str, ...] = ()
    external_viewers: tuple[str, ...] = ()
    actions: tuple[WorkflowActionSpec, ...] = ()


ORACLE_GUI_WINDOWS: tuple[WindowSpec, ...] = (
    WindowSpec(
        key="dashboard",
        title="ORACLE Project Dashboard",
        description="Project state, molecule preview, xyzin sections and workflow status.",
        category="project",
        capabilities=(
            "open/create ORACLE projects",
            "inspect xyzin sections",
            "show workflow readiness",
        ),
        actions=(
            WorkflowActionSpec(
                key="validate",
                label="Validate molecule",
                command="validate",
            ),
        ),
    ),
    WindowSpec(
        key="babel",
        title="ORACLE-Babel / Preprocessing",
        description="Import external geometries and materialize the shared enriched XYZ state.",
        category="project",
        produced_sections=("SOURCE", "BASIC", "SYMMETRY", "TOPOLOGY", "SYNTHONS"),
        capabilities=(
            "import XYZ, Z-matrix, QM geometry and SMILES/RDKit sources",
            "run symmetry detection once with explicit thresholds",
            "materialize topology and synthons for all downstream tools",
        ),
        actions=(
            WorkflowActionSpec(
                key="preprocess",
                label="Import geometry",
                command="babel preprocess",
            ),
        ),
    ),
    WindowSpec(
        key="avogadro",
        title="Molecule Editor / Avogadro Bridge",
        description="Open the XYZ block in Avogadro and import edited coordinates.",
        category="structure",
        required_sections=("SOURCE",),
        external_viewers=("Avogadro", "Avogadro2"),
        capabilities=(
            "view/edit the first XYZ block externally",
            "reimport edited Cartesian coordinates",
            "invalidate or refresh dependent xyzin sections after geometry edits",
        ),
        actions=(
            WorkflowActionSpec(
                key="replace_xyz",
                label="Import edited XYZ block",
                command="replace xyz block",
            ),
        ),
    ),
    WindowSpec(
        key="topology",
        title="Molecular Structure / Synthons",
        description="Inspect molecular structure, synthons, topology and fragment libraries.",
        category="structure",
        required_sections=("TOPOLOGY", "SYNTHONS"),
        produced_sections=("FRAGMENTS",),
        capabilities=(
            "inspect bonds, rings, charges and bond orders",
            "compare and edit synthon classes",
            "fragment molecules and build reusable fragment libraries",
            "prepare nano-lego fragment assembly workflows",
        ),
        publication_exports=("structure tables", "synthon tables", "fragment maps"),
        external_viewers=("Avogadro", "Avogadro2"),
        actions=(
            WorkflowActionSpec(
                key="fragment_plan",
                label="Plan fragments",
                command="fragments plan",
                required_sections=("TOPOLOGY", "SYNTHONS"),
                produced_sections=("FRAGMENTS",),
            ),
            WorkflowActionSpec(
                key="fragment_build",
                label="Build fragments",
                command="fragments build",
                required_sections=("TOPOLOGY", "SYNTHONS"),
                produced_sections=("FRAGMENTS",),
            ),
        ),
    ),
    WindowSpec(
        key="gicforge",
        title="GICForge",
        description="Build, symmetrize and diagnose frozen GICs and B matrices.",
        category="project",
        required_sections=("SYMMETRY", "TOPOLOGY", "SYNTHONS"),
        produced_sections=("GIC", "SYCART"),
        capabilities=(
            "build primitive and special internal coordinates",
            "symmetrize GICs without mixing coordinate families",
            "evaluate B matrices on the current geometry",
            "write Gaussian input using the frozen GIC contract",
        ),
        actions=(
            WorkflowActionSpec(
                key="gic_build",
                label="Build GICs",
                command="gicforge build",
                required_sections=("SYMMETRY", "TOPOLOGY", "SYNTHONS"),
                produced_sections=("GIC",),
            ),
            WorkflowActionSpec(
                key="gic_bmatrix",
                label="Evaluate B matrix",
                command="gicforge bmatrix",
                required_sections=("GIC",),
            ),
            WorkflowActionSpec(
                key="gic_gaussian_input",
                label="Write Gaussian input",
                command="gicforge gaussian-input",
                required_sections=("GIC",),
            ),
        ),
    ),
    WindowSpec(
        key="gf",
        title="GF / PED",
        description="Run harmonic GF/PED from a Cartesian Hessian and frozen GICs.",
        category="vibrational",
        required_sections=("GIC", "CARTESIAN_HESSIAN"),
        produced_sections=("GF_PED",),
        capabilities=(
            "solve harmonic vibrational GF/PED models",
            "scale internal force constants",
            "apply local force-field filtering and nonbonded subtractions",
            "export mode/PED diagnostics for the vibrational workbench",
        ),
        actions=(
            WorkflowActionSpec(
                key="gf_run",
                label="Run GF/PED",
                command="gf",
                required_sections=("GIC", "CARTESIAN_HESSIAN"),
                produced_sections=("GF_PED",),
            ),
        ),
    ),
    WindowSpec(
        key="sefit",
        title="SEFit / MORPHEUS",
        description="Fit semiexperimental structures and multi-molecule MORPHEUS models.",
        category="rotational",
        required_sections=("ISOTOPOLOGUES",),
        produced_sections=("MORPHEUS",),
        capabilities=(
            "fit semiexperimental structures from rotational constants",
            "compare observed, corrected and calculated rotational constants",
            "run ensemble/multistructure MORPHEUS refinements",
        ),
        actions=(
            WorkflowActionSpec(
                key="semiexp_fit",
                label="Run SEFit",
                command="semiexp",
                required_sections=("ISOTOPOLOGUES",),
                produced_sections=("MORPHEUS",),
            ),
            WorkflowActionSpec(
                key="semiexp_benchmark",
                label="Run paper benchmark",
                command="semiexp-benchmark",
            ),
        ),
    ),
    WindowSpec(
        key="rovib_thermo",
        title="Rovib / Thermo Utilities",
        description="Inspect rotational and vibrational sections and run thermochemistry utilities.",
        category="thermochemistry",
        required_sections=("BASIC", "ROTATIONAL"),
        produced_sections=("THERMO",),
        capabilities=(
            "summarize rotational and vibrational xyzin sections",
            "build vibrational and rovibrational density of states",
            "run thermochemistry from normalized ORACLE sections",
        ),
        actions=(
            WorkflowActionSpec(
                key="thermo_run",
                label="Run Thermo",
                command="thermo",
                required_sections=("BASIC", "ROTATIONAL"),
                produced_sections=("THERMO",),
            ),
        ),
    ),
    WindowSpec(
        key="anharmonic",
        title="Anharmonic: VPT2 / VCI / DVR",
        description="Prepare, run and collect VPT2/VCI and DVR workflow state.",
        category="vibrational",
        produced_sections=("VPT2_VCI", "DVR"),
        capabilities=(
            "run VPT2/VCI from normalized QFF data",
            "run DVR directly or collect post-run output",
            "compare anharmonic origins for spectral assignment",
        ),
        actions=(
            WorkflowActionSpec(
                key="vpt2_vci_run",
                label="Run VPT2/VCI",
                command="vpt2-vci",
                required_sections=("QFF",),
                produced_sections=("VPT2_VCI",),
            ),
            WorkflowActionSpec(
                key="dvr_run",
                label="Run DVR",
                command="dvr run",
                required_sections=("DVR",),
                produced_sections=("DVR",),
            ),
        ),
    ),
    WindowSpec(
        key="qm_jobs",
        title="QM Jobs",
        description="Generate inputs, monitor external jobs and normalize QM outputs.",
        category="project",
        produced_sections=("CARTESIAN_HESSIAN", "NORMAL_MODES", "QFF"),
        capabilities=(
            "write Gaussian inputs",
            "promote FCHK Hessian, normal-mode and QFF data",
            "promote rovibrational output sections from Gaussian logs",
            "monitor external QM jobs through service adapters",
        ),
        external_viewers=("Avogadro", "Molden"),
        actions=(
            WorkflowActionSpec(
                key="gaussian_from_gic",
                label="Gaussian from GICs",
                command="gicforge gaussian-input",
                required_sections=("GIC",),
            ),
        ),
    ),
    WindowSpec(
        key="diagnostics",
        title="Diagnostics / Regression",
        description="Run numerical audits against corpus, Fortran77 and benchmark fixtures.",
        category="diagnostics",
        capabilities=(
            "audit GICForge Python/Fortran equivalence",
            "audit demanding GIC regression corpora",
            "run paper benchmark commands",
        ),
        actions=(
            WorkflowActionSpec(
                key="gic_fortran_audit",
                label="GICForge Python/Fortran audit",
                command="gicforge fortran-audit",
            ),
            WorkflowActionSpec(
                key="gic_corpus_audit",
                label="GIC regression corpus audit",
                command="gicforge corpus-audit",
            ),
        ),
    ),
    WindowSpec(
        key="rotational_spectroscopy",
        title="Rotational Spectroscopy",
        description="Collect rotational calculations, inspect assignments and draw publication spectra.",
        category="spectroscopy",
        required_sections=("ROTATIONAL",),
        capabilities=(
            "collect rotational constants and rovibrational corrections",
            "compare isotopologues and semiexperimental residuals",
            "simulate rotational stick/envelope spectra",
            "scale or select rotational components for planar/asymmetric systems",
        ),
        publication_exports=("CSV line list", "SVG spectrum", "PDF spectrum", "LaTeX table"),
        actions=(
            WorkflowActionSpec(
                key="rotational_summary",
                label="Summarize rotational state",
                command="rovib summarize",
                required_sections=("ROTATIONAL",),
            ),
            WorkflowActionSpec(
                key="semiexp_fit",
                label="Run SEFit",
                command="semiexp",
                required_sections=("ISOTOPOLOGUES",),
                produced_sections=("MORPHEUS",),
            ),
        ),
    ),
    WindowSpec(
        key="vibrational_spectroscopy",
        title="Vibrational Spectroscopy",
        description="Collect harmonic/anharmonic modes, compare origins and draw spectra.",
        category="spectroscopy",
        required_sections=("VIBRATIONAL",),
        produced_sections=("GF_PED", "VPT2_VCI", "DVR"),
        capabilities=(
            "compare normal modes from GF, Gaussian, VPT2/VCI and DVR",
            "draw heat maps of normal-mode overlaps or contributions",
            "scale force constants and compare frequency shifts",
            "simulate IR/Raman-style stick and broadened spectra",
        ),
        publication_exports=("CSV peak table", "SVG spectrum", "PDF spectrum", "mode heat-map"),
        actions=(
            WorkflowActionSpec(
                key="vibrational_dos",
                label="Build vibrational DOS",
                command="rovib dos",
                required_sections=("VIBRATIONAL",),
            ),
            WorkflowActionSpec(
                key="gf_run",
                label="Run GF/PED",
                command="gf",
                required_sections=("GIC", "CARTESIAN_HESSIAN"),
                produced_sections=("GF_PED",),
            ),
            WorkflowActionSpec(
                key="vpt2_vci_run",
                label="Run VPT2/VCI",
                command="vpt2-vci",
                required_sections=("QFF",),
                produced_sections=("VPT2_VCI",),
            ),
        ),
    ),
    WindowSpec(
        key="electronic_spectroscopy",
        title="Electronic Spectroscopy",
        description="Collect electronic-state data, visualize orbitals and prepare spectra.",
        category="spectroscopy",
        capabilities=(
            "collect electronic-state and transition data from QM adapters",
            "visualize orbitals and densities with external viewers",
            "draw UV/visible or electronic stick/broadened spectra",
            "export publication-ready electronic spectra",
        ),
        publication_exports=("CSV transition table", "SVG spectrum", "PDF spectrum"),
        external_viewers=("Avogadro", "Avogadro2", "Molden"),
        actions=(
            WorkflowActionSpec(
                key="promote_fchk",
                label="Promote FCHK data",
                command="gaussian promote-fchk",
                produced_sections=("CARTESIAN_HESSIAN", "NORMAL_MODES", "QFF"),
            ),
        ),
    ),
    WindowSpec(
        key="thermochemistry_kinetics",
        title="Thermochemistry / Kinetics",
        description="Collect thermal functions, kinetic models and publication tables.",
        category="thermochemistry",
        required_sections=("BASIC", "ROTATIONAL"),
        produced_sections=("THERMO", "KINETICS"),
        capabilities=(
            "run thermochemistry from normalized ORACLE sections",
            "combine vibrational and rovibrational density of states",
            "prepare kinetic-model inputs and compare rates",
            "export publication thermochemistry and kinetics tables",
        ),
        publication_exports=("CSV thermo table", "LaTeX thermo table", "SVG plots", "PDF plots"),
        actions=(
            WorkflowActionSpec(
                key="thermo_run",
                label="Run Thermo",
                command="thermo",
                required_sections=("BASIC", "ROTATIONAL"),
                produced_sections=("THERMO",),
            ),
            WorkflowActionSpec(
                key="rovib_dos",
                label="Build rovibrational DOS",
                command="rovib dos-rovib",
                required_sections=("ROTATIONAL", "VIBRATIONAL"),
            ),
        ),
    ),
)


WINDOWS_BY_KEY: dict[str, WindowSpec] = {spec.key: spec for spec in ORACLE_GUI_WINDOWS}


def window_spec(key: str) -> WindowSpec:
    try:
        return WINDOWS_BY_KEY[key]
    except KeyError as exc:
        raise KeyError(f"unknown ORACLE GUI window: {key}") from exc


def all_known_sections() -> tuple[str, ...]:
    sections: set[str] = set()
    for spec in ORACLE_GUI_WINDOWS:
        sections.update(spec.required_sections)
        sections.update(spec.produced_sections)
        for action in spec.actions:
            sections.update(action.required_sections)
            sections.update(action.produced_sections)
    return tuple(sorted(sections))
