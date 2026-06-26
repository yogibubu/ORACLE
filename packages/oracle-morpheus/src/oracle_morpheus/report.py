from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path

import numpy as np

from merlino_core.numerics import rank_condition

from .contracts import IsotopologueObservation, ParameterClassConstraint, SemiexperimentalFitRequest
from .fit import (
    SemiexperimentalFitResult,
    _active_mask,
    _atomic_number,
    _build_measurement_model,
    _gic_model,
    _gic_fixed_patterns,
    _gicforge_a1_mask,
    _jacobian_constants_wrt_gics,
    _combined_fixed_parameters,
    _fixed_primitives_from_patterns,
    _hydrogen_fixed_primitives,
    _merge_primitives,
    _parameter_class_transform,
    _primitive_constrained_transform,
    _semiexp_warning_rows,
    _symmetry_expanded_fixed_primitives,
    fit_semiexperimental_geometry,
)
from .geometry_input import read_geometry_input


@dataclass(frozen=True)
class SemiexperimentalGICPreview:
    atoms: tuple[str, ...]
    gic_labels: tuple[str, ...]
    rows: tuple["SemiexperimentalGICPreviewRow", ...]
    suggested_classes: tuple[ParameterClassConstraint, ...]
    warnings: tuple[str, ...]

    @property
    def text(self) -> str:
        lines = [
            "ORACLE semiexperimental GIC preview",
            f"Atoms: {len(self.atoms)}",
            f"Non-redundant GICs: {len(self.gic_labels)}",
            "",
            "Suggested parameter classes:",
        ]
        lines.extend(
            f"  {item.name}:{item.mode}:{'|'.join(item.patterns)}" for item in self.suggested_classes
        )
        if not self.suggested_classes:
            lines.append("  none")
        lines.extend(["", "GIC labels:"])
        lines.extend(f"  {row.label} [{row.kind}] class={row.suggested_class or '-'}" for row in self.rows)
        if self.warnings:
            lines.extend(["", "Warnings:", *[f"  {item}" for item in self.warnings]])
        return "\n".join(lines)


@dataclass(frozen=True)
class SemiexperimentalGICPreviewRow:
    label: str
    kind: str
    atoms: tuple[int, ...]
    suggested_class: str
    state: str


@dataclass(frozen=True)
class SemiexperimentalValidationIssue:
    severity: str
    message: str


@dataclass(frozen=True)
class SemiexperimentalConditioningPreview:
    rank: int
    condition_number: float
    n_observations: int
    n_effective_parameters: int
    components: tuple[str, ...]
    warnings: tuple[str, ...]

    @property
    def text(self) -> str:
        lines = [
            "Semiexperimental conditioning preview",
            f"observations: {self.n_observations}",
            f"effective parameters: {self.n_effective_parameters}",
            f"rank: {self.rank}",
            f"condition number: {self.condition_number:.8g}",
            f"components: {','.join(self.components)}",
        ]
        if self.warnings:
            lines.extend(["Warnings:", *[f"  {item}" for item in self.warnings]])
        return "\n".join(lines)


@dataclass(frozen=True)
class SemiexperimentalBenchmarkCase:
    label: str
    request: SemiexperimentalFitRequest
    max_iter: int | None = None
    step: float = 1.0e-4
    damping: float = 1.0e-8
    max_step: float = 0.25
    prune_condition: float = 0.0


@dataclass(frozen=True)
class SemiexperimentalBenchmarkRow:
    label: str
    rms_MHz: float
    rotational_rms_MHz: float
    iterations: int
    rank: int
    incremental_rank: int
    condition_number: float
    stationary_point: str
    n_parameters: int
    n_kraitchman: int
    gicforge_calls: int
    coordinate_model_reuse_steps: int
    b_projector_secant_updates: int


def preview_semiexperimental_gics(
    xyz: Path,
    observations: tuple[IsotopologueObservation, ...] = (),
) -> SemiexperimentalGICPreview:
    geometry_input = read_geometry_input(Path(xyz))
    atoms = tuple(geometry_input.atoms)
    coords = geometry_input.coordinates_angstrom
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    _prims, _u_matrix, labels = _gic_model(np.asarray(coords, dtype=float), z_numbers)
    suggestions = suggest_parameter_classes(tuple(atoms), labels, observations)
    rows = _preview_rows(labels, suggestions, geometry_input.fixed_parameters)
    warnings = _preview_warnings(labels, suggestions)
    return SemiexperimentalGICPreview(tuple(atoms), labels, rows, suggestions, warnings)


def validate_semiexperimental_request(request: SemiexperimentalFitRequest) -> tuple[SemiexperimentalValidationIssue, ...]:
    issues: list[SemiexperimentalValidationIssue] = []
    try:
        geometry_input = read_geometry_input(Path(request.initial_geometry))
        atoms = geometry_input.atoms
    except Exception as exc:
        return (SemiexperimentalValidationIssue("error", f"Cannot read parent geometry: {exc}"),)
    labels = [obs.label for obs in request.observations]
    if len(labels) != len(set(labels)):
        issues.append(SemiexperimentalValidationIssue("error", "Duplicate isotopologue labels"))
    for obs in request.observations:
        if any(value <= 0.0 for value in obs.constants.as_tuple()):
            issues.append(SemiexperimentalValidationIssue("error", f"{obs.label}: rotational constants must be positive"))
        if obs.weights is not None and any(value <= 0.0 for value in obs.weights.as_tuple()):
            issues.append(SemiexperimentalValidationIssue("error", f"{obs.label}: sigma-derived weights must be positive"))
        seen_atoms = set()
        for atom_index, mass in obs.substitutions.items():
            if atom_index in seen_atoms:
                issues.append(SemiexperimentalValidationIssue("error", f"{obs.label}: duplicate substitution at atom {atom_index}"))
            seen_atoms.add(atom_index)
            if atom_index < 1 or atom_index > len(atoms):
                issues.append(SemiexperimentalValidationIssue("error", f"{obs.label}: substitution atom {atom_index} is out of range"))
            elif mass == 2 and atoms[atom_index - 1].upper() != "H":
                issues.append(SemiexperimentalValidationIssue("warning", f"{obs.label}: deuterium substitution on non-H atom {atom_index}"))
        if any(abs(value) > 0.25 * max(abs(base), 1.0) for value, base in zip(obs.correction.as_tuple(), obs.constants.as_tuple())):
            issues.append(SemiexperimentalValidationIssue("warning", f"{obs.label}: unusually large vibrational correction"))
    try:
        preview = preview_semiexperimental_gics(request.initial_geometry, request.observations)
    except Exception as exc:
        issues.append(SemiexperimentalValidationIssue("error", f"Cannot generate GIC preview: {exc}"))
        return tuple(issues)
    for parameter_class in request.parameter_classes:
        matches = [label for label in preview.gic_labels if any(pattern.lower() in label.lower() for pattern in parameter_class.patterns)]
        if not matches:
            issues.append(SemiexperimentalValidationIssue("error", f"Parameter class {parameter_class.name} matches no GIC"))
        kinds = {_gic_kind(label) for label in matches}
        if len(kinds) > 1:
            issues.append(SemiexperimentalValidationIssue("error", f"Parameter class {parameter_class.name} mixes coordinate types: {', '.join(sorted(kinds))}"))
    return tuple(issues)


def preview_semiexperimental_conditioning(
    request: SemiexperimentalFitRequest,
    *,
    step: float = 1.0e-4,
) -> SemiexperimentalConditioningPreview:
    geometry_input = read_geometry_input(Path(request.initial_geometry))
    atoms = geometry_input.atoms
    coords_arr = np.asarray(geometry_input.coordinates_angstrom, dtype=float)
    z_numbers = np.array([_atomic_number(symbol) for symbol in atoms], dtype=int)
    prims, u_matrix, labels = _gic_model(coords_arr, z_numbers)
    measurement = _build_measurement_model(request, atoms, coords_arr, prims, u_matrix, labels)
    fixed_parameters = _combined_fixed_parameters(request.fixed_parameters, geometry_input.fixed_parameters)
    fixed_gic_patterns = _gic_fixed_patterns(fixed_parameters)
    fixed_primitives = _merge_primitives(
        _fixed_primitives_from_patterns(fixed_parameters),
        _hydrogen_fixed_primitives(atoms, prims, fixed_parameters, coords=coords_arr),
    )
    fixed_primitives = _symmetry_expanded_fixed_primitives(atoms, coords_arr, prims, fixed_primitives)
    active = _active_mask(labels, fixed_gic_patterns, request.parameter_classes) & _gicforge_a1_mask(labels)
    jac_gic = _jacobian_constants_wrt_gics(atoms, coords_arr, request, prims, u_matrix, active, labels, measurement, step=step)
    transform, _names, _class_by_gic = _parameter_class_transform(labels, active, request.parameter_classes)
    transform, _names = _primitive_constrained_transform(
        coords_arr,
        prims,
        u_matrix,
        active,
        transform,
        _names,
        fixed_primitives,
    )
    jac = jac_gic @ transform
    weighted = jac * np.sqrt(measurement.weights)[:, None]
    conditioning = rank_condition(weighted)
    warnings = []
    if conditioning.rank < weighted.shape[1]:
        warnings.append("rank deficient for the current isotopologues/classes")
    if not np.isfinite(conditioning.condition_number) or conditioning.condition_number > 1.0e8:
        warnings.append("ill-conditioned parameter set")
    return SemiexperimentalConditioningPreview(
        conditioning.rank,
        conditioning.condition_number,
        int(weighted.shape[0]),
        int(weighted.shape[1]),
        measurement.components,
        tuple(warnings),
    )


def suggest_parameter_classes(
    atoms: tuple[str, ...],
    gic_labels: tuple[str, ...],
    observations: tuple[IsotopologueObservation, ...] = (),
) -> tuple[ParameterClassConstraint, ...]:
    suggestions: list[ParameterClassConstraint] = []
    h_substituted = _substituted_hydrogens(atoms, observations)
    h_atoms = tuple(idx + 1 for idx, atom in enumerate(atoms) if atom.upper() == "H")
    if h_atoms:
        for heavy in sorted({heavy for heavy, h in _heavy_h_bonds(atoms, gic_labels)}):
            patterns = tuple(
                f"R({heavy},{h})"
                for h in h_atoms
                if _matches_any(gic_labels, f"R({heavy},{h})") or _matches_any(gic_labels, f"R({h},{heavy})")
            )
            if len(patterns) >= 2:
                name = f"{atoms[heavy - 1].upper()}H_stretches"
                suggestions.append(ParameterClassConstraint(name, patterns, "shared"))
        angle_patterns = tuple(_first_gic_expression(label, "angle") for label in gic_labels if _gic_kind(label) == "angle" and _angle_has_h(label, atoms))
        angle_patterns = tuple(pattern for pattern in angle_patterns if pattern)
        if angle_patterns and not h_substituted:
            suggestions.append(ParameterClassConstraint("XH_angles", angle_patterns, "fixed"))
    return tuple(suggestions)


def write_semiexperimental_html_report(
    path: Path,
    result: SemiexperimentalFitResult,
    request: SemiexperimentalFitRequest,
) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    html = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>ORACLE semiexperimental geometry report</title>",
        "<style>body{font-family:Helvetica,Arial,sans-serif;margin:32px;line-height:1.35}"
        "table{border-collapse:collapse;margin:16px 0;width:100%}"
        "th,td{border:1px solid #ccc;padding:5px 7px;text-align:left;font-size:13px}"
        "th{background:#f3f3f3} code{background:#f7f7f7;padding:1px 3px}</style>",
        "</head><body>",
        "<h1>ORACLE Semiexperimental Geometry Report</h1>",
        "<h2>Diagnostics</h2>",
        "<table><tr><th>Quantity</th><th>Value</th></tr>",
        _row("RMS", f"{result.rms_MHz:.8g}"),
        _row("Iterations", str(result.iterations)),
        _row("Stationary point", result.stationary_point),
        _row("Convergence", result.diagnostics.convergence_reason),
        _row("Linear solver", result.diagnostics.linear_solver),
        _row("Robust loss", result.diagnostics.robust_loss),
        _row("Robust scale", f"{result.diagnostics.robust_scale:.8g}"),
        _row("Downweighted rows", str(result.diagnostics.robust_downweighted_observations)),
        _row("Downweighted isotopologues", str(result.diagnostics.robust_downweighted_isotopologues)),
        _row("Rank", str(result.diagnostics.rank)),
        _row("Condition number", f"{result.diagnostics.condition_number:.8g}"),
        _row("Observable", result.diagnostics.observable),
        _row("Components", ",".join(result.diagnostics.components)),
        "</table>",
        "<h2>Warnings</h2>",
        _diagnostic_warnings_table(result),
        "<h2>Parameter Classes</h2>",
        _classes_table(request.parameter_classes),
        "<h2>GIC Parameters</h2>",
        _parameters_table(result),
        "<h2>Final Cartesian Geometry Parameters</h2>",
        _geometry_parameters_table(result),
        "<h2>Rotational Constants</h2>",
        _rotational_constants_table(result),
        "<h2>Residuals</h2>",
        _residuals_table(result),
        "<h2>Kraitchman Comparison</h2>",
        _kraitchman_table(result),
        "</body></html>",
    ]
    target.write_text("\n".join(html) + "\n", encoding="utf-8")
    return target


def run_semiexperimental_benchmark(
    cases: tuple[SemiexperimentalBenchmarkCase, ...],
    *,
    outdir: Path | None = None,
    max_iter: int | None = None,
) -> tuple[SemiexperimentalBenchmarkRow, ...]:
    rows = []
    for case in cases:
        case_out = Path(outdir) / case.label if outdir is not None else None
        result = fit_semiexperimental_geometry(
            case.request,
            max_iter=max_iter if max_iter is not None else case.max_iter,
            step=case.step,
            damping=case.damping,
            max_step=case.max_step,
            prune_condition=case.prune_condition,
            outdir=case_out,
        )
        rot_diffs = [item.difference_MHz for item in result.rotational_constants]
        rotational_rms = float(np.sqrt(np.mean(np.asarray(rot_diffs, dtype=float) ** 2))) if rot_diffs else 0.0
        rows.append(
            SemiexperimentalBenchmarkRow(
                case.label,
                result.rms_MHz,
                rotational_rms,
                result.iterations,
                result.diagnostics.rank,
                result.diagnostics.incremental_rank,
                result.diagnostics.condition_number,
                result.stationary_point,
                len(result.parameters),
                len(result.kraitchman),
                result.diagnostics.gicforge_calls,
                result.diagnostics.coordinate_model_reuse_steps,
                result.diagnostics.b_projector_secant_updates,
            )
        )
    return tuple(rows)


def benchmark_csv(rows: tuple[SemiexperimentalBenchmarkRow, ...]) -> str:
    lines = [
        "label,rms,rotational_rms_MHz,iterations,rank,incremental_rank,condition_number,"
        "stationary_point,n_parameters,n_kraitchman,gicforge_calls,coordinate_model_reuse_steps,"
        "b_projector_secant_updates"
    ]
    for row in rows:
        lines.append(
            f"{row.label},{row.rms_MHz:.12g},{row.rotational_rms_MHz:.12g},{row.iterations},{row.rank},"
            f"{row.incremental_rank},{row.condition_number:.12g},{row.stationary_point},"
            f"{row.n_parameters},{row.n_kraitchman},{row.gicforge_calls},{row.coordinate_model_reuse_steps},"
            f"{row.b_projector_secant_updates}"
        )
    return "\n".join(lines) + "\n"


def semiexperimental_latex_tables(result: SemiexperimentalFitResult) -> dict[str, str]:
    return {
        "parameters": _latex_parameter_table(result),
        "rotational_constants": _latex_rotational_constants_table(result),
        "residuals": _latex_residual_table(result),
        "kraitchman": _latex_kraitchman_table(result),
    }


def _preview_warnings(
    labels: tuple[str, ...],
    suggestions: tuple[ParameterClassConstraint, ...],
) -> tuple[str, ...]:
    warnings = []
    if not labels:
        warnings.append("No non-redundant GIC labels were generated")
    if not suggestions:
        warnings.append("No automatic parameter-class suggestion was found")
    return tuple(warnings)


def _preview_rows(
    labels: tuple[str, ...],
    suggestions: tuple[ParameterClassConstraint, ...],
    fixed_parameters: tuple[str, ...] = (),
) -> tuple[SemiexperimentalGICPreviewRow, ...]:
    rows = []
    for label in labels:
        assigned = next((item.name for item in suggestions if any(pattern.lower() in label.lower() for pattern in item.patterns)), "")
        state = "fixed_by_input" if any(pattern.lower() in label.lower() for pattern in fixed_parameters) else "active"
        rows.append(SemiexperimentalGICPreviewRow(label, _gic_kind(label), _gic_atoms(label), assigned, state))
    return tuple(rows)


def _gic_kind(label: str) -> str:
    for kind, markers in _GIC_KIND_MARKERS.items():
        if any(marker in label for marker in markers):
            return kind
    if "ring" in label.lower():
        return "ring"
    return "mixed"


def _gic_atoms(label: str) -> tuple[int, ...]:
    atoms = []
    for markers in _GIC_KIND_MARKERS.values():
        for marker in markers:
            start = 0
            while True:
                pos = label.find(marker, start)
                if pos < 0:
                    break
                end = label.find(")", pos)
                if end < 0:
                    break
                atoms.extend(int(part.strip()) for part in label[pos + len(marker):end].split(",") if part.strip().lstrip("-").isdigit() and int(part.strip()) > 0)
                start = end + 1
    return tuple(sorted(set(atoms)))


_GIC_KIND_MARKERS = {
    "bond": ("R(", "B(", "Bond(", "Stretch(", "bond("),
    "angle": ("A(", "Angle(", "Bend(", "angle("),
    "dihedral": ("D(", "Dihedral(", "Torsion(", "dihedral("),
    "out_of_plane": ("U(", "out_of_plane("),
    "linear_bend": ("L(", "Linear(", "LinearBend(", "linear_bend("),
}


def _first_gic_expression(label: str, kind: str) -> str:
    for marker in _GIC_KIND_MARKERS.get(kind, ()):
        start = 0
        while True:
            pos = label.find(marker, start)
            if pos < 0:
                break
            end = label.find(")", pos)
            if end < 0:
                break
            return label[pos:end + 1]
    return ""


def _substituted_hydrogens(atoms: tuple[str, ...], observations: tuple[IsotopologueObservation, ...]) -> set[int]:
    result = set()
    for obs in observations:
        for atom_index in obs.substitutions:
            if 1 <= atom_index <= len(atoms) and atoms[atom_index - 1].upper() == "H":
                result.add(atom_index)
    return result


def _heavy_h_bonds(atoms: tuple[str, ...], labels: tuple[str, ...]) -> tuple[tuple[int, int], ...]:
    pairs = []
    for label in labels:
        for left, right in _label_atom_pairs(label, "bond"):
            if atoms[left - 1].upper() == "H" and atoms[right - 1].upper() != "H":
                pairs.append((right, left))
            elif atoms[right - 1].upper() == "H" and atoms[left - 1].upper() != "H":
                pairs.append((left, right))
    return tuple(pairs)


def _label_atom_pairs(label: str, kind: str) -> tuple[tuple[int, int], ...]:
    pairs = []
    for marker in _GIC_KIND_MARKERS.get(kind, (f"{kind}(",)):
        start = 0
        while True:
            pos = label.find(marker, start)
            if pos < 0:
                break
            end = label.find(")", pos)
            if end < 0:
                break
            parts = [int(part.strip()) for part in label[pos + len(marker):end].split(",") if part.strip().lstrip("-").isdigit()]
            parts = [part for part in parts if part > 0]
            if len(parts) >= 2:
                pairs.append((parts[0], parts[-1]))
            start = end + 1
    return tuple(pairs)


def _matches_any(labels: tuple[str, ...], pattern: str) -> bool:
    low = pattern.lower()
    return any(low in label.lower() for label in labels)


def _angle_has_h(label: str, atoms: tuple[str, ...]) -> bool:
    expr = _first_gic_expression(label, "angle")
    if not expr:
        return False
    atom_text = expr[expr.find("(") + 1 : -1]
    atom_ids = [int(part.strip()) for part in atom_text.split(",") if part.strip().lstrip("-").isdigit()]
    return any(1 <= idx <= len(atoms) and atoms[idx - 1].upper() == "H" for idx in atom_ids)


def _row(key: str, value: str) -> str:
    return f"<tr><th>{escape(key)}</th><td>{escape(value)}</td></tr>"


def _classes_table(classes: tuple[ParameterClassConstraint, ...]) -> str:
    if not classes:
        return "<p>No parameter classes.</p>"
    rows = ["<table><tr><th>Name</th><th>Mode</th><th>Patterns</th></tr>"]
    for item in classes:
        rows.append(
            f"<tr><td>{escape(item.name)}</td><td>{escape(item.mode)}</td>"
            f"<td><code>{escape('|'.join(item.patterns))}</code></td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _diagnostic_warnings_table(result: SemiexperimentalFitResult) -> str:
    warnings = _semiexp_warning_rows(
        result.diagnostics,
        (),
        result.parameters,
        result.geometry_parameters,
        None,
        None,
        None,
    )
    if not warnings:
        return "<p>No diagnostic warnings.</p>"
    rows = ["<table><tr><th>Severity</th><th>Code</th><th>Message</th><th>Context</th></tr>"]
    for item in warnings:
        rows.append(
            f"<tr><td>{escape(item.severity)}</td><td>{escape(item.code)}</td>"
            f"<td>{escape(item.message)}</td><td><code>{escape(item.context)}</code></td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _parameters_table(result: SemiexperimentalFitResult) -> str:
    rows = ["<table><tr><th>Name</th><th>Value</th><th>Sigma</th><th>Active</th><th>Class</th></tr>"]
    for item in result.parameters:
        rows.append(
            f"<tr><td>{escape(item.name)}</td><td>{item.value:.10g}</td><td>{item.sigma:.10g}</td>"
            f"<td>{int(item.active)}</td><td>{escape(item.parameter_class)}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _geometry_parameters_table(result: SemiexperimentalFitResult) -> str:
    if not result.geometry_parameters:
        return "<p>No final topological geometry parameter table available.</p>"
    rows = [
        "<table><tr><th>Kind</th><th>Label</th><th>Atoms</th><th>Symbols</th>"
        "<th>Bond length / Angstrom</th><th>Sigma / Angstrom</th>"
        "<th>Angle or dihedral / degree</th><th>Sigma / degree</th></tr>"
    ]
    for item in result.geometry_parameters:
        rows.append(
            f"<tr><td>{escape(item.kind)}</td><td>{escape(item.label)}</td>"
            f"<td>{'-'.join(str(idx) for idx in item.atom_indices)}</td>"
            f"<td>{escape('-'.join(item.atom_symbols))}</td>"
            f"<td>{'' if item.value_angstrom is None else f'{item.value_angstrom:.8f}'}</td>"
            f"<td>{'' if item.sigma_angstrom is None else f'{item.sigma_angstrom:.8f}'}</td>"
            f"<td>{'' if item.value_degree is None else f'{item.value_degree:.6f}'}</td>"
            f"<td>{'' if item.sigma_degree is None else f'{item.sigma_degree:.6f}'}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _residuals_table(result: SemiexperimentalFitResult) -> str:
    rows = ["<table><tr><th>Isotopologue</th><th>Observable</th><th>Observed</th><th>Calculated</th><th>Residual</th></tr>"]
    for item in result.residuals:
        rows.append(
            f"<tr><td>{escape(item.isotopologue)}</td><td>{escape(item.constant)}</td>"
            f"<td>{item.observed_equilibrium_MHz:.10g}</td><td>{item.calculated_MHz:.10g}</td>"
            f"<td>{item.residual_MHz:.10g}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _rotational_constants_table(result: SemiexperimentalFitResult) -> str:
    if not result.rotational_constants:
        return "<p>No rotational-constant comparison available.</p>"
    rows = [
        "<table><tr><th>Isotopologue</th><th>Component</th>"
        "<th>Corrected experimental / MHz</th><th>Calculated / MHz</th><th>Difference / MHz</th></tr>"
    ]
    for item in result.rotational_constants:
        rows.append(
            f"<tr><td>{escape(item.isotopologue)}</td><td>{escape(item.component)}</td>"
            f"<td>{item.corrected_experimental_MHz:.10g}</td><td>{item.calculated_MHz:.10g}</td>"
            f"<td>{item.difference_MHz:.10g}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _kraitchman_table(result: SemiexperimentalFitResult) -> str:
    if not result.kraitchman:
        return "<p>No single-substitution Kraitchman comparison available.</p>"
    rows = ["<table><tr><th>Isotopologue</th><th>Atom</th><th>Axis</th><th>Kraitchman abs A</th><th>Fit abs A</th><th>Difference A</th></tr>"]
    for item in result.kraitchman:
        rows.append(
            f"<tr><td>{escape(item.isotopologue)}</td><td>{item.atom_index} {escape(item.atom)}</td>"
            f"<td>{escape(item.coordinate)}</td><td>{item.kraitchman_abs_angstrom:.10g}</td>"
            f"<td>{item.fitted_abs_angstrom:.10g}</td><td>{item.difference_angstrom:.10g}</td></tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _latex_parameter_table(result: SemiexperimentalFitResult) -> str:
    lines = [
        "\\begin{tabular}{lrrrl}",
        "\\toprule",
        "Parameter & Value & Sigma & Active & Class \\\\",
        "\\midrule",
    ]
    for item in result.parameters:
        lines.append(f"{_tex(item.name)} & {item.value:.8g} & {item.sigma:.3g} & {int(item.active)} & {_tex(item.parameter_class)} \\\\")
    lines.extend(["\\bottomrule", "\\end{tabular}", ""])
    return "\n".join(lines)


def _latex_residual_table(result: SemiexperimentalFitResult) -> str:
    lines = [
        "\\begin{tabular}{llrrr}",
        "\\toprule",
        "Isotopologue & Observable & Observed & Calculated & Residual \\\\",
        "\\midrule",
    ]
    for item in result.residuals:
        lines.append(f"{_tex(item.isotopologue)} & {_tex(item.constant)} & {item.observed_equilibrium_MHz:.8g} & {item.calculated_MHz:.8g} & {item.residual_MHz:.3g} \\\\")
    lines.extend(["\\bottomrule", "\\end{tabular}", ""])
    return "\n".join(lines)


def _latex_rotational_constants_table(result: SemiexperimentalFitResult) -> str:
    lines = [
        "\\begin{tabular}{llrrr}",
        "\\toprule",
        "Isotopologue & Component & Corrected exp. & Calculated & Difference \\\\",
        "\\midrule",
    ]
    for item in result.rotational_constants:
        lines.append(
            f"{_tex(item.isotopologue)} & {_tex(item.component)} & "
            f"{item.corrected_experimental_MHz:.8g} & {item.calculated_MHz:.8g} & "
            f"{item.difference_MHz:.3g} \\\\"
        )
    lines.extend(["\\bottomrule", "\\end{tabular}", ""])
    return "\n".join(lines)


def _latex_kraitchman_table(result: SemiexperimentalFitResult) -> str:
    lines = [
        "\\begin{tabular}{lllrrr}",
        "\\toprule",
        "Isotopologue & Atom & Axis & Kraitchman & Fit & Difference \\\\",
        "\\midrule",
    ]
    for item in result.kraitchman:
        lines.append(f"{_tex(item.isotopologue)} & {item.atom_index} {_tex(item.atom)} & {_tex(item.coordinate)} & {item.kraitchman_abs_angstrom:.6g} & {item.fitted_abs_angstrom:.6g} & {item.difference_angstrom:.3g} \\\\")
    lines.extend(["\\bottomrule", "\\end{tabular}", ""])
    return "\n".join(lines)


def _tex(text: str) -> str:
    return str(text).replace("\\", "\\textbackslash{}").replace("_", "\\_").replace("&", "\\&").replace("%", "\\%")
