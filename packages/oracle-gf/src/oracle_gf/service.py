from __future__ import annotations

from dataclasses import dataclass
import csv
from io import StringIO
from pathlib import Path

import numpy as np

from oracle_gaussian import hessian_input_from_gaussian_fchk
from oracle_gicforge import read_gic_definition_from_xyzin

from .internal import (
    GFLocalOptions,
    InternalGFResult,
    gf_from_hessian_input_and_xyzin,
    gf_from_hessian_input_with_oracle_gics,
    topology_bonds_from_xyzin,
)
from .nonbonded import nonbonded_cartesian_hessian_correction, synthon_charges_from_xyzin


@dataclass(frozen=True)
class GFReport:
    fchk_path: Path
    result: InternalGFResult
    text: str
    xyzin_path: Path | None = None
    scale_path: Path | None = None


def run_gf_report_from_fchk(fchk_path: Path) -> GFReport:
    """Read an FCHK adapter and return a formatted quick GF/PED report."""
    path = Path(fchk_path)
    hessian_input = hessian_input_from_gaussian_fchk(path)
    result = gf_from_hessian_input_with_oracle_gics(hessian_input)
    return GFReport(path, result, format_gf_report(path, result))


def run_xyzin_gf_report_from_fchk(
    fchk_path: Path,
    xyzin_path: Path,
    *,
    scale_path: Path | None = None,
    scale_records: tuple[str, ...] = (),
    local: bool = False,
    force_threshold: float | None = None,
    block_by_irrep: bool = False,
    subtract_electrostatic: bool = False,
    subtract_uff_vdw: bool = False,
    nonbonded_14_scale: float = 0.5,
) -> GFReport:
    """Run the frozen-xyzin GF branch from a Cartesian Hessian FCHK adapter."""
    path = Path(fchk_path)
    xyzin = Path(xyzin_path)
    definition = read_gic_definition_from_xyzin(xyzin)
    hessian_input = hessian_input_from_gaussian_fchk(path)
    correction = None
    correction_label = "NONE"
    if subtract_electrostatic or subtract_uff_vdw:
        topology_bonds = topology_bonds_from_xyzin(xyzin)
        charges = None
        labels: list[str] = []
        if subtract_electrostatic:
            charges, charge_source = synthon_charges_from_xyzin(xyzin, len(hessian_input.atomic_numbers))
            labels.append(f"ELECTROSTATIC({charge_source})")
        if subtract_uff_vdw:
            labels.append("UFF_VDW")
        correction = nonbonded_cartesian_hessian_correction(
            hessian_input.cartesian_coordinates_bohr,
            hessian_input.atomic_numbers,
            topology_bonds,
            charges=charges,
            electrostatic=subtract_electrostatic,
            uff_vdw=subtract_uff_vdw,
            one_four_scale=nonbonded_14_scale,
        )
        correction_label = "+".join(labels) + f"; 1-4 scale={float(nonbonded_14_scale):g}"
    names = tuple(gic.name for gic in definition.gics)
    labels = tuple(_gic_display_label(gic.identifier, gic.gaussian_expression) for gic in definition.gics)
    scaling = pulay_scaling_factors(
        len(definition.gics),
        labels=labels,
        names=names,
        scale_path=scale_path,
        scale_records=scale_records,
    )
    result = gf_from_hessian_input_and_xyzin(
        hessian_input,
        xyzin,
        scaling_factors=scaling,
        local_options=GFLocalOptions(enabled=local),
        force_threshold=force_threshold,
        block_by_irrep=block_by_irrep,
        cartesian_hessian_correction=correction,
        cartesian_hessian_correction_label=correction_label,
    )
    return GFReport(
        path,
        result,
        format_gf_report(path, result, xyzin_path=xyzin, scale_path=scale_path),
        xyzin_path=xyzin,
        scale_path=scale_path,
    )


def format_gf_report(
    fchk_path: Path,
    result: InternalGFResult,
    *,
    xyzin_path: Path | None = None,
    scale_path: Path | None = None,
) -> str:
    lines = [
        "GF/PED from ORACLE non-redundant GICs",
        f"Source FCHK: {Path(fchk_path)}",
        f"Coordinate source: {result.coordinate_source}",
        f"Point group: {result.point_group}",
        f"Symmetrized GICs: {result.symmetrized_gics}",
        f"Matrix model: {result.matrix_model}",
        f"Cartesian Hessian correction: {result.hessian_correction}",
        f"GIC count: {len(result.gic_labels)}",
        "",
        "Frequencies (cm-1):",
    ]
    if xyzin_path is not None:
        lines.insert(2, f"Frozen xyzin: {Path(xyzin_path)}")
    if result.scaling_factors is not None:
        changed = int(np.sum(np.abs(result.scaling_factors - 1.0) > 1.0e-12))
        lines.insert(
            3 if xyzin_path is not None else 2,
            "Pulay Hessian scaling: "
            f"applied ({changed} factors != 1; file={Path(scale_path) if scale_path else 'inline/default'})",
        )
    if result.force_threshold is not None:
        lines.insert(
            3 if xyzin_path is not None else 2,
            f"Force-constant threshold: {result.force_threshold:g}",
        )
    if result.block_labels:
        counts = ", ".join(
            f"{label}:{result.block_labels.count(label)}"
            for label in dict.fromkeys(result.block_labels)
        )
        lines.insert(3 if xyzin_path is not None else 2, f"GF blocks: {counts}")
    for idx, freq in enumerate(result.frequencies_cm, start=1):
        lines.append(f"  mode {idx:3d}: {freq:12.3f}")

    lines.extend(["", "GIC labels:"])
    for idx, label in enumerate(result.gic_labels, start=1):
        name = result.gic_names[idx - 1] if idx <= len(result.gic_names) else f"GIC{idx:03d}"
        irrep = result.gic_irreps[idx - 1] if idx <= len(result.gic_irreps) else "UNK"
        lines.append(f"  GIC{idx:03d}: {name:12s} irrep={irrep:6s} {label}")

    lines.extend(["", "PED (%) rows=GIC cols=modes:"])
    header = "          " + " ".join(f"M{idx:02d}" for idx in range(1, len(result.frequencies_cm) + 1))
    lines.append(header)
    for idx, row in enumerate(result.ped.values, start=1):
        values = " ".join(f"{value:7.2f}" for value in row)
        lines.append(f"  GIC{idx:03d} {values}")
    return "\n".join(lines)


def gf_csv_tables(report: GFReport) -> dict[str, str]:
    """Return CSV tables for GF frequencies, GIC labels, matrices, normal modes and PED."""
    freq_rows = [["mode", "frequency_cm-1"]]
    freq_rows.extend(
        [[idx, f"{freq:.10g}"] for idx, freq in enumerate(report.result.frequencies_cm, start=1)]
    )

    label_rows = [["gic", "name", "irrep", "label"]]
    for idx, label in enumerate(report.result.gic_labels, start=1):
        name = report.result.gic_names[idx - 1] if idx <= len(report.result.gic_names) else f"GIC{idx:03d}"
        irrep = report.result.gic_irreps[idx - 1] if idx <= len(report.result.gic_irreps) else "UNK"
        label_rows.append([f"GIC{idx:03d}", name, irrep, label])

    ped_rows = [
        ["gic", "name", "irrep", *[f"mode_{idx}" for idx in range(1, len(report.result.frequencies_cm) + 1)]]
    ]
    for idx, row in enumerate(report.result.ped.values, start=1):
        name = report.result.gic_names[idx - 1] if idx <= len(report.result.gic_names) else f"GIC{idx:03d}"
        irrep = report.result.gic_irreps[idx - 1] if idx <= len(report.result.gic_irreps) else "UNK"
        ped_rows.append([f"GIC{idx:03d}", name, irrep, *[f"{value:.10g}" for value in row]])

    return {
        "frequencies.csv": _csv_text(freq_rows),
        "gic_labels.csv": _csv_text(label_rows),
        "ped.csv": _csv_text(ped_rows),
        "normal_modes.csv": _csv_text(_gic_mode_table(report.result.modes_internal, "mode", report.result)),
        "force_constants.csv": _csv_text(_square_gic_table(report.result.force_constants, report.result)),
        "g_matrix.csv": _csv_text(_square_gic_table(report.result.g_matrix, report.result)),
    }


def write_csv_tables(report: GFReport, outdir: Path, *, prefix: str = "gf") -> dict[str, Path]:
    target_dir = Path(outdir)
    target_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for name, text in gf_csv_tables(report).items():
        path = target_dir / f"{prefix}_{name}"
        path.write_text(text, encoding="utf-8")
        written[name] = path
    return written


def pulay_scaling_factors(
    n_gics: int,
    *,
    labels: tuple[str, ...],
    names: tuple[str, ...] = (),
    scale_path: Path | None = None,
    scale_records: tuple[str, ...] = (),
) -> np.ndarray | None:
    """Build diagonal Pulay scaling factors from file and inline records."""
    records: list[tuple[str, float]] = []
    if scale_path is not None:
        records.extend(_read_scaling_records(Path(scale_path)))
    for record in scale_records:
        records.append(_parse_scaling_record(record))
    if not records:
        return None
    factors = np.ones(n_gics, dtype=float)
    for selector, factor in records:
        if factor < 0.0:
            raise ValueError("Pulay scaling factors must be non-negative")
        for index in _resolve_scaling_selector(selector, labels=labels, names=names):
            factors[index] = factor
    return factors


def _read_scaling_records(path: Path) -> list[tuple[str, float]]:
    records: list[tuple[str, float]] = []
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.split("#", 1)[0].split("!", 1)[0].strip()
        if not line:
            continue
        lower = line.lower().replace(",", " ").split()
        if lower[:2] in (["selector", "factor"], ["gic", "factor"], ["name", "factor"]):
            continue
        records.append(_parse_scaling_record(line))
    return records


def _parse_scaling_record(record: str) -> tuple[str, float]:
    text = record.strip()
    if "=" in text:
        selector, value = text.split("=", 1)
    elif "," in text:
        selector, value = text.split(",", 1)
    else:
        parts = text.split()
        if len(parts) != 2:
            raise ValueError(f"Invalid Pulay scaling record: {record!r}")
        selector, value = parts
    return selector.strip(), float(value.strip())


def _resolve_scaling_selector(
    selector: str,
    *,
    labels: tuple[str, ...],
    names: tuple[str, ...],
) -> tuple[int, ...]:
    token = selector.strip()
    lower = token.lower()
    n_gics = len(labels)
    if lower in {"*", "all", "default"}:
        return tuple(range(n_gics))
    if lower.startswith("gic") and lower[3:].isdigit():
        index = int(lower[3:]) - 1
        if 0 <= index < n_gics:
            return (index,)
    if token.isdigit():
        index = int(token) - 1
        if 0 <= index < n_gics:
            return (index,)
    exact = [idx for idx, name in enumerate(names) if name == token]
    exact.extend(idx for idx, label in enumerate(labels) if label == token)
    if exact:
        return tuple(dict.fromkeys(exact))
    matches = [
        idx
        for idx, (name, label) in enumerate(zip(_padded_names(names, n_gics), labels))
        if token in name or token in label
    ]
    if len(matches) == 1:
        return (matches[0],)
    if matches:
        raise ValueError(f"Pulay scaling selector {token!r} is ambiguous")
    raise ValueError(f"Pulay scaling selector {token!r} did not match any GIC")


def _padded_names(names: tuple[str, ...], n_gics: int) -> tuple[str, ...]:
    if len(names) >= n_gics:
        return names[:n_gics]
    return names + tuple("" for _ in range(n_gics - len(names)))


def _gic_display_label(identifier: str, expression: str) -> str:
    return f"{identifier} {expression}" if expression and expression != "NONE" else identifier


def _csv_text(rows: list[list[object]]) -> str:
    stream = StringIO()
    writer = csv.writer(stream)
    writer.writerows(rows)
    return stream.getvalue()


def _gic_mode_table(matrix: np.ndarray, label: str, result: InternalGFResult) -> list[list[object]]:
    values = np.asarray(matrix, dtype=float)
    rows: list[list[object]] = [
        ["gic", "name", "irrep", *[f"{label}_{idx}" for idx in range(1, values.shape[1] + 1)]]
    ]
    for idx, row in enumerate(values, start=1):
        name = result.gic_names[idx - 1] if idx <= len(result.gic_names) else f"GIC{idx:03d}"
        irrep = result.gic_irreps[idx - 1] if idx <= len(result.gic_irreps) else "UNK"
        rows.append([f"GIC{idx:03d}", name, irrep, *[f"{value:.10g}" for value in row]])
    return rows


def _square_gic_table(matrix: np.ndarray, result: InternalGFResult) -> list[list[object]]:
    values = np.asarray(matrix, dtype=float)
    rows: list[list[object]] = [
        ["gic", "name", "irrep", *[f"GIC{idx:03d}" for idx in range(1, values.shape[1] + 1)]]
    ]
    for idx, row in enumerate(values, start=1):
        name = result.gic_names[idx - 1] if idx <= len(result.gic_names) else f"GIC{idx:03d}"
        irrep = result.gic_irreps[idx - 1] if idx <= len(result.gic_irreps) else "UNK"
        rows.append([f"GIC{idx:03d}", name, irrep, *[f"{value:.10g}" for value in row]])
    return rows
