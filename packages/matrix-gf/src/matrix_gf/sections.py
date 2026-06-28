from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from matrix_core import (
    key_value_section_lines,
    parse_key_value_section,
    read_sectioned_lines,
    replace_section,
    section_content,
)


ORACLE_XYZ_GF_PED_SCHEMA = "oracle.xyz.gf_ped.v1"


@dataclass(frozen=True)
class GFModeRow:
    index: int
    frequency_cm: float
    block: str = ""


@dataclass(frozen=True)
class GFGICRow:
    identifier: str
    name: str
    irrep: str
    label: str
    ped: tuple[float, ...] = ()
    scaling_factor: float | None = None


@dataclass(frozen=True)
class GFLargeAmplitudeCoordinateRow:
    identifier: str
    name: str
    irrep: str
    family: str
    label: str


@dataclass(frozen=True)
class GFLargeAmplitudeBlockRow:
    label: str
    family: str
    gics: tuple[str, ...]
    frequencies_cm: tuple[float, ...]
    max_f_coupling_to_rest: float = 0.0
    relative_f_coupling_to_rest: float = 0.0
    max_g_coupling_to_rest: float = 0.0
    relative_g_coupling_to_rest: float = 0.0


@dataclass(frozen=True)
class GFLargeAmplitudeModeRow:
    index: int
    frequency_cm: float
    ped_percent: float


@dataclass(frozen=True)
class GFPEDSection:
    source_kind: str = "xyzin"
    source_path: Path | None = None
    hessian_source: str = ""
    coordinate_source: str = ""
    report_path: Path | None = None
    csv_dir: Path | None = None
    status: str = "complete"
    point_group: str = "UNKNOWN"
    symmetrized_gics: bool = False
    matrix_model: str = "FULL"
    hessian_correction: str = "NONE"
    force_threshold: float | None = None
    modes: tuple[GFModeRow, ...] = ()
    gics: tuple[GFGICRow, ...] = ()
    large_amplitude_coordinates: tuple[GFLargeAmplitudeCoordinateRow, ...] = ()
    large_amplitude_blocks: tuple[GFLargeAmplitudeBlockRow, ...] = ()
    large_amplitude_modes: tuple[GFLargeAmplitudeModeRow, ...] = ()
    schema: str = ORACLE_XYZ_GF_PED_SCHEMA

    def __post_init__(self) -> None:
        for attr in ("source_path", "report_path", "csv_dir"):
            value = getattr(self, attr)
            if value is not None:
                object.__setattr__(self, attr, Path(value))
        object.__setattr__(self, "modes", tuple(self.modes))
        object.__setattr__(self, "gics", tuple(self.gics))
        object.__setattr__(
            self,
            "large_amplitude_coordinates",
            tuple(self.large_amplitude_coordinates),
        )
        object.__setattr__(self, "large_amplitude_blocks", tuple(self.large_amplitude_blocks))
        object.__setattr__(self, "large_amplitude_modes", tuple(self.large_amplitude_modes))


def gf_ped_section_from_report(
    report,
    *,
    source_kind: str | None = None,
    source_path: Path | str | None = None,
    report_path: Path | str | None = None,
    csv_dir: Path | str | None = None,
    status: str = "complete",
) -> GFPEDSection:
    """Build a normalized #GF_PED section from an ORACLE GF report object."""
    result = report.result
    frequencies = np.asarray(result.frequencies_cm, dtype=float).reshape(-1)
    block_labels = tuple(getattr(result, "block_labels", ()))
    modes = tuple(
        GFModeRow(
            index=idx,
            frequency_cm=float(frequency),
            block=block_labels[idx - 1] if idx - 1 < len(block_labels) else "",
        )
        for idx, frequency in enumerate(frequencies, start=1)
    )

    ped_values = np.asarray(result.ped.values, dtype=float)
    gic_labels = tuple(getattr(result, "gic_labels", ()))
    gic_names = tuple(getattr(result, "gic_names", ()))
    gic_irreps = tuple(getattr(result, "gic_irreps", ()))
    scaling = getattr(result, "scaling_factors", None)
    scaling_values = None if scaling is None else np.asarray(scaling, dtype=float).reshape(-1)
    gics: list[GFGICRow] = []
    for index, label in enumerate(gic_labels):
        ped = (
            tuple(float(value) for value in ped_values[index, :])
            if ped_values.ndim == 2 and index < ped_values.shape[0]
            else ()
        )
        scale = (
            float(scaling_values[index])
            if scaling_values is not None and index < scaling_values.size
            else None
        )
        gics.append(
            GFGICRow(
                identifier=f"GIC{index + 1:03d}",
                name=gic_names[index] if index < len(gic_names) else f"GIC{index + 1:03d}",
                irrep=gic_irreps[index] if index < len(gic_irreps) else "UNK",
                label=str(label),
                ped=ped,
                scaling_factor=scale,
            )
        )

    large = getattr(result, "large_amplitude", None)
    large_coordinates = ()
    large_blocks = ()
    large_modes = ()
    if large is not None:
        large_coordinates = tuple(
            GFLargeAmplitudeCoordinateRow(
                identifier=f"GIC{coordinate.index:03d}",
                name=coordinate.name,
                irrep=coordinate.irrep,
                family=coordinate.family,
                label=coordinate.label,
            )
            for coordinate in large.coordinates
        )
        large_blocks = tuple(
            GFLargeAmplitudeBlockRow(
                label=block.label,
                family=block.family,
                gics=tuple(f"GIC{index:03d}" for index in block.indices),
                frequencies_cm=block.frequencies_cm,
                max_f_coupling_to_rest=block.max_f_coupling_to_rest,
                relative_f_coupling_to_rest=block.relative_f_coupling_to_rest,
                max_g_coupling_to_rest=block.max_g_coupling_to_rest,
                relative_g_coupling_to_rest=block.relative_g_coupling_to_rest,
            )
            for block in large.blocks
        )
        large_modes = tuple(
            GFLargeAmplitudeModeRow(
                index=item.mode,
                frequency_cm=item.frequency_cm,
                ped_percent=item.ped_percent,
            )
            for item in large.mode_contributions
        )

    resolved_source_path = Path(source_path) if source_path is not None else Path(report.fchk_path)
    resolved_source_kind = source_kind or _infer_source_kind(report, resolved_source_path)
    return GFPEDSection(
        source_kind=resolved_source_kind,
        source_path=resolved_source_path,
        hessian_source=getattr(report, "hessian_source", "") or str(resolved_source_path),
        coordinate_source=getattr(result, "coordinate_source", ""),
        report_path=None if report_path is None else Path(report_path),
        csv_dir=None if csv_dir is None else Path(csv_dir),
        status=status,
        point_group=getattr(result, "point_group", "UNKNOWN"),
        symmetrized_gics=bool(getattr(result, "symmetrized_gics", False)),
        matrix_model=getattr(result, "matrix_model", "FULL"),
        hessian_correction=getattr(result, "hessian_correction", "NONE"),
        force_threshold=getattr(result, "force_threshold", None),
        modes=modes,
        gics=tuple(gics),
        large_amplitude_coordinates=large_coordinates,
        large_amplitude_blocks=large_blocks,
        large_amplitude_modes=large_modes,
    )


def gf_ped_section_lines(section: GFPEDSection) -> list[str]:
    values = {
        "STATUS": section.status,
        "SOURCE_KIND": section.source_kind,
        "SOURCE_PATH": section.source_path,
        "HESSIAN_SOURCE": section.hessian_source,
        "COORDINATE_SOURCE": section.coordinate_source,
        "REPORT": section.report_path,
        "CSV_DIR": section.csv_dir,
        "POINT_GROUP": section.point_group,
        "SYMMETRIZED_GICS": int(section.symmetrized_gics),
        "MATRIX_MODEL": section.matrix_model,
        "HESSIAN_CORRECTION": section.hessian_correction,
        "FORCE_THRESHOLD": None
        if section.force_threshold is None
        else _format_float(section.force_threshold),
        "MODE_COUNT": len(section.modes),
        "GIC_COUNT": len(section.gics),
    }
    lines = key_value_section_lines(
        ORACLE_XYZ_GF_PED_SCHEMA,
        values,
        key_order=(
            "STATUS",
            "SOURCE_KIND",
            "SOURCE_PATH",
            "HESSIAN_SOURCE",
            "COORDINATE_SOURCE",
            "REPORT",
            "CSV_DIR",
            "POINT_GROUP",
            "SYMMETRIZED_GICS",
            "MATRIX_MODEL",
            "HESSIAN_CORRECTION",
            "FORCE_THRESHOLD",
            "MODE_COUNT",
            "GIC_COUNT",
        ),
    )
    lines.append("[MODES]")
    if section.modes:
        for mode in section.modes:
            block = f" BLOCK={mode.block}" if mode.block else ""
            lines.append(f"{mode.index} FREQUENCY_CM-1={_format_float(mode.frequency_cm)}{block}")
    else:
        lines.append("NONE")

    lines.append("[GICS]")
    if section.gics:
        for gic in section.gics:
            scale = (
                "" if gic.scaling_factor is None else f" SCALE={_format_float(gic.scaling_factor)}"
            )
            lines.append(
                f"{gic.identifier} NAME={gic.name} IRREP={gic.irrep}{scale} LABEL={gic.label}"
            )
    else:
        lines.append("NONE")

    lines.append("[PED]")
    if section.gics and any(gic.ped for gic in section.gics):
        for gic in section.gics:
            lines.append(
                f"{gic.identifier} VALUES={','.join(_format_float(value) for value in gic.ped)}"
            )
    else:
        lines.append("NONE")
    lines.append("[LARGE_AMPLITUDE_COORDINATES]")
    if section.large_amplitude_coordinates:
        for coordinate in section.large_amplitude_coordinates:
            lines.append(
                f"{coordinate.identifier} NAME={coordinate.name} IRREP={coordinate.irrep} "
                f"FAMILY={coordinate.family} LABEL={coordinate.label}"
            )
    else:
        lines.append("NONE")

    lines.append("[LARGE_AMPLITUDE_BLOCKS]")
    if section.large_amplitude_blocks:
        for block in section.large_amplitude_blocks:
            lines.append(
                f"{block.label} FAMILY={block.family} DIM={len(block.gics)} "
                f"GICS={','.join(block.gics)} "
                f"FREQUENCIES_CM-1={','.join(_format_float(value) for value in block.frequencies_cm)} "
                f"F_COUPLE_MAX={_format_float(block.max_f_coupling_to_rest)} "
                f"F_COUPLE_REL={_format_float(block.relative_f_coupling_to_rest)} "
                f"G_COUPLE_MAX={_format_float(block.max_g_coupling_to_rest)} "
                f"G_COUPLE_REL={_format_float(block.relative_g_coupling_to_rest)}"
            )
    else:
        lines.append("NONE")

    lines.append("[LARGE_AMPLITUDE_MODE_PED]")
    if section.large_amplitude_modes:
        for mode in section.large_amplitude_modes:
            lines.append(
                f"{mode.index} FREQUENCY_CM-1={_format_float(mode.frequency_cm)} "
                f"PED_PERCENT={_format_float(mode.ped_percent)}"
            )
    else:
        lines.append("NONE")
    return lines


def parse_gf_ped_section(lines: Iterable[str]) -> GFPEDSection:
    raw_lines = list(lines)
    values = parse_key_value_section(_header_lines(raw_lines))
    schema = values.get("SCHEMA", ORACLE_XYZ_GF_PED_SCHEMA)
    if schema != ORACLE_XYZ_GF_PED_SCHEMA:
        raise ValueError(f"unsupported GF_PED schema: {schema}")
    modes = tuple(
        _parse_mode_line(line) for line in _subsection(raw_lines, "MODES") if _data_line(line)
    )
    ped_by_gic = {
        identifier: values
        for identifier, values in (
            _parse_ped_line(line) for line in _subsection(raw_lines, "PED") if _data_line(line)
        )
    }
    gics = tuple(
        _parse_gic_line(line, ped_by_gic=ped_by_gic)
        for line in _subsection(raw_lines, "GICS")
        if _data_line(line)
    )
    large_coordinates = tuple(
        _parse_large_amplitude_coordinate_line(line)
        for line in _subsection(raw_lines, "LARGE_AMPLITUDE_COORDINATES")
        if _data_line(line)
    )
    large_blocks = tuple(
        _parse_large_amplitude_block_line(line)
        for line in _subsection(raw_lines, "LARGE_AMPLITUDE_BLOCKS")
        if _data_line(line)
    )
    large_modes = tuple(
        _parse_large_amplitude_mode_line(line)
        for line in _subsection(raw_lines, "LARGE_AMPLITUDE_MODE_PED")
        if _data_line(line)
    )
    return GFPEDSection(
        source_kind=values.get("SOURCE_KIND", "xyzin"),
        source_path=_optional_path(values.get("SOURCE_PATH")),
        hessian_source=values.get("HESSIAN_SOURCE", ""),
        coordinate_source=values.get("COORDINATE_SOURCE", ""),
        report_path=_optional_path(values.get("REPORT")),
        csv_dir=_optional_path(values.get("CSV_DIR")),
        status=values.get("STATUS", "complete"),
        point_group=values.get("POINT_GROUP", "UNKNOWN"),
        symmetrized_gics=_bool_value(values.get("SYMMETRIZED_GICS"), default=False),
        matrix_model=values.get("MATRIX_MODEL", "FULL"),
        hessian_correction=values.get("HESSIAN_CORRECTION", "NONE"),
        force_threshold=_optional_float(values.get("FORCE_THRESHOLD")),
        modes=modes,
        gics=gics,
        large_amplitude_coordinates=large_coordinates,
        large_amplitude_blocks=large_blocks,
        large_amplitude_modes=large_modes,
        schema=schema,
    )


def read_gf_ped_section(path: Path | str) -> GFPEDSection:
    content = section_content(read_sectioned_lines(Path(path)), "GF_PED")
    if not content:
        raise ValueError("missing #GF_PED section")
    return parse_gf_ped_section(content)


def write_gf_ped_section(path: Path | str, section: GFPEDSection) -> None:
    replace_section(Path(path), "GF_PED", gf_ped_section_lines(section))


def write_gf_ped_section_from_report(
    path: Path | str,
    report,
    *,
    source_kind: str | None = None,
    source_path: Path | str | None = None,
    report_path: Path | str | None = None,
    csv_dir: Path | str | None = None,
    status: str = "complete",
) -> GFPEDSection:
    section = gf_ped_section_from_report(
        report,
        source_kind=source_kind,
        source_path=source_path,
        report_path=report_path,
        csv_dir=csv_dir,
        status=status,
    )
    write_gf_ped_section(path, section)
    return section


def _infer_source_kind(report, source_path: Path) -> str:
    xyzin_path = getattr(report, "xyzin_path", None)
    if xyzin_path is not None and Path(xyzin_path) == source_path:
        return "xyzin"
    return "fchk"


def _header_lines(lines: list[str]) -> list[str]:
    out: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            break
        out.append(line)
    return out


def _subsection(section_lines: list[str], name: str) -> list[str]:
    header = f"[{name.upper()}]"
    start = None
    for idx, line in enumerate(section_lines):
        if line.strip().upper() == header:
            start = idx + 1
            break
    if start is None:
        return []
    end = len(section_lines)
    for idx in range(start, len(section_lines)):
        text = section_lines[idx].strip()
        if text.startswith("[") and text.endswith("]"):
            end = idx
            break
    return list(section_lines[start:end])


def _parse_mode_line(line: str) -> GFModeRow:
    parts = line.split()
    if not parts:
        raise ValueError("empty GF_PED mode line")
    fields = _key_values(parts[1:])
    return GFModeRow(
        index=int(parts[0]),
        frequency_cm=float(fields["FREQUENCY_CM-1"]),
        block=fields.get("BLOCK", ""),
    )


def _parse_gic_line(line: str, *, ped_by_gic: Mapping[str, tuple[float, ...]]) -> GFGICRow:
    if " LABEL=" in line:
        prefix, label = line.split(" LABEL=", 1)
    else:
        prefix, label = line, ""
    parts = prefix.split()
    if not parts:
        raise ValueError("empty GF_PED GIC line")
    fields = _key_values(parts[1:])
    identifier = parts[0]
    return GFGICRow(
        identifier=identifier,
        name=fields.get("NAME", identifier),
        irrep=fields.get("IRREP", "UNK"),
        label=label,
        ped=ped_by_gic.get(identifier, ()),
        scaling_factor=_optional_float(fields.get("SCALE")),
    )


def _parse_large_amplitude_coordinate_line(line: str) -> GFLargeAmplitudeCoordinateRow:
    if " LABEL=" in line:
        prefix, label = line.split(" LABEL=", 1)
    else:
        prefix, label = line, ""
    parts = prefix.split()
    if not parts:
        raise ValueError("empty GF_PED large-amplitude coordinate line")
    fields = _key_values(parts[1:])
    identifier = parts[0]
    return GFLargeAmplitudeCoordinateRow(
        identifier=identifier,
        name=fields.get("NAME", identifier),
        irrep=fields.get("IRREP", "UNK"),
        family=fields.get("FAMILY", ""),
        label=label,
    )


def _parse_large_amplitude_block_line(line: str) -> GFLargeAmplitudeBlockRow:
    parts = line.split()
    if not parts:
        raise ValueError("empty GF_PED large-amplitude block line")
    fields = _key_values(parts[1:])
    return GFLargeAmplitudeBlockRow(
        label=parts[0],
        family=fields.get("FAMILY", ""),
        gics=tuple(item for item in fields.get("GICS", "").split(",") if item),
        frequencies_cm=_float_tuple(fields.get("FREQUENCIES_CM-1", "")),
        max_f_coupling_to_rest=_optional_float(fields.get("F_COUPLE_MAX")) or 0.0,
        relative_f_coupling_to_rest=_optional_float(fields.get("F_COUPLE_REL")) or 0.0,
        max_g_coupling_to_rest=_optional_float(fields.get("G_COUPLE_MAX")) or 0.0,
        relative_g_coupling_to_rest=_optional_float(fields.get("G_COUPLE_REL")) or 0.0,
    )


def _parse_large_amplitude_mode_line(line: str) -> GFLargeAmplitudeModeRow:
    parts = line.split()
    if not parts:
        raise ValueError("empty GF_PED large-amplitude mode line")
    fields = _key_values(parts[1:])
    return GFLargeAmplitudeModeRow(
        index=int(parts[0]),
        frequency_cm=float(fields["FREQUENCY_CM-1"]),
        ped_percent=float(fields["PED_PERCENT"]),
    )


def _parse_ped_line(line: str) -> tuple[str, tuple[float, ...]]:
    parts = line.split(maxsplit=1)
    if len(parts) != 2:
        raise ValueError(f"invalid GF_PED PED line: {line}")
    fields = _key_values([parts[1]])
    return parts[0], _float_tuple(fields.get("VALUES", ""))


def _key_values(parts: Iterable[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for part in parts:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        result[key.upper()] = value
    return result


def _float_tuple(text: str) -> tuple[float, ...]:
    if not text or text.upper() == "NONE":
        return ()
    return tuple(float(item) for item in text.split(",") if item)


def _optional_path(raw: str | None) -> Path | None:
    if raw is None or not raw.strip():
        return None
    return Path(raw)


def _optional_float(raw: str | None) -> float | None:
    if raw is None or not str(raw).strip():
        return None
    return float(str(raw).replace("D", "E").replace("d", "e"))


def _bool_value(raw: str | None, *, default: bool) -> bool:
    if raw is None:
        return default
    return str(raw).strip().upper() in {"1", "TRUE", "YES", "Y"}


def _data_line(line: str) -> bool:
    text = line.strip()
    return bool(text and text.upper() != "NONE")


def _format_float(value: float) -> str:
    return f"{float(value):.12g}"
