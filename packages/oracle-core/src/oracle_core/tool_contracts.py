from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from .sectioned_xyz import is_section_header_line, read_sectioned_lines, xyz_tail_start


PLANNED_FRAMEWORK_NAME = "MATRIX"
PLANNED_FRAMEWORK_EXPANSION = (
    "Molecular Analysis Toolkit for Reusable Integrated eXperiments"
)


@dataclass(frozen=True)
class ToolContract:
    key: str
    display_name: str
    current_package: str
    standalone_command: str
    required_sections: tuple[str, ...] = ()
    optional_sections: tuple[str, ...] = ()
    produced_sections: tuple[str, ...] = ()
    owned_sections: tuple[str, ...] = ()
    status: str = "implemented"
    planned_name: str = ""
    expanded_name: str = ""
    notes: str = ""

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ToolReadiness:
    contract: ToolContract
    xyzin_path: Path
    present_sections: tuple[str, ...]
    missing_required_sections: tuple[str, ...]

    @property
    def ready(self) -> bool:
        return not self.missing_required_sections

    def to_dict(self) -> dict[str, object]:
        return {
            "key": self.contract.key,
            "display_name": self.contract.display_name,
            "planned_name": self.contract.planned_name,
            "xyzin_path": str(self.xyzin_path),
            "ready": self.ready,
            "required_sections": self.contract.required_sections,
            "present_sections": self.present_sections,
            "missing_required_sections": self.missing_required_sections,
        }


TOOL_CONTRACTS: tuple[ToolContract, ...] = (
    ToolContract(
        key="babel",
        display_name="ORACLE-Babel",
        current_package="oracle-chem/oracle-babel",
        standalone_command="oracle babel preprocess SOURCE OUTPUT",
        produced_sections=("SOURCE", "BASIC", "SYMMETRY", "TOPOLOGY", "SYNTHONS"),
        notes="Normalizes external geometry/QM/SMILES sources into the shared xyzin state.",
    ),
    ToolContract(
        key="fragments",
        display_name="Fragments / nano-lego",
        current_package="oracle-fragments",
        standalone_command="oracle fragments build molecule.xyzin",
        required_sections=("TOPOLOGY", "SYNTHONS"),
        produced_sections=("FRAGMENTS",),
        owned_sections=("FRAGMENTS",),
        status="implemented-skeleton",
        notes="Fragmentation and future assembly are clients of saved topology/synthons.",
    ),
    ToolContract(
        key="gicforge",
        display_name="GICForge",
        current_package="oracle-gicforge",
        standalone_command="oracle gicforge build molecule.xyzin",
        required_sections=("VALIDATION", "SYMMETRY", "TOPOLOGY", "SYNTHONS"),
        optional_sections=("FRAGMENTS",),
        produced_sections=("GIC", "SYCART"),
        owned_sections=("GIC", "SYCART"),
        planned_name="NEO",
        expanded_name="Nonredundant Equivariant Orthogonalizer",
        notes="Runtime name stays GICForge until the refactory is stable; final tool name should be NEO.",
    ),
    ToolContract(
        key="gf",
        display_name="GF / PED",
        current_package="oracle-gf",
        standalone_command="oracle gf --xyzin molecule.xyzin",
        required_sections=("GIC", "CARTESIAN_HESSIAN"),
        optional_sections=("SYNTHONS", "GF_PED"),
        produced_sections=("GF_PED",),
        owned_sections=("GF_PED",),
        notes="May optionally import an FCHK through the Gaussian adapter before consuming xyzin state.",
    ),
    ToolContract(
        key="morpheus",
        display_name="SEFit / MORPHEUS",
        current_package="oracle-morpheus",
        standalone_command="oracle semiexp --xyzin molecule.xyzin --job job.toml --outdir run",
        required_sections=("ISOTOPOLOGUES",),
        optional_sections=("GIC", "SYCART", "MORPHEUS"),
        produced_sections=("MORPHEUS",),
        owned_sections=("MORPHEUS", "ISOTOPOLOGUES"),
        notes="Consumes frozen coordinate models or symmetry-Cartesian state; owns semiexperimental fit state.",
    ),
    ToolContract(
        key="trinity",
        display_name="TRINITY",
        current_package="oracle-trinity",
        standalone_command="oracle trinity prepare molecule.xyzin --run-dir run --engine-command CMD",
        required_sections=("BASIC",),
        optional_sections=("GIC", "SYCART", "TRINITY"),
        produced_sections=("TRINITY",),
        owned_sections=("TRINITY",),
        status="prepare-only",
        notes="External energy/gradient optimizer branch; runner loop is intentionally not implemented yet.",
    ),
    ToolContract(
        key="rovib",
        display_name="Rovib utilities",
        current_package="oracle-rovib",
        standalone_command="oracle rovib summarize molecule.xyzin",
        required_sections=("ROTATIONAL",),
        optional_sections=("VIBRATIONAL", "DELTABVIB", "CORIOLIS", "QCENT"),
        produced_sections=("CORIOLIS", "QCENT"),
        owned_sections=("ROTATIONAL", "VIBRATIONAL", "DELTABVIB", "CORIOLIS", "QCENT"),
    ),
    ToolContract(
        key="thermo",
        display_name="Thermo",
        current_package="oracle-thermo",
        standalone_command="oracle thermo molecule.xyzin",
        required_sections=("BASIC", "ROTATIONAL"),
        optional_sections=("VIBRATIONAL",),
        produced_sections=("THERMO",),
        owned_sections=("THERMO",),
    ),
    ToolContract(
        key="vpt2_vci",
        display_name="VPT2 / VCI",
        current_package="oracle-vpt2-vci",
        standalone_command="oracle vpt2-vci --xyzin molecule.xyzin --run-dir run",
        required_sections=("QFF",),
        optional_sections=("VPT2_VCI",),
        produced_sections=("VPT2_VCI",),
        owned_sections=("VPT2_VCI",),
        notes="FCHK and indexed QFF text are adapter entry points; normalized standalone input is #QFF.",
    ),
    ToolContract(
        key="dvr",
        display_name="DVR",
        current_package="oracle-dvr",
        standalone_command="oracle dvr run --xyzin molecule.xyzin",
        required_sections=("DVR",),
        produced_sections=("DVR",),
        owned_sections=("DVR",),
        notes="Gaussian scan logs are prepare-time adapter inputs; post-run state is collected into #DVR.",
    ),
    ToolContract(
        key="gui",
        display_name="Desktop GUI",
        current_package="oracle-gui",
        standalone_command="python -m oracle_gui molecule.xyzin",
        optional_sections=("BASIC", "GIC", "GF_PED", "MORPHEUS", "TRINITY", "VPT2_VCI", "DVR"),
        status="orchestrator",
        planned_name="ORACLE",
        expanded_name="Operator for Routing, Analysis, Control, Launch and Exploration",
        notes="The GUI remains the user-facing ORACLE application; it must not own scientific parsers.",
    ),
)


def tool_contracts(*, include_gui: bool = True) -> tuple[ToolContract, ...]:
    if include_gui:
        return TOOL_CONTRACTS
    return tuple(contract for contract in TOOL_CONTRACTS if contract.key != "gui")


def tool_contract(key: str) -> ToolContract:
    normalized = key.strip().lower().replace("-", "_")
    for contract in TOOL_CONTRACTS:
        if contract.key == normalized:
            return contract
        if contract.display_name.lower() == normalized:
            return contract
        if contract.planned_name and contract.planned_name.lower() == normalized:
            return contract
    raise KeyError(f"unknown ORACLE tool contract: {key}")


def tool_contract_lines(
    contracts: tuple[ToolContract, ...] | None = None,
    *,
    include_notes: bool = True,
) -> list[str]:
    rows = tool_contracts() if contracts is None else contracts
    lines: list[str] = []
    for contract in rows:
        lines.append(f"{contract.key}: {contract.display_name}")
        if contract.planned_name:
            lines.append(f"  planned_name: {contract.planned_name}")
        if contract.expanded_name:
            lines.append(f"  expanded_name: {contract.expanded_name}")
        lines.append(f"  package: {contract.current_package}")
        lines.append(f"  command: {contract.standalone_command}")
        lines.append(f"  required: {_join_sections(contract.required_sections)}")
        lines.append(f"  optional: {_join_sections(contract.optional_sections)}")
        lines.append(f"  produced: {_join_sections(contract.produced_sections)}")
        lines.append(f"  owned: {_join_sections(contract.owned_sections)}")
        lines.append(f"  status: {contract.status}")
        if include_notes and contract.notes:
            lines.append(f"  notes: {contract.notes}")
    return lines


def tool_contract_markdown_table(contracts: tuple[ToolContract, ...] | None = None) -> str:
    rows = tool_contracts() if contracts is None else contracts
    lines = [
        "| Key | Current name | Planned name | Expanded name | Package | "
        "Required sections | Produced sections | Status |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for contract in rows:
        lines.append(
            "| "
            + " | ".join(
                (
                    contract.key,
                    contract.display_name,
                    contract.planned_name or "",
                    contract.expanded_name or "",
                    contract.current_package,
                    ", ".join(contract.required_sections) or "none",
                    ", ".join(contract.produced_sections) or "none",
                    contract.status,
                )
            )
            + " |"
        )
    return "\n".join(lines)


def tool_contracts_json(contracts: tuple[ToolContract, ...] | None = None) -> str:
    rows = tool_contracts() if contracts is None else contracts
    return json.dumps([contract.to_dict() for contract in rows], indent=2, sort_keys=True)


def xyzin_section_names(path: Path | str) -> tuple[str, ...]:
    lines = read_sectioned_lines(Path(path))
    names: list[str] = []
    for raw in lines[xyz_tail_start(lines) :]:
        if is_section_header_line(raw):
            names.append(raw.strip()[1:].strip().upper())
    return tuple(names)


def tool_contract_readiness(path: Path | str, contract: ToolContract | str) -> ToolReadiness:
    target = Path(path)
    resolved = tool_contract(contract) if isinstance(contract, str) else contract
    present = xyzin_section_names(target)
    present_set = set(present)
    missing = tuple(section for section in resolved.required_sections if section not in present_set)
    return ToolReadiness(
        contract=resolved,
        xyzin_path=target,
        present_sections=present,
        missing_required_sections=missing,
    )


def tool_contract_readinesses(
    path: Path | str,
    contracts: tuple[ToolContract, ...] | None = None,
) -> tuple[ToolReadiness, ...]:
    rows = tool_contracts() if contracts is None else contracts
    return tuple(tool_contract_readiness(path, contract) for contract in rows)


def tool_readiness_lines(readinesses: tuple[ToolReadiness, ...]) -> list[str]:
    lines: list[str] = []
    for readiness in readinesses:
        missing = _join_sections(readiness.missing_required_sections)
        lines.append(
            f"{readiness.contract.key}: ready={int(readiness.ready)} "
            f"missing={missing}"
        )
    return lines


def tool_readiness_markdown_table(readinesses: tuple[ToolReadiness, ...]) -> str:
    lines = [
        "| Key | Current name | Planned name | Ready | Missing required sections |",
        "| --- | --- | --- | --- | --- |",
    ]
    for readiness in readinesses:
        lines.append(
            "| "
            + " | ".join(
                (
                    readiness.contract.key,
                    readiness.contract.display_name,
                    readiness.contract.planned_name or "",
                    str(int(readiness.ready)),
                    ", ".join(readiness.missing_required_sections) or "none",
                )
            )
            + " |"
        )
    return "\n".join(lines)


def tool_readiness_json(readinesses: tuple[ToolReadiness, ...]) -> str:
    return json.dumps(
        [readiness.to_dict() for readiness in readinesses],
        indent=2,
        sort_keys=True,
    )


def _join_sections(sections: tuple[str, ...]) -> str:
    return ", ".join(sections) if sections else "none"
