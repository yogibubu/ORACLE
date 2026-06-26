"""Point-group irrep labels used by frozen GIC definitions."""

from __future__ import annotations

import re


_IRREP_ORDER_BY_POINT_GROUP = {
    "C1": ("A",),
    "CS": ("A'", "A''"),
    "CI": ("Ag", "Au"),
    "C2": ("A", "B"),
    "C2V": ("A1", "A2", "B1", "B2"),
    "D2": ("A", "B1", "B2", "B3"),
    "C2H": ("Ag", "Bg", "Au", "Bu"),
    "D2H": ("Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"),
    "C3V": ("A1", "A2", "E"),
    "C4V": ("A1", "A2", "B1", "B2", "E"),
    "C6V": ("A1", "A2", "B1", "B2", "E1", "E2"),
    "D3H": ("A1'", "A2'", "E'", "A1''", "A2''", "E''"),
    "D4H": (
        "A1g",
        "A2g",
        "B1g",
        "B2g",
        "Eg",
        "A1u",
        "A2u",
        "B1u",
        "B2u",
        "Eu",
    ),
    "TD": ("A1", "A2", "E", "T1", "T2"),
    "OH": ("A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu", "T1u", "T2u"),
    "IH": ("Ag", "T1g", "T2g", "Gg", "Hg", "Au", "T1u", "T2u", "Gu", "Hu"),
}


def normalized_point_group(point_group: str | None) -> str:
    text = (point_group or "C1").strip()
    return text or "C1"


def irrep_sequence(point_group: str | None) -> tuple[str, ...]:
    """Return a deterministic irrep order with the totally symmetric irrep first."""
    group = normalized_point_group(point_group)
    exact = _IRREP_ORDER_BY_POINT_GROUP.get(group.upper())
    if exact:
        return exact

    match = re.fullmatch(r"([CD])(\d+)([VHD]?)", group.upper())
    if not match:
        return ("A",)
    family, n_text, suffix = match.groups()
    n = int(n_text)
    if family == "C":
        if suffix == "V":
            base = ["A1", "A2"]
        elif suffix == "H":
            base = ["Ag", "Au"]
        else:
            base = ["A"]
        if n % 2 == 0:
            base.extend(["B1", "B2"] if suffix == "V" else ["B"])
        base.extend(f"E{idx}" for idx in range(1, max(1, n // 2)))
        return tuple(dict.fromkeys(base))

    base = ["A1", "A2"]
    if n % 2 == 0:
        base.extend(["B1", "B2"])
    base.extend(f"E{idx}" for idx in range(1, max(1, n // 2)))
    if suffix == "H":
        return tuple(
            dict.fromkeys(
                [f"{name}g" for name in base if not name.startswith("E")]
                + [name + "g" for name in base if name.startswith("E")]
                + [f"{name}u" for name in base if not name.startswith("E")]
                + [name + "u" for name in base if name.startswith("E")]
            )
        )
    return tuple(dict.fromkeys(base))


def total_symmetric_irrep(point_group: str | None) -> str:
    return irrep_sequence(point_group)[0]


def non_total_irrep_sequence(point_group: str | None) -> tuple[str, ...]:
    irreps = irrep_sequence(point_group)
    return tuple(irrep for irrep in irreps if irrep != irreps[0])


def is_total_symmetric_irrep(point_group: str | None, irrep: str | None) -> bool:
    return (irrep or "").strip() == total_symmetric_irrep(point_group)


def irrep_name_prefix(irrep: str | None) -> str:
    """Convert an irrep label to a compact coordinate-name prefix."""
    text = (irrep or "X").strip() or "X"
    return text.replace("'", "p").replace('"', "pp")
