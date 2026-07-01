from __future__ import annotations

import re
from dataclasses import dataclass

from .contracts import ParameterClassConstraint


@dataclass(frozen=True)
class PrimitiveClassSpec:
    """User-facing chemical class defined by primitive-coordinate patterns."""

    name: str
    patterns: tuple[str, ...]

    def validate(self) -> None:
        if not self.name.strip():
            raise ValueError("Primitive class name cannot be empty")
        if not self.patterns or any(not pattern.strip() for pattern in self.patterns):
            raise ValueError("Primitive class patterns cannot be empty")


@dataclass(frozen=True)
class DerivedPrimitiveClassPlan:
    """Disjoint GIC classes and blocked complement derived from primitive classes."""

    fixed_patterns: tuple[str, ...]
    parameter_classes: tuple[ParameterClassConstraint, ...]
    rejected_labels: tuple[str, ...]
    budget: int


def parse_primitive_class_spec(raw: str) -> PrimitiveClassSpec:
    """Parse `name:primitive[|primitive...]` into a primitive class specification."""

    parts = str(raw).split(":", 1)
    if len(parts) != 2:
        raise ValueError("--primitive-class must be name:primitive[|primitive...]")
    spec = PrimitiveClassSpec(
        name=parts[0].strip(),
        patterns=tuple(part.strip() for part in parts[1].split("|") if part.strip()),
    )
    spec.validate()
    return spec


def derive_primitive_class_plan(
    gic_labels: tuple[str, ...],
    primitive_classes: tuple[PrimitiveClassSpec, ...],
    *,
    min_fraction: float = 0.70,
    cross_fraction_max: float = 0.20,
    max_classes: int | None = None,
) -> DerivedPrimitiveClassPlan:
    """Map primitive-class definitions onto a disjoint reduced GIC class model.

    Classes are interpreted as a priority cascade: broader classes should be
    declared first and more specific classes later. A GIC is assigned to the
    last class whose primitive coefficient reaches `min_fraction`; earlier
    classes are used only when all later classes stay below
    `cross_fraction_max`. Unsupported GICs are frozen. If `max_classes` is set,
    only the best-supported classes are kept.
    """

    if not primitive_classes:
        return DerivedPrimitiveClassPlan((), (), (), 0)
    if min_fraction <= 0.0:
        raise ValueError("min_fraction must be positive")
    if cross_fraction_max < 0.0:
        raise ValueError("cross_fraction_max must be non-negative")
    for spec in primitive_classes:
        spec.validate()

    class_patterns = {
        spec.name: tuple(_canonical_primitive(pattern) for pattern in spec.patterns)
        for spec in primitive_classes
    }
    assignments: dict[str, list[tuple[str, float]]] = {spec.name: [] for spec in primitive_classes}
    rejected: list[str] = []
    fixed: list[str] = []

    for label in gic_labels:
        gid = _gic_id(label)
        coeffs = _primitive_coefficients(label)
        scores = {
            spec.name: max(
                (coeffs.get(pattern, 0.0) for pattern in class_patterns[spec.name]),
                default=0.0,
            )
            for spec in primitive_classes
        }
        assigned_name = ""
        assigned_score = 0.0
        for index, spec in enumerate(primitive_classes):
            score = scores[spec.name]
            later_scores = [scores[item.name] for item in primitive_classes[index + 1 :]]
            if score >= min_fraction and all(item <= cross_fraction_max for item in later_scores):
                assigned_name = spec.name
                assigned_score = score
        if assigned_name:
            assignments[assigned_name].append((gid, assigned_score))
        else:
            fixed.append(gid)
            if any(score > 0.0 for score in scores.values()):
                rejected.append(gid)

    class_order = sorted(
        primitive_classes,
        key=lambda spec: (
            -len(assignments[spec.name]),
            -sum(score for _gid, score in assignments[spec.name]),
            spec.name,
        ),
    )
    if max_classes is not None:
        if max_classes < 0:
            raise ValueError("max_classes must be non-negative")
        class_order = class_order[:max_classes]
    kept_names = {spec.name for spec in class_order if assignments[spec.name]}

    constraints: list[ParameterClassConstraint] = []
    for spec in class_order:
        labels = tuple(gid for gid, _score in assignments[spec.name])
        if labels:
            constraints.append(ParameterClassConstraint(spec.name, labels, "shared"))

    for name, labels in assignments.items():
        if name not in kept_names:
            fixed.extend(gid for gid, _score in labels)

    selected = {pattern for constraint in constraints for pattern in constraint.patterns}
    fixed.extend(
        _gic_id(label)
        for label in gic_labels
        if _gic_id(label) not in selected and _gic_id(label) not in fixed
    )
    return DerivedPrimitiveClassPlan(
        fixed_patterns=tuple(dict.fromkeys(fixed)),
        parameter_classes=tuple(constraints),
        rejected_labels=tuple(dict.fromkeys(rejected)),
        budget=len(constraints) if max_classes is None else max_classes,
    )


def _gic_id(label: str) -> str:
    return str(label).split(None, 1)[0]


def _canonical_primitive(pattern: str) -> str:
    text = re.sub(r"\s+", "", str(pattern))
    match = re.fullmatch(r"([A-Za-z][A-Za-z0-9_]*)\(([^)]*)\)", text)
    if match is None:
        return text
    kind = match.group(1).upper()
    atoms = tuple(int(item) for item in match.group(2).split(",") if item)
    return _canonical_internal(kind, atoms)


def _canonical_internal(kind: str, atoms: tuple[int, ...]) -> str:
    if kind == "R" and len(atoms) == 2:
        atoms = tuple(sorted(atoms))
    elif kind in {"A", "B"} and len(atoms) == 3:
        atoms = min(atoms, tuple(reversed(atoms)))
    elif kind in {"D", "T"} and len(atoms) == 4:
        atoms = min(atoms, tuple(reversed(atoms)))
    args = ",".join(str(item) for item in atoms)
    return f"{kind}({args})"


def _primitive_coefficients(label: str) -> dict[str, float]:
    coeffs: dict[str, float] = {}
    for value, primitive in re.findall(
        r"([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*\*\s*"
        r"([A-Za-z][A-Za-z0-9_]*\(\s*\d+(?:\s*,\s*\d+)*\s*\))",
        str(label),
    ):
        key = _canonical_primitive(primitive)
        coeffs[key] = max(coeffs.get(key, 0.0), abs(float(value)))
    return coeffs
