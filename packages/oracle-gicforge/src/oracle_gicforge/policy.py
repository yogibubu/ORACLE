"""Shared GICForge policy constants.

This module is the Python source of truth for the frozen ORACLE GIC contract.
Fortran kernels mirror these names and policies explicitly; downstream tools
should import these constants instead of hard-coding family names locally.
"""

from __future__ import annotations

from dataclasses import dataclass


GIC_BACKEND = "oracle-native-primitive.v1"
SYCART_BACKEND = "oracle-native-cartesian-nullspace.v1"
B_MATRIX_BACKEND = "oracle-native-analytic-bmatrix.v1"
RANK_METHOD = "analytic_b_matrix_mgs_greedy"
RANK_TOLERANCE = 1.0e-7
DIAGNOSTIC_FINITE_DIFFERENCE_STEP = 1.0e-5
LINEAR_ANGLE_DEGREES = 175.0
ORDINARY_REDUCTION_CLASS = "ORDINARY"
SPECIAL_REDUCTION_CLASS = "SPECIAL_PROTECTED"
REDUCTION_POLICY = "SPECIAL_PROTECTED_FIRST_THEN_ORDINARY_ANALYTIC_RANK"
LOCAL_SYMMETRIZATION_METHOD = "LOCAL_BLOCK_SALC"
POINT_GROUP_PROJECTOR_METHOD = "POINT_GROUP_PROJECTOR"
SYMMETRIZATION_POLICY = "MERLINO_TYPE_BLOCK_SUM_AND_DIFFERENCE"
PROJECTOR_SYMMETRIZATION_POLICY = "HOMOGENEOUS_TYPE_BLOCK_POINT_GROUP_PROJECTOR"


@dataclass(frozen=True)
class PrimitiveFamilyPolicy:
    family: str
    function: str
    prefix: str
    reduction_class: str
    symmetry_block: str


PRIMITIVE_FAMILY_POLICIES = (
    PrimitiveFamilyPolicy("STRETCH", "R", "Str", ORDINARY_REDUCTION_CLASS, "STRETCH"),
    PrimitiveFamilyPolicy("BEND", "A", "Bend", ORDINARY_REDUCTION_CLASS, "BEND"),
    PrimitiveFamilyPolicy(
        "LINEAR_BEND",
        "L",
        "LinB",
        ORDINARY_REDUCTION_CLASS,
        "LINEAR_BEND",
    ),
    PrimitiveFamilyPolicy("TORSION", "D", "Tors", ORDINARY_REDUCTION_CLASS, "TORSION"),
    PrimitiveFamilyPolicy(
        "OUT_OF_PLANE",
        "U",
        "OuPl",
        ORDINARY_REDUCTION_CLASS,
        "OUT_OF_PLANE",
    ),
    PrimitiveFamilyPolicy(
        "FRAG_DISTANCE",
        "FC_DIST",
        "FCDi",
        SPECIAL_REDUCTION_CLASS,
        "SPECIAL_FRAGMENT_DISTANCE",
    ),
    PrimitiveFamilyPolicy(
        "FRAG_CENTER_ATOM_DISTANCE",
        "FCA_DIST",
        "FCAt",
        SPECIAL_REDUCTION_CLASS,
        "SPECIAL_FRAGMENT_CENTER_ATOM",
    ),
    PrimitiveFamilyPolicy(
        "FRAG_TRANSLATION",
        "FTRANS",
        "FTrn",
        SPECIAL_REDUCTION_CLASS,
        "SPECIAL_FRAGMENT_TRANSLATION",
    ),
    PrimitiveFamilyPolicy(
        "FRAG_ORIENTATION",
        "FROT",
        "FRot",
        SPECIAL_REDUCTION_CLASS,
        "SPECIAL_FRAGMENT_ORIENTATION",
    ),
    PrimitiveFamilyPolicy(
        "CENTER_ATOM_DISTANCE",
        "CENTER_ATOM_DIST",
        "CnAt",
        SPECIAL_REDUCTION_CLASS,
        "SPECIAL_CENTER_ATOM",
    ),
)

PRIMITIVE_POLICY_BY_FAMILY = {
    policy.family: policy for policy in PRIMITIVE_FAMILY_POLICIES
}
PRIMITIVE_POLICY_BY_FUNCTION = {
    policy.function: policy for policy in PRIMITIVE_FAMILY_POLICIES
}
PRIMITIVE_FAMILY_ORDER = tuple(policy.family for policy in PRIMITIVE_FAMILY_POLICIES)
SPECIAL_PRIMITIVE_FAMILIES = frozenset(
    policy.family
    for policy in PRIMITIVE_FAMILY_POLICIES
    if policy.reduction_class == SPECIAL_REDUCTION_CLASS
)
SYMMETRY_BLOCK_BY_FAMILY = {
    policy.family: policy.symmetry_block for policy in PRIMITIVE_FAMILY_POLICIES
}


def primitive_reduction_class(family: str) -> str:
    policy = PRIMITIVE_POLICY_BY_FAMILY.get(family)
    if policy is None:
        return ORDINARY_REDUCTION_CLASS
    return policy.reduction_class


def primitive_prefix(family: str) -> str:
    return PRIMITIVE_POLICY_BY_FAMILY[family].prefix


def primitive_symmetry_block(family: str) -> str:
    return SYMMETRY_BLOCK_BY_FAMILY.get(family, family)
