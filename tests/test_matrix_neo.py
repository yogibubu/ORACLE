from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

import numpy as np
import pytest

from matrix_core import section_content
from matrix_chem import (
    MolecularGeometry,
    analyze_molecular_symmetry,
    preprocess_to_enriched_xyz,
    read_enriched_xyz,
    write_validation_section,
)
from matrix_fragments import write_fragment_build_section, write_interaction_center_section
from matrix_gaussian import (
    GaussianWriteError,
    ensure_gaussian_output_pickett,
    read_gaussian_cartesian_input,
)
from matrix_neo import (
    STAGE_COORDINATE_GENERATION,
    FrozenGIC,
    GICForgeContractError,
    GICDefinition,
    GICPointGroupOperation,
    GICPrimitive,
    build_gic_b_matrix_from_xyzin,
    build_gic_definition_from_xyzin,
    coordinate_generator_by_family,
    default_coordinate_generator_registry,
    generate_stretch_coordinates,
    gaussian_gic_lines_from_xyzin,
    gic_report_from_xyzin,
    gic_b_matrix_lines,
    gic_definition_section_lines,
    special_symmetry_source_blocks,
    symmetrize_gic_definition,
    irrep_characters_for_operations,
    total_symmetric_gic_names,
    validate_coordinate_generator_registry,
    write_gicforge_build_sections,
    write_gicforge_gaussian_input,
    write_gicforge_plan_sections,
)
from matrix_neo.definition import (
    _analytic_b_row,
    _encode_ring_pucker_term,
    _finite_difference_b_row,
    _primitive_candidates,
    _primitive_value,
    _ring_pucker_component_terms,
    _select_ranked_primitives,
    _topology_bonds,
)


def _tetrahedral_methane_coordinates() -> np.ndarray:
    return np.array(
        [
            (0.0, 0.0, 0.0),
            (1.0, 1.0, 1.0),
            (1.0, -1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
        ],
        dtype=float,
    )


def _test_molecule_path(name: str) -> Path:
    return (
        Path(__file__).resolve().parents[1]
        / "tests"
        / "fixtures"
        / "test_molecules"
        / "molecules"
        / name
    )


def _skip_if_smiles_requires_rdkit(path: Path) -> None:
    header = "\n".join(path.read_text(encoding="utf-8").splitlines()[:3]).upper()
    if "SMILES" in header:
        pytest.importorskip("rdkit")


def _symmetry_operation_records(
    symmetry_section: list[str],
) -> tuple[tuple[str, tuple[int, ...], np.ndarray], ...]:
    records = []
    for line in symmetry_section:
        if " LABEL=" not in line or " PERMUTATION=" not in line or " MATRIX=" not in line:
            continue
        label = line.split(" LABEL=", 1)[1].split()[0]
        permutation_text = line.split(" PERMUTATION=", 1)[1].split()[0]
        matrix_text = line.split(" MATRIX=", 1)[1].split()[0]
        permutation = tuple(int(item) for item in permutation_text.split(","))
        matrix_values = [float(item) for item in matrix_text.split(",")]
        records.append((label, permutation, np.asarray(matrix_values, dtype=float).reshape(3, 3)))
    return tuple(records)


def _assert_symmetry_operations_closed(symmetry_section: list[str]) -> None:
    operations = _symmetry_operation_records(symmetry_section)
    keys = {
        (
            permutation,
            tuple(round(float(value), 8) for value in matrix.reshape(-1)),
        )
        for _label, permutation, matrix in operations
    }
    for _left_label, left_permutation, left_matrix in operations:
        for _right_label, right_permutation, right_matrix in operations:
            product_permutation = tuple(right_permutation[atom - 1] for atom in left_permutation)
            product_matrix = left_matrix @ right_matrix
            key = (
                product_permutation,
                tuple(round(float(value), 8) for value in product_matrix.reshape(-1)),
            )
            assert key in keys


def _test_dihedral_value(coords: np.ndarray, atoms: tuple[int, int, int, int]) -> float:
    i, j, k, ell = (atom - 1 for atom in atoms)
    p0, p1, p2, p3 = coords[i], coords[j], coords[k], coords[ell]
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def _manual_rpck_value(primitive: GICPrimitive, coords: np.ndarray) -> float:
    total = 0.0
    for ref in primitive.refs:
        coefficient_text, atoms_text = ref.split(":", 1)
        atoms = tuple(int(atom) for atom in atoms_text.split("-"))
        assert len(atoms) == 4
        total += float(coefficient_text) * _test_dihedral_value(coords, atoms)
    return float(total)


def _rpck_test_primitive(
    identifier: str,
    name: str,
    ring: tuple[int, ...],
    component_index: int,
) -> GICPrimitive:
    terms = _ring_pucker_component_terms(ring)[component_index]
    return GICPrimitive(
        identifier,
        name,
        "RING_PUCKER_COMPONENT",
        "RPCK",
        ring,
        refs=tuple(_encode_ring_pucker_term(coefficient, atoms) for coefficient, atoms in terms),
    )


def _tetrahedral_operations() -> tuple[GICPointGroupOperation, ...]:
    from itertools import permutations, product

    vertices = _tetrahedral_methane_coordinates()[1:]
    operations: list[GICPointGroupOperation] = []
    counters: Counter[str] = Counter()
    for permutation in permutations(range(3)):
        for signs in product((-1.0, 1.0), repeat=3):
            if signs[0] * signs[1] * signs[2] != 1.0:
                continue
            matrix = np.zeros((3, 3), dtype=float)
            for row, column in enumerate(permutation):
                matrix[row, column] = signs[row]
            operation_class = _tetrahedral_operation_class(matrix)
            counters[operation_class] += 1
            label = (
                "E" if operation_class == "E" else f"{operation_class}_{counters[operation_class]}"
            )
            operations.append(
                GICPointGroupOperation(
                    label,
                    tuple(tuple(float(value) for value in row) for row in matrix),
                    _tetrahedral_permutation(matrix, vertices),
                )
            )
    return tuple(sorted(operations, key=lambda operation: operation.label != "E"))


def _octahedral_operation_matrices() -> tuple[tuple[tuple[float, ...], ...], ...]:
    from itertools import permutations, product

    matrices: list[tuple[tuple[float, ...], ...]] = []
    for permutation in permutations(range(3)):
        for signs in product((-1.0, 1.0), repeat=3):
            matrix = np.zeros((3, 3), dtype=float)
            for row, column in enumerate(permutation):
                matrix[row, column] = signs[row]
            matrices.append(tuple(tuple(float(value) for value in row) for row in matrix))
    return tuple(matrices)


def _dnd_operations(n: int) -> tuple[GICPointGroupOperation, ...]:
    sd = np.array(((0.0, 1.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)))
    operations = [
        GICPointGroupOperation(
            "E",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (),
        )
    ]
    for power in range(1, n):
        matrix = _rotation_matrix((0.0, 0.0, 1.0), 2.0 * np.pi * power / n)
        operations.append(_operation_without_permutation(f"C{n}z^{power}", matrix))
    for k in range(n):
        matrix = _rotation_matrix((np.cos(np.pi * k / n), np.sin(np.pi * k / n), 0.0), np.pi)
        operations.append(_operation_without_permutation(f"C2_xy_{n}_{k}", matrix))
    for power in range(n):
        matrix = sd @ _rotation_matrix((0.0, 0.0, 1.0), 2.0 * np.pi * power / n)
        operations.append(
            _operation_without_permutation(
                f"Dnd_C2_xy_{2 * n}_{(2 * power + 1) % (2 * n)}",
                matrix,
            )
        )
    for k in range(n):
        c2 = _rotation_matrix((np.cos(np.pi * k / n), np.sin(np.pi * k / n), 0.0), np.pi)
        operations.append(
            _operation_without_permutation(
                f"Dnd_C{2 * n}z^{(2 * k + 1) % (2 * n)}",
                sd @ c2,
            )
        )
    return tuple(operations)


def _operation_without_permutation(label: str, matrix: np.ndarray) -> GICPointGroupOperation:
    return GICPointGroupOperation(
        label,
        tuple(tuple(float(value) for value in row) for row in matrix),
        (),
    )


def _rotation_matrix(axis: tuple[float, float, float], theta: float) -> np.ndarray:
    axis_array = np.array(axis, dtype=float)
    axis_array /= np.linalg.norm(axis_array)
    x, y, z = axis_array
    c = np.cos(theta)
    s = np.sin(theta)
    one_c = 1.0 - c
    return np.array(
        [
            [c + x * x * one_c, x * y * one_c - z * s, x * z * one_c + y * s],
            [y * x * one_c + z * s, c + y * y * one_c, y * z * one_c - x * s],
            [z * x * one_c - y * s, z * y * one_c + x * s, c + z * z * one_c],
        ],
        dtype=float,
    )


def _d2d_stretch_operations() -> tuple[GICPointGroupOperation, ...]:
    ligand_coords = _tetrahedral_methane_coordinates()[1:]
    operations = []
    for operation in _dnd_operations(2):
        matrix = np.asarray(operation.rotation, dtype=float)
        transformed = ligand_coords @ matrix.T
        permutation = [1]
        for position in transformed:
            matches = np.where(np.all(np.isclose(ligand_coords, position, atol=1.0e-8), axis=1))[0]
            assert len(matches) == 1
            permutation.append(int(matches[0]) + 2)
        operations.append(
            GICPointGroupOperation(operation.label, operation.rotation, tuple(permutation))
        )
    return tuple(operations)


def _icosahedral_vertices() -> np.ndarray:
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    vertices = []
    for y in (-1.0, 1.0):
        for z in (-phi, phi):
            vertices.append((0.0, y, z))
    for x in (-1.0, 1.0):
        for y in (-phi, phi):
            vertices.append((x, y, 0.0))
    for x in (-phi, phi):
        for z in (-1.0, 1.0):
            vertices.append((x, 0.0, z))
    return np.array(vertices, dtype=float)


def _tetrahedral_operation_class(matrix: np.ndarray) -> str:
    if np.allclose(matrix, np.eye(3)):
        return "E"
    det = float(np.linalg.det(matrix))
    trace = float(np.trace(matrix))
    if det > 0.0 and abs(trace) <= 1.0e-8:
        return "C3"
    if det > 0.0 and abs(trace + 1.0) <= 1.0e-8:
        return "C2"
    if det < 0.0 and abs(trace + 1.0) <= 1.0e-8:
        return "S4"
    if det < 0.0 and abs(trace - 1.0) <= 1.0e-8:
        return "sigma_d"
    raise AssertionError(f"unexpected Td operation:\n{matrix}")


def _tetrahedral_permutation(matrix: np.ndarray, vertices: np.ndarray) -> tuple[int, ...]:
    transformed = vertices @ matrix.T
    permutation = [1]
    for position in transformed:
        matches = np.where(np.all(np.isclose(vertices, position, atol=1.0e-8), axis=1))[0]
        assert len(matches) == 1
        permutation.append(int(matches[0]) + 2)
    return tuple(permutation)


def _assert_character_rows_orthonormal(
    irreps: dict[str, tuple[float, ...]],
    *,
    group_order: int,
) -> None:
    names = tuple(irreps)
    for left in names:
        for right in names:
            overlap = sum(float(a) * float(b) for a, b in zip(irreps[left], irreps[right])) / float(
                group_order
            )
            assert overlap == pytest.approx(1.0 if left == right else 0.0, abs=1.0e-8)


def _assert_merlino_character_rows(
    point_group: str,
    labels: tuple[str, ...],
    expected: tuple[tuple[str, tuple[float, ...]], ...],
) -> None:
    rows = irrep_characters_for_operations(labels, point_group)
    assert tuple(name for name, _chars in rows) == tuple(name for name, _chars in expected)
    for (name, chars), (expected_name, expected_chars) in zip(rows, expected):
        assert name == expected_name
        assert chars == pytest.approx(expected_chars, abs=1.0e-8)


def test_coordinate_generator_registry_is_deterministic_and_stage_separated():
    registry = default_coordinate_generator_registry()
    by_family = coordinate_generator_by_family()

    assert validate_coordinate_generator_registry() == ()
    assert all(entry.stage == STAGE_COORDINATE_GENERATION for entry in registry)
    assert len({entry.name for entry in registry}) == len(registry)
    assert len({entry.coordinate_family for entry in registry}) == len(registry)
    assert {
        "STRETCH",
        "LOCAL_XH_STRETCH",
        "BEND",
        "LOCAL_ANGLE_SALC",
        "TORSION",
        "RING_PUCKER_COMPONENT",
        "BUTTERFLY",
        "SPECIAL_COORDINATE",
    } <= set(by_family)
    assert by_family["LOCAL_XH_STRETCH"].implemented_by == (
        "matrix_neo.generators.generate_stretch_coordinates"
    )
    assert "through 9" in by_family["LOCAL_ANGLE_SALC"].notes


def _stretch_generator_signature(records):
    return tuple((item.family, item.function, item.atoms) for item in records)


def _stretch_candidate_signature(candidates):
    return tuple(
        (item.family, item.function, item.atoms)
        for item in candidates
        if item.function == "R" and item.family in {"STRETCH", "LOCAL_XH_STRETCH"}
    )


def _assert_stretch_generator_matches_primitive_candidates(
    bonds,
    coords,
    atom_symbols,
    *,
    xh_stretch_policy="SYMMETRIZE",
    local_xh_classes=(),
):
    generated = generate_stretch_coordinates(
        tuple(bonds),
        atom_symbols=tuple(atom_symbols),
        xh_stretch_policy=xh_stretch_policy,
        local_xh_classes=tuple(local_xh_classes),
    )
    candidates = _primitive_candidates(
        tuple(bonds),
        coords=np.asarray(coords, dtype=float),
        natoms=len(atom_symbols),
        atom_symbols=tuple(atom_symbols),
        xh_stretch_policy=xh_stretch_policy,
        local_xh_classes=tuple(local_xh_classes),
    )

    assert _stretch_generator_signature(generated) == _stretch_candidate_signature(candidates)


def test_stretch_generator_matches_primitive_candidates_for_small_xh_cases():
    _assert_stretch_generator_matches_primitive_candidates(
        ((1, 2), (1, 3)),
        np.asarray([(0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0)]),
        ("O", "H", "H"),
        xh_stretch_policy="LOCAL_SELECTED",
        local_xh_classes=("XH2",),
    )
    _assert_stretch_generator_matches_primitive_candidates(
        ((1, 2), (1, 3), (1, 4)),
        np.asarray(
            [
                (0.0, 0.0, 0.0),
                (0.94, 0.0, 0.0),
                (-0.313, 0.885, 0.0),
                (-0.313, -0.442, 0.766),
            ]
        ),
        ("N", "H", "H", "H"),
        xh_stretch_policy="LOCAL_SELECTED",
        local_xh_classes=("XH3",),
    )
    _assert_stretch_generator_matches_primitive_candidates(
        ((1, 2), (1, 3), (1, 4), (1, 5)),
        _tetrahedral_methane_coordinates(),
        ("C", "H", "H", "H", "H"),
        xh_stretch_policy="SYMMETRIZE",
    )


def test_stretch_generator_matches_primitive_candidates_for_ring_fixtures(tmp_path):
    for filename in ("benzene.inp", "pyrrole.inp"):
        source = _test_molecule_path(filename)
        xyzin = tmp_path / f"{filename}.xyzin"
        preprocess_to_enriched_xyz(source, xyzin)
        write_validation_section(xyzin)
        lines = xyzin.read_text(encoding="utf-8").splitlines()
        geometry = read_enriched_xyz(xyzin)
        bonds = _topology_bonds(lines, natoms=geometry.natoms)

        _assert_stretch_generator_matches_primitive_candidates(
            bonds,
            geometry.coordinates_angstrom,
            geometry.atoms,
            xh_stretch_policy="SYMMETRIZE",
        )


def test_gicforge_plan_requires_validation_pass(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#VALIDATION",
                "SCHEMA oracle.xyz.validation.v1",
                "STATUS FAIL",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(GICForgeContractError, match="status must be PASS"):
        write_gicforge_plan_sections(path)


def test_gicforge_plan_writes_gic_and_sycart_sections(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#VALIDATION",
                "SCHEMA oracle.xyz.validation.v1",
                "STATUS PASS",
                "",
                "#TOPOLOGY",
                "SCHEMA oracle.xyz.topology.v1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    write_gicforge_plan_sections(path, symmetrize=True, sycart=True)
    lines = path.read_text(encoding="utf-8").splitlines()

    gic = section_content(lines, "GIC")
    sycart = section_content(lines, "SYCART")
    assert gic[0] == "SCHEMA oracle.xyz.gic.v1"
    assert "STATUS PLANNED" in gic
    assert "SYMMETRIZE TRUE" in gic
    assert "OUT_OF_PLANE_MODE OUT_OF_PLANE" in gic
    assert "FRAGMENT_MODE SPECIAL_COORDINATES" in gic
    assert "PENDING GICFORGE_IMPLEMENTATION" in gic
    assert sycart[0] == "SCHEMA oracle.xyz.sycart.v1"


def test_gicforge_plan_can_request_improper_dihedrals(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#VALIDATION",
                "SCHEMA oracle.xyz.validation.v1",
                "STATUS PASS",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    write_gicforge_plan_sections(path, improper_dihedrals=True)
    gic = section_content(path.read_text(encoding="utf-8").splitlines(), "GIC")

    assert "OUT_OF_PLANE_MODE IMPROPER_DIHEDRAL" in gic


def test_gicforge_plan_can_request_pseudo_bond_fragment_mode(tmp_path):
    path = tmp_path / "molecule.xyz"
    path.write_text(
        "\n".join(
            [
                "1",
                "h",
                "H 0 0 0",
                "",
                "#VALIDATION",
                "SCHEMA oracle.xyz.validation.v1",
                "STATUS PASS",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    write_gicforge_plan_sections(path, fragment_mode="h-bonds")
    gic = section_content(path.read_text(encoding="utf-8").splitlines(), "GIC")

    assert "FRAGMENT_MODE PSEUDO_BONDS" in gic


def test_preprocess_validate_then_gicforge_plan_pipeline(tmp_path):
    source = tmp_path / "water.xyz"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    target = tmp_path / "water.xyzin"

    preprocess_to_enriched_xyz(source, target)
    validation = write_validation_section(target)
    write_gicforge_plan_sections(target, symmetrize=True)
    lines = target.read_text(encoding="utf-8").splitlines()

    assert validation.status == "PASS"
    assert "STATUS PASS" in section_content(lines, "VALIDATION")
    assert "SYMMETRIZE TRUE" in section_content(lines, "GIC")


def test_gicforge_writes_gaussian_input_after_gic_plan(tmp_path):
    source = tmp_path / "water.xyz"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "water.xyzin"
    gjf = tmp_path / "water.gjf"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_gicforge_plan_sections(xyzin, symmetrize=True)
    write_gicforge_gaussian_input(
        xyzin,
        gjf,
        route="#p b3lyp/sto-3g opt",
        title="water from ORACLE",
    )

    text = gjf.read_text(encoding="utf-8")
    geometry = read_gaussian_cartesian_input(gjf)
    assert "#p b3lyp/sto-3g opt=readallgic output=pickett" in text
    assert "! ORACLE_SCHEMA" not in text
    assert geometry.atoms == ("O", "H", "H")
    assert geometry.comment == "water from ORACLE"


def test_gaussian_route_normalization_enforces_pickett_output():
    assert ensure_gaussian_output_pickett("#p hf/sto-3g") == "#p hf/sto-3g output=pickett"
    assert (
        ensure_gaussian_output_pickett("#p hf/sto-3g output=wfn")
        == "#p hf/sto-3g output=pickett"
    )


def test_fragment_delta_gaussian_inputs_enforce_pickett_output(tmp_path):
    from matrix_neo.survibfit.fragment_delta_correction import prepare_low_level_gaussian_inputs

    xyz = tmp_path / "fragment.xyz"
    xyz.write_text("1\nfragment\nN 0.0 0.0 0.0\n", encoding="utf-8")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "entries": [
                    {
                        "query_fragment_id": "frag1",
                        "suggested_low_level_input_xyz": str(xyz),
                    }
                ]
            }
        )
        + "\n",
        encoding="utf-8",
    )

    result = prepare_low_level_gaussian_inputs(
        manifest,
        route="#p hf/cc-pVTZ output=wfn",
        nproc=1,
    )

    assert result["generated_inputs"] == 1
    gjf = Path(result["inputs"][0])
    assert "#p hf/cc-pVTZ output=pickett" in gjf.read_text(encoding="utf-8")


def test_gicforge_build_writes_frozen_gics_and_sycart(tmp_path):
    source = tmp_path / "water.xyz"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "water.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True, sycart=True)
    lines = xyzin.read_text(encoding="utf-8").splitlines()
    gic = section_content(lines, "GIC")
    sycart = section_content(lines, "SYCART")

    assert definition.target_rank == 3
    assert definition.rank == 3
    assert gic[0] == "SCHEMA oracle.xyz.gic.v1"
    assert "STATUS BUILT" in gic
    assert "BACKEND merlino-python-gicforge.v1" in gic
    assert "POINT_GROUP C2v" in gic
    assert "SYMMETRY_MODE POINT_GROUP_PROJECTOR" in gic
    assert "SYMMETRY_GROUP C2v" in gic
    assert "TOTAL_SYMMETRIC_IRREP A1" in gic
    assert "TOTAL_SYMMETRIC_GIC_COUNT 2" in gic
    assert "TOTAL_SYMMETRIC_GICS A1Str001,A1Bend001" in gic
    assert "GIC_COUNT 3" in gic
    assert not any("PENDING" in line for line in gic)
    assert "[SYMMETRY_DIAGNOSTICS]" in gic
    assert "METHOD POINT_GROUP_PROJECTOR" in gic
    assert "STATUS APPLIED" in gic
    assert any(
        line.startswith(
            "GROUP 1 BLOCK=STRETCH FAMILY=STRETCH "
            "SIGNATURE=OPS=E,sigma_yz,sigma_xy,C2y^1 "
            "SOURCES=Stre0001,Stre0002 OUTPUTS=A1Str001,B2Str001"
        )
        for line in gic
    )
    assert any(
        line.startswith(
            "GIC001 NAME=A1Str001 FAMILY=STRETCH IRREP=A1 "
            "COEFFS=P001:0.707106781187,P002:0.707106781187"
        )
        for line in gic
    )
    assert "A1Str001 = 0.707106781187*(R(1,2))+0.707106781187*(R(1,3))" in gic
    assert "B2Str001 = 0.707106781187*(R(1,2))-0.707106781187*(R(1,3))" in gic
    assert "A1Bend001 = A(2,1,3)" in gic
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.status == "APPLIED"
    assert definition.symmetry_diagnostics.symmetry_group == "C2v"
    assert definition.symmetry_diagnostics.total_symmetric_irrep == "A1"
    assert definition.symmetry_diagnostics.total_symmetric_gics == (
        "A1Str001",
        "A1Bend001",
    )
    assert total_symmetric_gic_names(definition) == (
        "A1Str001",
        "A1Bend001",
    )
    assert definition.gics[0].name == "A1Str001"
    assert definition.gics[0].coefficients[0][1] == pytest.approx(1.0 / np.sqrt(2.0))
    matrix = build_gic_b_matrix_from_xyzin(xyzin)
    assert matrix.coordinate_names[:2] == ("A1Str001", "B2Str001")
    assert matrix.rows[0][5] == pytest.approx(1.0 / np.sqrt(2.0))
    assert matrix.rows[0][7] == pytest.approx(1.0 / np.sqrt(2.0))
    assert matrix.rows[1][5] == pytest.approx(1.0 / np.sqrt(2.0))
    assert matrix.rows[1][7] == pytest.approx(-1.0 / np.sqrt(2.0))
    report_lines = gic_report_from_xyzin(xyzin)
    assert "Method: POINT_GROUP_PROJECTOR" in report_lines
    assert "Total irrep: A1" in report_lines
    assert "Target rank rationale: 3N-6 non-linear vibrational rank for N=3" in report_lines
    assert "Mode: NONE" in report_lines
    assert "Policy: no built fragments were consumed by this GIC definition." in report_lines
    assert (
        "STRETCH OPS=E,sigma_yz,sigma_xy,C2y^1: Stre0001,Stre0002 -> A1Str001,B2Str001"
    ) in report_lines
    assert "Selected by family: BEND:1, STRETCH:2" in report_lines
    assert any(
        line.startswith("A1Str001 irrep=A1 family=STRETCH")
        and "P001:+0.707106781187" in line
        and "P002:+0.707106781187" in line
        for line in report_lines
    )
    assert sycart[0] == "SCHEMA oracle.xyz.sycart.v1"
    assert "STATUS BUILT" in sycart
    assert "COORD_COUNT 3" in sycart
    assert any(line.startswith("SYC001 IRREP=A COMPONENTS=") for line in sycart)


def test_gicforge_build_definition_uses_saved_topology(tmp_path):
    source = tmp_path / "water.xyz"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "water.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = build_gic_definition_from_xyzin(xyzin)

    assert [primitive.gaussian_expression() for primitive in definition.primitives] == [
        "R(1,2)",
        "R(1,3)",
        "A(2,1,3)",
    ]


def test_xh_stretch_policy_defaults_to_ordinary_stretches(tmp_path):
    source = tmp_path / "water.xyz"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "water.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    gic = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "GIC")

    assert "LOCAL_XH_STRETCH" not in {primitive.family for primitive in definition.primitives}
    assert "XH_STRETCH_POLICY SYMMETRIZE" in gic


def test_xh_stretch_policy_can_keep_xh2_class_local(tmp_path):
    source = tmp_path / "water.xyz"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "water.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(
        xyzin,
        xh_stretch_policy="local-selected",
        local_xh_classes=("XH2",),
    )
    gic = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "GIC")
    families = Counter(primitive.family for primitive in definition.primitives)

    assert families["LOCAL_XH_STRETCH"] == 2
    assert families["BEND"] == 1
    assert "LOCAL_XH_CLASSES XH2" in gic


def test_xh_stretch_policy_can_keep_xh3_class_local(tmp_path):
    source = tmp_path / "ammonia.xyz"
    source.write_text(
        "\n".join(
            [
                "4",
                "ammonia",
                "N 0.00000000 0.00000000 0.00000000",
                "H 0.94000000 0.00000000 0.00000000",
                "H -0.31300000 0.88500000 0.00000000",
                "H -0.31300000 -0.44200000 0.76600000",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "ammonia.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(
        xyzin,
        xh_stretch_policy="local-selected",
        local_xh_classes=("XH3",),
    )

    assert Counter(primitive.family for primitive in definition.primitives)["LOCAL_XH_STRETCH"] == 3


def test_gicforge_build_gaussian_input_includes_readgic_block(tmp_path):
    source = tmp_path / "water.xyz"
    source.write_text(
        "\n".join(
            [
                "3",
                "water",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "water.xyzin"
    gjf = tmp_path / "water.gjf"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_gicforge_build_sections(xyzin)
    write_gicforge_gaussian_input(xyzin, gjf, route="#p hf/sto-3g")

    text = gjf.read_text(encoding="utf-8")
    assert gaussian_gic_lines_from_xyzin(xyzin) == [
        "GIC001 = R(1,2)",
        "GIC002 = R(1,3)",
        "GIC003 = A(2,1,3)",
    ]
    assert "#p hf/sto-3g geom=readallgic output=pickett" in text
    assert "GIC001 = R(1,2)" in text
    assert "GIC003 = A(2,1,3)" in text

    write_gicforge_gaussian_input(xyzin, gjf, route="#p hf/sto-3g opt")
    text = gjf.read_text(encoding="utf-8")
    assert "#p hf/sto-3g opt=readallgic output=pickett" in text

    write_gicforge_gaussian_input(
        xyzin,
        gjf,
        route="#p hf/sto-3g opt=(MaxCycle=1) geom=readallgic",
    )
    text = gjf.read_text(encoding="utf-8")
    assert "#p hf/sto-3g opt=(readallgic,MaxCycle=1) output=pickett" in text
    assert "geom=readallgic" not in text.lower()

    with pytest.raises(GaussianWriteError):
        write_gicforge_gaussian_input(
            xyzin,
            gjf,
            route="#p hf/sto-3g geom=(readallgic,gicallsym)",
        )
    assert "gicallsym" not in text.lower()


@pytest.mark.parametrize(
    ("molecule", "expected_point_group", "expected_rank", "expected_tokens"),
    [
        ("pyrrole.inp", "C2v", 24, ("A2RPck001", "B1RPck001", "D(")),
        ("azulene.inp", "C2v", 48, ("BtFl", "RPck", "D(")),
        ("spiro.inp", "D2", 87, ("ASpir001", "ARPck001", "D(")),
        ("cubane.inp", "Oh", 42, ("A2uBtFl001", "D(")),
        ("ferrocene.inp", "D5h", 57, ("A1pCnAt001", "CxC031(Inactive)", "SQRT(")),
        (
            "ferrocene_staggered.inp",
            "D5d",
            57,
            ("A1gCnAt001", "A2uCnAt001", "CxC031(Inactive)", "SQRT("),
        ),
    ],
)
def test_gicforge_gaussian_input_smoke_for_symmetrized_nontrivial_cases(
    tmp_path,
    molecule,
    expected_point_group,
    expected_rank,
    expected_tokens,
):
    source = _test_molecule_path(molecule)
    _skip_if_smiles_requires_rdkit(source)
    xyzin = tmp_path / f"{Path(molecule).stem}.xyzin"
    gjf = tmp_path / f"{Path(molecule).stem}.gjf"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    if molecule.startswith("ferrocene"):
        write_interaction_center_section(xyzin)
    definition = write_gicforge_build_sections(
        xyzin,
        symmetrize=True,
        sycart=True,
        improper_dihedrals=True,
    )
    write_gicforge_gaussian_input(
        xyzin,
        gjf,
        route="#p hf/sto-3g opt=(MaxCycle=1) scf=tight",
        title=f"MATRIX symmetrized GIC smoke: {Path(molecule).stem}",
    )
    text = gjf.read_text(encoding="utf-8")
    definition_lines = [
        line
        for line in text.splitlines()
        if "=" in line and not line.lstrip().startswith(("#", "!", "%"))
    ]
    definition_names = [
        line.split("=", 1)[0].strip().replace("(Inactive)", "") for line in definition_lines
    ]

    assert definition.point_group == expected_point_group
    assert definition.rank == expected_rank
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.status == "APPLIED"
    assert "#p hf/sto-3g opt=(readallgic,MaxCycle=1) scf=tight output=pickett" in text
    assert "geom=readallgic" not in text.lower()
    assert "! ORACLE_SCHEMA" not in text
    assert len(definition_lines) >= len(definition.gics)
    assert "U(" not in text
    assert "QPck" not in text
    assert "PhiP" not in text
    assert "gicallsym" not in text.lower()
    assert "PENDING" not in text
    assert "UNASSIGNED" not in text
    assert not any("'" in name for name in definition_names)
    for token in expected_tokens:
        assert token in text


def test_gicforge_native_analytic_b_rows_match_diagnostic_finite_difference():
    coords = np.asarray(
        [
            [0.00, 0.00, 0.00],
            [1.20, 0.00, 0.00],
            [1.80, 1.05, 0.20],
            [2.55, 1.20, 1.05],
        ],
        dtype=float,
    )
    linear_coords = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.001, 0.0],
        ],
        dtype=float,
    )
    primitives = [
        (GICPrimitive("P001", "Str", "STRETCH", "R", (1, 2)), coords),
        (GICPrimitive("P002", "Bend", "BEND", "A", (1, 2, 3)), coords),
        (GICPrimitive("P003", "Tors", "TORSION", "D", (1, 2, 3, 4)), coords),
        (GICPrimitive("P004", "OuPl", "OUT_OF_PLANE", "U", (2, 1, 3, 4)), coords),
        (
            GICPrimitive("P005", "LinB", "LINEAR_BEND", "L", (1, 2, 3), mode=-1),
            linear_coords,
        ),
    ]

    for primitive, primitive_coords in primitives:
        analytic = _analytic_b_row(primitive, primitive_coords)
        diagnostic_fd = _finite_difference_b_row(
            primitive,
            primitive_coords,
            step_angstrom=1.0e-6,
        )
        assert np.max(np.abs(analytic - diagnostic_fd)) < 1.0e-5


def test_gicforge_generates_out_of_plane_or_improper_dihedral_mode():
    coords = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.2],
            [0.0, 1.0, -0.1],
            [0.0, 0.0, 1.1],
        ],
        dtype=float,
    )
    bonds = ((1, 2), (1, 3), (1, 4))

    oop_candidates = _primitive_candidates(bonds, coords=coords, natoms=4)
    oop = next(primitive for primitive in oop_candidates if primitive.family == "OUT_OF_PLANE")
    assert oop.function == "U"
    assert oop.name == "OuPl0001"
    assert oop.gaussian_expression() == "U(1,2,3,4)"

    improper_candidates = _primitive_candidates(
        bonds,
        coords=coords,
        natoms=4,
        improper_dihedrals=True,
    )
    improper = next(
        primitive for primitive in improper_candidates if primitive.family == "IMPROPER_DIHEDRAL"
    )
    assert improper.function == "IMPD"
    assert improper.name == "ImpD0001"
    assert improper.gaussian_expression() == "D(2,1,4,3)"
    assert _primitive_value(improper, coords) == pytest.approx(
        _test_dihedral_value(coords, (2, 1, 4, 3))
    )
    assert (
        np.max(
            np.abs(_analytic_b_row(improper, coords) - _finite_difference_b_row(improper, coords))
        )
        < 1.0e-5
    )


def test_gicforge_skips_all_cyclic_out_of_plane_candidates_like_merlino():
    coords = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-0.4, -0.3, 0.8],
        ],
        dtype=float,
    )
    bonds = ((1, 2), (1, 3), (1, 4))
    rings = ((1, (1, 2, 3, 4)),)

    oop_candidates = _primitive_candidates(
        bonds,
        rings=rings,
        coords=coords,
        natoms=4,
    )
    improper_candidates = _primitive_candidates(
        bonds,
        rings=rings,
        coords=coords,
        natoms=4,
        improper_dihedrals=True,
    )

    assert "OUT_OF_PLANE" not in {primitive.family for primitive in oop_candidates}
    assert "IMPROPER_DIHEDRAL" not in {primitive.family for primitive in improper_candidates}


def test_gicforge_build_handles_corpus_zmatrix_case(tmp_path):
    source = (
        Path(__file__).resolve().parents[1]
        / "tests"
        / "fixtures"
        / "test_molecules"
        / "molecules"
        / "h2cozmat.inp"
    )
    xyzin = tmp_path / "h2co.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)

    assert definition.target_rank == 6
    assert definition.rank == 6
    assert len(definition.gics) == 6


def test_gicforge_corpus_fused_ring_case_keeps_merlino_ring_families(tmp_path):
    source = _test_molecule_path("naphtalene.inp")
    xyzin = tmp_path / "naphtalene.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    families = Counter(primitive.family for primitive in definition.primitives)

    assert definition.target_rank == 48
    assert definition.rank == 48
    assert families["CYCLIC_BEND"] > 0
    assert families["RING_PUCKER_COMPONENT"] > 0
    assert families["BUTTERFLY"] > 0


def test_gicforge_ring_puckering_coefficients_match_merlino_six_ring():
    ring = (6, 4, 1, 3, 5, 2)
    components = _ring_pucker_component_terms(ring)

    assert len(components) == 3
    assert [term[1] for term in components[0]] == [
        (2, 6, 4, 1),
        (6, 4, 1, 3),
        (4, 1, 3, 5),
        (1, 3, 5, 2),
        (3, 5, 2, 6),
        (5, 2, 6, 4),
    ]
    np.testing.assert_allclose(
        [[coefficient for coefficient, _atoms in component] for component in components],
        [
            [
                0.5773502691896257,
                -0.28867513459481275,
                -0.2886751345948132,
                0.5773502691896257,
                -0.2886751345948125,
                -0.28867513459481325,
            ],
            [
                0.0,
                0.5,
                -0.5,
                0.0,
                0.5,
                -0.5,
            ],
            [
                0.408248290463863,
                -0.408248290463863,
                0.408248290463863,
                -0.408248290463863,
                0.408248290463863,
                -0.408248290463863,
            ],
        ],
        atol=1.0e-14,
    )


@pytest.mark.parametrize(
    ("molecule", "expected_rank", "expected_rpck_count"),
    [
        ("naphtalene.inp", 48, 6),
        ("phenantrene.inp", 66, 9),
        ("pyrene.inp", 72, 8),
    ],
)
def test_gicforge_ring_puckering_numeric_corpus_fused(
    tmp_path,
    molecule,
    expected_rank,
    expected_rpck_count,
):
    source = _test_molecule_path(molecule)
    xyzin = tmp_path / f"{molecule}.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    rpck_gics = [gic for gic in definition.gics if gic.family == "RING_PUCKER_COMPONENT"]
    gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin)

    assert definition.target_rank == expected_rank
    assert definition.rank == expected_rank
    assert len(rpck_gics) == expected_rpck_count
    if expected_rpck_count == 0:
        assert not any(line.startswith("RPck") for line in gaussian_lines)
        return
    assert any(line.startswith("RPck") for line in gaussian_lines)
    assert not any(line.startswith("QPck") for line in gaussian_lines)
    assert not any(line.startswith("PhiP") for line in gaussian_lines)
    primitive_by_id = {primitive.identifier: primitive for primitive in definition.primitives}
    for gic in rpck_gics:
        assert len(gic.coefficients) > 1
        assert all(
            primitive_by_id[primitive_id].family == "RING_PUCKER_COMPONENT"
            and primitive_by_id[primitive_id].function == "D"
            for primitive_id, _coefficient in gic.coefficients
        )


def test_gicforge_build_uses_built_fragments_for_relative_coordinates(tmp_path):
    source = tmp_path / "dimer.xyz"
    source.write_text(
        "\n".join(
            [
                "6",
                "water dimer fragments",
                "O 0.00 0.00 0.00",
                "H 0.96 0.00 0.00",
                "H -0.24 0.93 0.00",
                "O 0.00 0.00 3.20",
                "H 0.96 0.00 3.20",
                "H -0.24 0.93 3.20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "dimer.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_fragment_build_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    gic = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "GIC")
    primitive_lines = [line for line in gic if line.startswith("P")]
    gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin)

    assert definition.target_rank == 12
    assert definition.rank == 12
    assert len(definition.gics) == 12
    assert "RANK_METHOD analytic_b_matrix_mgs_greedy" in gic
    assert "REDUCTION_POLICY SPECIAL_PROTECTED_FIRST_THEN_ORDINARY_ANALYTIC_RANK" in gic
    assert "[REDUCTION_DIAGNOSTICS]" in gic
    assert any(line.startswith("SELECTED P") for line in gic)
    assert any(
        "FAMILY=FRAG_DISTANCE CLASS=SPECIAL_PROTECTED FUNCTION=FC_DIST" in line
        for line in primitive_lines
    )
    assert any(
        "FAMILY=FRAG_CENTER_ATOM_DISTANCE CLASS=SPECIAL_PROTECTED FUNCTION=FCA_DIST" in line
        for line in primitive_lines
    )
    assert any(
        "FAMILY=FRAG_ORIENTATION CLASS=SPECIAL_PROTECTED FUNCTION=FROT" in line
        for line in primitive_lines
    )
    assert any("REFS=F002,F001" in line for line in primitive_lines)
    assert "F001=Fragment(1-3)" in gaussian_lines
    assert "F002=Fragment(4-6)" in gaussian_lines
    assert "CxF001(Inactive)=XCntr(F001)" in gaussian_lines
    assert any("= SQRT((CxF001-CxF002)**2" in line for line in gaussian_lines)
    assert any("= SQRT((CxF002-X(1))**2" in line for line in gaussian_lines)
    assert any("KxF002F001" in line for line in gaussian_lines)
    assert any("ExF002F001" in line for line in gaussian_lines)
    assert any("1.0D-24" in line and line.startswith("KnF002F001") for line in gaussian_lines)
    assert any(" = ExF002F001" in line for line in gaussian_lines)
    report_lines = gic_report_from_xyzin(xyzin)
    assert "Mode: SPECIAL_COORDINATES" in report_lines
    assert (
        "Policy: automatic weak-complex default; keep fragments as protected bodies."
        in report_lines
    )
    assert any("family=FRAG_DISTANCE" in line for line in report_lines)


def test_gicforge_pseudo_bond_fragment_mode_uses_standard_internal_coordinates(tmp_path):
    source = tmp_path / "formic_acid_water.xyz"
    source.write_text(
        "\n".join(
            [
                "8",
                "formic acid-water non-covalent probe",
                "6    -1.171727   -0.018999   -0.001370",
                "1    -2.256869    0.130369   -0.015107",
                "8    -0.651376   -1.113045    0.004964",
                "8    -0.533544    1.143228    0.007603",
                "1     0.435203    0.960633    0.019343",
                "8     1.930145   -0.025976   -0.052269",
                "1     2.594643   -0.128917    0.633345",
                "1     1.351317   -0.802634    0.008831",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "formic_acid_water.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_fragment_build_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, fragment_mode="pseudo-bonds")
    matrix = build_gic_b_matrix_from_xyzin(xyzin)
    rows = np.asarray(matrix.rows, dtype=float)
    gic = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "GIC")

    assert definition.fragment_mode == "PSEUDO_BONDS"
    assert definition.pseudo_bonds == ((5, 6),)
    assert definition.pseudo_bond_kinds == ("HBOND",)
    assert definition.target_rank == 18
    assert definition.rank == 18
    assert np.linalg.matrix_rank(rows, tol=1.0e-8) == 18
    assert {gic.family for gic in definition.gics} == {"STRETCH", "BEND", "TORSION"}
    assert not any(gic.family.startswith("FRAG_") for gic in definition.gics)
    assert any(
        primitive.family == "STRETCH" and primitive.atoms == (5, 6)
        for primitive in definition.primitives
    )
    assert "FRAGMENT_MODE PSEUDO_BONDS" in gic
    assert "PSEUDO_BOND_COUNT 1" in gic
    assert "[PSEUDO_BONDS]" in gic
    assert "1 5 6 KIND=HBOND" in gic
    report_lines = gic_report_from_xyzin(xyzin)
    assert "Mode: PSEUDO_BONDS" in report_lines
    assert (
        "Policy: explicit graph-joining mode; do not build protected fragment coordinates."
        in report_lines
    )
    assert "Pseudo-bonds: 5-6:HBOND" in report_lines


def test_gicforge_non_covalent_probe_runs_fragment_and_hbond_coordinate_models(tmp_path):
    xyz = "\n".join(
        [
            "8",
            "formic acid-water non-covalent probe",
            "6    -1.171727   -0.018999   -0.001370",
            "1    -2.256869    0.130369   -0.015107",
            "8    -0.651376   -1.113045    0.004964",
            "8    -0.533544    1.143228    0.007603",
            "1     0.435203    0.960633    0.019343",
            "8     1.930145   -0.025976   -0.052269",
            "1     2.594643   -0.128917    0.633345",
            "1     1.351317   -0.802634    0.008831",
        ]
    )
    results = {}

    for mode in ("special-coordinates", "pseudo-bonds"):
        source = tmp_path / f"{mode}.xyz"
        xyzin = tmp_path / f"{mode}.xyzin"
        source.write_text(xyz + "\n", encoding="utf-8")
        preprocess_to_enriched_xyz(source, xyzin)
        write_validation_section(xyzin)
        write_fragment_build_section(xyzin)
        definition = write_gicforge_build_sections(xyzin, fragment_mode=mode)
        matrix = build_gic_b_matrix_from_xyzin(xyzin)
        rows = np.asarray(matrix.rows, dtype=float)
        results[mode] = definition

        assert definition.target_rank == 18
        assert definition.rank == 18
        assert np.linalg.matrix_rank(rows, tol=1.0e-8) == 18

    fragment_definition = results["special-coordinates"]
    fragment_families = {gic.family for gic in fragment_definition.gics}
    assert fragment_definition.fragment_mode == "SPECIAL_COORDINATES"
    assert {
        "FRAG_DISTANCE",
        "FRAG_CENTER_ATOM_DISTANCE",
        "FRAG_TRANSLATION",
        "FRAG_ORIENTATION",
    }.issubset(fragment_families)

    hbond_definition = results["pseudo-bonds"]
    hbond_families = {gic.family for gic in hbond_definition.gics}
    assert hbond_definition.fragment_mode == "PSEUDO_BONDS"
    assert hbond_definition.pseudo_bonds == ((5, 6),)
    assert hbond_definition.pseudo_bond_kinds == ("HBOND",)
    assert hbond_families == {"STRETCH", "BEND", "TORSION"}
    assert not any(
        primitive.family.startswith(("FRAG_", "RING_", "CYCLIC_"))
        for primitive in hbond_definition.primitives
    )


def test_gicforge_classifies_ring_and_butterfly_primitives_like_merlino():
    coords = np.asarray(
        [
            (-1.2, 1.4, 0.0),
            (-2.0, 0.7, 0.0),
            (-2.0, -0.7, 0.0),
            (-1.2, -1.4, 0.0),
            (0.0, -0.7, 0.0),
            (0.0, 0.7, 0.0),
            (1.2, 1.4, 0.0),
            (2.0, 0.7, 0.0),
            (2.0, -0.7, 0.0),
            (1.2, -1.4, 0.0),
        ],
        dtype=float,
    )
    bonds = (
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (1, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 10),
        (5, 10),
    )
    rings = (
        (1, (1, 2, 3, 4, 5, 6)),
        (2, (5, 6, 7, 8, 9, 10)),
    )

    candidates = _primitive_candidates(bonds, rings=rings, coords=coords, natoms=10)
    families = Counter(primitive.family for primitive in candidates)
    by_atoms = {primitive.atoms: primitive.family for primitive in candidates}

    assert families["CYCLIC_BEND"] > 0
    assert families["RING_PUCKER_COMPONENT"] == 6
    assert families["CONDENSED_RING_TORSION"] > 0
    assert families["BUTTERFLY"] > 0
    assert by_atoms[(4, 5, 6, 7)] == "BUTTERFLY"
    assert by_atoms[(4, 5, 10, 9)] == "CONDENSED_RING_TORSION"
    assert by_atoms[(1, 2, 3, 4, 5, 6)] == "RING_PUCKER_COMPONENT"
    assert by_atoms[(5, 6, 7, 8, 9, 10)] == "RING_PUCKER_COMPONENT"


def test_gicforge_symmetrizes_special_fragment_coordinates(tmp_path):
    source = tmp_path / "dimer.xyz"
    source.write_text(
        "\n".join(
            [
                "6",
                "water dimer fragments",
                "O 0.00 0.00 0.00",
                "H 0.96 0.00 0.00",
                "H -0.24 0.93 0.00",
                "O 0.00 0.00 3.20",
                "H 0.96 0.00 3.20",
                "H -0.24 0.93 3.20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "dimer.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_fragment_build_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    matrix = build_gic_b_matrix_from_xyzin(xyzin)
    gic = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "GIC")

    assert definition.symmetry_diagnostics is not None
    special_groups = [
        group
        for group in definition.symmetry_diagnostics.groups
        if group.block.startswith("SPECIAL_")
    ]
    assert special_groups
    assert special_groups[0].block == "SPECIAL_FRAGMENT_CENTER_ATOM"
    assert special_groups[0].output_gics == ("A1FCAtS001", "A2FCAtD001")
    assert "TOTAL_SYMMETRIC_IRREP A1" in gic
    assert "A1FCAtS001" in total_symmetric_gic_names(definition)
    assert "A2FCAtD001" not in total_symmetric_gic_names(definition)
    assert any(
        line.startswith(
            "GIC003 NAME=A1FCAtS001 FAMILY=FRAG_CENTER_ATOM_DISTANCE IRREP=A1 "
            "COEFFS=P007:0.707106781187,P008:0.707106781187"
        )
        for line in gic
    )
    assert any(line.startswith("A1FCAtS001 = 0.707106781187*(") for line in gic)
    assert matrix.coordinate_names[2:4] == ("A1FCAtS001", "A2FCAtD001")
    assert len(matrix.rows) == definition.rank


def test_gicforge_b_matrix_includes_fragment_coordinate_rows(tmp_path):
    source = tmp_path / "dimer.xyz"
    source.write_text(
        "\n".join(
            [
                "6",
                "water dimer fragments",
                "O 0.00 0.00 0.00",
                "H 0.96 0.00 0.00",
                "H -0.24 0.93 0.00",
                "O 0.00 0.00 3.20",
                "H 0.96 0.00 3.20",
                "H -0.24 0.93 3.20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    xyzin = tmp_path / "dimer.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    write_fragment_build_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    matrix = build_gic_b_matrix_from_xyzin(xyzin)

    assert matrix.coordinate_labels == tuple(f"GIC{idx:03d}" for idx in range(1, 13))
    assert matrix.cartesian_columns[:6] == ("1:X", "1:Y", "1:Z", "2:X", "2:Y", "2:Z")
    assert len(matrix.rows) == 12
    assert len(matrix.rows[0]) == 18
    assert all(value == value for row in matrix.rows for value in row)

    center_distance_row = matrix.rows[matrix.coordinate_names.index("FCDi0001")]
    assert center_distance_row[2] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[5] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[8] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[11] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[14] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert center_distance_row[17] == pytest.approx(1.0 / 3.0, abs=1.0e-8)

    translation_row = matrix.rows[matrix.coordinate_names.index("FTrn0001")]
    assert translation_row[0] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert translation_row[3] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert translation_row[6] == pytest.approx(-1.0 / 3.0, abs=1.0e-8)
    assert translation_row[9] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert translation_row[12] == pytest.approx(1.0 / 3.0, abs=1.0e-8)
    assert translation_row[15] == pytest.approx(1.0 / 3.0, abs=1.0e-8)

    lines = gic_b_matrix_lines(matrix)
    assert lines[0] == "SCHEMA oracle.gic.bmatrix.v1"
    assert "BACKEND oracle-native-analytic-bmatrix.v1" in lines
    assert "DERIVATIVE_MODE ANALYTIC" in lines
    assert "ROW_COUNT 12" in lines
    assert "COLUMN_COUNT 18" in lines
    center_distance_label = matrix.coordinate_labels[matrix.coordinate_names.index("FCDi0001")]
    assert any(
        line.startswith(f"{center_distance_label} NAME=FCDi0001 IRREP=UNASSIGNED VALUES=")
        for line in lines
    )

    coords = np.asarray(definition.reference_coordinates_angstrom, dtype=float)
    frot = next(primitive for primitive in definition.primitives if primitive.function == "FROT")
    analytic = _analytic_b_row(frot, coords)
    diagnostic_fd = _finite_difference_b_row(frot, coords, step_angstrom=1.0e-6)
    assert np.max(np.abs(analytic - diagnostic_fd)) < 1.0e-5

    blocks = special_symmetry_source_blocks(definition)
    assert {block.block for block in blocks} >= {
        "SPECIAL_FRAGMENT_DISTANCE",
        "SPECIAL_FRAGMENT_CENTER_ATOM",
        "SPECIAL_FRAGMENT_TRANSLATION",
        "SPECIAL_FRAGMENT_ORIENTATION",
    }

    report_lines = gic_report_from_xyzin(xyzin)
    assert "ORACLE GICForge Report" in report_lines
    assert "Rank method: analytic_b_matrix_mgs_greedy" in report_lines
    assert any("SPECIAL_FRAGMENT_DISTANCE" in line for line in report_lines)


def test_gicforge_center_atom_distance_uses_analytic_b_and_gaussian_center():
    primitive = GICPrimitive(
        identifier="P001",
        name="CnAt0001",
        family="CENTER_ATOM_DISTANCE",
        function="CENTER_ATOM_DIST",
        atoms=(1, 2),
        ref_atoms=(3,),
        refs=("C001", "A3"),
    )
    coords = np.asarray(
        [
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 2.0],
        ],
        dtype=float,
    )

    analytic = _analytic_b_row(primitive, coords)
    diagnostic_fd = _finite_difference_b_row(primitive, coords, step_angstrom=1.0e-6)
    assert np.max(np.abs(analytic - diagnostic_fd)) < 1.0e-8
    assert analytic[2] == pytest.approx(-0.5)
    assert analytic[5] == pytest.approx(-0.5)
    assert analytic[8] == pytest.approx(1.0)

    definition = GICDefinition(
        backend="test",
        point_group="C1",
        symmetrize=False,
        target_rank=1,
        rank=1,
        candidate_count=1,
        reference_coordinates_angstrom=tuple(tuple(row) for row in coords),
        primitives=(primitive,),
        gics=(
            FrozenGIC(
                identifier="GIC001",
                name="CnAt0001",
                family="CENTER_ATOM_DISTANCE",
                irrep="A",
                primitive_id="P001",
                gaussian_expression="NONE",
                coefficients=(("P001", 1.0),),
            ),
        ),
    )
    gaussian_lines = gic_definition_section_lines(definition)
    assert "CxC001(Inactive)=(X(1)+X(2))/2" in gaussian_lines
    assert any(line.startswith("GIC001 = SQRT((CxC001-X(3))**2") for line in gaussian_lines)


def test_gicforge_rank_reduction_protects_special_coordinates():
    ordinary = GICPrimitive(
        identifier="P001",
        name="Str0001",
        family="STRETCH",
        function="R",
        atoms=(1, 2),
    )
    special = GICPrimitive(
        identifier="P002",
        name="CnAt0001",
        family="CENTER_ATOM_DISTANCE",
        function="CENTER_ATOM_DIST",
        atoms=(1,),
        ref_atoms=(2,),
        refs=("C001", "A2"),
    )
    coords = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )

    selected, rank = _select_ranked_primitives(
        (ordinary, special),
        coords,
        target_rank=1,
        rank_tolerance=1.0e-7,
    )

    assert rank == 1
    assert selected == (special,)
    assert selected[0].reduction_class == "SPECIAL_PROTECTED"


def test_gicforge_symmetry_labels_use_point_group_irreps():
    primitives = tuple(
        GICPrimitive(
            identifier=f"P{idx:03d}",
            name=f"Str{idx:04d}",
            family="STRETCH",
            function="R",
            atoms=(1, idx + 1),
        )
        for idx in range(1, 5)
    )
    definition = GICDefinition(
        backend="test",
        point_group="C2v",
        symmetrize=False,
        target_rank=4,
        rank=4,
        candidate_count=4,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 5,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("O", "H", "H", "H", "H"),
    )

    assert [gic.irrep for gic in symmetrized.gics] == ["A1", "A2", "B1", "B2"]
    assert [gic.name for gic in symmetrized.gics] == [
        "A1StrS001",
        "A2StrD001",
        "B1StrD001",
        "B2StrD001",
    ]
    assert total_symmetric_gic_names(symmetrized) == ("A1StrS001",)


def test_gicforge_irrep_characters_cover_merlino_family_groups():
    c3v = dict(
        irrep_characters_for_operations(
            ("E", "C3z^1", "C3z^2", "sigma_v_3_0"),
            "C3v",
        )
    )
    assert c3v["A1"] == pytest.approx((1.0, 1.0, 1.0, 1.0))
    assert c3v["A2"] == pytest.approx((1.0, 1.0, 1.0, -1.0))
    assert c3v["E"] == pytest.approx((2.0, -1.0, -1.0, 0.0))

    c4v = dict(
        irrep_characters_for_operations(
            ("E", "C4z^1", "C4z^2", "C4z^3", "sigma_v_4_0", "sigma_v_4_1"),
            "C4v",
        )
    )
    assert c4v["B1"] == pytest.approx((1.0, -1.0, 1.0, -1.0, 1.0, -1.0))
    assert c4v["B2"] == pytest.approx((1.0, -1.0, 1.0, -1.0, -1.0, 1.0))
    assert c4v["E"] == pytest.approx((2.0, 0.0, -2.0, 0.0, 0.0, 0.0), abs=1.0e-12)

    c2h = dict(irrep_characters_for_operations(("E", "C2z", "i", "sigma_xy"), "C2h"))
    assert c2h["Ag"] == pytest.approx((1.0, 1.0, 1.0, 1.0))
    assert c2h["Bu"] == pytest.approx((1.0, -1.0, -1.0, 1.0))

    d2h = dict(
        irrep_characters_for_operations(
            ("E", "C2z", "C2y", "C2x", "i", "sigma_xy", "sigma_xz", "sigma_yz"),
            "D2h",
        )
    )
    assert d2h["B1g"] == pytest.approx((1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0))
    assert d2h["B3u"] == pytest.approx((1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0))

    d3h = dict(
        irrep_characters_for_operations(
            (
                "E",
                "C3z^1",
                "C3z^2",
                "C2_xy_3_0",
                "sigma_xy",
                "sigma_h*C3z^1",
                "sigma_v_3_0",
            ),
            "D3h",
        )
    )
    assert d3h["A1'"] == pytest.approx((1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
    assert d3h["A2'"] == pytest.approx((1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0))
    assert d3h["E'"] == pytest.approx(
        (2.0, -1.0, -1.0, 0.0, 2.0, -1.0, 0.0),
        abs=1.0e-12,
    )

    d6h_labels = (
        "E",
        "i",
        "sigma_yz",
        "sigma_xz",
        "sigma_xy",
        "C2z^1",
        "C2x^1",
        "C2y^1",
        "C3z^1",
        "sigma_h*C3z^1",
        "C3z^2",
        "sigma_h*C3z^2",
        "C6z^1",
        "sigma_h*C6z^1",
        "C6z^5",
        "sigma_h*C6z^5",
        "C2_xy_3_1",
        "C2_xy_3_2",
        "C2_xy_6_1",
        "C2_xy_6_5",
        "sigma_v_3_1",
        "sigma_v_3_2",
        "sigma_v_6_1",
        "sigma_v_6_5",
    )
    d6h = dict(irrep_characters_for_operations(d6h_labels, "D6h"))
    assert d6h["A1g"] == pytest.approx((1.0,) * 24)
    assert d6h["E2g"][0] == pytest.approx(2.0)
    assert d6h["E1u"][0] == pytest.approx(2.0)
    _assert_character_rows_orthonormal(d6h, group_order=24)

    d2d_ops = _dnd_operations(2)
    d2d = dict(
        irrep_characters_for_operations(
            tuple(operation.label for operation in d2d_ops),
            "D2d",
            operation_matrices=tuple(operation.rotation for operation in d2d_ops),
        )
    )
    assert tuple(d2d) == ("A1", "A2", "B1", "B2", "E")
    _assert_character_rows_orthonormal(d2d, group_order=8)

    d3d_ops = _dnd_operations(3)
    d3d = dict(
        irrep_characters_for_operations(
            tuple(operation.label for operation in d3d_ops),
            "D3d",
            operation_matrices=tuple(operation.rotation for operation in d3d_ops),
        )
    )
    assert tuple(d3d) == ("A1g", "A2g", "Eg", "A1u", "A2u", "Eu")
    _assert_character_rows_orthonormal(d3d, group_order=12)


def test_gicforge_label_only_characters_match_merlino3():
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    phi_bar = (1.0 - np.sqrt(5.0)) / 2.0
    cases = (
        (
            "Td",
            ("E", "C3_t", "C3_t2", "C2_t", "sigma_yz"),
            (
                ("A1", (1.0, 1.0, 1.0, 1.0, 1.0)),
                ("A2", (1.0, 1.0, 1.0, 1.0, -1.0)),
                ("E", (2.0, -1.0, -1.0, 2.0, 0.0)),
                ("T1", (3.0, 0.0, 0.0, -1.0, -1.0)),
                ("T2", (3.0, 0.0, 0.0, -1.0, 1.0)),
            ),
        ),
        (
            "Oh",
            ("E", "C3_o", "C2_o", "C4_o", "C4_o2", "i", "sigma_yz", "S4"),
            (
                ("A1g", (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
                ("A1u", (1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0)),
                ("A2g", (1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0)),
                ("A2u", (1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0)),
                ("Eg", (2.0, -1.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0)),
                ("Eu", (2.0, -1.0, 2.0, 0.0, 0.0, -2.0, 0.0, 0.0)),
                ("T1g", (3.0, 0.0, -1.0, 1.0, 1.0, 3.0, -1.0, 1.0)),
                ("T1u", (3.0, 0.0, -1.0, 1.0, 1.0, -3.0, 1.0, -1.0)),
                ("T2g", (3.0, 0.0, -1.0, -1.0, -1.0, 3.0, 1.0, -1.0)),
                ("T2u", (3.0, 0.0, -1.0, -1.0, -1.0, -3.0, -1.0, 1.0)),
            ),
        ),
        (
            "Ih",
            ("E", "C5_i", "C5_i2", "C5_i3", "C5_i4", "C3_i", "C2_i", "i"),
            (
                ("Ag", (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
                ("T1g", (3.0, phi, phi_bar, phi_bar, phi, 0.0, -1.0, 3.0)),
                ("T2g", (3.0, phi_bar, phi, phi, phi_bar, 0.0, -1.0, 3.0)),
                ("Gg", (4.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.0, 4.0)),
                ("Hg", (5.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 5.0)),
                ("Au", (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0)),
                ("T1u", (3.0, phi, phi_bar, phi_bar, phi, 0.0, -1.0, -3.0)),
                ("T2u", (3.0, phi_bar, phi, phi, phi_bar, 0.0, -1.0, -3.0)),
                ("Gu", (4.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.0, -4.0)),
                ("Hu", (5.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, -5.0)),
            ),
        ),
        (
            "D2d",
            ("E", "C2z", "C2_xy_2_0", "C2_xy_2_1", "S4"),
            (
                ("A1'", (1.0, 1.0, 1.0, 1.0, 1.0)),
                ("A1''", (1.0, 1.0, 1.0, 1.0, -1.0)),
                ("A2'", (1.0, 1.0, -1.0, -1.0, 1.0)),
                ("A2''", (1.0, 1.0, -1.0, -1.0, -1.0)),
                ("B1'", (1.0, -1.0, -1.0, -1.0, 1.0)),
                ("B1''", (1.0, -1.0, -1.0, -1.0, -1.0)),
                ("B2'", (1.0, -1.0, 1.0, 1.0, 1.0)),
                ("B2''", (1.0, -1.0, 1.0, 1.0, -1.0)),
            ),
        ),
        (
            "D3d",
            ("E", "C3z^1", "C3z^2", "C2_xy_3_0", "i"),
            (
                ("A1g", (1.0, 1.0, 1.0, 1.0, 1.0)),
                ("A1u", (1.0, 1.0, 1.0, 1.0, -1.0)),
                ("A2g", (1.0, 1.0, 1.0, -1.0, 1.0)),
                ("A2u", (1.0, 1.0, 1.0, -1.0, -1.0)),
                ("E1g", (4.0, -2.0, -2.0, 0.0, 4.0)),
                ("E1u", (4.0, -2.0, -2.0, 0.0, -4.0)),
            ),
        ),
    )

    for point_group, labels, expected in cases:
        _assert_merlino_character_rows(point_group, labels, expected)


def test_oracle_symmetry_detects_tetrahedral_operations():
    coords = _tetrahedral_methane_coordinates()
    symmetry = analyze_molecular_symmetry(
        MolecularGeometry(
            atoms=("C", "H", "H", "H", "H"),
            coordinates_angstrom=coords,
            comment="methane td",
        ),
        distance_tolerance=1.0e-5,
        inertia_tolerance=1.0e-5,
        max_rotation_order=6,
    )

    assert symmetry.point_group == "Td"
    assert len(symmetry.operations) == 24


def test_oracle_symmetry_detects_icosahedral_operations():
    vertices = _icosahedral_vertices()
    symmetry = analyze_molecular_symmetry(
        MolecularGeometry(
            atoms=tuple("B" for _ in range(len(vertices))),
            coordinates_angstrom=vertices,
            comment="icosahedron ih",
        ),
        distance_tolerance=1.0e-5,
        inertia_tolerance=1.0e-5,
        max_rotation_order=6,
    )

    assert symmetry.point_group == "Ih"
    assert len(symmetry.operations) == 120


def test_gicforge_polyhedral_irrep_characters_use_operation_matrices():
    operations = _tetrahedral_operations()
    irreps = dict(
        irrep_characters_for_operations(
            tuple(operation.label for operation in operations),
            "Td",
            operation_matrices=tuple(operation.rotation for operation in operations),
        )
    )

    assert tuple(irreps) == ("A1", "A2", "E", "T1", "T2")
    assert Counter(round(value) for value in irreps["E"]) == Counter({2: 4, -1: 8, 0: 12})
    assert Counter(round(value) for value in irreps["T2"]) == Counter({3: 1, 0: 8, -1: 9, 1: 6})


def test_gicforge_octahedral_irrep_characters_use_operation_matrices():
    matrices = _octahedral_operation_matrices()
    labels = tuple(
        "E" if np.allclose(matrix, np.eye(3)) else f"Oh{idx}" for idx, matrix in enumerate(matrices)
    )
    irreps = dict(
        irrep_characters_for_operations(
            labels,
            "Oh",
            operation_matrices=matrices,
        )
    )

    assert tuple(irreps) == (
        "A1g",
        "A2g",
        "Eg",
        "T1g",
        "T2g",
        "A1u",
        "A2u",
        "Eu",
        "T1u",
        "T2u",
    )
    assert Counter(round(value) for value in irreps["A1g"]) == Counter({1: 48})
    assert Counter(round(value) for value in irreps["A1u"]) == Counter({1: 24, -1: 24})
    assert irreps["T1u"] == pytest.approx(
        tuple(float(np.trace(np.asarray(matrix, dtype=float))) for matrix in matrices)
    )


def test_gicforge_icosahedral_irrep_characters_use_operation_matrices():
    vertices = _icosahedral_vertices()
    symmetry = analyze_molecular_symmetry(
        MolecularGeometry(
            atoms=tuple("B" for _ in range(len(vertices))),
            coordinates_angstrom=vertices,
            comment="icosahedron ih",
        ),
        distance_tolerance=1.0e-5,
        inertia_tolerance=1.0e-5,
        max_rotation_order=6,
    )
    irreps = dict(
        irrep_characters_for_operations(
            tuple(operation.label for operation in symmetry.operations),
            "Ih",
            operation_matrices=tuple(operation.rotation for operation in symmetry.operations),
        )
    )

    assert tuple(irreps) == (
        "Ag",
        "T1g",
        "T2g",
        "Gg",
        "Hg",
        "Au",
        "T1u",
        "T2u",
        "Gu",
        "Hu",
    )
    _assert_character_rows_orthonormal(irreps, group_order=120)
    assert irreps["T1u"] == pytest.approx(
        tuple(
            float(np.trace(np.asarray(operation.rotation, dtype=float)))
            for operation in symmetry.operations
        )
    )


def test_gicforge_point_group_projector_uses_operations_without_type_mixing():
    primitives = (
        GICPrimitive("P001", "Str0001", "STRETCH", "R", (1, 2)),
        GICPrimitive("P002", "Str0002", "STRETCH", "R", (1, 3)),
        GICPrimitive("P003", "Bend0001", "BEND", "A", (2, 1, 3)),
    )
    definition = GICDefinition(
        backend="test",
        point_group="C2v",
        symmetrize=False,
        target_rank=3,
        rank=3,
        candidate_count=3,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 3,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )
    operations = (
        GICPointGroupOperation(
            "E",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (1, 2, 3),
        ),
        GICPointGroupOperation(
            "sigma_yz",
            ((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (1, 3, 2),
        ),
        GICPointGroupOperation(
            "sigma_xy",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)),
            (1, 2, 3),
        ),
        GICPointGroupOperation(
            "C2y^1",
            ((-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)),
            (1, 3, 2),
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("O", "H", "H"),
        symmetry_operations=operations,
    )

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert [gic.name for gic in symmetrized.gics] == [
        "A1Str001",
        "B2Str001",
        "A1Bend001",
    ]
    assert [gic.irrep for gic in symmetrized.gics] == ["A1", "B2", "A1"]
    assert [
        (group.block, group.family, group.output_gics)
        for group in symmetrized.symmetry_diagnostics.groups
    ] == [
        ("STRETCH", "STRETCH", ("A1Str001", "B2Str001")),
        ("BEND", "BEND", ("A1Bend001",)),
    ]
    assert total_symmetric_gic_names(symmetrized) == ("A1Str001", "A1Bend001")


def test_gicforge_point_group_projector_symmetrizes_ring_puckering_components():
    ring_left = (1, 2, 3, 4, 5, 6)
    ring_right = (7, 8, 9, 10, 11, 12)
    primitives = (
        _rpck_test_primitive("P001", "RPck0001", ring_left, 0),
        _rpck_test_primitive("P002", "RPck0002", ring_left, 1),
        _rpck_test_primitive("P003", "RPck0003", ring_right, 0),
        _rpck_test_primitive("P004", "RPck0004", ring_right, 1),
    )
    definition = GICDefinition(
        backend="test",
        point_group="C2",
        symmetrize=False,
        target_rank=4,
        rank=4,
        candidate_count=4,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 12,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression="NONE",
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )
    operations = (
        GICPointGroupOperation(
            "E",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            tuple(range(1, 13)),
        ),
        GICPointGroupOperation(
            "C2z^1",
            ((-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
            (7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6),
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("C",) * 12,
        symmetry_operations=operations,
    )
    gaussian_lines = gic_definition_section_lines(symmetrized)

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert [gic.name for gic in symmetrized.gics] == [
        "ARPck001",
        "ARPck002",
        "BRPck001",
        "BRPck002",
    ]
    assert [gic.irrep for gic in symmetrized.gics] == ["A", "A", "B", "B"]
    assert symmetrized.symmetry_diagnostics.groups[0].block == "RING_PUCKER_COMPONENT"
    assert symmetrized.symmetry_diagnostics.groups[0].family == "RING_PUCKER_COMPONENT"
    assert total_symmetric_gic_names(symmetrized) == ("ARPck001", "ARPck002")
    assert any(line.startswith("ARPck001 =") for line in gaussian_lines)
    assert any(line.startswith("ARPck002 =") for line in gaussian_lines)
    assert any(line.startswith("BRPck001 =") for line in gaussian_lines)
    assert any(line.startswith("BRPck002 =") for line in gaussian_lines)
    assert not any("RPck" in line and "(Inactive)" in line for line in gaussian_lines)
    assert not any(line.startswith("QPck") for line in gaussian_lines)
    assert not any(line.startswith("PhiP") for line in gaussian_lines)


def test_gicforge_benzene_d6h_projector_symmetrizes_ring_puckering(tmp_path):
    source = _test_molecule_path("benzene.inp")
    xyzin = tmp_path / "benzene.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin)
    gic = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "GIC")
    rpck = [gic for gic in definition.gics if gic.family == "RING_PUCKER_COMPONENT"]

    assert definition.point_group == "D6h"
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert "POINT_GROUP D6h" in gic
    assert "SYMMETRY_MODE POINT_GROUP_PROJECTOR" in gic
    assert [gic.irrep for gic in rpck] == ["B2g", "E2u", "E2u"]
    assert [gic.name for gic in rpck] == ["B2gRPck001", "E2uRPck001", "E2uRPck002"]
    assert not any(name.startswith("A1gRPck") for name in total_symmetric_gic_names(definition))
    assert any(line.startswith("B2gRPck001 =") for line in gaussian_lines)
    assert any(line.startswith("E2uRPck001 =") for line in gaussian_lines)
    assert any(line.startswith("E2uRPck002 =") for line in gaussian_lines)
    assert not any("RPck" in line and "(Inactive)" in line for line in gaussian_lines)
    assert not any(line.startswith("QPck") for line in gaussian_lines)
    assert not any(line.startswith("PhiP") for line in gaussian_lines)


def test_gicforge_pyrrole_c2v_projector_keeps_ring_coordinates(tmp_path):
    source = _test_molecule_path("pyrrole.inp")
    xyzin = tmp_path / "pyrrole.xyzin"
    gjf = tmp_path / "pyrrole.gjf"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin)
    total_gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin, total_symmetric_only=True)
    gic = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "GIC")
    rpck = [gic for gic in definition.gics if gic.family == "RING_PUCKER_COMPONENT"]
    oop = [gic for gic in definition.gics if gic.family == "OUT_OF_PLANE"]
    cyclic_bends = [gic for gic in definition.gics if gic.family == "CYCLIC_BEND"]

    assert definition.point_group == "C2v"
    assert definition.target_rank == 24
    assert definition.rank == 24
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert definition.symmetry_diagnostics.status == "APPLIED"
    assert "POINT_GROUP C2v" in gic
    assert "SYMMETRY_MODE POINT_GROUP_PROJECTOR" in gic
    assert [gic.name for gic in rpck] == ["A2RPck001", "B1RPck001"]
    assert [gic.irrep for gic in rpck] == ["A2", "B1"]
    assert [gic.irrep for gic in oop] == ["A2", "A2", "B1", "B1", "B1"]
    assert [gic.name for gic in cyclic_bends] == ["A1CyBe001", "B2CyBe001"]
    assert [gic.irrep for gic in cyclic_bends] == ["A1", "B2"]
    assert not any(
        any(token in name for token in ("RPck", "OuPl"))
        for name in total_symmetric_gic_names(definition)
    )
    assert any(line.startswith("A2RPck001 =") for line in gaussian_lines)
    assert any(line.startswith("B1RPck001 =") for line in gaussian_lines)
    assert any(
        line.startswith("A2RPck001(Frozen) =")
        for line in gaussian_gic_lines_from_xyzin(
            xyzin,
            freeze_non_total=True,
        )
    )
    assert not any("RPck" in line and "(Inactive)" in line for line in gaussian_lines)
    assert not any(line.startswith("QPck") for line in gaussian_lines)
    assert not any(line.startswith("PhiP") for line in gaussian_lines)
    assert not any("RPck" in line or "OuPl" in line for line in total_gaussian_lines)

    write_gicforge_gaussian_input(xyzin, gjf, route="#p hf/sto-3g opt=readallgic")
    text = gjf.read_text(encoding="utf-8")
    assert "#p hf/sto-3g opt=readallgic output=pickett" in text
    assert "A1Str001 = R(1,2)" in text
    assert "A1CyBe001 = " in text
    assert "A2RPck001(Frozen) = " in text
    assert "B1RPck001(Frozen) = " in text
    assert "A2OuPl001(Frozen) = " in text
    assert "B1OuPl003(Frozen) = " in text
    assert "A1Str001(Frozen)" not in text
    assert "A1CyBe001(Frozen)" not in text


def test_gicforge_projector_handles_large_noncondensed_ring_puckering(tmp_path):
    source = _test_molecule_path("cyclottane.inp")
    xyzin = tmp_path / "cyclottane.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin)
    matrix = build_gic_b_matrix_from_xyzin(xyzin)
    rows = np.asarray(matrix.rows, dtype=float)
    rpck = [gic for gic in definition.gics if gic.family == "RING_PUCKER_COMPONENT"]
    cyclic_bends = [gic for gic in definition.gics if gic.family == "CYCLIC_BEND"]

    assert definition.point_group == "D2"
    assert definition.target_rank == 66
    assert definition.rank == 66
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert definition.symmetry_diagnostics.status == "APPLIED"
    assert np.linalg.matrix_rank(rows, tol=1.0e-8) == 66
    assert len(rpck) == 5
    assert len(cyclic_bends) == 5
    assert not any(line.startswith("QPck") for line in gaussian_lines)
    assert not any(line.startswith("PhiP") for line in gaussian_lines)
    assert all(group.block == group.family for group in definition.symmetry_diagnostics.groups)


@pytest.mark.parametrize(
    (
        "molecule",
        "expected_group",
        "expected_rank",
        "expected_rpck_count",
        "expected_symmetry_method",
    ),
    [
        ("naphtalene.inp", "D2h", 48, 6, "POINT_GROUP_PROJECTOR"),
        ("phenantrene.inp", "C2v", 66, 9, "POINT_GROUP_PROJECTOR"),
        ("anthracene.inp", "D2h", 66, 9, "POINT_GROUP_PROJECTOR"),
        ("pyrene.inp", "D2h", 72, 8, "POINT_GROUP_PROJECTOR"),
        ("fluorene.inp", "C2v", 63, 8, "POINT_GROUP_PROJECTOR"),
    ],
)
def test_gicforge_projector_handles_fused_ring_corpus_without_fallback(
    tmp_path,
    molecule,
    expected_group,
    expected_rank,
    expected_rpck_count,
    expected_symmetry_method,
):
    source = _test_molecule_path(molecule)
    xyzin = tmp_path / f"{Path(molecule).stem}.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    rpck = [gic for gic in definition.gics if gic.family == "RING_PUCKER_COMPONENT"]

    assert definition.point_group == expected_group
    assert definition.target_rank == expected_rank
    assert definition.rank == expected_rank
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.method == expected_symmetry_method
    assert definition.symmetry_diagnostics.status == "APPLIED"
    assert len(rpck) == expected_rpck_count
    assert all(gic.name.startswith(gic.irrep) for gic in rpck)
    if expected_symmetry_method == "POINT_GROUP_PROJECTOR":
        assert all(group.block == group.family for group in definition.symmetry_diagnostics.groups)


@pytest.mark.parametrize(
    (
        "molecule",
        "expected_group",
        "expected_rank",
        "expected_rpck_count",
        "expected_btfl_count",
        "expected_symmetry_method",
    ),
    [
        ("azulene.inp", "C2v", 48, 6, 1, "POINT_GROUP_PROJECTOR"),
        ("norbornane.inp", "C2v", 51, 3, 0, "POINT_GROUP_PROJECTOR"),
        ("norbornene.inp", "Cs", 45, 3, 0, "POINT_GROUP_PROJECTOR"),
        ("norbornadiene.inp", "C2v", 39, 3, 0, "POINT_GROUP_PROJECTOR"),
        ("norcamphor.inp", "C1", 48, 3, 0, "LOCAL_BLOCK_SALC"),
    ],
)
def test_gicforge_projector_handles_bridged_and_nonbenzenoid_ring_probe(
    tmp_path,
    molecule,
    expected_group,
    expected_rank,
    expected_rpck_count,
    expected_btfl_count,
    expected_symmetry_method,
):
    source = _test_molecule_path(molecule)
    _skip_if_smiles_requires_rdkit(source)
    xyzin = tmp_path / f"{Path(molecule).stem}.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin)
    rpck = [gic for gic in definition.gics if gic.family == "RING_PUCKER_COMPONENT"]
    btfl = [gic for gic in definition.gics if gic.family == "BUTTERFLY"]

    assert definition.point_group == expected_group
    assert definition.target_rank == expected_rank
    assert definition.rank == expected_rank
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.method == expected_symmetry_method
    assert definition.symmetry_diagnostics.status == "APPLIED"
    assert len(rpck) == expected_rpck_count
    assert len(btfl) == expected_btfl_count
    assert not any(line.startswith("QPck") for line in gaussian_lines)
    assert not any(line.startswith("PhiP") for line in gaussian_lines)
    if expected_symmetry_method == "POINT_GROUP_PROJECTOR":
        assert all(group.block == group.family for group in definition.symmetry_diagnostics.groups)


def test_gicforge_projector_handles_spiro_ring_probe(tmp_path):
    source = _test_molecule_path("spiro.inp")
    xyzin = tmp_path / "spiro.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin)
    spiro = [gic for gic in definition.gics if gic.family == "SPIRO_BEND"]
    rpck = [gic for gic in definition.gics if gic.family == "RING_PUCKER_COMPONENT"]

    assert definition.point_group == "D2"
    assert definition.target_rank == 87
    assert definition.rank == 87
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert definition.symmetry_diagnostics.status == "APPLIED"
    assert len(spiro) == 3
    assert len(rpck) == 6
    assert all(gic.name.startswith(gic.irrep) for gic in spiro + rpck)
    assert any("Spir" in line for line in gaussian_lines)
    assert not any(line.startswith("QPck") for line in gaussian_lines)
    assert not any(line.startswith("PhiP") for line in gaussian_lines)
    assert all(group.block == group.family for group in definition.symmetry_diagnostics.groups)


def test_gicforge_projector_handles_cubane_oh_without_rank_loss(tmp_path):
    source = _test_molecule_path("cubane.inp")
    xyzin = tmp_path / "cubane.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    matrix = build_gic_b_matrix_from_xyzin(xyzin)
    rows = np.asarray(matrix.rows, dtype=float)

    assert definition.point_group == "Oh"
    assert definition.target_rank == 42
    assert definition.rank == 42
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert definition.symmetry_diagnostics.status == "APPLIED"
    assert np.linalg.matrix_rank(rows, tol=1.0e-8) == 42
    families = {group.family for group in definition.symmetry_diagnostics.groups}
    assert {"BUTTERFLY", "CYCLIC_BEND"}.issubset(families)


def test_gicforge_ferrocene_uses_metal_ring_center_coordinates(tmp_path):
    source = _test_molecule_path("ferrocene.inp")
    xyzin = tmp_path / "ferrocene.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    centers = write_interaction_center_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    matrix = build_gic_b_matrix_from_xyzin(xyzin)
    gaussian_lines = gaussian_gic_lines_from_xyzin(xyzin)
    symmetry_section = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "SYMMETRY")
    rows = np.asarray(matrix.rows, dtype=float)

    ring_center_ids = {
        center.identifier for center in centers.centers if center.kind == "RING_CENTER"
    }
    center_atom_primitives = tuple(
        primitive
        for primitive in definition.primitives
        if primitive.family == "CENTER_ATOM_DISTANCE"
    )
    center_atom_gics = tuple(gic for gic in definition.gics if gic.family == "CENTER_ATOM_DISTANCE")

    assert definition.point_group == "D5h"
    assert any(line == "POINT_GROUP D5h" for line in symmetry_section)
    assert sum(1 for line in symmetry_section if " LABEL=" in line and "PERMUTATION=" in line) == 20
    _assert_symmetry_operations_closed(symmetry_section)
    assert definition.target_rank == 57
    assert definition.rank == 57
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert np.linalg.matrix_rank(rows, tol=1.0e-8) == 57
    assert len(center_atom_primitives) == 2
    assert len(center_atom_gics) == 2
    assert {primitive.refs[0] for primitive in center_atom_primitives} == ring_center_ids
    assert all(primitive.ref_atoms == (2,) for primitive in center_atom_primitives)
    assert {gic.irrep for gic in center_atom_gics} == {"A1'", "A2''"}
    assert all(gic.name.startswith(gic.irrep.replace("'", "p")) for gic in center_atom_gics)
    assert any("(X(1)+X(3)+X(6)+X(5)+X(4))/5" in line for line in gaussian_lines)
    assert any("(X(7)+X(8)+X(9)+X(10)+X(11))/5" in line for line in gaussian_lines)
    assert any("X(2)" in line and "CnAt" in line for line in gaussian_lines)


def test_gicforge_staggered_ferrocene_d5d_projector_closes_b_rank(tmp_path):
    source = _test_molecule_path("ferrocene_staggered.inp")
    xyzin = tmp_path / "ferrocene_staggered.xyzin"

    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    centers = write_interaction_center_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    matrix = build_gic_b_matrix_from_xyzin(xyzin)
    symmetry_section = section_content(xyzin.read_text(encoding="utf-8").splitlines(), "SYMMETRY")
    rows = np.asarray(matrix.rows, dtype=float)
    ring_centers = tuple(center for center in centers.centers if center.kind == "RING_CENTER")
    center_atom_gics = tuple(gic for gic in definition.gics if gic.family == "CENTER_ATOM_DISTANCE")

    assert definition.point_group == "D5d"
    assert any(line == "POINT_GROUP D5d" for line in symmetry_section)
    assert sum(1 for line in symmetry_section if " LABEL=" in line and "PERMUTATION=" in line) == 20
    _assert_symmetry_operations_closed(symmetry_section)
    assert definition.target_rank == 57
    assert definition.rank == 57
    assert definition.symmetry_diagnostics is not None
    assert definition.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert np.linalg.matrix_rank(rows, tol=1.0e-8) == 57
    assert [center.atoms for center in ring_centers] == [(1, 3, 6, 5, 4), (7, 8, 9, 10, 11)]
    assert {gic.irrep for gic in center_atom_gics} == {"A1g", "A2u"}


def test_gicforge_point_group_projector_handles_c3v_degenerate_irrep():
    primitives = tuple(
        GICPrimitive(
            identifier=f"P{idx:03d}",
            name=f"Str{idx:04d}",
            family="STRETCH",
            function="R",
            atoms=(1, idx + 1),
        )
        for idx in range(1, 4)
    )
    definition = GICDefinition(
        backend="test",
        point_group="C3v",
        symmetrize=False,
        target_rank=3,
        rank=3,
        candidate_count=3,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 4,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )
    identity = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    operations = (
        GICPointGroupOperation("E", identity, (1, 2, 3, 4)),
        GICPointGroupOperation("C3z^1", identity, (1, 3, 4, 2)),
        GICPointGroupOperation("C3z^2", identity, (1, 4, 2, 3)),
        GICPointGroupOperation("sigma_v_3_0", identity, (1, 2, 4, 3)),
        GICPointGroupOperation("sigma_v_3_1", identity, (1, 4, 3, 2)),
        GICPointGroupOperation("sigma_v_3_2", identity, (1, 3, 2, 4)),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("N", "H", "H", "H"),
        symmetry_operations=operations,
    )

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert [gic.name for gic in symmetrized.gics] == [
        "A1Str001",
        "EStr001",
        "EStr002",
    ]
    assert [gic.irrep for gic in symmetrized.gics] == ["A1", "E", "E"]
    assert total_symmetric_gic_names(symmetrized) == ("A1Str001",)


def test_gicforge_point_group_projector_handles_dnd_even_group():
    primitives = tuple(
        GICPrimitive(
            identifier=f"P{idx:03d}",
            name=f"Str{idx:04d}",
            family="STRETCH",
            function="R",
            atoms=(1, idx + 1),
        )
        for idx in range(1, 5)
    )
    definition = GICDefinition(
        backend="test",
        point_group="D2d",
        symmetrize=False,
        target_rank=4,
        rank=4,
        candidate_count=4,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 5,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("C", "H", "H", "H", "H"),
        symmetry_operations=_d2d_stretch_operations(),
    )

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert Counter(gic.irrep for gic in symmetrized.gics) == Counter({"A1": 1, "B2": 1, "E": 2})
    assert total_symmetric_gic_names(symmetrized) == ("A1Str001",)


def test_gicforge_point_group_projector_handles_tetrahedral_stretches():
    primitives = tuple(
        GICPrimitive(
            identifier=f"P{idx:03d}",
            name=f"Str{idx:04d}",
            family="STRETCH",
            function="R",
            atoms=(1, idx + 1),
        )
        for idx in range(1, 5)
    )
    definition = GICDefinition(
        backend="test",
        point_group="Td",
        symmetrize=False,
        target_rank=4,
        rank=4,
        candidate_count=4,
        reference_coordinates_angstrom=tuple(map(tuple, _tetrahedral_methane_coordinates())),
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("C", "H", "H", "H", "H"),
        symmetry_operations=_tetrahedral_operations(),
    )

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert [gic.name for gic in symmetrized.gics] == [
        "A1Str001",
        "T2Str001",
        "T2Str002",
        "T2Str003",
    ]
    assert [gic.irrep for gic in symmetrized.gics] == ["A1", "T2", "T2", "T2"]
    assert total_symmetric_gic_names(symmetrized) == ("A1Str001",)


def test_gicforge_projector_symmetrizes_special_fragment_center_atom_coordinates():
    primitives = (
        GICPrimitive(
            "P001",
            "FCAt0001",
            "FRAG_CENTER_ATOM_DISTANCE",
            "FCA_DIST",
            (1, 2),
            ref_atoms=(3,),
        ),
        GICPrimitive(
            "P002",
            "FCAt0002",
            "FRAG_CENTER_ATOM_DISTANCE",
            "FCA_DIST",
            (4, 5),
            ref_atoms=(6,),
        ),
    )
    definition = GICDefinition(
        backend="test",
        point_group="C2",
        symmetrize=False,
        target_rank=2,
        rank=2,
        candidate_count=2,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 6,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )
    operations = (
        GICPointGroupOperation(
            "E",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (1, 2, 3, 4, 5, 6),
        ),
        GICPointGroupOperation(
            "C2z^1",
            ((-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
            (4, 5, 6, 1, 2, 3),
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("O", "H", "H", "O", "H", "H"),
        symmetry_operations=operations,
    )

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert [gic.name for gic in symmetrized.gics] == ["AFCAt001", "BFCAt001"]
    assert [gic.irrep for gic in symmetrized.gics] == ["A", "B"]
    assert symmetrized.symmetry_diagnostics.groups[0].block == ("SPECIAL_FRAGMENT_CENTER_ATOM")
    assert total_symmetric_gic_names(symmetrized) == ("AFCAt001",)


def test_gicforge_projector_symmetrizes_butterfly_without_torsion_mixing():
    primitives = (
        GICPrimitive("P001", "BtFl0001", "BUTTERFLY", "D", (1, 2, 3, 4)),
        GICPrimitive("P002", "BtFl0002", "BUTTERFLY", "D", (5, 6, 7, 8)),
        GICPrimitive("P003", "Tors0001", "TORSION", "D", (1, 3, 4, 2)),
        GICPrimitive("P004", "Tors0002", "TORSION", "D", (5, 7, 8, 6)),
    )
    definition = GICDefinition(
        backend="test",
        point_group="C2",
        symmetrize=False,
        target_rank=4,
        rank=4,
        candidate_count=4,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 8,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )
    operations = (
        GICPointGroupOperation(
            "E",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            tuple(range(1, 9)),
        ),
        GICPointGroupOperation(
            "C2z^1",
            ((-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
            (5, 6, 7, 8, 1, 2, 3, 4),
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("C",) * 8,
        symmetry_operations=operations,
    )

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert [
        (group.block, group.family, group.output_gics)
        for group in symmetrized.symmetry_diagnostics.groups
    ] == [
        ("BUTTERFLY", "BUTTERFLY", ("ABtFl001", "BBtFl001")),
        ("TORSION", "TORSION", ("ATors001", "BTors001")),
    ]
    assert [gic.irrep for gic in symmetrized.gics] == ["A", "B", "A", "B"]
    assert total_symmetric_gic_names(symmetrized) == ("ABtFl001", "ATors001")


def test_gicforge_projector_symmetrizes_special_fragment_orientation_components():
    primitives = (
        GICPrimitive(
            "P001",
            "FRot0001",
            "FRAG_ORIENTATION",
            "FROT",
            (1, 2),
            mode=0,
            ref_atoms=(3, 4),
            frame_atoms=(1, 2),
            ref_frame_atoms=(3, 4),
        ),
        GICPrimitive(
            "P002",
            "FRot0002",
            "FRAG_ORIENTATION",
            "FROT",
            (1, 2),
            mode=2,
            ref_atoms=(3, 4),
            frame_atoms=(1, 2),
            ref_frame_atoms=(3, 4),
        ),
    )
    definition = GICDefinition(
        backend="test",
        point_group="C2",
        symmetrize=False,
        target_rank=2,
        rank=2,
        candidate_count=2,
        reference_coordinates_angstrom=((0.0, 0.0, 0.0),) * 4,
        primitives=primitives,
        gics=tuple(
            FrozenGIC(
                identifier=f"GIC{idx:03d}",
                name=primitive.name,
                family=primitive.family,
                irrep="UNASSIGNED",
                primitive_id=primitive.identifier,
                gaussian_expression=primitive.gaussian_expression(),
                coefficients=((primitive.identifier, 1.0),),
            )
            for idx, primitive in enumerate(primitives, start=1)
        ),
    )
    operations = (
        GICPointGroupOperation(
            "E",
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (1, 2, 3, 4),
        ),
        GICPointGroupOperation(
            "C2z^1",
            ((-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
            (1, 2, 3, 4),
        ),
    )

    symmetrized = symmetrize_gic_definition(
        definition,
        atom_symbols=("C", "H", "C", "H"),
        symmetry_operations=operations,
    )

    assert symmetrized.symmetry_diagnostics is not None
    assert symmetrized.symmetry_diagnostics.method == "POINT_GROUP_PROJECTOR"
    assert symmetrized.symmetry_diagnostics.groups[0].block == "SPECIAL_FRAGMENT_ORIENTATION"
    assert [gic.name for gic in symmetrized.gics] == ["AFRot001", "BFRot001"]
    assert [gic.irrep for gic in symmetrized.gics] == ["A", "B"]
    assert total_symmetric_gic_names(symmetrized) == ("AFRot001",)
