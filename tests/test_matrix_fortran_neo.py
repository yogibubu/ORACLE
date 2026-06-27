from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from matrix_chem import (
    preprocess_to_enriched_xyz,
    read_enriched_xyz,
    write_validation_section,
)
from matrix_engines import (
    LEGACY_GICFORGE_FILES,
    gicforge_fortran_layout,
    run_legacy_gicforge,
    validate_legacy_gicforge_sources,
)
from matrix_fragments import write_fragment_build_section
from matrix_neo import (
    DEFAULT_FORTRAN_AUDIT_MOLECULES,
    GICForgeFortranAudit,
    GICForgeFortranAuditResult,
    audit_gicforge_fortran_corpus,
    build_gic_b_matrix,
    format_gicforge_fortran_audit_cases,
    format_gicforge_fortran_audit_summary,
    write_gicforge_build_sections,
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


ORIGINAL_MERLINO_GICFORGE = Path("/Users/vincenzobarone/merlino3.0/bin/gicforge.x")
GOLDEN_CORPUS = (
    Path(__file__).resolve().parent / "fixtures" / "golden_corpus" / "neo_gic_golden.json"
)


def _row_space_basis(matrix: np.ndarray) -> tuple[np.ndarray, int]:
    singular_values: np.ndarray
    _u, singular_values, vh = np.linalg.svd(matrix, full_matrices=False)
    tolerance = max(matrix.shape) * float(singular_values[0]) * np.finfo(float).eps * 10.0
    rank = int(np.sum(singular_values > tolerance))
    return vh[:rank], rank


def _row_space_residual(left: np.ndarray, right: np.ndarray) -> float:
    left_basis, left_rank = _row_space_basis(left)
    right_basis, right_rank = _row_space_basis(right)
    if left_rank != right_rank:
        return float("inf")
    projected = left_basis @ right_basis.T @ right_basis
    return float(np.linalg.norm(left_basis - projected) / np.sqrt(float(left_rank)))


def test_legacy_merlino_gicforge_sources_are_vendored():
    missing = validate_legacy_gicforge_sources(Path(__file__).resolve().parents[1])

    assert missing == ()
    assert "dina25.f" in LEGACY_GICFORGE_FILES
    assert "gicprune.f" in LEGACY_GICFORGE_FILES


def test_fortran_audit_report_includes_projector_diagnostics(tmp_path):
    audit = GICForgeFortranAudit(
        root=tmp_path,
        workdir=tmp_path / "audit",
        tolerance=2.0e-8,
        results=(
            GICForgeFortranAuditResult(
                molecule="case.inp",
                source=tmp_path / "case.inp",
                status="PASS",
                oracle_rank=6,
                fortran_rank=6,
                oracle_row_rank=6,
                fortran_row_rank=6,
                row_space_residual=1.0e-12,
                projector_status="POINT_GROUP_PROJECTOR",
                symmetry_group_count=4,
                special_symmetry_group_count=1,
                mixed_symmetry_group_count=0,
                total_symmetric_gic_count=2,
                salc_coefficient_gic_count=3,
                salc_coefficient_max_norm_error=2.0e-16,
            ),
        ),
    )

    summary = "\n".join(format_gicforge_fortran_audit_summary(audit))
    cases = "\n".join(format_gicforge_fortran_audit_cases(audit))

    assert "MIXED_SYMMETRY_GROUPS 0" in summary
    assert "SALC_COEFFICIENT_GICS 3" in summary
    assert "MAX_SALC_COEFFICIENT_NORM_ERROR 2e-16" in summary
    assert "projector_status=POINT_GROUP_PROJECTOR" in cases
    assert "symmetry_groups=4" in cases
    assert "special_symmetry_groups=1" in cases
    assert "total_symmetric_gics=2" in cases
    assert "salc_coefficient_gics=3" in cases
    assert "salc_norm_error=2e-16" in cases


def test_default_fortran_audit_covers_official_golden_parity_roles():
    registry = json.loads(GOLDEN_CORPUS.read_text(encoding="utf-8"))
    by_role = {
        role: {Path(entry["path"]).name for entry in registry["entries"] if role in entry["roles"]}
        for role in {
            "ring",
            "fused_ring",
            "bridged_ring",
            "spiro_ring",
            "python_fortran_parity",
        }
    }
    default = set(DEFAULT_FORTRAN_AUDIT_MOLECULES)

    assert by_role["python_fortran_parity"] <= default
    assert by_role["fused_ring"] <= default
    assert by_role["bridged_ring"] <= default
    assert by_role["spiro_ring"] <= default
    assert {"benzene.inp", "pyrrole.inp"} <= default
    assert {
        "azulene.inp",
        "pyrene.inp",
        "norbornane.inp",
        "norbornene.inp",
        "norbornadiene.inp",
        "norcamphor.inp",
        "thujone.inp",
        "ribose.inp",
        "cubane.inp",
        "spiro.inp",
        "cyclottane.inp",
    } <= default


def test_legacy_merlino_group_dispatch_includes_ih_and_dnd_extensions():
    root = Path(__file__).resolve().parents[1]
    symm = (
        root / "engines" / "fortran" / "gicforge" / "legacy_merlino" / "symm.f"
    ).read_text(encoding="utf-8")

    assert "FAM .EQ. 'Dnd'" in symm
    assert "CALL OPS_DND" in symm
    assert "GROUP .EQ. 'Ih'" in symm
    assert "CALL OPS_IH" in symm
    assert "SUBROUTINE OPS_IH" in symm


def test_legacy_merlino_ring_and_butterfly_blocks_remain_reference():
    root = Path(__file__).resolve().parents[1]
    legacy = root / "engines" / "fortran" / "gicforge" / "legacy_merlino"
    mksalc = (legacy / "mksalc.f").read_text(encoding="utf-8")
    mkcyc = (legacy / "mkcyc.f").read_text(encoding="utf-8")
    gicprune = (legacy / "gicprune.f").read_text(encoding="utf-8")

    assert "Subroutine BtFly" in mksalc
    assert "Subroutine CySalc" in mksalc
    assert "Subroutine PrtPckQP" in mksalc
    assert "SQRT(RPck" in mksalc
    assert "ATAN2(RPck" in mksalc
    assert "Subroutine CycAng" in mkcyc
    assert "Subroutine CyGNSVD" in mkcyc
    assert "Subroutine PruneGICBlocks" in gicprune


def test_python_pseudo_bond_hbond_mode_tracks_legacy_fortran_hbond_contract(tmp_path):
    root = Path(__file__).resolve().parents[1]
    legacy = root / "engines" / "fortran" / "gicforge" / "legacy_merlino"
    mkprim = (legacy / "mkprim.f").read_text(encoding="utf-8")
    coord = (legacy / "coord.f").read_text(encoding="utf-8")
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
    definition = write_gicforge_build_sections(xyzin, fragment_mode="h-bonds")

    assert definition.pseudo_bonds == ((5, 6),)
    assert definition.pseudo_bond_kinds == ("HBOND",)
    assert "Subroutine FindHBnd" in mkprim
    assert "Subroutine MkHBnd" in mkprim
    assert "BDPCS3_HB_ANGLE_MIN" in mkprim
    assert "NBond(JAt)=NBond(JAt)+1" in mkprim
    assert "IBond(NBond(JAt),JAt)=KAt" in mkprim
    assert "IBond(NBond(KAt),KAt)=JAt" in mkprim
    assert "DoHBnd = KWd(17)" in coord


def test_legacy_merlino_gicforge_backend_compiles():
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran is not available")

    root = Path(__file__).resolve().parents[1]
    layout = gicforge_fortran_layout(root)

    result = subprocess.run(
        [str(layout.legacy_compile_script)],
        check=True,
        cwd=root,
        capture_output=True,
        text=True,
    )

    assert layout.legacy_executable.is_file()
    assert str(layout.legacy_executable) in result.stdout


def test_legacy_merlino_benzene_fortran_baseline_is_unchanged(tmp_path):
    root = Path(__file__).resolve().parents[1]
    layout = gicforge_fortran_layout(root)
    if shutil.which("gfortran") is None and not layout.legacy_executable.is_file():
        pytest.skip("gfortran is not available")

    source = _test_molecule_path("benzene.inp")
    xyzin = tmp_path / "benzene.xyzin"
    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    geometry = read_enriched_xyz(xyzin)

    legacy = run_legacy_gicforge(
        tmp_path / "benzene-legacy",
        atoms=geometry.atoms,
        coordinates_angstrom=geometry.coordinates_angstrom,
        point_group="D6h",
        title="benzene",
        keywords=("GNIC", "BMAT"),
        repo_root=root,
    )

    assert "Point Group from symm.f: C2v" in legacy.provout


@pytest.mark.parametrize(
    ("molecule", "expected_rank"),
    [
        ("pyrrole.inp", 24),
        ("benzene.inp", 30),
        ("naphtalene.inp", 48),
        ("phenantrene.inp", 66),
        ("anthracene.inp", 66),
        ("pyrene.inp", 72),
        ("fluorene.inp", 63),
    ],
)
def test_legacy_merlino_executable_bmatrix_span_matches_oracle_corpus(
    tmp_path,
    molecule,
    expected_rank,
):
    root = Path(__file__).resolve().parents[1]
    layout = gicforge_fortran_layout(root)
    if shutil.which("gfortran") is None and not layout.legacy_executable.is_file():
        pytest.skip("gfortran is not available")

    source = _test_molecule_path(molecule)
    xyzin = tmp_path / f"{Path(molecule).stem}.xyzin"
    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin)
    geometry = read_enriched_xyz(xyzin)
    oracle_b_matrix = np.asarray(build_gic_b_matrix(definition).rows, dtype=float)

    legacy = run_legacy_gicforge(
        tmp_path / f"{Path(molecule).stem}-legacy",
        atoms=geometry.atoms,
        coordinates_angstrom=geometry.coordinates_angstrom,
        point_group="C1",
        title=Path(molecule).stem,
        keywords=("GNIC", "BMAT"),
        repo_root=root,
    )
    legacy_b_matrix = np.asarray(legacy.b_matrix_rows, dtype=float)
    legacy_prefixes = {label[:4] for label in legacy.gic_labels}
    oracle_rpck_labels = {
        gic.name for gic in definition.gics if gic.family == "RING_PUCKER_COMPONENT"
    }

    assert definition.target_rank == expected_rank
    assert definition.rank == expected_rank
    assert legacy.final_counts[-1] == expected_rank
    assert legacy_b_matrix.shape == oracle_b_matrix.shape
    assert _row_space_basis(legacy_b_matrix)[1] == expected_rank
    assert _row_space_basis(oracle_b_matrix)[1] == expected_rank
    assert _row_space_residual(oracle_b_matrix, legacy_b_matrix) < 2.0e-8
    fused = molecule in {
        "naphtalene.inp",
        "phenantrene.inp",
        "anthracene.inp",
        "pyrene.inp",
        "fluorene.inp",
    }
    if fused:
        assert {"RPck"}.issubset(legacy_prefixes)
        assert "QPck" not in legacy_prefixes
        assert "PhiP" not in legacy_prefixes
        assert oracle_rpck_labels
    else:
        assert {"RPck", "QPck", "PhiP"}.issubset(legacy_prefixes)
        assert oracle_rpck_labels


@pytest.mark.parametrize(
    ("molecule", "expected_rank"),
    [
        ("azulene.inp", 48),
        ("norbornane.inp", 51),
        ("norbornene.inp", 45),
        ("norbornadiene.inp", 39),
        ("norcamphor.inp", 48),
        ("spiro.inp", 87),
    ],
)
def test_legacy_merlino_executable_bmatrix_span_matches_bridged_ring_probe(
    tmp_path,
    molecule,
    expected_rank,
):
    root = Path(__file__).resolve().parents[1]
    layout = gicforge_fortran_layout(root)
    if shutil.which("gfortran") is None and not layout.legacy_executable.is_file():
        pytest.skip("gfortran is not available")

    source = _test_molecule_path(molecule)
    _skip_if_smiles_requires_rdkit(source)
    xyzin = tmp_path / f"{Path(molecule).stem}.xyzin"
    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    definition = write_gicforge_build_sections(xyzin, symmetrize=True)
    geometry = read_enriched_xyz(xyzin)
    oracle_b_matrix = np.asarray(build_gic_b_matrix(definition).rows, dtype=float)

    legacy = run_legacy_gicforge(
        tmp_path / f"{Path(molecule).stem}-legacy",
        atoms=geometry.atoms,
        coordinates_angstrom=geometry.coordinates_angstrom,
        point_group="C1",
        title=Path(molecule).stem,
        keywords=("GNIC", "BMAT"),
        repo_root=root,
    )
    legacy_b_matrix = np.asarray(legacy.b_matrix_rows, dtype=float)
    legacy_prefixes = {label[:4] for label in legacy.gic_labels}

    assert definition.target_rank == expected_rank
    assert definition.rank == expected_rank
    assert legacy.final_counts[-1] == expected_rank
    assert "Post-pruning GIC count below" not in legacy.provout
    assert "RPck" in legacy_prefixes
    if molecule == "spiro.inp":
        assert "Spir" in legacy_prefixes
        assert any(gic.family == "SPIRO_BEND" for gic in definition.gics)
    assert _row_space_basis(legacy_b_matrix)[1] == expected_rank
    assert _row_space_basis(oracle_b_matrix)[1] == expected_rank
    assert _row_space_residual(oracle_b_matrix, legacy_b_matrix) < 2.0e-7
    if molecule == "azulene.inp":
        assert "BtFl" in legacy_prefixes
        assert "QPck" not in legacy_prefixes
        assert "PhiP" not in legacy_prefixes


@pytest.mark.parametrize(
    ("molecule", "expected_rank"),
    [
        ("c2h2.inp", 7),
        ("c4s.inp", 10),
        ("thujone.inp", 75),
        ("ribose.inp", 51),
        ("cubane.inp", 42),
        ("cyclottane.inp", 66),
    ],
)
def test_fortran_audit_handles_added_gic_regressions(
    tmp_path,
    molecule,
    expected_rank,
):
    root = Path(__file__).resolve().parents[1]
    layout = gicforge_fortran_layout(root)
    if shutil.which("gfortran") is None and not layout.legacy_executable.is_file():
        pytest.skip("gfortran is not available")

    audit = audit_gicforge_fortran_corpus(
        root=_test_molecule_path("benzene.inp").parent,
        molecules=(molecule,),
        workdir=tmp_path / "audit",
        repo_root=root,
    )
    result = audit.results[0]

    assert result.status == "PASS", result.message
    assert result.oracle_rank == expected_rank
    assert result.fortran_rank == expected_rank
    assert result.oracle_row_rank == expected_rank
    assert result.fortran_row_rank == expected_rank


@pytest.mark.skipif(
    not ORIGINAL_MERLINO_GICFORGE.is_file(),
    reason="original merlino3.0 GICForge executable is not available",
)
@pytest.mark.parametrize(
    "molecule",
    [
        "pyrrole.inp",
        "benzene.inp",
        "naphtalene.inp",
        "phenantrene.inp",
        "anthracene.inp",
        "pyrene.inp",
        "fluorene.inp",
    ],
)
def test_vendored_gicforge_matches_original_merlino_executable(tmp_path, molecule):
    root = Path(__file__).resolve().parents[1]
    layout = gicforge_fortran_layout(root)
    if shutil.which("gfortran") is None and not layout.legacy_executable.is_file():
        pytest.skip("gfortran is not available")

    source = _test_molecule_path(molecule)
    xyzin = tmp_path / f"{Path(molecule).stem}.xyzin"
    preprocess_to_enriched_xyz(source, xyzin)
    write_validation_section(xyzin)
    geometry = read_enriched_xyz(xyzin)
    run_kwargs = {
        "atoms": geometry.atoms,
        "coordinates_angstrom": geometry.coordinates_angstrom,
        "point_group": "C1",
        "title": Path(molecule).stem,
        "keywords": ("GNIC", "BMAT"),
    }

    vendored = run_legacy_gicforge(
        tmp_path / f"{Path(molecule).stem}-vendored",
        repo_root=root,
        **run_kwargs,
    )
    original = run_legacy_gicforge(
        tmp_path / f"{Path(molecule).stem}-merlino",
        repo_root=root,
        executable=ORIGINAL_MERLINO_GICFORGE,
        **run_kwargs,
    )

    fused = molecule in {
        "naphtalene.inp",
        "phenantrene.inp",
        "anthracene.inp",
        "pyrene.inp",
        "fluorene.inp",
    }
    vendored_b = np.asarray(vendored.b_matrix_rows, dtype=float)
    original_b = np.asarray(original.b_matrix_rows, dtype=float)
    if fused:
        vendored_prefixes = {label[:4] for label in vendored.gic_labels}
        assert vendored.final_counts[-1] == original.final_counts[-1]
        assert {"RPck"}.issubset(vendored_prefixes)
        assert "QPck" not in vendored_prefixes
        assert "PhiP" not in vendored_prefixes
        assert _row_space_residual(vendored_b, original_b) < 2.0e-8
    else:
        assert vendored.final_counts == original.final_counts
        assert vendored.gic_labels == original.gic_labels
        np.testing.assert_allclose(
            vendored_b,
            original_b,
            rtol=0.0,
            atol=1.0e-12,
        )


def test_merlino_group_builders_cover_ih_and_dnd(tmp_path):
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran is not available")

    root = Path(__file__).resolve().parents[1]
    source = root / "engines" / "fortran" / "gicforge" / "legacy_merlino" / "symm.f"
    driver = tmp_path / "test_groups.f"
    executable = tmp_path / "test_groups"
    driver.write_text(
        """
      Program TGRP
      Double Precision R(3,3,200)
      Integer NOPS
      Call BUILD_GROUP_OPS('I',' ',0,R,NOPS)
      If (NOPS .NE. 60) Stop 1
      Call BUILD_GROUP_OPS('Ih',' ',0,R,NOPS)
      If (NOPS .NE. 120) Stop 2
      Call BUILD_GROUP_OPS('D2d','Dnd',2,R,NOPS)
      If (NOPS .NE. 8) Stop 3
      Call BUILD_GROUP_OPS('D3d','Dnd',3,R,NOPS)
      If (NOPS .NE. 12) Stop 4
      End
""",
        encoding="utf-8",
    )

    subprocess.run(
        [gfortran, "-std=legacy", str(driver), str(source), "-o", str(executable)],
        check=True,
        cwd=root,
        capture_output=True,
        text=True,
    )
    subprocess.run([str(executable)], check=True, capture_output=True, text=True)


def test_fortran_fragment_tric_bmat_compiles_and_runs(tmp_path):
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran is not available")

    root = Path(__file__).resolve().parents[1]
    source = root / "engines" / "fortran" / "gicforge" / "frag_tric_bmat.f"
    driver = tmp_path / "driver.f"
    executable = tmp_path / "driver"
    driver.write_text(
        """
      Program TFRAG
      Implicit Real*8(A-H,O-Z)
      Integer FAT(3),RAT(3),IPROT(3),ISEL(3)
      Dimension C(3,6),B(18),BMAT(3,3),Q(3,3),WORK(3)
      Data FAT /4,5,6/
      Data RAT /1,2,3/
C
      Do 10 I=1,18
       B(I)=0.0D0
   10 Continue
      Do 20 I=1,6
       Do 30 J=1,3
        C(J,I)=0.0D0
   30  Continue
   20 Continue
C
      C(1,1)=0.00D0
      C(2,1)=0.00D0
      C(3,1)=0.00D0
      C(1,2)=0.96D0
      C(2,2)=0.00D0
      C(3,2)=0.00D0
      C(1,3)=-0.24D0
      C(2,3)=0.93D0
      C(3,3)=0.00D0
      C(1,4)=0.00D0
      C(2,4)=0.00D0
      C(3,4)=3.20D0
      C(1,5)=0.96D0
      C(2,5)=0.00D0
      C(3,5)=3.20D0
      C(1,6)=-0.24D0
      C(2,6)=0.93D0
      C(3,6)=3.20D0
C
      Call ORCFTRN(6,3,FAT,3,RAT,1,B)
      If(DAbs(B(10)-1.0D0/3.0D0).gt.1.0D-10) Stop 1
      If(DAbs(B(1)+1.0D0/3.0D0).gt.1.0D-10) Stop 2
C
      Call ORCFCDI(6,3,FAT,3,RAT,C,B,VAL,IFAIL)
      If(IFAIL.ne.0) Stop 3
      If(DAbs(VAL-3.20D0).gt.1.0D-10) Stop 4
      If(DAbs(B(12)-1.0D0/3.0D0).gt.1.0D-10) Stop 5
      If(DAbs(B(3)+1.0D0/3.0D0).gt.1.0D-10) Stop 6
C
      Call ORCFROT(6,3,FAT,3,RAT,5,6,2,3,1,C,B,VAL,IFAIL)
      If(IFAIL.ne.0) Stop 7
      If(DAbs(VAL).gt.1.0D-10) Stop 8
      SUM=0.0D0
      Do 40 I=1,18
       If(B(I).ne.B(I)) Stop 9
       SUM=SUM+DAbs(B(I))
   40 Continue
      If(SUM.le.1.0D-8) Stop 10
C
      Call ORCGSPC('FRAG_DISTANCE',ISPEC)
      If(ISPEC.ne.1) Stop 11
      Call ORCGSPC('CENTER_ATOM_DISTANCE',ISPEC)
      If(ISPEC.ne.1) Stop 12
      Call ORCGSPC('STRETCH',ISPEC)
      If(ISPEC.ne.0) Stop 13
C
      Do 50 I=1,3
       IPROT(I)=0
       ISEL(I)=0
       WORK(I)=0.0D0
       Do 60 J=1,3
        BMAT(J,I)=0.0D0
        Q(J,I)=0.0D0
   60  Continue
   50 Continue
      BMAT(1,1)=1.0D0
      BMAT(1,2)=1.0D0
      BMAT(2,3)=1.0D0
      IPROT(2)=1
      Call ORCGSEL(3,3,BMAT,IPROT,1,1.0D-7,ISEL,NSEL,IRANK,
     $  Q,WORK,IFAIL)
      If(IFAIL.ne.0) Stop 14
      If(NSEL.ne.1) Stop 15
      If(IRANK.ne.1) Stop 16
      If(ISEL(1).ne.2) Stop 17
C
      IPROT(1)=1
      IPROT(2)=1
      IPROT(3)=0
      BMAT(1,1)=1.0D0
      BMAT(2,1)=0.0D0
      BMAT(3,1)=0.0D0
      BMAT(1,2)=0.0D0
      BMAT(2,2)=1.0D0
      BMAT(3,2)=0.0D0
      BMAT(1,3)=0.0D0
      BMAT(2,3)=0.0D0
      BMAT(3,3)=1.0D0
      Call ORCGSEL(3,3,BMAT,IPROT,1,1.0D-7,ISEL,NSEL,IRANK,
     $  Q,WORK,IFAIL)
      If(IFAIL.ne.2) Stop 18
      End
""",
        encoding="ascii",
    )

    subprocess.run(
        [
            gfortran,
            "-std=legacy",
            "-Wall",
            "-Wextra",
            "-fcheck=all",
            str(source),
            str(driver),
            "-o",
            str(executable),
        ],
        check=True,
        cwd=root,
    )
    subprocess.run([str(executable)], check=True, cwd=root)
