from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from oracle_chem import (
    preprocess_to_enriched_xyz,
    read_enriched_xyz,
    write_validation_section,
)
from oracle_engines import (
    LEGACY_GICFORGE_FILES,
    gicforge_fortran_layout,
    run_legacy_gicforge,
    validate_legacy_gicforge_sources,
)
from oracle_gicforge import (
    GICForgeFortranAudit,
    GICForgeFortranAuditResult,
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
            ),
        ),
    )

    summary = "\n".join(format_gicforge_fortran_audit_summary(audit))
    cases = "\n".join(format_gicforge_fortran_audit_cases(audit))

    assert "MIXED_SYMMETRY_GROUPS 0" in summary
    assert "projector_status=POINT_GROUP_PROJECTOR" in cases
    assert "symmetry_groups=4" in cases
    assert "special_symmetry_groups=1" in cases
    assert "total_symmetric_gics=2" in cases


def test_legacy_merlino_group_dispatch_matches_classifier_families():
    root = Path(__file__).resolve().parents[1]
    symm = (
        root / "engines" / "fortran" / "gicforge" / "legacy_merlino" / "symm.f"
    ).read_text(encoding="utf-8")

    assert "FAM .EQ. 'Cnv'" in symm
    assert "FAM .EQ. 'Cnh'" in symm
    assert "FAM .EQ. 'Dn'" in symm
    assert "GROUP .EQ. 'td'" in symm
    assert "GROUP .EQ. 'oh'" in symm
    assert "GROUP .EQ. 'ih'" in symm
    assert "CALL OPS_IH" in symm


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


@pytest.mark.parametrize(
    ("molecule", "expected_rank"),
    [
        ("naphtalene.inp", 48),
        ("phenantrene.inp", 66),
        ("pyrene.inp", 72),
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
    if molecule == "pyrene.inp":
        assert "Dihe" in legacy_prefixes
    else:
        assert {"RPck", "QPck", "PhiP"}.issubset(legacy_prefixes)
    assert oracle_rpck_labels


def test_legacy_merlino_icosahedral_group_builder_runs(tmp_path):
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran is not available")

    root = Path(__file__).resolve().parents[1]
    source = root / "engines" / "fortran" / "gicforge" / "legacy_merlino" / "symm.f"
    driver = tmp_path / "test_ih.f"
    executable = tmp_path / "test_ih"
    driver.write_text(
        """
      Program TIH
      Double Precision R(3,3,200)
      Integer NOPS
      Call OPS_IH(R,NOPS)
      If (NOPS .NE. 120) Stop 1
      Call OPS_I(R,NOPS)
      If (NOPS .NE. 60) Stop 2
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
