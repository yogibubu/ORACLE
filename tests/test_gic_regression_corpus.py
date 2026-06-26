from __future__ import annotations

from pathlib import Path

from oracle_babel import rdkit_available
from oracle_gicforge import (
    audit_gic_corpus_geometry,
    discover_gic_corpus,
    summarize_gic_corpus,
)


CORPUS = Path(__file__).resolve().parent / "fixtures" / "test_molecules" / "molecules"
ENV_HELPERS = Path(__file__).resolve().parents[1] / "scripts" / "oracle_env.sh"


def test_gic_regression_corpus_is_available():
    inputs = sorted(CORPUS.glob("*.inp"))
    names = {path.name for path in inputs}

    assert len(inputs) >= 100
    assert {
        "benzene.inp",
        "cyclopentadiene_B3.inp",
        "nitrobenzene.inp",
        "azulene.inp",
        "norbornane.inp",
    } <= names


def test_gic_regression_corpus_keeps_qm_adapter_outputs():
    assert (CORPUS / "apinene.log").is_file()
    assert (CORPUS / "apinene.out").is_file()
    assert (CORPUS / "c6h5.fchk").is_file()
    assert (CORPUS / "c6h5.gjf").is_file()


def test_gic_regression_corpus_inventory_classifies_files():
    summary = summarize_gic_corpus(CORPUS)
    inp_entries = discover_gic_corpus(CORPUS, suffixes=["inp"])

    assert summary.total_files == 153
    assert summary.suffix_counts[".inp"] == 126
    assert summary.role_counts["legacy_gic_input"] == 126
    assert len(inp_entries) == 126
    assert {entry.role for entry in inp_entries} == {"legacy_gic_input"}


def test_gic_regression_corpus_geometry_audit_tracks_parser_budget():
    audit = audit_gic_corpus_geometry(CORPUS)
    failures = {entry.name for entry in audit.entries if not entry.passed}

    assert audit.total_files == 129
    assert audit.passed_files >= 114
    assert audit.failed_files <= 15
    assert audit.source_format_counts["gaussian_cartesian_input"] == 99
    if rdkit_available():
        assert audit.source_format_counts["gaussian_zmatrix_input"] <= 15
    else:
        assert audit.source_format_counts["gaussian_zmatrix_input"] == 15
    if rdkit_available():
        assert "testvib.inp" in failures
    else:
        assert audit.passed_files == 114
        assert audit.failed_files == 15
        assert audit.error_counts == {"GeometryParseError": 1, "RDKitUnavailableError": 14}
        assert {"azulene.inp", "pyrrole_smile1.inp", "testvib.inp"} <= failures


def test_oracle_environment_helpers_define_oracle_style_commands():
    text = ENV_HELPERS.read_text(encoding="utf-8")

    for name in (
        "oracle-set()",
        "oracle-run()",
        "oracle-run-bg()",
        "oracle-run-check()",
        "oracle-test-all()",
        "oracle-clean()",
        "oracle-create-venv()",
        "oracle-install-runtime-deps()",
        "oracle-gic-corpus-list()",
        "oracle-gic-corpus-summary()",
        "oracle-gic-corpus-audit()",
    ):
        assert name in text
    assert "ORACLE_AUTO_CREATE_VENV" in text
    assert "ORACLE_AUTO_INSTALL_RUNTIME_DEPS" in text
    assert "rdkit" in text
