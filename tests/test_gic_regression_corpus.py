from __future__ import annotations

from pathlib import Path


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


def test_oracle_environment_helpers_define_merlino_style_commands():
    text = ENV_HELPERS.read_text(encoding="utf-8")

    for name in (
        "oracle-set()",
        "oracle-run()",
        "oracle-run-bg()",
        "oracle-run-check()",
        "oracle-test-all()",
        "oracle-clean()",
    ):
        assert name in text
