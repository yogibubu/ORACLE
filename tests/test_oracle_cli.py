from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from tools import oracle_run


def test_validate_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyz"

    def fake_validate(target, *, require_fragments=False):
        calls["target"] = target
        calls["require_fragments"] = require_fragments
        return SimpleNamespace(status="PASS")

    monkeypatch.setattr("oracle_chem.write_validation_section", fake_validate)

    rc = oracle_run.main(["validate", str(path), "--require-fragments"])

    assert rc == 0
    assert calls == {"target": path, "require_fragments": True}
    assert "(PASS)" in capsys.readouterr().out


def test_lcb25_fetch_cli_calls_sync(tmp_path, monkeypatch, capsys):
    calls = {}

    def fake_sync(root, *, datasets=None, force=False):
        calls["root"] = root
        calls["datasets"] = datasets
        calls["force"] = force
        return Path(root) / "manifest.json"

    monkeypatch.setattr("oracle_babel.sync_lcb25_library", fake_sync)

    rc = oracle_run.main(
        ["lcb25", "fetch", "--root", str(tmp_path / "cache"), "--dataset", "se", "--force"]
    )

    assert rc == 0
    assert calls == {"root": tmp_path / "cache", "datasets": ["se"], "force": True}
    assert "Synced LCB25 library" in capsys.readouterr().out


def test_fragments_plan_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyz"

    def fake_write(target):
        calls["target"] = target

    monkeypatch.setattr("oracle_fragments.write_fragment_plan_section", fake_write)

    rc = oracle_run.main(["fragments", "plan", str(path)])

    assert rc == 0
    assert calls == {"target": path}
    assert "Planned ORACLE fragment workflow" in capsys.readouterr().out


def test_gicforge_plan_cli_calls_writer(tmp_path, monkeypatch, capsys):
    calls = {}
    path = tmp_path / "molecule.xyz"

    def fake_write(target, *, symmetrize=False, sycart=False):
        calls["target"] = target
        calls["symmetrize"] = symmetrize
        calls["sycart"] = sycart

    monkeypatch.setattr("oracle_gicforge.write_gicforge_plan_sections", fake_write)

    rc = oracle_run.main(["gicforge", "plan", str(path), "--symmetrize", "--sycart"])

    assert rc == 0
    assert calls == {"target": path, "symmetrize": True, "sycart": True}
    assert "Planned GICForge workflow" in capsys.readouterr().out
