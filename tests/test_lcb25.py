from __future__ import annotations

import json
import zipfile

import pytest

from oracle_babel import (
    LCB25_DATASETS,
    extract_lcb25_archive,
    lcb25_dataset_url,
    lcb25_download_plan,
    sync_lcb25_library,
)


def test_lcb25_download_plan_uses_official_zip_names():
    plan = lcb25_download_plan()

    assert tuple(item.label for item in plan) == LCB25_DATASETS
    assert lcb25_dataset_url("PCS2").endswith("/PCS2.zip")
    assert lcb25_dataset_url("SE").endswith("/SE.zip")
    assert lcb25_dataset_url("HPCS2").endswith("/HPCS2.zip")


def test_lcb25_archive_extraction_returns_xyz_files(tmp_path):
    archive = tmp_path / "SE.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("water.xyz", "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n")
        zf.writestr("notes.txt", "ignore")

    extracted = extract_lcb25_archive(archive, tmp_path / "out")

    assert len(extracted) == 1
    assert extracted[0].name == "water.xyz"


def test_lcb25_archive_extraction_rejects_path_traversal(tmp_path):
    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("../escape.xyz", "1\nbad\nH 0 0 0\n")

    with pytest.raises(ValueError, match="unsafe path"):
        extract_lcb25_archive(archive, tmp_path / "out")


def test_lcb25_sync_writes_manifest_without_network(tmp_path, monkeypatch):
    def fake_download(label, target_dir):
        archive = target_dir / f"{label}.zip"
        with zipfile.ZipFile(archive, "w") as zf:
            zf.writestr(f"{label.lower()}/water.xyz", "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n")
        return archive

    monkeypatch.setattr("oracle_babel.lcb25.download_lcb25_dataset", fake_download)

    manifest_path = sync_lcb25_library(tmp_path / "lcb25", datasets=("se",))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert manifest["schema"] == "oracle.lcb25.cache.v1"
    assert manifest["datasets"][0]["label"] == "SE"
    assert manifest["datasets"][0]["xyz_count"] == 1
    assert len(manifest["datasets"][0]["archive_sha256"]) == 64
    assert manifest["datasets"][0]["xyz_files"] == ["xyz/SE/se/water.xyz"]
