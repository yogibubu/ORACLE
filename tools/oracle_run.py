from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _add_repo_packages_to_path() -> None:
    for src in sorted((REPO_ROOT / "packages").glob("*/src")):
        text = str(src)
        if text not in sys.path:
            sys.path.insert(0, text)


def main(argv: list[str] | None = None) -> int:
    _add_repo_packages_to_path()
    from oracle_core.cli import main as oracle_main

    return oracle_main(argv, repo_root=REPO_ROOT)


if __name__ == "__main__":
    raise SystemExit(main())
