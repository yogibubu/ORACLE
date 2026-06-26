from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

for src in sorted((ROOT / "packages").glob("*/src")):
    text = str(src)
    if text not in sys.path:
        sys.path.insert(0, text)

from oracle_core.cli import main


if __name__ == "__main__":
    raise SystemExit(main(repo_root=ROOT))
