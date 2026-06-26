from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _add_repo_packages_to_path() -> None:
    for src in sorted((REPO_ROOT / "packages").glob("*/src")):
        text = str(src)
        if text not in sys.path:
            sys.path.insert(0, text)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="oracle", description="ORACLE workflow CLI")
    sub = parser.add_subparsers(dest="command")

    init = sub.add_parser("init", help="Create an ORACLE project workspace")
    init.add_argument("workdir", type=Path)

    lcb25 = sub.add_parser("lcb25", help="Manage the local ORACLE LCB25 geometry cache")
    lcb25_sub = lcb25.add_subparsers(dest="lcb25_command")
    fetch = lcb25_sub.add_parser("fetch", help="Download/extract LCB25 geometries once")
    fetch.add_argument("--root", type=Path, default=REPO_ROOT / "data" / "lcb25")
    fetch.add_argument("--dataset", action="append", help="PCS2, SE or HPCS2; repeatable")
    fetch.add_argument("--force", action="store_true")

    sub.add_parser("merlino", help="Delegate to the current Merlino CLI during migration")
    return parser


def main(argv: list[str] | None = None) -> int:
    _add_repo_packages_to_path()
    args, remainder = build_parser().parse_known_args(argv)
    if args.command == "init":
        from oracle_core.workspace import ensure_workspace

        ensure_workspace(args.workdir)
        print(f"Created ORACLE workspace: {args.workdir}")
        return 0
    if args.command == "lcb25" and args.lcb25_command == "fetch":
        from oracle_babel import sync_lcb25_library

        manifest = sync_lcb25_library(args.root, datasets=args.dataset, force=args.force)
        print(f"Synced LCB25 library: {manifest}")
        return 0
    if args.command == "merlino":
        sys.argv = ["merlino", *remainder]
        runpy.run_module("merlino", run_name="__main__")
        return 0
    build_parser().print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
