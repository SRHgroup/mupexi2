"""Command line interface wrapper for MuPeXI2.

This module intentionally keeps the legacy CLI behaviour from the original MuPeXI script
(getopt-based args). It simply routes argv to the core implementation.
"""
from __future__ import annotations

import sys
from .core import main as core_main


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    core_main(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
