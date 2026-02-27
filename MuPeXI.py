#!/usr/bin/env python3
"""Backwards-compatible entrypoint.

Keeps the historical MuPeXI.py invocation style, but routes into the packaged code.
"""
from mupexi2.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
