#!/usr/bin/env python3
"""Canonical CLI entrypoint for the endpoint modeling pipeline."""

from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from version1.endpoint_modeling.pipeline import cli_main


if __name__ == "__main__":
    cli_main()
