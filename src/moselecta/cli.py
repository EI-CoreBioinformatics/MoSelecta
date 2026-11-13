"""Command-line entry point for MoSelecta."""

from __future__ import annotations

import sys
from typing import Optional

from .pipeline import main as _pipeline_main


def main(argv: Optional[list[str]] = None) -> int:
    """Run the CLI and present expected filesystem failures concisely."""

    try:
        return _pipeline_main(argv)
    except OSError as exc:
        print(f"MoSelecta input/output error: {exc}", file=sys.stderr)
        return 2

__all__ = ["main"]
