"""Backward-compatible import location for the MoSelecta CLI pipeline.

The implementation was moved to :mod:`moselecta.pipeline` in v0.1.8 so that
``core`` no longer acts as the main orchestration module.  This wrapper is kept
for compatibility with older imports such as ``from moselecta.core import main``.
"""

from __future__ import annotations

from .pipeline import main

__all__ = ["main"]
