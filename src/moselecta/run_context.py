"""Runtime context helpers for MoSelecta pipeline execution.

This module intentionally keeps only lightweight run state: parsed CLI args,
progress/timing helpers, enabled output keys, and the list of outputs written.
It is a refactor-only layer that makes future evidence/genome/parallel pipeline
stages easier to add without passing the same bookkeeping objects everywhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .config import expand_name_list, resolve_enabled_outputs
from .log_utils import ProgressLogger, RunTimer


@dataclass
class RunContext:
    """Mutable state shared by one MoSelecta command-line run."""

    args: Any
    timer: RunTimer
    progress: ProgressLogger
    enabled_outputs: set[str]
    explicit_skips: set[str]
    skip_all_outputs: bool
    written_output_items: list[tuple[str, Path]] = field(default_factory=list)

    @classmethod
    def from_args(cls, args: Any) -> "RunContext":
        """Create a context after pruning/output presets have been resolved."""
        timer = RunTimer()
        progress = ProgressLogger(args.progress_log_level, args.progress_interval)
        explicit_skips = expand_name_list(args.skip_output)
        enabled_outputs = resolve_enabled_outputs(args.output_profile, args.write_output, args.skip_output)
        return cls(
            args=args,
            timer=timer,
            progress=progress,
            enabled_outputs=enabled_outputs,
            explicit_skips=explicit_skips,
            skip_all_outputs="all" in explicit_skips,
        )

    def output_enabled(self, key: str) -> bool:
        """Return true when an output key should be written for this run."""
        return key in self.enabled_outputs

    def ensure_output(self, key: str) -> None:
        """Enable an output key unless it was explicitly skipped by the user."""
        if not self.skip_all_outputs and key not in self.explicit_skips:
            self.enabled_outputs.add(key)

    def note_output(self, key: str, label: str, path: Path) -> None:
        """Record a written output if its output key is enabled."""
        if self.output_enabled(key):
            self.written_output_items.append((label, path))
