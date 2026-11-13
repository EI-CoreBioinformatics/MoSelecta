"""Progress logging and timing helpers."""

from __future__ import annotations

import os
import sys
import time
from functools import lru_cache
from typing import Optional

try:
    import resource
except ImportError:  # pragma: no cover - non-POSIX fallback
    resource = None

from .constants import PROGRESS_LEVELS


def format_duration(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, sec = divmod(int(round(seconds)), 60)
    if minutes < 60:
        return f"{minutes}m {sec:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes:02d}m {sec:02d}s"


class RunTimer:
    """Records coarse stage timings for the end-of-run log."""

    def __init__(self) -> None:
        self.start_time = time.monotonic()
        self.stage_start_time = self.start_time
        self.rows: list[tuple[str, float, float, float, float]] = []

    def mark(self, stage_name: str) -> None:
        now = time.monotonic()
        duration = now - self.stage_start_time
        cumulative = now - self.start_time
        self.rows.append((stage_name, duration, cumulative, current_rss_mb(), peak_rss_mb()))
        self.stage_start_time = now

    def total_elapsed(self) -> float:
        return time.monotonic() - self.start_time


def peak_rss_mb() -> float:
    """Return process peak resident memory in MiB when supported."""
    if resource is None:
        return _linux_current_rss_mb()
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    # Linux reports KiB; macOS/BSD commonly report bytes.
    if sys.platform == "darwin":
        return value / (1024.0 * 1024.0)
    return value / 1024.0


def _linux_current_rss_mb() -> float:
    """Return current resident memory from Linux ``/proc`` or zero."""
    try:
        with open("/proc/self/statm", "r", encoding="ascii") as handle:
            fields = handle.read().split()
        resident_pages = int(fields[1])
        return resident_pages * os.sysconf("SC_PAGE_SIZE") / (1024.0 * 1024.0)
    except (OSError, ValueError, IndexError):
        return 0.0


@lru_cache(maxsize=1)
def _darwin_proc_pidinfo():
    """Load the native current-RSS query only on macOS."""
    import ctypes

    function = ctypes.CDLL("/usr/lib/libproc.dylib").proc_pidinfo
    function.argtypes = (
        ctypes.c_int, ctypes.c_int, ctypes.c_uint64,
        ctypes.c_void_p, ctypes.c_int,
    )
    function.restype = ctypes.c_int
    return function


def _darwin_current_rss_mb() -> float:
    """Return proc_taskinfo.pti_resident_size, not the monotonic peak."""
    try:
        import ctypes

        # PROC_PIDTASKINFO (4) begins with virtual_size and resident_size.
        # The kernel requires room for the complete proc_taskinfo structure
        # (96 bytes on current macOS), although only its first two fields are
        # needed here.
        fields = (ctypes.c_uint64 * 16)()
        size = ctypes.sizeof(fields)
        if _darwin_proc_pidinfo()(os.getpid(), 4, 0, fields, size) < 16:
            return 0.0
        return fields[1] / (1024.0 * 1024.0)
    except (OSError, ValueError):
        return 0.0


def current_rss_mb() -> float:
    """Return current process resident memory in MiB where supported.

    Current RSS may not fall immediately after objects are released because
    Python's allocator can retain arenas for reuse. Fall back to peak RSS only
    on platforms without a native current-RSS query.
    """
    current = (
        _darwin_current_rss_mb() if sys.platform == "darwin"
        else _linux_current_rss_mb()
    )
    return current if current > 0 else peak_rss_mb()


class ProgressLogger:
    """Small stderr progress logger for long genome-scale runs."""

    def __init__(self, level: str = "basic", interval: int = 10000) -> None:
        if level not in PROGRESS_LEVELS:
            level = "basic"
        self.level_name = level
        self.level = PROGRESS_LEVELS[level]
        self.interval = max(1, int(interval))
        self.start_time = time.monotonic()
        self.last_time = self.start_time

    def enabled_for(self, level: str) -> bool:
        return self.level >= PROGRESS_LEVELS.get(level, 1)

    def log(self, message: str, *, level: str = "basic") -> None:
        if not self.enabled_for(level):
            return
        now = time.monotonic()
        elapsed = now - self.start_time
        delta = now - self.last_time
        self.last_time = now
        print(f"[progress {elapsed:,.1f}s +{delta:,.1f}s] {message}", file=sys.stderr, flush=True)

    def every(self, count: int, total: Optional[int], message: str, *, level: str = "basic", interval: Optional[int] = None) -> None:
        if not self.enabled_for(level):
            return
        step = max(1, int(interval or self.interval))
        if count <= 0 or (count % step != 0 and (total is None or count != total)):
            return
        if total:
            pct = 100.0 * count / total
            self.log(f"{message}: {count:,}/{total:,} ({pct:.1f}%)", level=level)
        else:
            self.log(f"{message}: {count:,}", level=level)
