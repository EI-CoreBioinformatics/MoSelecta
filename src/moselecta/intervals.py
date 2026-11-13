"""Generic interval and precision/recall helper functions."""

from __future__ import annotations

from bisect import bisect_left
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Callable, Generic, Hashable, Iterable, Optional, TypeVar

from .scoring import fraction_or_none


T = TypeVar("T")
K = TypeVar("K", bound=Hashable)


@dataclass(frozen=True, slots=True)
class _IntervalIndexEntry(Generic[T]):
    start0: int
    end0: int
    order: int
    item: T


@dataclass(frozen=True, slots=True)
class _IntervalIndexBucket(Generic[T]):
    entries: tuple[_IntervalIndexEntry[T], ...]
    starts: tuple[int, ...]
    prefix_max_ends: tuple[int, ...]


class IntervalOverlapIndex(Generic[K, T]):
    """Small in-memory interval index for immutable genomic records.

    Records are grouped by an arbitrary hashable key, commonly chromosome or
    ``(chromosome, strand)``. Queries return overlapping records in their
    original input order, preserving the deterministic behaviour of the
    previous linear scans while avoiding comparisons with distant loci.
    """

    def __init__(
        self,
        records: Iterable[T],
        *,
        key: Callable[[T], K],
        start: Callable[[T], int],
        end: Callable[[T], int],
    ) -> None:
        grouped: dict[K, list[_IntervalIndexEntry[T]]] = defaultdict(list)
        for order, item in enumerate(records):
            start0 = int(start(item))
            end0 = int(end(item))
            if end0 <= start0:
                continue
            grouped[key(item)].append(_IntervalIndexEntry(start0, end0, order, item))

        buckets: dict[K, _IntervalIndexBucket[T]] = {}
        for group_key, entries in grouped.items():
            entries.sort(key=lambda entry: (entry.start0, entry.end0, entry.order))
            prefix_max_ends: list[int] = []
            current_max = -1
            for entry in entries:
                current_max = max(current_max, entry.end0)
                prefix_max_ends.append(current_max)
            buckets[group_key] = _IntervalIndexBucket(
                entries=tuple(entries),
                starts=tuple(entry.start0 for entry in entries),
                prefix_max_ends=tuple(prefix_max_ends),
            )
        self._buckets = buckets

    def query(self, group_key: K, start0: int, end0: int) -> list[T]:
        """Return records overlapping ``[start0, end0)`` in input order."""
        if end0 <= start0:
            return []
        bucket = self._buckets.get(group_key)
        if bucket is None:
            return []

        # All possible hits begin before the query end. Walk backwards only
        # while an earlier interval can still extend beyond the query start.
        idx = bisect_left(bucket.starts, int(end0)) - 1
        hits: list[_IntervalIndexEntry[T]] = []
        while idx >= 0 and bucket.prefix_max_ends[idx] > start0:
            entry = bucket.entries[idx]
            if entry.end0 > start0:
                hits.append(entry)
            idx -= 1
        hits.sort(key=lambda entry: entry.order)
        return [entry.item for entry in hits]


def harmonic_mean(precision: Optional[float], sensitivity: Optional[float]) -> Optional[float]:
    if precision is None or sensitivity is None:
        return None
    if precision + sensitivity <= 0:
        return 0.0
    return 2.0 * precision * sensitivity / (precision + sensitivity)


def prf_from_counts(query_supported: int, query_total: int, target_supported: int, target_total: int) -> dict[str, Optional[float]]:
    precision = fraction_or_none(query_supported, query_total)
    sensitivity = fraction_or_none(target_supported, target_total)
    return {
        "precision": precision,
        "sensitivity": sensitivity,
        "f1": harmonic_mean(precision, sensitivity),
    }


def exact_feature_prf(query_features: list[tuple] | set[tuple], target_features: list[tuple] | set[tuple]) -> dict[str, Optional[float]]:
    query_set = set(query_features)
    target_set = set(target_features)
    shared = len(query_set & target_set)
    return prf_from_counts(shared, len(query_set), shared, len(target_set))


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    cleaned = sorted((int(s), int(e)) for s, e in intervals if int(e) > int(s))
    if not cleaned:
        return []
    merged: list[tuple[int, int]] = []
    cur_s, cur_e = cleaned[0]
    for s, e in cleaned[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged


def interval_total_bp(intervals: list[tuple[int, int]]) -> int:
    return sum(e - s for s, e in merge_intervals(intervals))


def interval_list_overlap_bp(a_intervals: list[tuple[int, int]], b_intervals: list[tuple[int, int]]) -> int:
    a = merge_intervals(a_intervals)
    b = merge_intervals(b_intervals)
    i = j = 0
    overlap = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if e > s:
            overlap += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return overlap


def interval_list_overlap_bp_sorted(
    a_intervals: list[tuple[int, int]] | tuple[tuple[int, int], ...],
    b_intervals: list[tuple[int, int]] | tuple[tuple[int, int], ...],
) -> int:
    """Return overlap for already sorted, internally non-overlapping intervals.

    This avoids the sorting and allocation performed by
    :func:`interval_list_overlap_bp`. It is intended for immutable CDS/exon
    chains that are normalised during parsing and collapsing.
    """
    i = j = 0
    overlap = 0
    while i < len(a_intervals) and j < len(b_intervals):
        a_start, a_end = a_intervals[i]
        b_start, b_end = b_intervals[j]
        start0 = max(a_start, b_start)
        end0 = min(a_end, b_end)
        if end0 > start0:
            overlap += end0 - start0
        if a_end < b_end:
            i += 1
        else:
            j += 1
    return overlap


def interval_prf_sorted(
    query_intervals: list[tuple[int, int]] | tuple[tuple[int, int], ...],
    target_intervals: list[tuple[int, int]] | tuple[tuple[int, int], ...],
) -> dict[str, Optional[float]]:
    """Precision/recall/F1 for sorted, internally non-overlapping intervals."""
    query_total = sum(end0 - start0 for start0, end0 in query_intervals)
    target_total = sum(end0 - start0 for start0, end0 in target_intervals)
    shared = interval_list_overlap_bp_sorted(query_intervals, target_intervals)
    return prf_from_counts(shared, query_total, shared, target_total)


def interval_prf(query_intervals: list[tuple[int, int]], target_intervals: list[tuple[int, int]]) -> dict[str, Optional[float]]:
    q_total = interval_total_bp(query_intervals)
    t_total = interval_total_bp(target_intervals)
    shared = interval_list_overlap_bp(query_intervals, target_intervals)
    return prf_from_counts(shared, q_total, shared, t_total)


def subtract_intervals(query_intervals: list[tuple[int, int]], subtracting_intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return parts of query_intervals not covered by subtracting_intervals."""
    query = merge_intervals(query_intervals)
    subtracting = merge_intervals(subtracting_intervals)
    if not query:
        return []
    if not subtracting:
        return query
    out: list[tuple[int, int]] = []
    j = 0
    for qs, qe in query:
        cur = qs
        while j < len(subtracting) and subtracting[j][1] <= qs:
            j += 1
        k = j
        while k < len(subtracting) and subtracting[k][0] < qe:
            ss, se = subtracting[k]
            if ss > cur:
                out.append((cur, min(ss, qe)))
            cur = max(cur, se)
            if cur >= qe:
                break
            k += 1
        if cur < qe:
            out.append((cur, qe))
    return merge_intervals(out)


def classify_intervals_relative_to_span(intervals: list[tuple[int, int]], span_start: int, span_end: int) -> tuple[int, int, int, int]:
    """Classify interval bases as before, inside, or after a reference span.

    Returns (before_bp, inside_bp, after_bp, total_bp). Intervals that straddle
    a boundary are split by base count.
    """
    before = inside = after = 0
    for s, e in merge_intervals(intervals):
        if e <= s:
            continue
        if s < span_start:
            before += max(0, min(e, span_start) - s)
        inside_s = max(s, span_start)
        inside_e = min(e, span_end)
        if inside_e > inside_s:
            inside += inside_e - inside_s
        if e > span_end:
            after += e - max(s, span_end)
    return before, inside, after, before + inside + after


def five_three_from_before_after(strand: str, before_bp: int, after_bp: int) -> tuple[int, int]:
    """Map genomic before/after-span bases to biological 5'/3' side."""
    if strand == "+":
        return before_bp, after_bp
    return after_bp, before_bp


def consensus_intervals_from_source_intervals(source_intervals: dict[str, list[tuple[int, int]]], min_source_support: int = 2) -> list[tuple[int, int]]:
    """Return intervals covered by at least min_source_support distinct sources.

    Source-level interval counts are maintained, so multiple overlapping models
    from the same source count once for that source at each position.
    """
    events: list[tuple[int, int, str]] = []
    for source, intervals in source_intervals.items():
        for s, e in merge_intervals(intervals):
            if e <= s:
                continue
            events.append((s, 1, source))
            events.append((e, -1, source))
    if not events:
        return []
    events.sort(key=lambda x: (x[0], -x[1]))
    active: Counter[str] = Counter()
    out: list[tuple[int, int]] = []
    prev = events[0][0]
    idx = 0
    while idx < len(events):
        pos = events[idx][0]
        if pos > prev and sum(1 for count in active.values() if count > 0) >= min_source_support:
            out.append((prev, pos))
        while idx < len(events) and events[idx][0] == pos:
            _pos, delta, source = events[idx]
            active[source] += delta
            if active[source] <= 0:
                active.pop(source, None)
            idx += 1
        prev = pos
    return merge_intervals(out)
