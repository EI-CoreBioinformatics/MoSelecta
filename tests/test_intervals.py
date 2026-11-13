from __future__ import annotations

from moselecta.intervals import (
    classify_intervals_relative_to_span,
    interval_list_overlap_bp,
    interval_prf,
    merge_intervals,
    subtract_intervals,
)


def test_merge_intervals_sorts_merges_and_drops_empty_intervals() -> None:
    assert merge_intervals([(10, 20), (5, 7), (7, 10), (20, 20), (18, 25)]) == [(5, 25)]


def test_interval_overlap_uses_merged_half_open_intervals() -> None:
    assert interval_list_overlap_bp([(0, 10), (20, 30)], [(5, 25)]) == 10


def test_interval_prf_reports_base_precision_sensitivity_f1() -> None:
    prf = interval_prf([(0, 10)], [(5, 15)])
    assert prf["precision"] == 0.5
    assert prf["sensitivity"] == 0.5
    assert prf["f1"] == 0.5


def test_subtract_intervals_splits_query_by_subtracting_coverage() -> None:
    assert subtract_intervals([(0, 20)], [(3, 5), (10, 15)]) == [(0, 3), (5, 10), (15, 20)]


def test_classify_intervals_relative_to_span_splits_boundary_crossing_bases() -> None:
    assert classify_intervals_relative_to_span([(0, 12), (18, 25)], 10, 20) == (10, 4, 5, 19)
