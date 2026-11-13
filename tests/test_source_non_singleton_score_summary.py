from __future__ import annotations

import csv
from pathlib import Path

import pytest

from moselecta.console_reporting import (
    summarise_score_rows_by_source,
    summarise_source_scaled_scores,
)
from moselecta.structural import write_model_selection_score_summary_tsv


def score_row(
    source: str,
    score: float,
    sublocus_model_count: int,
    *,
    raw_score: float | None = None,
    exon_count: int = 2,
    corroborating_other_source_count: int = 1,
    final_pruned: bool = False,
) -> dict[str, object]:
    return {
        "source_name": source,
        "model_selection_score_scaled": score,
        "model_selection_score": score if raw_score is None else raw_score,
        "sublocus_model_count": sublocus_model_count,
        "sublocus_id": f"sub_{source}_{score}",
        "cds_exon_count": exon_count,
        "corroborating_other_source_count": corroborating_other_source_count,
        "final_pruned": final_pruned,
    }


def source_summary_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["section"] == "source_scaled_score"
        ]


def test_no_singletons_matches_existing_mean() -> None:
    rows = [score_row("A", 0.2, 2), score_row("A", 0.8, 2)]
    stats = summarise_source_scaled_scores(rows)["A"]
    assert stats["count"] == stats["non_singleton_count"] == 2
    assert stats["singleton_count"] == 0
    assert stats["scaled_sum"] / stats["count"] == pytest.approx(0.5)
    assert stats["non_singleton_scaled_sum"] / stats["non_singleton_count"] == pytest.approx(0.5)


def test_all_singletons_has_undefined_structured_mean(tmp_path: Path) -> None:
    rows = [score_row("A", 0.2, 1), score_row("A", 0.8, 1)]
    stats = summarise_source_scaled_scores(rows)["A"]
    assert stats["count"] == stats["singleton_count"] == 2
    assert stats["non_singleton_count"] == 0
    path = tmp_path / "summary.tsv"
    write_model_selection_score_summary_tsv(rows, path)
    row = source_summary_rows(path)[0]
    assert row["mean_scaled_score_all"] == "0.500000"
    assert row["mean_scaled_score_non_singletons"] == "."
    assert row["n_scored_non_singletons"] == "0"
    assert row["n_scored_singletons"] == "2"
    assert row["singleton_fraction"] == "1.000000"


def test_mixed_population_uses_exact_subset() -> None:
    rows = [
        score_row("A", 0.1, 1),
        score_row("A", 0.7, 3),
        score_row("A", 0.9, 3),
    ]
    stats = summarise_source_scaled_scores(rows)["A"]
    assert stats["count"] == 3
    assert stats["non_singleton_count"] == 2
    assert stats["singleton_count"] == 1
    assert stats["non_singleton_scaled_sum"] / stats["non_singleton_count"] == pytest.approx(0.8)


def test_same_source_models_in_one_sublocus_are_non_singletons() -> None:
    rows = [score_row("A", 0.3, 2), score_row("A", 0.9, 2)]
    stats = summarise_source_scaled_scores(rows)["A"]
    assert stats["non_singleton_count"] == 2
    assert stats["singleton_count"] == 0


def test_different_source_models_in_one_sublocus_are_non_singletons() -> None:
    rows = [score_row("A", 0.3, 2), score_row("B", 0.9, 2)]
    stats = summarise_source_scaled_scores(rows)
    assert stats["A"]["non_singleton_count"] == 1
    assert stats["B"]["non_singleton_count"] == 1


def test_exon_count_does_not_classify_singletons() -> None:
    rows = [
        score_row("A", 0.4, 2, exon_count=1),
        score_row("A", 0.6, 1, exon_count=5),
    ]
    stats = summarise_source_scaled_scores(rows)["A"]
    assert stats["non_singleton_scaled_sum"] == pytest.approx(0.4)
    assert stats["singleton_count"] == 1


def test_corroboration_count_does_not_classify_singletons() -> None:
    rows = [
        score_row("A", 0.4, 2, corroborating_other_source_count=0),
        score_row("A", 0.6, 1, corroborating_other_source_count=3),
    ]
    stats = summarise_source_scaled_scores(rows)["A"]
    assert stats["non_singleton_scaled_sum"] == pytest.approx(0.4)
    assert stats["singleton_count"] == 1


def test_pruning_status_does_not_change_population() -> None:
    rows = [
        score_row("A", 0.2, 2, final_pruned=False),
        score_row("A", 0.8, 2, final_pruned=True),
    ]
    stats = summarise_source_scaled_scores(rows)["A"]
    assert stats["non_singleton_count"] == 2
    assert stats["non_singleton_scaled_sum"] == pytest.approx(1.0)


def test_source_order_precision_and_output_are_deterministic(tmp_path: Path) -> None:
    rows = [
        score_row("B", 0.1234567, 2),
        score_row("A", 0.1234567, 1),
        score_row("C", 0.9, 2),
    ]
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    write_model_selection_score_summary_tsv(rows, first)
    write_model_selection_score_summary_tsv(rows, second)
    assert first.read_bytes() == second.read_bytes()
    source_rows = source_summary_rows(first)
    assert [row["group"] for row in source_rows] == ["C", "A", "B"]
    assert source_rows[1]["mean_scaled_score_all"] == "0.123457"


def test_existing_summary_aggregation_is_unchanged() -> None:
    rows = [
        score_row("A", 0.2, 1, raw_score=-0.4),
        score_row("A", 0.8, 2, raw_score=0.6),
    ]
    assert summarise_score_rows_by_source(rows) == {
        "A": {"count": 2.0, "scaled_sum": 1.0, "raw_sum": pytest.approx(0.2)}
    }


def test_invalid_sublocus_model_count_is_rejected() -> None:
    with pytest.raises(ValueError, match="sublocus_model_count"):
        summarise_source_scaled_scores([score_row("A", 0.5, 0)])
