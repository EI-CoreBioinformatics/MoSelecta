"""Regression tests for selector-consistent model-selection order ranks."""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import pytest

from moselecta.scoring import (
    _assign_model_selection_score_ranks,
    _assign_score_ranks,
    primary_model_selection_order_key,
)
from moselecta.selection import model_selection_row_sort_key, select_sublocus_winners


@dataclass(frozen=True)
class SyntheticModel:
    output_transcript_id: str


@dataclass(frozen=True)
class SyntheticSublocus:
    sublocus_id: str
    chrom: str
    start0: int
    end0: int
    models: tuple[SyntheticModel, ...]


def row(model_id: str, score: float = 1.0, **updates: object) -> dict[str, object]:
    result: dict[str, object] = {
        "sublocus_id": "synthetic.sublocus001",
        "output_transcript_id": model_id,
        "model_selection_score": score,
        "model_completeness_score": 0.8,
        "intrinsic_penalty_score": 0.1,
        "unsupported_short_cds_exon_penalty_score": 0.0,
        "sequence_integrity_penalty_contribution": 0.0,
        "consensus_similarity_score": 0.6,
        "structural_support_score": 0.6,
        "exact_structure_support_score": 0.0,
        "source_prior_score": 0.5,
        "cds_length": 900,
        "model_selection_score_rank": ".",
        "model_selection_order_rank": ".",
        "model_selection_tie_group": ".",
        "model_selection_tie_break_reason": ".",
    }
    result.update(updates)
    return result


def assign(rows: list[dict[str, object]]) -> None:
    _assign_model_selection_score_ranks(
        {"synthetic.sublocus001": rows},
        tie_epsilon=1e-6,
    )


def ids_by_rank(rows: list[dict[str, object]]) -> list[str]:
    return [
        str(item["output_transcript_id"])
        for item in sorted(rows, key=lambda item: int(item["model_selection_order_rank"]))
    ]


def selected_winner(rows: list[dict[str, object]]) -> str:
    by_id = {str(item["output_transcript_id"]): item for item in rows}
    models = tuple(SyntheticModel(model_id) for model_id in by_id)
    sublocus = SyntheticSublocus("synthetic.sublocus001", "chrSynthetic", 0, 1000, models)
    return select_sublocus_winners([sublocus], by_id)[0].output_transcript_id


@pytest.mark.parametrize(
    ("report_only_updates", "selector_updates", "expected_reason"),
    [
        (
            {"unsupported_short_cds_exon_penalty_score": 0.0, "consensus_similarity_score": 0.2},
            {"unsupported_short_cds_exon_penalty_score": 0.8, "consensus_similarity_score": 0.9},
            "consensus_similarity_score",
        ),
        (
            {"sequence_integrity_penalty_contribution": -0.8, "consensus_similarity_score": 0.2},
            {"sequence_integrity_penalty_contribution": 0.0, "consensus_similarity_score": 0.9},
            "consensus_similarity_score",
        ),
        (
            {"exact_structure_support_score": 0.9, "source_prior_score": 0.1},
            {"exact_structure_support_score": 0.0, "source_prior_score": 0.9},
            "source_prior_score",
        ),
    ],
)
@pytest.mark.parametrize("score", [1.0, 0.0, -1.0])
def test_confirmed_report_only_reversals_now_follow_selector(
    report_only_updates: dict[str, object],
    selector_updates: dict[str, object],
    expected_reason: str,
    score: float,
) -> None:
    rows = [
        row("old_report_first", score, **report_only_updates),
        row("selector_first", score, **selector_updates),
    ]
    assign(rows)

    assert ids_by_rank(rows) == ["selector_first", "old_report_first"]
    assert selected_winner(rows) == "selector_first"
    assert {item["model_selection_tie_break_reason"] for item in rows} == {expected_reason}
    assert {item["model_selection_score_rank"] for item in rows} == {1}


@pytest.mark.parametrize(
    ("field", "first_value", "second_value", "expected_first"),
    [
        ("model_completeness_score", 0.9, 0.7, "a"),
        ("intrinsic_penalty_score", 0.0, 0.2, "a"),
        ("consensus_similarity_score", 0.9, 0.7, "a"),
        ("structural_support_score", 0.9, 0.7, "a"),
        ("source_prior_score", 0.9, 0.7, "a"),
        ("cds_length", 1200, 900, "a"),
    ],
)
def test_each_numeric_secondary_key_drives_rank_and_reason(
    field: str,
    first_value: object,
    second_value: object,
    expected_first: str,
) -> None:
    rows = [row("a", **{field: first_value}), row("b", **{field: second_value})]
    assign(rows)

    assert ids_by_rank(rows)[0] == expected_first
    assert {item["model_selection_tie_break_reason"] for item in rows} == {field}


def test_lexical_output_id_is_final_fallback() -> None:
    rows = [row("id_b"), row("id_a")]
    assign(rows)

    assert ids_by_rank(rows) == ["id_a", "id_b"]
    assert {item["model_selection_tie_break_reason"] for item in rows} == {"output_transcript_id"}


@pytest.mark.parametrize(
    ("higher_score", "expected_shared_score_rank", "expected_reason"),
    [
        (0.0, True, "consensus_similarity_score"),
        (0.5e-6, True, "model_selection_score"),
        (1.0e-6, True, "model_selection_score"),
        (1.01e-6, False, "."),
        (1.0, False, "."),
    ],
)
def test_exact_epsilon_and_unequal_score_semantics(
    higher_score: float,
    expected_shared_score_rank: bool,
    expected_reason: str,
) -> None:
    rows = [
        row("higher", higher_score, consensus_similarity_score=0.9),
        row("lower", 0.0, consensus_similarity_score=0.1),
    ]
    assign(rows)

    assert ids_by_rank(rows)[0] == "higher"
    ranks = [int(item["model_selection_score_rank"]) for item in rows]
    assert (len(set(ranks)) == 1) is expected_shared_score_rank
    assert {item["model_selection_tie_break_reason"] for item in rows} == {expected_reason}


def test_score_rank_remains_competition_rank() -> None:
    rows = [row("a", 1.0), row("b", 1.0), row("c", 0.5)]
    assign(rows)
    ranks = {str(item["output_transcript_id"]): item["model_selection_score_rank"] for item in rows}

    assert ranks == {"a": 1, "b": 1, "c": 3}


def test_three_way_complete_order_is_permutation_invariant() -> None:
    original = [
        row("low", 1.0, consensus_similarity_score=0.1),
        row("high", 1.0, consensus_similarity_score=0.9),
        row("middle", 1.0, consensus_similarity_score=0.5),
    ]
    observed: set[tuple[str, ...]] = set()
    for permutation in itertools.permutations(original):
        rows = [dict(item) for item in permutation]
        assign(rows)
        observed.add(tuple(ids_by_rank(rows)))
        assert ids_by_rank(rows) == [
            str(item["output_transcript_id"])
            for item in sorted(rows, key=primary_model_selection_order_key)
        ]
        assert selected_winner(rows) == ids_by_rank(rows)[0]

    assert observed == {("high", "middle", "low")}


def test_selection_wrapper_delegates_to_authoritative_helper() -> None:
    rows = [
        row("a", 1.0, consensus_similarity_score=0.4),
        row("b", 1.0, consensus_similarity_score=0.8),
        row("c", 2.0),
    ]
    assert [model_selection_row_sort_key(item) for item in rows] == [
        primary_model_selection_order_key(item) for item in rows
    ]


def test_missing_and_nan_values_keep_existing_numeric_fallback() -> None:
    rows = [
        row("missing", model_completeness_score=None, consensus_similarity_score=0.2),
        row("nan", model_completeness_score=math.nan, consensus_similarity_score=0.9),
    ]
    assign(rows)

    assert ids_by_rank(rows) == ["nan", "missing"]
    assert {item["model_selection_tie_break_reason"] for item in rows} == {
        "consensus_similarity_score"
    }


def test_structural_support_order_retains_legacy_report_only_tie_breakers() -> None:
    rows = [
        row(
            "clean",
            model_structural_support_score=1.0,
            unsupported_short_cds_exon_penalty_score=0.0,
            consensus_similarity_score=0.2,
        ),
        row(
            "higher_consensus",
            model_structural_support_score=1.0,
            unsupported_short_cds_exon_penalty_score=0.8,
            consensus_similarity_score=0.9,
        ),
    ]
    _assign_score_ranks(
        {"synthetic.sublocus001": rows},
        score_field="model_structural_support_score",
        score_rank_field="model_structural_support_score_rank",
        order_rank_field="model_structural_support_order_rank",
        tie_group_field="model_structural_support_tie_group",
        tie_break_reason_field="model_structural_support_tie_break_reason",
        tie_epsilon=1e-6,
    )
    by_id = {str(item["output_transcript_id"]): item for item in rows}

    assert by_id["clean"]["model_structural_support_order_rank"] == 1
    assert by_id["higher_consensus"]["model_structural_support_order_rank"] == 2
    assert {item["model_structural_support_tie_break_reason"] for item in rows} == {
        "unsupported_short_cds_exon_penalty_score"
    }


def test_scaled_score_clipping_does_not_replace_raw_total_order() -> None:
    rows = [
        row("less_negative", -0.1, model_selection_score_scaled=0.0),
        row("more_negative", -1.0, model_selection_score_scaled=0.0),
    ]
    assign(rows)

    assert ids_by_rank(rows) == ["less_negative", "more_negative"]
