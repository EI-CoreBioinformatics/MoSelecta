from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.models import CodingModel
from moselecta.structural import (
    exact_interval_feature_key,
    make_exact_support_map,
    model_internal_cds_exons,
    model_terminal_cds_exons_with_side,
    terminal_cds_exons_match,
    unsupported_short_cds_exon_penalty_for_model,
)


def _model(name: str, source: str, parts: list[tuple[int, int]], strand: str = "+") -> CodingModel:
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"g_{name}",
        input_transcript_id=name,
        output_gene_id=f"out_g_{name}",
        output_transcript_id=name,
        chrom="chr1",
        strand=strand,
        cds_parts=parts,
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def _penalty(
    query: CodingModel,
    others: list[CodingModel] | None = None,
    *,
    terminal_intrinsic: int = 3,
    internal_intrinsic: int = 6,
    soft_5: int = 15,
    soft_3: int = 15,
    soft_internal: int = 15,
    wobble: int = 3,
) -> dict[str, object]:
    models = [query] + list(others or [])
    internal_support = make_exact_support_map(
        models,
        lambda m: [exact_interval_feature_key(m, s, e) for s, e in model_internal_cds_exons(m)],
    )
    terminal_features = [feature for model in models for feature in model_terminal_cds_exons_with_side(model)]

    def terminal_support_getter(feature):
        model_ids = set()
        sources = set()
        for other in terminal_features:
            if terminal_cds_exons_match(feature, other, wobble):
                model_ids.add(str(other["model_id"]))
                sources.add(str(other["source"]))
        return model_ids, sources

    return unsupported_short_cds_exon_penalty_for_model(
        query,
        internal_exon_support=internal_support,
        terminal_support_getter=terminal_support_getter,
        min_terminal_cds_exon_bp=terminal_intrinsic,
        min_internal_cds_exon_bp=internal_intrinsic,
        soft_5prime_terminal_cds_exon_bp=soft_5,
        soft_3prime_terminal_cds_exon_bp=soft_3,
        soft_internal_cds_exon_bp=soft_internal,
        component_weights={"5prime_terminal": 0.25, "3prime_terminal": 0.25, "internal": 0.50},
    )


def test_unsupported_eight_bp_internal_exon_receives_graded_penalty() -> None:
    row = _penalty(_model("q", "A", [(0, 50), (100, 108), (200, 260)]))
    assert row["unsupported_moderately_short_internal_cds_exon_count"] == 1
    assert row["unsupported_short_internal_cds_exon_severity"] == pytest.approx(7 / 9)
    assert row["unsupported_short_cds_exon_penalty_score"] == pytest.approx(0.5 * 7 / 9)


def test_same_internal_exon_supported_by_another_source_has_no_penalty() -> None:
    query = _model("q", "A", [(0, 50), (100, 108), (200, 260)])
    support = _model("s", "B", [(0, 50), (100, 108), (200, 260)])
    row = _penalty(query, [support])
    assert row["supported_moderately_short_internal_cds_exon_count"] == 1
    assert row["unsupported_moderately_short_internal_cds_exon_count"] == 0
    assert row["unsupported_short_cds_exon_penalty_score"] == 0.0


def test_same_source_model_does_not_count_as_independent_support() -> None:
    query = _model("q", "A", [(0, 50), (100, 108), (200, 260)])
    duplicate = _model("dup", "A", [(0, 50), (100, 108), (200, 260)])
    row = _penalty(query, [duplicate])
    assert row["unsupported_moderately_short_internal_cds_exon_count"] == 1
    assert row["unsupported_short_cds_exon_penalty_score"] > 0


def test_penalty_is_graduated_and_absent_at_or_above_soft_threshold() -> None:
    row_8 = _penalty(_model("q8", "A", [(0, 50), (100, 108), (200, 260)]))
    row_14 = _penalty(_model("q14", "A", [(0, 50), (100, 114), (200, 260)]))
    row_20 = _penalty(_model("q20", "A", [(0, 50), (100, 120), (200, 260)]))
    assert row_8["unsupported_short_cds_exon_penalty_score"] > row_14["unsupported_short_cds_exon_penalty_score"] > 0
    assert row_20["unsupported_short_cds_exon_penalty_score"] == 0.0


def test_intrinsically_short_internal_exon_is_not_double_counted() -> None:
    row = _penalty(_model("q", "A", [(0, 50), (100, 105), (200, 260)]))
    assert row["moderately_short_internal_cds_exon_count"] == 0
    assert row["unsupported_short_cds_exon_penalty_score"] == 0.0


def test_biological_terminal_side_is_strand_aware() -> None:
    plus = _penalty(_model("plus", "A", [(0, 8), (100, 150)], "+"))
    minus = _penalty(_model("minus", "A", [(0, 8), (100, 150)], "-"))
    assert plus["unsupported_short_5prime_terminal_cds_exon_count"] == 1
    assert plus["unsupported_short_3prime_terminal_cds_exon_count"] == 0
    assert minus["unsupported_short_5prime_terminal_cds_exon_count"] == 0
    assert minus["unsupported_short_3prime_terminal_cds_exon_count"] == 1


def test_three_prime_wobble_can_supply_other_source_support() -> None:
    query = _model("q", "A", [(0, 50), (200, 208)], "+")
    support = _model("s", "B", [(0, 50), (200, 210)], "+")
    row = _penalty(query, [support], wobble=3)
    assert row["unsupported_short_3prime_terminal_cds_exon_supported_by_other_source"] == "true"
    assert row["unsupported_short_3prime_terminal_cds_exon_count"] == 0
    assert row["unsupported_short_cds_exon_penalty_score"] == 0.0


def test_soft_interval_disables_cleanly_when_intrinsic_threshold_is_higher() -> None:
    row = _penalty(
        _model("q", "A", [(0, 50), (100, 114), (200, 260)]),
        internal_intrinsic=20,
        soft_internal=15,
    )
    assert row["moderately_short_internal_cds_exon_count"] == 0
    assert row["unsupported_short_cds_exon_penalty_score"] == 0.0


def test_model_only_score_uses_penalty_without_changing_intrinsic_score() -> None:
    from moselecta.scoring import apply_model_selection_scores

    def score_row(name: str, penalty: float) -> dict[str, object]:
        return {
            "sublocus_id": "sub1",
            "source_name": name,
            "output_transcript_id": name,
            "cds_length": 300,
            "sublocus_source_count": 2,
            "sublocus_model_count": 2,
            "pairwise_other_source_count": 0,
            "structural_support_score": 0.0,
            "consensus_similarity_score": 0.0,
            "pairwise_similarity_score": 0.0,
            "model_completeness_score": 0.0,
            "terminal_extension_score": 0.0,
            "exact_structure_support_score": 0.0,
            "intrinsic_penalty_score": 0.0,
            "partial_structure_penalty": 0.0,
            "internal_extra_fraction": 0.0,
            "terminal_extra_fraction": 0.0,
            "unsupported_short_cds_exon_penalty_score": penalty,
        }

    weights = {
        "structural_support": 0.0,
        "consensus_similarity": 0.0,
        "pairwise_similarity": 0.0,
        "completeness": 0.0,
        "source_prior": 0.0,
        "terminal_extension": 0.0,
        "exact_structure_support": 0.0,
        "unsupported_short_cds_exon_penalty": 0.5,
        "intrinsic_penalty": 0.0,
        "partial_structure_penalty": 0.0,
        "internal_extra_penalty": 0.0,
        "terminal_extra_penalty": 0.0,
        "protein_evidence": 0.0,
        "transcript_evidence": 0.0,
    }
    rows = [score_row("clean", 0.0), score_row("short", 0.5)]
    apply_model_selection_scores(rows, {}, weights, 1e-12)
    by_id = {row["output_transcript_id"]: row for row in rows}
    assert by_id["clean"]["model_structural_support_score"] == 0.0
    assert by_id["short"]["model_structural_support_score"] == pytest.approx(-0.25)
    assert by_id["short"]["intrinsic_penalty_score"] == 0.0
    assert by_id["short"]["unsupported_short_cds_exon_penalty_in_selection"] == "true"
    assert by_id["clean"]["model_structural_support_order_rank"] == 1
