from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.models import CodingModel, Sublocus
from moselecta.selection import (
    final_prune_low_score_nonsingletons,
    model_confidence_nonsingleton_threshold,
)


@pytest.mark.parametrize(
    ("sublocus_sources", "run_sources", "expected_threshold", "expected_tier"),
    [
        (2, 3, 0.25, "base"),
        (3, 3, 0.25, "base"),
        (2, 4, 0.25, "base"),
        (2, 5, 0.35, "low_depth_5plus"),
        (3, 5, 0.25, "base"),
        (2, 6, 0.35, "low_depth_5plus"),
        (3, 6, 0.35, "low_depth_5plus"),
        (2, 7, 0.45, "very_low_depth_7plus"),
        (3, 7, 0.35, "low_depth_5plus"),
        (4, 7, 0.25, "base"),
        (24, 83, 0.25, "base"),
        (3, 83, 0.35, "low_depth_5plus"),
        (2, 83, 0.45, "very_low_depth_7plus"),
    ],
)
def test_high_model_confidence_threshold_matrix(
    sublocus_sources: int,
    run_sources: int,
    expected_threshold: float,
    expected_tier: str,
) -> None:
    threshold, tier, fraction = model_confidence_nonsingleton_threshold(
        "high", run_sources, sublocus_sources
    )
    assert threshold == pytest.approx(expected_threshold)
    assert tier == expected_tier
    assert fraction == pytest.approx(sublocus_sources / run_sources)


def test_medium_never_depth_escalates_and_one_source_uses_base() -> None:
    assert model_confidence_nonsingleton_threshold("medium", 7, 2)[:2] == (0.20, "base")
    assert model_confidence_nonsingleton_threshold("high", 7, 1)[:2] == (0.25, "base")


@pytest.mark.parametrize("level", ["medium", "high", "very-high"])
def test_fallback_requires_at_least_three_run_model_sources(level: str) -> None:
    threshold, tier, fraction = model_confidence_nonsingleton_threshold(level, 2, 2)
    assert threshold is None
    assert tier == "insufficient_model_sources"
    assert fraction == pytest.approx(1.0)


def _model(transcript_id: str, source_name: str) -> CodingModel:
    return CodingModel(
        source_name=source_name,
        input_file=Path(f"{source_name}.gff3"),
        input_gene_id=f"gene_{transcript_id}",
        input_transcript_id=transcript_id,
        output_gene_id=f"MoSelecta_gene_{transcript_id}",
        output_transcript_id=transcript_id,
        chrom="chr1",
        strand="+",
        cds_parts=[(99, 150), (249, 400)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def _fallback_run(
    *,
    level: str = "high",
    run_sources: int = 7,
    model_score: float = 0.10,
    protein_supplied: bool = False,
    protein_score: float = 0.0,
    transcript_supplied: bool = False,
    transcript_score: float = 0.0,
    model_selection_score: float = 1.0,
    context_concern: bool = False,
    use_transcript_cdna_context: bool = True,
) -> dict[str, object]:
    candidate = _model("candidate", "A")
    partner = _model("partner", "B")
    sub = Sublocus("sub001", "sl001", "chr1", "+", 99, 400, [candidate, partner])
    row = {
        "output_transcript_id": "candidate",
        "source_name": "A",
        "model_structural_support_score_scaled": model_score,
        "protein_evidence_score_scaled": protein_score,
        "transcript_evidence_score_scaled": transcript_score,
        "model_selection_score_scaled": model_selection_score,
        "final_protein_singleton_rescue_applied": "true",
        "transcript_cdna_context_pruning_concern": "true" if context_concern else "false",
        "transcript_cdna_context_pruning_concern_class": "low_cdna_fraction" if context_concern else ".",
        "transcript_cdna_context_pruning_supporting_flags": "cds_chain_contained_in_longer_cdna" if context_concern else ".",
    }
    result = {
        "initial_retained_models": [candidate],
        "initial_retained_rows": [row],
        "retained_models": [candidate],
        "retained_rows": [row],
        "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }
    return final_prune_low_score_nonsingletons(
        result,
        [sub],
        {},
        enabled=True,
        scope="all",
        protein_evidence_supplied=protein_supplied,
        transcript_evidence_supplied=transcript_supplied,
        protein_adequacy_required=True,
        transcript_adequacy_required=True,
        protein_overlap_model_fraction=0.0,
        protein_supported_model_fraction=0.0,
        transcript_scoring_eligible_model_fraction=0.0,
        transcript_supported_model_fraction=0.0,
        min_protein_overlap_model_fraction=0.75,
        min_protein_supported_model_fraction=0.50,
        min_transcript_scoring_eligible_model_fraction=0.75,
        min_transcript_supported_model_fraction=0.50,
        max_model_structural_support_score_scaled=0.99,
        max_protein_evidence_score_scaled=0.35,
        max_transcript_evidence_score_scaled=0.50,
        max_model_selection_score_scaled=0.99,
        pruning_level_name=level,
        run_model_source_count=run_sources,
        use_transcript_cdna_context=use_transcript_cdna_context,
    )


def test_no_evidence_fallback_uses_model_score_only_and_writes_audit() -> None:
    out = _fallback_run(model_selection_score=1.0)
    assert len(out["removed_rows"]) == 1
    row = out["removed_rows"][0]
    assert row["final_low_score_nonsingleton_policy"] == "model_confidence_fallback"
    assert row["final_low_score_nonsingleton_prune_route"] == "model_confidence_score_thresholds"
    assert row["final_low_score_nonsingleton_run_model_source_count"] == 7
    assert row["final_low_score_nonsingleton_sublocus_source_count"] == 2
    assert row["final_low_score_nonsingleton_source_depth_fraction"] == pytest.approx(2 / 7)
    assert row["final_low_score_nonsingleton_effective_model_confidence_threshold"] == 0.45
    assert row["final_low_score_nonsingleton_model_confidence_tier"] == "very_low_depth_7plus"
    assert row["final_prune_classes"] == "low_score_nonsingleton"


@pytest.mark.parametrize(
    ("protein_supplied", "protein_score", "transcript_supplied", "transcript_score"),
    [
        (True, 0.90, False, 0.0),
        (False, 0.0, True, 0.90),
        (True, 0.90, True, 0.0),
    ],
)
def test_inadequate_supplied_high_local_evidence_vetoes_fallback_score_route(
    protein_supplied: bool,
    protein_score: float,
    transcript_supplied: bool,
    transcript_score: float,
) -> None:
    out = _fallback_run(
        protein_supplied=protein_supplied,
        protein_score=protein_score,
        transcript_supplied=transcript_supplied,
        transcript_score=transcript_score,
    )
    assert len(out["removed_rows"]) == 0
    assert out["retained_rows"][0]["final_low_score_nonsingleton_policy"] == "model_confidence_fallback"


def test_fallback_cdna_context_route_ignores_high_transcript_score() -> None:
    out = _fallback_run(
        transcript_supplied=True,
        transcript_score=0.95,
        context_concern=True,
    )
    row = out["removed_rows"][0]
    assert row["final_low_score_nonsingleton_prune_route"] == "model_confidence_transcript_cdna_context"
    assert "transcript_cdna_context_nonsingleton" in row["final_prune_classes"]
    assert out["low_score_nonsingleton_removed_model_confidence_context_route"] == 1


def test_fallback_cdna_concern_retains_when_context_route_is_disabled() -> None:
    out = _fallback_run(
        transcript_supplied=True,
        transcript_score=0.95,
        context_concern=True,
        use_transcript_cdna_context=False,
    )
    assert len(out["removed_rows"]) == 0
    row = out["retained_rows"][0]
    assert row["final_low_score_nonsingleton_policy"] == "model_confidence_fallback"
    assert row["final_low_score_nonsingleton_use_transcript_cdna_context"] == "false"
    assert row["final_low_score_nonsingleton_prune_reason"] == "score_and_context_thresholds_not_met"


def test_high_local_protein_support_vetoes_fallback_cdna_context_route() -> None:
    out = _fallback_run(
        protein_supplied=True,
        protein_score=0.90,
        transcript_supplied=True,
        transcript_score=0.95,
        context_concern=True,
    )
    assert len(out["removed_rows"]) == 0
    row = out["retained_rows"][0]
    assert row["final_low_score_nonsingleton_policy"] == "model_confidence_fallback"
    assert row["final_low_score_nonsingleton_cdna_context_pruning_concern"] == "true"
    assert row["final_low_score_nonsingleton_prune_reason"] == "score_and_context_thresholds_not_met"


def test_missing_cdna_context_is_neutral_in_fallback() -> None:
    out = _fallback_run(transcript_supplied=True, transcript_score=0.95, context_concern=False)
    assert len(out["removed_rows"]) == 0
    assert out["retained_rows"][0]["final_low_score_nonsingleton_prune_reason"] == "score_and_context_thresholds_not_met"


def test_two_source_guard_blocks_fallback_removal_and_is_audited() -> None:
    out = _fallback_run(run_sources=2)
    assert len(out["removed_rows"]) == 0
    row = out["retained_rows"][0]
    assert row["final_low_score_nonsingleton_model_confidence_tier"] == "insufficient_model_sources"
    assert row["final_low_score_nonsingleton_effective_model_confidence_threshold"] == "."
    assert out["low_score_nonsingleton_skipped_model_source_guard"] == 1


def test_evidence_aware_model_override_does_not_change_canonical_fallback_threshold() -> None:
    out = _fallback_run(model_score=0.50)
    assert len(out["removed_rows"]) == 0
    row = out["retained_rows"][0]
    assert row["final_low_score_nonsingleton_max_model_structural_support_score_scaled"] == 0.99
    assert row["final_low_score_nonsingleton_effective_model_confidence_threshold"] == 0.45


def test_fallback_levels_are_monotonic_for_boundary_scores() -> None:
    scores = (0.19, 0.23, 0.34, 0.44, 0.46)
    removed_by_level = {}
    for level in ("medium", "high", "very-high"):
        removed_by_level[level] = {
            score for score in scores if len(_fallback_run(level=level, model_score=score)["removed_rows"]) == 1
        }
    assert removed_by_level["medium"] <= removed_by_level["high"]
    assert removed_by_level["high"] == removed_by_level["very-high"]
