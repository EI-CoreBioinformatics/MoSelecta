from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from moselecta.models import CodingModel, Sublocus


def _selection_weights(transcript_weight: float) -> dict[str, float]:
    return {
        "structural_support": 0.0,
        "consensus_similarity": 0.0,
        "pairwise_similarity": 0.0,
        "completeness": 0.0,
        "source_prior": 0.0,
        "terminal_extension": 0.0,
        "exact_structure_support": 0.0,
        "intrinsic_penalty": 0.0,
        "partial_structure_penalty": 0.0,
        "internal_extra_penalty": 0.0,
        "terminal_extra_penalty": 0.0,
        "protein_evidence": 0.0,
        "transcript_evidence": transcript_weight,
    }


def _score_row(transcript_id: str) -> dict[str, object]:
    return {
        "sublocus_id": "sub001",
        "source_name": transcript_id,
        "output_transcript_id": transcript_id,
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
    }


def _model(transcript_id: str, start: int) -> CodingModel:
    return CodingModel(
        source_name=transcript_id,
        input_file=Path(f"{transcript_id}.gff3"),
        input_gene_id=f"g_{transcript_id}",
        input_transcript_id=transcript_id,
        output_gene_id=f"out_g_{transcript_id}",
        output_transcript_id=transcript_id,
        chrom="chr1",
        strand="+",
        cds_parts=[(start, start + 50), (start + 100, start + 200)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def test_positive_transcript_selection_weight_changes_model_selection_ranking() -> None:
    from moselecta.scoring import apply_model_selection_scores

    base_rows = [_score_row("txA"), _score_row("txB")]
    evidence = {
        "txA": {"transcript_evidence_score": 0.0, "transcript_evidence_score_scaled": 0.0},
        "txB": {"transcript_evidence_score": 1.0, "transcript_evidence_score_scaled": 1.0},
    }

    inactive_rows = deepcopy(base_rows)
    apply_model_selection_scores(
        inactive_rows,
        source_prior_scores={},
        selection_weights=_selection_weights(1.0),
        tie_epsilon=1e-12,
        transcript_evidence_rows_by_transcript=evidence,
        transcript_evidence_active=False,
    )
    inactive = {row["output_transcript_id"]: row for row in inactive_rows}
    assert inactive["txA"]["model_selection_order_rank"] == 1
    assert inactive["txB"]["model_selection_order_rank"] == 2
    assert inactive["txB"]["transcript_evidence_score_in_selection"] == "false"

    active_rows = deepcopy(base_rows)
    apply_model_selection_scores(
        active_rows,
        source_prior_scores={},
        selection_weights=_selection_weights(1.0),
        tie_epsilon=1e-12,
        transcript_evidence_rows_by_transcript=evidence,
        transcript_evidence_active=True,
    )
    active = {row["output_transcript_id"]: row for row in active_rows}
    assert active["txB"]["model_selection_order_rank"] == 1
    assert active["txA"]["model_selection_order_rank"] == 2
    assert active["txB"]["transcript_evidence_score_in_selection"] == "true"
    assert active["txB"]["score_delta_from_evidence"] == 1.0


def test_transcript_evidence_rescues_an_otherwise_pruned_singleton() -> None:
    from moselecta.selection import final_prune_singleton_selected_models

    model = _model("tx1", 100)
    sub = Sublocus(
        sublocus_id="sub001",
        superlocus_id="sl001",
        chrom="chr1",
        strand="+",
        start0=model.cds_start0,
        end0=model.cds_end0,
        models=[model],
    )
    base_row = {
        "output_transcript_id": "tx1",
        "source_name": "tx1",
        "final_prune_classes": ".",
    }
    prune_result = {
        "initial_retained_models": [model],
        "initial_retained_rows": [base_row],
        "retained_models": [model],
        "retained_rows": [base_row],
        "removed_rows": [],
        "retained_selection_info_by_transcript": {"tx1": base_row},
    }
    transcript_rows = {
        "tx1": {
            "transcript_evidence_score_scaled": 0.90,
            "transcript_supporting_source_count": 2,
            "transcript_scoring_eligible_count": 2,
            "transcript_stranded_scoring_eligible_count": 2,
            "transcript_unstranded_scoring_eligible_count": 0,
            "transcript_evidence_status": "strong_transcript_support",
            "transcript_evidence_flags": "strong_transcript_support",
            "transcript_cdna_context_evaluable": "true",
            "transcript_cdna_context_cds_to_min_cdna_fraction": 0.50,
        }
    }

    out = final_prune_singleton_selected_models(
        prune_result,
        [sub],
        {"tx1": {"intrinsic_penalty_components": "short_cds", "intrinsic_penalty_score": 1.0}},
        prune_level="all",
        prune_sources=set(),
        pruning_level_name="high",
        rescue_transcript_supported_singletons=True,
        transcript_evidence_rows_by_transcript=transcript_rows,
        rescue_min_transcript_evidence_score_scaled=0.75,
        rescue_min_transcript_supporting_source_count=1,
        rescue_min_transcript_scoring_eligible_count=1,
    )

    assert out["retained_models"] == [model]
    assert out["removed_rows"] == []
    row = out["retained_rows"][0]
    assert row["final_transcript_singleton_rescue_criteria_met"] == "true"
    assert row["final_transcript_singleton_rescue_applied"] == "true"
    assert row["final_transcript_singleton_rescued"] == "true"
    assert row["final_transcript_singleton_rescue_cdna_context_fraction"] == 0.50


def test_transcript_singleton_rescue_requires_minimum_cdna_context_fraction() -> None:
    from moselecta.selection import final_prune_singleton_selected_models

    model = _model("tx_low_context", 100)
    sub = Sublocus("sub001", "sl001", "chr1", "+", model.cds_start0, model.cds_end0, [model])
    base_row = {
        "output_transcript_id": model.output_transcript_id,
        "source_name": model.source_name,
        "final_prune_classes": ".",
    }
    prune_result = {
        "initial_retained_models": [model],
        "initial_retained_rows": [base_row],
        "retained_models": [model],
        "retained_rows": [base_row],
        "removed_rows": [],
        "retained_selection_info_by_transcript": {model.output_transcript_id: base_row},
    }
    transcript_rows = {
        model.output_transcript_id: {
            "transcript_evidence_score_scaled": 0.90,
            "transcript_supporting_source_count": 2,
            "transcript_scoring_eligible_count": 2,
            "transcript_stranded_scoring_eligible_count": 2,
            "transcript_unstranded_scoring_eligible_count": 0,
            "transcript_cdna_context_evaluable": "true",
            "transcript_cdna_context_cds_to_min_cdna_fraction": 0.299,
            "transcript_cdna_context_pruning_concern": "false",
        }
    }
    out = final_prune_singleton_selected_models(
        prune_result,
        [sub],
        {model.output_transcript_id: {"intrinsic_penalty_components": "short_cds", "intrinsic_penalty_score": 1.0}},
        prune_level="all",
        prune_sources=set(),
        pruning_level_name="high",
        rescue_transcript_supported_singletons=True,
        transcript_evidence_rows_by_transcript=transcript_rows,
    )
    assert out["retained_models"] == []
    row = out["removed_rows"][0]
    assert row["final_transcript_singleton_rescue_criteria_met"] == "false"
    assert row["final_transcript_singleton_rescue_cdna_context_reason"] == "context_fraction_failed_below_minimum"
    assert row["final_transcript_singleton_rescue_cdna_context_fraction_state"] == "failed_below_minimum"


def test_transcript_singleton_rescue_accepts_boundary_and_neutral_not_evaluable_context() -> None:
    from moselecta.selection import singleton_transcript_rescue_fields

    common = {
        "transcript_evidence_score_scaled": 0.90,
        "transcript_supporting_source_count": 2,
        "transcript_scoring_eligible_count": 2,
        "transcript_stranded_scoring_eligible_count": 2,
        "transcript_unstranded_scoring_eligible_count": 0,
        "transcript_cdna_context_pruning_concern": "false",
    }
    rows = {
        "boundary": {
            **common,
            "transcript_cdna_context_evaluable": "true",
            "transcript_cdna_context_cds_to_min_cdna_fraction": 0.30,
        },
        "not_evaluable": {
            **common,
            "transcript_cdna_context_evaluable": "false",
            "transcript_cdna_context_cds_to_min_cdna_fraction": ".",
        },
    }

    boundary = singleton_transcript_rescue_fields(
        "boundary",
        rows,
        rescue_transcript_supported_singletons=True,
        min_score_scaled=0.75,
        min_supporting_source_count=1,
        min_scoring_eligible_count=1,
        min_cdna_context_fraction=0.30,
        require_stranded_support=False,
    )
    neutral = singleton_transcript_rescue_fields(
        "not_evaluable",
        rows,
        rescue_transcript_supported_singletons=True,
        min_score_scaled=0.75,
        min_supporting_source_count=1,
        min_scoring_eligible_count=1,
        min_cdna_context_fraction=0.30,
        require_stranded_support=False,
    )

    assert boundary["final_transcript_singleton_rescue_criteria_met"] == "true"
    assert boundary["final_transcript_singleton_rescue_cdna_context_reason"] == "context_fraction_passed"
    assert boundary["final_transcript_singleton_rescue_cdna_context_fraction_state"] == "passed"
    assert neutral["final_transcript_singleton_rescue_criteria_met"] == "true"
    assert neutral["final_transcript_singleton_rescue_cdna_context_reason"] == "context_fraction_not_evaluable_neutral"
    assert neutral["final_transcript_singleton_rescue_cdna_context_fraction_state"] == "not_evaluable_neutral"


def test_high_transcript_score_protects_singleton_from_low_score_pruning() -> None:
    from moselecta.selection import final_prune_low_score_singletons

    low = _model("low", 100)
    high = _model("high", 1000)
    sub_low = Sublocus("sub_low", "sl_low", "chr1", "+", low.cds_start0, low.cds_end0, [low])
    sub_high = Sublocus("sub_high", "sl_high", "chr1", "+", high.cds_start0, high.cds_end0, [high])
    common = {
        "model_structural_support_score_scaled": 0.0,
        "protein_evidence_score_scaled": 0.0,
        "model_selection_score_scaled": 0.0,
        "final_protein_singleton_rescue_applied": "false",
        "final_transcript_singleton_rescue_applied": "false",
    }
    rows = [
        {"output_transcript_id": "low", "source_name": "low", "transcript_evidence_score_scaled": 0.10, **common},
        {"output_transcript_id": "high", "source_name": "high", "transcript_evidence_score_scaled": 0.90, **common},
    ]
    prune_result = {
        "initial_retained_models": [low, high],
        "initial_retained_rows": rows,
        "retained_models": [low, high],
        "retained_rows": rows,
        "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }

    out = final_prune_low_score_singletons(
        prune_result,
        [sub_low, sub_high],
        {},
        enabled=True,
        protein_evidence_supplied=False,
        transcript_evidence_supplied=True,
        protein_adequacy_required=True,
        transcript_adequacy_required=True,
        protein_overlap_model_fraction=0.0,
        protein_supported_model_fraction=0.0,
        transcript_scoring_eligible_model_fraction=1.0,
        transcript_supported_model_fraction=1.0,
        min_protein_overlap_model_fraction=0.75,
        min_protein_supported_model_fraction=0.50,
        min_transcript_scoring_eligible_model_fraction=0.75,
        min_transcript_supported_model_fraction=0.50,
        max_model_structural_support_score_scaled=0.25,
        max_protein_evidence_score_scaled=0.25,
        max_transcript_evidence_score_scaled=0.25,
        max_model_selection_score_scaled=0.25,
        pruning_level_name="mid",
    )

    assert [m.output_transcript_id for m in out["retained_models"]] == ["high"]
    assert [row["output_transcript_id"] for row in out["removed_rows"]] == ["low"]
    assert out["retained_rows"][0]["final_low_score_singleton_prune_reason"] == "score_and_context_thresholds_not_met"
    assert out["removed_rows"][0]["final_low_score_singleton_pruned"] == "true"


def test_high_transcript_score_protects_nonsingleton_from_low_score_pruning() -> None:
    from moselecta.selection import final_prune_low_score_nonsingletons

    low = _model("low", 100)
    high = _model("high", 100)
    sub = Sublocus("sub001", "sl001", "chr1", "+", 100, 300, [low, high])
    common = {
        "model_structural_support_score_scaled": 0.0,
        "protein_evidence_score_scaled": 0.0,
        "model_selection_score_scaled": 0.0,
        "final_protein_singleton_rescue_applied": "false",
        "final_transcript_singleton_rescue_applied": "false",
        "final_selection_stage": "backfill",
        "backfill_selected": "true",
    }
    rows = [
        {"output_transcript_id": "low", "source_name": "low", "transcript_evidence_score_scaled": 0.10, **common},
        {"output_transcript_id": "high", "source_name": "high", "transcript_evidence_score_scaled": 0.90, **common},
    ]
    prune_result = {
        "initial_retained_models": [low, high],
        "initial_retained_rows": rows,
        "retained_models": [low, high],
        "retained_rows": rows,
        "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }

    out = final_prune_low_score_nonsingletons(
        prune_result,
        [sub],
        {},
        enabled=True,
        scope="all",
        protein_evidence_supplied=False,
        transcript_evidence_supplied=True,
        protein_adequacy_required=True,
        transcript_adequacy_required=True,
        protein_overlap_model_fraction=0.0,
        protein_supported_model_fraction=0.0,
        transcript_scoring_eligible_model_fraction=1.0,
        transcript_supported_model_fraction=1.0,
        min_protein_overlap_model_fraction=0.75,
        min_protein_supported_model_fraction=0.50,
        min_transcript_scoring_eligible_model_fraction=0.75,
        min_transcript_supported_model_fraction=0.50,
        max_model_structural_support_score_scaled=0.20,
        max_protein_evidence_score_scaled=0.20,
        max_transcript_evidence_score_scaled=0.20,
        max_model_selection_score_scaled=0.20,
        pruning_level_name="mid",
    )

    assert [m.output_transcript_id for m in out["retained_models"]] == ["high"]
    assert [row["output_transcript_id"] for row in out["removed_rows"]] == ["low"]
    assert out["retained_rows"][0]["final_low_score_nonsingleton_prune_reason"] == "score_and_context_thresholds_not_met"
    assert out["removed_rows"][0]["final_low_score_nonsingleton_pruned"] == "true"


def test_uncalculated_cdna_context_is_neutral_in_general_scoring() -> None:
    from moselecta.scoring import apply_model_selection_scores

    rows = [_score_row("tx1")]
    evidence = {
        "tx1": {
            "transcript_evidence_score": 0.85,
            "transcript_evidence_score_scaled": 0.85,
            "transcript_cdna_context_status": "not_calculated_no_compatible_transcript",
            "transcript_cdna_context_cds_to_min_cdna_fraction": ".",
            "transcript_cdna_context_flags": ".",
        }
    }
    apply_model_selection_scores(
        rows,
        source_prior_scores={},
        selection_weights=_selection_weights(1.0),
        tie_epsilon=1e-12,
        transcript_evidence_rows_by_transcript=evidence,
        transcript_evidence_active=True,
        transcript_cdna_context_penalty_weight=0.15,
    )
    row = rows[0]
    assert row["transcript_cdna_context_evaluable"] == "false"
    assert row["transcript_cdna_context_concern_score"] == 0.0
    assert row["transcript_cdna_context_penalty"] == 0.0
    assert row["transcript_evidence_score_context_adjusted"] == 0.85
    assert row["model_selection_score"] == 0.85


def test_low_cdna_fraction_modestly_adjusts_general_transcript_score() -> None:
    from moselecta.scoring import apply_model_selection_scores

    rows = [_score_row("tx1")]
    evidence = {
        "tx1": {
            "transcript_evidence_score": 1.0,
            "transcript_evidence_score_scaled": 1.0,
            "transcript_cdna_context_status": "calculated",
            "transcript_cdna_context_cds_to_min_cdna_fraction": 0.05,
            "transcript_cdna_context_flags": ".",
        }
    }
    apply_model_selection_scores(
        rows,
        source_prior_scores={},
        selection_weights=_selection_weights(1.0),
        tie_epsilon=1e-12,
        transcript_evidence_rows_by_transcript=evidence,
        transcript_evidence_active=True,
        transcript_cdna_context_neutral_fraction=0.30,
        transcript_cdna_context_severe_fraction=0.15,
        transcript_cdna_context_penalty_weight=0.15,
    )
    row = rows[0]
    assert abs(float(row["transcript_cdna_context_concern_score"]) - (5.0 / 6.0)) < 1e-12
    assert abs(float(row["transcript_cdna_context_penalty"]) - 0.125) < 1e-12
    assert abs(float(row["transcript_evidence_score_context_adjusted"]) - 0.875) < 1e-12
    assert abs(float(row["model_selection_score"]) - 0.875) < 1e-12
    assert row["transcript_cdna_context_pruning_concern"] == "true"
    assert row["transcript_cdna_context_pruning_concern_class"] == "severe_low_fraction"


def test_moderate_cdna_fraction_needs_structural_flag_for_pruning_concern() -> None:
    from moselecta.scoring import transcript_cdna_context_concern_fields

    no_flag = transcript_cdna_context_concern_fields({
        "transcript_cdna_context_status": "calculated",
        "transcript_cdna_context_cds_to_min_cdna_fraction": 0.20,
        "transcript_cdna_context_flags": ".",
    })
    with_flag = transcript_cdna_context_concern_fields({
        "transcript_cdna_context_status": "calculated",
        "transcript_cdna_context_cds_to_min_cdna_fraction": 0.20,
        "transcript_cdna_context_flags": "cds_chain_contained_in_longer_cdna",
    })
    assert no_flag["transcript_cdna_context_concern_score"] > 0
    assert no_flag["transcript_cdna_context_pruning_concern"] == "false"
    assert no_flag["transcript_cdna_context_pruning_concern_class"] == "low_fraction_without_structural_context"
    assert with_flag["transcript_cdna_context_pruning_concern"] == "true"
    assert with_flag["transcript_cdna_context_pruning_concern_class"] == "moderate_low_fraction_with_structural_context"


def test_severe_cdna_context_blocks_transcript_only_singleton_rescue() -> None:
    from moselecta.selection import final_prune_singleton_selected_models

    model = _model("tx1", 100)
    sub = Sublocus("sub001", "sl001", "chr1", "+", model.cds_start0, model.cds_end0, [model])
    base_row = {"output_transcript_id": "tx1", "source_name": "tx1", "final_prune_classes": "."}
    prune_result = {
        "initial_retained_models": [model],
        "initial_retained_rows": [base_row],
        "retained_models": [model],
        "retained_rows": [base_row],
        "removed_rows": [],
        "retained_selection_info_by_transcript": {"tx1": base_row},
    }
    transcript_rows = {
        "tx1": {
            "transcript_evidence_score_scaled": 1.0,
            "transcript_supporting_source_count": 2,
            "transcript_scoring_eligible_count": 2,
            "transcript_stranded_scoring_eligible_count": 2,
            "transcript_unstranded_scoring_eligible_count": 0,
            "transcript_evidence_status": "strong_transcript_support",
            "transcript_evidence_flags": "strong_transcript_support",
            "transcript_cdna_context_evaluable": "true",
            "transcript_cdna_context_cds_to_min_cdna_fraction": 0.05,
            "transcript_cdna_context_pruning_concern": "true",
            "transcript_cdna_context_pruning_concern_class": "severe_low_fraction",
        }
    }
    out = final_prune_singleton_selected_models(
        prune_result,
        [sub],
        {"tx1": {"intrinsic_penalty_components": "short_cds", "intrinsic_penalty_score": 1.0}},
        prune_level="all",
        prune_sources=set(),
        pruning_level_name="high",
        rescue_transcript_supported_singletons=True,
        transcript_evidence_rows_by_transcript=transcript_rows,
        transcript_rescue_block_on_cdna_context_concern=True,
    )
    assert out["retained_models"] == []
    assert len(out["removed_rows"]) == 1
    row = out["removed_rows"][0]
    assert row["final_transcript_singleton_rescue_criteria_met"] == "false"
    assert row["final_transcript_singleton_rescue_cdna_context_eligible"] == "false"
    assert row["final_transcript_singleton_rescue_cdna_context_reason"] == "blocked_by_severe_low_fraction"


def test_high_transcript_score_can_be_pruned_by_singleton_cdna_context_route() -> None:
    from moselecta.selection import final_prune_low_score_singletons

    model = _model("tx1", 100)
    sub = Sublocus("sub001", "sl001", "chr1", "+", model.cds_start0, model.cds_end0, [model])
    row = {
        "output_transcript_id": "tx1",
        "source_name": "tx1",
        "model_structural_support_score_scaled": 0.10,
        "protein_evidence_score_scaled": 0.0,
        "transcript_evidence_score_scaled": 1.0,
        "model_selection_score_scaled": 1.0,
        "final_protein_singleton_rescue_applied": "false",
        "final_transcript_singleton_rescue_applied": "false",
        "transcript_cdna_context_pruning_concern": "true",
        "transcript_cdna_context_pruning_concern_class": "severe_low_fraction",
        "transcript_cdna_context_pruning_supporting_flags": ".",
    }
    prune_result = {
        "initial_retained_models": [model], "initial_retained_rows": [row],
        "retained_models": [model], "retained_rows": [row], "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }
    kwargs = dict(
        enabled=True,
        protein_evidence_supplied=False,
        transcript_evidence_supplied=True,
        protein_adequacy_required=True,
        transcript_adequacy_required=True,
        protein_overlap_model_fraction=0.0,
        protein_supported_model_fraction=0.0,
        transcript_scoring_eligible_model_fraction=1.0,
        transcript_supported_model_fraction=1.0,
        min_protein_overlap_model_fraction=0.75,
        min_protein_supported_model_fraction=0.50,
        min_transcript_scoring_eligible_model_fraction=0.75,
        min_transcript_supported_model_fraction=0.50,
        max_model_structural_support_score_scaled=0.25,
        max_protein_evidence_score_scaled=0.25,
        max_transcript_evidence_score_scaled=0.40,
        max_model_selection_score_scaled=0.25,
        pruning_level_name="mid",
    )
    out = final_prune_low_score_singletons(prune_result, [sub], {}, **kwargs)
    assert out["retained_models"] == []
    removed = out["removed_rows"][0]
    assert removed["final_low_score_singleton_prune_route"] == "transcript_cdna_context"
    assert "transcript_cdna_context_singleton" in removed["final_prune_classes"]
    assert out["low_score_singleton_removed_context_route"] == 1

    disabled = final_prune_low_score_singletons(
        prune_result, [sub], {}, use_transcript_cdna_context=False, **kwargs
    )
    assert disabled["retained_models"] == [model]
    assert disabled["retained_rows"][0]["final_low_score_singleton_prune_reason"] == "score_and_context_thresholds_not_met"


def test_cdna_context_route_requires_transcript_adequacy() -> None:
    from moselecta.selection import final_prune_low_score_singletons

    model = _model("tx1", 100)
    sub = Sublocus("sub001", "sl001", "chr1", "+", model.cds_start0, model.cds_end0, [model])
    row = {
        "output_transcript_id": "tx1", "source_name": "tx1",
        "model_structural_support_score_scaled": 0.10,
        "protein_evidence_score_scaled": 0.0,
        "transcript_evidence_score_scaled": 1.0,
        "model_selection_score_scaled": 1.0,
        "final_protein_singleton_rescue_applied": "false",
        "final_transcript_singleton_rescue_applied": "false",
        "transcript_cdna_context_pruning_concern": "true",
        "transcript_cdna_context_pruning_concern_class": "severe_low_fraction",
    }
    prune_result = {
        "initial_retained_models": [model], "initial_retained_rows": [row],
        "retained_models": [model], "retained_rows": [row], "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }
    out = final_prune_low_score_singletons(
        prune_result, [sub], {}, enabled=True,
        protein_evidence_supplied=False, transcript_evidence_supplied=True,
        protein_adequacy_required=True, transcript_adequacy_required=True,
        protein_overlap_model_fraction=0.0, protein_supported_model_fraction=0.0,
        transcript_scoring_eligible_model_fraction=0.20, transcript_supported_model_fraction=0.20,
        min_protein_overlap_model_fraction=0.75, min_protein_supported_model_fraction=0.50,
        min_transcript_scoring_eligible_model_fraction=0.75, min_transcript_supported_model_fraction=0.50,
        max_model_structural_support_score_scaled=0.25, max_protein_evidence_score_scaled=0.25,
        max_transcript_evidence_score_scaled=0.40, max_model_selection_score_scaled=0.25,
        pruning_level_name="mid",
    )
    assert out["retained_models"] == [model]
    assert out["retained_rows"][0]["final_low_score_singleton_prune_reason"] == "external_evidence_adequacy_failed"


def test_high_transcript_score_backfill_nonsingleton_can_be_pruned_by_cdna_context_route() -> None:
    from moselecta.selection import final_prune_low_score_nonsingletons

    model = _model("tx1", 100)
    partner = _model("tx2", 100)
    sub = Sublocus("sub001", "sl001", "chr1", "+", 100, 300, [model, partner])
    row = {
        "output_transcript_id": "tx1", "source_name": "tx1",
        "model_structural_support_score_scaled": 0.10,
        "protein_evidence_score_scaled": 0.0,
        "transcript_evidence_score_scaled": 1.0,
        "model_selection_score_scaled": 1.0,
        "final_protein_singleton_rescue_applied": "false",
        "final_transcript_singleton_rescue_applied": "false",
        "final_selection_stage": "backfill", "backfill_selected": "true",
        "transcript_cdna_context_pruning_concern": "true",
        "transcript_cdna_context_pruning_concern_class": "moderate_low_fraction_with_structural_context",
        "transcript_cdna_context_pruning_supporting_flags": "cds_chain_contained_in_longer_cdna",
    }
    prune_result = {
        "initial_retained_models": [model], "initial_retained_rows": [row],
        "retained_models": [model], "retained_rows": [row], "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }
    out = final_prune_low_score_nonsingletons(
        prune_result, [sub], {}, enabled=True, scope="backfill",
        protein_evidence_supplied=False, transcript_evidence_supplied=True,
        protein_adequacy_required=True, transcript_adequacy_required=True,
        protein_overlap_model_fraction=0.0, protein_supported_model_fraction=0.0,
        transcript_scoring_eligible_model_fraction=1.0, transcript_supported_model_fraction=1.0,
        min_protein_overlap_model_fraction=0.75, min_protein_supported_model_fraction=0.50,
        min_transcript_scoring_eligible_model_fraction=0.75, min_transcript_supported_model_fraction=0.50,
        max_model_structural_support_score_scaled=0.20, max_protein_evidence_score_scaled=0.20,
        max_transcript_evidence_score_scaled=0.40, max_model_selection_score_scaled=0.20,
        pruning_level_name="mid",
    )
    assert out["retained_models"] == []
    removed = out["removed_rows"][0]
    assert removed["final_low_score_nonsingleton_prune_route"] == "evidence_aware_transcript_cdna_context"
    assert "transcript_cdna_context_nonsingleton" in removed["final_prune_classes"]
    assert out["low_score_nonsingleton_removed_context_route"] == 1
