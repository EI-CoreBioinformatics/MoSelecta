from __future__ import annotations

from pathlib import Path

from moselecta.auto_pruning import (
    ScenarioSummary,
    add_incremental_metrics,
    apply_monotonic_pruning_closure,
    removal_provenance_from_result,
    validate_scenarios,
)
from moselecta.models import CodingModel, Sublocus
from moselecta.selection import final_prune_low_score_singletons


def _model(model_id: str, source: str = "A", start: int = 100) -> CodingModel:
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"g_{model_id}",
        input_transcript_id=model_id,
        output_gene_id=f"out_g_{model_id}",
        output_transcript_id=model_id,
        chrom="chr1",
        strand="+",
        cds_parts=[(start, start + 150)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def _low_score_kwargs() -> dict[str, object]:
    return {
        "enabled": True,
        "protein_evidence_supplied": True,
        "transcript_evidence_supplied": False,
        "protein_adequacy_required": True,
        "transcript_adequacy_required": True,
        "protein_overlap_model_fraction": 1.0,
        "protein_supported_model_fraction": 1.0,
        "transcript_scoring_eligible_model_fraction": 0.0,
        "transcript_supported_model_fraction": 0.0,
        "min_protein_overlap_model_fraction": 0.75,
        "min_protein_supported_model_fraction": 0.50,
        "min_transcript_scoring_eligible_model_fraction": 0.75,
        "min_transcript_supported_model_fraction": 0.50,
        "max_model_structural_support_score_scaled": 0.25,
        "max_protein_evidence_score_scaled": 0.25,
        "max_transcript_evidence_score_scaled": 0.40,
        "max_model_selection_score_scaled": 0.25,
        "pruning_level_name": "medium",
    }


def test_low_score_singleton_uses_protein_rescue_criteria_not_prior_route_application() -> None:
    model = _model("m1")
    sub = Sublocus("sub1", "sl1", "chr1", "+", 100, 250, [model])
    row = {
        "output_transcript_id": "m1",
        "source_name": "A",
        "model_structural_support_score_scaled": 0.10,
        "protein_evidence_score_scaled": 0.10,
        "transcript_evidence_score_scaled": 0.0,
        "model_selection_score_scaled": 0.10,
        "final_protein_singleton_rescue_criteria_met": "true",
        "final_protein_singleton_rescue_criteria_reason": "protein_exact_cds_chain_support",
        "final_protein_singleton_rescue_applied": "false",
        "final_transcript_singleton_rescue_criteria_met": "false",
        "final_transcript_singleton_rescue_applied": "false",
    }
    result = {
        "initial_retained_models": [model],
        "initial_retained_rows": [row],
        "retained_models": [model],
        "retained_rows": [row],
        "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }
    out = final_prune_low_score_singletons(result, [sub], {}, **_low_score_kwargs())
    assert out["retained_models"] == [model]
    retained = out["retained_rows"][0]
    assert retained["final_low_score_singleton_prune_reason"] == "protected_by_protein_evidence_rescue"
    assert retained["final_protein_singleton_rescue_applied"] == "true"
    assert retained["final_protein_singleton_rescued"] == "true"
    assert out["singleton_rescued_by_protein_evidence_count"] == 1


def test_low_score_singleton_uses_transcript_rescue_criteria_not_prior_route_application() -> None:
    model = _model("m1")
    sub = Sublocus("sub1", "sl1", "chr1", "+", 100, 250, [model])
    row = {
        "output_transcript_id": "m1",
        "source_name": "A",
        "model_structural_support_score_scaled": 0.10,
        "protein_evidence_score_scaled": 0.0,
        "transcript_evidence_score_scaled": 0.10,
        "model_selection_score_scaled": 0.10,
        "final_protein_singleton_rescue_criteria_met": "false",
        "final_protein_singleton_rescue_applied": "false",
        "final_transcript_singleton_rescue_criteria_met": "true",
        "final_transcript_singleton_rescue_criteria_reason": (
            "transcript_evidence_score_scaled+transcript_supporting_source_count+"
            "transcript_scoring_eligible_count+transcript_cdna_context_not_evaluable_neutral"
        ),
        "final_transcript_singleton_rescue_applied": "false",
    }
    kwargs = _low_score_kwargs()
    kwargs.update({
        "protein_evidence_supplied": False,
        "transcript_evidence_supplied": True,
        "protein_overlap_model_fraction": 0.0,
        "protein_supported_model_fraction": 0.0,
        "transcript_scoring_eligible_model_fraction": 1.0,
        "transcript_supported_model_fraction": 1.0,
    })
    result = {
        "initial_retained_models": [model],
        "initial_retained_rows": [row],
        "retained_models": [model],
        "retained_rows": [row],
        "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }
    out = final_prune_low_score_singletons(result, [sub], {}, **kwargs)
    assert out["retained_models"] == [model]
    retained = out["retained_rows"][0]
    assert retained["final_low_score_singleton_prune_reason"] == "protected_by_transcript_evidence_rescue"
    assert retained["final_transcript_singleton_rescue_applied"] == "true"
    assert retained["final_transcript_singleton_rescued"] == "true"
    assert out["singleton_rescued_by_transcript_evidence_count"] == 1


def test_monotonic_closure_inherits_shadow_provenance_without_claiming_native_shadow_removal() -> None:
    restored = _model("restored", source="B")
    lower_result = {
        "removed_rows": [{
            "output_transcript_id": "restored",
            "source_name": "B",
            "final_prune_level": "medium",
            "final_prune_classes": "opposite_strand_shadow",
            "final_shadow_pruned": "true",
        }]
    }
    high_result = {
        "retained_models": [restored],
        "retained_rows": [{
            "output_transcript_id": "restored",
            "source_name": "B",
            "final_prune_level": "high",
            "final_prune_classes": ".",
            "final_shadow_pruned": "false",
        }],
        "removed_rows": [{
            "output_transcript_id": "native",
            "source_name": "A",
            "final_prune_level": "high",
            "final_prune_classes": "singleton_sublocus",
        }],
        "retained_selection_info_by_transcript": {"restored": {}},
    }
    closed = apply_monotonic_pruning_closure(
        high_result,
        required_removed_ids={"restored"},
        lower_level="medium",
        current_level="high",
        lower_provenance_by_id=removal_provenance_from_result(lower_result, level="medium"),
    )
    assert closed["retained_models"] == []
    rows = {row["output_transcript_id"]: row for row in closed["removed_rows"]}
    inherited = rows["restored"]
    assert inherited["final_monotonic_inherited_removal"] == "true"
    assert inherited["final_monotonic_inherited_from_level"] == "medium"
    assert inherited["final_monotonic_inherited_original_primary_reason"] == "opposite_strand_shadow"
    assert inherited["final_shadow_pruned"] == "false"
    assert inherited["final_prune_classes"] == "inherited_lower_level_removal"


def _scenario(
    level: str,
    removed: set[str],
    *,
    inherited: set[str] | None = None,
    source_by_id: dict[str, str] | None = None,
) -> ScenarioSummary:
    all_ids = {"a", "b", "c"}
    inherited = inherited or set()
    source_by_id = source_by_id or {"a": "A", "b": "B", "c": "C"}
    return ScenarioSummary(
        level=level,
        starting_models=3,
        candidate_models=len(removed),
        protected_models=0,
        removed_models=len(removed),
        retained_models=3 - len(removed),
        removed_fraction=len(removed) / 3,
        affected_subloci=len(removed),
        retained_subloci=3 - len(removed),
        removed_ids=removed,
        retained_ids=all_ids - removed,
        primary_reason_by_id={model_id: "singleton_sublocus" for model_id in removed},
        all_triggered_by_id={model_id: ("singleton_sublocus",) for model_id in removed},
        native_removed_ids=removed - inherited,
        inherited_removed_ids=inherited,
        source_by_id=source_by_id,
        removal_by_source={},
        source_selected_counts={"A": 1, "B": 1, "C": 1},
    )


def test_incremental_source_counts_use_actual_model_set_and_validation_reports_native_closure() -> None:
    scenarios = {
        "medium": _scenario("medium", {"a"}),
        "high": _scenario("high", {"a", "b"}, inherited={"a"}),
        "very-high": _scenario("very-high", {"a", "b", "c"}),
    }
    add_incremental_metrics(scenarios)
    assert scenarios["high"].incremental_removed_ids == {"b"}
    assert scenarios["high"].incremental_by_source == {"B": 1}
    validation = validate_scenarios(scenarios)
    assert validation.passed
    assert validation.native_restored_medium_to_high == 1
