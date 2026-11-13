from __future__ import annotations

from pathlib import Path

from moselecta.models import CodingModel, Sublocus
from moselecta.selection import (
    final_prune_low_score_singletons,
    final_prune_singleton_selected_models,
)


def _model(transcript_id: str, source: str) -> CodingModel:
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"g_{transcript_id}",
        input_transcript_id=transcript_id,
        output_gene_id=f"out_g_{transcript_id}",
        output_transcript_id=transcript_id,
        chrom="chr1",
        strand="+",
        cds_parts=[(100, 250)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def _stage4(candidate: CodingModel) -> dict[str, object]:
    return {
        "monosublocus_rows": [{"monosublocus_id": "mono1", "member_count": 2}],
        "monosublocus_member_rows": [
            {"monosublocus_id": "mono1", "output_transcript_id": candidate.output_transcript_id},
            {"monosublocus_id": "mono1", "output_transcript_id": "partner"},
        ],
        "monosublocus_edge_rows": [
            {
                "monosublocus_id": "mono1",
                "model_a_id": candidate.output_transcript_id,
                "model_b_id": "partner",
                "source_a": candidate.source_name,
                "source_b": "Annevo",
                "cds_exon_overlap_bp": 150,
                "overlap_fraction_of_a": 1.0,
                "overlap_fraction_of_b": 0.85,
            }
        ],
    }


def _prune(
    *,
    protein_adequacy: bool,
    transcript_adequacy: bool,
    protein_score: float = 0.0,
    transcript_score: float = 0.0,
    transcript_fraction: float = 0.50,
    policy: str = "adequacy-conditioned",
    sequence_integrity_flags: str = ".",
    allow_internal_stop_evidence_rescue: bool = False,
) -> dict[str, object]:
    model = _model("candidate", "Galba")
    sub = Sublocus("sub1", "sl1", "chr1", "+", 100, 250, [model])
    row = {
        "output_transcript_id": model.output_transcript_id,
        "source_name": model.source_name,
        "final_prune_classes": ".",
    }
    return final_prune_singleton_selected_models(
        {
            "initial_retained_models": [model],
            "initial_retained_rows": [row],
            "retained_models": [model],
            "retained_rows": [row],
            "removed_rows": [],
            "retained_selection_info_by_transcript": {model.output_transcript_id: row},
        },
        [sub],
        {model.output_transcript_id: {
            "intrinsic_penalty_components": ".",
            "intrinsic_penalty_score": 0.0,
            "sequence_integrity_flags": sequence_integrity_flags,
        }},
        prune_level="all",
        prune_sources=set(),
        pruning_level_name="high",
        rescue_monosublocus_overlap=True,
        stage4_result=_stage4(model),
        monosublocus_overlap_rescue_evidence_policy=policy,
        protein_adequacy_passed=protein_adequacy,
        transcript_adequacy_passed=transcript_adequacy,
        rescue_protein_supported_singletons=True,
        protein_evidence_rows_by_transcript={
            model.output_transcript_id: {
                "protein_evidence_score_scaled": protein_score,
                "protein_exact_cds_chain_source_count": 0,
            }
        },
        rescue_transcript_supported_singletons=True,
        transcript_evidence_rows_by_transcript={
            model.output_transcript_id: {
                "transcript_evidence_score_scaled": transcript_score,
                "transcript_supporting_source_count": 1,
                "transcript_scoring_eligible_count": 1,
                "transcript_stranded_scoring_eligible_count": 1,
                "transcript_unstranded_scoring_eligible_count": 0,
                "transcript_cdna_context_evaluable": "true",
                "transcript_cdna_context_cds_to_min_cdna_fraction": transcript_fraction,
                "transcript_cdna_context_pruning_concern": "false",
            }
        },
        allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
    )


def test_overlap_only_rescue_is_allowed_when_no_evidence_class_is_adequate() -> None:
    result = _prune(protein_adequacy=False, transcript_adequacy=False)
    assert len(result["retained_models"]) == 1
    row = result["retained_rows"][0]
    assert row["final_singleton_overlap_rescue_decision"] == "rescued_no_adequate_evidence"
    assert row["final_singleton_prune_reason"] == "singleton_sublocus_rescued_by_monosublocus_overlap"


def test_overlap_only_rescue_is_denied_when_adequate_evidence_lacks_strong_local_support() -> None:
    result = _prune(protein_adequacy=True, transcript_adequacy=False, protein_score=0.40)
    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_singleton_overlap_rescue_geometry_passed"] == "true"
    assert row["final_singleton_overlap_rescue_decision"] == "denied"
    assert row["final_singleton_overlap_rescue_denial_reason"] == "adequate_evidence_without_strong_local_support"


def test_strong_local_support_rescues_even_when_its_evidence_class_is_not_adequate() -> None:
    result = _prune(
        protein_adequacy=True,
        transcript_adequacy=False,
        protein_score=0.40,
        transcript_score=0.82,
    )
    assert len(result["retained_models"]) == 1
    row = result["retained_rows"][0]
    assert row["final_singleton_overlap_rescue_strong_local_transcript"] == "true"
    assert row["final_singleton_overlap_rescue_decision"] == "rescued_strong_local_evidence"


def test_internal_stop_blocks_group_c_strong_local_evidence_rescue_by_default() -> None:
    result = _prune(
        protein_adequacy=True,
        transcript_adequacy=False,
        transcript_score=0.82,
        sequence_integrity_flags="internal_stop_codon",
    )
    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_internal_stop_evidence_criteria_met"] == "true"
    assert row["final_internal_stop_evidence_rescue_blocked"] == "true"
    assert row["final_singleton_overlap_rescue_strong_local_transcript"] == "false"
    assert row["final_singleton_overlap_rescue_decision"] == "denied"


def test_internal_stop_override_restores_group_c_strong_local_evidence_rescue() -> None:
    result = _prune(
        protein_adequacy=True,
        transcript_adequacy=False,
        transcript_score=0.82,
        sequence_integrity_flags="internal_stop_codon",
        allow_internal_stop_evidence_rescue=True,
    )
    assert len(result["retained_models"]) == 1
    row = result["retained_rows"][0]
    assert row["final_singleton_overlap_rescue_strong_local_transcript"] == "true"
    assert row["final_internal_stop_evidence_rescue_applied"] == "true"
    assert (
        row["final_internal_stop_evidence_rescue_applied_route"]
        == "monosublocus_overlap_transcript_evidence"
    )


def test_unconditional_policy_preserves_previous_geometry_only_rescue() -> None:
    result = _prune(
        protein_adequacy=True,
        transcript_adequacy=True,
        policy="unconditional",
    )
    assert len(result["retained_models"]) == 1
    assert result["retained_rows"][0]["final_singleton_overlap_rescue_decision"] == "rescued_unconditional"


def test_overlap_policy_summary_counts_survive_low_score_singleton_stage() -> None:
    result = _prune(protein_adequacy=False, transcript_adequacy=False)
    model = result["retained_models"][0]
    sub = Sublocus("sub1", "sl1", "chr1", "+", 100, 250, [model])
    chained = final_prune_low_score_singletons(
        result,
        [sub],
        {model.output_transcript_id: {}},
        enabled=True,
        protein_evidence_supplied=False,
        transcript_evidence_supplied=False,
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
        max_model_structural_support_score_scaled=0.25,
        max_protein_evidence_score_scaled=0.25,
        max_transcript_evidence_score_scaled=0.40,
        max_model_selection_score_scaled=0.25,
        pruning_level_name="medium",
    )
    assert chained["singleton_overlap_rescue_geometry_candidate_count"] == 1
    assert chained["singleton_overlap_rescue_overlap_bp_geometry_candidate_count"] == 1
    assert chained["singleton_overlap_rescue_reciprocal_rejected_count"] == 0
    assert chained["singleton_overlap_rescue_models_with_same_source_edge_count"] == 0
    assert chained["singleton_overlap_rescued_without_adequate_evidence_count"] == 1
    assert chained["singleton_overlap_rescue_policy_denied_count"] == 0
