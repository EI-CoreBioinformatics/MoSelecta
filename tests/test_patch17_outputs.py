from __future__ import annotations

from moselecta.outputs import model_gff_scoring_attrs
from moselecta.transcript_evidence import transcript_summary_tsv_rows


def test_transcript_summary_includes_model_association_and_context_metrics() -> None:
    rows = transcript_summary_tsv_rows(
        transcript_paths=[],
        transcript_source_labels=[],
        raw_transcripts=[],
        superloci=[],
        collapsed_transcripts=[],
        subloci=[],
        raw_to_collapsed_mapping_requested=False,
        model_support_summary={
            "model_support_row_count": 12,
            "models_with_scoring_eligible_transcript_evidence": 8,
            "models_with_cdna_context_calculated": 5,
            "models_with_cdna_context_pruning_concern": 2,
            "model_status_counts": "strong_transcript_support:4",
        },
    )
    values = {
        (row["summary_scope"], row["metric"]): row["value"]
        for row in rows
    }
    assert values[("model_support", "model_support_row_count")] == 12
    assert values[("model_support", "models_with_scoring_eligible_transcript_evidence")] == 8
    assert values[("model_support", "models_with_cdna_context_calculated")] == 5
    assert values[("model_support", "models_with_cdna_context_pruning_concern")] == 2
    assert values[("model_support", "model_status_counts")] == "strong_transcript_support:4"


def test_gff_scoring_attrs_include_compact_transcript_fields() -> None:
    attrs = model_gff_scoring_attrs({
        "transcript_evidence_status": "strong_transcript_support",
        "transcript_evidence_score": 0.9,
        "transcript_evidence_score_scaled": 0.9,
        "transcript_evidence_score_context_adjusted": 0.8,
        "transcript_evidence_score_context_adjusted_scaled": 0.8,
        "transcript_evidence_score_weight": 2.0,
        "transcript_evidence_score_in_selection": "true",
        "transcript_cdna_context_status": "calculated",
        "transcript_cdna_context_cds_to_min_cdna_fraction": 0.12,
        "transcript_cdna_context_concern_class": "severe_low_fraction",
        "transcript_cdna_context_pruning_concern": "true",
    })
    assert attrs["transcript_evidence_score"] == "0.900000"
    assert attrs["transcript_evidence_score_context_adjusted"] == "0.800000"
    assert attrs["transcript_evidence_status"] == "strong_transcript_support"
    assert attrs["transcript_cdna_context_cds_to_min_cdna_fraction"] == "0.120000"
    assert attrs["transcript_cdna_context_pruning_concern"] == "true"


def test_gff_scoring_attrs_omit_transcript_fields_when_not_configured() -> None:
    attrs = model_gff_scoring_attrs({
        "transcript_evidence_status": "transcript_evidence_not_configured",
        "transcript_evidence_score": 0.0,
    })
    assert not any(key.startswith("transcript_") for key in attrs)
