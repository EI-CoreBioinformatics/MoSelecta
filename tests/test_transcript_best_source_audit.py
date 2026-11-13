from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from moselecta.models import CodingModel
from moselecta.scoring import apply_model_selection_scores
from moselecta.transcript_evidence import (
    TranscriptEvidenceCollapsedTranscript,
    TranscriptEvidenceSuperlocus,
    _best_per_source_pairwise_components,
    calculate_transcript_model_support,
    transcript_model_pair_metrics,
)


COMPONENT_FIELDS = {
    "mean_best_source_model_cds_intron_support": "pairwise_model_cds_intron_support",
    "mean_best_source_splice_site_f1": "pairwise_splice_site_f1",
    "mean_best_source_full_intron_chain_f1": "pairwise_full_intron_chain_f1",
    "mean_best_source_internal_cds_exon_compatibility": "pairwise_internal_cds_exon_compatibility",
    "mean_best_source_terminal_cds_exon_compatibility": "pairwise_terminal_cds_exon_compatibility",
    "mean_best_source_model_cds_exonic_coverage": "pairwise_model_cds_exonic_coverage",
}


def make_model(parts: tuple[tuple[int, int], ...], model_id: str = "model") -> CodingModel:
    return CodingModel(
        source_name="model_source",
        input_file=Path("model.gff3"),
        input_gene_id=f"{model_id}_gene",
        input_transcript_id=model_id,
        output_gene_id=f"{model_id}_gene",
        output_transcript_id=model_id,
        chrom="chr1",
        strand="+",
        cds_parts=list(parts),
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def make_transcript(
    parts: tuple[tuple[int, int], ...],
    transcript_id: str,
    sources: tuple[str, ...],
) -> TranscriptEvidenceCollapsedTranscript:
    return TranscriptEvidenceCollapsedTranscript(
        transcript_evidence_id=transcript_id,
        chrom="chr1",
        strand="+",
        exon_parts=parts,
        evidence_superlocus_id="transcript_superlocus",
        collapse_key_type="exact_intron_chain",
        source_member_counts={source: 1 for source in sources},
        member_count=len(sources),
        span_start0=parts[0][0],
        span_end0=parts[-1][1],
        left_terminal_outer_min0=parts[0][0],
        left_terminal_outer_max0=parts[0][0],
        right_terminal_outer_min0=parts[-1][1],
        right_terminal_outer_max0=parts[-1][1],
    )


def make_superlocus(
    transcripts: list[TranscriptEvidenceCollapsedTranscript],
) -> TranscriptEvidenceSuperlocus:
    return TranscriptEvidenceSuperlocus(
        evidence_superlocus_id="transcript_superlocus",
        chrom="chr1",
        strand="+",
        start0=min(tx.start0 for tx in transcripts),
        end0=max(tx.end0 for tx in transcripts),
        collapsed_transcripts=transcripts,
    )


def test_one_source_mean_best_audit_fields_match_pairwise_components() -> None:
    parts = ((100, 150), (250, 300), (400, 450))
    model = make_model(parts)
    tx = make_transcript(parts, "tx_exact", ("source_a",))
    pair_metrics = transcript_model_pair_metrics(model, tx)

    result = _best_per_source_pairwise_components([(tx, pair_metrics)])

    assert result["source_count"] == 1
    assert result["mean_best_source_pairwise_score"] == pytest.approx(
        pair_metrics["pairwise_transcript_score"]
    )
    for output_key, pairwise_key in COMPONENT_FIELDS.items():
        assert result[output_key] == pytest.approx(pair_metrics[pairwise_key])


def test_multiple_sources_use_each_sources_best_pair_for_every_component() -> None:
    model = make_model(((100, 150), (250, 300), (400, 450)))
    tx_a_low = make_transcript(
        ((100, 150), (250, 300), (420, 450)), "tx_a_low", ("source_a",)
    )
    tx_a_best = make_transcript(
        ((100, 150), (250, 300), (400, 450)), "tx_a_best", ("source_a",)
    )
    tx_b_best = make_transcript(
        ((100, 150), (260, 300), (400, 450)), "tx_b_best", ("source_b",)
    )
    associations = [
        (tx, transcript_model_pair_metrics(model, tx))
        for tx in (tx_a_low, tx_a_best, tx_b_best)
    ]
    selected_metrics = [associations[1][1], associations[2][1]]

    result = _best_per_source_pairwise_components(associations)

    assert associations[1][1]["pairwise_transcript_score"] > associations[0][1]["pairwise_transcript_score"]
    assert result["source_count"] == 2
    assert result["mean_best_source_pairwise_score"] == pytest.approx(
        sum(float(metrics["pairwise_transcript_score"]) for metrics in selected_metrics) / 2
    )
    for output_key, pairwise_key in COMPONENT_FIELDS.items():
        expected = sum(float(metrics[pairwise_key]) for metrics in selected_metrics) / 2
        assert result[output_key] == pytest.approx(expected)
    assert result["mean_best_source_splice_site_f1"] > 0.0
    assert result["mean_best_source_full_intron_chain_f1"] > 0.0
    assert result["mean_best_source_model_cds_intron_support"] > 0.0


def test_internal_component_excludes_na_but_preserves_true_zero() -> None:
    two_exon_model = make_model(((100, 150), (250, 300)), "two_exon")
    two_exon_tx_a = make_transcript(
        ((100, 150), (250, 300)), "two_a", ("source_a",)
    )
    two_exon_tx_b = make_transcript(
        ((100, 150), (250, 300)), "two_b", ("source_b",)
    )
    two_result = _best_per_source_pairwise_components([
        (tx, transcript_model_pair_metrics(two_exon_model, tx))
        for tx in (two_exon_tx_a, two_exon_tx_b)
    ])

    assert two_result["source_count"] == 2
    assert two_result["mean_best_source_internal_cds_exon_compatibility"] is None

    three_exon_model = make_model(
        ((1000, 1050), (1150, 1200), (1300, 1350)), "three_exon"
    )
    unsupported_internal_tx = make_transcript(
        ((1000, 1050), (1150, 1190), (1300, 1350)), "three_zero", ("source_a",)
    )
    three_result = _best_per_source_pairwise_components([
        (
            unsupported_internal_tx,
            transcript_model_pair_metrics(three_exon_model, unsupported_internal_tx),
        )
    ])

    assert three_result["source_count"] == 1
    assert three_result["mean_best_source_internal_cds_exon_compatibility"] == pytest.approx(0.0)


def _legacy_best_source_components(scoring_associations):
    """Reproduce the pre-fix component lookup for score-invariance testing."""
    best_by_source = {}
    for tx, metrics in scoring_associations:
        score = float(metrics.get("pairwise_transcript_score") or 0.0)
        for source in tx.sources:
            previous = best_by_source.get(source)
            if previous is None or score > float(previous[1].get("pairwise_transcript_score") or 0.0):
                best_by_source[source] = (tx, metrics)
    best_items = list(best_by_source.values())
    output = {
        "source_count": len(best_by_source),
        "mean_best_source_pairwise_score": (
            sum(float(metrics.get("pairwise_transcript_score") or 0.0) for _tx, metrics in best_items)
            / len(best_items)
            if best_items
            else 0.0
        ),
    }
    for key in (
        "model_cds_intron_support",
        "splice_site_f1",
        "full_intron_chain_f1",
        "terminal_cds_exon_compatibility",
        "model_cds_exonic_coverage",
    ):
        output[f"mean_best_source_{key}"] = (
            sum(float(metrics.get(key) or 0.0) for _tx, metrics in best_items) / len(best_items)
            if best_items
            else 0.0
        )
    internal_values = [
        float(metrics["pairwise_internal_cds_exon_compatibility"])
        for _tx, metrics in best_items
        if metrics.get("pairwise_internal_cds_exon_compatibility") is not None
    ]
    output["mean_best_source_internal_cds_exon_compatibility"] = (
        sum(internal_values) / len(internal_values) if internal_values else None
    )
    return output


def _selection_weights() -> dict[str, float]:
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
        "transcript_evidence": 2.0,
    }


def _selection_row(model_id: str) -> dict[str, object]:
    return {
        "sublocus_id": "sub001",
        "source_name": model_id,
        "output_transcript_id": model_id,
        "cds_length": 150,
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


def test_audit_component_fix_preserves_scores_best_source_and_selection(monkeypatch) -> None:
    import moselecta.transcript_evidence as transcript_evidence

    parts = ((100, 150), (250, 300), (400, 450))
    model = make_model(parts, "supported")
    tx_a = make_transcript(parts, "tx_a", ("source_a",))
    tx_b = make_transcript(
        ((100, 150), (260, 300), (400, 450)), "tx_b", ("source_b",)
    )
    superlocus = make_superlocus([tx_a, tx_b])

    fixed_row = calculate_transcript_model_support(
        [model], [superlocus], ["source_a", "source_b"]
    )[0]
    monkeypatch.setattr(
        transcript_evidence,
        "_best_per_source_pairwise_components",
        _legacy_best_source_components,
    )
    legacy_row = calculate_transcript_model_support(
        [model], [superlocus], ["source_a", "source_b"]
    )[0]

    invariant_fields = (
        "transcript_best_pairwise_id",
        "transcript_best_pairwise_score",
        "transcript_mean_best_source_pairwise_score",
        "transcript_weighted_mean_scoring_eligible_pairwise_score",
        "transcript_supporting_source_count",
        "transcript_supporting_source_fraction",
        "transcript_supporting_sources",
        "transcript_chain_source_score",
        "transcript_support_confidence",
        "transcript_evidence_score",
        "transcript_evidence_score_scaled",
        "transcript_evidence_score_type",
        "transcript_evidence_score_components",
        "transcript_evidence_status",
        "transcript_evidence_flags",
    )
    assert {key: fixed_row[key] for key in invariant_fields} == {
        key: legacy_row[key] for key in invariant_fields
    }
    for field in (
        "transcript_mean_best_source_model_cds_intron_support",
        "transcript_mean_best_source_full_intron_chain_f1",
        "transcript_mean_best_source_terminal_cds_exon_compatibility",
        "transcript_mean_best_source_model_cds_exonic_coverage",
    ):
        assert legacy_row[field] == 0.0
        assert fixed_row[field] > 0.0
    assert legacy_row["transcript_mean_best_source_splice_site_f1"] == pytest.approx(
        fixed_row["transcript_mean_best_source_splice_site_f1"]
    )
    assert fixed_row["transcript_mean_best_source_splice_site_f1"] > 0.0

    fixed_evidence = {
        "supported": fixed_row,
        "unsupported": {
            "transcript_evidence_score": 0.0,
            "transcript_evidence_score_scaled": 0.0,
        },
    }
    legacy_evidence = {
        "supported": legacy_row,
        "unsupported": {
            "transcript_evidence_score": 0.0,
            "transcript_evidence_score_scaled": 0.0,
        },
    }
    fixed_selection = [_selection_row("supported"), _selection_row("unsupported")]
    legacy_selection = deepcopy(fixed_selection)
    apply_model_selection_scores(
        fixed_selection,
        source_prior_scores={},
        selection_weights=_selection_weights(),
        tie_epsilon=1e-12,
        transcript_evidence_rows_by_transcript=fixed_evidence,
        transcript_evidence_active=True,
    )
    apply_model_selection_scores(
        legacy_selection,
        source_prior_scores={},
        selection_weights=_selection_weights(),
        tie_epsilon=1e-12,
        transcript_evidence_rows_by_transcript=legacy_evidence,
        transcript_evidence_active=True,
    )
    assert fixed_selection == legacy_selection
