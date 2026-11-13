from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

from moselecta.alternative_splice import (
    ALTERNATIVE_SPLICE_CANDIDATE_COLUMNS,
    ALTERNATIVE_SPLICE_GENE_COLUMNS,
    ALTERNATIVE_SPLICE_RUN_COLUMNS,
    AlternativeSpliceReportPolicy,
    analyze_alternative_splice_report,
)
from moselecta.intervals import IntervalOverlapIndex
from moselecta.models import CodingModel, FeatureRecord, SequenceIntegrityMetrics
from moselecta.protein_evidence import ProteinEvidenceCollapsedAlignment
from moselecta.structural import (
    cds_exon_chain_containment_relation,
    cds_exon_chains_match,
    cds_frame_concordance,
    cds_intron_chains_match,
)
from moselecta.transcript_evidence import (
    TranscriptEvidenceCollapsedTranscript,
    TranscriptEvidenceSuperlocus,
    associated_transcripts_for_model,
    evaluate_coherent_exact_transcript_support,
)


def _record(
    start0: int,
    end0: int,
    phase: str,
    *,
    strand: str,
    feature_type: str = "CDS",
) -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1", source="test", feature_type=feature_type,
        start=start0 + 1, end=end0, score=".", strand=strand, phase=phase,
        attrs={"Parent": "tx"}, raw_attrs="Parent=tx", input_file=Path("test.gff3"),
        source_name="test", line_number=1, record_index=1, raw_line="",
    )


def _model(
    name: str,
    parts: list[tuple[int, int]],
    *,
    source: str = "A",
    strand: str = "+",
    phases: list[str] | None = None,
    superlocus: str = "super1",
    sublocus: str = "sub1",
    non_cds_parts: list[tuple[str, int, int]] | None = None,
) -> CodingModel:
    phase_values = phases if phases is not None else ["0"] * len(parts)
    children = [
        _record(start0, end0, phase, strand=strand)
        for (start0, end0), phase in zip(parts, phase_values)
    ]
    for feature_type, start0, end0 in non_cds_parts or []:
        children.append(_record(start0, end0, ".", strand=strand, feature_type=feature_type))
    model = CodingModel(
        source_name=source, input_file=Path(f"{source}.gff3"), input_gene_id=f"g_{name}",
        input_transcript_id=f"input_{name}", output_gene_id=f"gene_{name}",
        output_transcript_id=name, chrom="chr1", strand=strand, cds_parts=parts,
        gene_record=None, transcript_record=None, child_records=children,
    )
    model.superlocus_id = superlocus
    model.sublocus_id = sublocus
    return model


def _score_row(model: CodingModel, score: float = 1.0, **overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "output_transcript_id": model.output_transcript_id,
        "model_selection_score": score,
        "model_selection_score_scaled": max(0.0, min(score, 1.0)),
        "model_selection_order_rank": 1,
        "model_completeness_score": 1.0,
        "intrinsic_penalty_score": 0.0,
        "consensus_similarity_score": 1.0,
        "structural_support_score": 1.0,
        "source_prior_score": 1.0,
        "cds_length": model.cds_length,
        "model_structural_support_score": 1.0,
        "model_structural_support_score_scaled": 1.0,
        "pairwise_similarity_score": 1.0,
        "protein_evidence_score_scaled": 0.0,
        "transcript_evidence_score_scaled": 0.0,
        "cds_introns_total": len(model.cds_introns),
        "cds_introns_supported_by_other_source": 0,
        "cds_intron_support_source_min": 1 if model.cds_introns else ".",
        "cds_splice_sites_total": len(model.cds_splice_sites),
        "cds_splice_sites_supported_by_other_source": 0,
        "cds_splice_site_support_source_min": (
            1 if model.cds_splice_sites else "."
        ),
        "internal_cds_exons_total": len(model.internal_cds_exons),
        "internal_cds_exons_supported_by_other_source": 0,
        "internal_cds_exon_support_source_min": (
            1 if model.internal_cds_exons else "."
        ),
        "terminal_cds_exons_total": len(model.terminal_cds_exons),
        "terminal_cds_exons_supported_by_other_source": 0,
        "terminal_cds_exon_support_source_min": (
            1 if model.terminal_cds_exons else "."
        ),
        "cds_base_consensus_precision": 0.0,
        "intrinsic_penalty_cds_length_not_multiple_of_3_count": 0,
        "intrinsic_penalty_short_cds_count": 0,
        "intrinsic_penalty_short_cds_intron_count": 0,
        "intrinsic_penalty_short_internal_cds_exon_count": 0,
        "intrinsic_penalty_short_terminal_cds_exon_count": 0,
        "intrinsic_penalty_long_cds_intron_count": 0,
    }
    row.update(overrides)
    return row


def _policy(**overrides: object) -> AlternativeSpliceReportPolicy:
    values: dict[str, object] = {
        "terminal_cds_3prime_wobble_bp": 3,
        "protein_terminal_cds_3prime_wobble_bp": 3,
        "transcript_min_cdna_context_fraction": 0.30,
        "transcript_require_stranded_support": False,
        "transcript_block_on_cdna_context_concern": True,
        "transcript_cdna_context_neutral_fraction": 0.30,
        "transcript_cdna_context_severe_fraction": 0.15,
    }
    values.update(overrides)
    return AlternativeSpliceReportPolicy(**values)  # type: ignore[arg-type]


def _analyse(
    anchors: list[CodingModel],
    candidates: list[CodingModel],
    *,
    rows: dict[str, dict[str, object]] | None = None,
    **kwargs: object,
) -> dict[str, list[dict[str, object]]]:
    models = anchors + candidates
    score_rows = rows or {model.output_transcript_id: _score_row(model) for model in models}
    return analyze_alternative_splice_report(
        anchor_models=anchors,
        all_scored_models=models,
        model_score_rows_by_transcript=score_rows,
        policy=kwargs.pop("policy", _policy()),  # type: ignore[arg-type]
        **kwargs,  # type: ignore[arg-type]
    )


def test_multi_exon_same_sublocus_unique_ownership() -> None:
    primary = _model("primary", [(0, 90), (180, 270)])
    candidate = _model("candidate", [(0, 90), (180, 240), (300, 360)], source="B")
    row = _analyse([primary], [candidate])["candidate_rows"][0]
    assert row["ownership_status"] == "unique_owner"
    assert row["structural_association_route"] == "multi_exon_same_sublocus_shared_cds_splice_site"
    assert row["shared_cds_splice_site_count"] == 2


def test_model_source_depth_counts_distinct_prediction_sources_only() -> None:
    primary = _model(
        "primary", [(0, 90), (180, 270), (360, 450)], source="P"
    )
    candidate = _model(
        "candidate", [(0, 90), (180, 240), (330, 450)], source="A"
    )
    same_b1 = _model("same_b1", list(candidate.cds_parts), source="B")
    same_b2 = _model("same_b2", list(candidate.cds_parts), source="B")
    same_c = _model("same_c", list(candidate.cds_parts), source="C")
    report = _analyse(
        [primary],
        [candidate, same_b1, same_b2, same_c],
        run_metadata={
            "protein_evidence_source_count": 17,
            "transcript_evidence_source_count": 23,
        },
    )
    row = next(
        item for item in report["candidate_rows"]
        if item["candidate_output_transcript_id"] == "candidate"
    )
    assert row["alternative_splice_exact_model_source_count"] == 3
    assert row["alternative_splice_exact_model_source_fraction"] == 0.75
    assert report["run_rows"][0]["prediction_source_count"] == 4


def test_mosaic_minimum_source_depth_allows_different_sources_by_component() -> None:
    primary = _model(
        "primary", [(0, 90), (180, 270), (360, 450)], source="P"
    )
    candidate = _model(
        "candidate", [(0, 90), (180, 240), (330, 450)], source="A"
    )
    support_b = _model("support_b", [(0, 90), (180, 240)], source="B")
    support_c = _model("support_c", [(180, 240), (330, 450)], source="C")
    support_d = _model("support_d", [(0, 90), (330, 450)], source="D")
    models = [primary, candidate, support_b, support_c, support_d]
    rows = {model.output_transcript_id: _score_row(model) for model in models}
    rows["candidate"].update({
        "cds_introns_supported_by_other_source": len(candidate.cds_introns),
        "cds_intron_support_source_min": 3,
        "cds_splice_sites_supported_by_other_source": len(candidate.cds_splice_sites),
        "cds_splice_site_support_source_min": 3,
        "internal_cds_exons_supported_by_other_source": len(candidate.internal_cds_exons),
        "internal_cds_exon_support_source_min": 3,
        "terminal_cds_exons_supported_by_other_source": len(candidate.terminal_cds_exons),
        "terminal_cds_exon_support_source_min": 3,
        "cds_base_consensus_precision": 1.0,
    })
    report = _analyse(
        [primary],
        [candidate, support_b, support_c, support_d],
        rows=rows,
    )
    row = next(
        item for item in report["candidate_rows"]
        if item["candidate_output_transcript_id"] == "candidate"
    )
    assert row["alternative_splice_exact_model_source_count"] == 1
    assert row["alternative_splice_mosaic_min_component_source_count"] == 3
    assert row["alternative_splice_mosaic_min_component_source_fraction"] == 0.6
    assert row["full_structural_mosaic"] is True


def test_transitive_sublocus_with_multiple_primaries_is_ambiguous() -> None:
    left = _model("left", [(0, 90), (180, 270)], source="A")
    middle = _model("middle", [(0, 90), (180, 240), (330, 420)], source="B")
    right = _model("right", [(180, 240), (330, 420)], source="C")
    result = _analyse([left, right], [middle])
    rows = result["candidate_rows"]
    assert len(rows) == 2
    assert {row["primary_output_transcript_id"] for row in rows} == {"left", "right"}
    assert {row["ownership_status"] for row in rows} == {"ambiguous_multiple_owners"}
    assert {row["fixed_rule_potential_candidate"] for row in rows} == {False}


@pytest.mark.parametrize(
    ("candidate_parts", "expected_owner"),
    [
        ([(0, 50), (100, 150)], True),
        ([(0, 50), (100, 149)], False),
    ],
)
def test_single_exon_primary_two_exon_gate1_threshold(
    candidate_parts: list[tuple[int, int]], expected_owner: bool
) -> None:
    primary = _model("primary", [(50, 150)], sublocus="single")
    candidate = _model("candidate", candidate_parts, source="B", sublocus="multi")
    row = _analyse([primary], [candidate])["candidate_rows"][0]
    assert (row["ownership_status"] == "unique_owner") is expected_owner
    if expected_owner:
        assert row["primary_cds_overlap_fraction"] == pytest.approx(0.50)
        assert row["candidate_cds_overlap_fraction"] == pytest.approx(0.50)


def test_single_exon_primary_route_still_checks_other_retained_primaries() -> None:
    single = _model("single", [(50, 150)], sublocus="single")
    multi = _model("multi", [(0, 50), (100, 180)], source="C", sublocus="multi")
    candidate = _model("candidate", [(0, 50), (100, 150)], source="B", sublocus="multi")
    rows = _analyse([single, multi], [candidate])["candidate_rows"]
    assert len(rows) == 2
    assert {row["ownership_status"] for row in rows} == {"ambiguous_multiple_owners"}


@pytest.mark.parametrize(
    "parts",
    [
        [(50, 150)],
        [(0, 40), (60, 100), (120, 160)],
    ],
)
def test_single_exon_primary_route_requires_exactly_two_candidate_cds_exons(
    parts: list[tuple[int, int]],
) -> None:
    primary = _model("primary", [(50, 150)], sublocus="single")
    candidate = _model("candidate", parts, source="B", sublocus="other")
    row = _analyse([primary], [candidate])["candidate_rows"][0]
    assert row["ownership_status"] == "zero_owner"


def test_same_sublocus_without_direct_shared_splice_has_no_owner() -> None:
    primary = _model("primary", [(0, 90), (180, 270)])
    candidate = _model("candidate", [(20, 110), (200, 290)], source="B")
    row = _analyse([primary], [candidate])["candidate_rows"][0]
    assert row["ownership_status"] == "zero_owner"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_reciprocal_overlap_is_strand_independent(strand: str) -> None:
    primary = _model("primary", [(50, 150)], strand=strand, sublocus="single")
    candidate = _model("candidate", [(0, 50), (100, 150)], strand=strand, source="B", sublocus="multi")
    row = _analyse([primary], [candidate])["candidate_rows"][0]
    assert row["primary_cds_overlap_fraction"] == pytest.approx(0.50)
    assert row["candidate_cds_overlap_fraction"] == pytest.approx(0.50)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_frame_concordant_and_shifted_models(strand: str) -> None:
    primary = _model("primary", [(0, 90), (180, 270)], strand=strand, phases=["0", "0"])
    same = _model("same", [(0, 90), (180, 270)], strand=strand, phases=["0", "0"])
    shifted = _model("shifted", [(0, 90), (180, 270)], strand=strand, phases=["1", "1"])
    assert cds_frame_concordance(primary, same)["same_frame_overlap_fraction"] == 1.0
    shifted_metrics = cds_frame_concordance(primary, shifted)
    assert shifted_metrics["frame_concordance_status"] == "calculated"
    assert shifted_metrics["same_frame_overlap_fraction"] == 0.0


def test_missing_and_partial_phase_are_explicit() -> None:
    primary = _model("primary", [(0, 90), (180, 270)], phases=["0", "0"])
    missing = _model("missing", [(0, 90), (180, 270)], phases=[".", "."])
    partial = _model("partial", [(0, 90), (180, 270)], phases=["0", "."])
    assert cds_frame_concordance(primary, missing)["frame_concordance_status"] == "unavailable_missing_phase"
    partial_metrics = cds_frame_concordance(primary, partial)
    assert partial_metrics["frame_concordance_status"] == "partially_available"
    assert partial_metrics["same_frame_overlap_fraction"] is None
    assert partial_metrics["frame_comparable_overlap_bp"] == 90


def test_missing_first_phase_is_not_silently_assumed_zero() -> None:
    primary = _model("primary", [(0, 90), (180, 270)], phases=["1", "1"])
    partial = _model("partial", [(0, 90), (180, 270)], phases=[".", "1"])
    metrics = cds_frame_concordance(primary, partial)
    assert metrics["candidate_phase_status"] == "partially_available"
    assert metrics["candidate_phase_mismatch_count"] == 0
    assert metrics["frame_concordance_status"] == "partially_available"


def test_utr_only_difference_does_not_create_distinct_cds() -> None:
    primary = _model("primary", [(0, 90), (180, 270)], non_cds_parts=[("exon", 0, 270)])
    candidate = _model("candidate", [(0, 90), (180, 270)], source="B", non_cds_parts=[("exon", 0, 300)])
    row = _analyse([primary], [candidate])["candidate_rows"][0]
    assert row["exact_cds_equivalent_to_primary"] is True
    assert row["utr_only_difference_from_primary"] is True
    assert row["fixed_rule_potential_candidate"] is False


@pytest.mark.parametrize(
    ("strand", "left", "right"),
    [
        ("+", [(0, 90), (180, 270)], [(0, 90), (180, 273)]),
        ("-", [(0, 90), (180, 270)], [(-3, 90), (180, 270)]),
    ],
)
def test_terminal_3prime_wobble_both_strands(
    strand: str, left: list[tuple[int, int]], right: list[tuple[int, int]]
) -> None:
    assert cds_exon_chains_match(
        _model("left", left, strand=strand),
        _model("right", right, strand=strand),
        3,
    )


def test_terminal_wobble_is_not_transitive() -> None:
    a = _model("a", [(0, 90), (180, 270)])
    b = _model("b", [(0, 90), (180, 273)])
    c = _model("c", [(0, 90), (180, 276)])
    assert cds_exon_chains_match(a, b, 3)
    assert cds_exon_chains_match(b, c, 3)
    assert not cds_exon_chains_match(a, c, 3)


def test_exact_model_support_is_representative_centred_not_transitive() -> None:
    primary = _model("primary", [(0, 90), (180, 360)], source="P")
    representative = _model("representative", [(0, 90), (180, 270)], source="A")
    direct = _model("direct", [(0, 90), (180, 273)], source="B")
    transitive_only = _model("transitive_only", [(0, 90), (180, 276)], source="C")
    rows = _analyse([primary], [representative, direct, transitive_only])["candidate_rows"]
    representative_row = next(row for row in rows if row["candidate_output_transcript_id"] == "representative")
    assert representative_row["exact_model_supporting_model_ids"] == "direct"
    assert representative_row["exact_model_other_sources"] == "B"


def test_same_intron_chain_terminal_variant_is_redundant() -> None:
    primary = _model("primary", [(0, 90), (180, 270)])
    candidate = _model("candidate", [(-20, 90), (180, 300)], source="B")
    assert cds_intron_chains_match(primary, candidate)
    row = _analyse([primary], [candidate])["candidate_rows"][0]
    assert row["same_intron_chain_terminal_variant"] is True
    assert row["fixed_rule_potential_candidate"] is False


def test_bidirectional_true_cds_chain_containment() -> None:
    short = _model("short", [(100, 150), (200, 250)])
    long = _model("long", [(0, 50), (100, 150), (200, 250), (300, 350)])
    relation = cds_exon_chain_containment_relation(short, long)
    reverse = cds_exon_chain_containment_relation(long, short)
    assert relation == {"a_contained_in_b": True, "b_contained_in_a": False}
    assert reverse == {"a_contained_in_b": False, "b_contained_in_a": True}


def test_genomic_containment_without_true_cds_chain_containment() -> None:
    inner = _model("inner", [(100, 160), (200, 250)])
    outer = _model("outer", [(0, 50), (100, 150), (200, 250), (300, 350)])
    assert outer.cds_start0 <= inner.cds_start0 and inner.cds_end0 <= outer.cds_end0
    assert not any(cds_exon_chain_containment_relation(inner, outer).values())


def test_same_source_duplicates_do_not_satisfy_multisource_support() -> None:
    primary = _model("primary", [(0, 90), (180, 270)], source="P")
    candidate = _model("candidate", [(0, 90), (180, 240), (300, 360)], source="A")
    duplicate = _model("duplicate", list(candidate.cds_parts), source="A")
    corroborator = _model("corroborator", list(candidate.cds_parts), source="B")
    same_source = _analyse([primary], [candidate, duplicate])["candidate_rows"][0]
    assert same_source["exact_model_other_source_count"] == 0
    assert same_source["exact_multisource_model_support"] is False
    multi_source = _analyse([primary], [candidate, duplicate, corroborator])["candidate_rows"][0]
    assert multi_source["exact_model_other_source_count"] == 1
    assert multi_source["exact_multisource_model_support"] is True


def test_full_structural_mosaic_requires_every_applicable_raw_component() -> None:
    primary = _model("primary", [(0, 90), (180, 270)], source="P")
    candidate = _model("candidate", [(0, 90), (180, 240), (300, 360)], source="A")
    candidate_row = _score_row(
        candidate,
        cds_introns_supported_by_other_source=len(candidate.cds_introns),
        cds_splice_sites_supported_by_other_source=len(candidate.cds_splice_sites),
        internal_cds_exons_supported_by_other_source=len(candidate.internal_cds_exons),
        terminal_cds_exons_supported_by_other_source=len(candidate.terminal_cds_exons),
        cds_base_consensus_precision=1.0,
    )
    rows = {"primary": _score_row(primary), "candidate": candidate_row}
    complete = _analyse([primary], [candidate], rows=rows)["candidate_rows"][0]
    assert complete["full_structural_mosaic"] is True
    assert complete["mosaic_cds_introns_status"] == "complete"
    candidate_row["terminal_cds_exons_supported_by_other_source"] = 1
    incomplete = _analyse([primary], [candidate], rows=rows)["candidate_rows"][0]
    assert incomplete["full_structural_mosaic"] is False
    assert incomplete["mosaic_terminal_cds_exons_status"] == "incomplete"


def test_exact_equivalent_candidates_keep_model_specific_quality_states() -> None:
    primary = _model("primary", [(0, 90), (180, 270)], source="P")
    first = _model("first", [(0, 90), (180, 240), (300, 360)], source="A")
    second = _model("second", list(first.cds_parts), source="B")
    rows = {model.output_transcript_id: _score_row(model) for model in [primary, first, second]}
    rows["second"]["intrinsic_penalty_short_cds_count"] = 1
    result = _analyse([primary], [first, second], rows=rows)
    by_id = {row["candidate_output_transcript_id"]: row for row in result["candidate_rows"]}
    assert by_id["first"]["intrinsic_default_exclusion"] is False
    assert by_id["second"]["intrinsic_default_exclusion"] is True


def test_local_exact_protein_support_ignores_global_adequacy() -> None:
    primary = _model("primary", [(0, 90), (180, 270)], source="P")
    candidate = _model("candidate", [(0, 90), (180, 240), (300, 360)], source="A")
    protein = ProteinEvidenceCollapsedAlignment(
        collapsed_id="protein_exact", chrom="chr1", strand="+",
        cds_parts=tuple(candidate.cds_parts), evidence_superlocus_id="psl1",
        source_member_counts={"protein_source": 1}, member_count=1,
    )
    row = _analyse(
        [primary], [candidate], protein_collapsed_alignments=[protein]
    )["candidate_rows"][0]
    assert row["protein_exact_cds_chain_support"] is True
    assert row["protein_exact_supporting_sources"] == "protein_source"
    assert "exact_protein_support" in row["fixed_support_routes"]


def _transcript(
    name: str,
    parts: tuple[tuple[int, int], ...],
    *,
    strand: str = "+",
    source: str = "tx_source",
) -> TranscriptEvidenceCollapsedTranscript:
    return TranscriptEvidenceCollapsedTranscript(
        transcript_evidence_id=name, chrom="chr1", strand=strand, exon_parts=parts,
        evidence_superlocus_id="tsl1", collapse_key_type="exact_intron_chain",
        source_member_counts={source: 1}, member_count=1, representative_raw_id=f"raw_{name}",
        span_start0=parts[0][0], span_end0=parts[-1][1],
        left_terminal_outer_min0=parts[0][0], left_terminal_outer_max0=parts[0][0],
        right_terminal_outer_min0=parts[-1][1], right_terminal_outer_max0=parts[-1][1],
    )


def _transcript_index(items: list[TranscriptEvidenceCollapsedTranscript]) -> IntervalOverlapIndex[str, TranscriptEvidenceSuperlocus]:
    locus = TranscriptEvidenceSuperlocus(
        evidence_superlocus_id="tsl1", chrom="chr1", strand="+",
        start0=min(item.start0 for item in items), end0=max(item.end0 for item in items),
        collapsed_transcripts=items,
    )
    return IntervalOverlapIndex([locus], key=lambda sl: sl.chrom, start=lambda sl: sl.start0, end=lambda sl: sl.end0)


def test_coherent_transcript_support_conjunction() -> None:
    candidate = _model("candidate", [(0, 90), (180, 270)])
    exact = _transcript("exact", ((0, 90), (180, 270)))
    associations = associated_transcripts_for_model(candidate, _transcript_index([exact]))
    result = evaluate_coherent_exact_transcript_support(
        candidate, associations, min_cdna_context_fraction=0.30,
        require_stranded_support=True, block_on_cdna_context_concern=True,
        cdna_context_neutral_fraction=0.30, cdna_context_severe_fraction=0.15,
    )
    assert result["exact_coherent_transcript_support"] is True
    assert result["qualifying_exact_transcript_group_ids"] == "exact"


def test_transcript_summary_facts_from_different_records_cannot_be_combined() -> None:
    candidate = _model("candidate", [(0, 90), (180, 270)])
    exact_unstranded = _transcript("exact_unstranded", ((0, 90), (180, 270)), strand=".")
    stranded_partial = _transcript("stranded_partial", ((0, 90), (180, 240), (330, 390)), strand="+")
    associations = associated_transcripts_for_model(
        candidate, _transcript_index([exact_unstranded, stranded_partial])
    )
    result = evaluate_coherent_exact_transcript_support(
        candidate, associations, min_cdna_context_fraction=0.30,
        require_stranded_support=True, block_on_cdna_context_concern=True,
        cdna_context_neutral_fraction=0.30, cdna_context_severe_fraction=0.15,
    )
    assert result["exact_transcript_chain_group_count"] == 1
    assert result["exact_coherent_transcript_support"] is False


def test_unavailable_genome_is_neutral_and_boundary_state_remains_neutral() -> None:
    primary = _model("primary", [(0, 90), (180, 270)], source="P")
    candidate = _model("candidate", [(0, 90), (180, 240), (300, 360)], source="A")
    unavailable = _analyse([primary], [candidate])["candidate_rows"][0]
    assert unavailable["sequence_integrity_status"] == "unavailable"
    assert unavailable["sequence_integrity_default_exclusion"] is False
    candidate.sequence_integrity = SequenceIntegrityMetrics(
        status="evaluated_no_concern", start_codon=".",
        start_codon_status="not_evaluable_at_contig_boundary", stop_codon="TAA",
        stop_codon_status="included_in_cds", internal_stop_count=0,
        phase_status="consistent", phase_mismatch_count=0,
        noncanonical_splice_intron_count=0, noncanonical_splice_intron_fraction=0.0,
        noncanonical_splice_details=".", has_ambiguous_bases=True, penalty_score=0.0,
        penalty_components=".", flags="ambiguous_sequence,five_prime_contig_boundary",
    )
    boundary = _analyse([primary], [candidate])["candidate_rows"][0]
    assert boundary["sequence_integrity_default_exclusion"] is False
    assert "five_prime_contig_boundary" in boundary["sequence_neutral_boundary_or_ambiguity_state"]


@pytest.mark.parametrize(
    ("updates", "field"),
    [
        ({"internal_stop_count": 1}, "sequence_internal_stop_concern"),
        ({"start_codon_status": "missing"}, "sequence_missing_start_concern"),
        ({"stop_codon_status": "missing"}, "sequence_missing_stop_concern"),
        ({"phase_status": "inconsistent", "phase_mismatch_count": 1}, "sequence_phase_inconsistency_concern"),
        ({"noncanonical_splice_intron_count": 1, "noncanonical_splice_intron_fraction": 0.5}, "sequence_noncanonical_splice_concern"),
    ],
)
def test_each_contract_sequence_concern_is_an_exclusion(
    updates: dict[str, object], field: str
) -> None:
    primary = _model("primary", [(0, 90), (180, 270)], source="P")
    candidate = _model("candidate", [(0, 90), (180, 240), (300, 360)], source="A")
    values: dict[str, object] = {
        "status": "evaluated_with_concern", "start_codon": "ATG",
        "start_codon_status": "present", "stop_codon": "TAA",
        "stop_codon_status": "included_in_cds", "internal_stop_count": 0,
        "phase_status": "consistent", "phase_mismatch_count": 0,
        "noncanonical_splice_intron_count": 0, "noncanonical_splice_intron_fraction": 0.0,
        "noncanonical_splice_details": ".", "has_ambiguous_bases": False,
        "penalty_score": 1.0, "penalty_components": "test", "flags": "test_concern",
    }
    values.update(updates)
    candidate.sequence_integrity = SequenceIntegrityMetrics(**values)  # type: ignore[arg-type]
    row = _analyse([primary], [candidate])["candidate_rows"][0]
    assert row[field] is True
    assert row["sequence_integrity_default_exclusion"] is True


@pytest.mark.parametrize(
    "field",
    [
        "intrinsic_penalty_cds_length_not_multiple_of_3_count",
        "intrinsic_penalty_short_cds_count",
        "intrinsic_penalty_short_cds_intron_count",
        "intrinsic_penalty_short_internal_cds_exon_count",
        "intrinsic_penalty_short_terminal_cds_exon_count",
        "intrinsic_penalty_long_cds_intron_count",
    ],
)
def test_each_contract_intrinsic_concern_is_an_exclusion(field: str) -> None:
    primary = _model("primary", [(0, 90), (180, 270)], source="P")
    candidate = _model("candidate", [(0, 90), (180, 240), (300, 360)], source="A")
    rows = {model.output_transcript_id: _score_row(model) for model in [primary, candidate]}
    rows["candidate"][field] = 1
    row = _analyse([primary], [candidate], rows=rows)["candidate_rows"][0]
    assert row["intrinsic_default_exclusion"] is True


def test_deterministic_per_gene_rank_across_original_subloci() -> None:
    primary = _model("primary", [(50, 150)], source="P", sublocus="single")
    lower = _model("lower", [(0, 50), (100, 150)], source="A", sublocus="subA")
    higher = _model("higher", [(0, 50), (100, 150)], source="B", sublocus="subB")
    rows = {
        "primary": _score_row(primary, 2.0),
        "lower": _score_row(lower, 0.5, model_selection_order_rank=1),
        "higher": _score_row(higher, 0.8, model_selection_order_rank=1),
    }
    result = _analyse([primary], [lower, higher], rows=rows)
    by_id = {row["candidate_output_transcript_id"]: row for row in result["candidate_rows"]}
    assert by_id["higher"]["alternative_candidate_order_rank"] == 1
    assert by_id["lower"]["alternative_candidate_order_rank"] == 2
    assert by_id["higher"]["historical_model_selection_order_rank"] == 1
    assert by_id["lower"]["historical_model_selection_order_rank"] == 1


def test_diagnostic_schema_contracts_have_no_duplicates() -> None:
    assert len(ALTERNATIVE_SPLICE_CANDIDATE_COLUMNS) == len(set(ALTERNATIVE_SPLICE_CANDIDATE_COLUMNS))
    assert len(ALTERNATIVE_SPLICE_GENE_COLUMNS) == len(set(ALTERNATIVE_SPLICE_GENE_COLUMNS))
    assert len(ALTERNATIVE_SPLICE_RUN_COLUMNS) == len(set(ALTERNATIVE_SPLICE_RUN_COLUMNS))
    forbidden_public_terms = (
        "stage2", "stage_2", "contract", "preview", "pending_calibration",
        "gate1", "gate2", "gate3", "gate4",
    )
    for column in (
        *ALTERNATIVE_SPLICE_CANDIDATE_COLUMNS,
        *ALTERNATIVE_SPLICE_GENE_COLUMNS,
        *ALTERNATIVE_SPLICE_RUN_COLUMNS,
    ):
        assert not any(term in column.lower() for term in forbidden_public_terms), column
    primary = _model("primary", [(0, 90), (180, 270)])
    candidate = _model("candidate", [(0, 90), (180, 240), (300, 360)], source="B")
    result = _analyse([primary], [candidate])
    assert not (set(ALTERNATIVE_SPLICE_CANDIDATE_COLUMNS) - set(result["candidate_rows"][0]))
    assert not (set(ALTERNATIVE_SPLICE_GENE_COLUMNS) - set(result["gene_rows"][0]))
    assert not (set(ALTERNATIVE_SPLICE_RUN_COLUMNS) - set(result["run_rows"][0]))


@pytest.mark.parametrize("pruning_level", ("none", "auto"))
def test_cli_report_only_outputs_do_not_change_existing_scientific_files(
    tmp_path: Path,
    pruning_level: str,
) -> None:
    def write_gff(path: Path, source: str, transcript: str) -> None:
        path.write_text(
            "##gff-version 3\n"
            f"chr1\t{source}\tgene\t1\t270\t.\t+\t.\tID=g_{transcript}\n"
            f"chr1\t{source}\tmRNA\t1\t270\t.\t+\t.\tID={transcript};Parent=g_{transcript}\n"
            f"chr1\t{source}\texon\t1\t90\t.\t+\t.\tParent={transcript}\n"
            f"chr1\t{source}\texon\t181\t270\t.\t+\t.\tParent={transcript}\n"
            f"chr1\t{source}\tCDS\t1\t90\t.\t+\t0\tParent={transcript}\n"
            f"chr1\t{source}\tCDS\t181\t270\t.\t+\t0\tParent={transcript}\n"
        )

    first = tmp_path / "first.gff3"
    second = tmp_path / "second.gff3"
    write_gff(first, "A", "txA")
    write_gff(second, "B", "txB")
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")

    def run(prefix: Path, report: bool) -> None:
        command = [
            sys.executable, "-m", "moselecta", "--inputs", str(first), str(second),
            "--source-name", "A", "B", "--out-prefix", str(prefix),
            "--pruning-level", pruning_level, "--output-profile", "minimal",
            "--progress-log-level", "none",
        ]
        if report:
            command.append("--alternative-splice-report")
        subprocess.run(command, cwd=root, env=env, check=True, text=True, capture_output=True)

    baseline = tmp_path / "baseline"
    report = tmp_path / "report"
    run(baseline, False)
    run(report, True)
    baseline_files = {path.name.removeprefix("baseline."): path for path in tmp_path.glob("baseline.*")}
    report_files = {path.name.removeprefix("report."): path for path in tmp_path.glob("report.*")}
    for suffix, path in baseline_files.items():
        assert suffix in report_files
        if suffix == "run_configuration.json":
            continue
        assert path.read_bytes() == report_files[suffix].read_bytes(), suffix
    assert {name for name in report_files if name.startswith("alternative_splice_")} == {
        "alternative_splice_candidates.tsv",
        "alternative_splice_genes.tsv",
        "alternative_splice_run_summary.tsv",
    }
    with report_files["alternative_splice_candidates.tsv"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["alternative_splice_decision"] == "rejected_exact_equivalent"


def test_cli_report_requires_final_stop_point(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    completed = subprocess.run(
        [
            sys.executable, "-m", "moselecta",
            "--inputs", str(tmp_path / "not-read.gff3"),
            "--out-prefix", str(tmp_path / "out"),
            "--stop-after", "backfill",
            "--alternative-splice-report",
        ],
        cwd=root,
        env=env,
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode != 0
    assert "Alternative-splice diagnostics require --stop-after final" in completed.stderr


def test_cli_default_final_report_and_audit_outputs_share_one_selector_result(
    tmp_path: Path,
) -> None:
    first = tmp_path / "first.gff3"
    second = tmp_path / "second.gff3"
    first.write_text(
        "##gff-version 3\n"
        "chr1\tA\tgene\t1\t450\t.\t+\t.\tID=gA\n"
        "chr1\tA\tmRNA\t1\t450\t.\t+\t.\tID=txA;Parent=gA\n"
        "chr1\tA\texon\t1\t90\t.\t+\t.\tParent=txA\n"
        "chr1\tA\texon\t181\t270\t.\t+\t.\tParent=txA\n"
        "chr1\tA\texon\t361\t450\t.\t+\t.\tParent=txA\n"
        "chr1\tA\tCDS\t1\t90\t.\t+\t0\tParent=txA\n"
        "chr1\tA\tCDS\t181\t270\t.\t+\t0\tParent=txA\n"
        "chr1\tA\tCDS\t361\t450\t.\t+\t0\tParent=txA\n"
    )
    second.write_text(
        "##gff-version 3\n"
        "chr1\tB\tgene\t1\t450\t.\t+\t.\tID=gB\n"
        "chr1\tB\tmRNA\t1\t450\t.\t+\t.\tID=txB;Parent=gB\n"
        "chr1\tB\texon\t1\t90\t.\t+\t.\tParent=txB\n"
        "chr1\tB\texon\t181\t240\t.\t+\t.\tParent=txB\n"
        "chr1\tB\texon\t331\t450\t.\t+\t.\tParent=txB\n"
        "chr1\tB\tCDS\t1\t90\t.\t+\t0\tParent=txB\n"
        "chr1\tB\tCDS\t181\t240\t.\t+\t0\tParent=txB\n"
        "chr1\tB\tCDS\t331\t450\t.\t+\t0\tParent=txB\n"
    )
    first_replica = tmp_path / "first_replica.gff3"
    second_replica = tmp_path / "second_replica.gff3"
    first_replica.write_text(first.read_text())
    second_replica.write_text(second.read_text())
    unrelated_protein = tmp_path / "unrelated_protein.gff3"
    unrelated_protein.write_text(
        "##gff-version 3\n"
        "chr9\tprotein\tCDS\t1\t90\t.\t+\t0\tID=p1.cds1;Parent=p1\n"
    )
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")

    def run(
        prefix: Path, extra: list[str] | None = None
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable, "-m", "moselecta",
                "--inputs", str(first), str(first_replica), str(second), str(second_replica),
                "--source-name", "A", "A2", "B", "B2",
                "--out-prefix", str(prefix),
                "--pruning-level", "none",
                "--output-profile", "standard",
                "--progress-log-level", "none",
                *(extra or []),
            ],
            cwd=root, env=env, check=True, text=True, capture_output=True,
        )

    default_prefix = tmp_path / "default"
    report_prefix = tmp_path / "report"
    audit_prefix = tmp_path / "audit"
    primary_prefix = tmp_path / "primary"
    primary_report_prefix = tmp_path / "primary_report"
    evidence_required_prefix = tmp_path / "evidence_required"
    progress_prefix = tmp_path / "progress"
    default_completed = run(default_prefix)
    report_completed = run(report_prefix, ["--alternative-splice-report"])
    run(audit_prefix, [
        "--alternative-splice-report",
        "--write-output", "alternative_splice_candidate_audit_gff3",
    ])
    primary_completed = run(primary_prefix, ["--no-alternative-splice"])
    primary_report_completed = run(primary_report_prefix, [
        "--no-alternative-splice", "--alternative-splice-report",
    ])
    run(evidence_required_prefix, [
        "--evidence-proteins", str(unrelated_protein),
        "--alternative-splice-evidence-requirement", "protein",
        "--alternative-splice-report",
    ])
    progress_completed = run(
        progress_prefix, ["--progress-log-level", "basic"]
    )
    assert "Alternative-splice selection\n============================" in default_completed.stdout
    assert "Threshold preset: intermediate" in default_completed.stdout
    assert "Support policy: exact" in default_completed.stdout
    assert "External evidence requirement: none" in default_completed.stdout
    assert "Alternatives emitted into final annotation: true" in default_completed.stdout
    assert "Accepted alternative transcripts (selector result): 1" in default_completed.stdout
    assert "Genes with alternative transcripts (selector result): 1" in default_completed.stdout
    assert "Final annotation: genes=1; transcripts=2" in default_completed.stdout
    assert "Final annotation genes: 1" in default_completed.stdout
    assert "Final annotation transcripts: 2" in default_completed.stdout
    assert "  primary transcripts: 1" in default_completed.stdout
    assert "  alternative transcripts: 1" in default_completed.stdout
    assert "Alternative-splice candidates:" not in default_completed.stdout
    assert "standardised primary annotation GFF3:" in default_completed.stdout
    assert "final annotation GFF3:" in default_completed.stdout
    assert "canonical final public annotation" not in default_completed.stdout.lower()
    assert "Final annotation\t" in default_completed.stdout
    assert "Final public annotation" not in default_completed.stdout

    assert "Alternative-splice augmentation: disabled" in primary_completed.stdout
    assert "Threshold preset:" not in primary_completed.stdout
    assert "Final annotation transcripts: 1" in primary_completed.stdout
    assert "  alternative transcripts: 0" in primary_completed.stdout

    assert "Threshold preset: intermediate" in primary_report_completed.stdout
    assert (
        "Alternatives emitted into final annotation: false"
        in primary_report_completed.stdout
    )
    assert (
        "Accepted alternative transcripts (selector result): 1"
        in primary_report_completed.stdout
    )
    assert "Final annotation: genes=1; transcripts=1" in primary_report_completed.stdout
    assert "  alternative transcripts: 0" in primary_report_completed.stdout
    assert report_completed.returncode == 0

    assert "Writing standardised primary annotation" in progress_completed.stderr
    assert "Finished standardised primary annotation" in progress_completed.stderr
    assert "Starting alternative-splice selection" in progress_completed.stderr
    assert "Finished alternative-splice selection: accepted=1" in progress_completed.stderr
    assert "Writing final annotation" in progress_completed.stderr
    assert "Finished final annotation: genes=1; transcripts=2" in progress_completed.stderr
    assert "retained primaries were frozen" not in progress_completed.stderr
    assert "canonical final public annotation" not in progress_completed.stderr.lower()
    report_files = {
        path.name.removeprefix("report."): path for path in tmp_path.glob("report.*")
    }
    audit_files = {
        path.name.removeprefix("audit."): path for path in tmp_path.glob("audit.*")
    }
    for suffix, report_path in report_files.items():
        if suffix in {
            "alternative_splice_candidate_audit.gff3",
            "run_configuration.json",
        }:
            continue
        assert suffix in audit_files
        assert report_path.read_bytes() == audit_files[suffix].read_bytes(), suffix
    for suffix in (
        ".alternative_splice_candidates.tsv",
        ".alternative_splice_genes.tsv",
        ".alternative_splice_run_summary.tsv",
    ):
        assert Path(str(report_prefix) + suffix).read_bytes() == Path(
            str(audit_prefix) + suffix
        ).read_bytes()
    standardised_primary_bytes = Path(
        str(report_prefix) + ".standardised_primary_models.gff3"
    ).read_bytes()
    for prefix in (
        default_prefix, audit_prefix, primary_prefix, primary_report_prefix,
    ):
        assert Path(
            str(prefix) + ".standardised_primary_models.gff3"
        ).read_bytes() == standardised_primary_bytes
    assert Path(str(default_prefix) + ".final_annotation.gff3").exists()
    default_final_lines = Path(
        str(default_prefix) + ".final_annotation.gff3"
    ).read_text().splitlines()
    primary_final_lines = Path(
        str(primary_prefix) + ".final_annotation.gff3"
    ).read_text().splitlines()
    assert any("ID=MoSelecta_000001.2" in line for line in default_final_lines)
    assert not any("ID=MoSelecta_000001.2" in line for line in primary_final_lines)
    assert [
        line for line in default_final_lines
        if (
            "MoSelecta_000001.1" in line
            and not line.startswith("#")
            and line.split("\t")[2] != "gene"
        )
    ] == [
        line for line in primary_final_lines
        if (
            "MoSelecta_000001.1" in line
            and not line.startswith("#")
            and line.split("\t")[2] != "gene"
        )
    ]
    assert "# moselecta_alternative_splice_enabled=true" in default_final_lines
    assert "# moselecta_alternative_splice_selector_evaluated=true" in default_final_lines
    assert "# moselecta_alternative_splice_alternatives_emitted=true" in default_final_lines
    assert "# moselecta_alternative_splice_enabled=false" in primary_final_lines
    assert "# moselecta_alternative_splice_selector_evaluated=false" in primary_final_lines
    assert "# moselecta_alternative_splice_alternatives_emitted=false" in primary_final_lines
    primary_report_final = Path(
        str(primary_report_prefix) + ".final_annotation.gff3"
    ).read_text().splitlines()
    assert "# moselecta_alternative_splice_enabled=false" in primary_report_final
    assert "# moselecta_alternative_splice_selector_evaluated=true" in primary_report_final
    assert "# moselecta_alternative_splice_alternatives_emitted=false" in primary_report_final
    evidence_required_final = Path(
        str(evidence_required_prefix) + ".final_annotation.gff3"
    ).read_text()
    assert "ID=MoSelecta_000001.2" not in evidence_required_final
    with Path(
        str(evidence_required_prefix) + ".alternative_splice_candidates.tsv"
    ).open() as handle:
        evidence_required_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert any(
        row["alternative_splice_decision"]
        == "rejected_external_evidence_requirement"
        for row in evidence_required_rows
    )
    assert not Path(str(default_prefix) + ".alternative_splice_candidates.tsv").exists()
    assert not Path(str(report_prefix) + ".alternative_splice_candidate_audit.gff3").exists()
    assert Path(str(audit_prefix) + ".alternative_splice_candidate_audit.gff3").exists()
    public_diagnostic_paths = [
        Path(str(audit_prefix) + suffix)
        for suffix in (
            ".alternative_splice_candidates.tsv",
            ".alternative_splice_genes.tsv",
            ".alternative_splice_run_summary.tsv",
            ".alternative_splice_candidate_audit.gff3",
            ".final_annotation.gff3",
        )
    ]
    for path in public_diagnostic_paths:
        text = path.read_text().lower()
        for obsolete_term in (
            "stage2", "stage_2", "contract-v1", "preview",
            "pending_calibration", "gate3", "unset_pending_calibration",
        ):
            assert obsolete_term not in text, (path, obsolete_term)

    for name, options in (
        ("permissive", ["--alternative-splice-preset", "permissive"]),
        ("stringent", ["--alternative-splice-preset", "stringent"]),
        (
            "mosaic",
            ["--alternative-splice-support-policy", "exact-plus-mosaic"],
        ),
    ):
        configured_prefix = tmp_path / name
        run(configured_prefix, options)
        assert Path(
            str(configured_prefix) + ".standardised_primary_models.gff3"
        ).read_bytes() == Path(
            str(report_prefix) + ".standardised_primary_models.gff3"
        ).read_bytes()

def test_full_output_profile_writes_reports_but_not_large_audit(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    prefix = tmp_path / "full"
    subprocess.run(
        [
            sys.executable, "-m", "moselecta",
            "--inputs",
            str(root / "tests" / "data" / "alternative_splice_smoke_a.gff3"),
            str(root / "tests" / "data" / "alternative_splice_smoke_b.gff3"),
            "--source-name", "A", "B",
            "--out-prefix", str(prefix),
            "--pruning-level", "none",
            "--output-profile", "full",
            "--progress-log-level", "none",
        ],
        cwd=root, env=env, check=True, text=True, capture_output=True,
    )
    assert Path(str(prefix) + ".alternative_splice_candidates.tsv").exists()
    assert Path(str(prefix) + ".alternative_splice_genes.tsv").exists()
    assert Path(str(prefix) + ".alternative_splice_run_summary.tsv").exists()
    assert not Path(str(prefix) + ".alternative_splice_candidate_audit.gff3").exists()


def test_earlier_stop_full_profile_does_not_activate_default_selector(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    prefix = tmp_path / "backfill"
    completed = subprocess.run(
        [
            sys.executable, "-m", "moselecta",
            "--inputs",
            str(root / "tests/data/alternative_splice_smoke_a.gff3"),
            str(root / "tests/data/alternative_splice_smoke_b.gff3"),
            "--source-name", "A", "B",
            "--out-prefix", str(prefix),
            "--stop-after", "backfill",
            "--pruning-level", "none",
            "--output-profile", "full",
            "--progress-log-level", "none",
        ],
        cwd=root,
        env=env,
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0
    assert not Path(str(prefix) + ".alternative_splice_candidates.tsv").exists()
    assert not Path(str(prefix) + ".final_annotation.gff3").exists()


@pytest.mark.parametrize(
    "requirement",
    ("protein-or-transcript", "protein", "transcript", "protein-and-transcript"),
)
def test_cli_rejects_unavailable_external_evidence_requirement(
    tmp_path: Path,
    requirement: str,
) -> None:
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    completed = subprocess.run(
        [
            sys.executable, "-m", "moselecta",
            "--inputs", str(root / "tests/data/alternative_splice_smoke_a.gff3"),
            "--out-prefix", str(tmp_path / requirement),
            "--alternative-splice-evidence-requirement", requirement,
        ],
        cwd=root,
        env=env,
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode != 0
    assert "requires compatible configured evidence inputs" in completed.stderr


def test_primary_only_unused_selector_options_warn(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    completed = subprocess.run(
        [
            sys.executable, "-m", "moselecta",
            "--inputs", str(root / "tests/data/alternative_splice_smoke_a.gff3"),
            "--out-prefix", str(tmp_path / "primary"),
            "--pruning-level", "none",
            "--output-profile", "minimal",
            "--progress-log-level", "none",
            "--no-alternative-splice",
            "--alternative-splice-preset", "intermediate",
            "--alternative-splice-evidence-requirement", "protein",
        ],
        cwd=root,
        env=env,
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0
    assert "selector options have no effect" in completed.stderr
    assert "--alternative-splice-evidence-requirement" in completed.stderr
