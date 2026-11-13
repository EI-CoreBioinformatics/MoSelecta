"""Regression coverage for the O5a internal-stop and incomplete-end policy."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from moselecta.final_annotation import FinalAnnotation, FinalGene, FinalTranscript
from moselecta.models import CodingModel, SequenceIntegrityMetrics, Sublocus
from moselecta.selection import (
    count_retained_internal_stop_evidence_overrides,
    final_prune_sequence_integrity_nonsingletons,
    final_prune_singleton_selected_models,
    write_final_pruned_outputs,
)
from moselecta.utr_augmentation import (
    UTREvidence,
    UTRSettings,
    augment_final_annotation,
)


def _model(model_id: str, source: str = "SourceA", start0: int = 100) -> CodingModel:
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"gene_{model_id}",
        input_transcript_id=model_id,
        output_gene_id=f"out_gene_{model_id}",
        output_transcript_id=model_id,
        chrom="chr1",
        strand="+",
        cds_parts=[(start0, start0 + 51), (start0 + 150, start0 + 301)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def _prune_result(models: list[CodingModel], rows: list[dict[str, object]]) -> dict[str, object]:
    return {
        "initial_retained_models": list(models),
        "initial_retained_rows": list(rows),
        "retained_models": list(models),
        "retained_rows": list(rows),
        "removed_rows": [],
        "retained_selection_info_by_transcript": {
            str(row["output_transcript_id"]): dict(row) for row in rows
        },
    }


def _protein_row(score: float, exact_chain_sources: int = 0) -> dict[str, object]:
    return {
        "protein_evidence_score_scaled": score,
        "protein_exact_cds_chain_source_count": exact_chain_sources,
        "protein_supporting_source_count": max(exact_chain_sources, int(score > 0)),
    }


def _transcript_row(score: float = 0.90) -> dict[str, object]:
    return {
        "transcript_evidence_score_scaled": score,
        "transcript_supporting_source_count": 2,
        "transcript_scoring_eligible_count": 2,
        "transcript_stranded_scoring_eligible_count": 2,
        "transcript_unstranded_scoring_eligible_count": 0,
        "transcript_cdna_context_evaluable": "true",
        "transcript_cdna_context_cds_to_min_cdna_fraction": 0.50,
        "transcript_cdna_context_pruning_concern": "false",
    }


def _run_singleton(
    flags: str,
    *,
    protein: dict[str, object] | None = None,
    transcript: dict[str, object] | None = None,
    allow_internal_stop_evidence_rescue: bool = False,
) -> dict[str, object]:
    model = _model("candidate")
    sub = Sublocus("sub1", "sl1", "chr1", "+", model.cds_start0, model.cds_end0, [model])
    row = {
        "output_transcript_id": model.output_transcript_id,
        "source_name": model.source_name,
        "final_prune_classes": ".",
    }
    structural = {
        model.output_transcript_id: {
            "intrinsic_penalty_components": ".",
            "intrinsic_penalty_score": 0.0,
            "sequence_integrity_flags": flags,
        }
    }
    return final_prune_singleton_selected_models(
        _prune_result([model], [row]),
        [sub],
        structural,
        prune_level="sequence_integrity",
        prune_sources=set(),
        pruning_level_name="high",
        rescue_protein_supported_singletons=True,
        protein_evidence_rows_by_transcript=(
            {model.output_transcript_id: protein} if protein is not None else {}
        ),
        rescue_transcript_supported_singletons=True,
        transcript_evidence_rows_by_transcript=(
            {model.output_transcript_id: transcript} if transcript is not None else {}
        ),
        allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
    )


@pytest.mark.parametrize(
    ("protein", "transcript"),
    [
        pytest.param(None, None, id="no-evidence"),
        pytest.param(_protein_row(0.749), None, id="protein-below-threshold"),
        pytest.param(_protein_row(0.750), None, id="protein-at-threshold"),
        pytest.param(_protein_row(0.10, 1), None, id="protein-exact-chain-only"),
        pytest.param(None, _transcript_row(), id="strong-transcript"),
    ],
)
def test_singleton_internal_stop_is_not_evidence_rescued_by_default(
    protein: dict[str, object] | None,
    transcript: dict[str, object] | None,
) -> None:
    result = _run_singleton(
        "internal_stop_codon", protein=protein, transcript=transcript
    )
    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_sequence_integrity_singleton_candidate"] == "true"
    assert row["final_sequence_integrity_singleton_rescued"] == "false"
    assert row["final_singleton_rescue_reason"] == "."
    assert row["final_allow_internal_stop_evidence_rescue"] == "false"
    criteria_met = protein in (_protein_row(0.750), _protein_row(0.10, 1)) or transcript is not None
    assert row["final_internal_stop_evidence_criteria_met"] == str(criteria_met).lower()
    assert row["final_internal_stop_evidence_rescue_blocked"] == str(criteria_met).lower()


@pytest.mark.parametrize(
    ("protein", "transcript", "route"),
    [
        pytest.param(_protein_row(0.750), None, "protein_evidence", id="protein-at-threshold"),
        pytest.param(_protein_row(0.10, 1), None, "protein_evidence", id="protein-exact-chain-only"),
        pytest.param(None, _transcript_row(), "transcript_evidence", id="strong-transcript"),
    ],
)
def test_singleton_internal_stop_override_restores_evidence_rescue(
    protein: dict[str, object] | None,
    transcript: dict[str, object] | None,
    route: str,
) -> None:
    result = _run_singleton(
        "internal_stop_codon",
        protein=protein,
        transcript=transcript,
        allow_internal_stop_evidence_rescue=True,
    )
    assert len(result["retained_models"]) == 1
    row = result["retained_rows"][0]
    assert row["final_singleton_rescue_reason"] == route
    assert row["final_allow_internal_stop_evidence_rescue"] == "true"
    assert row["final_internal_stop_evidence_rescue_applied"] == "true"
    assert row["final_internal_stop_evidence_rescue_applied_route"] == route


def _run_internal_stop_cross_stage_provenance(
    *, allow_internal_stop_evidence_rescue: bool
) -> dict[str, object]:
    protein = _protein_row(0.80)
    singleton_result = _run_singleton(
        "internal_stop_codon",
        protein=protein,
        allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
    )
    model = singleton_result["initial_retained_models"][0]
    sub = Sublocus(
        "sub1", "sl1", "chr1", "+", model.cds_start0, model.cds_end0, [model]
    )
    return final_prune_sequence_integrity_nonsingletons(
        singleton_result,
        [sub],
        {model.output_transcript_id: {
            "sequence_integrity_flags": "internal_stop_codon",
        }},
        enabled=True,
        scope="all",
        pruning_level_name="high",
        protect_protein_supported=True,
        protein_evidence_rows_by_transcript={model.output_transcript_id: protein},
        allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
    )


def test_later_nonsingleton_pass_preserves_singleton_override_provenance_and_tsv(
    tmp_path: Path,
) -> None:
    result = _run_internal_stop_cross_stage_provenance(
        allow_internal_stop_evidence_rescue=True
    )
    assert len(result["retained_models"]) == 1
    row = result["retained_rows"][0]
    assert row["final_internal_stop_evidence_criteria_met"] == "true"
    assert row["final_internal_stop_evidence_criteria_routes"] == "protein_evidence"
    assert row["final_internal_stop_evidence_rescue_applied"] == "true"
    assert row["final_internal_stop_evidence_rescue_applied_route"] == "protein_evidence"
    assert count_retained_internal_stop_evidence_overrides(result["retained_rows"]) == 1
    assert result["internal_stop_evidence_rescue_applied_count"] == 1

    retained_path = tmp_path / "retained.tsv"
    write_final_pruned_outputs(
        result,
        final_pruned_selected_models_tsv=retained_path,
        final_pruned_selected_removed_tsv=tmp_path / "removed.tsv",
    )
    with retained_path.open(newline="") as handle:
        written_row = next(csv.DictReader(handle, delimiter="\t"))
    assert written_row["final_internal_stop_evidence_criteria_met"] == "true"
    assert written_row["final_internal_stop_evidence_rescue_applied"] == "true"
    assert written_row["final_internal_stop_evidence_rescue_applied_route"] == "protein_evidence"


def test_later_nonsingleton_pass_preserves_default_blocked_provenance() -> None:
    result = _run_internal_stop_cross_stage_provenance(
        allow_internal_stop_evidence_rescue=False
    )
    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_internal_stop_evidence_criteria_met"] == "true"
    assert row["final_internal_stop_evidence_criteria_routes"] == "protein_evidence"
    assert row["final_internal_stop_evidence_rescue_blocked"] == "true"
    assert (
        row["final_internal_stop_evidence_rescue_block_reason"]
        == "internal_stop_codon_override_disabled"
    )
    assert row["final_internal_stop_evidence_rescue_applied"] == "false"


def _run_nonsingleton(
    flags: str,
    *,
    protein: dict[str, object] | None = None,
    transcript: dict[str, object] | None = None,
    backfill: bool = False,
    allow_internal_stop_evidence_rescue: bool = False,
) -> dict[str, object]:
    candidate = _model("candidate", start0=100)
    companion = _model("companion", source="SourceB", start0=600)
    sub = Sublocus("sub1", "sl1", "chr1", "+", 100, 901, [candidate, companion])
    candidate_row = {
        "output_transcript_id": candidate.output_transcript_id,
        "source_name": candidate.source_name,
        "final_selection_stage": "backfill" if backfill else "monosublocus",
        "backfill_selected": "true" if backfill else "false",
        "sequence_integrity_flags": flags,
    }
    companion_row = {
        "output_transcript_id": companion.output_transcript_id,
        "source_name": companion.source_name,
        "final_selection_stage": "monosublocus",
        "sequence_integrity_flags": ".",
    }
    structural = {
        candidate.output_transcript_id: {"sequence_integrity_flags": flags},
        companion.output_transcript_id: {"sequence_integrity_flags": "."},
    }
    return final_prune_sequence_integrity_nonsingletons(
        _prune_result([candidate, companion], [candidate_row, companion_row]),
        [sub],
        structural,
        enabled=True,
        scope="all",
        pruning_level_name="high",
        protect_protein_supported=True,
        protein_evidence_rows_by_transcript=(
            {candidate.output_transcript_id: protein} if protein is not None else {}
        ),
        protect_transcript_supported=True,
        transcript_evidence_rows_by_transcript=(
            {candidate.output_transcript_id: transcript} if transcript is not None else {}
        ),
        allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
    )


@pytest.mark.parametrize(
    ("protein", "transcript"),
    [
        pytest.param(None, None, id="no-evidence"),
        pytest.param(_protein_row(0.80), None, id="high-protein"),
        pytest.param(_protein_row(0.10, 1), None, id="protein-exact-chain-only"),
        pytest.param(None, _transcript_row(), id="strong-transcript"),
    ],
)
def test_nonsingleton_internal_stop_is_not_evidence_protected_by_default(
    protein: dict[str, object] | None,
    transcript: dict[str, object] | None,
) -> None:
    result = _run_nonsingleton(
        "internal_stop_codon", protein=protein, transcript=transcript
    )
    retained_ids = {
        model.output_transcript_id for model in result["retained_models"]
    }
    assert retained_ids == {"companion"}
    row = next(
        row for row in result["removed_rows"]
        if row["output_transcript_id"] == "candidate"
    )
    assert row["final_sequence_integrity_nonsingleton_pruned"] == "true"
    criteria_met = protein in (_protein_row(0.80), _protein_row(0.10, 1)) or transcript is not None
    assert row["final_internal_stop_evidence_criteria_met"] == str(criteria_met).lower()
    assert row["final_internal_stop_evidence_rescue_blocked"] == str(criteria_met).lower()


@pytest.mark.parametrize(
    ("protein", "transcript", "protected_by", "route"),
    [
        pytest.param(_protein_row(0.80), None, "protein", "protein_evidence", id="high-protein"),
        pytest.param(_protein_row(0.10, 1), None, "protein", "protein_evidence", id="protein-exact-chain-only"),
        pytest.param(None, _transcript_row(), "transcript", "transcript_evidence", id="strong-transcript"),
    ],
)
def test_nonsingleton_internal_stop_override_restores_evidence_protection(
    protein: dict[str, object] | None,
    transcript: dict[str, object] | None,
    protected_by: str,
    route: str,
) -> None:
    result = _run_nonsingleton(
        "internal_stop_codon",
        protein=protein,
        transcript=transcript,
        allow_internal_stop_evidence_rescue=True,
    )
    retained_ids = {model.output_transcript_id for model in result["retained_models"]}
    assert retained_ids == {"candidate", "companion"}
    row = next(row for row in result["retained_rows"] if row["output_transcript_id"] == "candidate")
    assert row[f"final_sequence_integrity_nonsingleton_{protected_by}_protected"] == "true"
    assert row["final_internal_stop_evidence_rescue_applied"] == "true"
    assert row["final_internal_stop_evidence_rescue_applied_route"] == route


@pytest.mark.parametrize("flags", ["missing_start_codon", "missing_stop_codon"])
def test_incomplete_singleton_remains_protein_rescuable(flags: str) -> None:
    result = _run_singleton(flags, protein=_protein_row(0.80))
    assert len(result["retained_models"]) == 1
    assert result["retained_rows"][0]["final_protein_singleton_rescue_applied"] == "true"


@pytest.mark.parametrize(
    ("protein", "transcript", "route"),
    [
        pytest.param(_protein_row(0.80), None, "protein_evidence", id="protein"),
        pytest.param(None, _transcript_row(), "transcript_evidence", id="transcript"),
    ],
)
def test_ordinary_singleton_evidence_rescue_is_unchanged(
    protein: dict[str, object] | None,
    transcript: dict[str, object] | None,
    route: str,
) -> None:
    result = _run_singleton(
        "missing_start_codon", protein=protein, transcript=transcript
    )
    assert len(result["retained_models"]) == 1
    row = result["retained_rows"][0]
    assert row["final_singleton_rescue_reason"] == route
    assert row["final_model_has_internal_stop_codon"] == "false"


@pytest.mark.parametrize("flags", ["missing_start_codon", "missing_stop_codon"])
def test_incomplete_backfill_nonsingleton_remains_protein_protectable(flags: str) -> None:
    result = _run_nonsingleton(flags, protein=_protein_row(0.80), backfill=True)
    retained_ids = {model.output_transcript_id for model in result["retained_models"]}
    assert retained_ids == {"candidate", "companion"}
    row = next(
        row for row in result["retained_rows"]
        if row["output_transcript_id"] == "candidate"
    )
    assert row["final_sequence_integrity_nonsingleton_protein_protected"] == "true"


def _utr_model(start_status: str, stop_status: str) -> CodingModel:
    flags = []
    if start_status == "missing":
        flags.append("missing_start_codon")
    if stop_status == "missing":
        flags.append("missing_stop_codon")
    model = CodingModel(
        source_name="prediction",
        input_file=Path("prediction.gff3"),
        input_gene_id="gene_utr",
        input_transcript_id="utr",
        output_gene_id="out_gene_utr",
        output_transcript_id="utr",
        chrom="chr1",
        strand="+",
        cds_parts=[(100, 150), (200, 250)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )
    model.sequence_integrity = SequenceIntegrityMetrics(
        status="evaluated_with_concern",
        start_codon="ATG" if start_status == "present" else "CCC",
        start_codon_status=start_status,
        stop_codon="TAA" if stop_status == "included_in_cds" else "CCC",
        stop_codon_status=stop_status,
        internal_stop_count=0,
        phase_status="consistent",
        phase_mismatch_count=0,
        noncanonical_splice_intron_count=0,
        noncanonical_splice_intron_fraction=0.0,
        noncanonical_splice_details=".",
        has_ambiguous_bases=False,
        penalty_score=0.30 + (0.50 if stop_status == "missing" else 0.0),
        penalty_components=",".join(flags),
        flags=",".join(flags),
    )
    return model


@pytest.mark.parametrize(
    ("start_status", "stop_status", "expected_five", "expected_three"),
    [
        pytest.param("missing", "included_in_cds", False, True, id="missing-start-blocks-five-prime"),
        pytest.param("present", "missing", True, False, id="missing-stop-blocks-three-prime"),
        pytest.param("missing", "missing", False, False, id="both-incomplete-block-both-sides"),
        pytest.param("present", "included_in_cds", True, True, id="complete-allows-both-sides"),
        pytest.param("boundary_neutral", "included_in_cds", True, True, id="boundary-neutral-start-is-not-blocked"),
        pytest.param("present", "boundary_neutral", True, True, id="boundary-neutral-stop-is-not-blocked"),
    ],
)
def test_utr_augmentation_blocks_only_explicitly_incomplete_cds_ends(
    start_status: str,
    stop_status: str,
    expected_five: bool,
    expected_three: bool,
) -> None:
    model = _utr_model(start_status, stop_status)
    annotation = FinalAnnotation((
        FinalGene("MoSelecta_000001", (
            FinalTranscript("MoSelecta_000001.1", model),
        )),
    ))
    evidence = UTREvidence(
        evidence_id="rna",
        kind="transcript",
        chrom="chr1",
        strand="+",
        exon_parts=((80, 150), (200, 270)),
        source_member_counts=(("rna", 1),),
    )
    settings = UTRSettings(
        five_prime_exonic_override_bp=100_000,
        three_prime_exonic_override_bp=100_000,
        five_prime_genomic_override_bp=100_000,
        three_prime_genomic_override_bp=100_000,
    )
    result = augment_final_annotation(annotation, [evidence], [], settings)
    augmented = result.annotation.genes[0].transcripts[0]
    assert bool(augmented.five_prime_utr_parts) is expected_five
    assert bool(augmented.three_prime_utr_parts) is expected_three
    assert augmented.model.cds_parts == [(100, 150), (200, 250)]
    audit = result.audit_rows[0]
    assert audit["five_prime_anchor_block_reason"] == (
        "missing_start_codon" if start_status == "missing" else "."
    )
    assert audit["three_prime_anchor_block_reason"] == (
        "missing_stop_codon" if stop_status == "missing" else "."
    )
    if start_status == "missing":
        assert audit["five_prime_quality_status"] == "blocked_incomplete_cds_anchor"
    if stop_status == "missing":
        assert audit["three_prime_quality_status"] == "blocked_incomplete_cds_anchor"


def test_blocked_incomplete_utr_side_does_not_contribute_to_limit_learning() -> None:
    model = _utr_model("missing", "included_in_cds")
    annotation = FinalAnnotation((
        FinalGene("MoSelecta_000001", (
            FinalTranscript("MoSelecta_000001.1", model),
        )),
    ))
    evidence = UTREvidence(
        evidence_id="five-only",
        kind="transcript",
        chrom="chr1",
        strand="+",
        exon_parts=((80, 150), (200, 250)),
        source_member_counts=(("rna", 1),),
    )
    settings = UTRSettings(limit_min_observations=1)
    result = augment_final_annotation(annotation, [evidence], [], settings)
    assert result.annotation.genes[0].transcripts[0].five_prime_utr_parts == ()
    assert result.limits.five_prime_exonic_source == "sparse_data_exonic_failsafe"
    assert result.limits.five_prime_genomic_source == "sparse_data_unavailable"
    assert (
        "incomplete_cds_anchor_missing_start_codon", 1
    ) in result.summary.rejection_counts
