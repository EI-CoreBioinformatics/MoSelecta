from __future__ import annotations

import csv
from pathlib import Path

import pytest

from moselecta.models import CodingModel
from moselecta.protein_evidence import (
    ProteinEvidenceCollapsedAlignment,
    ProteinEvidenceSuperlocus,
    calculate_protein_model_support,
    pairwise_model_protein_metrics,
    write_protein_model_support_tsv,
)
from moselecta.transcript_evidence import (
    TranscriptEvidenceCollapsedTranscript,
    TranscriptEvidenceSuperlocus,
    calculate_transcript_model_support,
    transcript_model_pair_metrics,
    write_transcript_model_support_tsv,
)


def cds_parts(exon_count: int, *, offset: int = 0) -> tuple[tuple[int, int], ...]:
    return tuple((offset + 100 + 150 * index, offset + 150 + 150 * index) for index in range(exon_count))


def make_model(parts: tuple[tuple[int, int], ...], model_id: str) -> CodingModel:
    return CodingModel(
        source_name="model",
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


def make_protein_alignment(
    parts: tuple[tuple[int, int], ...],
    alignment_id: str,
) -> ProteinEvidenceCollapsedAlignment:
    return ProteinEvidenceCollapsedAlignment(
        collapsed_id=alignment_id,
        chrom="chr1",
        strand="+",
        cds_parts=parts,
        evidence_superlocus_id="protein_superlocus",
        source_member_counts={"protein": 1},
        member_count=1,
    )


def make_transcript(
    parts: tuple[tuple[int, int], ...],
    transcript_id: str,
    *,
    member_count: int = 1,
) -> TranscriptEvidenceCollapsedTranscript:
    return TranscriptEvidenceCollapsedTranscript(
        transcript_evidence_id=transcript_id,
        chrom="chr1",
        strand="+",
        exon_parts=parts,
        evidence_superlocus_id="transcript_superlocus",
        collapse_key_type="exact_intron_chain",
        source_member_counts={"transcript": member_count},
        member_count=member_count,
        span_start0=parts[0][0],
        span_end0=parts[-1][1],
        left_terminal_outer_min0=parts[0][0],
        left_terminal_outer_max0=parts[0][0],
        right_terminal_outer_min0=parts[-1][1],
        right_terminal_outer_max0=parts[-1][1],
    )


def read_tsv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["output_transcript_id"]: row for row in rows}


@pytest.mark.parametrize("exon_count", [1, 2, 3, 4])
def test_perfect_protein_pairwise_score_across_cds_topologies(exon_count: int) -> None:
    parts = cds_parts(exon_count)
    metrics = pairwise_model_protein_metrics(
        make_model(parts, f"model_{exon_count}"),
        make_protein_alignment(parts, f"protein_{exon_count}"),
        wobble_bp=3,
    )

    assert metrics["score"] == pytest.approx(1.0)
    if exon_count <= 2:
        assert metrics["internal_cds_exon_f1"] is None
    else:
        assert metrics["internal_cds_exon_f1"] == pytest.approx(1.0)


def test_protein_internal_exon_na_is_renormalized_and_serialized(tmp_path: Path) -> None:
    two_parts = cds_parts(2)
    three_parts = cds_parts(3, offset=2000)
    unsupported_three_parts = (three_parts[0], (three_parts[1][0], three_parts[1][1] - 10), three_parts[2])
    two_model = make_model(two_parts, "two_exon")
    three_model = make_model(three_parts, "three_exon_unsupported")
    two_alignment = make_protein_alignment(two_parts, "protein_two")
    three_alignment = make_protein_alignment(unsupported_three_parts, "protein_three")
    superlocus = ProteinEvidenceSuperlocus(
        evidence_superlocus_id="protein_superlocus",
        chrom="chr1",
        strand="+",
        start0=two_parts[0][0],
        end0=unsupported_three_parts[-1][1],
        collapsed_alignments=[two_alignment, three_alignment],
    )

    rows = calculate_protein_model_support(
        [two_model, three_model],
        [superlocus],
        ["protein"],
        terminal_cds_3prime_wobble_bp=3,
    )
    by_id = {row["output_transcript_id"]: row for row in rows}
    two_row = by_id["two_exon"]
    three_row = by_id["three_exon_unsupported"]

    assert two_row["protein_internal_cds_exon_f1"] is None
    assert two_row["protein_supported_internal_cds_exon_fraction"] is None
    assert two_row["protein_evidence_score"] == pytest.approx(1.0)
    assert "+protein_internal_cds_exon_f1:" not in two_row["protein_evidence_score_components"]
    assert "/available_weight:0.800" in two_row["protein_evidence_score_components"]

    assert three_row["protein_internal_cds_exon_f1"] == pytest.approx(0.0)
    assert three_row["protein_supported_internal_cds_exon_fraction"] == pytest.approx(0.0)
    assert "+protein_internal_cds_exon_f1:0.000000*0.200" in three_row["protein_evidence_score_components"]

    output_path = tmp_path / "protein_support.tsv"
    write_protein_model_support_tsv(rows, output_path)
    output = read_tsv(output_path)
    assert output["two_exon"]["protein_internal_cds_exon_f1"] == "."
    assert output["two_exon"]["protein_supported_internal_cds_exon_fraction"] == "."
    assert output["three_exon_unsupported"]["protein_internal_cds_exon_f1"] == "0.000000"
    assert output["three_exon_unsupported"]["protein_supported_internal_cds_exon_fraction"] == "0.000000"


def test_applicable_imperfect_protein_pairwise_score_is_unchanged() -> None:
    model_parts = cds_parts(3)
    alignment_parts = (model_parts[0], (model_parts[1][0], model_parts[1][1] - 10), model_parts[2])
    metrics = pairwise_model_protein_metrics(
        make_model(model_parts, "three_exon"),
        make_protein_alignment(alignment_parts, "protein_three"),
        wobble_bp=3,
    )

    assert metrics["internal_cds_exon_f1"] == pytest.approx(0.0)
    assert metrics["score"] == pytest.approx(0.5806034482758621)


@pytest.mark.parametrize("exon_count", [1, 2, 3, 4])
def test_perfect_transcript_pairwise_score_across_cds_topologies(exon_count: int) -> None:
    parts = cds_parts(exon_count)
    metrics = transcript_model_pair_metrics(
        make_model(parts, f"model_{exon_count}"),
        make_transcript(parts, f"transcript_{exon_count}"),
    )

    assert metrics["pairwise_transcript_score"] == pytest.approx(1.0)
    if exon_count <= 2:
        assert metrics["pairwise_internal_cds_exon_compatibility"] is None
    else:
        assert metrics["pairwise_internal_cds_exon_compatibility"] == pytest.approx(1.0)


def test_transcript_internal_exon_na_is_renormalized_and_serialized(tmp_path: Path) -> None:
    two_parts = cds_parts(2)
    three_parts = cds_parts(3, offset=2000)
    unsupported_three_parts = (three_parts[0], (three_parts[1][0], three_parts[1][1] - 10), three_parts[2])
    two_model = make_model(two_parts, "two_exon")
    three_model = make_model(three_parts, "three_exon_unsupported")
    two_tx = make_transcript(two_parts, "transcript_two", member_count=10)
    three_tx = make_transcript(unsupported_three_parts, "transcript_three", member_count=10)
    superlocus = TranscriptEvidenceSuperlocus(
        evidence_superlocus_id="transcript_superlocus",
        chrom="chr1",
        strand="+",
        start0=two_parts[0][0],
        end0=unsupported_three_parts[-1][1],
        collapsed_transcripts=[two_tx, three_tx],
    )

    rows = calculate_transcript_model_support(
        [two_model, three_model],
        [superlocus],
        ["transcript"],
    )
    by_id = {row["output_transcript_id"]: row for row in rows}
    two_row = by_id["two_exon"]
    three_row = by_id["three_exon_unsupported"]

    assert two_row["transcript_best_pairwise_score"] == pytest.approx(1.0)
    assert two_row["transcript_evidence_score"] == pytest.approx(1.0)
    assert two_row["transcript_mean_best_source_internal_cds_exon_compatibility"] is None
    assert two_row["transcript_supported_internal_cds_exon_fraction"] is None
    assert two_row["transcript_scoring_supported_internal_cds_exon_fraction"] is None

    assert three_row["transcript_mean_best_source_internal_cds_exon_compatibility"] == pytest.approx(0.0)
    assert three_row["transcript_supported_internal_cds_exon_fraction"] == pytest.approx(0.0)
    assert three_row["transcript_scoring_supported_internal_cds_exon_fraction"] == pytest.approx(0.0)

    output_path = tmp_path / "transcript_support.tsv"
    write_transcript_model_support_tsv(rows, output_path)
    output = read_tsv(output_path)
    assert output["two_exon"]["transcript_mean_best_source_internal_cds_exon_compatibility"] == "."
    assert output["two_exon"]["transcript_supported_internal_cds_exon_fraction"] == "."
    assert output["two_exon"]["transcript_scoring_supported_internal_cds_exon_fraction"] == "."
    assert output["three_exon_unsupported"]["transcript_mean_best_source_internal_cds_exon_compatibility"] == "0.000000"
    assert output["three_exon_unsupported"]["transcript_supported_internal_cds_exon_fraction"] == "0.000000"
    assert output["three_exon_unsupported"]["transcript_scoring_supported_internal_cds_exon_fraction"] == "0.000000"


def test_applicable_imperfect_transcript_pairwise_score_is_unchanged() -> None:
    model_parts = cds_parts(3)
    transcript_parts = (model_parts[0], (model_parts[1][0], model_parts[1][1] - 10), model_parts[2])
    metrics = transcript_model_pair_metrics(
        make_model(model_parts, "three_exon"),
        make_transcript(transcript_parts, "transcript_three"),
    )

    assert metrics["pairwise_internal_cds_exon_compatibility"] == pytest.approx(0.0)
    assert metrics["pairwise_transcript_score"] == pytest.approx(0.5466666666666666)
