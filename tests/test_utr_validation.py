from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.genome_sequence import CANONICAL_SPLICE_MOTIFS, IndexedFasta
from moselecta.utr_validation import (
    UTRValidationRecorder,
    validation_audit_rows,
    validate_utr_side,
)


def _fasta(tmp_path: Path, *, motif: str, strand: str = "+") -> IndexedFasta:
    sequence = list("A" * 100)
    donor, acceptor = motif.split("-")
    if strand == "+":
        sequence[20:22] = donor
        sequence[38:40] = acceptor
    else:
        complement = str.maketrans("ACGT", "TGCA")
        sequence[38:40] = donor.translate(complement)[::-1]
        sequence[20:22] = acceptor.translate(complement)[::-1]
    path = tmp_path / f"{strand}.fa"
    path.write_text(">chr1\n" + "".join(sequence) + "\n")
    return IndexedFasta(path)


def _validate(*, fasta: IndexedFasta | None = None, strand: str = "+", introns=((20, 40),)):
    return validate_utr_side(
        transcript_id="tx1",
        evidence_origin="transcript",
        evidence_id="evidence1",
        chrom="chr1",
        strand=strand,
        side="five_prime",
        introns=introns,
        minimum_intron_bp=20,
        fasta=fasta,
    )


@pytest.mark.parametrize("motif", sorted(CANONICAL_SPLICE_MOTIFS))
@pytest.mark.parametrize("strand", ["+", "-"])
def test_authoritative_motifs_pass_on_both_strands(
    tmp_path: Path, motif: str, strand: str
) -> None:
    with _fasta(tmp_path, motif=motif, strand=strand) as fasta:
        result = _validate(fasta=fasta, strand=strand)
    assert result.status == "passed"
    assert result.observations[0].motif == motif
    assert result.observations[0].motif_status == "passed"


def test_resolved_unrecognized_motif_rejects(tmp_path: Path) -> None:
    with _fasta(tmp_path, motif="GG-TT") as fasta:
        result = _validate(fasta=fasta)
    assert result.status == "rejected"
    assert result.reason_codes == ("resolved_unrecognized_motif",)


@pytest.mark.parametrize("length", [1, 3, 5, 19])
def test_all_subminimum_introns_use_one_reason(length: int) -> None:
    result = _validate(introns=((20, 20 + length),))
    assert result.status == "rejected"
    assert "short_utr_intron" in result.reason_codes
    assert result.observations[0].intron_length == length


@pytest.mark.parametrize("length", [20, 21])
def test_minimum_is_inclusive_without_genome(length: int) -> None:
    result = _validate(introns=((20, 20 + length),))
    assert result.status == "not_fully_evaluated"
    assert "short_utr_intron" not in result.reason_codes
    assert result.reason_codes == ("no_genome_fasta",)


def test_unspliced_side_has_no_motif_test() -> None:
    result = _validate(introns=())
    assert result.status == "passed"
    assert result.observations == ()
    assert result.reason_codes == ()


def test_confirmed_failure_precedes_unevaluable_context(tmp_path: Path) -> None:
    with _fasta(tmp_path, motif="GG-TT") as fasta:
        result = _validate(fasta=fasta, introns=((20, 40), (80, 120)))
    assert result.status == "rejected"
    assert "resolved_unrecognized_motif" in result.reason_codes
    assert "splice_context_out_of_range" in result.reason_codes


def test_ambiguous_motif_is_retained_as_unevaluated(tmp_path: Path) -> None:
    with _fasta(tmp_path, motif="GN-AG") as fasta:
        result = _validate(fasta=fasta)
    assert result.status == "not_fully_evaluated"
    assert result.reason_codes == ("ambiguous_splice_sequence",)


def test_missing_contig_is_retained_as_unevaluated(tmp_path: Path) -> None:
    fasta = _fasta(tmp_path, motif="GT-AG")
    with fasta:
        result = validate_utr_side(
            transcript_id="tx", evidence_origin="prediction", evidence_id="p",
            chrom="missing", strand="+", side="three_prime",
            introns=((20, 40),), minimum_intron_bp=20, fasta=fasta,
        )
    assert result.status == "not_fully_evaluated"
    assert result.reason_codes == ("genome_contig_missing",)


def test_invalid_geometry_is_deterministically_rejected() -> None:
    result = _validate(introns=((40, 20),))
    assert result.status == "rejected"
    assert result.reason_codes == (
        "invalid_utr_geometry",
        "no_genome_fasta",
        "short_utr_intron",
    )


def test_audit_is_exception_only_and_omits_routine_no_fasta() -> None:
    assert validation_audit_rows(
        _validate(introns=((20, 40),)), rejected_action="exclude_candidate_side"
    ) == []
    assert validation_audit_rows(
        _validate(introns=()), rejected_action="exclude_candidate_side"
    ) == []

    rejected = _validate(introns=((20, 25), (40, 60)))
    rows = validation_audit_rows(
        rejected, rejected_action="exclude_candidate_side"
    )
    assert len(rows) == 1
    assert rows[0]["intron_start"] == 21
    assert rows[0]["intron_end"] == 25
    assert rows[0]["side_reason_codes"] == "no_genome_fasta,short_utr_intron"
    assert rows[0]["intron_reason_codes"] == "no_genome_fasta,short_utr_intron"
    assert rows[0]["final_action"] == "exclude_candidate_side"


def test_side_geometry_reason_is_distinct_from_junction_reasons() -> None:
    validation = validate_utr_side(
        transcript_id="tx1",
        evidence_origin="prediction",
        evidence_id="p1",
        chrom="chr1",
        strand="+",
        side="three_prime",
        introns=(),
        minimum_intron_bp=20,
        geometry_valid=False,
    )
    rows = validation_audit_rows(
        validation, rejected_action="suppress_source_side"
    )
    assert len(rows) == 1
    assert rows[0]["exception_scope"] == "side"
    assert rows[0]["side_reason_codes"] == "invalid_utr_geometry"
    assert rows[0]["intron_reason_codes"] == "."
    assert rows[0]["motif_status"] == "not_applicable"


def test_recorder_reconciles_summary_and_streams_only_exceptions(
    tmp_path: Path,
) -> None:
    detail_path = tmp_path / "exceptions.tsv"
    recorder = UTRValidationRecorder(detail_path)
    recorder.record(
        _validate(introns=()), rejected_action="exclude_candidate_side"
    )
    recorder.record(
        _validate(introns=((20, 40),)),
        rejected_action="exclude_candidate_side",
    )
    recorder.record(
        _validate(introns=((20, 25), (40, 60))),
        rejected_action="exclude_candidate_side",
    )
    summary = recorder.summary_rows()[0]
    recorder.close()

    assert summary["sides_assessed"] == 3
    assert summary["sides_accepted"] == 1
    assert summary["sides_rejected"] == 1
    assert summary["sides_not_fully_evaluated"] == 1
    assert summary["sides_without_splice_junction"] == 1
    assert summary["junctions_assessed"] == 3
    assert summary["individual_rejection_violations"] == 1
    assert summary["rejection_reason_counts"] == "short_utr_intron:1"
    assert summary["unevaluable_reason_counts"] == "no_genome_fasta:3"
    assert recorder.detail_rows_written == 1
    assert len(detail_path.read_text().splitlines()) == 2
    assert not any(isinstance(value, list) for value in recorder.__dict__.values())


def test_exception_observation_mode_does_not_retain_passing_junctions(
    tmp_path: Path,
) -> None:
    with _fasta(tmp_path, motif="GT-AG") as fasta:
        result = validate_utr_side(
            transcript_id="tx1",
            evidence_origin="transcript",
            evidence_id="evidence1",
            chrom="chr1",
            strand="+",
            side="five_prime",
            introns=((20, 40),),
            minimum_intron_bp=20,
            fasta=fasta,
            observation_mode="exceptions",
        )
    assert result.status == "passed"
    assert result.observations == ()
    assert result.junction_count == 1
    assert result.motif_passed_count == 1
