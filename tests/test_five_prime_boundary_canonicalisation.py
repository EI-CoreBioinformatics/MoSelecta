from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from moselecta.annotation_outputs import public_transcript_attributes
from moselecta.cds_phase import (
    CANONICALISED_FIVE_PRIME_BOUNDARY,
    calculate_cds_phase_chain,
    cds_phase_output_status,
)
from moselecta.five_prime_boundary import (
    AUDIT_COLUMNS,
    canonicalise_model_five_prime_boundary,
    compact_five_prime_boundary_provenance,
    five_prime_boundary_audit_row,
    write_five_prime_boundary_actions_tsv,
    write_five_prime_boundary_audit_tsv,
    write_five_prime_boundary_summary_tsv,
)
from moselecta.genome_sequence import calculate_model_sequence_integrity, reverse_complement
from moselecta.final_annotation import AnnotationSummary, FinalGene, FinalTranscript
from moselecta.models import CodingModel, FeatureRecord
from moselecta.phase_recovery import phase_recovery_audit_row, recover_model_cds_phase
from moselecta.terminal_stop import canonicalise_model_terminal_stop


class MemoryFasta:
    def __init__(self, sequence: str):
        self.sequence = sequence

    def contig_length(self, name: str) -> int | None:
        return len(self.sequence) if name == "chr1" else None

    def fetch(self, name: str, start0: int, end0: int) -> str:
        if name != "chr1" or start0 < 0 or end0 > len(self.sequence):
            raise ValueError("out of range")
        return self.sequence[start0:end0]


def _record(
    kind: str, part: tuple[int, int], *, strand: str, phase: str = ".", index: int = 1
) -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1", source="test", feature_type=kind,
        start=part[0] + 1, end=part[1], score=".", strand=strand, phase=phase,
        attrs={"ID": f"t.{kind}.{index}", "Parent": "t"}, raw_attrs=".",
        input_file=Path("input.gff3"), source_name="test", line_number=index,
        record_index=index, raw_line="",
    )


def _case(
    offset: int = 1,
    *,
    strand: str = "+",
    phase: str = ".",
    adjacent_stop: bool = False,
    exon_inside: bool = True,
    utr: str = "none",
    part_lengths: tuple[int, ...] | None = None,
    sequence: str | None = None,
) -> tuple[CodingModel, MemoryFasta]:
    coding = "ATGCCC" + ("" if adjacent_stop else "TAA")
    raw = "C" * offset + coding if sequence is None else sequence
    if part_lengths is None:
        part_lengths = (len(raw),)
    assert sum(part_lengths) == len(raw)
    chars = list("C" * 400)
    chunks: list[tuple[int, int]] = []
    cursor = 0
    if strand == "+":
        start = 100
        for index, length in enumerate(part_lengths):
            part = (start + index * 30, start + index * 30 + length)
            chars[part[0]:part[1]] = raw[cursor:cursor + length]
            chunks.append(part)
            cursor += length
        if adjacent_stop:
            chars[chunks[-1][1]:chunks[-1][1] + 3] = "TAA"
        first = chunks[0]
        exon = (first[0] - (5 if exon_inside else 0), first[1])
        if len(chunks) == 1 and adjacent_stop:
            exon = (exon[0], first[1] + 3)
        utr_part = (first[0] - 5, first[0]) if utr == "adjacent" else (first[0] - 8, first[0] - 3)
    else:
        high = 300
        for index, length in enumerate(part_lengths):
            end = high - index * 30
            part = (end - length, end)
            chars[part[0]:part[1]] = reverse_complement(raw[cursor:cursor + length])
            chunks.append(part)
            cursor += length
        if adjacent_stop:
            chars[chunks[-1][0] - 3:chunks[-1][0]] = reverse_complement("TAA")
        first = chunks[0]
        exon = (first[0], first[1] + (5 if exon_inside else 0))
        if len(chunks) == 1 and adjacent_stop:
            exon = (first[0] - 3, exon[1])
        utr_part = (first[1], first[1] + 5) if utr == "adjacent" else (first[1] + 3, first[1] + 8)
    children = [_record("exon", exon, strand=strand)]
    for index, part in enumerate(chunks, start=1):
        children.append(_record("CDS", part, strand=strand, phase=phase, index=index))
    if utr != "none":
        children.append(_record("five_prime_UTR", utr_part, strand=strand))
    model = CodingModel(
        source_name="test", input_file=Path("input.gff3"), input_gene_id="g",
        input_transcript_id="t", output_gene_id="g.out", output_transcript_id="t.out",
        chrom="chr1", strand=strand, cds_parts=sorted(chunks), gene_record=None,
        transcript_record=None, child_records=children,
    )
    return model, MemoryFasta("".join(chars))


def _cds_records(model: CodingModel) -> list[FeatureRecord]:
    return [record for record in model.child_records if record.feature_type == "CDS"]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [1, 2])
def test_class_b_canonicalises_both_trims_and_strands(strand: str, offset: int) -> None:
    model, fasta = _case(offset, strand=strand)
    old_parts = list(model.cds_parts)
    old_exons = [(r.start0, r.end0) for r in model.child_records if r.feature_type == "exon"]
    result = canonicalise_model_five_prime_boundary(model, fasta)
    audit = result.five_prime_boundary_canonicalisation
    assert audit is not None and audit.status == "canonicalised"
    assert audit.bases_trimmed == offset
    assert result.cds_length == model.cds_length - offset
    assert [(r.start0, r.end0) for r in result.child_records if r.feature_type == "exon"] == old_exons
    expected_first = (
        (old_parts[0][0] + offset, old_parts[0][1])
        if strand == "+" else (old_parts[-1][0], old_parts[-1][1] - offset)
    )
    assert (expected_first in result.cds_parts)
    assert all(record.phase == "0" for record in _cds_records(result))
    assert cds_phase_output_status(result) == CANONICALISED_FIVE_PRIME_BOUNDARY
    assert "source_phase:." in compact_five_prime_boundary_provenance(result)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_multi_exon_only_moves_first_part_and_propagates_phase(strand: str) -> None:
    model, fasta = _case(1, strand=strand, part_lengths=(5, 5))
    old_biological = sorted(model.cds_parts, reverse=strand == "-")
    result = canonicalise_model_five_prime_boundary(model, fasta)
    new_biological = sorted(result.cds_parts, reverse=strand == "-")
    assert new_biological[1:] == old_biological[1:]
    assert new_biological[0][1] - new_biological[0][0] == 4
    expected = calculate_cds_phase_chain(result.cds_parts, strand, 0)
    assert {(r.start0, r.end0): r.phase for r in _cds_records(result)} == expected


def test_trim_preserves_source_gene_transcript_and_exon_geometry() -> None:
    model, fasta = _case(1)
    exon = next(record for record in model.child_records if record.feature_type == "exon")
    model.gene_record = _record("gene", (exon.start0, exon.end0 + 20), strand="+")
    model.transcript_record = _record("mRNA", (exon.start0, exon.end0 + 10), strand="+")
    original_gene = (model.gene_record.start0, model.gene_record.end0)
    original_transcript = (model.transcript_record.start0, model.transcript_record.end0)
    original_exons = [(exon.start0, exon.end0)]
    result = canonicalise_model_five_prime_boundary(model, fasta)
    assert (result.gene_record.start0, result.gene_record.end0) == original_gene
    assert (result.transcript_record.start0, result.transcript_record.end0) == original_transcript
    assert [
        (record.start0, record.end0)
        for record in result.child_records if record.feature_type == "exon"
    ] == original_exons


@pytest.mark.parametrize("offset", [1, 2])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_adjacent_five_prime_utr_extends_exactly_trim(strand: str, offset: int) -> None:
    model, fasta = _case(offset, strand=strand, utr="adjacent")
    old = next(r for r in model.child_records if r.feature_type == "five_prime_UTR")
    result = canonicalise_model_five_prime_boundary(model, fasta)
    new = next(r for r in result.child_records if r.feature_type == "five_prime_UTR")
    assert new.end0 - new.start0 == old.end0 - old.start0 + offset
    assert result.five_prime_boundary_canonicalisation.five_prime_utr_extended
    first = sorted(result.cds_parts, reverse=strand == "-")[0]
    assert (new.end0 == first[0]) if strand == "+" else (new.start0 == first[1])


@pytest.mark.parametrize("utr", ["none", "gap"])
def test_missing_or_gapped_utr_is_not_synthesised_or_changed(utr: str) -> None:
    model, fasta = _case(1, utr=utr)
    old_utrs = [(r.start0, r.end0) for r in model.child_records if r.feature_type == "five_prime_UTR"]
    result = canonicalise_model_five_prime_boundary(model, fasta)
    assert [(r.start0, r.end0) for r in result.child_records if r.feature_type == "five_prime_UTR"] == old_utrs
    assert not result.five_prime_boundary_canonicalisation.five_prime_utr_extended


def test_adjacent_stop_then_uses_existing_terminal_stop_plus_three() -> None:
    model, fasta = _case(1, adjacent_stop=True)
    corrected = canonicalise_model_five_prime_boundary(model, fasta)
    assert corrected.five_prime_boundary_canonicalisation.stop_representation == "adjacent"
    final = canonicalise_model_terminal_stop(corrected, fasta)
    assert final.terminal_stop_canonicalisation.status == "canonicalised"
    assert final.terminal_stop_canonicalisation.offset_bp == 3
    assert final.five_prime_boundary_canonicalisation.status == "canonicalised"


def test_included_stop_becomes_already_canonical_and_integrity_start_present() -> None:
    model, fasta = _case(2)
    corrected = canonicalise_model_five_prime_boundary(model, fasta)
    terminal = canonicalise_model_terminal_stop(corrected, fasta)
    assert terminal.terminal_stop_canonicalisation.status == "already_canonical"
    row = calculate_model_sequence_integrity(
        terminal, window=fasta.sequence, window_start0=0,
        contig_length=len(fasta.sequence),
    )
    assert row["sequence_start_codon_status"] == "present"
    assert row["sequence_stop_codon_status"] == "included_in_cds"
    assert row["sequence_internal_stop_count"] == 0


def test_correction_removes_phase_zero_false_internal_stop() -> None:
    model, fasta = _case(1, sequence="CATGCCTAACTAA")
    before = calculate_model_sequence_integrity(
        model, window=fasta.sequence, window_start0=0,
        contig_length=len(fasta.sequence),
    )
    corrected = canonicalise_model_five_prime_boundary(model, fasta)
    after = calculate_model_sequence_integrity(
        corrected, window=fasta.sequence, window_start0=0,
        contig_length=len(fasta.sequence),
    )
    assert before["sequence_internal_stop_count"] == 1
    assert after["sequence_internal_stop_count"] == 0


@pytest.mark.parametrize(
    ("kwargs", "reason", "status"),
    [
        ({"exon_inside": False}, "cds_starts_at_exon_boundary", "skipped"),
        ({"phase": "1"}, "numeric_source_phase", "not_eligible_numeric_source_phase"),
        ({"sequence": "ATGCCCTAA"}, "phase_zero_positively_supported", "skipped"),
        ({"sequence": "CCCCCCCCCC"}, "no_unique_start_stop_offset", "skipped"),
        ({"sequence": "CATGCCCCCC"}, "no_unique_start_stop_offset", "skipped"),
        ({"sequence": "CATGTAACCCTAA"}, "offset_frame_internal_stop", "skipped"),
        ({"sequence": "CATGN CCTAA".replace(" ", "")}, "ambiguous_sequence", "skipped"),
    ],
)
def test_negative_eligibility_is_stable(kwargs: dict[str, object], reason: str, status: str) -> None:
    model, fasta = _case(1, **kwargs)
    old_parts = list(model.cds_parts)
    result = canonicalise_model_five_prime_boundary(model, fasta)
    audit = result.five_prime_boundary_canonicalisation
    assert audit.status == status
    assert audit.skip_reason == reason
    assert result.cds_parts == old_parts


def test_start_codon_cannot_cross_first_cds_splice_boundary() -> None:
    model, fasta = _case(2, part_lengths=(4, 7))
    result = canonicalise_model_five_prime_boundary(model, fasta)
    assert result.five_prime_boundary_canonicalisation.skip_reason == "start_codon_not_in_first_cds_part"


def test_boundary_canonicalised_model_is_not_phase_recovered() -> None:
    model, fasta = _case(1)
    corrected = canonicalise_model_five_prime_boundary(model, fasta)
    result = recover_model_cds_phase(corrected, fasta)
    assert result is corrected
    assert result.cds_phase_recovery is None
    row = phase_recovery_audit_row(result)
    assert row["status"] == CANONICALISED_FIVE_PRIME_BOUNDARY
    assert row["original_source_first_phase"] == "."


def test_final_transcript_attributes_keep_dedicated_boundary_provenance() -> None:
    model, fasta = _case(1)
    corrected = canonicalise_model_five_prime_boundary(model, fasta)
    transcript = FinalTranscript(
        "MoSelecta_000001.1", corrected,
        annotation_summary=AnnotationSummary(),
    )
    gene = FinalGene("MoSelecta_000001", (transcript,))
    attrs = public_transcript_attributes(
        gene, transcript, include_alternative_metadata=False
    )
    assert attrs["cds_phase_output_status"] == CANONICALISED_FIVE_PRIME_BOUNDARY
    assert attrs["cds_5prime_boundary_canonicalised"].startswith("trim:1|")
    assert "cds_phase_recovered" not in attrs


def test_failed_post_trim_validation_rolls_back_with_stable_reason(monkeypatch: pytest.MonkeyPatch) -> None:
    model, fasta = _case(1)
    old_parts = list(model.cds_parts)
    monkeypatch.setattr(
        "moselecta.five_prime_boundary._validate_canonical_model",
        lambda *args, **kwargs: (None, "postcheck_test_failure"),
    )
    result = canonicalise_model_five_prime_boundary(model, fasta)
    assert result.cds_parts == old_parts
    assert result.five_prime_boundary_canonicalisation.status == "skipped"
    assert result.five_prime_boundary_canonicalisation.skip_reason == "postcheck_test_failure"


def test_audit_and_provenance_are_deterministic(tmp_path: Path) -> None:
    model, fasta = _case(1, utr="adjacent")
    result = canonicalise_model_five_prime_boundary(model, fasta)
    row = five_prime_boundary_audit_row(result)
    assert row["status"] == "canonicalised"
    assert row["bases_trimmed"] == 1
    assert row["exon_context"] == "strictly_inside_source_exon"
    assert row["original_source_first_phase"] == "."
    path = tmp_path / "audit.tsv"
    write_five_prime_boundary_audit_tsv([result], path)
    lines = path.read_text().splitlines()
    assert lines[0].split("\t") == AUDIT_COLUMNS
    assert len(lines) == 2


def test_summary_actions_and_exhaustive_audit_reconcile(tmp_path: Path) -> None:
    changed_model, changed_fasta = _case(1, utr="adjacent")
    skipped_model, skipped_fasta = _case(1, exon_inside=False)
    changed = canonicalise_model_five_prime_boundary(changed_model, changed_fasta)
    skipped = canonicalise_model_five_prime_boundary(skipped_model, skipped_fasta)
    models = [changed, skipped]
    audit_path = tmp_path / "all.tsv"
    actions_path = tmp_path / "actions.tsv"
    summary_path = tmp_path / "summary.tsv"
    write_five_prime_boundary_audit_tsv(models, audit_path)
    write_five_prime_boundary_actions_tsv(models, actions_path)
    write_five_prime_boundary_summary_tsv(
        models, summary_path, execution_status="performed"
    )
    with audit_path.open(newline="") as handle:
        audit = list(csv.DictReader(handle, delimiter="\t"))
    with actions_path.open(newline="") as handle:
        actions = list(csv.DictReader(handle, delimiter="\t"))
    with summary_path.open(newline="") as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    assert actions == [row for row in audit if row["status"] == "canonicalised"]
    assert len(actions) == 1
    assert actions[0]["five_prime_utr_extended"] == "true"
    assert summary["canonicalised_models"] == "1"
    assert summary["skipped_models"] == "1"
    assert summary["canonicalised_with_utr_extension"] == "1"
    assert json.loads(summary["canonicalised_bases_trimmed_counts"]) == {"1": 1}
    assert json.loads(summary["skip_reason_counts"]) == {
        "cds_starts_at_exon_boundary": 1
    }
