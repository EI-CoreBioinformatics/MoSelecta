from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from moselecta.annotation_outputs import (
    FinalAnnotationOutputMetadata,
    public_transcript_attributes,
    write_final_annotation_gff3,
)
from moselecta.cds_refinement import CDSTerminusRefinementSettings, refine_cds_termini
from moselecta.final_annotation import AnnotationSummary, FinalAnnotation, FinalGene, FinalTranscript
from moselecta.genome_sequence import (
    calculate_model_sequence_integrity,
    reverse_complement,
    sequence_integrity_metrics_from_row,
)
from moselecta.models import CodingModel, FeatureRecord
from moselecta.structural import cds_exon_chains_match
from moselecta.terminal_stop import (
    annotate_new_same_strand_overlaps,
    canonicalise_model_terminal_stop,
    terminal_stop_audit_row,
    write_terminal_stop_actions_tsv,
    write_terminal_stop_audit_tsv,
    write_terminal_stop_summary_tsv,
)


class MemoryFasta:
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()

    def contig_length(self, name: str) -> int | None:
        return len(self.sequence) if name == "chr1" else None

    def fetch(self, name: str, start0: int, end0: int) -> str:
        if name != "chr1" or start0 < 0 or end0 > len(self.sequence):
            raise ValueError("out of range")
        return self.sequence[start0:end0]


def _record(
    feature_type: str,
    start0: int,
    end0: int,
    *,
    strand: str,
    phase: str = ".",
    name: str = "m",
) -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1", source="test", feature_type=feature_type,
        start=start0 + 1, end=end0, score=".", strand=strand, phase=phase,
        attrs={"ID": f"{name}.{feature_type}", "Parent": f"{name}.t1"},
        raw_attrs=".", input_file=Path("input.gff3"), source_name="test",
        line_number=1, record_index=1, raw_line="",
    )


def _case(
    offset: int,
    strand: str,
    *,
    extension_mode: str = "within",
    phase: str = "0",
    coding: str = "ATGCCCGGGTAA",
    name: str = "m",
    adjacent_override: str | None = None,
) -> tuple[CodingModel, MemoryFasta, tuple[int, int], tuple[int, int]]:
    phase_value = int(phase) if phase in {"0", "1", "2"} else 0
    prefix = "C" * phase_value
    translated_source = (
        coding[:-offset] if offset > 0 else coding + "C" * (-offset)
    )
    # For +3/+2/+1 this expression removes 3/2/1 terminal bases respectively.
    raw_source = prefix + translated_source
    if offset > 0:
        tail = coding[len(translated_source):]
    else:
        tail = "CCC"
    if adjacent_override is not None:
        tail = adjacent_override + tail[len(adjacent_override):]
    transcript_space = raw_source + tail + "CCCCCC"
    genomic_payload = transcript_space if strand == "+" else reverse_complement(transcript_space)
    anchor = 100
    chars = list("C" * 300)
    payload_start = anchor if strand == "+" else 200 - len(genomic_payload)
    chars[payload_start:payload_start + len(genomic_payload)] = genomic_payload
    fasta = MemoryFasta("".join(chars))

    if strand == "+":
        cds = (anchor, anchor + len(raw_source))
        canonical = (anchor, anchor + len(prefix) + len(coding))
    else:
        high = 200
        cds = (high - len(raw_source), high)
        canonical = (high - len(prefix) - len(coding), high)
    if offset > 0 and extension_mode == "within":
        exon = canonical
    else:
        exon = cds
    children = [
        _record("exon", *exon, strand=strand, name=name),
        _record("CDS", *cds, strand=strand, phase=phase, name=name),
    ]
    model = CodingModel(
        source_name="test", input_file=Path("input.gff3"),
        input_gene_id=name, input_transcript_id=f"{name}.t1",
        output_gene_id=f"{name}.gene", output_transcript_id=f"{name}.mrna",
        chrom="chr1", strand=strand, cds_parts=[cds], gene_record=None,
        transcript_record=None, child_records=children,
    )
    return model, fasta, canonical, exon


def _integrity(model: CodingModel, fasta: MemoryFasta) -> dict[str, object]:
    return calculate_model_sequence_integrity(
        model, window=fasta.sequence, window_start0=0,
        contig_length=len(fasta.sequence),
    )


def _set_exon(model: CodingModel, part: tuple[int, int]) -> None:
    model.child_records = [
        _record("exon", *part, strand=model.strand)
        if record.feature_type == "exon"
        else record
        for record in model.child_records
    ]


def _add_adjacent_three_prime_utr(
    model: CodingModel,
    canonical: tuple[int, int],
    *,
    completely_consumed: bool = False,
) -> tuple[int, int]:
    old = model.cds_parts[0]
    if model.strand == "+":
        utr = (old[1], canonical[1] if completely_consumed else canonical[1] + 6)
        exon = (old[0], utr[1])
    else:
        utr = (canonical[0] if completely_consumed else canonical[0] - 6, old[0])
        exon = (utr[0], old[1])
    _set_exon(model, exon)
    model.child_records.append(_record("three_prime_UTR", *utr, strand=model.strand))
    return utr


def _three_prime_utrs(model: CodingModel) -> list[tuple[int, int]]:
    return [
        (record.start0, record.end0)
        for record in model.child_records
        if record.feature_type.lower() == "three_prime_utr"
    ]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_already_included_stop_is_unchanged(strand: str) -> None:
    model, fasta, canonical, _ = _case(0, strand)
    result = canonicalise_model_terminal_stop(model, fasta)
    assert result is model
    assert result.cds_parts == [canonical]
    assert result.terminal_stop_canonicalisation.status == "already_canonical"
    assert _integrity(result, fasta)["sequence_stop_codon_status"] == "included_in_cds"


def test_summary_actions_and_exhaustive_audit_have_stable_semantics(
    tmp_path: Path,
) -> None:
    changed, changed_fasta, _canonical, _ = _case(3, "+", name="changed")
    unchanged, unchanged_fasta, _canonical, _ = _case(0, "+", name="unchanged")
    changed = canonicalise_model_terminal_stop(changed, changed_fasta)
    unchanged = canonicalise_model_terminal_stop(unchanged, unchanged_fasta)
    models = [changed, unchanged]
    audit_path = tmp_path / "audit.tsv"
    actions_path = tmp_path / "actions.tsv"
    summary_path = tmp_path / "summary.tsv"
    write_terminal_stop_audit_tsv(models, audit_path)
    write_terminal_stop_actions_tsv(models, actions_path)
    write_terminal_stop_summary_tsv(
        models, summary_path, execution_status="performed"
    )
    with audit_path.open(newline="") as handle:
        audit = list(csv.DictReader(handle, delimiter="\t"))
    with actions_path.open(newline="") as handle:
        actions = list(csv.DictReader(handle, delimiter="\t"))
    with summary_path.open(newline="") as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    assert [row for row in audit if row["status"] == "canonicalised"] == actions
    assert {row["status"] for row in audit} == {"canonicalised", "already_canonical"}
    assert summary["total_models"] == "2"
    assert summary["canonicalised_models"] == "1"
    assert summary["already_canonical_models"] == "1"
    assert json.loads(summary["canonicalised_offset_counts"]) == {"+3": 1}


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [3, 2, 1])
@pytest.mark.parametrize("extension_mode", ["within", "extend_exon"])
def test_positive_offset_matrix_both_strands(
    strand: str, offset: int, extension_mode: str
) -> None:
    model, fasta, canonical, original_exon = _case(
        offset, strand, extension_mode=extension_mode
    )
    original_phase = next(r.phase for r in model.child_records if r.feature_type == "CDS")
    result = canonicalise_model_terminal_stop(model, fasta)
    audit = result.terminal_stop_canonicalisation
    assert result.cds_parts == [canonical]
    exons = [(r.start0, r.end0) for r in result.child_records if r.feature_type == "exon"]
    assert exons == [canonical]
    assert result.model_start0 == canonical[0]
    assert result.model_end0 == canonical[1]
    assert next(r.phase for r in result.child_records if r.feature_type == "CDS") == original_phase
    assert audit.status == "canonicalised"
    assert audit.offset_bp == offset
    assert audit.stop_codon == "TAA"
    assert audit.exon_adjusted is (extension_mode == "extend_exon")
    metrics = _integrity(result, fasta)
    assert metrics["sequence_stop_codon_status"] == "included_in_cds"
    assert metrics["sequence_internal_stop_count"] == 0
    row = terminal_stop_audit_row(result)
    assert row["offset"] == f"+{offset}"
    assert row["original_biological_3prime_cds_boundary"] != row[
        "canonical_biological_3prime_cds_boundary"
    ]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [-1, -2])
def test_negative_offset_trim_matrix_both_strands(strand: str, offset: int) -> None:
    model, fasta, canonical, original_exon = _case(offset, strand)
    result = canonicalise_model_terminal_stop(model, fasta)
    assert result.cds_parts == [canonical]
    assert [(r.start0, r.end0) for r in result.child_records if r.feature_type == "exon"] == [original_exon]
    assert result.model_start0 == original_exon[0]
    assert result.model_end0 == original_exon[1]
    assert result.terminal_stop_canonicalisation.offset_bp == offset
    assert _integrity(result, fasta)["sequence_stop_codon_status"] == "included_in_cds"
    assert _integrity(result, fasta)["sequence_internal_stop_count"] == 0


@pytest.mark.parametrize("phase", ["0", "1", "2"])
def test_numeric_phase_zero_one_two_establishes_frame(phase: str) -> None:
    model, fasta, canonical, _ = _case(3, "+", phase=phase)
    result = canonicalise_model_terminal_stop(model, fasta)
    assert result.cds_parts == [canonical]
    assert result.terminal_stop_canonicalisation.frame_confidence_basis == f"numeric_phase_{phase}"
    assert next(r.phase for r in result.child_records if r.feature_type == "CDS") == phase


def test_missing_phase_requires_confirmed_zero_frame_start() -> None:
    allowed, fasta, canonical, _ = _case(3, "+", phase=".")
    result = canonicalise_model_terminal_stop(allowed, fasta)
    assert result.cds_parts == [canonical]
    assert result.terminal_stop_canonicalisation.frame_confidence_basis == "start_confirmed_zero_frame"

    missing, missing_fasta, _, _ = _case(3, "+", phase=".", coding="CCCCCCGGGTAA")
    skipped = canonicalise_model_terminal_stop(missing, missing_fasta)
    assert skipped.cds_parts == missing.cds_parts
    assert skipped.terminal_stop_canonicalisation.skip_reason == "uncertain_frame_missing_start"

    ambiguous, ambiguous_fasta, _, _ = _case(3, "+", phase=".", coding="NTGCCCGGGTAA")
    skipped = canonicalise_model_terminal_stop(ambiguous, ambiguous_fasta)
    assert skipped.terminal_stop_canonicalisation.skip_reason == "uncertain_frame_ambiguous_start"


def test_ambiguous_completion_and_trim_is_skipped() -> None:
    # Existing terminal TAA plus the first T of a second TAA; adjacent AA makes
    # both -1 trimming and +2 completion locally valid.
    model, fasta, _, _ = _case(-1, "+", adjacent_override="AA")
    chars = list(fasta.sequence)
    chars[model.cds_end0 - 1] = "T"
    chars[model.cds_end0:model.cds_end0 + 2] = "AA"
    fasta = MemoryFasta("".join(chars))
    result = canonicalise_model_terminal_stop(model, fasta)
    assert result is model
    assert result.terminal_stop_canonicalisation.skip_reason == "ambiguous_partial_completion_vs_overrun_trim"


def test_adjacent_n_and_contig_boundary_are_skipped() -> None:
    model, fasta, _, _ = _case(3, "+")
    chars = list(fasta.sequence)
    chars[model.cds_end0:model.cds_end0 + 3] = "TNA"
    result = canonicalise_model_terminal_stop(model, MemoryFasta("".join(chars)))
    assert result.terminal_stop_canonicalisation.skip_reason == "ambiguous_genomic_bases"

    boundary_model, boundary_fasta, _, _ = _case(3, "+")
    boundary_fasta.sequence = boundary_fasta.sequence[:boundary_model.cds_end0]
    result = canonicalise_model_terminal_stop(boundary_model, boundary_fasta)
    assert result.terminal_stop_canonicalisation.skip_reason == "contig_boundary"


def test_splice_crossing_downstream_exon_and_missing_exon_are_skipped() -> None:
    crossing, fasta, _, _ = _case(3, "+", extension_mode="extend_exon")
    exon = next(r for r in crossing.child_records if r.feature_type == "exon")
    exon.end += 1
    result = canonicalise_model_terminal_stop(crossing, fasta)
    assert result.terminal_stop_canonicalisation.skip_reason == "splice_crossing_or_structurally_nonlocal"

    downstream, fasta, _, _ = _case(3, "+", extension_mode="extend_exon")
    downstream.child_records.append(
        _record("exon", downstream.cds_end0 + 10, downstream.cds_end0 + 20, strand="+", name="m")
    )
    result = canonicalise_model_terminal_stop(downstream, fasta)
    assert result.terminal_stop_canonicalisation.skip_reason == "downstream_exon_or_structurally_nonlocal"

    no_exon, fasta, _, _ = _case(3, "+")
    no_exon.child_records = [r for r in no_exon.child_records if r.feature_type == "CDS"]
    result = canonicalise_model_terminal_stop(no_exon, fasta)
    assert result.terminal_stop_canonicalisation.skip_reason == "missing_terminal_exon"


def test_internal_stop_does_not_veto_canonicalisation() -> None:
    model, fasta, canonical, _ = _case(3, "+", coding="ATGTAACCCTAA")
    result = canonicalise_model_terminal_stop(model, fasta)
    assert result.cds_parts == [canonical]
    metrics = _integrity(result, fasta)
    assert metrics["sequence_stop_codon_status"] == "included_in_cds"
    assert metrics["sequence_internal_stop_count"] == 1

    overrun, fasta, canonical, _ = _case(-2, "+", coding="ATGTAACCCTAA")
    result = canonicalise_model_terminal_stop(overrun, fasta)
    assert result.cds_parts == [canonical]
    assert _integrity(result, fasta)["sequence_internal_stop_count"] == 1


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [3, 2, 1, -1, -2])
def test_canonicalised_pairs_are_exact_without_wobble(strand: str, offset: int) -> None:
    variant, fasta, _, _ = _case(offset, strand, name="variant")
    canonical, _, _, _ = _case(0, strand, name="canonical")
    result = canonicalise_model_terminal_stop(variant, fasta)
    assert cds_exon_chains_match(result, canonical, wobble_bp=0)


def test_new_candidate_overlap_is_audited_but_not_vetoed() -> None:
    model, fasta, _, _ = _case(3, "+", name="left")
    result = canonicalise_model_terminal_stop(model, fasta)
    other = CodingModel(
        source_name="test", input_file=Path("other.gff3"), input_gene_id="right",
        input_transcript_id="right.t1", output_gene_id="right.gene",
        output_transcript_id="right.mrna", chrom="chr1", strand="+",
        cds_parts=[(result.cds_end0 - 2, result.cds_end0 + 10)],
        gene_record=None, transcript_record=None, child_records=[],
    )
    annotate_new_same_strand_overlaps([result, other])
    audit = result.terminal_stop_canonicalisation
    assert audit.status == "canonicalised"
    assert audit.newly_introduced_same_strand_cds_overlap_count == 1
    assert audit.newly_introduced_same_strand_cds_overlap_models == "right.mrna"


def test_canonicalised_stop_is_not_a_postselection_refinement_target() -> None:
    model, fasta, _, _ = _case(3, "+")
    result = canonicalise_model_terminal_stop(model, fasta)
    result.sequence_integrity = sequence_integrity_metrics_from_row(
        _integrity(result, fasta)
    )
    transcript = FinalTranscript(
        "MoSelecta_000001.1", result,
        annotation_summary=AnnotationSummary(
            sequence_integrity_status=result.sequence_integrity.status
        ),
    )
    gene = FinalGene("MoSelecta_000001", (transcript,))
    annotation = FinalAnnotation((gene,))
    refinement = refine_cds_termini(
        annotation, [], [], fasta, CDSTerminusRefinementSettings()
    )
    assert refinement.summary.execution_status == "zero_incomplete_models"
    assert refinement.annotation is annotation

    attrs = public_transcript_attributes(
        gene, transcript, include_alternative_metadata=False
    )
    assert attrs["terminal_stop_canonicalised"].startswith("offset:+3|")


def test_explicit_stop_codon_child_is_preserved_when_it_agrees() -> None:
    model, fasta, canonical, _ = _case(3, "+", extension_mode="extend_exon")
    model.child_records.append(
        _record("stop_codon", canonical[1] - 3, canonical[1], strand="+", name="m")
    )
    result = canonicalise_model_terminal_stop(model, fasta)
    stop_children = [
        (record.start0, record.end0)
        for record in result.child_records
        if record.feature_type == "stop_codon"
    ]
    assert stop_children == [(canonical[1] - 3, canonical[1])]
    assert result.cds_parts == [canonical]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_multiexon_terminal_part_is_changed_without_changing_splice_structure(
    strand: str,
) -> None:
    model, fasta, canonical_terminal, _ = _case(3, strand)
    if strand == "+":
        upstream = (50, 56)
        upstream_sequence = "ATGCCC"
    else:
        upstream = (244, 250)
        upstream_sequence = reverse_complement("ATGCCC")
    chars = list(fasta.sequence)
    chars[upstream[0]:upstream[1]] = upstream_sequence
    fasta = MemoryFasta("".join(chars))
    terminal = model.cds_parts[0]
    model.cds_parts = sorted([upstream, terminal])
    model.child_records.extend(
        [
            _record("exon", *upstream, strand=strand, name="m"),
            _record("CDS", *upstream, strand=strand, phase="0", name="m"),
        ]
    )
    model.__post_init__()
    result = canonicalise_model_terminal_stop(model, fasta)
    assert result.cds_parts == sorted([upstream, canonical_terminal])
    assert len(result.cds_introns) == 1
    assert result.cds_introns == model.cds_introns


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [3, 2, 1])
def test_extension_clips_only_overlapped_three_prime_utr(
    strand: str, offset: int
) -> None:
    model, fasta, canonical, _ = _case(offset, strand)
    original_utr = _add_adjacent_three_prime_utr(model, canonical)
    result = canonicalise_model_terminal_stop(model, fasta)
    expected = (
        (canonical[1], original_utr[1])
        if strand == "+"
        else (original_utr[0], canonical[0])
    )
    assert result.cds_parts == [canonical]
    assert _three_prime_utrs(result) == [expected]
    assert all(
        min(cds_end, utr_end) <= max(cds_start, utr_start)
        for cds_start, cds_end in result.cds_parts
        for utr_start, utr_end in _three_prime_utrs(result)
    )


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [3, 2, 1])
def test_extension_removes_fully_consumed_three_prime_utr(
    strand: str, offset: int
) -> None:
    model, fasta, canonical, _ = _case(offset, strand)
    _add_adjacent_three_prime_utr(model, canonical, completely_consumed=True)
    result = canonicalise_model_terminal_stop(model, fasta)
    assert result.cds_parts == [canonical]
    assert _three_prime_utrs(result) == []


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [-1, -2])
def test_trim_extends_directly_adjacent_three_prime_utr(
    strand: str, offset: int
) -> None:
    model, fasta, canonical, _ = _case(offset, strand)
    original_utr = _add_adjacent_three_prime_utr(model, canonical)
    result = canonicalise_model_terminal_stop(model, fasta)
    expected = (
        (canonical[1], original_utr[1])
        if strand == "+"
        else (original_utr[0], canonical[0])
    )
    assert result.cds_parts == [canonical]
    assert _three_prime_utrs(result) == [expected]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [-1, -2])
def test_trim_without_source_utr_does_not_synthesise_one(
    strand: str, offset: int
) -> None:
    model, fasta, canonical, original_exon = _case(offset, strand)
    result = canonicalise_model_terminal_stop(model, fasta)
    assert result.cds_parts == [canonical]
    assert _three_prime_utrs(result) == []
    assert [(r.start0, r.end0) for r in result.child_records if r.feature_type == "exon"] == [original_exon]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_nonadjacent_unaffected_three_prime_utr_and_five_prime_utr_are_unchanged(
    strand: str,
) -> None:
    model, fasta, canonical, _ = _case(3, strand)
    old = model.cds_parts[0]
    if strand == "+":
        three_utr = (canonical[1] + 3, canonical[1] + 9)
        five_utr = (old[0] - 9, old[0])
        exon = (five_utr[0], three_utr[1])
    else:
        three_utr = (canonical[0] - 9, canonical[0] - 3)
        five_utr = (old[1], old[1] + 9)
        exon = (three_utr[0], five_utr[1])
    _set_exon(model, exon)
    model.child_records.extend([
        _record("three_prime_UTR", *three_utr, strand=strand),
        _record("five_prime_UTR", *five_utr, strand=strand),
    ])
    result = canonicalise_model_terminal_stop(model, fasta)
    assert _three_prime_utrs(result) == [three_utr]
    assert [
        (record.start0, record.end0)
        for record in result.child_records
        if record.feature_type.lower() == "five_prime_utr"
    ] == [five_utr]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [3, -1])
def test_ordinary_final_gff_has_no_canonicalisation_induced_cds_utr_overlap(
    tmp_path: Path, strand: str, offset: int
) -> None:
    model, fasta, canonical, _ = _case(offset, strand)
    _add_adjacent_three_prime_utr(model, canonical)
    result = canonicalise_model_terminal_stop(model, fasta)
    result.sequence_integrity = sequence_integrity_metrics_from_row(
        _integrity(result, fasta)
    )
    transcript = FinalTranscript("MoSelecta_000001.1", result)
    annotation = FinalAnnotation((FinalGene("MoSelecta_000001", (transcript,)),))
    path = tmp_path / "final.gff3"
    write_final_annotation_gff3(
        annotation,
        path,
        FinalAnnotationOutputMetadata(
            alternative_splice_enabled=False,
            selector_evaluated=False,
            alternatives_emitted=False,
        ),
    )
    rows = [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    cds = [(int(row[3]) - 1, int(row[4])) for row in rows if row[2] == "CDS"]
    utrs = [
        (int(row[3]) - 1, int(row[4]))
        for row in rows
        if row[2].lower() == "three_prime_utr"
    ]
    assert cds == [canonical]
    assert utrs
    assert all(
        min(cds_end, utr_end) <= max(cds_start, utr_start)
        for cds_start, cds_end in cds
        for utr_start, utr_end in utrs
    )


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [3, 2, 1, -1, -2])
def test_source_utr_reconciliation_does_not_change_canonicalisation_decision(
    strand: str, offset: int
) -> None:
    plain, plain_fasta, canonical, _ = _case(offset, strand)
    with_utr, utr_fasta, _, _ = _case(offset, strand)
    _add_adjacent_three_prime_utr(with_utr, canonical)
    plain_result = canonicalise_model_terminal_stop(plain, plain_fasta)
    utr_result = canonicalise_model_terminal_stop(with_utr, utr_fasta)
    plain_audit = plain_result.terminal_stop_canonicalisation
    utr_audit = utr_result.terminal_stop_canonicalisation
    assert utr_result.cds_parts == plain_result.cds_parts == [canonical]
    assert (
        utr_audit.status,
        utr_audit.offset_bp,
        utr_audit.original_boundary,
        utr_audit.canonical_boundary,
        utr_audit.stop_codon,
        utr_audit.frame_confidence_basis,
    ) == (
        plain_audit.status,
        plain_audit.offset_bp,
        plain_audit.original_boundary,
        plain_audit.canonical_boundary,
        plain_audit.stop_codon,
        plain_audit.frame_confidence_basis,
    )
