from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from moselecta.intervals import (
    IntervalOverlapIndex,
    interval_list_overlap_bp,
    interval_list_overlap_bp_sorted,
    interval_prf,
    interval_prf_sorted,
)
from moselecta.models import CodingModel
from moselecta.protein_evidence import (
    ProteinEvidenceCollapsedAlignment,
    ProteinEvidenceSuperlocus,
    associated_evidence_for_model,
)
from moselecta.transcript_evidence import (
    TranscriptEvidenceCollapsedTranscript,
    TranscriptEvidenceSuperlocus,
    associated_transcripts_for_model,
    transcript_model_pair_metrics,
)


@dataclass
class Record:
    name: str
    chrom: str
    start0: int
    end0: int


def test_interval_overlap_index_matches_linear_scan_and_preserves_input_order() -> None:
    records = [
        Record("late_input", "chr1", 90, 160),
        Record("early_input", "chr1", 0, 100),
        Record("nested", "chr1", 25, 75),
        Record("other_chrom", "chr2", 40, 120),
        Record("distant", "chr1", 300, 400),
    ]
    index = IntervalOverlapIndex(
        records,
        key=lambda record: record.chrom,
        start=lambda record: record.start0,
        end=lambda record: record.end0,
    )
    expected = [
        record
        for record in records
        if record.chrom == "chr1" and record.end0 > 50 and record.start0 < 120
    ]
    assert index.query("chr1", 50, 120) == expected
    assert index.query("chr1", 160, 300) == []
    assert index.query("chrX", 0, 1000) == []


def test_fast_sorted_interval_helpers_match_general_helpers() -> None:
    query = ((0, 10), (20, 30), (40, 50))
    target = ((5, 25), (35, 45))
    assert interval_list_overlap_bp_sorted(query, target) == interval_list_overlap_bp(list(query), list(target))
    assert interval_prf_sorted(query, target) == interval_prf(list(query), list(target))


def make_model() -> CodingModel:
    return CodingModel(
        source_name="model",
        input_file=Path("model.gff3"),
        input_gene_id="g1",
        input_transcript_id="m1",
        output_gene_id="g1",
        output_transcript_id="m1",
        chrom="chr1",
        strand="+",
        cds_parts=[(100, 150), (250, 320)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def test_indexed_protein_association_matches_previous_linear_logic() -> None:
    model = make_model()
    aln_hit_1 = ProteinEvidenceCollapsedAlignment(
        "p1", "chr1", "+", ((100, 150), (250, 320)), "psl1", {"p": 1}, 1
    )
    aln_no_exonic_overlap = ProteinEvidenceCollapsedAlignment(
        "p2", "chr1", "+", ((160, 190),), "psl1", {"p": 1}, 1
    )
    aln_hit_2 = ProteinEvidenceCollapsedAlignment(
        "p3", "chr1", "+", ((300, 340),), "psl2", {"p": 1}, 1
    )
    opposite = ProteinEvidenceCollapsedAlignment(
        "p4", "chr1", "-", ((100, 150),), "psl3", {"p": 1}, 1
    )
    superloci = [
        ProteinEvidenceSuperlocus("psl1", "chr1", "+", 90, 330, collapsed_alignments=[aln_hit_1, aln_no_exonic_overlap]),
        ProteinEvidenceSuperlocus("psl2", "chr1", "+", 300, 350, collapsed_alignments=[aln_hit_2]),
        ProteinEvidenceSuperlocus("psl3", "chr1", "-", 90, 160, collapsed_alignments=[opposite]),
    ]
    index = IntervalOverlapIndex(
        superloci,
        key=lambda sl: (sl.chrom, sl.strand),
        start=lambda sl: sl.start0,
        end=lambda sl: sl.end0,
    )
    observed = associated_evidence_for_model(model, index)
    expected = []
    for sl in superloci:
        if (sl.chrom, sl.strand) != (model.chrom, model.strand):
            continue
        if sl.end0 <= model.cds_start0 or sl.start0 >= model.cds_end0:
            continue
        for aln in sl.collapsed_alignments:
            if interval_list_overlap_bp(model.cds_parts, list(aln.cds_parts)) > 0:
                expected.append(aln)
    assert [aln.collapsed_id for aln in observed] == [aln.collapsed_id for aln in expected]


def make_tx(tx_id: str, strand: str, parts: tuple[tuple[int, int], ...], super_id: str) -> TranscriptEvidenceCollapsedTranscript:
    return TranscriptEvidenceCollapsedTranscript(
        transcript_evidence_id=tx_id,
        chrom="chr1",
        strand=strand,
        exon_parts=parts,
        evidence_superlocus_id=super_id,
        collapse_key_type="exact_intron_chain" if len(parts) > 1 else "exact_single_exon_interval",
        source_member_counts={"t": 1},
        member_count=1,
        span_start0=parts[0][0],
        span_end0=parts[-1][1],
        left_terminal_outer_min0=parts[0][0],
        left_terminal_outer_max0=parts[0][0],
        right_terminal_outer_min0=parts[-1][1],
        right_terminal_outer_max0=parts[-1][1],
    )


def test_indexed_transcript_association_matches_previous_linear_logic() -> None:
    model = make_model()
    tx_hit = make_tx("t1", "+", ((90, 150), (250, 340)), "tsl1")
    tx_broad = make_tx("t2", "+", ((120, 145),), "tsl1")
    tx_opposite = make_tx("t3", "-", ((90, 150), (250, 340)), "tsl2")
    tx_unstranded = make_tx("t4", ".", ((100, 150), (250, 320)), "tsl3")
    superloci = [
        TranscriptEvidenceSuperlocus("tsl1", "chr1", "+", 90, 340, collapsed_transcripts=[tx_hit, tx_broad]),
        TranscriptEvidenceSuperlocus("tsl2", "chr1", "-", 90, 340, collapsed_transcripts=[tx_opposite]),
        TranscriptEvidenceSuperlocus("tsl3", "chr1", ".", 100, 320, collapsed_transcripts=[tx_unstranded]),
    ]
    index = IntervalOverlapIndex(
        superloci,
        key=lambda sl: sl.chrom,
        start=lambda sl: sl.start0,
        end=lambda sl: sl.end0,
    )
    observed = associated_transcripts_for_model(model, index)

    expected = []
    for sl in superloci:
        if sl.chrom != model.chrom or sl.end0 <= model.cds_start0 or sl.start0 >= model.cds_end0:
            continue
        for tx in sl.collapsed_transcripts:
            if tx.end0 <= model.cds_start0 or tx.start0 >= model.cds_end0:
                continue
            if tx.strand not in {model.strand, "."}:
                continue
            metrics = transcript_model_pair_metrics(model, tx)
            if metrics["association_class"] != ".":
                expected.append((tx, metrics))

    assert [tx.transcript_evidence_id for tx, _ in observed] == [tx.transcript_evidence_id for tx, _ in expected]
    assert [metrics for _, metrics in observed] == [metrics for _, metrics in expected]
