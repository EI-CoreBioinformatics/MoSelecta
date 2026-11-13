from __future__ import annotations

from pathlib import Path

from moselecta.log_utils import RunTimer, current_rss_mb, peak_rss_mb
from moselecta.models import CodingModel, FeatureRecord, RelationshipEdge
from moselecta.protein_evidence import (
    ProteinEvidenceRawAlignment,
    build_protein_superloci,
    collapse_protein_alignments,
    release_protein_raw_alignments,
    write_protein_superloci_tsv,
)
from moselecta.transcript_evidence import (
    TranscriptEvidenceRawTranscript,
    build_transcript_superloci,
    collapse_transcript_evidence,
    release_transcript_raw_transcripts,
    write_transcript_superloci_tsv,
)


def test_high_count_core_and_evidence_records_use_slots() -> None:
    model = CodingModel(
        source_name="A",
        input_file=Path("a.gff3"),
        input_gene_id="g",
        input_transcript_id="t",
        output_gene_id="g",
        output_transcript_id="t",
        chrom="chr1",
        strand="+",
        cds_parts=[(0, 10), (20, 30)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )
    raw_protein = ProteinEvidenceRawAlignment(
        "p", "p1", 1, Path("p.gff3"), "p", "chr1", "+", ((0, 10),), 1
    )
    raw_tx = TranscriptEvidenceRawTranscript(
        "r", "t1", 1, Path("t.gtf"), "r", "chr1", "+", ((0, 10),), 1
    )
    edge = RelationshipEdge("sl", "sub", model, model, "test")
    record = FeatureRecord(
        "chr1", "src", "CDS", 1, 10, ".", "+", "0", {}, "", Path("a.gff3"), "A", 1, 1, ""
    )
    for obj in (model, raw_protein, raw_tx, edge, record):
        assert not hasattr(obj, "__dict__")


def test_release_protein_raw_records_preserves_superlocus_summary(tmp_path: Path) -> None:
    raw = [
        ProteinEvidenceRawAlignment(
            "p1", "sourceA", 1, Path("p.gff3"), "a1", "chr1", "+", ((0, 10),), 1
        ),
        ProteinEvidenceRawAlignment(
            "p2", "sourceB", 2, Path("p.gff3"), "a2", "chr1", "+", ((0, 10),), 1
        ),
    ]
    superloci = build_protein_superloci(raw)
    collapse_protein_alignments(superloci)
    assert release_protein_raw_alignments(raw, superloci) == 2
    assert raw == []
    assert superloci[0].raw_alignments == []
    out = tmp_path / "protein_superloci.tsv"
    write_protein_superloci_tsv(superloci, out)
    fields = out.read_text().splitlines()[1].split("\t")
    assert fields[6] == "2"
    assert fields[8] == "2"
    assert fields[9] == "sourceA:1,sourceB:1"


def test_release_transcript_raw_records_preserves_superlocus_summary(tmp_path: Path) -> None:
    raw = [
        TranscriptEvidenceRawTranscript(
            "t1", "sourceA", 1, Path("t.gtf"), "t1", "chr1", "+", ((0, 10), (20, 30)), 2
        ),
        TranscriptEvidenceRawTranscript(
            "t2", "sourceB", 2, Path("t.gtf"), "t2", "chr1", "+", ((0, 10), (20, 30)), 2
        ),
    ]
    superloci = build_transcript_superloci(raw)
    collapse_transcript_evidence(superloci)
    assert release_transcript_raw_transcripts(raw, superloci) == 2
    assert raw == []
    assert superloci[0].raw_transcripts == []
    out = tmp_path / "transcript_superloci.tsv"
    write_transcript_superloci_tsv(superloci, out)
    fields = out.read_text().splitlines()[1].split("\t")
    assert fields[6] == "2"
    assert fields[8] == "2"
    assert fields[9] == "sourceA:1,sourceB:1"


def test_runtime_timer_records_memory_columns() -> None:
    timer = RunTimer()
    timer.mark("stage")
    stage, duration, cumulative, current_rss, peak_rss = timer.rows[0]
    assert stage == "stage"
    assert duration >= 0
    assert cumulative >= 0
    assert current_rss > 0
    assert peak_rss > 0
    assert current_rss_mb() > 0
    assert peak_rss_mb() > 0
    assert current_rss_mb() <= peak_rss_mb() + 2.0
