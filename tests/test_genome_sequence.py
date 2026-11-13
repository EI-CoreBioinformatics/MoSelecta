from __future__ import annotations

import csv
import os
from pathlib import Path

import pytest

from moselecta.genome_sequence import (
    IndexedFasta,
    calculate_model_sequence_integrity,
    calculate_sequence_integrity_rows,
    sequence_integrity_log_lines,
    sequence_integrity_summary_rows,
)
from moselecta.models import CodingModel, FeatureRecord, Superlocus


def _cds_record(start0: int, end0: int, phase: str, strand: str = "+") -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1",
        source="test",
        feature_type="CDS",
        start=start0 + 1,
        end=end0,
        score=".",
        strand=strand,
        phase=phase,
        attrs={"Parent": "tx"},
        raw_attrs="Parent=tx",
        input_file=Path("test.gff3"),
        source_name="test",
        line_number=1,
        record_index=1,
        raw_line="",
    )


def _model(name: str, parts: list[tuple[int, int]], phases: list[str], strand: str = "+") -> CodingModel:
    records = [_cds_record(start, end, phase, strand) for (start, end), phase in zip(parts, phases)]
    return CodingModel(
        source_name="test",
        input_file=Path("test.gff3"),
        input_gene_id=f"g_{name}",
        input_transcript_id=name,
        output_gene_id=f"out_g_{name}",
        output_transcript_id=name,
        chrom="chr1",
        strand=strand,
        cds_parts=parts,
        gene_record=None,
        transcript_record=None,
        child_records=records,
    )


def _place(chars: list[str], start0: int, sequence: str) -> None:
    chars[start0 : start0 + len(sequence)] = list(sequence)


def test_fasta_index_is_built_and_random_access_works(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1 description\nACGTAC\nGTACGT\n>chr2\nTTTT\n")
    index = Path(str(fasta) + ".fai")
    assert not index.exists()
    with IndexedFasta(fasta) as reader:
        assert reader.index_built is True
        assert reader.fetch("chr1", 2, 10) == "GTACGTAC"
        assert reader.contig_length("chr2") == 4
    assert index.exists()


def test_sequence_metrics_start_stop_internal_stop_and_splice_classes() -> None:
    chars = list("C" * 130)
    _place(chars, 0, "ATGAAATAA")
    _place(chars, 20, "ATGAAA")
    _place(chars, 26, "TAA")
    _place(chars, 40, "ATGTAATAA")
    _place(chars, 60, "ATGAAA")
    _place(chars, 66, "GTCAG")
    _place(chars, 71, "CCCTAA")
    _place(chars, 80, "ATGAAA")
    _place(chars, 86, "AAAAA")
    _place(chars, 91, "CCCTAA")
    window = "".join(chars)

    included = calculate_model_sequence_integrity(
        _model("included", [(0, 9)], ["0"]),
        window=window,
        window_start0=0,
        contig_length=len(window),
    )
    assert included["sequence_start_codon_status"] == "present"
    assert included["sequence_stop_codon_status"] == "included_in_cds"
    assert included["sequence_integrity_penalty_score"] == 0.0

    adjacent = calculate_model_sequence_integrity(
        _model("adjacent", [(20, 26)], ["0"]),
        window=window,
        window_start0=0,
        contig_length=len(window),
    )
    assert adjacent["sequence_stop_codon_status"] == "adjacent_to_cds"
    assert adjacent["sequence_stop_codon"] == "TAA"
    assert adjacent["sequence_integrity_penalty_score"] == 0.0

    internal = calculate_model_sequence_integrity(
        _model("internal", [(40, 49)], ["0"]),
        window=window,
        window_start0=0,
        contig_length=len(window),
    )
    assert internal["sequence_internal_stop_count"] == 1
    assert internal["sequence_integrity_penalty_score"] == pytest.approx(1.0)

    canonical = calculate_model_sequence_integrity(
        _model("canonical", [(60, 66), (71, 77)], ["0", "0"]),
        window=window,
        window_start0=0,
        contig_length=len(window),
    )
    assert canonical["sequence_phase_status"] == "consistent"
    assert canonical["sequence_noncanonical_splice_intron_count"] == 0
    assert canonical["sequence_integrity_penalty_score"] == 0.0

    noncanonical = calculate_model_sequence_integrity(
        _model("noncanonical", [(80, 86), (91, 97)], ["0", "0"]),
        window=window,
        window_start0=0,
        contig_length=len(window),
    )
    assert noncanonical["sequence_noncanonical_splice_intron_count"] == 1
    assert noncanonical["sequence_noncanonical_splice_intron_fraction"] == 1.0
    assert "all_splices_noncanonical" in str(noncanonical["sequence_integrity_flags"])
    assert noncanonical["sequence_integrity_penalty_score"] == pytest.approx(1.0)


def test_minus_strand_sequence_and_splice_orientation() -> None:
    chars = list("C" * 60)
    # Transcript-oriented CDS is ATGAAA + CCCTAA. On the minus strand the high
    # genomic segment is the biological 5-prime segment.
    _place(chars, 10, "TTAGGG")
    _place(chars, 16, "CTAAC")  # reverse-oriented canonical GT-AG intron
    _place(chars, 21, "TTTCAT")
    window = "".join(chars)
    model = _model("minus", [(10, 16), (21, 27)], ["0", "0"], strand="-")
    row = calculate_model_sequence_integrity(
        model,
        window=window,
        window_start0=0,
        contig_length=len(window),
    )
    assert row["sequence_start_codon_status"] == "present"
    assert row["sequence_stop_codon_status"] == "included_in_cds"
    assert row["sequence_noncanonical_splice_intron_count"] == 0
    assert row["sequence_integrity_penalty_score"] == 0.0


def test_ambiguous_splice_is_noncanonical_and_boundary_missing_start_is_neutral() -> None:
    chars = list("C" * 40)
    _place(chars, 0, "CCCAAA")
    _place(chars, 6, "NNCNN")
    _place(chars, 11, "CCCTAA")
    row = calculate_model_sequence_integrity(
        _model("boundary", [(0, 6), (11, 17)], ["0", "0"]),
        window="".join(chars),
        window_start0=0,
        contig_length=len(chars),
    )
    assert row["sequence_start_codon_status"] == "not_evaluable_at_contig_boundary"
    assert row["sequence_noncanonical_splice_intron_count"] == 1
    assert row["sequence_has_ambiguous_bases"] == "true"
    assert row["sequence_integrity_penalty_score"] == pytest.approx(1.0)


def test_superlocus_processing_handles_missing_contig_non_strict(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">other\n" + "A" * 20 + "\n")
    model = _model("missing", [(0, 9)], ["0"])
    sl = Superlocus("sl1", "chr1", "+", 0, 9, [model])
    with IndexedFasta(fasta) as reader:
        rows = calculate_sequence_integrity_rows([sl], reader, strict=False)
    assert rows[0]["sequence_integrity_status"] == "not_evaluable_missing_contig"
    assert rows[0]["sequence_integrity_penalty_score"] == 0.0

def test_cli_genome_sequence_outputs_and_model_only_penalty(tmp_path: Path) -> None:
    import os
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[1]
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "CCCAAACCC" + "A" * 91 + "\n")
    gff = tmp_path / "model.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
    )
    prefix = tmp_path / "out"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    subprocess.run(
        [
            sys.executable,
            "-m",
            "moselecta",
            "--inputs",
            str(gff),
            "--source-name",
            "A",
            "--genome-fasta",
            str(fasta),
            "--out-prefix",
            str(prefix),
            "--output-profile",
            "standard",
            "--minimum-coding-length-bp",
            "1",
            "--write-output",
            "sequence_integrity_model_metrics_tsv,model_selection_scores_tsv",
            "--pruning-level",
            "none",
            "--progress-log-level",
            "none",
        ],
        cwd=root,
        env=env,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert Path(str(fasta) + ".fai").exists()
    assert Path(str(prefix) + ".sequence_integrity_summary.tsv").exists()
    assert Path(str(prefix) + ".sequence_integrity_model_metrics.tsv").exists()
    with Path(str(prefix) + ".sequence_integrity_model_metrics.tsv").open() as handle:
        detailed_header = next(csv.reader(handle, delimiter="\t"))
    assert "sequence_start_codon_status" in detailed_header
    assert "sequence_noncanonical_splice_details" in detailed_header
    with Path(str(prefix) + ".model_selection_scores.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert "sequence_integrity_status" in (reader.fieldnames or [])
        assert "sequence_integrity_penalty_score" in (reader.fieldnames or [])
        assert "sequence_integrity_flags" in (reader.fieldnames or [])
        assert "sequence_start_codon_status" not in (reader.fieldnames or [])
        assert "sequence_noncanonical_splice_details" not in (reader.fieldnames or [])
        row = next(reader)
    assert row["sequence_integrity_status"] == "evaluated_with_concern"
    assert float(row["sequence_integrity_penalty_score"]) == pytest.approx(0.75)
    assert row["sequence_integrity_penalty_in_selection"] == "true"
    assert float(row["model_structural_support_score_before_sequence_integrity"]) > float(
        row["model_structural_support_score"]
    )

@pytest.mark.parametrize("motif", ["GC-AG", "AT-AC"])
def test_recognised_minor_splice_motifs_are_canonical(motif: str) -> None:
    chars = list("C" * 50)
    _place(chars, 0, "ATGAAA")
    donor, acceptor = motif.split("-")
    intron = donor + "C" + acceptor
    _place(chars, 6, intron)
    _place(chars, 11, "CCCTAA")
    row = calculate_model_sequence_integrity(
        _model(motif, [(0, 6), (11, 17)], ["0", "0"]),
        window="".join(chars),
        window_start0=0,
        contig_length=len(chars),
    )
    assert row["sequence_noncanonical_splice_intron_count"] == 0
    assert row["sequence_integrity_penalty_score"] == 0.0


def test_phase_inconsistency_is_penalised_without_trimming_downstream_cds() -> None:
    chars = list("C" * 40)
    _place(chars, 0, "ATGAA")  # five coding bases; downstream expected phase is 1
    _place(chars, 5, "GTCAG")
    _place(chars, 10, "ACCCTAA")
    row = calculate_model_sequence_integrity(
        _model("phase", [(0, 5), (10, 17)], ["0", "0"]),
        window="".join(chars),
        window_start0=0,
        contig_length=len(chars),
    )
    assert row["sequence_phase_status"] == "inconsistent"
    assert row["sequence_phase_mismatch_count"] == 1
    assert float(row["sequence_integrity_penalty_score"]) >= 0.75


def test_diagnostic_only_sequence_mode_retains_metrics_without_score_penalty(tmp_path: Path) -> None:
    import os
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[1]
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "CCCAAACCC" + "A" * 91 + "\n")
    gff = tmp_path / "model.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
    )
    prefix = tmp_path / "diag"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    subprocess.run(
        [
            sys.executable,
            "-m",
            "moselecta",
            "--inputs",
            str(gff),
            "--source-name",
            "A",
            "--genome-fasta",
            str(fasta),
            "--no-sequence-integrity-scoring",
            "--out-prefix",
            str(prefix),
            "--output-profile",
            "standard",
            "--minimum-coding-length-bp",
            "1",
            "--write-output",
            "model_selection_scores_tsv",
            "--pruning-level",
            "none",
            "--progress-log-level",
            "none",
        ],
        cwd=root,
        env=env,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    with Path(str(prefix) + ".model_selection_scores.tsv").open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert float(row["sequence_integrity_penalty_score"]) > 0
    assert float(row["sequence_integrity_penalty_weight"]) == 0.0
    assert row["sequence_integrity_penalty_in_selection"] == "false"
    assert row["model_structural_support_score_before_sequence_integrity"] == row["model_structural_support_score"]


def test_attached_metrics_writer_matches_row_writer(tmp_path: Path) -> None:
    from moselecta.genome_sequence import (
        calculate_and_attach_sequence_integrity,
        sequence_integrity_row_from_model,
        write_sequence_integrity_model_metrics_from_models,
        write_sequence_integrity_model_metrics_tsv,
    )

    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nATGAAATAA" + "C" * 31 + "\n")
    model = _model("attached", [(0, 9)], ["0"])
    sl = Superlocus("sl1", "chr1", "+", 0, 9, [model])
    with IndexedFasta(fasta) as reader:
        summary, log_lines = calculate_and_attach_sequence_integrity([sl], reader)
    assert model.sequence_integrity is not None
    assert model.cds_phase_output is not None
    assert model.cds_phase_output.status == "preserved_numeric_source"
    assert model.cds_phase_output.first_phase == 0
    assert any(row["metric"] == "models_evaluated" and row["value"] == 1 for row in summary)
    assert "evaluated=1" in log_lines[0]

    attached_path = tmp_path / "attached.tsv"
    legacy_path = tmp_path / "legacy.tsv"
    write_sequence_integrity_model_metrics_from_models([model], attached_path)
    write_sequence_integrity_model_metrics_tsv([sequence_integrity_row_from_model(model)], legacy_path)
    assert attached_path.read_bytes() == legacy_path.read_bytes()
    row = next(csv.DictReader(attached_path.open(), delimiter="\t"))
    assert row["cds_phase_output_status"] == "preserved_numeric_source"
    assert row["cds_phase_output_first_phase"] == "0"
    assert row["cds_phase_output_reason"] == "numeric_source_phase"


def test_diagnostic_only_sequence_metrics_do_not_break_score_ties() -> None:
    from moselecta.constants import MODEL_SELECTION_DEFAULT_WEIGHTS
    from moselecta.scoring import apply_model_selection_scores

    def row(transcript_id: str, sequence_penalty: float) -> dict[str, object]:
        return {
            "sublocus_id": "sl1.sub1",
            "sublocus_source_count": 1,
            "sublocus_model_count": 2,
            "pairwise_other_source_count": 0,
            "source_name": "A",
            "output_transcript_id": transcript_id,
            "structural_support_score": 0.0,
            "consensus_similarity_score": 0.0,
            "pairwise_similarity_score": 0.0,
            "model_completeness_score": 0.0,
            "terminal_extension_score": 0.0,
            "exact_structure_support_score": 0.0,
            "unsupported_short_cds_exon_penalty_score": 0.0,
            "sequence_integrity_penalty_score": sequence_penalty,
            "intrinsic_penalty_score": 0.0,
            "partial_structure_penalty": 0.0,
            "internal_extra_fraction": 0.0,
            "terminal_extra_fraction": 0.0,
            "cds_length": 100,
        }

    rows = [row("a_penalized", 1.0), row("z_clean", 0.0)]
    weights = dict(MODEL_SELECTION_DEFAULT_WEIGHTS)
    weights["sequence_integrity_penalty"] = 0.0
    apply_model_selection_scores(
        rows,
        source_prior_scores={"A": 0.0},
        selection_weights=weights,
        tie_epsilon=1e-6,
    )
    by_id = {str(item["output_transcript_id"]): item for item in rows}
    assert by_id["a_penalized"]["model_structural_support_order_rank"] == 1
    assert by_id["a_penalized"]["sequence_integrity_penalty_in_selection"] == "false"



def test_sequence_integrity_summary_and_log_include_source_breakdown(tmp_path: Path) -> None:
    rows = [
        {
            "source_name": "Galba",
            "sequence_integrity_status": "evaluated_with_concern",
            "sequence_start_codon_status": "missing",
            "sequence_stop_codon_status": "present_in_cds",
            "sequence_internal_stop_count": 0,
            "sequence_phase_status": "consistent",
            "sequence_noncanonical_splice_intron_count": 0,
            "cds_intron_count": 1,
            "sequence_integrity_penalty_score": 0.30,
            "sequence_integrity_flags": "missing_start_codon",
        },
        {
            "source_name": "Helixer",
            "sequence_integrity_status": "evaluated_with_concern",
            "sequence_start_codon_status": "present",
            "sequence_stop_codon_status": "present_in_cds",
            "sequence_internal_stop_count": 0,
            "sequence_phase_status": "consistent",
            "sequence_noncanonical_splice_intron_count": 2,
            "cds_intron_count": 2,
            "sequence_integrity_penalty_score": 1.0,
            "sequence_integrity_flags": "noncanonical_splice,all_splices_noncanonical",
        },
        {
            "source_name": "Helixer",
            "sequence_integrity_status": "evaluated_no_concern",
            "sequence_start_codon_status": "present",
            "sequence_stop_codon_status": "adjacent_to_cds",
            "sequence_internal_stop_count": 0,
            "sequence_phase_status": "consistent",
            "sequence_noncanonical_splice_intron_count": 0,
            "cds_intron_count": 3,
            "sequence_integrity_penalty_score": 0.0,
            "sequence_integrity_flags": ".",
        },
    ]
    summary = sequence_integrity_summary_rows(
        rows,
        fasta_path=tmp_path / "genome.fa",
        index_path=tmp_path / "genome.fa.fai",
        index_built=False,
        genetic_code=1,
        penalty_active=True,
    )
    by_key = {
        (str(row["scope"]), str(row["source_name"]), str(row["metric"])): row
        for row in summary
    }
    helixer_concern = by_key[("source", "Helixer", "models_evaluated_with_concern")]
    assert helixer_concern["value"] == 1
    assert float(helixer_concern["fraction"]) == pytest.approx(0.5)
    helixer_introns = by_key[("source", "Helixer", "noncanonical_splice_intron_count")]
    assert helixer_introns["value"] == 2
    assert float(helixer_introns["fraction"]) == pytest.approx(2 / 5)
    assert by_key[("source", "Galba", "models_with_missing_start")]["value"] == 1
    assert by_key[("all_sources", ".", "models_total")]["value"] == 3

    lines = sequence_integrity_log_lines(rows)
    assert any(line.startswith("source\tmodels\tevaluated") for line in lines)
    assert any(line.startswith("Galba\t1\t1(100.00%)") for line in lines)
    assert any("Helixer\t2\t2(100.00%)\t1(50.00%)\t1(50.00%)" in line for line in lines)


def test_stale_fasta_index_is_rebuilt_automatically(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nAAAA\n")
    with IndexedFasta(fasta) as reader:
        assert reader.index_action == "created"
        assert reader.fetch("chr1", 0, 4) == "AAAA"
    index = Path(str(fasta) + ".fai")
    old_index_mtime = index.stat().st_mtime_ns

    fasta.write_text(">chr1\nCCCCCC\n>chr2\nTTTT\n")
    os.utime(fasta, ns=(old_index_mtime + 10_000_000, old_index_mtime + 10_000_000))
    with IndexedFasta(fasta) as reader:
        assert reader.index_action == "rebuilt"
        assert reader.index_rebuild_reason == "genome_fasta_newer_than_index"
        assert reader.fetch("chr1", 0, 6) == "CCCCCC"
        assert reader.fetch("chr2", 0, 4) == "TTTT"
    with IndexedFasta(fasta) as reader:
        assert reader.index_action == "reused"


def test_structurally_invalid_fasta_index_is_rebuilt_automatically(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGTACGT\n")
    index = Path(str(fasta) + ".fai")
    index.write_text("chr1\t999\t999999\t8\t9\n")
    # Keep the index timestamp newer so geometry validation, not mtime, triggers.
    fasta_mtime = fasta.stat().st_mtime_ns
    os.utime(index, ns=(fasta_mtime + 10_000_000, fasta_mtime + 10_000_000))
    with IndexedFasta(fasta) as reader:
        assert reader.index_action == "rebuilt"
        assert "out_of_range" in reader.index_rebuild_reason or "exceeds_fasta_size" in reader.index_rebuild_reason
        assert reader.fetch("chr1", 0, 8) == "ACGTACGT"
