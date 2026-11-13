from __future__ import annotations

from pathlib import Path
import csv
import json

import pytest

from moselecta.cds_phase import (
    CDSFrameEvidence,
    RECOVERED_MISSING_SOURCE_PHASE,
    REFINED_NUMERIC_PHASE,
    calculate_cds_phase_chain,
    cds_phase_output_status,
)
from moselecta.genome_sequence import (
    calculate_model_sequence_integrity,
    reverse_complement,
)
from moselecta.models import CodingModel, FeatureRecord
from moselecta.phase_recovery import (
    choose_recovered_phase,
    phase_recovery_audit_row,
    recover_model_cds_phase,
    write_phase_recovery_actions_tsv,
    write_phase_recovery_audit_tsv,
    write_phase_recovery_summary_tsv,
)
from moselecta import pipeline as pipeline_module


class MemoryFasta:
    def __init__(self, sequence: str):
        self.sequence = sequence

    def contig_length(self, name: str) -> int | None:
        return len(self.sequence) if name == "chr1" else None

    def fetch(self, name: str, start0: int, end0: int) -> str:
        if name != "chr1" or start0 < 0 or end0 > len(self.sequence):
            raise ValueError("out of range")
        return self.sequence[start0:end0]


def _record(start0: int, end0: int, phase: str, strand: str) -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1", source="test", feature_type="CDS",
        start=start0 + 1, end=end0, score=".", strand=strand, phase=phase,
        attrs={"Parent": "tx"}, raw_attrs="Parent=tx",
        input_file=Path("test.gff3"), source_name="test",
        line_number=1, record_index=1, raw_line="",
    )


def _case(
    sequence: str,
    *,
    strand: str = "+",
    phase: str = ".",
    adjacent: str = "CCC",
    part_lengths: tuple[int, ...] | None = None,
) -> tuple[CodingModel, MemoryFasta]:
    lengths = part_lengths or (len(sequence),)
    assert sum(lengths) == len(sequence)
    chars = list("C" * 260)
    if strand == "+":
        biological_parts = []
        cursor = 40
        for length in lengths:
            biological_parts.append((cursor, cursor + length))
            cursor += length + 20
    else:
        biological_parts = []
        cursor = 210
        for length in lengths:
            biological_parts.append((cursor - length, cursor))
            cursor -= length + 20
    offset = 0
    for part, length in zip(biological_parts, lengths):
        chunk = sequence[offset : offset + length]
        genomic = chunk if strand == "+" else reverse_complement(chunk)
        chars[part[0] : part[1]] = genomic
        offset += length
    parts = sorted(biological_parts)
    if strand == "+":
        chars[parts[-1][1] : parts[-1][1] + 3] = adjacent
    else:
        chars[parts[0][0] - 3 : parts[0][0]] = reverse_complement(adjacent)
    model = CodingModel(
        source_name="test", input_file=Path("test.gff3"),
        input_gene_id="gene", input_transcript_id="tx",
        output_gene_id="out_gene", output_transcript_id="out_tx",
        chrom="chr1", strand=strand, cds_parts=parts,
        gene_record=None, transcript_record=None,
        child_records=[_record(start, end, phase, strand) for start, end in parts],
    )
    return model, MemoryFasta("".join(chars))


def _integrity(model: CodingModel, fasta: MemoryFasta) -> dict[str, object]:
    return calculate_model_sequence_integrity(
        model,
        window=fasta.sequence,
        window_start0=0,
        contig_length=len(fasta.sequence),
    )


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    ("phase", "sequence"),
    [(1, "CCCCAAATAA"), (2, "CCCCCAAATAA")],
)
def test_recovers_phase_from_included_stop_on_both_strands(
    strand: str, phase: int, sequence: str
) -> None:
    model, fasta = _case(sequence, strand=strand)
    before = _integrity(model, fasta)
    recovered = recover_model_cds_phase(model, fasta)
    after = _integrity(recovered, fasta)
    audit = recovered.cds_phase_recovery
    assert audit is not None
    assert audit.status == RECOVERED_MISSING_SOURCE_PHASE
    assert audit.recovered_first_phase == phase
    assert audit.stop_representation == "included"
    assert audit.original_effective_phase == 0
    assert after["sequence_start_codon_status"] == "missing"
    assert after["sequence_stop_codon_status"] == "included_in_cds"
    assert tuple(recovered.cds_parts) == tuple(model.cds_parts)
    assert before["sequence_stop_codon_status"] == "missing"


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    ("phase", "sequence"),
    [(1, "CCCCCCCCCC"), (2, "CCCCCCCCCCC")],
)
def test_recovers_phase_from_immediately_adjacent_stop(
    strand: str, phase: int, sequence: str
) -> None:
    model, fasta = _case(sequence, strand=strand, adjacent="TGA")
    recovered = recover_model_cds_phase(model, fasta)
    audit = recovered.cds_phase_recovery
    assert audit is not None
    assert audit.recovered_first_phase == phase
    assert audit.stop_representation == "adjacent"
    assert audit.stop_codon == "TGA"
    assert _integrity(recovered, fasta)["sequence_stop_codon_status"] == "adjacent_to_cds"


def test_summary_actions_and_exhaustive_audit_reconcile(tmp_path: Path) -> None:
    recovered_model, recovered_fasta = _case("CCCCAAATAA", phase=".")
    numeric_model, numeric_fasta = _case("ATGCCCTAA", phase="0")
    skipped_model, skipped_fasta = _case("CCCCCCCCCC", phase=".")
    recovered = recover_model_cds_phase(recovered_model, recovered_fasta)
    numeric = recover_model_cds_phase(numeric_model, numeric_fasta)
    skipped = recover_model_cds_phase(skipped_model, skipped_fasta)
    models = [recovered, numeric, skipped]
    audit_path = tmp_path / "audit.tsv"
    actions_path = tmp_path / "actions.tsv"
    summary_path = tmp_path / "summary.tsv"
    write_phase_recovery_audit_tsv(models, audit_path)
    write_phase_recovery_actions_tsv(models, actions_path)
    write_phase_recovery_summary_tsv(
        models, summary_path, execution_status="performed"
    )
    with audit_path.open(newline="") as handle:
        audit = list(csv.DictReader(handle, delimiter="\t"))
    with actions_path.open(newline="") as handle:
        actions = list(csv.DictReader(handle, delimiter="\t"))
    with summary_path.open(newline="") as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    assert actions == [
        row for row in audit
        if row["status"] == RECOVERED_MISSING_SOURCE_PHASE
    ]
    assert len(actions) == 1
    assert summary["total_models"] == "3"
    assert summary["recovered_models"] == "1"
    assert summary["numeric_source_phase_models"] == "1"
    assert summary["skipped_models"] == "1"
    assert json.loads(summary["recovered_phase_counts"]) == {"1": 1}
    assert json.loads(summary["skip_reason_counts"]) == {
        "no_unique_nonzero_stop_supported_frame": 1
    }


@pytest.mark.parametrize("strand", ["+", "-"])
def test_multi_exon_recovery_propagates_phase_across_mod3_junctions(
    strand: str,
) -> None:
    sequence = "CCCCCCAAAGTAA"  # phase 1: CCC CCA AAG TAA
    model, fasta = _case(
        sequence, strand=strand, part_lengths=(5, 8)
    )
    recovered = recover_model_cds_phase(model, fasta)
    expected = calculate_cds_phase_chain(recovered.cds_parts, strand, 1)
    observed = {
        (record.start0, record.end0): record.phase
        for record in recovered.child_records
        if record.feature_type == "CDS"
    }
    assert observed == expected
    assert _integrity(recovered, fasta)["sequence_phase_status"] == "consistent"


def test_recovery_removes_false_phase_zero_internal_stop_and_changes_normal_integrity_input() -> None:
    model, fasta = _case("TAACCCCTAA")
    before = _integrity(model, fasta)
    recovered = recover_model_cds_phase(model, fasta)
    after = _integrity(recovered, fasta)
    assert before["sequence_internal_stop_count"] == 1
    assert after["sequence_internal_stop_count"] == 0
    assert before["sequence_integrity_penalty_score"] > after["sequence_integrity_penalty_score"]
    assert after["sequence_start_codon_status"] == "missing"
    assert after["sequence_stop_codon_status"] == "included_in_cds"


@pytest.mark.parametrize("phase", ["0", "1", "2"])
def test_numeric_source_phase_is_never_overridden(phase: str) -> None:
    model, fasta = _case("CCCCAAATAA", phase=phase)
    recovered = recover_model_cds_phase(model, fasta)
    assert recovered is model
    assert recovered.child_records[0].phase == phase
    assert recovered.cds_phase_recovery is None
    assert phase_recovery_audit_row(recovered)["status"] == "not_eligible_numeric_source_phase"


@pytest.mark.parametrize(
    ("sequence", "adjacent", "reason"),
    [
        ("ATGCCCCTAA", "CCC", "phase_zero_positively_supported"),
        ("CATGCCCTAA", "CCC", "five_prime_leading_base_offset"),
        ("CCCCTAACCCTAA", "CCC", "recovered_frame_internal_stop"),
        ("CCCCCCCCCC", "CCC", "no_unique_nonzero_stop_supported_frame"),
        ("CATGCCCCCC", "CCC", "no_unique_nonzero_stop_supported_frame"),
        ("CCCCCCCCCC", "TAR", "ambiguous_terminal_or_frame_evidence"),
        ("CCCCNAAATAA", "CCC", "ambiguous_terminal_or_frame_evidence"),
    ],
)
def test_recovery_skip_classes(sequence: str, adjacent: str, reason: str) -> None:
    model, fasta = _case(sequence, adjacent=adjacent)
    unchanged = recover_model_cds_phase(model, fasta)
    assert unchanged.cds_phase_recovery is not None
    assert unchanged.cds_phase_recovery.status == "skipped"
    assert unchanged.cds_phase_recovery.reason == reason
    assert all(record.phase == "." for record in unchanged.child_records)


def test_p0_supported_even_with_internal_stop_is_not_recovered() -> None:
    model, fasta = _case("ATGTAACCTAA")
    unchanged = recover_model_cds_phase(model, fasta)
    assert unchanged.cds_phase_recovery is not None
    assert unchanged.cds_phase_recovery.reason == "phase_zero_positively_supported"


def test_leading_atg_label_is_only_used_after_offset_internal_stop_check() -> None:
    model, fasta = _case("CATGTAACCCTAA")
    unchanged = recover_model_cds_phase(model, fasta)
    assert unchanged.cds_phase_recovery is not None
    assert unchanged.cds_phase_recovery.reason == "recovered_frame_internal_stop"


def test_multiple_nonzero_candidates_are_rejected_by_policy_decision() -> None:
    def frame(phase: int, *, stop: bool = False) -> CDSFrameEvidence:
        return CDSFrameEvidence(
            phase=phase, definite_start=False, possible_start=False,
            definite_stop=stop, possible_stop=stop,
            stop_representation="included" if stop else ".",
            terminal_codon="TAA" if stop else ".",
            definite_internal_stop_count=0, possible_internal_stop_count=0,
        )
    phase, reason = choose_recovered_phase(
        {0: frame(0), 1: frame(1, stop=True), 2: frame(2, stop=True)}
    )
    assert phase is None
    assert reason == "no_unique_nonzero_stop_supported_frame"


def test_audit_and_final_phase_status_preserve_recovery_provenance() -> None:
    model, fasta = _case("CCCCAAATAA")
    recovered = recover_model_cds_phase(model, fasta)
    row = phase_recovery_audit_row(recovered)
    assert row["status"] == RECOVERED_MISSING_SOURCE_PHASE
    assert row["recovered_first_phase"] == 1
    assert row["original_stop_status"] == "missing"
    assert row["recovered_stop_status"] == "included_in_cds"
    assert cds_phase_output_status(recovered) == RECOVERED_MISSING_SOURCE_PHASE
    recovered.cds_terminus_refined = "upstream_start"
    assert cds_phase_output_status(recovered) == REFINED_NUMERIC_PHASE
    assert recovered.cds_phase_recovery.recovered_first_phase == 1


def test_pipeline_recovers_after_stop_canonicalisation_before_integrity_and_scoring(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "models.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t21\t30\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t21\t30\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t21\t30\t.\t+\t.\tParent=t1\n"
        "chr1\ttest\tCDS\t21\t30\t.\t+\t.\tParent=t1\n"
    )
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "C" * 20 + "CCCCAAATAA" + "C" * 20 + "\n")
    prefix = tmp_path / "out"
    assert pipeline_module.main([
        "--inputs", str(gff), "--source-name", "Test",
        "--genome-fasta", str(fasta), "--out-prefix", str(prefix),
        "--pruning-level", "none", "--no-alternative-splice",
        "--output-profile", "standard", "--progress-log-level", "none",
        "--minimum-coding-length-bp", "1",
        "--write-output", "sequence_integrity_model_metrics_tsv",
    ]) == 0

    with Path(str(prefix) + ".cds_phase_recovery.actions.tsv").open() as handle:
        audit = next(csv.DictReader(handle, delimiter="\t"))
    assert audit["status"] == RECOVERED_MISSING_SOURCE_PHASE
    assert audit["recovered_first_phase"] == "1"
    assert audit["stop_representation"] == "included"

    with Path(str(prefix) + ".sequence_integrity_model_metrics.tsv").open() as handle:
        integrity = next(csv.DictReader(handle, delimiter="\t"))
    assert integrity["sequence_start_codon_status"] == "missing"
    assert integrity["sequence_stop_codon_status"] == "included_in_cds"
    assert integrity["sequence_internal_stop_count"] == "0"
    assert integrity["cds_phase_output_status"] == RECOVERED_MISSING_SOURCE_PHASE

    final_gff = Path(str(prefix) + ".final_annotation.gff3").read_text()
    assert "cds_phase_recovered=1" in final_gff
    assert "cds_phase_output_status=recovered_missing_source_phase" in final_gff
    assert "\tCDS\t21\t30\t.\t+\t1\t" in final_gff
