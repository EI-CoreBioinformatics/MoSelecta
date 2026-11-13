from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

import pytest

from moselecta.extreme_short_cds import (
    EXCLUSION_CLASS,
    exclude_extreme_short_cds,
    write_extreme_short_cds_exclusions_tsv,
    write_extreme_short_cds_summary_tsv,
)
from moselecta.genome_sequence import reverse_complement
from moselecta.loci import build_superloci
from moselecta.models import CodingModel, FeatureRecord, ParsedInput
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


def _record(
    start0: int,
    end0: int,
    *,
    strand: str,
    phase: str,
    transcript_id: str,
) -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1", source="test", feature_type="CDS",
        start=start0 + 1, end=end0, score=".", strand=strand, phase=phase,
        attrs={"Parent": transcript_id}, raw_attrs=f"Parent={transcript_id}",
        input_file=Path("test.gff3"), source_name="test",
        line_number=1, record_index=1, raw_line="",
    )


def _case(
    sequence: str,
    *,
    strand: str = "+",
    phase: str = "0",
    adjacent: str = "CCC",
    part_lengths: tuple[int, ...] | None = None,
    source: str = "test",
    transcript_id: str = "tx",
    at_contig_edge: bool = False,
) -> tuple[ParsedInput, MemoryFasta]:
    lengths = part_lengths or (len(sequence),)
    assert sum(lengths) == len(sequence)
    genome = list("C" * 360)
    biological_parts: list[tuple[int, int]] = []
    if strand == "+":
        cursor = 0 if at_contig_edge else 40
        for length in lengths:
            biological_parts.append((cursor, cursor + length))
            cursor += length + 30
    else:
        cursor = len(genome) if at_contig_edge else 300
        for length in lengths:
            biological_parts.append((cursor - length, cursor))
            cursor -= length + 30
    offset = 0
    for part, length in zip(biological_parts, lengths):
        chunk = sequence[offset : offset + length]
        genome[part[0] : part[1]] = (
            chunk if strand == "+" else reverse_complement(chunk)
        )
        offset += length
    parts = sorted(biological_parts)
    if not at_contig_edge:
        if strand == "+":
            genome[parts[-1][1] : parts[-1][1] + 3] = adjacent
        else:
            genome[parts[0][0] - 3 : parts[0][0]] = reverse_complement(adjacent)
    model = CodingModel(
        source_name=source, input_file=Path(f"{source}.gff3"),
        input_gene_id=f"gene_{transcript_id}",
        input_transcript_id=transcript_id,
        output_gene_id=f"out_gene_{transcript_id}",
        output_transcript_id=f"out_{transcript_id}",
        chrom="chr1", strand=strand, cds_parts=parts,
        gene_record=None, transcript_record=None,
        child_records=[
            _record(
                start, end, strand=strand, phase=phase,
                transcript_id=transcript_id,
            )
            for start, end in parts
        ],
    )
    parsed = ParsedInput(
        input_file=Path(f"{source}.gff3"), source_name=source, records=[],
        coding_models=[model], excluded_models=[], used_record_indices=set(),
    )
    return parsed, MemoryFasta("".join(genome))


@pytest.mark.parametrize(
    ("length", "excluded"),
    [(29, True), (30, False), (31, False)],
)
def test_boundary_without_established_included_stop(
    length: int, excluded: bool
) -> None:
    parsed, _ = _case("A" * length)
    result = exclude_extreme_short_cds(
        [parsed], minimum_coding_length_bp=30
    )
    assert bool(result.excluded_models) is excluded
    assert bool(parsed.coding_models) is not excluded


@pytest.mark.parametrize(
    ("sequence", "phase", "expected_coding_length", "excluded"),
    [
        ("A" * 27 + "TAA", "0", 27, True),
        ("C" + "A" * 27 + "TAA", "1", 28, True),
        ("CC" + "A" * 27 + "TAA", "2", 29, True),
        ("A" * 30 + "TAA", "0", 30, False),
        ("A" * 33 + "TAA", "0", 33, False),
    ],
)
def test_included_canonical_stop_is_the_only_three_base_subtraction(
    sequence: str,
    phase: str,
    expected_coding_length: int,
    excluded: bool,
) -> None:
    parsed, fasta = _case(sequence, phase=phase)
    result = exclude_extreme_short_cds(
        [parsed], minimum_coding_length_bp=30, fasta=fasta
    )
    assert bool(result.excluded_models) is excluded
    if excluded:
        row = result.exclusions[0]
        assert row.terminal_stop_status == "included_in_cds"
        assert row.included_terminal_stop_subtracted is True
        assert row.terminal_stop_subtraction_bp == 3
        assert row.coding_length_bp == expected_coding_length


def test_adjacent_stop_and_ambiguous_terminal_stop_are_not_subtracted() -> None:
    adjacent, adjacent_fasta = _case("A" * 30, adjacent="TAA")
    adjacent_result = exclude_extreme_short_cds(
        [adjacent], minimum_coding_length_bp=30, fasta=adjacent_fasta
    )
    assert adjacent_result.excluded_models == 0

    ambiguous, ambiguous_fasta = _case("A" * 27 + "TAN")
    ambiguous_result = exclude_extreme_short_cds(
        [ambiguous], minimum_coding_length_bp=30, fasta=ambiguous_fasta
    )
    assert ambiguous_result.excluded_models == 0


@pytest.mark.parametrize(
    ("sequence", "expected_coding_length"),
    [("TAA", 0), ("A" * 12, 12), ("A" * 18, 18)],
)
def test_pathologically_small_models_are_excluded(
    sequence: str, expected_coding_length: int
) -> None:
    parsed, fasta = _case(sequence)
    result = exclude_extreme_short_cds(
        [parsed], minimum_coding_length_bp=30, fasta=fasta
    )
    assert result.excluded_models == 1
    assert result.exclusions[0].coding_length_bp == expected_coding_length


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("part_lengths", [(29,), (14, 15), (7, 9, 13)])
def test_spliced_cds_length_not_genomic_span_controls_exclusion(
    strand: str, part_lengths: tuple[int, ...]
) -> None:
    parsed, _ = _case(
        "A" * 29, strand=strand, part_lengths=part_lengths
    )
    result = exclude_extreme_short_cds(
        [parsed], minimum_coding_length_bp=30
    )
    assert result.excluded_models == 1
    assert result.exclusions[0].total_cds_length_bp == 29
    assert result.exclusions[0].cds_exon_count == len(part_lengths)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_subthreshold_contig_edge_model_is_not_rescued(strand: str) -> None:
    parsed, fasta = _case("A" * 18, strand=strand, at_contig_edge=True)
    result = exclude_extreme_short_cds(
        [parsed], minimum_coding_length_bp=30, fasta=fasta
    )
    assert result.excluded_models == 1


def test_source_identity_and_duplicate_support_do_not_rescue() -> None:
    first, _ = _case("A" * 18, source="Original", transcript_id="a")
    second, _ = _case("A" * 18, source="Predictor", transcript_id="b")
    second.coding_models[0] = replace(
        second.coding_models[0], cds_parts=list(first.coding_models[0].cds_parts)
    )
    result = exclude_extreme_short_cds(
        [first, second], minimum_coding_length_bp=30
    )
    assert result.excluded_models == 2
    assert not first.coding_models and not second.coding_models
    assert result.excluded_source_counts == {"Original": 1, "Predictor": 1}


def test_summary_and_ledger_reconcile_and_are_deterministic(
    tmp_path: Path,
) -> None:
    excluded, _ = _case("A" * 18, source="B", transcript_id="bad")
    retained, _ = _case("A" * 30, source="A", transcript_id="good")
    result = exclude_extreme_short_cds(
        [excluded, retained], minimum_coding_length_bp=30
    )
    summary_path = tmp_path / "summary.tsv"
    ledger_path = tmp_path / "ledger.tsv"
    write_extreme_short_cds_summary_tsv(result, summary_path)
    write_extreme_short_cds_exclusions_tsv(result, ledger_path)
    with summary_path.open(newline="") as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    with ledger_path.open(newline="") as handle:
        ledger = list(csv.DictReader(handle, delimiter="\t"))
    assert int(summary["candidates_evaluated"]) == 2
    assert int(summary["excluded_models"]) == len(ledger) == 1
    assert int(summary["retained_models"]) == 1
    assert json.loads(summary["excluded_source_counts"]) == {"B": 1}
    assert ledger[0]["exclusion_class"] == EXCLUSION_CLASS
    assert ledger[0]["coding_length_bp"] == "18"


def test_filter_removes_a_tiny_locus_connector_before_clustering() -> None:
    left, _ = _case("A" * 30, transcript_id="left")
    right, _ = _case("A" * 30, transcript_id="right")
    connector, _ = _case("A" * 23, transcript_id="connector")
    left_model = replace(left.coding_models[0], cds_parts=[(0, 30)])
    right_model = replace(right.coding_models[0], cds_parts=[(49, 79)])
    connector_model = replace(connector.coding_models[0], cds_parts=[(28, 51)])
    parsed = ParsedInput(
        input_file=Path("test.gff3"), source_name="test", records=[],
        coding_models=[left_model, connector_model, right_model],
        excluded_models=[], used_record_indices=set(),
    )
    result = exclude_extreme_short_cds(
        [parsed], minimum_coding_length_bp=30
    )
    superloci = build_superloci(parsed.coding_models)
    assert result.excluded_models == 1
    assert [model.output_transcript_id for model in parsed.coding_models] == [
        "out_left", "out_right"
    ]
    assert len(superloci) == 2


def test_pipeline_reports_default_on_no_genome_and_excludes_before_outputs(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "models.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t60\t.\t+\t.\tID=g_good\n"
        "chr1\ttest\tmRNA\t1\t60\t.\t+\t.\tID=t_good;Parent=g_good\n"
        "chr1\ttest\texon\t1\t60\t.\t+\t.\tParent=t_good\n"
        "chr1\ttest\tCDS\t1\t60\t.\t+\t0\tParent=t_good\n"
        "chr1\ttest\tgene\t100\t117\t.\t+\t.\tID=g_short\n"
        "chr1\ttest\tmRNA\t100\t117\t.\t+\t.\tID=t_short;Parent=g_short\n"
        "chr1\ttest\texon\t100\t117\t.\t+\t.\tParent=t_short\n"
        "chr1\ttest\tCDS\t100\t117\t.\t+\t0\tParent=t_short\n"
    )
    prefix = tmp_path / "out"
    assert pipeline_module.main([
        "--inputs", str(gff), "--source-name", "Test",
        "--out-prefix", str(prefix), "--pruning-level", "none",
        "--output-profile", "full", "--progress-log-level", "none",
    ]) == 0
    with Path(str(prefix) + ".extreme_short_cds_exclusions.tsv").open(
        newline=""
    ) as handle:
        ledger = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["original_transcript_id"] for row in ledger] == ["t_short"]
    for suffix in (
        "superlocus_members.tsv", "model_structural_support.tsv",
        "model_selection_scores.tsv", "selected_models.tsv",
        "backfill_candidates.tsv", "final_selected_models.tsv",
        "final_selected_models.pruned.tsv", "alternative_splice_candidates.tsv",
    ):
        path = Path(str(prefix) + "." + suffix)
        if path.exists():
            assert "t_short" not in path.read_text(), suffix
    assert "t_short" not in Path(str(prefix) + ".final_annotation.gff3").read_text()
    configuration = json.loads(
        Path(str(prefix) + ".run_configuration.json").read_text()
    )
    assert configuration["extreme_short_cds_exclusion"] == {
        "default_enabled": True,
        "effective_enabled": True,
        "evidence_rescue": False,
        "included_terminal_stop_detection": "unavailable_without_genome_fasta",
        "included_terminal_stop_subtraction_bp": 3,
        "minimum_coding_length_bp": 30,
        "no_genome_policy": "subtract_nothing_unless_included_stop_established",
    }


def test_h3_precedes_short_filter_and_each_model_has_one_early_exclusion(
    tmp_path: Path,
) -> None:
    sequences = {
        "h3": (1, "ATGTAACCCTAA"),
        "short": (41, "ATGAAAAAAAAAAAATAA"),
        "eligible": (81, "ATG" + "AAA" * 9),
    }
    lines = ["##gff-version 3"]
    genome = list("C" * 180)
    for name, (start, sequence) in sequences.items():
        end = start + len(sequence) - 1
        genome[start - 1 : end] = sequence
        lines.extend([
            f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID=g_{name}",
            f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID=t_{name};Parent=g_{name}",
            f"chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\tParent=t_{name}",
            f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tParent=t_{name}",
        ])
    gff = tmp_path / "early_filters.gff3"
    gff.write_text("\n".join(lines) + "\n")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "".join(genome) + "\n")
    prefix = tmp_path / "early"
    assert pipeline_module.main([
        "--inputs", str(gff), "--source-name", "Test",
        "--genome-fasta", str(fasta), "--out-prefix", str(prefix),
        "--pruning-level", "none", "--output-profile", "standard",
        "--progress-log-level", "none",
    ]) == 0
    with Path(str(prefix) + ".definite_internal_stop_exclusions.tsv").open(
        newline=""
    ) as handle:
        h3_rows = list(csv.DictReader(handle, delimiter="\t"))
    with Path(str(prefix) + ".extreme_short_cds_exclusions.tsv").open(
        newline=""
    ) as handle:
        short_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["original_transcript_id"] for row in h3_rows] == ["t_h3"]
    assert [row["original_transcript_id"] for row in short_rows] == ["t_short"]
    final_gff = Path(str(prefix) + ".final_annotation.gff3").read_text()
    assert "t_h3" not in final_gff
    assert "t_short" not in final_gff
    assert "t_eligible" in final_gff
    configuration = json.loads(
        Path(str(prefix) + ".run_configuration.json").read_text()
    )
    assert configuration["extreme_short_cds_exclusion"][
        "included_terminal_stop_detection"
    ] == "genome_sequence_available"
