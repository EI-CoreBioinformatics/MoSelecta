from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

import pytest

from moselecta.cds_phase import RECOVERED_MISSING_SOURCE_PHASE, analyse_cds_frame
from moselecta.definite_internal_stop import (
    INVALID_CDS_GEOMETRY,
    RESOLVED_FRAME_INTERNAL_STOP,
    UNRESOLVED_ALL_FRAMES_INTERNAL_STOP,
    empty_definite_internal_stop_result,
    recover_phases_and_exclude_definite_internal_stops,
    write_definite_internal_stop_exclusions_tsv,
    write_definite_internal_stop_summary_tsv,
)
from moselecta.genome_sequence import reverse_complement
from moselecta.loci import build_superloci
from moselecta.models import (
    CDSPhaseRecovery,
    CodingModel,
    FeatureRecord,
    FivePrimeBoundaryCanonicalisation,
    ParsedInput,
)
from moselecta.phase_recovery import model_sequence_context
from moselecta import pipeline as pipeline_module


class CountingMemoryFasta:
    def __init__(self, sequence: str):
        self.sequence = sequence
        self.fetch_count = 0

    def contig_length(self, name: str) -> int | None:
        return len(self.sequence) if name == "chr1" else None

    def fetch(self, name: str, start0: int, end0: int) -> str:
        self.fetch_count += 1
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
) -> tuple[ParsedInput, CountingMemoryFasta]:
    lengths = part_lengths or (len(sequence),)
    assert sum(lengths) == len(sequence)
    genome = list("C" * 320)
    if strand == "+":
        biological_parts = []
        cursor = 40
        for length in lengths:
            biological_parts.append((cursor, cursor + length))
            cursor += length + 20
    else:
        biological_parts = []
        cursor = 270
        for length in lengths:
            biological_parts.append((cursor - length, cursor))
            cursor -= length + 20
    offset = 0
    for part, length in zip(biological_parts, lengths):
        chunk = sequence[offset : offset + length]
        genome[part[0] : part[1]] = (
            chunk if strand == "+" else reverse_complement(chunk)
        )
        offset += length
    parts = sorted(biological_parts)
    if strand == "+":
        genome[parts[-1][1] : parts[-1][1] + 3] = adjacent
    else:
        genome[parts[0][0] - 3 : parts[0][0]] = reverse_complement(adjacent)
    model = CodingModel(
        source_name="test", input_file=Path("test.gff3"),
        input_gene_id="gene", input_transcript_id="tx",
        output_gene_id="out_gene", output_transcript_id="out_tx",
        chrom="chr1", strand=strand, cds_parts=parts,
        gene_record=None, transcript_record=None,
        child_records=[_record(start, end, phase, strand) for start, end in parts],
    )
    parsed = ParsedInput(
        input_file=Path("test.gff3"), source_name="test", records=[],
        coding_models=[model], excluded_models=[], used_record_indices=set(),
    )
    return parsed, CountingMemoryFasta("".join(genome))


def _overlapping_case(
    *,
    strand: str,
) -> tuple[ParsedInput, CountingMemoryFasta, str]:
    """Build a 1-bp CDS overlap whose literal concatenation stops in all frames."""
    literal_sequence = "TAGTACTTGATTGAT"
    first_chunk = literal_sequence[:7]
    second_chunk = literal_sequence[7:]
    assert first_chunk[-1] == second_chunk[0] == "T"
    genome = list("C" * 160)
    if strand == "+":
        parts = [(40, 47), (46, 54)]
        genome[40:47] = first_chunk
        genome[46:54] = second_chunk
    else:
        parts = [(40, 48), (47, 54)]
        genome[40:48] = reverse_complement(second_chunk)
        genome[47:54] = reverse_complement(first_chunk)
    if strand == "+":
        genome[54:57] = "CCC"
    else:
        genome[37:40] = reverse_complement("CCC")
    model = CodingModel(
        source_name="test", input_file=Path("overlap.gff3"),
        input_gene_id="gene_overlap", input_transcript_id="tx_overlap",
        output_gene_id="out_gene_overlap",
        output_transcript_id="out_tx_overlap",
        chrom="chr1", strand=strand, cds_parts=parts,
        gene_record=None, transcript_record=None,
        child_records=[_record(start, end, ".", strand) for start, end in parts],
    )
    parsed = ParsedInput(
        input_file=Path("overlap.gff3"), source_name="test", records=[],
        coding_models=[model], excluded_models=[], used_record_indices=set(),
    )
    return parsed, CountingMemoryFasta("".join(genome)), literal_sequence


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    ("phase", "sequence"),
    [
        ("0", "ATGTAACCCTAA"),
        ("1", "CTAACCCTAA"),
        ("2", "CCTAACCCTAA"),
    ],
)
def test_resolved_active_frame_internal_stop_is_excluded(
    strand: str, phase: str, sequence: str
) -> None:
    parsed, fasta = _case(sequence, strand=strand, phase=phase)
    _, audit_models, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert len(audit_models) == 1
    assert parsed.coding_models == []
    assert result.excluded_models == 1
    assert result.exclusions[0].exclusion_class == RESOLVED_FRAME_INTERNAL_STOP
    assert result.exclusions[0].active_first_cds_phase == int(phase)
    assert result.excluded_resolved_frame_internal_stop == 1


@pytest.mark.parametrize("strand", ["+", "-"])
def test_unresolved_all_three_frames_internal_stop_is_excluded_multiexon(
    strand: str,
) -> None:
    parsed, fasta = _case(
        "TAGTACTTGATTGAT", strand=strand, phase=".", part_lengths=(7, 8)
    )
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert parsed.coding_models == []
    exclusion = result.exclusions[0]
    assert exclusion.exclusion_class == UNRESOLVED_ALL_FRAMES_INTERNAL_STOP
    assert (
        exclusion.frame0_internal_stop_count,
        exclusion.frame1_internal_stop_count,
        exclusion.frame2_internal_stop_count,
    ) == (1, 1, 1)
    assert result.excluded_invalid_cds_geometry == 0


@pytest.mark.parametrize("strand", ["+", "-"])
def test_one_base_overlapping_cds_is_explicit_invalid_geometry(
    strand: str,
) -> None:
    parsed, fasta, _ = _overlapping_case(strand=strand)
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert parsed.coding_models == []
    assert result.excluded_models == 1
    assert result.excluded_invalid_cds_geometry == 1
    assert result.excluded_unresolved_all_frames_internal_stop == 0
    exclusion = result.exclusions[0]
    assert exclusion.exclusion_class == INVALID_CDS_GEOMETRY
    assert exclusion.sequence_evaluability == INVALID_CDS_GEOMETRY
    assert (
        exclusion.frame0_internal_stop_count,
        exclusion.frame1_internal_stop_count,
        exclusion.frame2_internal_stop_count,
    ) == (0, 0, 0)


def test_overlap_literal_concatenation_would_have_stops_in_all_frames() -> None:
    parsed, fasta, expected_literal_sequence = _overlapping_case(strand="+")
    sequence, adjacent = model_sequence_context(parsed.coding_models[0], fasta)
    assert sequence == expected_literal_sequence
    assert all(
        analyse_cds_frame(
            sequence, adjacent, phase, cds_is_acgt=True
        ).definite_internal_stop_count > 0
        for phase in (0, 1, 2)
    )
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert result.exclusions[0].exclusion_class == INVALID_CDS_GEOMETRY
    assert result.excluded_unresolved_all_frames_internal_stop == 0


def test_valid_nonoverlapping_multiexon_clean_frame_remains_eligible() -> None:
    parsed, fasta = _case(
        "CATGCCCTAA", phase=".", part_lengths=(4, 6)
    )
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert len(parsed.coding_models) == 1
    assert result.excluded_models == 0
    assert result.excluded_invalid_cds_geometry == 0


def test_unresolved_model_with_one_clean_frame_is_retained() -> None:
    parsed, fasta = _case("ATGTAACCCTAA", phase=".")
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert len(parsed.coding_models) == 1
    assert result.excluded_models == 0
    assert result.retained_unresolved_with_at_least_one_clean_frame == 1


def test_ambiguous_stop_like_codon_does_not_count_as_literal() -> None:
    parsed, fasta = _case("TARTACTTGATTGAT", phase=".")
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert len(parsed.coding_models) == 1
    assert result.excluded_models == 0
    assert result.retained_unresolved_with_at_least_one_clean_frame == 1


def test_resolved_ambiguous_possible_stop_is_retained() -> None:
    parsed, fasta = _case("ATGTARCCCTAA", phase="0")
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert len(parsed.coding_models) == 1
    assert result.excluded_models == 0


@pytest.mark.parametrize(
    ("sequence", "expected_counts", "excluded"),
    [
        ("TAACTAACCCCCCCC", (1, 1, 0), False),
        ("TAACTAACTAACCCC", (1, 1, 1), True),
        ("TAACTAACTAACTAACTAACTAACCCC", (2, 2, 2), True),
    ],
)
def test_unresolved_two_frame_and_all_frame_stop_populations(
    sequence: str, expected_counts: tuple[int, int, int], excluded: bool
) -> None:
    parsed, fasta = _case(sequence, phase=".")
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert (len(parsed.coding_models) == 0) is excluded
    if excluded:
        row = result.exclusions[0]
        assert (
            row.frame0_internal_stop_count,
            row.frame1_internal_stop_count,
            row.frame2_internal_stop_count,
        ) == expected_counts
    else:
        assert result.retained_unresolved_with_at_least_one_clean_frame == 1


def test_recovered_and_class_b_resolved_provenance_use_active_frame_rule() -> None:
    recovered_parsed, recovered_fasta = _case("CTAACCCTAA", phase="1")
    recovered_parsed.coding_models[0] = replace(
        recovered_parsed.coding_models[0],
        cds_phase_recovery=CDSPhaseRecovery(
            status=RECOVERED_MISSING_SOURCE_PHASE,
            recovered_first_phase=1,
        ),
    )
    _, _, recovered_result = recover_phases_and_exclude_definite_internal_stops(
        [recovered_parsed], recovered_fasta
    )
    assert recovered_parsed.coding_models == []
    assert recovered_result.exclusions[0].phase_provenance == RECOVERED_MISSING_SOURCE_PHASE

    boundary_parsed, boundary_fasta = _case("ATGTAACCCTAA", phase="0")
    boundary_parsed.coding_models[0].five_prime_boundary_canonicalisation = (
        FivePrimeBoundaryCanonicalisation(status="canonicalised")
    )
    _, _, boundary_result = recover_phases_and_exclude_definite_internal_stops(
        [boundary_parsed], boundary_fasta
    )
    assert boundary_parsed.coding_models == []
    assert (
        boundary_result.exclusions[0].phase_provenance
        == "canonicalised_five_prime_boundary"
    )


@pytest.mark.parametrize(
    ("sequence", "adjacent"),
    [
        ("ATGCCCTAA", "CCC"),  # included terminal stop
        ("ATGCCC", "TAA"),  # immediately adjacent terminal stop
        ("ATGCCCT", "CCC"),  # trailing partial codon without a stop
    ],
)
def test_terminal_adjacent_and_nonstop_partial_codons_are_not_internal(
    sequence: str, adjacent: str
) -> None:
    parsed, fasta = _case(sequence, phase="0", adjacent=adjacent)
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert len(parsed.coding_models) == 1
    assert result.excluded_models == 0


def test_stop_before_trailing_partial_is_internal_and_excluded() -> None:
    parsed, fasta = _case("ATGTAACC", phase="0")
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert parsed.coding_models == []
    assert result.excluded_models == 1


def test_phase_recovery_bridge_reuses_one_fasta_context_and_retains_clean_frame() -> None:
    parsed, fasta = _case("CCCCAAATAA", phase=".", part_lengths=(4, 6))
    phase_counts, audit_models, result = (
        recover_phases_and_exclude_definite_internal_stops([parsed], fasta)
    )
    assert len(parsed.coding_models) == 1
    assert audit_models[0].cds_phase_recovery.status == RECOVERED_MISSING_SOURCE_PHASE
    assert phase_counts == {"recovered_phase_1": 1}
    assert result.resolved_frame_models == 1
    assert result.excluded_models == 0
    # Two CDS chunks plus one adjacent codon: phase recovery and H3 share them.
    assert fasta.fetch_count == 3


def test_h3_bridge_is_removed_before_superlocus_construction() -> None:
    genome = list("C" * 140)
    genome[20:29] = "ATGCCCTAA"
    genome[50:77] = "TAACTAACTAACTAACTAACTAACCCC"
    genome[100:109] = "ATGCCCTAA"
    fasta = CountingMemoryFasta("".join(genome))

    def model(name: str, parts: list[tuple[int, int]]) -> CodingModel:
        return CodingModel(
            source_name="test", input_file=Path("bridge.gff3"),
            input_gene_id=f"g_{name}", input_transcript_id=name,
            output_gene_id=f"og_{name}", output_transcript_id=f"ot_{name}",
            chrom="chr1", strand="+", cds_parts=parts,
            gene_record=None, transcript_record=None,
            child_records=[_record(start, end, ".", "+") for start, end in parts],
        )

    left = model("left", [(20, 29)])
    right = model("right", [(100, 109)])
    bridge = model("bridge", [(20, 109)])
    assert len(build_superloci([left, bridge, right])) == 1
    parsed = ParsedInput(
        input_file=Path("bridge.gff3"), source_name="test", records=[],
        coding_models=[left, bridge, right], excluded_models=[],
        used_record_indices=set(),
    )
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert [item.original_transcript_id for item in result.exclusions] == ["bridge"]
    assert len(build_superloci(parsed.coding_models)) == 2


def test_class_a_leading_base_case_is_not_broadly_excluded() -> None:
    parsed, fasta = _case("CATGCCCTAA", phase=".")
    _, audit_models, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], fasta
    )
    assert len(parsed.coding_models) == 1
    assert audit_models[0].cds_phase_recovery.reason == "five_prime_leading_base_offset"
    assert result.excluded_models == 0


def test_summary_and_ledger_reconcile_with_sorted_source_json(tmp_path: Path) -> None:
    parsed_a, fasta_a = _case("ATGTAACCCTAA", phase="0")
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed_a], fasta_a
    )
    summary_path = tmp_path / "summary.tsv"
    ledger_path = tmp_path / "ledger.tsv"
    write_definite_internal_stop_summary_tsv(result, summary_path)
    write_definite_internal_stop_exclusions_tsv(result, ledger_path)
    with summary_path.open(newline="") as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    with ledger_path.open(newline="") as handle:
        ledger = list(csv.DictReader(handle, delimiter="\t"))
    assert summary["total_models"] == "1"
    assert summary["evaluable_models"] == "1"
    assert summary["excluded_models"] == str(len(ledger)) == "1"
    assert summary["excluded_invalid_cds_geometry"] == "0"
    assert json.loads(summary["excluded_source_counts"]) == {"test": 1}
    assert ledger[0]["exclusion_class"] == RESOLVED_FRAME_INTERNAL_STOP


def test_unavailable_sequence_is_retained_as_not_evaluable() -> None:
    parsed, _ = _case("ATGTAACCCTAA", phase="0")
    unavailable = CountingMemoryFasta("C" * 10)
    _, _, result = recover_phases_and_exclude_definite_internal_stops(
        [parsed], unavailable
    )
    assert len(parsed.coding_models) == 1
    assert result.evaluable_models == 0
    assert result.not_evaluable_models == 1
    assert result.excluded_models == 0


def test_no_genome_result_is_explicitly_inactive() -> None:
    result = empty_definite_internal_stop_result(7)
    assert result.execution_status == "prerequisite_unavailable:no_genome_fasta"
    assert result.total_models == 7
    assert result.not_evaluable_models == 7
    assert result.excluded_models == 0


def test_pipeline_hard_exclusion_precedes_scoring_and_ignores_evidence_rescue_switch(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "models.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t21\t32\t.\t+\t.\tID=g_bad\n"
        "chr1\ttest\tmRNA\t21\t32\t.\t+\t.\tID=t_bad;Parent=g_bad\n"
        "chr1\ttest\texon\t21\t32\t.\t+\t.\tParent=t_bad\n"
        "chr1\ttest\tCDS\t21\t32\t.\t+\t0\tParent=t_bad\n"
        "chr1\ttest\tgene\t21\t32\t.\t+\t.\tID=g_good\n"
        "chr1\ttest\tmRNA\t21\t32\t.\t+\t.\tID=t_good;Parent=g_good\n"
        "chr1\ttest\texon\t21\t32\t.\t+\t.\tParent=t_good\n"
        "chr1\ttest\tCDS\t21\t32\t.\t+\t1\tParent=t_good\n"
    )
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "C" * 20 + "ATGTAACCCTAA" + "C" * 20 + "\n")
    protein = tmp_path / "protein.gff3"
    protein.write_text(
        "##gff-version 3\n"
        "chr1\tprotein\tCDS\t21\t32\t.\t+\t0\tID=p1.cds;Parent=p1\n"
    )
    prefix = tmp_path / "out"
    assert pipeline_module.main([
        "--inputs", str(gff), "--source-name", "Test",
        "--genome-fasta", str(fasta), "--out-prefix", str(prefix),
        "--evidence-proteins", str(protein),
        "--pruning-level", "none", "--no-alternative-splice",
        "--allow-internal-stop-evidence-rescue",
        "--minimum-coding-length-bp", "1",
        "--output-profile", "standard", "--progress-log-level", "none",
        "--stop-after", "scored-candidates",
    ]) == 0
    scores = Path(str(prefix) + ".model_selection_scores.tsv").read_text()
    assert "t_bad" not in scores
    assert "t_good" in scores
    with Path(str(prefix) + ".definite_internal_stop_exclusions.tsv").open(
        newline=""
    ) as handle:
        ledger = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["original_transcript_id"] for row in ledger] == ["t_bad"]
    assert Path(
        str(prefix) + ".definite_internal_stop_exclusion.summary.tsv"
    ).exists()
    active_configuration = json.loads(
        Path(str(prefix) + ".run_configuration.json").read_text()
    )
    assert active_configuration["definite_internal_stop_exclusion"]["requested_enabled"] is True
    assert active_configuration["definite_internal_stop_exclusion"]["effective_enabled"] is True

    legacy_prefix = tmp_path / "legacy"
    assert pipeline_module.main([
        "--inputs", str(gff), "--source-name", "Test",
        "--genome-fasta", str(fasta), "--out-prefix", str(legacy_prefix),
        "--evidence-proteins", str(protein),
        "--pruning-level", "none", "--no-alternative-splice",
        "--no-definite-internal-stop-exclusion",
        "--minimum-coding-length-bp", "1",
        "--output-profile", "standard", "--progress-log-level", "none",
        "--stop-after", "scored-candidates",
    ]) == 0
    legacy_scores = Path(
        str(legacy_prefix) + ".model_selection_scores.tsv"
    ).read_text()
    assert "t_bad" in legacy_scores
    assert "t_good" in legacy_scores
    with Path(
        str(legacy_prefix) + ".definite_internal_stop_exclusion.summary.tsv"
    ).open(newline="") as handle:
        legacy_summary = next(csv.DictReader(handle, delimiter="\t"))
    assert legacy_summary["execution_status"] == "disabled_by_user"
    assert legacy_summary["excluded_models"] == "0"
    configuration = json.loads(
        Path(str(legacy_prefix) + ".run_configuration.json").read_text()
    )
    assert configuration["definite_internal_stop_exclusion"]["requested_enabled"] is False
    assert configuration["definite_internal_stop_exclusion"]["effective_enabled"] is False

    final_prefix = tmp_path / "final"
    assert pipeline_module.main([
        "--inputs", str(gff), "--source-name", "Test",
        "--genome-fasta", str(fasta), "--out-prefix", str(final_prefix),
        "--evidence-proteins", str(protein),
        "--pruning-level", "none", "--no-refine-cds-termini",
        "--minimum-coding-length-bp", "1",
        "--alternative-splice-report", "--output-profile", "full",
        "--progress-log-level", "none",
    ]) == 0
    for suffix in (
        "superlocus_members.tsv",
        "sublocus_members.tsv",
        "model_structural_support.tsv",
        "model_selection_scores.tsv",
        "selected_models.tsv",
        "backfill_candidates.tsv",
        "final_selected_models.tsv",
        "alternative_splice_candidates.tsv",
        "final_annotation.gff3",
    ):
        text = Path(str(final_prefix) + "." + suffix).read_text()
        assert "t_bad" not in text, suffix


def test_invalid_overlap_is_excluded_before_all_candidate_stages_and_not_rescued(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "overlap_models.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t21\t34\t.\t+\t.\tID=g_overlap\n"
        "chr1\ttest\tmRNA\t21\t34\t.\t+\t.\tID=t_overlap;Parent=g_overlap\n"
        "chr1\ttest\texon\t21\t27\t.\t+\t.\tParent=t_overlap\n"
        "chr1\ttest\tCDS\t21\t27\t.\t+\t.\tParent=t_overlap\n"
        "chr1\ttest\texon\t27\t34\t.\t+\t.\tParent=t_overlap\n"
        "chr1\ttest\tCDS\t27\t34\t.\t+\t.\tParent=t_overlap\n"
        "chr1\ttest\tgene\t81\t89\t.\t+\t.\tID=g_good\n"
        "chr1\ttest\tmRNA\t81\t89\t.\t+\t.\tID=t_good;Parent=g_good\n"
        "chr1\ttest\texon\t81\t89\t.\t+\t.\tParent=t_good\n"
        "chr1\ttest\tCDS\t81\t89\t.\t+\t0\tParent=t_good\n"
    )
    genome = list("C" * 140)
    genome[20:34] = "TAGTACTGATTGAT"
    genome[80:89] = "ATGCCCTAA"
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "".join(genome) + "\n")
    protein = tmp_path / "protein.gff3"
    protein.write_text(
        "##gff-version 3\n"
        "chr1\tprotein\tCDS\t21\t34\t.\t+\t0\tID=p1.cds;Parent=p1\n"
    )
    prefix = tmp_path / "overlap_out"
    assert pipeline_module.main([
        "--inputs", str(gff), "--source-name", "Test",
        "--genome-fasta", str(fasta), "--out-prefix", str(prefix),
        "--evidence-proteins", str(protein),
        "--allow-internal-stop-evidence-rescue",
        "--minimum-coding-length-bp", "1",
        "--pruning-level", "none", "--no-refine-cds-termini",
        "--no-alternative-splice", "--output-profile", "full",
        "--progress-log-level", "none",
    ]) == 0

    with Path(
        str(prefix) + ".definite_internal_stop_exclusions.tsv"
    ).open(newline="") as handle:
        ledger = list(csv.DictReader(handle, delimiter="\t"))
    assert len(ledger) == 1
    assert ledger[0]["original_transcript_id"] == "t_overlap"
    assert ledger[0]["exclusion_class"] == INVALID_CDS_GEOMETRY
    with Path(
        str(prefix) + ".definite_internal_stop_exclusion.summary.tsv"
    ).open(newline="") as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    assert summary["excluded_models"] == "1"
    assert summary["excluded_invalid_cds_geometry"] == "1"
    assert summary["excluded_unresolved_all_frames_internal_stop"] == "0"

    for suffix in (
        "superlocus_members.tsv",
        "sublocus_members.tsv",
        "model_structural_support.tsv",
        "model_selection_scores.tsv",
        "selected_models.tsv",
        "backfill_candidates.tsv",
        "final_selected_models.tsv",
    ):
        with Path(str(prefix) + "." + suffix).open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert all("t_overlap" not in row.values() for row in rows), suffix
    final_gff = Path(str(prefix) + ".final_annotation.gff3").read_text()
    assert "t_overlap" not in final_gff
    assert "t_good" in final_gff


def test_all_excluded_run_writes_h3_reports_before_no_candidate_exit(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "only_bad.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t12\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t1\t12\t.\t+\t.\tParent=t1\n"
        "chr1\ttest\tCDS\t1\t12\t.\t+\t0\tParent=t1\n"
    )
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nATGTAACCCTAA\n")
    prefix = tmp_path / "only_bad"
    with pytest.raises(SystemExit, match="No valid coding models remained"):
        pipeline_module.main([
            "--inputs", str(gff), "--source-name", "Test",
            "--genome-fasta", str(fasta), "--out-prefix", str(prefix),
            "--output-profile", "minimal", "--progress-log-level", "none",
            "--write-output", "definite_internal_stop_exclusions_tsv",
        ])
    summary_path = Path(
        str(prefix) + ".definite_internal_stop_exclusion.summary.tsv"
    )
    ledger_path = Path(str(prefix) + ".definite_internal_stop_exclusions.tsv")
    assert summary_path.exists()
    assert ledger_path.exists()
    with summary_path.open(newline="") as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    assert summary["excluded_models"] == "1"
