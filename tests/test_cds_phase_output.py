from __future__ import annotations

from copy import deepcopy
from io import StringIO
from pathlib import Path

import pytest

from moselecta.cds_phase import (
    PRESERVED_NUMERIC_SOURCE,
    STANDARDISED_SOURCE_MISSING_TO_ZERO,
    UNRESOLVED_SOURCE_MISSING,
    calculate_cds_phase_chain,
    emitted_cds_phase_by_part,
)
from moselecta.genome_sequence import (
    calculate_model_sequence_integrity,
    infer_cds_phase_output,
    reverse_complement,
)
from moselecta.constants import MODEL_SELECTION_DEFAULT_WEIGHTS
from moselecta.annotation_outputs import (
    write_cds_only_transcript_children,
    write_source_model_children,
)
from moselecta.models import CDSPhaseOutput, CodingModel, FeatureRecord, Sublocus
from moselecta.scoring import apply_model_selection_scores
from moselecta.selection import select_sublocus_winners


def _record(start0: int, end0: int, phase: str, strand: str) -> FeatureRecord:
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


def _case(
    cds_sequence: str,
    *,
    phase: str = ".",
    strand: str = "+",
    adjacent: str = "CCC",
) -> tuple[CodingModel, str]:
    chars = list("C" * 80)
    if strand == "+":
        start0 = 20
        chars[start0 : start0 + len(cds_sequence)] = cds_sequence
        chars[start0 + len(cds_sequence) : start0 + len(cds_sequence) + 3] = adjacent
    else:
        start0 = 40
        genomic = reverse_complement(cds_sequence)
        chars[start0 : start0 + len(genomic)] = genomic
        chars[start0 - 3 : start0] = reverse_complement(adjacent)
    end0 = start0 + len(cds_sequence)
    model = CodingModel(
        source_name="test",
        input_file=Path("test.gff3"),
        input_gene_id="gene",
        input_transcript_id="tx",
        output_gene_id="out_gene",
        output_transcript_id="out_tx",
        chrom="chr1",
        strand=strand,
        cds_parts=[(start0, end0)],
        gene_record=None,
        transcript_record=None,
        child_records=[_record(start0, end0, phase, strand)],
    )
    return model, "".join(chars)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("phase", ["0", "1", "2"])
def test_numeric_source_first_phase_is_authoritative(strand: str, phase: str) -> None:
    model, window = _case("CATGAAATAA", phase=phase, strand=strand)
    decision = infer_cds_phase_output(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    model.cds_phase_output = decision
    assert decision.status == PRESERVED_NUMERIC_SOURCE
    assert decision.first_phase == int(phase)
    assert emitted_cds_phase_by_part(model)[model.cds_parts[0]] == phase


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    ("sequence", "adjacent", "reason"),
    [
        ("ATGAAATAA", "CCC", "phase_zero_start_and_stop_clean"),
        ("ATGAAAAAA", "CCC", "phase_zero_start_only_clean"),
        ("CCCAAATAA", "CCC", "phase_zero_stop_only_clean"),
        ("ATGCTAACCC", "TAA", "phase_zero_uniquely_clean_among_supported_frames"),
    ],
)
def test_missing_source_phase_standardises_only_clean_supported_zero_frame(
    strand: str,
    sequence: str,
    adjacent: str,
    reason: str,
) -> None:
    model, window = _case(sequence, strand=strand, adjacent=adjacent)
    decision = infer_cds_phase_output(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    model.cds_phase_output = decision
    assert decision.status == STANDARDISED_SOURCE_MISSING_TO_ZERO
    assert decision.first_phase == 0
    assert decision.reason == reason
    assert emitted_cds_phase_by_part(model)[model.cds_parts[0]] == "0"


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    ("sequence", "adjacent", "reason"),
    [
        ("CATGAAAAAA", "CCC", "only_nonzero_frame_supported"),
        ("CCATGAAAAAA", "CCC", "only_nonzero_frame_supported"),
        ("ATGCCCCCCC", "TAA", "competing_supported_clean_frame"),
        ("CCCCCCCCC", "CCC", "no_informative_terminus"),
        ("ATGTAATAA", "CCC", "phase_zero_internal_stop"),
    ],
)
def test_missing_source_phase_remains_missing_when_numeric_output_would_reinterpret(
    strand: str,
    sequence: str,
    adjacent: str,
    reason: str,
) -> None:
    model, window = _case(sequence, strand=strand, adjacent=adjacent)
    decision = infer_cds_phase_output(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    model.cds_phase_output = decision
    assert decision.status == UNRESOLVED_SOURCE_MISSING
    assert decision.first_phase is None
    assert decision.reason == reason
    assert emitted_cds_phase_by_part(model)[model.cds_parts[0]] == "."


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("adjacent", ["TAN", "TAR"])
def test_ambiguous_competing_adjacent_stop_does_not_veto_definite_phase_zero(
    strand: str,
    adjacent: str,
) -> None:
    model, window = _case(
        "ATGCCCCCCC", strand=strand, adjacent=adjacent
    )
    decision = infer_cds_phase_output(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    assert decision.status == STANDARDISED_SOURCE_MISSING_TO_ZERO
    assert decision.first_phase == 0
    assert decision.reason == "phase_zero_start_only_clean"


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("adjacent", ["TAA", "TAG", "TGA"])
def test_literal_competing_adjacent_stop_still_vetoes_phase_zero(
    strand: str,
    adjacent: str,
) -> None:
    model, window = _case(
        "ATGCCCCCCC", strand=strand, adjacent=adjacent
    )
    decision = infer_cds_phase_output(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    assert decision.status == UNRESOLVED_SOURCE_MISSING
    assert decision.first_phase is None
    assert decision.reason == "competing_supported_clean_frame"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_ambiguity_with_definite_included_competing_stop_remains_unresolved(
    strand: str,
) -> None:
    model, window = _case("ATGCCCCTAR", strand=strand)
    decision = infer_cds_phase_output(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    assert decision.status == UNRESOLVED_SOURCE_MISSING
    assert decision.first_phase is None
    assert decision.reason == "competing_supported_clean_frame"


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    ("sequence", "adjacent"),
    [
        ("ATGAAANAA", "CCC"),
        ("ATGCCCNCC", "CCC"),
        ("ATG" + "N" * 30 + "---CCC", "CCC"),
    ],
)
def test_ambiguity_away_from_decisive_evidence_allows_unique_phase_zero_output(
    strand: str, sequence: str, adjacent: str
) -> None:
    model, window = _case(sequence, strand=strand, adjacent=adjacent)
    before = calculate_model_sequence_integrity(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    decision = infer_cds_phase_output(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    model.cds_phase_output = decision
    after = calculate_model_sequence_integrity(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    assert decision.status == STANDARDISED_SOURCE_MISSING_TO_ZERO
    assert decision.first_phase == 0
    assert before == after


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    ("sequence", "reason"),
    [
        ("CANGCCTAA", "competing_supported_clean_frame"),
        ("NTGCCCCCC", "no_informative_terminus"),
        ("CCCCCCTAN", "no_informative_terminus"),
    ],
)
def test_ambiguity_that_can_change_frame_choice_remains_unresolved(
    strand: str, sequence: str, reason: str
) -> None:
    model, window = _case(sequence, strand=strand)
    decision = infer_cds_phase_output(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    assert decision.status == UNRESOLVED_SOURCE_MISSING
    assert decision.reason == reason


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("first_phase", [0, 1, 2])
def test_downstream_phase_chain_is_deterministic_on_both_strands(
    strand: str, first_phase: int
) -> None:
    parts = [(10, 12), (20, 22), (30, 32)]
    phases = calculate_cds_phase_chain(parts, strand, first_phase)
    ordered = sorted(parts, reverse=strand == "-")
    expected: list[str] = []
    cumulative = -first_phase
    for index, part in enumerate(ordered):
        expected.append(
            str(first_phase if index == 0 else (3 - cumulative % 3) % 3)
        )
        cumulative += part[1] - part[0]
    assert [phases[part] for part in ordered] == expected
    assert set(expected) == {"0", "1", "2"}


@pytest.mark.parametrize("cds_only", [False, True])
def test_gff_child_writers_share_the_approved_phase_chain(cds_only: bool) -> None:
    parts = [(10, 14), (20, 25), (30, 36)]
    model = CodingModel(
        source_name="test",
        input_file=Path("test.gff3"),
        input_gene_id="gene",
        input_transcript_id="tx",
        output_gene_id="out_gene",
        output_transcript_id="out_tx",
        chrom="chr1",
        strand="+",
        cds_parts=parts,
        gene_record=None,
        transcript_record=None,
        child_records=[_record(start, end, ".", "+") for start, end in parts],
        cds_phase_output=CDSPhaseOutput(
            status=STANDARDISED_SOURCE_MISSING_TO_ZERO,
            first_phase=0,
            reason="phase_zero_start_and_stop_clean",
        ),
    )
    out = StringIO()
    if cds_only:
        write_cds_only_transcript_children(out, model, "public_tx")
    else:
        write_source_model_children(out, model, "public_tx")
    phases = [
        fields[7]
        for line in out.getvalue().splitlines()
        if (fields := line.split("\t"))[2] == "CDS"
    ]
    assert phases == ["0", "2", "0"]


@pytest.mark.parametrize(
    ("sequence", "expected_reason"),
    [
        ("CATGAAAAAA", "only_nonzero_frame_supported"),
        ("CCATGAAAAAA", "only_nonzero_frame_supported"),
    ],
)
def test_output_phase_inference_does_not_change_integrity_interpretation(
    sequence: str, expected_reason: str
) -> None:
    model, window = _case(sequence)
    before = calculate_model_sequence_integrity(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    geometry = (tuple(model.cds_parts), tuple((r.start0, r.end0) for r in model.child_records))
    decision = infer_cds_phase_output(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    model.cds_phase_output = decision
    after = calculate_model_sequence_integrity(
        model, window=window, window_start0=0, contig_length=len(window)
    )
    assert decision.reason == expected_reason
    assert decision.first_phase is None
    assert before == after
    assert geometry == (
        tuple(model.cds_parts),
        tuple((r.start0, r.end0) for r in model.child_records),
    )


def test_output_phase_state_does_not_change_score_rank_or_selection() -> None:
    def score_row(transcript_id: str, completeness: float) -> dict[str, object]:
        return {
            "sublocus_id": "sub1",
            "sublocus_source_count": 1,
            "sublocus_model_count": 2,
            "pairwise_other_source_count": 0,
            "source_name": "test",
            "output_transcript_id": transcript_id,
            "structural_support_score": 0.0,
            "consensus_similarity_score": 0.0,
            "pairwise_similarity_score": 0.0,
            "model_completeness_score": completeness,
            "terminal_extension_score": 0.0,
            "exact_structure_support_score": 0.0,
            "unsupported_short_cds_exon_penalty_score": 0.0,
            "sequence_integrity_penalty_score": 0.0,
            "intrinsic_penalty_score": 0.0,
            "partial_structure_penalty": 0.0,
            "internal_extra_fraction": 0.0,
            "terminal_extra_fraction": 0.0,
            "cds_length": 9,
        }

    baseline = [score_row("winner", 1.0), score_row("loser", 0.0)]
    with_phase_audit = deepcopy(baseline)
    for row in with_phase_audit:
        row["cds_phase_output_status"] = STANDARDISED_SOURCE_MISSING_TO_ZERO
        row["cds_phase_output_first_phase"] = 0
        row["cds_phase_output_reason"] = "phase_zero_start_and_stop_clean"
    weights = dict(MODEL_SELECTION_DEFAULT_WEIGHTS)
    apply_model_selection_scores(baseline, {"test": 0.0}, weights, 1e-6)
    apply_model_selection_scores(with_phase_audit, {"test": 0.0}, weights, 1e-6)
    for plain, audited in zip(baseline, with_phase_audit):
        for key in (
            "model_structural_support_score",
            "model_structural_support_score_rank",
            "model_structural_support_order_rank",
            "model_selection_score",
            "model_selection_score_rank",
            "model_selection_order_rank",
        ):
            assert audited[key] == plain[key]

    winner, winner_window = _case("ATGCCCCCCC", adjacent="TAN")
    loser, _loser_window = _case("ATGCCCCCCC", adjacent="TAN")
    winner.output_transcript_id = "winner"
    loser.output_transcript_id = "loser"
    sublocus = Sublocus("sub1", "sl1", "chr1", "+", 20, 30, [winner, loser])
    rows_by_id = {str(row["output_transcript_id"]): row for row in baseline}
    selected_before = select_sublocus_winners([sublocus], rows_by_id)
    integrity_before = calculate_model_sequence_integrity(
        winner,
        window=winner_window,
        window_start0=0,
        contig_length=len(winner_window),
    )
    geometry_before = (
        tuple(winner.cds_parts),
        tuple((record.start0, record.end0) for record in winner.child_records),
    )
    winner.cds_phase_output = infer_cds_phase_output(
        winner,
        window=winner_window,
        window_start0=0,
        contig_length=len(winner_window),
    )
    integrity_after = calculate_model_sequence_integrity(
        winner,
        window=winner_window,
        window_start0=0,
        contig_length=len(winner_window),
    )
    selected_after = select_sublocus_winners([sublocus], rows_by_id)
    assert winner.cds_phase_output.status == STANDARDISED_SOURCE_MISSING_TO_ZERO
    assert integrity_after == integrity_before
    assert geometry_before == (
        tuple(winner.cds_parts),
        tuple((record.start0, record.end0) for record in winner.child_records),
    )
    assert [model.output_transcript_id for model in selected_before] == ["winner"]
    assert selected_after == selected_before
