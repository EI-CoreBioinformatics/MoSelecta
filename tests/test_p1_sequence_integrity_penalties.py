"""Focused tests for the opt-in P1 terminal-codon penalty configuration."""

from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.config import build_arg_parser
from moselecta.constants import (
    MODEL_SELECTION_DEFAULT_WEIGHTS,
    SEQUENCE_INTEGRITY_P0_MISSING_START_PENALTY,
    SEQUENCE_INTEGRITY_P0_MISSING_STOP_PENALTY,
    SEQUENCE_INTEGRITY_P1_MISSING_START_PENALTY,
    SEQUENCE_INTEGRITY_P1_MISSING_STOP_PENALTY,
)
from moselecta.genome_sequence import (
    calculate_model_sequence_integrity,
    unavailable_sequence_integrity_row,
)
from moselecta.models import CodingModel, FeatureRecord, Sublocus
from moselecta.run_configuration import build_run_configuration
from moselecta.scoring import (
    apply_model_selection_scores,
    model_selection_score_positive_scale_max,
)
from moselecta.selection import stage5_backfill_selection


def cds_record(start0: int, end0: int, phase: str = "0") -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1", source="test", feature_type="CDS", start=start0 + 1,
        end=end0, score=".", strand="+", phase=phase, attrs={"Parent": "tx"},
        raw_attrs="Parent=tx", input_file=Path("test.gff3"), source_name="test",
        line_number=1, record_index=1, raw_line="",
    )


def model(
    model_id: str,
    parts: tuple[tuple[int, int], ...],
    *,
    source: str = "A",
    phases: tuple[str, ...] | None = None,
) -> CodingModel:
    phases = phases or tuple("0" for _ in parts)
    item = CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"{model_id}.gene",
        input_transcript_id=f"{model_id}.input",
        output_gene_id=f"{model_id}.gene.out",
        output_transcript_id=model_id,
        chrom="chr1",
        strand="+",
        cds_parts=list(parts),
        gene_record=None,
        transcript_record=None,
        child_records=[
            cds_record(start0, end0, phase)
            for (start0, end0), phase in zip(parts, phases)
        ],
    )
    item.superlocus_id = "super1"
    item.sublocus_id = f"super1.{model_id}"
    return item


def sequence_row(
    coding_sequence: str,
    *,
    p1: bool,
    adjacent_sequence: str = "CCC",
    start0: int = 10,
    contig_length: int = 80,
) -> dict[str, object]:
    chars = list("C" * contig_length)
    chars[start0 : start0 + len(coding_sequence)] = coding_sequence
    adjacent_start = start0 + len(coding_sequence)
    chars[adjacent_start : adjacent_start + len(adjacent_sequence)] = adjacent_sequence
    return calculate_model_sequence_integrity(
        model("sequence", ((start0, adjacent_start),)),
        window="".join(chars),
        window_start0=0,
        contig_length=contig_length,
        missing_start_penalty=(
            SEQUENCE_INTEGRITY_P1_MISSING_START_PENALTY
            if p1 else SEQUENCE_INTEGRITY_P0_MISSING_START_PENALTY
        ),
        missing_stop_penalty=(
            SEQUENCE_INTEGRITY_P1_MISSING_STOP_PENALTY
            if p1 else SEQUENCE_INTEGRITY_P0_MISSING_STOP_PENALTY
        ),
    )


@pytest.mark.parametrize(
    ("sequence", "p0", "p1", "flags"),
    [
        ("CCCAAATAA", 0.30, 0.50, "missing_start_codon"),
        ("ATGAAACCC", 0.50, 0.75, "missing_stop_codon"),
        ("CCCAAACCC", 0.80, 1.25, "missing_start_codon,missing_stop_codon"),
        ("ATGAAATAA", 0.00, 0.00, "."),
    ],
)
def test_p1_changes_only_missing_terminal_codon_magnitudes(
    sequence: str, p0: float, p1: float, flags: str
) -> None:
    p0_row = sequence_row(sequence, p1=False)
    p1_row = sequence_row(sequence, p1=True)
    assert p0_row["sequence_integrity_penalty_score"] == pytest.approx(p0)
    assert p1_row["sequence_integrity_penalty_score"] == pytest.approx(p1)
    assert p0_row["sequence_integrity_flags"] == flags
    assert p1_row["sequence_integrity_flags"] == flags


def test_p1_leaves_ambiguous_and_boundary_terminal_codons_neutral() -> None:
    ambiguous_start = sequence_row("NNNAAATAA", p1=True)
    ambiguous_stop = sequence_row("ATGAAA", p1=True, adjacent_sequence="NNN")
    boundary_start = sequence_row("CCCAAATAA", p1=True, start0=0)
    boundary_stop = sequence_row(
        "ATGAAACCC", p1=True, start0=71, contig_length=80, adjacent_sequence=""
    )
    assert ambiguous_start["sequence_start_codon_status"] == "ambiguous"
    assert ambiguous_stop["sequence_stop_codon_status"] == "ambiguous"
    assert boundary_start["sequence_start_codon_status"] == "not_evaluable_at_contig_boundary"
    assert boundary_stop["sequence_stop_codon_status"] == "not_evaluable_at_contig_boundary"
    for row in (ambiguous_start, ambiguous_stop, boundary_start, boundary_stop):
        assert row["sequence_integrity_penalty_score"] == 0.0


def test_p1_leaves_internal_stop_phase_splice_and_unavailable_cases_unchanged() -> None:
    internal_p0 = sequence_row("ATGTAACCCTAA", p1=False)
    internal_p1 = sequence_row("ATGTAACCCTAA", p1=True)
    assert internal_p0["sequence_integrity_penalty_score"] == 1.0
    assert internal_p1["sequence_integrity_penalty_score"] == 1.0

    chars = list("C" * 80)
    chars[10:16] = "ATGAAA"
    chars[16:18] = "GT"
    chars[28:30] = "AG"
    chars[30:36] = "CCCTAA"
    phase_model = model("phase", ((10, 16), (30, 36)), phases=("0", "1"))
    phase_rows = [
        calculate_model_sequence_integrity(
            phase_model, window="".join(chars), window_start0=0,
            contig_length=80, missing_start_penalty=start, missing_stop_penalty=stop,
        )
        for start, stop in ((0.30, 0.50), (0.50, 0.75))
    ]
    assert [row["sequence_integrity_penalty_score"] for row in phase_rows] == [0.75, 0.75]

    chars[16:18] = "CC"
    chars[28:30] = "CC"
    splice_model = model("splice", ((10, 16), (30, 36)), phases=("0", "0"))
    splice_rows = [
        calculate_model_sequence_integrity(
            splice_model, window="".join(chars), window_start0=0,
            contig_length=80, missing_start_penalty=start, missing_stop_penalty=stop,
        )
        for start, stop in ((0.30, 0.50), (0.50, 0.75))
    ]
    assert [row["sequence_integrity_penalty_score"] for row in splice_rows] == [1.0, 1.0]

    unavailable = unavailable_sequence_integrity_row(
        model("unavailable", ((10, 19),)), "not_evaluated", "genome_unavailable"
    )
    assert unavailable["sequence_integrity_penalty_score"] == 0.0


def scoring_row(
    model_id: str,
    *,
    structural_support: float,
    pairwise_similarity: float,
    sequence_penalty: float,
    source: str = "A",
) -> dict[str, object]:
    return {
        "sublocus_id": "super1.rank",
        "sublocus_source_count": 3,
        "sublocus_model_count": 2,
        "pairwise_other_source_count": 1,
        "source_name": source,
        "output_transcript_id": model_id,
        "structural_support_score": structural_support,
        "consensus_similarity_score": 0.0,
        "pairwise_similarity_score": pairwise_similarity,
        "model_completeness_score": 0.0,
        "terminal_extension_score": 0.0,
        "exact_structure_support_score": 0.0,
        "corroborating_other_source_count": 1 if model_id == "clean" else 0,
        "run_possible_other_source_count": 4,
        "unsupported_short_cds_exon_penalty_score": 0.0,
        "structural_coherence_raw_score": 0.0,
        "sequence_integrity_penalty_score": sequence_penalty,
        "intrinsic_penalty_score": 0.25,
        "partial_structure_penalty": 0.0,
        "internal_extra_fraction": 0.0,
        "terminal_extra_fraction": 0.0,
        "cds_length": 100,
    }


def score_pair(missing_penalty: float) -> dict[str, dict[str, object]]:
    rows = [
        scoring_row(
            "clean", structural_support=1.0, pairwise_similarity=0.0,
            sequence_penalty=0.0, source="C",
        ),
        scoring_row(
            "missing", structural_support=0.7, pairwise_similarity=0.7,
            sequence_penalty=missing_penalty, source="M",
        ),
    ]
    weights = dict(MODEL_SELECTION_DEFAULT_WEIGHTS)
    weights["structural_corroboration"] = 0.0
    apply_model_selection_scores(
        rows,
        source_prior_scores={"C": 0.0, "M": 0.0},
        selection_weights=weights,
        tie_epsilon=1e-6,
    )
    return {str(row["output_transcript_id"]): row for row in rows}


def test_p1_score_difference_reporting_ranking_and_positive_maximum() -> None:
    p0 = score_pair(0.30)
    p1 = score_pair(0.50)
    assert p1["missing"]["model_selection_score"] == pytest.approx(
        float(p0["missing"]["model_selection_score"]) - 0.20
    )
    assert p1["missing"]["sequence_integrity_penalty_contribution"] == 0.50
    assert "-sequence_integrity_penalty_score:0.500000*1.000" in str(
        p1["missing"]["model_selection_score_components"]
    )
    assert p0["missing"]["model_selection_order_rank"] == 1
    assert p1["clean"]["model_selection_order_rank"] == 1
    assert p0["missing"]["intrinsic_penalty_score"] == p1["missing"]["intrinsic_penalty_score"] == 0.25

    weights = dict(MODEL_SELECTION_DEFAULT_WEIGHTS)
    weights["structural_corroboration"] = 0.5
    assert model_selection_score_positive_scale_max(weights) == pytest.approx(9.85)


def as_sublocus(item: CodingModel) -> Sublocus:
    return Sublocus(
        item.sublocus_id, "super1", item.chrom, item.strand,
        item.cds_start0, item.cds_end0, [item],
        "single_exon" if item.is_single_exon else "multi_exon",
    )


def test_p1_changes_backfill_order_through_existing_model_selection_scores() -> None:
    clean = model("clean", ((10, 19),), source="C")
    missing = model("missing", ((30, 39),), source="M")
    orders: list[list[str]] = []
    for penalty in (0.30, 0.50):
        rows = score_pair(penalty)
        result = stage5_backfill_selection(
            [as_sublocus(clean), as_sublocus(missing)],
            rows,
            {"selected_models": [], "suppressed_rows": [], "selection_info_by_transcript": {}},
            min_cds_exon_overlap_fraction=0.10,
        )
        orders.append([
            str(row["output_transcript_id"])
            for row in sorted(
                result["backfilled_rows"],  # type: ignore[arg-type]
                key=lambda row: int(row["backfill_selection_order"]),
            )
        ])
    assert orders == [["missing", "clean"], ["clean", "missing"]]


def test_r1me_comparison_uses_p1_revised_model_selection_score() -> None:
    candidate = model("clean", ((100, 200), (300, 400)), source="C")
    blocker = model("missing", ((120, 140), (160, 180)), source="M")
    outcomes: list[set[str]] = []
    for penalty in (0.30, 0.50):
        rows = score_pair(penalty)
        rows["clean"]["corroborating_other_source_count"] = 1
        rows["missing"]["corroborating_other_source_count"] = 0
        result = stage5_backfill_selection(
            [as_sublocus(candidate), as_sublocus(blocker)],
            rows,
            {
                "selected_models": [blocker],
                "suppressed_rows": [],
                "selection_info_by_transcript": {"missing": {"monosublocus_selected": "true"}},
            },
            min_cds_exon_overlap_fraction=0.10,
            strict_r1_backfill_replacement=True,
            strict_r1_allow_multiexon_blockers=True,
            sequence_penalty_missing_start=penalty,
        )
        outcomes.append({item.output_transcript_id for item in result["final_models"]})  # type: ignore[index,union-attr]
    assert outcomes == [{"missing"}, {"clean"}]


def test_p1_defaults_p0_override_and_serialization() -> None:
    parser = build_arg_parser(full_help=True)
    required = ["--inputs", "input.gff3", "--out-prefix", "out"]
    p1 = parser.parse_args(required)
    p0 = parser.parse_args([
        *required,
        "--sequence-penalty-missing-start", "0.30",
        "--sequence-penalty-missing-stop", "0.50",
    ])
    assert (p1.sequence_penalty_missing_start, p1.sequence_penalty_missing_stop) == (0.50, 0.75)
    assert (p0.sequence_penalty_missing_start, p0.sequence_penalty_missing_stop) == (0.30, 0.50)

    configuration = build_run_configuration(
        p1,
        raw_arguments=[],
        enabled_output_keys=set(),
        effective_pruning_level=None,
    )
    assert configuration["scoring"]["sequence_penalty_missing_start"] == 0.50
    assert configuration["scoring"]["sequence_penalty_missing_stop"] == 0.75
