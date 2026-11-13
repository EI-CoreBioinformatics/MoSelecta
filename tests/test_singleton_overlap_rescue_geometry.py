from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.config import build_arg_parser
from moselecta.models import CodingModel, Sublocus
from moselecta.selection import final_prune_singleton_selected_models


def _model(transcript_id: str = "candidate", source: str = "SourceA") -> CodingModel:
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"gene_{transcript_id}",
        input_transcript_id=transcript_id,
        output_gene_id=f"out_gene_{transcript_id}",
        output_transcript_id=transcript_id,
        chrom="chr1",
        strand="+",
        cds_parts=[(100, 250)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def _edge(
    partner_id: str,
    partner_source: str,
    self_fraction: float,
    partner_fraction: float,
    overlap_bp: int = 150,
) -> dict[str, object]:
    return {
        "monosublocus_id": "mono1",
        "model_a_id": "candidate",
        "model_b_id": partner_id,
        "source_a": "SourceA",
        "source_b": partner_source,
        "cds_exon_overlap_bp": overlap_bp,
        "overlap_fraction_of_a": self_fraction,
        "overlap_fraction_of_b": partner_fraction,
    }


def _run(
    edges: list[dict[str, object]],
    *,
    holder_member_count: int | None = None,
    min_overlap_bp: int = 150,
    min_self_fraction: float = 0.80,
    min_reciprocal_fraction: float = 0.80,
    policy: str = "adequacy-conditioned",
    protein_adequacy: bool = False,
    transcript_adequacy: bool = False,
    protein_score: float = 0.0,
    transcript_score: float = 0.0,
) -> dict[str, object]:
    candidate = _model()
    sub = Sublocus("sub1", "sl1", "chr1", "+", 100, 250, [candidate])
    member_ids = ["candidate", *[str(edge["model_b_id"]) for edge in edges]]
    member_count = holder_member_count if holder_member_count is not None else len(set(member_ids))
    stage4_result = {
        "monosublocus_rows": [{"monosublocus_id": "mono1", "member_count": member_count}],
        "monosublocus_member_rows": [
            {"monosublocus_id": "mono1", "output_transcript_id": model_id}
            for model_id in dict.fromkeys(member_ids)
        ],
        "monosublocus_edge_rows": edges,
    }
    row = {
        "output_transcript_id": "candidate",
        "source_name": "SourceA",
        "final_prune_classes": ".",
    }
    return final_prune_singleton_selected_models(
        {
            "initial_retained_models": [candidate],
            "initial_retained_rows": [row],
            "retained_models": [candidate],
            "retained_rows": [row],
            "removed_rows": [],
            "retained_selection_info_by_transcript": {"candidate": row},
        },
        [sub],
        {"candidate": {"intrinsic_penalty_components": ".", "intrinsic_penalty_score": 0.0}},
        prune_level="all",
        prune_sources=set(),
        pruning_level_name="high",
        rescue_monosublocus_overlap=True,
        stage4_result=stage4_result,
        rescue_min_holder_member_count=2,
        rescue_min_overlap_bp=min_overlap_bp,
        rescue_min_self_overlap_fraction=min_self_fraction,
        rescue_min_reciprocal_overlap_fraction=min_reciprocal_fraction,
        monosublocus_overlap_rescue_evidence_policy=policy,
        protein_adequacy_passed=protein_adequacy,
        transcript_adequacy_passed=transcript_adequacy,
        rescue_protein_supported_singletons=protein_score > 0,
        protein_evidence_rows_by_transcript={
            "candidate": {
                "protein_evidence_score_scaled": protein_score,
                "protein_exact_cds_chain_source_count": 0,
            }
        },
        rescue_transcript_supported_singletons=transcript_score > 0,
        transcript_evidence_rows_by_transcript={
            "candidate": {
                "transcript_evidence_score_scaled": transcript_score,
                "transcript_supporting_source_count": 1,
                "transcript_scoring_eligible_count": 1,
                "transcript_stranded_scoring_eligible_count": 1,
                "transcript_unstranded_scoring_eligible_count": 0,
                "transcript_cdna_context_evaluable": "true",
                "transcript_cdna_context_cds_to_min_cdna_fraction": 0.80,
                "transcript_cdna_context_pruning_concern": "false",
            }
        },
    )


def _result_row(result: dict[str, object]) -> dict[str, object]:
    rows = result["retained_rows"] or result["removed_rows"]
    return rows[0]


def test_one_sided_containment_fails_reciprocal_geometry() -> None:
    result = _run([_edge("partner", "SourceB", 0.95, 0.30)])
    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_singleton_overlap_rescue_geometry_reason"] == "reciprocal_overlap_below_threshold"
    assert row["final_singleton_rescue_reciprocal_overlap_fraction"] == pytest.approx(0.30)
    assert result["singleton_overlap_rescue_reciprocal_rejected_count"] == 1


def test_reciprocal_085_different_source_geometry_rescues() -> None:
    result = _run([_edge("partner", "SourceB", 0.85, 0.85)])
    row = result["retained_rows"][0]
    assert row["final_singleton_rescue_reason"] == "monosublocus_overlap"
    assert row["final_singleton_overlap_rescue_geometry_passed"] == "true"
    assert row["final_singleton_rescue_reciprocal_overlap_fraction"] == pytest.approx(0.85)


def test_self_overlap_079_fails_existing_and_reciprocal_requirements() -> None:
    result = _run([_edge("partner", "SourceB", 0.79, 1.00)])
    assert result["retained_models"] == []
    assert result["removed_rows"][0]["final_singleton_overlap_rescue_geometry_reason"] == "self_overlap_below_threshold"


@pytest.mark.parametrize(
    ("fraction", "rescued"),
    [(0.80, True), (0.799999, False)],
)
def test_reciprocal_boundary_uses_normal_tolerance(fraction: float, rescued: bool) -> None:
    result = _run([_edge("partner", "SourceB", 0.95, fraction)])
    assert bool(result["retained_models"]) is rescued


def test_existing_overlap_bp_and_holder_member_requirements_remain_effective() -> None:
    low_bp = _run([_edge("partner", "SourceB", 0.95, 0.95, overlap_bp=149)])
    assert low_bp["retained_models"] == []
    assert _result_row(low_bp)["final_singleton_overlap_rescue_geometry_reason"] == "overlap_bp_below_threshold"

    small_holder = _run(
        [_edge("partner", "SourceB", 0.95, 0.95)],
        holder_member_count=1,
    )
    assert small_holder["retained_models"] == []
    assert _result_row(small_holder)["final_singleton_overlap_rescue_geometry_reason"] == "holder_member_count_below_threshold"


def test_same_source_partner_cannot_rescue() -> None:
    result = _run([_edge("partner", "SourceA", 0.95, 0.95)])
    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_singleton_overlap_rescue_geometry_reason"] == "no_qualifying_different_source_partner"
    assert row["final_singleton_overlap_rescue_same_source_edge"] == "true"


def test_qualifying_different_source_partner_wins_over_same_source_partner() -> None:
    result = _run([
        _edge("same", "SourceA", 1.00, 1.00, overlap_bp=200),
        _edge("different", "SourceB", 0.85, 0.85),
    ])
    row = result["retained_rows"][0]
    assert row["final_singleton_rescue_partner_model_id"] == "different"
    assert row["final_singleton_rescue_partner_source"] == "SourceB"
    assert row["final_singleton_overlap_rescue_same_source_edge"] == "true"
    assert row["final_singleton_overlap_rescue_geometry_reason"] == "passed"
    assert result["singleton_overlap_rescue_models_with_same_source_edge_count"] == 1


def test_no_partner_meeting_both_new_requirements_does_not_rescue() -> None:
    result = _run([
        _edge("same", "SourceA", 0.95, 0.95),
        _edge("low_reciprocal", "SourceB", 0.95, 0.30),
    ])
    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_singleton_overlap_rescue_reciprocal_rejected"] == "true"
    assert row["final_singleton_overlap_rescue_same_source_edge"] == "true"


def test_failed_diagnostic_partner_and_reason_come_from_same_furthest_edge() -> None:
    result = _run([
        _edge("high_fraction_low_bp", "SourceB", 0.99, 0.98, overlap_bp=149),
        _edge("reciprocal_failure", "SourceC", 0.90, 0.40, overlap_bp=175),
    ])

    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_singleton_rescue_partner_model_id"] == "reciprocal_failure"
    assert row["final_singleton_rescue_self_overlap_fraction"] == pytest.approx(0.90)
    assert row["final_singleton_rescue_partner_overlap_fraction"] == pytest.approx(0.40)
    assert row["final_singleton_rescue_reciprocal_overlap_fraction"] == pytest.approx(0.40)
    assert row["final_singleton_overlap_rescue_geometry_reason"] == "reciprocal_overlap_below_threshold"


def test_failed_diagnostic_tie_uses_existing_deterministic_edge_order() -> None:
    result = _run([
        _edge("lower_self", "SourceB", 0.90, 0.40, overlap_bp=190),
        _edge("higher_self", "SourceC", 0.95, 0.30, overlap_bp=150),
    ])

    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_singleton_rescue_partner_model_id"] == "higher_self"
    assert row["final_singleton_rescue_self_overlap_fraction"] == pytest.approx(0.95)
    assert row["final_singleton_rescue_partner_overlap_fraction"] == pytest.approx(0.30)
    assert row["final_singleton_rescue_reciprocal_overlap_fraction"] == pytest.approx(0.30)
    assert row["final_singleton_overlap_rescue_geometry_reason"] == "reciprocal_overlap_below_threshold"


def test_same_source_rejection_is_recorded_even_when_reciprocal_also_fails() -> None:
    result = _run([_edge("same_low_reciprocal", "SourceA", 0.95, 0.30)])

    assert result["retained_models"] == []
    row = result["removed_rows"][0]
    assert row["final_singleton_rescue_partner_model_id"] == "same_low_reciprocal"
    assert row["final_singleton_overlap_rescue_geometry_reason"] == "reciprocal_overlap_below_threshold"
    assert row["final_singleton_overlap_rescue_reciprocal_rejected"] == "true"
    assert row["final_singleton_overlap_rescue_same_source_edge"] == "true"
    assert result["singleton_overlap_rescue_reciprocal_rejected_count"] == 1
    assert result["singleton_overlap_rescue_models_with_same_source_edge_count"] == 1


def test_single_exon_vs_multiexon_structural_disagreement_is_not_a_geometry_gate() -> None:
    # The candidate is single-exon; the edge may represent a multi-exon partner.
    # Rescue uses only CDS territory geometry and raw source identity.
    result = _run([_edge("multiexon_partner", "SourceB", 0.90, 0.90)])
    assert len(result["retained_models"]) == 1


def test_adequacy_conditioned_policy_interactions_remain_unchanged() -> None:
    qualifying = [_edge("partner", "SourceB", 0.90, 0.90)]
    no_adequacy = _run(qualifying)
    assert len(no_adequacy["retained_models"]) == 1

    adequate_without_local = _run(qualifying, protein_adequacy=True)
    assert adequate_without_local["retained_models"] == []
    assert _result_row(adequate_without_local)["final_singleton_overlap_rescue_denial_reason"] == "adequate_evidence_without_strong_local_support"


def test_strong_protein_rescue_is_independent_of_failed_reciprocal_geometry() -> None:
    result = _run(
        [_edge("partner", "SourceB", 0.95, 0.30)],
        protein_adequacy=True,
        protein_score=0.80,
    )
    row = result["retained_rows"][0]
    assert row["final_singleton_rescue_reason"] == "protein_evidence"
    assert row["final_protein_singleton_rescue_applied"] == "true"
    assert row["final_singleton_overlap_rescue_geometry_passed"] == "false"


def test_strong_transcript_rescue_is_independent_of_failed_reciprocal_geometry() -> None:
    result = _run(
        [_edge("partner", "SourceB", 0.95, 0.30)],
        transcript_adequacy=True,
        transcript_score=0.80,
    )
    row = result["retained_rows"][0]
    assert row["final_singleton_rescue_reason"] == "transcript_evidence"
    assert row["final_transcript_singleton_rescue_applied"] == "true"
    assert row["final_singleton_overlap_rescue_geometry_passed"] == "false"


def test_unconditional_policy_bypasses_only_evidence_adequacy() -> None:
    failed_geometry = _run(
        [_edge("partner", "SourceB", 0.95, 0.30)],
        policy="unconditional",
        protein_adequacy=True,
    )
    assert failed_geometry["retained_models"] == []

    same_source = _run(
        [_edge("partner", "SourceA", 0.95, 0.95)],
        policy="unconditional",
        protein_adequacy=True,
    )
    assert same_source["retained_models"] == []

    qualifying_geometry = _run(
        [_edge("partner", "SourceB", 0.90, 0.90)],
        policy="unconditional",
        protein_adequacy=True,
    )
    assert len(qualifying_geometry["retained_models"]) == 1
    assert _result_row(qualifying_geometry)["final_singleton_overlap_rescue_decision"] == "rescued_unconditional"


def test_reciprocal_threshold_default_and_override_are_honoured() -> None:
    args = build_arg_parser(full_help=True).parse_args(
        ["--inputs", "input.gff3", "--out-prefix", "out"]
    )
    overridden_args = build_arg_parser(full_help=True).parse_args([
        "--inputs", "input.gff3", "--out-prefix", "out",
        "--singleton-rescue-min-reciprocal-overlap-fraction", "0.75",
    ])
    assert args.singleton_rescue_min_reciprocal_overlap_fraction == 0.80
    assert overridden_args.singleton_rescue_min_reciprocal_overlap_fraction == 0.75

    default_result = _run([_edge("partner", "SourceB", 0.85, 0.85)])
    strict_result = _run(
        [_edge("partner", "SourceB", 0.85, 0.85)],
        min_reciprocal_fraction=0.90,
    )
    assert len(default_result["retained_models"]) == 1
    assert strict_result["retained_models"] == []


def test_existing_self_overlap_override_remains_effective() -> None:
    result = _run(
        [_edge("partner", "SourceB", 0.85, 0.95)],
        min_self_fraction=0.90,
        min_reciprocal_fraction=0.80,
    )
    assert result["retained_models"] == []
    assert _result_row(result)["final_singleton_overlap_rescue_geometry_reason"] == "self_overlap_below_threshold"


def test_success_audit_identifies_selected_different_source_edge_and_threshold() -> None:
    result = _run([_edge("partner", "SourceB", 0.85, 0.90)])
    row = result["retained_rows"][0]
    assert row["final_singleton_rescue_self_overlap_fraction"] == pytest.approx(0.85)
    assert row["final_singleton_rescue_partner_overlap_fraction"] == pytest.approx(0.90)
    assert row["final_singleton_rescue_reciprocal_overlap_fraction"] == pytest.approx(0.85)
    assert row["final_singleton_rescue_min_reciprocal_overlap_fraction"] == pytest.approx(0.80)
    assert row["final_singleton_rescue_partner_model_id"] == "partner"
    assert row["final_singleton_rescue_partner_source"] == "SourceB"
    assert row["final_singleton_overlap_rescue_geometry_reason"] == "passed"
