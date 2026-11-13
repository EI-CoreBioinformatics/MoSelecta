"""Focused synthetic tests for strict R1 backfill replacement."""

from pathlib import Path

import pytest

from moselecta.models import CodingModel, Sublocus
from moselecta.selection import (
    make_initial_final_prune_result,
    selected_model_conflicts_for_backfill,
    stage5_backfill_selection,
    strict_r1_replacement_eligibility,
    write_r1_replacement_outputs,
)


def model(model_id: str, source: str, parts: tuple[tuple[int, int], ...]) -> CodingModel:
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
        child_records=[],
    )
    item.superlocus_id = "super1"
    item.sublocus_id = f"super1.{model_id}"
    return item


def sublocus(item: CodingModel) -> Sublocus:
    return Sublocus(
        item.sublocus_id,
        "super1",
        item.chrom,
        item.strand,
        item.cds_start0,
        item.cds_end0,
        [item],
        "single_exon" if item.is_single_exon else "multi_exon",
    )


def score_row(item: CodingModel, score: float, k: int) -> dict[str, object]:
    return {
        "output_transcript_id": item.output_transcript_id,
        "source_name": item.source_name,
        "sublocus_id": item.sublocus_id,
        "model_selection_score": score,
        "model_selection_score_scaled": max(0.0, score / 4.85),
        "corroborating_other_source_count": k,
        "run_possible_other_source_count": 6,
        "structural_corroboration_score": 0.0,
        "structural_corroboration_alpha": 16.0,
        "structural_corroboration_weight": 0.0,
        "structural_corroboration_contribution": 0.0,
        "structural_corroboration_in_selection": "false",
        "cds_length": item.cds_length,
    }


def stage4(blockers: list[CodingModel], *, suppressed: list[CodingModel] | None = None) -> dict[str, object]:
    suppressed = suppressed or []
    return {
        "selected_models": blockers,
        "suppressed_rows": [
            {"output_transcript_id": item.output_transcript_id} for item in suppressed
        ],
        "selection_info_by_transcript": {
            item.output_transcript_id: {
                "monosublocus_id": f"mono.{index}",
                "monosublocus_selected": "true",
            }
            for index, item in enumerate(blockers, start=1)
        },
    }


def run_r1(
    candidates: list[CodingModel],
    blockers: list[CodingModel],
    rows: dict[str, dict[str, object]],
    *,
    enabled: bool = True,
    allow_multiexon_blockers: bool = False,
    suppressed: list[CodingModel] | None = None,
    weight: float = 0.0,
    reverse_inputs: bool = False,
) -> tuple[dict[str, object], list[Sublocus]]:
    items = [*blockers, *candidates]
    subs = [sublocus(item) for item in items]
    if reverse_inputs:
        subs.reverse()
        blockers = list(reversed(blockers))
    result = stage5_backfill_selection(
        subs,
        rows,
        stage4(blockers, suppressed=suppressed),
        min_cds_exon_overlap_fraction=0.10,
        strict_r1_backfill_replacement=enabled,
        strict_r1_allow_multiexon_blockers=allow_multiexon_blockers,
        structural_corroboration_weight=weight,
        structural_corroboration_alpha=16.0,
    )
    return result, subs


def standard_models() -> tuple[CodingModel, CodingModel, CodingModel]:
    candidate = model("candidate", "C", ((100, 200), (300, 400)))
    blocker1 = model("blocker1", "A", ((120, 180),))
    blocker2 = model("blocker2", "B", ((320, 380),))
    return candidate, blocker1, blocker2


def rows_for(
    candidate: CodingModel,
    blockers: list[CodingModel],
    *,
    candidate_score: float = 3.0,
    candidate_k: int = 1,
    blocker_scores: list[float] | None = None,
    blocker_ks: list[int] | None = None,
) -> dict[str, dict[str, object]]:
    blocker_scores = blocker_scores or [2.0] * len(blockers)
    blocker_ks = blocker_ks or [0] * len(blockers)
    result = {candidate.output_transcript_id: score_row(candidate, candidate_score, candidate_k)}
    result.update({
        blocker.output_transcript_id: score_row(blocker, score, k)
        for blocker, score, k in zip(blockers, blocker_scores, blocker_ks)
    })
    return result


def final_ids(result: dict[str, object]) -> set[str]:
    return {item.output_transcript_id for item in result["final_models"]}  # type: ignore[index,union-attr]


def test_r1_replaces_one_uncorroborated_single_exon_winner() -> None:
    candidate, blocker, _ = standard_models()
    result, _subs = run_r1([candidate], [blocker], rows_for(candidate, [blocker]))
    assert final_ids(result) == {"candidate"}
    assert len(result["r1_replacement_rows"]) == 1  # type: ignore[arg-type]


def test_r1_replaces_two_blockers_atomically() -> None:
    candidate, blocker1, blocker2 = standard_models()
    result, _subs = run_r1(
        [candidate], [blocker1, blocker2], rows_for(candidate, [blocker1, blocker2])
    )
    assert final_ids(result) == {"candidate"}
    audit = result["r1_replacement_rows"][0]  # type: ignore[index]
    assert audit["replaced_model_count"] == 2
    assert audit["replaced_model_ids"] == "blocker1,blocker2"


@pytest.mark.parametrize("weight", [0.0, 0.5])
def test_r1_works_with_zero_or_active_corroboration_weight(weight: float) -> None:
    candidate, blocker, _ = standard_models()
    rows = rows_for(candidate, [blocker])
    rows["candidate"]["structural_corroboration_weight"] = weight
    result, _subs = run_r1([candidate], [blocker], rows, weight=weight)
    assert final_ids(result) == {"candidate"}
    assert result["r1_replacement_rows"][0]["structural_corroboration_weight"] == weight  # type: ignore[index]


@pytest.mark.parametrize(
    ("candidate_k", "candidate_score", "blocker_scores", "blocker_ks", "reason"),
    [
        (0, 3.0, [2.0], [0], "candidate_corroboration_count_lt_1"),
        (1, 3.0, [2.0], [1], "blocking_model_corroboration_count_ge_1"),
        (1, 2.0, [2.0], [0], "candidate_score_not_strictly_greater_than_every_blocker"),
        (1, 1.5, [2.0], [0], "candidate_score_not_strictly_greater_than_every_blocker"),
    ],
)
def test_r1_rejects_count_support_and_score_failures(
    candidate_k: int,
    candidate_score: float,
    blocker_scores: list[float],
    blocker_ks: list[int],
    reason: str,
) -> None:
    candidate, blocker, _ = standard_models()
    rows = rows_for(
        candidate,
        [blocker],
        candidate_k=candidate_k,
        candidate_score=candidate_score,
        blocker_scores=blocker_scores,
        blocker_ks=blocker_ks,
    )
    result, _subs = run_r1([candidate], [blocker], rows)
    assert final_ids(result) == {"blocker1"}
    assert result["backfill_rejected_rows"][0]["r1_replacement_reason"] == reason  # type: ignore[index]


def test_r1_rejects_selected_multi_exon_blocker() -> None:
    candidate, _blocker, _ = standard_models()
    blocker = model("multi_blocker", "A", ((120, 180), (300, 350)))
    result, _subs = run_r1([candidate], [blocker], rows_for(candidate, [blocker]))
    assert final_ids(result) == {"multi_blocker"}
    assert result["backfill_rejected_rows"][0]["r1_replacement_reason"] == "blocking_model_not_single_exon"  # type: ignore[index]


def test_r1_does_not_replace_with_single_exon_candidate() -> None:
    candidate = model("single_candidate", "C", ((100, 200),))
    blocker = model("single_blocker", "A", ((120, 180),))
    result, _subs = run_r1([candidate], [blocker], rows_for(candidate, [blocker]))
    assert final_ids(result) == {"single_blocker"}
    assert result["backfill_rejected_rows"][0]["r1_replacement_reason"] == "candidate_not_multi_exon"  # type: ignore[index]


def test_r1_rejects_backfilled_single_exon_blocker() -> None:
    candidate, blocker, _ = standard_models()
    rows = rows_for(candidate, [blocker])
    conflicts = selected_model_conflicts_for_backfill(
        candidate, [blocker], min_cds_exon_overlap_fraction=0.10
    )
    eligible, reason, _ = strict_r1_replacement_eligibility(
        candidate,
        conflicts,
        model_by_id={candidate.output_transcript_id: candidate, blocker.output_transcript_id: blocker},
        model_structural_rows_by_transcript=rows,
        monosublocus_selected_ids=set(),
    )
    assert not eligible
    assert reason == "blocking_model_not_monosublocus_selected"


def test_r1_rejects_three_blockers() -> None:
    candidate, blocker1, blocker2 = standard_models()
    blocker3 = model("blocker3", "D", ((150, 170),))
    blockers = [blocker1, blocker2, blocker3]
    result, _subs = run_r1([candidate], blockers, rows_for(candidate, blockers))
    assert final_ids(result) == {item.output_transcript_id for item in blockers}
    assert result["backfill_rejected_rows"][0]["r1_replacement_reason"] == "blocking_model_count_gt_2"  # type: ignore[index]


def test_r1_rejects_candidate_with_another_nonreplaceable_selected_conflict() -> None:
    candidate, blocker, _ = standard_models()
    other = model("other_multi", "D", ((100, 200), (300, 350)))
    blockers = [blocker, other]
    result, _subs = run_r1([candidate], blockers, rows_for(candidate, blockers))
    assert final_ids(result) == {"blocker1", "other_multi"}
    assert result["backfill_rejected_rows"][0]["r1_replacement_reason"] == "blocking_model_not_single_exon"  # type: ignore[index]


def test_r1_respects_shared_splice_site_conflict_with_selected_multi_exon() -> None:
    candidate, blocker, _ = standard_models()
    other = model("splice_conflict", "D", ((50, 200), (300, 450)))
    blockers = [blocker, other]
    result, _subs = run_r1([candidate], blockers, rows_for(candidate, blockers))
    rejected = result["backfill_rejected_rows"][0]  # type: ignore[index]
    assert "shared_cds_splice_site" in str(rejected["backfill_rejection_reason"])
    assert rejected["r1_replacement_reason"] == "blocking_model_not_single_exon"


def test_two_blocker_atomic_failure_preserves_both_blockers() -> None:
    candidate, blocker1, blocker2 = standard_models()
    blockers = [blocker1, blocker2]
    rows = rows_for(candidate, blockers, blocker_ks=[0, 1])
    result, _subs = run_r1([candidate], blockers, rows)
    assert final_ids(result) == {"blocker1", "blocker2"}
    assert result["r1_replacement_rows"] == []


def test_noneligible_suppressed_candidate_is_not_considered() -> None:
    candidate, blocker, _ = standard_models()
    result, _subs = run_r1(
        [candidate], [blocker], rows_for(candidate, [blocker]), suppressed=[candidate]
    )
    assert final_ids(result) == {"blocker1"}
    assert result["backfill_candidate_rows"] == []


def test_r1_disabled_preserves_existing_backfill_rejection() -> None:
    candidate, blocker, _ = standard_models()
    result, _subs = run_r1(
        [candidate], [blocker], rows_for(candidate, [blocker]), enabled=False
    )
    assert final_ids(result) == {"blocker1"}
    rejected = result["backfill_rejected_rows"][0]  # type: ignore[index]
    assert "r1_replacement_attempted" not in rejected


def test_r1_does_not_change_ordinary_nonconflicting_backfill_acceptance() -> None:
    candidate = model("candidate", "C", ((500, 600), (700, 800)))
    selected = model("selected", "A", ((100, 200),))
    rows = rows_for(candidate, [selected])
    result, _subs = run_r1([candidate], [selected], rows)
    assert final_ids(result) == {"candidate", "selected"}
    accepted = result["backfilled_rows"][0]  # type: ignore[index]
    assert accepted["final_selection_stage"] == "backfill"
    assert result["r1_replacement_rows"] == []


def test_successful_replacement_updates_selection_state_and_prevents_reintroduction() -> None:
    candidate, blocker, _ = standard_models()
    result, _subs = run_r1([candidate], [blocker], rows_for(candidate, [blocker]))
    assert "blocker1" not in final_ids(result)
    assert "blocker1" not in result["final_selection_info_by_transcript"]  # type: ignore[operator]
    assert result["r1_replaced_model_ids"] == ["blocker1"]


def test_later_backfill_decision_uses_revised_selected_population() -> None:
    candidate, blocker, _ = standard_models()
    later = model("later", "D", ((130, 170),))
    rows = rows_for(candidate, [blocker])
    rows[later.output_transcript_id] = score_row(later, 1.0, 0)
    result, _subs = run_r1([candidate, later], [blocker], rows)
    assert final_ids(result) == {"candidate"}
    assert any(row["output_transcript_id"] == "later" for row in result["backfill_rejected_rows"])  # type: ignore[index]


def test_r1_is_deterministic_under_input_iteration_order() -> None:
    candidate, blocker1, blocker2 = standard_models()
    blockers = [blocker1, blocker2]
    rows = rows_for(candidate, blockers)
    first, _ = run_r1([candidate], blockers, rows)
    second, _ = run_r1([candidate], blockers, rows, reverse_inputs=True)
    assert final_ids(first) == final_ids(second)
    assert first["r1_replacement_rows"] == second["r1_replacement_rows"]


def test_r1_provenance_and_summary_output(tmp_path: Path) -> None:
    candidate, blocker1, blocker2 = standard_models()
    result, _ = run_r1(
        [candidate], [blocker1, blocker2], rows_for(candidate, [blocker1, blocker2]), weight=0.5
    )
    replacements = tmp_path / "replacements.tsv"
    summary = tmp_path / "summary.tsv"
    write_r1_replacement_outputs(result, replacements_tsv=replacements, summary_tsv=summary)
    replacement_text = replacements.read_text()
    summary_text = summary.read_text()
    assert "candidate_model_selection_score" in replacement_text
    assert "blocker1,blocker2" in replacement_text
    assert "strict_r1_supported_multiexon_replaces_uncorroborated_single_exon_winner" in replacement_text
    assert "two_blocker_replacements\t1" in summary_text


def test_final_pruning_receives_revised_selected_population() -> None:
    candidate, blocker, _ = standard_models()
    rows = rows_for(candidate, [blocker])
    result, subs = run_r1([candidate], [blocker], rows)
    prune_input = make_initial_final_prune_result(result, subs, rows)
    retained = {item.output_transcript_id for item in prune_input["retained_models"]}  # type: ignore[index,union-attr]
    assert retained == {"candidate"}


def test_r1_summary_counts_rejection_conditions() -> None:
    candidate, blocker, _ = standard_models()
    result, _ = run_r1(
        [candidate], [blocker], rows_for(candidate, [blocker], candidate_k=0)
    )
    summary = result["r1_summary_rows"]  # type: ignore[assignment]
    assert {row["details"]: row["value"] for row in summary if row["metric"] == "rejections"} == {
        "candidate_corroboration_count_lt_1": 1
    }


def test_r1me_replaces_one_two_exon_holder() -> None:
    candidate, _blocker, _ = standard_models()
    blocker = model("two_exon_blocker", "A", ((120, 140), (160, 180)))
    result, _ = run_r1(
        [candidate], [blocker], rows_for(candidate, [blocker]),
        allow_multiexon_blockers=True,
    )
    assert final_ids(result) == {"candidate"}
    audit = result["r1_replacement_rows"][0]  # type: ignore[index]
    assert audit["candidate_cds_exon_count"] == 2
    assert audit["replaced_model_cds_exon_counts"] == "2"
    assert audit["r1_eligibility_class"] == "extended_multiexon_blocker_eligibility"


def test_r1me_replaces_one_three_exon_holder() -> None:
    candidate, _blocker, _ = standard_models()
    blocker = model("three_exon_blocker", "A", ((120, 135), (145, 160), (170, 185)))
    result, _ = run_r1(
        [candidate], [blocker], rows_for(candidate, [blocker]),
        allow_multiexon_blockers=True,
    )
    assert final_ids(result) == {"candidate"}
    assert result["r1_replacement_rows"][0]["replaced_model_cds_exon_counts"] == "3"  # type: ignore[index]


def test_r1me_replaces_two_multiexon_blockers_atomically() -> None:
    candidate, _blocker, _ = standard_models()
    left = model("left_multi", "A", ((120, 140), (160, 180)))
    right = model("right_multi", "B", ((320, 340), (360, 380)))
    blockers = [left, right]
    result, _ = run_r1(
        [candidate], blockers, rows_for(candidate, blockers),
        allow_multiexon_blockers=True,
    )
    assert final_ids(result) == {"candidate"}
    assert result["r1_replacement_rows"][0]["replaced_model_exon_classes"] == "multi_exon,multi_exon"  # type: ignore[index]


def test_r1me_replaces_mixed_single_and_multiexon_blockers() -> None:
    candidate, single, _ = standard_models()
    multi = model("right_multi", "B", ((320, 340), (360, 380)))
    blockers = [single, multi]
    result, _ = run_r1(
        [candidate], blockers, rows_for(candidate, blockers),
        allow_multiexon_blockers=True,
    )
    assert final_ids(result) == {"candidate"}
    assert result["r1_replacement_rows"][0]["replaced_model_exon_classes"] == "single_exon,multi_exon"  # type: ignore[index]


def test_r1me_does_not_require_blocker_cds_containment() -> None:
    candidate, _blocker, _ = standard_models()
    blocker = model("extending_blocker", "A", ((50, 150), (350, 450)))
    result, _ = run_r1(
        [candidate], [blocker], rows_for(candidate, [blocker]),
        allow_multiexon_blockers=True,
    )
    assert blocker.cds_start0 < candidate.cds_start0
    assert blocker.cds_end0 > candidate.cds_end0
    assert final_ids(result) == {"candidate"}


def test_r1me_separate_subloci_are_eligible() -> None:
    candidate, _blocker, _ = standard_models()
    blocker = model("multi_blocker", "A", ((120, 140), (160, 180)))
    assert candidate.sublocus_id != blocker.sublocus_id
    result, _ = run_r1(
        [candidate], [blocker], rows_for(candidate, [blocker]),
        allow_multiexon_blockers=True,
    )
    assert final_ids(result) == {"candidate"}


def test_r1me_retains_original_single_exon_eligibility_class() -> None:
    candidate, blocker, _ = standard_models()
    result, _ = run_r1(
        [candidate], [blocker], rows_for(candidate, [blocker]),
        allow_multiexon_blockers=True,
    )
    audit = result["r1_replacement_rows"][0]  # type: ignore[index]
    assert audit["r1_eligibility_class"] == "original_single_exon_blocker_eligibility"
    assert audit["r1_replacement_reason"] == (
        "strict_r1_supported_multiexon_replaces_uncorroborated_single_exon_winner"
    )


@pytest.mark.parametrize(
    ("candidate_k", "candidate_score", "blocker_score", "blocker_k", "reason"),
    [
        (0, 3.0, 2.0, 0, "candidate_corroboration_count_lt_1"),
        (1, 3.0, 2.0, 1, "blocking_model_corroboration_count_ge_1"),
        (1, 2.0, 2.0, 0, "candidate_score_not_strictly_greater_than_every_blocker"),
        (1, 1.9, 2.0, 0, "candidate_score_not_strictly_greater_than_every_blocker"),
    ],
)
def test_r1me_rejects_support_and_score_failures(
    candidate_k: int,
    candidate_score: float,
    blocker_score: float,
    blocker_k: int,
    reason: str,
) -> None:
    candidate, _blocker, _ = standard_models()
    blocker = model("multi_blocker", "A", ((120, 140), (160, 180)))
    rows = rows_for(
        candidate, [blocker], candidate_k=candidate_k, candidate_score=candidate_score,
        blocker_scores=[blocker_score], blocker_ks=[blocker_k],
    )
    result, _ = run_r1(
        [candidate], [blocker], rows, allow_multiexon_blockers=True,
    )
    assert final_ids(result) == {"multi_blocker"}
    assert result["backfill_rejected_rows"][0]["r1_replacement_reason"] == reason  # type: ignore[index]


def test_r1me_rejects_backfilled_multiexon_blocker() -> None:
    candidate, _blocker, _ = standard_models()
    blocker = model("backfilled_multi", "A", ((120, 140), (160, 180)))
    rows = rows_for(candidate, [blocker])
    conflicts = selected_model_conflicts_for_backfill(
        candidate, [blocker], min_cds_exon_overlap_fraction=0.10
    )
    eligible, reason, _ = strict_r1_replacement_eligibility(
        candidate,
        conflicts,
        model_by_id={candidate.output_transcript_id: candidate, blocker.output_transcript_id: blocker},
        model_structural_rows_by_transcript=rows,
        monosublocus_selected_ids=set(),
        allow_multiexon_blockers=True,
    )
    assert not eligible
    assert reason == "blocking_model_not_monosublocus_selected"


def test_r1me_rejects_more_than_two_mixed_blockers() -> None:
    candidate, single, _ = standard_models()
    left = model("left_multi", "A", ((120, 140), (160, 180)))
    right = model("right_multi", "B", ((320, 340), (360, 380)))
    blockers = [single, left, right]
    result, _ = run_r1(
        [candidate], blockers, rows_for(candidate, blockers),
        allow_multiexon_blockers=True,
    )
    assert final_ids(result) == {item.output_transcript_id for item in blockers}
    assert result["backfill_rejected_rows"][0]["r1_replacement_reason"] == "blocking_model_count_gt_2"  # type: ignore[index]


def test_r1me_atomic_failure_on_one_ineligible_multiexon_blocker() -> None:
    candidate, _blocker, _ = standard_models()
    left = model("left_multi", "A", ((120, 140), (160, 180)))
    right = model("right_multi", "B", ((320, 340), (360, 380)))
    blockers = [left, right]
    result, _ = run_r1(
        [candidate], blockers, rows_for(candidate, blockers, blocker_ks=[0, 1]),
        allow_multiexon_blockers=True,
    )
    assert final_ids(result) == {"left_multi", "right_multi"}
    assert result["r1_replacement_rows"] == []


def test_r1me_shared_splice_conflict_with_corroborated_model_is_nonreplaceable() -> None:
    candidate, single, _ = standard_models()
    splice_conflict = model("splice_conflict", "D", ((50, 200), (300, 450)))
    blockers = [single, splice_conflict]
    result, _ = run_r1(
        [candidate], blockers, rows_for(candidate, blockers, blocker_ks=[0, 1]),
        allow_multiexon_blockers=True,
    )
    assert final_ids(result) == {"blocker1", "splice_conflict"}
    assert result["backfill_rejected_rows"][0]["r1_replacement_reason"] == (  # type: ignore[index]
        "blocking_model_corroboration_count_ge_1"
    )


def test_r1me_suppressed_candidate_is_not_in_eligible_backfill_population() -> None:
    candidate, _blocker, _ = standard_models()
    blocker = model("multi_blocker", "A", ((120, 140), (160, 180)))
    result, _ = run_r1(
        [candidate], [blocker], rows_for(candidate, [blocker]),
        allow_multiexon_blockers=True, suppressed=[candidate],
    )
    assert final_ids(result) == {"multi_blocker"}
    assert result["backfill_candidate_rows"] == []


def test_r1me_is_deterministic_under_input_iteration_order() -> None:
    candidate, _blocker, _ = standard_models()
    left = model("left_multi", "A", ((120, 140), (160, 180)))
    right = model("right_multi", "B", ((320, 340), (360, 380)))
    blockers = [left, right]
    rows = rows_for(candidate, blockers)
    first, _ = run_r1(
        [candidate], blockers, rows, allow_multiexon_blockers=True,
    )
    second, _ = run_r1(
        [candidate], blockers, rows, allow_multiexon_blockers=True, reverse_inputs=True,
    )
    assert first["r1_replacement_rows"] == second["r1_replacement_rows"]
    assert final_ids(first) == final_ids(second)


def test_r1me_provenance_and_summary_distinguish_exon_classes(tmp_path: Path) -> None:
    candidate, single, _ = standard_models()
    multi = model("right_multi", "B", ((320, 340), (360, 380)))
    blockers = [single, multi]
    result, _ = run_r1(
        [candidate], blockers, rows_for(candidate, blockers),
        allow_multiexon_blockers=True, weight=0.5,
    )
    replacements = tmp_path / "replacements.tsv"
    summary = tmp_path / "summary.tsv"
    write_r1_replacement_outputs(result, replacements_tsv=replacements, summary_tsv=summary)
    replacement_text = replacements.read_text()
    summary_text = summary.read_text()
    assert "candidate_cds_exon_count" in replacement_text
    assert "replaced_model_cds_exon_counts" in replacement_text
    assert "single_exon,multi_exon" in replacement_text
    assert "extended_multiexon_blocker_eligibility" in replacement_text
    assert "multi_exon_blockers_replaced\t1" in summary_text
    assert "single_exon_blockers_replaced\t1" in summary_text


def test_r1me_disabled_reproduces_original_multiexon_rejection_and_schema(tmp_path: Path) -> None:
    candidate, _blocker, _ = standard_models()
    blocker = model("multi_blocker", "A", ((120, 140), (160, 180)))
    result, _ = run_r1([candidate], [blocker], rows_for(candidate, [blocker]))
    assert final_ids(result) == {"multi_blocker"}
    assert result["backfill_rejected_rows"][0]["r1_replacement_reason"] == "blocking_model_not_single_exon"  # type: ignore[index]

    original_candidate, original_blocker, _ = standard_models()
    original, _ = run_r1(
        [original_candidate], [original_blocker], rows_for(original_candidate, [original_blocker])
    )
    replacements = tmp_path / "replacements.tsv"
    write_r1_replacement_outputs(original, replacements_tsv=replacements, summary_tsv=None)
    header = replacements.read_text().splitlines()[0]
    assert "candidate_cds_exon_count" not in header
    assert "r1_eligibility_class" not in header
