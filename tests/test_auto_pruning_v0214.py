from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.auto_pruning import (
    AutoDecision,
    ScenarioSummary,
    add_incremental_metrics,
    assess_flood_adjusted_consensus,
    automatic_diagnostics_rows,
    calculate_source_consensus_metrics,
    calculate_source_input_count_metrics,
    choose_auto_pruning_level,
    classify_evidence_adequacy,
    consensus_metrics_from_distribution,
    explain_auto_pruning_decision,
    flood_adjusted_consensus_log_lines,
    strong_deep_consensus_incidence_threshold,
    validate_scenarios,
)
from moselecta.models import CodingModel, Sublocus


def _scenario(
    level: str,
    starting: int,
    removed: set[str],
    source_by_id: dict[str, str] | None = None,
) -> ScenarioSummary:
    all_ids = {f"m{i}" for i in range(starting)}
    retained = all_ids - removed
    source_by_id = source_by_id or {model_id: "A" for model_id in all_ids}
    removal_by_source: dict[str, int] = {}
    for model_id in removed:
        source = source_by_id[model_id]
        removal_by_source[source] = removal_by_source.get(source, 0) + 1
    return ScenarioSummary(
        level=level,
        starting_models=starting,
        candidate_models=len(removed),
        protected_models=0,
        removed_models=len(removed),
        retained_models=len(retained),
        removed_fraction=len(removed) / starting,
        affected_subloci=len(removed),
        retained_subloci=len(retained),
        removed_ids=removed,
        retained_ids=retained,
        primary_reason_by_id={item: "singleton_sublocus" for item in removed},
        all_triggered_by_id={item: ("singleton_sublocus",) for item in removed},
        source_by_id=source_by_id,
        removal_by_source=removal_by_source,
    )


def _scenarios(*, flood_increment: bool = False) -> dict[str, ScenarioSummary]:
    source_by_id = {f"m{i}": "A" for i in range(1000)}
    if flood_increment:
        source_by_id.update({f"m{i}": "Flood" for i in range(100, 150)})
    scenarios = {
        "medium": _scenario("medium", 1000, {f"m{i}" for i in range(100)}, source_by_id),
        "high": _scenario("high", 1000, {f"m{i}" for i in range(150)}, source_by_id),
        "very-high": _scenario("very-high", 1000, {f"m{i}" for i in range(170)}, source_by_id),
    }
    add_incremental_metrics(scenarios)
    return scenarios


def _flooded_consensus():
    return consensus_metrics_from_distribution(
        {1: 2060, 2: 200, 3: 300},
        raw_source_count=3,
        per_source_singleton_count={"A": 30, "B": 30, "Flood": 2000},
        per_source_total_subloci={"A": 530, "B": 530, "Flood": 2300},
    )


@pytest.mark.parametrize(
    ("source_count", "below", "at", "above", "threshold"),
    (
        (3, {1: 70, 2: 100, 3: 109}, {1: 70, 2: 100, 3: 110}, {1: 70, 2: 100, 3: 111}, 0.55),
        (4, {1: 40, 2: 80, 3: 99}, {1: 40, 2: 80, 3: 100}, {1: 40, 2: 80, 3: 101}, 0.60),
        (5, {1: 40, 2: 80, 3: 99}, {1: 40, 2: 80, 3: 100}, {1: 40, 2: 80, 3: 101}, 0.60),
    ),
)
def test_deep_consensus_thresholds_immediately_below_at_and_above(
    source_count: int,
    below: dict[int, int],
    at: dict[int, int],
    above: dict[int, int],
    threshold: float,
) -> None:
    assert strong_deep_consensus_incidence_threshold(source_count) == threshold
    assert consensus_metrics_from_distribution(
        below, raw_source_count=source_count
    ).backbone_class == "intermediate"
    assert consensus_metrics_from_distribution(
        at, raw_source_count=source_count
    ).backbone_class == "strong"
    assert consensus_metrics_from_distribution(
        above, raw_source_count=source_count
    ).backbone_class == "strong"


def test_other_strong_criteria_and_intermediate_thresholds_still_bind() -> None:
    too_few_multisource = consensus_metrics_from_distribution(
        {1: 1, 3: 99}, raw_source_count=3
    )
    assert too_few_multisource.deep_consensus_incidence_fraction > 0.55
    assert too_few_multisource.backbone_class == "weak_or_discordant"
    intermediate = consensus_metrics_from_distribution(
        {1: 76, 2: 100, 3: 100}, raw_source_count=3
    )
    assert intermediate.deep_consensus_incidence_fraction == pytest.approx(0.5208333333)
    assert intermediate.backbone_class == "intermediate"


def test_flood_discount_recalculates_denominators_and_preserves_source_count() -> None:
    raw = _flooded_consensus()
    assessment = assess_flood_adjusted_consensus(raw)
    assert raw.backbone_class == "weak_or_discordant"
    assert raw.strong_source_excess_sources == ("Flood",)
    assert raw.per_source_singleton_share["Flood"] >= 0.60
    assert raw.per_source_singleton_rate["Flood"] >= 0.20
    assert raw.per_source_singleton_rate_ratio["Flood"] >= 2.0
    assert assessment.qualifying_excess_source == "Flood"
    assert assessment.discounted_source_singleton_subloci == 2000
    adjusted = assessment.adjusted
    assert adjusted is not None
    assert adjusted.source_count_distribution == {1: 60, 2: 200, 3: 300}
    assert adjusted.total_subloci == raw.total_subloci - 2000
    assert adjusted.total_source_incidences == raw.total_source_incidences - 2000
    assert adjusted.raw_source_count == adjusted.effective_source_count == 3
    assert adjusted.multisource_subloci == raw.multisource_subloci
    assert adjusted.backbone_class == "strong"


def test_flood_discount_is_deterministic_and_source_order_independent() -> None:
    first = _flooded_consensus()
    second = consensus_metrics_from_distribution(
        {3: 300, 1: 2060, 2: 200},
        raw_source_count=3,
        per_source_singleton_count={"Flood": 2000, "B": 30, "A": 30},
        per_source_total_subloci={"Flood": 2300, "B": 530, "A": 530},
    )
    first_adjusted = assess_flood_adjusted_consensus(first).adjusted
    second_adjusted = assess_flood_adjusted_consensus(second).adjusted
    assert first_adjusted is not None and second_adjusted is not None
    assert first_adjusted.source_count_distribution == second_adjusted.source_count_distribution
    assert first_adjusted.backbone_class == second_adjusted.backbone_class


def test_distributed_disagreement_has_no_qualifying_excess_source() -> None:
    metrics = consensus_metrics_from_distribution(
        {1: 900, 2: 300, 3: 200},
        raw_source_count=3,
        per_source_singleton_count={"A": 300, "B": 300, "C": 300},
        per_source_total_subloci={"A": 800, "B": 800, "C": 800},
    )
    assessment = assess_flood_adjusted_consensus(metrics)
    assert metrics.backbone_class == "weak_or_discordant"
    assert metrics.strong_source_excess_sources == ()
    assert assessment.adjusted is None


def test_two_substantial_singleton_sources_do_not_qualify_for_adjustment() -> None:
    metrics = consensus_metrics_from_distribution(
        {1: 1000, 2: 200, 3: 200},
        raw_source_count=3,
        per_source_singleton_count={"A": 500, "B": 500, "C": 0},
        per_source_total_subloci={"A": 700, "B": 700, "C": 400},
    )
    # Neither reaches the required 0.60 share; this is deliberately
    # distributed disagreement rather than one unambiguous flooding source.
    assert metrics.strong_source_excess_sources == ()
    assert assess_flood_adjusted_consensus(metrics).adjusted is None


def test_incomplete_source_pattern_is_not_promoted_by_another_sources_isolates() -> None:
    metrics = consensus_metrics_from_distribution(
        {1: 10, 2: 1000, 3: 5},
        raw_source_count=3,
        per_source_singleton_count={"A": 5, "B": 5, "Truncated": 0},
        per_source_total_subloci={"A": 1005, "B": 1005, "Truncated": 5},
    )
    assert metrics.backbone_class == "weak_or_discordant"
    assert assess_flood_adjusted_consensus(metrics).adjusted is None


def test_truncated_source_negative_control_adjusts_but_remains_weak() -> None:
    raw = consensus_metrics_from_distribution(
        {1: 1050, 2: 900, 3: 10},
        raw_source_count=3,
        per_source_singleton_count={"A": 45, "Flood": 1000, "Truncated": 5},
        per_source_total_subloci={"A": 955, "Flood": 1910, "Truncated": 15},
    )
    original_distribution = dict(raw.source_count_distribution)
    original_singletons = dict(raw.per_source_singleton_count)
    original_total_subloci = raw.total_subloci
    assessment = assess_flood_adjusted_consensus(raw)

    assert raw.backbone_class == "weak_or_discordant"
    assert raw.strong_source_excess_sources == ("Flood",)
    assert assessment.qualifying_excess_source == "Flood"
    assert assessment.discounted_source_singleton_subloci == 1000
    assert assessment.unadjusted is raw
    assert assessment.adjusted is not None
    assert assessment.adjusted.source_count_distribution == {1: 50, 2: 900, 3: 10}
    assert assessment.adjusted.total_subloci == 960
    assert assessment.adjusted.total_source_incidences == 1880
    assert assessment.adjusted.backbone_class == "weak_or_discordant"

    decision = choose_auto_pruning_level(
        "auto",
        scenarios=_scenarios(),
        validation=validate_scenarios(_scenarios()),
        consensus=raw,
        adjusted_consensus=assessment.adjusted,
        adequacy_pass_count=0,
        selected_singleton_model_fraction=0.5,
    )
    assert decision.unconstrained_level == decision.effective_level == "medium"
    assert decision.used_flood_adjustment is True
    assert raw.source_count_distribution == original_distribution
    assert raw.per_source_singleton_count == original_singletons
    assert raw.total_subloci == original_total_subloci


def test_disabled_and_degenerate_adjustments_are_conservative() -> None:
    raw = _flooded_consensus()
    disabled = assess_flood_adjusted_consensus(raw, enabled=False)
    assert disabled.qualifying_excess_source == "Flood"
    assert disabled.adjusted is None
    empty = consensus_metrics_from_distribution({}, raw_source_count=3)
    assert empty.backbone_class == "weak_or_discordant"
    assert assess_flood_adjusted_consensus(empty).adjusted is None


def _model(model_id: str, source: str, start: int) -> CodingModel:
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"g_{model_id}",
        input_transcript_id=model_id,
        output_gene_id=f"out_g_{model_id}",
        output_transcript_id=model_id,
        chrom="chr1",
        strand="+",
        cds_parts=[(start, start + 30)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def test_normal_and_distribution_fixture_interfaces_have_equivalent_semantics() -> None:
    distribution = {1: 70, 2: 100, 3: 110}
    subloci: list[Sublocus] = []
    index = 0
    sources = ("A", "B", "C")
    for source_count, number in distribution.items():
        for _ in range(number):
            models = [
                _model(f"m{index}_{source}", source, index * 40)
                for source in sources[:source_count]
            ]
            subloci.append(
                Sublocus(
                    f"sub{index}", "sl1", "chr1", "+", index * 40,
                    index * 40 + 30, models,
                )
            )
            index += 1
    normal = calculate_source_consensus_metrics(subloci)
    fixture = consensus_metrics_from_distribution(
        distribution,
        raw_source_count=3,
        per_source_singleton_count=normal.per_source_singleton_count,
        per_source_total_subloci=normal.per_source_total_subloci,
    )
    assert normal.source_count_distribution == fixture.source_count_distribution
    assert normal.total_source_incidences == fixture.total_source_incidences
    assert normal.deep_consensus_incidence_fraction == fixture.deep_consensus_incidence_fraction
    assert normal.backbone_class == fixture.backbone_class == "strong"


def test_auto_uses_adjusted_backbone_but_conservative_uses_unadjusted() -> None:
    scenarios = _scenarios()
    raw = _flooded_consensus()
    adjusted = assess_flood_adjusted_consensus(raw).adjusted
    validation = validate_scenarios(scenarios)
    standard = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validation, consensus=raw,
        adjusted_consensus=adjusted, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.5,
    )
    conservative = choose_auto_pruning_level(
        "auto-conservative", scenarios=scenarios, validation=validation,
        consensus=raw, adjusted_consensus=adjusted, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.5, max_level="high",
    )
    assert standard.effective_level == "high"
    assert standard.used_flood_adjustment is True
    assert standard.backbone_class_used == "strong"
    assert conservative.effective_level == "medium"
    assert conservative.used_flood_adjustment is False
    assert conservative.backbone_class_used == "weak_or_discordant"


def test_source_excess_disabled_keeps_standard_auto_on_raw_backbone() -> None:
    scenarios = _scenarios()
    raw = _flooded_consensus()
    adjusted = assess_flood_adjusted_consensus(raw).adjusted
    decision = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=raw, adjusted_consensus=adjusted, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.5, use_source_excess=False,
    )
    assert decision.effective_level == "medium"
    assert decision.used_flood_adjustment is False


def test_revised_three_source_threshold_promotes_strong_but_not_0521_case() -> None:
    scenarios = _scenarios()
    validation = validate_scenarios(scenarios)
    strong = consensus_metrics_from_distribution(
        {1: 70, 2: 100, 3: 110}, raw_source_count=3
    )
    conservative_case = consensus_metrics_from_distribution(
        {1: 76, 2: 100, 3: 100}, raw_source_count=3
    )
    assert choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validation, consensus=strong,
        adequacy_pass_count=0, selected_singleton_model_fraction=0.2,
    ).effective_level == "high"
    assert choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validation,
        consensus=conservative_case, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2,
    ).effective_level == "medium"


def test_existing_intermediate_source_excess_escalation_is_preserved() -> None:
    raw = consensus_metrics_from_distribution(
        {1: 1060, 2: 200, 3: 300},
        raw_source_count=3,
        per_source_singleton_count={"A": 30, "B": 30, "Flood": 1000},
        per_source_total_subloci={"A": 530, "B": 530, "Flood": 1300},
    )
    assert raw.backbone_class == "intermediate"
    scenarios = _scenarios(flood_increment=True)
    decision = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=raw, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2,
    )
    assert decision.effective_level == "high"
    assert "INTERMEDIATE_BACKBONE_SINGLE_STRONG_SOURCE_EXCESS" in decision.reason_codes


def test_adjustment_does_not_change_adequate_evidence_path_or_user_bounds() -> None:
    scenarios = _scenarios()
    raw = _flooded_consensus()
    adjusted = assess_flood_adjusted_consensus(raw).adjusted
    common = dict(
        mode="auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=raw, adequacy_pass_count=1,
        selected_singleton_model_fraction=0.2, min_level="low", max_level="medium",
    )
    without = choose_auto_pruning_level(**common)
    with_adjustment = choose_auto_pruning_level(
        **common, adjusted_consensus=adjusted
    )
    assert without == with_adjustment
    assert with_adjustment.unconstrained_level in {"high", "very-high"}
    assert with_adjustment.effective_level == "medium"
    assert with_adjustment.used_flood_adjustment is False
    explanation = explain_auto_pruning_decision(
        with_adjustment,
        consensus=raw,
        flood_assessment=assess_flood_adjusted_consensus(raw),
    )
    assert "adequate-evidence pathway" in explanation
    assert "does not alter this pathway" in explanation


def test_input_count_warnings_use_run_median_without_collateral_lycorma_flags() -> None:
    metrics = calculate_source_input_count_metrics(
        {"Helixer": 18_722, "Annevo": 19_022, "Tiberius": 89_580}
    )
    assert {item.run_median_input_count for item in metrics.values()} == {19_022.0}
    assert metrics["Tiberius"].ratio_to_run_median == pytest.approx(4.709284)
    assert metrics["Tiberius"].warning == "INPUT_CODING_MODEL_COUNT_HIGH"
    assert metrics["Helixer"].warning == "."
    assert metrics["Annevo"].warning == "."
    # The historical leave-one-out values remain available but no longer drive
    # warning codes.
    assert metrics["Helixer"].other_source_median_count == 54_301.0
    assert metrics["Helixer"].ratio_to_other_source_median < 0.5


def test_input_count_warnings_detect_truncated_melospiza_without_false_high() -> None:
    metrics = calculate_source_input_count_metrics(
        {"Helixer": 1_643, "Annevo": 19_016, "Tiberius": 28_236}
    )
    assert {item.run_median_input_count for item in metrics.values()} == {19_016.0}
    assert metrics["Helixer"].ratio_to_run_median == pytest.approx(0.0864, rel=1e-3)
    assert metrics["Helixer"].warning == "INPUT_CODING_MODEL_COUNT_LOW"
    assert metrics["Annevo"].warning == "."
    assert metrics["Tiberius"].ratio_to_run_median == pytest.approx(1.485, rel=1e-3)
    assert metrics["Tiberius"].warning == "."


@pytest.mark.parametrize(
    ("counts", "expected_median"),
    (
        ({"A": 100, "B": 110, "C": 120}, 110.0),
        ({"A": 100, "B": 110, "C": 120, "D": 130}, 115.0),
        ({"A": 80, "B": 90, "C": 100, "D": 110, "E": 120}, 100.0),
        ({"A": 100, "B": 100, "C": 100}, 100.0),
    ),
)
def test_input_count_run_median_normal_even_odd_and_equal_cases(
    counts: dict[str, int], expected_median: float
) -> None:
    metrics = calculate_source_input_count_metrics(counts)
    assert {item.run_median_input_count for item in metrics.values()} == {
        expected_median
    }
    assert {item.warning for item in metrics.values()} == {"."}


def test_input_count_warning_boundaries_are_inclusive_and_borderlines_are_clear() -> None:
    exact = calculate_source_input_count_metrics({"Low": 50, "Middle": 100, "High": 200})
    assert exact["Low"].ratio_to_run_median == 0.5
    assert exact["Low"].warning == "INPUT_CODING_MODEL_COUNT_LOW"
    assert exact["High"].ratio_to_run_median == 2.0
    assert exact["High"].warning == "INPUT_CODING_MODEL_COUNT_HIGH"
    borderline = calculate_source_input_count_metrics(
        {"Low": 51, "Middle": 100, "High": 199}
    )
    assert {item.warning for item in borderline.values()} == {"."}


def test_input_count_warning_zero_cases_are_safe_and_deterministic() -> None:
    defined = calculate_source_input_count_metrics({"Zero": 0, "B": 100, "C": 100})
    assert defined["Zero"].ratio_to_run_median == 0.0
    assert defined["Zero"].warning == "INPUT_CODING_MODEL_COUNT_LOW"
    undefined = calculate_source_input_count_metrics({"A": 10, "B": 0, "C": 0})
    assert {item.run_median_input_count for item in undefined.values()} == {0.0}
    assert all(item.ratio_to_run_median is None for item in undefined.values())
    assert {item.warning for item in undefined.values()} == {"."}
    assert undefined["A"].ratio_to_other_source_median == pytest.approx(float("inf"))


def test_input_count_warning_results_are_source_order_independent() -> None:
    forward = calculate_source_input_count_metrics(
        {"Helixer": 18_722, "Annevo": 19_022, "Tiberius": 89_580}
    )
    reverse = calculate_source_input_count_metrics(
        {"Tiberius": 89_580, "Annevo": 19_022, "Helixer": 18_722}
    )
    assert forward == reverse


def test_input_count_warnings_do_not_affect_consensus_or_decisions() -> None:
    scenarios = _scenarios()
    consensus = _flooded_consensus()
    assessment_before = assess_flood_adjusted_consensus(consensus)
    before = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=consensus, adequacy_pass_count=0,
        adjusted_consensus=assessment_before.adjusted,
        selected_singleton_model_fraction=0.2,
    )
    calculate_source_input_count_metrics(
        {"Helixer": 18_722, "Annevo": 19_022, "Tiberius": 89_580}
    )
    assessment_after = assess_flood_adjusted_consensus(consensus)
    after = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=consensus, adjusted_consensus=assessment_after.adjusted,
        adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2,
    )
    assert consensus.strong_source_excess_sources == ("Flood",)
    assert assessment_before == assessment_after
    assert before == after


def test_diagnostics_preserve_raw_fields_and_append_adjusted_decision_details() -> None:
    scenarios = _scenarios(flood_increment=True)
    raw = _flooded_consensus()
    assessment = assess_flood_adjusted_consensus(raw)
    validation = validate_scenarios(scenarios)
    standard = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validation, consensus=raw,
        adjusted_consensus=assessment.adjusted, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2,
    )
    conservative = choose_auto_pruning_level(
        "auto-conservative", scenarios=scenarios, validation=validation,
        consensus=raw, adjusted_consensus=assessment.adjusted,
        adequacy_pass_count=0, selected_singleton_model_fraction=0.2,
        max_level="high",
    )
    absent = classify_evidence_adequacy(
        "protein", supplied=False, first_value=0, second_value=0,
        first_threshold=0.75, second_threshold=0.5,
    )
    transcript = classify_evidence_adequacy(
        "transcript", supplied=False, first_value=0, second_value=0,
        first_threshold=0.75, second_threshold=0.5,
    )
    rows = automatic_diagnostics_rows(
        protein=absent, transcript=transcript, consensus=raw,
        selected_singleton_model_fraction=0.2, validation=validation,
        standard_decision=standard, conservative_decision=conservative,
        flood_assessment=assessment,
        input_count_metrics=calculate_source_input_count_metrics(
            {"A": 100, "B": 100, "Flood": 400}
        ),
        scenarios=scenarios,
    )
    values = {(row["section"], row["metric"]): row["value"] for row in rows}
    assert values[("source_consensus", "backbone_class")] == "weak_or_discordant"
    assert values[("source_consensus_v0.2.14", "unadjusted_backbone_class")] == "weak_or_discordant"
    assert values[("source_consensus_adjusted", "backbone_class")] == "strong"
    assert values[("source_consensus_adjustment", "discounted_source_singleton_subloci")] == 2000
    assert values[("decision", "auto.backbone_class_used")] == "strong"
    assert values[("decision", "auto-conservative.backbone_class_used")] == "weak_or_discordant"
    assert values[("alternative_pruning", "high_vs_medium_additional_models_removed")] == 50
    assert values[("input_count_summary", "run_median_input_count")] == "100.000000"
    assert values[("per_source_input", "Flood.ratio_to_run_median")] == "4.000000"
    assert values[("per_source_input", "Flood.other_source_median_count")] == "100.000000"
    assert values[("per_source_input", "Flood.ratio_to_other_source_median")] == "4.000000"
    assert values[("per_source_input", "Flood.warning_code")] == "INPUT_CODING_MODEL_COUNT_HIGH"
    explanation = values[("decision", "auto.explanation")]
    assert "discounted only for backbone assessment" in explanation
    assert "No prediction models were removed" in explanation


def test_conservative_medium_explanation_identifies_threshold_failure() -> None:
    metrics = consensus_metrics_from_distribution(
        {1: 76, 2: 100, 3: 100}, raw_source_count=3
    )
    decision = AutoDecision(
        "auto-conservative", "medium", "medium", "moderate", (), (),
        "medium", "high", metrics.backbone_class, False,
    )
    explanation = explain_auto_pruning_decision(decision, consensus=metrics)
    assert "deep-consensus incidence 0.521" in explanation
    assert ">= 0.550" in explanation
    assert "No single source" in explanation


def test_explanation_distinguishes_recommendation_from_user_maximum() -> None:
    scenarios = _scenarios()
    consensus = consensus_metrics_from_distribution(
        {1: 70, 2: 100, 3: 110}, raw_source_count=3
    )
    decision = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=consensus, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2,
        min_level="medium", max_level="medium",
    )
    explanation = explain_auto_pruning_decision(decision, consensus=consensus)
    assert decision.unconstrained_level == "high"
    assert decision.effective_level == "medium"
    assert "recommended high" in explanation
    assert "maximum pruning level is medium" in explanation
    assert "effective level is medium" in explanation


def test_explanation_distinguishes_recommendation_from_user_minimum() -> None:
    scenarios = _scenarios()
    consensus = consensus_metrics_from_distribution(
        {1: 76, 2: 100, 3: 100}, raw_source_count=3
    )
    decision = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=consensus, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2,
        min_level="high", max_level="very-high",
    )
    explanation = explain_auto_pruning_decision(decision, consensus=consensus)
    assert decision.unconstrained_level == "medium"
    assert decision.effective_level == "high"
    assert "recommended medium" in explanation
    assert "minimum pruning level is high" in explanation
    assert "effective level is high" in explanation


def test_explanation_has_no_constraint_text_when_level_is_unchanged() -> None:
    scenarios = _scenarios()
    consensus = consensus_metrics_from_distribution(
        {1: 70, 2: 100, 3: 110}, raw_source_count=3
    )
    decision = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=consensus, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2,
    )
    explanation = explain_auto_pruning_decision(decision, consensus=consensus)
    assert decision.unconstrained_level == decision.effective_level == "high"
    assert "recommended high" in explanation
    assert "user-specified" not in explanation


def test_conservative_explanation_uses_raw_weak_backbone_despite_adjusted_strong() -> None:
    scenarios = _scenarios()
    raw = _flooded_consensus()
    assessment = assess_flood_adjusted_consensus(raw)
    decision = choose_auto_pruning_level(
        "auto-conservative", scenarios=scenarios,
        validation=validate_scenarios(scenarios), consensus=raw,
        adjusted_consensus=assessment.adjusted, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2, max_level="high",
    )
    explanation = explain_auto_pruning_decision(
        decision, consensus=raw, flood_assessment=assessment
    )
    assert decision.unconstrained_level == decision.effective_level == "medium"
    assert "Flood-adjusted consensus is strong" in explanation
    assert "uses the unadjusted weak_or_discordant backbone" in explanation
    assert "Medium pruning was recommended" in explanation


def test_conservative_raw_strong_explanation_still_recommends_high() -> None:
    scenarios = _scenarios()
    raw = consensus_metrics_from_distribution(
        {1: 160, 2: 300, 3: 500},
        raw_source_count=3,
        per_source_singleton_count={"A": 30, "B": 30, "Flood": 100},
        per_source_total_subloci={"A": 905, "B": 905, "Flood": 450},
    )
    assessment = assess_flood_adjusted_consensus(raw)
    assert raw.backbone_class == "strong"
    assert assessment.adjusted is not None
    assert assessment.adjusted.backbone_class == "strong"
    decision = choose_auto_pruning_level(
        "auto-conservative", scenarios=scenarios,
        validation=validate_scenarios(scenarios), consensus=raw,
        adjusted_consensus=assessment.adjusted, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2, max_level="high",
    )
    explanation = explain_auto_pruning_decision(
        decision, consensus=raw, flood_assessment=assessment
    )
    assert decision.unconstrained_level == decision.effective_level == "high"
    assert "uses the unadjusted strong backbone" in explanation
    assert "High pruning was recommended" in explanation


def test_adjustment_log_reports_metrics_thresholds_basis_and_scope() -> None:
    scenarios = _scenarios()
    assessment = assess_flood_adjusted_consensus(_flooded_consensus())
    assert assessment.adjusted is not None
    standard = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=assessment.unadjusted,
        adjusted_consensus=assessment.adjusted, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2,
    )
    conservative = choose_auto_pruning_level(
        "auto-conservative", scenarios=scenarios,
        validation=validate_scenarios(scenarios),
        consensus=assessment.unadjusted,
        adjusted_consensus=assessment.adjusted, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2, max_level="high",
    )
    lines = flood_adjusted_consensus_log_lines(
        assessment, requested_mode="auto", standard_decision=standard,
        conservative_decision=conservative,
    )
    report = "\n".join(lines)
    assert "singleton incidence" in report
    assert "corroboration" in report
    assert "deep-consensus incidence" in report
    assert "deep-consensus sublocus fraction" in report
    assert "<= 0.3500" in report
    assert ">= 0.5500" in report
    assert "weak_or_discordant -> strong" in report
    assert "excessive_source=Flood" in report
    assert "discounted_source_singleton_subloci=2,000" in report
    assert "standard auto uses the adjusted backbone" in report
    assert "class_changed=true" in report
    assert "no subloci or prediction models were removed" in report

    evidence_decision = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=assessment.unadjusted,
        adjusted_consensus=assessment.adjusted, adequacy_pass_count=1,
        selected_singleton_model_fraction=0.2,
    )
    evidence_report = "\n".join(flood_adjusted_consensus_log_lines(
        assessment, requested_mode="auto",
        standard_decision=evidence_decision,
        conservative_decision=conservative,
    ))
    assert evidence_decision.used_flood_adjustment is False
    assert "adjusted assessment did not determine" in evidence_report


def test_adjustment_log_shows_weak_to_weak_and_conservative_basis() -> None:
    scenarios = _scenarios()
    raw = consensus_metrics_from_distribution(
        {1: 1050, 2: 900, 3: 10},
        raw_source_count=3,
        per_source_singleton_count={"A": 45, "Flood": 1000, "Truncated": 5},
        per_source_total_subloci={"A": 955, "Flood": 1910, "Truncated": 15},
    )
    assessment = assess_flood_adjusted_consensus(raw)
    assert assessment.adjusted is not None
    standard = choose_auto_pruning_level(
        "auto", scenarios=scenarios, validation=validate_scenarios(scenarios),
        consensus=raw, adjusted_consensus=assessment.adjusted,
        adequacy_pass_count=0, selected_singleton_model_fraction=0.2,
    )
    conservative = choose_auto_pruning_level(
        "auto-conservative", scenarios=scenarios,
        validation=validate_scenarios(scenarios), consensus=raw,
        adjusted_consensus=assessment.adjusted, adequacy_pass_count=0,
        selected_singleton_model_fraction=0.2, max_level="high",
    )
    report = "\n".join(flood_adjusted_consensus_log_lines(
        assessment, requested_mode="auto-conservative",
        standard_decision=standard, conservative_decision=conservative,
    ))
    assert "weak_or_discordant -> weak_or_discordant" in report
    assert "auto-conservative uses the unadjusted backbone" in report
    assert "adjusted class cannot promote" in report
