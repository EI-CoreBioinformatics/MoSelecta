from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.auto_pruning import (
    ScenarioSummary,
    add_incremental_metrics,
    choose_auto_pruning_level,
    classify_evidence_adequacy,
    consensus_metrics_from_distribution,
    validate_scenarios,
    write_scenario_summary,
)


def _scenario(level: str, starting: int, removed: set[str]) -> ScenarioSummary:
    all_ids = {f"m{i}" for i in range(starting)}
    retained = all_ids - removed
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
        primary_reason_by_id={model_id: "singleton_sublocus" for model_id in removed},
        all_triggered_by_id={model_id: ("singleton_sublocus",) for model_id in removed},
        removal_by_source={"A": len(removed)},
        source_selected_counts={"A": starting},
        source_removed_fractions={"A": len(removed) / starting},
        source_share_all_removals={"A": 1.0},
    )


def _scenarios(
    starting: int = 1000,
    medium_removed: int = 100,
    high_removed: int = 150,
    very_high_removed: int = 170,
) -> dict[str, ScenarioSummary]:
    scenarios = {
        "medium": _scenario("medium", starting, {f"m{i}" for i in range(medium_removed)}),
        "high": _scenario("high", starting, {f"m{i}" for i in range(high_removed)}),
        "very-high": _scenario("very-high", starting, {f"m{i}" for i in range(very_high_removed)}),
    }
    add_incremental_metrics(scenarios)
    return scenarios


def test_evidence_adequacy_states_and_marginal_band() -> None:
    assert classify_evidence_adequacy(
        "protein", supplied=False, first_value=0, second_value=0,
        first_threshold=0.75, second_threshold=0.50,
    ).state == "absent"
    assert classify_evidence_adequacy(
        "protein", supplied=True, first_value=0.80, second_value=0.50,
        first_threshold=0.75, second_threshold=0.50,
    ).state == "pass"
    marginal = classify_evidence_adequacy(
        "protein", supplied=True, first_value=0.72, second_value=0.60,
        first_threshold=0.75, second_threshold=0.50,
    )
    assert marginal.state == "marginal"
    assert marginal.limiting_margin == pytest.approx(-0.03)
    assert classify_evidence_adequacy(
        "protein", supplied=True, first_value=0.69, second_value=0.60,
        first_threshold=0.75, second_threshold=0.50,
    ).state == "fail"


def test_mock_discordant_three_source_distribution_is_weak() -> None:
    metrics = consensus_metrics_from_distribution({1: 4883, 2: 445, 3: 153})
    assert metrics.backbone_class == "weak_or_discordant"
    assert metrics.singleton_incidence_fraction == pytest.approx(4883 / 6232)
    assert metrics.corroboration_fraction == pytest.approx((6232 - 5481) / 6232)


HAT_AND_HATBG_DISTRIBUTIONS = (
    {1: 3717, 2: 2634, 3: 21983},
    {1: 8265, 2: 1797, 3: 1452, 4: 2608, 5: 20359},
    {1: 1755, 2: 2000, 3: 8180},
    {1: 6119, 2: 625, 3: 568, 4: 1969, 5: 7458},
    {1: 4930, 2: 5655, 3: 16468},
    {1: 33269, 2: 2584, 3: 2968, 4: 4533, 5: 14264},
    {1: 14753, 2: 7265, 3: 25081},
    {1: 26346, 2: 4064, 3: 4212, 4: 6604, 5: 19091},
    {1: 4883, 2: 4456, 3: 15364},
    {1: 31453, 2: 1987, 3: 2661, 4: 5033, 5: 11738},
)


@pytest.mark.parametrize("distribution", HAT_AND_HATBG_DISTRIBUTIONS)
def test_supplied_no_evidence_calibration_distributions_are_strong(distribution: dict[int, int]) -> None:
    metrics = consensus_metrics_from_distribution(distribution)
    assert metrics.backbone_class == "strong"
    scenarios = _scenarios()
    validation = validate_scenarios(scenarios)
    for mode in ("auto", "auto-conservative"):
        decision = choose_auto_pruning_level(
            mode,
            scenarios=scenarios,
            validation=validation,
            consensus=metrics,
            adequacy_pass_count=0,
            selected_singleton_model_fraction=0.5,
            max_level="very-high" if mode == "auto" else "high",
        )
        assert decision.unconstrained_level == "high"
        assert decision.effective_level == "high"
        assert "NO_ADEQUATE_EVIDENCE_AUTOMATIC_MAXIMUM_HIGH" in decision.reason_codes


def test_two_source_no_evidence_and_mock_distribution_both_choose_medium() -> None:
    scenarios = _scenarios()
    validation = validate_scenarios(scenarios)
    for metrics in (
        consensus_metrics_from_distribution({1: 100, 2: 900}, raw_source_count=2),
        consensus_metrics_from_distribution({1: 4883, 2: 445, 3: 153}),
    ):
        for mode in ("auto", "auto-conservative"):
            decision = choose_auto_pruning_level(
                mode,
                scenarios=scenarios,
                validation=validation,
                consensus=metrics,
                adequacy_pass_count=0,
                selected_singleton_model_fraction=0.8,
                max_level="very-high" if mode == "auto" else "high",
            )
            assert decision.effective_level == "medium"


def test_adequate_evidence_standard_auto_selects_very_high_inside_five_percent_cap() -> None:
    scenarios = _scenarios(starting=1000, medium_removed=100, high_removed=150, very_high_removed=170)
    decision = choose_auto_pruning_level(
        "auto",
        scenarios=scenarios,
        validation=validate_scenarios(scenarios),
        consensus=consensus_metrics_from_distribution({1: 50, 2: 100, 3: 500}),
        adequacy_pass_count=1,
        selected_singleton_model_fraction=0.20,
        max_level="very-high",
    )
    assert decision.effective_level == "very-high"
    assert "VERY_HIGH_INCREMENT_WITHIN_SAFETY_LIMIT" in decision.reason_codes


def test_adequate_evidence_standard_auto_stays_high_above_five_percent_cap() -> None:
    scenarios = _scenarios(starting=1000, medium_removed=100, high_removed=150, very_high_removed=250)
    decision = choose_auto_pruning_level(
        "auto",
        scenarios=scenarios,
        validation=validate_scenarios(scenarios),
        consensus=consensus_metrics_from_distribution({1: 50, 2: 100, 3: 500}),
        adequacy_pass_count=1,
        selected_singleton_model_fraction=0.20,
        max_level="very-high",
    )
    assert decision.effective_level == "high"
    assert "VERY_HIGH_EXTRA_REMOVAL_FRACTION_EXCEEDS_LIMIT" in decision.warning_codes


def test_conservative_uses_stricter_high_trigger_and_requires_explicit_very_high() -> None:
    scenarios = _scenarios(starting=1000, medium_removed=100, high_removed=125, very_high_removed=135)
    common = dict(
        scenarios=scenarios,
        validation=validate_scenarios(scenarios),
        consensus=consensus_metrics_from_distribution({1: 50, 2: 100, 3: 500}),
        adequacy_pass_count=1,
        selected_singleton_model_fraction=0.12,
    )
    decision = choose_auto_pruning_level(
        "auto-conservative", **common, max_level="high",
    )
    assert decision.effective_level == "medium"
    explicit = choose_auto_pruning_level(
        "auto-conservative",
        **{**common, "selected_singleton_model_fraction": 0.20},
        max_level="very-high",
        max_level_explicit=True,
    )
    assert explicit.effective_level == "very-high"


def test_validation_failure_falls_back_to_medium_low_confidence() -> None:
    scenarios = _scenarios()
    scenarios["high"].retained_ids.add("m0")
    decision = choose_auto_pruning_level(
        "auto",
        scenarios=scenarios,
        validation=validate_scenarios(scenarios),
        consensus=consensus_metrics_from_distribution({1: 50, 2: 100, 3: 500}),
        adequacy_pass_count=2,
        selected_singleton_model_fraction=0.5,
    )
    assert decision.effective_level == "medium"
    assert decision.confidence == "low"


def test_scenario_tsv_order_is_stable(tmp_path: Path) -> None:
    scenarios = _scenarios()
    out = tmp_path / "scenarios.tsv"
    write_scenario_summary(out, scenarios)
    lines = out.read_text().splitlines()
    assert [line.split("\t", 1)[0] for line in lines[1:]] == [
        "medium", "high", "very-high",
    ]
