from moselecta.alternative_splice_reporting import (
    AlternativeSpliceConsoleSummary,
    alternative_splice_console_lines,
)


def _summary(*, fraction: float | None, effective_minimum: int) -> AlternativeSpliceConsoleSummary:
    return AlternativeSpliceConsoleSummary(
        preset="intermediate",
        support_policy="exact-plus-mosaic",
        evidence_requirement="none",
        prediction_source_count=5,
        min_model_support_sources=2,
        min_model_support_source_fraction=fraction,
        effective_min_model_support_sources=effective_minimum,
        max_alternatives_per_gene=50,
        alternatives_emitted=True,
        uniquely_owned_candidates=0,
        fixed_rule_eligible_candidates=0,
        support_policy_eligible_candidates=0,
        external_evidence_requirement_failures=0,
        prethreshold_eligible_candidates=0,
        numerical_threshold_survivors=0,
        accepted_alternative_transcripts=0,
        genes_with_alternative_transcripts=0,
        maximum_accepted_alternatives_in_one_gene=0,
        final_annotation_genes=0,
        final_annotation_transcripts=0,
        final_annotation_primary_transcripts=0,
        final_annotation_alternative_transcripts=0,
    )


def test_console_reports_absolute_model_source_requirement() -> None:
    lines = alternative_splice_console_lines(
        _summary(fraction=None, effective_minimum=2)
    )
    assert "Prediction-model support requirement: >=2 of 5 sources" in lines


def test_console_reports_fractional_model_source_requirement() -> None:
    lines = alternative_splice_console_lines(
        _summary(fraction=0.5, effective_minimum=3)
    )
    assert (
        "Prediction-model support requirement: >=3 of 5 sources "
        "(configured fraction >=0.500)"
    ) in lines
