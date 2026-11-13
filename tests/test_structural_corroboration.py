"""Focused synthetic tests for diagnostic run-aware corroboration."""

from pathlib import Path

import pytest

from moselecta.corroboration import (
    CorroborationComparison,
    compare_structural_corroboration,
    corroboration_normalizations,
    intron_agreement_requirement_met,
    make_structural_corroboration_row,
)
from moselecta.constants import MODEL_SELECTION_DEFAULT_WEIGHTS
from moselecta.models import CodingModel, Sublocus
from moselecta.scoring import (
    apply_model_selection_scores,
    model_selection_score_positive_scale_max,
    saturating_corroboration_score,
)


def model(
    model_id: str,
    source: str,
    parts: tuple[tuple[int, int], ...],
    *,
    strand: str = "+",
) -> CodingModel:
    item = CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"{model_id}.gene",
        input_transcript_id=f"{model_id}.input",
        output_gene_id=f"{model_id}.gene.out",
        output_transcript_id=model_id,
        chrom="chr1",
        strand=strand,
        cds_parts=list(parts),
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )
    item.superlocus_id = "super1"
    item.sublocus_id = "super1.sublocus001"
    return item


def sublocus(models: list[CodingModel]) -> Sublocus:
    return Sublocus(
        "super1.sublocus001",
        "super1",
        "chr1",
        models[0].strand,
        min(item.cds_start0 for item in models),
        max(item.cds_end0 for item in models),
        models,
        "single_exon" if models[0].is_single_exon else "multi_exon",
    )


def assert_symmetric(a: CodingModel, b: CodingModel, expected: bool) -> None:
    assert compare_structural_corroboration(a, b).corroborates is expected
    assert compare_structural_corroboration(b, a).corroborates is expected


def test_exact_single_exon_models_corroborate_symmetrically() -> None:
    assert_symmetric(model("a", "A", ((100, 400),)), model("b", "B", ((100, 400),)), True)


def test_positive_strand_three_prime_wobble_is_inclusive_at_three_bp() -> None:
    candidate = model("a", "A", ((100, 400),))
    assert_symmetric(candidate, model("b", "B", ((100, 403),)), True)
    result = compare_structural_corroboration(candidate, model("c", "C", ((100, 404),)))
    assert not result.corroborates
    assert result.reason == "single_exon_three_prime_boundary"


def test_positive_strand_five_prime_change_is_rejected() -> None:
    result = compare_structural_corroboration(
        model("a", "A", ((100, 400),)),
        model("b", "B", ((101, 400),)),
    )
    assert not result.corroborates
    assert result.reason == "single_exon_five_prime_boundary"


def test_negative_strand_reverses_five_and_three_prime_boundaries() -> None:
    candidate = model("a", "A", ((100, 400),), strand="-")
    # Genomic start is the biological 3-prime boundary on the negative strand.
    assert_symmetric(candidate, model("b", "B", ((97, 400),), strand="-"), True)
    result = compare_structural_corroboration(
        candidate,
        model("c", "C", ((100, 399),), strand="-"),
    )
    assert not result.corroborates
    assert result.reason == "single_exon_five_prime_boundary"


@pytest.mark.parametrize("parts", [((100, 350),), ((100, 451),), ((99, 400),)])
def test_single_exon_substantial_truncations_extensions_and_five_prime_extension_fail(parts: tuple[tuple[int, int], ...]) -> None:
    assert not compare_structural_corroboration(
        model("a", "A", ((100, 400),)), model("b", "B", parts)
    ).corroborates


def test_exact_multi_exon_intron_chain_and_coverage_pass() -> None:
    parts = ((0, 100), (200, 300), (400, 500))
    assert_symmetric(model("a", "A", parts), model("b", "B", parts), True)


@pytest.mark.parametrize(
    ("candidate_parts", "other_parts"),
    [
        (
            ((0, 100), (200, 300), (400, 500), (600, 700)),
            ((0, 100), (200, 300), (400, 501), (601, 700)),
        ),
        (
            ((0, 100), (200, 300), (400, 500), (600, 700), (800, 900)),
            ((0, 100), (200, 300), (400, 500), (600, 701), (801, 900)),
        ),
    ],
)
def test_one_differing_intron_is_allowed_for_three_and_four_intron_models(
    candidate_parts: tuple[tuple[int, int], ...],
    other_parts: tuple[tuple[int, int], ...],
) -> None:
    assert_symmetric(model("a", "A", candidate_parts), model("b", "B", other_parts), True)


def test_intron_threshold_is_applied_reciprocally() -> None:
    four_introns = model(
        "short",
        "A",
        ((0, 100), (200, 300), (400, 500), (600, 700), (800, 900)),
    )
    six_introns = model(
        "long",
        "B",
        ((0, 100), (200, 300), (400, 500), (600, 650), (675, 700), (725, 750), (800, 900)),
    )
    result = compare_structural_corroboration(four_introns, six_introns)
    assert len(four_introns.cds_intron_set & six_introns.cds_intron_set) == 3
    assert not result.corroborates
    assert result.reason == "intron_agreement"
    assert_symmetric(four_introns, six_introns, False)


def test_intron_threshold_boundaries() -> None:
    assert intron_agreement_requirement_met(2, 2)
    assert not intron_agreement_requirement_met(1, 2)
    assert intron_agreement_requirement_met(2, 3)
    assert intron_agreement_requirement_met(3, 4)
    assert intron_agreement_requirement_met(4, 5)
    assert not intron_agreement_requirement_met(3, 5)
    assert intron_agreement_requirement_met(8, 10)
    assert not intron_agreement_requirement_met(7, 10)


def test_reciprocal_spliced_cds_coverage_passes_at_exactly_85_percent() -> None:
    candidate = model("a", "A", ((0, 100), (200, 300), (400, 500)))
    other = model("b", "B", ((45, 100), (200, 300), (400, 500)))
    result = compare_structural_corroboration(candidate, other)
    assert result.corroborates
    assert result.candidate_cds_coverage == pytest.approx(0.85)
    assert result.other_cds_coverage == pytest.approx(1.0)
    assert_symmetric(candidate, other, True)


def test_reciprocal_spliced_cds_coverage_rejects_below_85_percent() -> None:
    candidate = model("a", "A", ((0, 100), (200, 300), (400, 500)))
    other = model("b", "B", ((46, 100), (200, 300), (400, 500)))
    result = compare_structural_corroboration(candidate, other)
    assert not result.corroborates
    assert result.reason == "reciprocal_cds_coverage"
    assert_symmetric(candidate, other, False)


def test_partial_model_sharing_allowed_intron_count_still_fails_coverage() -> None:
    complete = model(
        "complete",
        "A",
        ((0, 100), (200, 300), (400, 500), (600, 700), (800, 900)),
    )
    partial = model(
        "partial",
        "B",
        ((0, 100), (200, 300), (400, 500), (600, 700), (800, 805)),
    )
    result = compare_structural_corroboration(complete, partial)
    assert result.shared_intron_count == 4
    assert not result.corroborates
    assert result.reason == "reciprocal_cds_coverage"


def test_fused_model_supported_only_in_one_portion_is_rejected() -> None:
    fusion = model(
        "fusion",
        "A",
        ((0, 100), (200, 300), (400, 500), (1000, 1100), (1200, 1300), (1400, 1500)),
    )
    left_gene = model("left", "B", ((0, 100), (200, 300), (400, 500)))
    result = compare_structural_corroboration(fusion, left_gene)
    assert not result.corroborates
    assert result.reason == "intron_agreement"
    assert_symmetric(fusion, left_gene, False)


def test_single_and_multi_exon_models_never_corroborate() -> None:
    result = compare_structural_corroboration(
        model("single", "A", ((0, 100),)),
        model("multi", "B", ((0, 40), (60, 100))),
    )
    assert not result.corroborates
    assert result.reason == "exon_class_mismatch"


def test_multiple_isoforms_from_one_source_contribute_one_vote() -> None:
    candidate = model("candidate", "A", ((100, 400),))
    b1 = model("b1", "B", ((100, 400),))
    b2 = model("b2", "B", ((100, 402),))
    c1 = model("c1", "C", ((101, 400),))
    models = [candidate, b1, b2, c1]
    row = make_structural_corroboration_row(
        candidate,
        sublocus(models),
        [
            (b1, compare_structural_corroboration(candidate, b1)),
            (b2, compare_structural_corroboration(candidate, b2)),
            (c1, compare_structural_corroboration(candidate, c1)),
        ],
        run_source_count=4,
        exact_cds_exon_chain_support_source_count=1,
    )
    assert row["corroborating_other_source_count"] == 1
    assert row["corroborating_other_sources"] == "B"
    assert row["corroborating_model_count"] == 2
    assert row["sublocus_other_source_count"] == 2
    assert row["run_possible_other_source_count"] == 3


def test_no_other_run_source_is_deterministically_zero() -> None:
    candidate = model("candidate", "A", ((100, 400),))
    row = make_structural_corroboration_row(
        candidate,
        sublocus([candidate]),
        [],
        run_source_count=1,
        exact_cds_exon_chain_support_source_count=0,
    )
    assert row["corroborating_other_source_count"] == 0
    assert row["structural_corroboration_linear_score"] == 0.0
    assert row["structural_corroboration_log_score"] == 0.0
    assert row["structural_corroboration_saturating_alpha_8_score"] == 0.0


@pytest.mark.parametrize("n", [2, 3, 6, 20, 40, 80])
def test_all_normalizations_are_bounded_and_reach_one_at_full_support(n: int) -> None:
    for k in range(n + 1):
        values = corroboration_normalizations(k, n)
        assert all(0.0 <= value <= 1.0 for value in values.values())
    full = corroboration_normalizations(n, n)
    assert all(value == pytest.approx(1.0) for value in full.values())


def test_normalizations_separate_two_from_twenty_and_saturate() -> None:
    two = corroboration_normalizations(2, 82)
    twenty = corroboration_normalizations(20, 82)
    forty = corroboration_normalizations(40, 82)
    eighty = corroboration_normalizations(80, 82)
    for key in two:
        assert two[key] < twenty[key] < forty[key] < eighty[key]
    # Alpha 8 gives diminishing absolute gains as support becomes extensive.
    key = "structural_corroboration_saturating_alpha_8_score"
    assert forty[key] - twenty[key] > eighty[key] - forty[key]


def test_diagnostic_source_votes_are_stable_under_input_ordering() -> None:
    candidate = model("candidate", "A", ((100, 400),))
    b = model("b", "B", ((100, 400),))
    c = model("c", "C", ((100, 403),))
    comparisons = [
        (b, CorroborationComparison(True, "corroborates")),
        (c, CorroborationComparison(True, "corroborates")),
    ]
    first = make_structural_corroboration_row(
        candidate,
        sublocus([candidate, b, c]),
        comparisons,
        run_source_count=4,
        exact_cds_exon_chain_support_source_count=1,
    )
    second = make_structural_corroboration_row(
        candidate,
        sublocus([c, candidate, b]),
        list(reversed(comparisons)),
        run_source_count=4,
        exact_cds_exon_chain_support_source_count=1,
    )
    assert first == second


def test_invalid_normalization_counts_and_alphas_are_rejected() -> None:
    with pytest.raises(ValueError):
        corroboration_normalizations(2, 1)
    with pytest.raises(ValueError):
        corroboration_normalizations(0, 1, alphas=(0.0,))


def scoring_row(model_id: str, *, k: int, n: int) -> dict[str, object]:
    return {
        "output_transcript_id": model_id,
        "source_name": "source_a",
        "sublocus_id": "sub1",
        "cds_length": 300,
        "corroborating_other_source_count": k,
        "run_possible_other_source_count": n,
    }


def scoring_weights(weight: float) -> dict[str, float]:
    weights = {key: 0.0 for key in MODEL_SELECTION_DEFAULT_WEIGHTS}
    weights["structural_corroboration"] = weight
    return weights


@pytest.mark.parametrize(("source_count", "k"), [(3, 1), (7, 2), (83, 20)])
def test_active_scoring_uses_run_depth_and_alpha_16(source_count: int, k: int) -> None:
    row = scoring_row("candidate", k=k, n=source_count - 1)
    apply_model_selection_scores(
        [row], {}, scoring_weights(0.5), 1e-12,
        structural_corroboration_alpha=16.0,
    )
    expected = saturating_corroboration_score(k, source_count - 1, 16.0)
    assert row["structural_corroboration_score"] == pytest.approx(expected)
    assert row["structural_corroboration_contribution"] == pytest.approx(0.5 * expected)
    assert row["model_selection_score"] == pytest.approx(0.5 * expected)
    assert row["structural_corroboration_in_selection"] == "true"


@pytest.mark.parametrize(("k", "n", "expected"), [(0, 0, 0.0), (0, 6, 0.0), (6, 6, 1.0)])
def test_zero_and_maximum_corroboration_scoring(k: int, n: int, expected: float) -> None:
    row = scoring_row("candidate", k=k, n=n)
    apply_model_selection_scores(
        [row], {}, scoring_weights(0.5), 1e-12,
        structural_corroboration_alpha=16.0,
    )
    assert row["structural_corroboration_score"] == pytest.approx(expected)
    assert row["structural_corroboration_contribution"] == pytest.approx(0.5 * expected)


def test_weight_zero_preserves_score_and_omits_component_text() -> None:
    without_counts = scoring_row("without", k=0, n=0)
    with_counts = scoring_row("with", k=5, n=6)
    weights = scoring_weights(0.0)
    apply_model_selection_scores([without_counts], {}, weights, 1e-12)
    apply_model_selection_scores([with_counts], {}, weights, 1e-12)
    assert with_counts["model_selection_score"] == without_counts["model_selection_score"]
    assert with_counts["model_selection_score_scaled"] == without_counts["model_selection_score_scaled"]
    assert "structural_corroboration_score" not in str(with_counts["model_selection_score_components"])
    assert with_counts["structural_corroboration_in_selection"] == "false"


def test_alpha_8_and_alternative_weight_are_configurable() -> None:
    row = scoring_row("candidate", k=2, n=6)
    apply_model_selection_scores(
        [row], {}, scoring_weights(0.75), 1e-12,
        structural_corroboration_alpha=8.0,
    )
    expected = saturating_corroboration_score(2, 6, 8.0)
    assert row["structural_corroboration_alpha"] == 8.0
    assert row["structural_corroboration_weight"] == 0.75
    assert row["model_selection_score"] == pytest.approx(0.75 * expected)


def test_negative_original_raw_score_remains_zero_on_scaled_axis() -> None:
    row = scoring_row("candidate", k=1, n=6)
    row["intrinsic_penalty_score"] = 1.0
    weights = scoring_weights(0.5)
    weights["intrinsic_penalty"] = 1.0
    apply_model_selection_scores([row], {}, weights, 1e-12)
    assert float(row["model_selection_score"]) < 0
    assert row["model_selection_score_scaled"] == 0.0


def test_alpha_16_weight_half_increases_no_evidence_max_to_4_85() -> None:
    weights = dict(MODEL_SELECTION_DEFAULT_WEIGHTS)
    weights["structural_corroboration"] = 0.5
    assert model_selection_score_positive_scale_max(
        weights,
        include_protein_evidence=False,
        include_transcript_evidence=False,
    ) == pytest.approx(4.85)


def test_external_evidence_positive_max_includes_corroboration_once() -> None:
    weights = dict(MODEL_SELECTION_DEFAULT_WEIGHTS)
    weights["structural_corroboration"] = 0.5
    assert model_selection_score_positive_scale_max(
        weights,
        include_protein_evidence=True,
        include_transcript_evidence=True,
    ) == pytest.approx(9.85)


def test_equal_corroboration_scores_keep_deterministic_identifier_order() -> None:
    rows = [scoring_row("model_b", k=2, n=6), scoring_row("model_a", k=2, n=6)]
    apply_model_selection_scores(
        rows, {}, scoring_weights(0.5), 1e-12,
        structural_corroboration_alpha=16.0,
    )
    ranks = {str(row["output_transcript_id"]): row["model_selection_order_rank"] for row in rows}
    assert ranks == {"model_b": 2, "model_a": 1}
