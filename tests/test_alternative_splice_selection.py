from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.alternative_splice_selection import (
    ADJACENT_CANONICAL_PRESET_COMPARISONS,
    CANONICAL_ALTERNATIVE_SPLICE_PRESETS,
    AlternativeSpliceEvidenceRequirement,
    AlternativeSplicePreset,
    AlternativeSpliceSelectionConfig,
    AlternativeSpliceSupportPolicy,
    cds_base_containment_relation,
    cds_interval_union_is_subset,
    cds_union_overlap_metrics,
    effective_min_model_support_sources,
    extend_alternative_splice_report,
    relaxed_contiguous_cds_chain_containment_relation,
    select_alternative_splice_candidates,
    union_cds_intervals,
)
from moselecta.final_annotation import augment_annotation, build_primary_annotation
from moselecta.models import CodingModel, FeatureRecord


def _record(start0: int, end0: int, *, strand: str = "+") -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1",
        source="test",
        feature_type="CDS",
        start=start0 + 1,
        end=end0,
        score=".",
        strand=strand,
        phase="0",
        attrs={"Parent": "tx"},
        raw_attrs="Parent=tx",
        input_file=Path("test.gff3"),
        source_name="test",
        line_number=1,
        record_index=1,
        raw_line="",
    )


def _model(
    name: str,
    parts: list[tuple[int, int]],
    *,
    source: str,
    strand: str = "+",
) -> CodingModel:
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"input_gene_{name}",
        input_transcript_id=f"input_tx_{name}",
        output_gene_id=f"internal_gene_{name}",
        output_transcript_id=name,
        chrom="chr1",
        strand=strand,
        cds_parts=parts,
        gene_record=None,
        transcript_record=None,
        child_records=[_record(start0, end0, strand=strand) for start0, end0 in parts],
        superlocus_id="super1",
        sublocus_id="sub1",
    )


def _score_row(raw: float, scaled: float | None = None) -> dict[str, object]:
    return {
        "model_selection_score": raw,
        "model_selection_score_scaled": max(0.0, raw) if scaled is None else scaled,
        "model_structural_support_score_scaled": 0.8,
        "protein_evidence_status": "protein_evidence_not_configured",
        "transcript_evidence_status": "transcript_evidence_not_configured",
        "sequence_integrity_status": "genome_not_configured",
    }


def _candidate_row(
    candidate: CodingModel,
    primary: CodingModel,
    *,
    candidate_raw: float = 0.8,
    candidate_scaled: float = 0.8,
    primary_raw: float = 1.0,
    rank: int = 1,
    exact_support: bool = True,
    protein_support: bool = False,
    transcript_support: bool = False,
    full_mosaic: bool = False,
    exact_model_source_count: int | None = None,
    mosaic_min_component_source_count: int | None = None,
    protein_strong_partial: bool = False,
    transcript_strong_partial: bool = False,
) -> dict[str, object]:
    historical_ratio: object = (
        candidate_raw / primary_raw if primary_raw != 0 else "."
    )
    row = {
        "alternative_candidate_assessment_id": f"assessment:{candidate.output_transcript_id}",
        "candidate_output_transcript_id": candidate.output_transcript_id,
        "primary_output_transcript_id": primary.output_transcript_id,
        "ownership_status": "unique_owner",
        "structural_candidate": True,
        "alternative_candidate_order_rank": rank,
        "exact_cds_equivalent_to_primary": False,
        "same_intron_chain_terminal_variant": False,
        "frame_concordance_status": "calculated",
        "sequence_integrity_default_exclusion": False,
        "intrinsic_default_exclusion": False,
        "exact_multisource_model_support": exact_support,
        "protein_exact_cds_chain_support": protein_support,
        "exact_coherent_transcript_support": transcript_support,
        "full_structural_mosaic": full_mosaic,
        "protein_strong_partial": protein_strong_partial,
        "transcript_strong_partial": transcript_strong_partial,
        "mosaic_cds_introns_status": "complete",
        "mosaic_cds_splice_sites_status": "complete",
        "mosaic_internal_cds_exons_status": "complete",
        "mosaic_terminal_cds_exons_status": "complete",
        "mosaic_cds_base_consensus_status": "complete",
        "reciprocal_cds_overlap": 0.90,
        "same_frame_overlap_fraction": 0.95,
        "candidate_model_selection_score": candidate_raw,
        "candidate_model_selection_score_scaled": candidate_scaled,
        "primary_model_selection_score": primary_raw,
        "candidate_primary_model_selection_score_ratio": historical_ratio,
    }
    if exact_model_source_count is not None:
        row["alternative_splice_exact_model_source_count"] = (
            exact_model_source_count
        )
    if mosaic_min_component_source_count is not None:
        row["alternative_splice_mosaic_min_component_source_count"] = (
            mosaic_min_component_source_count
        )
    return row


def _selection(
    primary: CodingModel,
    candidates: list[CodingModel],
    rows: list[dict[str, object]],
    *,
    config: AlternativeSpliceSelectionConfig | None = None,
    extra_models: list[CodingModel] | None = None,
):
    score_rows = {
        primary.output_transcript_id: _score_row(
            float(rows[0]["primary_model_selection_score"])
        )
    }
    score_rows.update({
        candidate.output_transcript_id: _score_row(
            float(row["candidate_model_selection_score"]),
            float(row["candidate_model_selection_score_scaled"]),
        )
        for candidate, row in zip(candidates, rows)
    })
    annotation = build_primary_annotation([primary], score_rows)
    result = select_alternative_splice_candidates(
        candidate_rows=rows,
        all_scored_models=[primary, *candidates, *(extra_models or [])],
        primary_annotation=annotation,
        config=config or AlternativeSpliceSelectionConfig.resolve(),
        terminal_cds_3prime_wobble_bp=3,
    )
    return result, annotation, score_rows


def _distinct_pair() -> tuple[CodingModel, CodingModel]:
    primary = _model(
        "primary", [(0, 90), (180, 270), (360, 450)], source="P"
    )
    candidate = _model(
        "candidate", [(0, 90), (180, 240), (330, 450)], source="A"
    )
    return primary, candidate


def test_canonical_preset_values_and_defaults() -> None:
    permissive = CANONICAL_ALTERNATIVE_SPLICE_PRESETS[AlternativeSplicePreset.PERMISSIVE]
    intermediate = CANONICAL_ALTERNATIVE_SPLICE_PRESETS[AlternativeSplicePreset.INTERMEDIATE]
    stringent = CANONICAL_ALTERNATIVE_SPLICE_PRESETS[AlternativeSplicePreset.STRINGENT]
    assert (
        permissive.min_reciprocal_cds_overlap,
        permissive.min_frame_concordance,
        permissive.min_model_selection_score_scaled,
        permissive.min_primary_score_ratio,
    ) == (0.80, 0.95, 0.50, 0.70)
    assert (
        intermediate.min_reciprocal_cds_overlap,
        intermediate.min_frame_concordance,
        intermediate.min_model_selection_score_scaled,
        intermediate.min_primary_score_ratio,
    ) == (0.85, 0.95, 0.60, 0.70)
    assert (
        stringent.min_reciprocal_cds_overlap,
        stringent.min_frame_concordance,
        stringent.min_model_selection_score_scaled,
        stringent.min_primary_score_ratio,
    ) == (0.85, 0.95, 0.60, 0.80)
    config = AlternativeSpliceSelectionConfig.resolve()
    assert config.requested_preset is AlternativeSplicePreset.INTERMEDIATE
    assert config.support_policy is AlternativeSpliceSupportPolicy.EXACT
    assert config.evidence_requirement is AlternativeSpliceEvidenceRequirement.NONE
    assert config.max_alternatives_per_gene == 50
    assert config.effective_thresholds == intermediate


def test_adjacent_canonical_preset_comparisons_are_unique_and_complete() -> None:
    assert ADJACENT_CANONICAL_PRESET_COMPARISONS == (
        (AlternativeSplicePreset.PERMISSIVE, AlternativeSplicePreset.INTERMEDIATE),
        (AlternativeSplicePreset.INTERMEDIATE, AlternativeSplicePreset.STRINGENT),
    )
    assert len(ADJACENT_CANONICAL_PRESET_COMPARISONS) == len(
        set(ADJACENT_CANONICAL_PRESET_COMPARISONS)
    )


def test_explicit_overrides_change_only_active_effective_values() -> None:
    config = AlternativeSpliceSelectionConfig.resolve(
        preset="stringent",
        min_reciprocal_cds_overlap=0.81,
        min_frame_concordance=0.91,
        min_model_selection_score_scaled=0.55,
        min_primary_score_ratio=0.61,
    )
    assert config.canonical_thresholds == CANONICAL_ALTERNATIVE_SPLICE_PRESETS[
        AlternativeSplicePreset.STRINGENT
    ]
    assert config.effective_thresholds.min_primary_score_ratio == 0.61
    assert set(config.overridden_thresholds) == {
        "min_reciprocal_cds_overlap",
        "min_frame_concordance",
        "min_model_selection_score_scaled",
        "min_primary_score_ratio",
    }


@pytest.mark.parametrize(
    ("source_count", "fraction", "expected"),
    [
        (5, 0.10, 2),
        (5, 0.50, 3),
        (10, 0.30, 3),
        (83, 0.10, 9),
        (83, 0.25, 21),
        (30, 0.10, 3),
    ],
)
def test_fractional_model_source_requirement_uses_robust_global_ceiling(
    source_count: int, fraction: float, expected: int
) -> None:
    assert effective_min_model_support_sources(
        prediction_source_count=source_count,
        configured_minimum=2,
        configured_fraction=fraction,
    ) == expected


def test_stricter_absolute_or_fractional_model_source_requirement_wins() -> None:
    assert effective_min_model_support_sources(
        prediction_source_count=5,
        configured_minimum=4,
        configured_fraction=0.50,
    ) == 4
    assert effective_min_model_support_sources(
        prediction_source_count=10,
        configured_minimum=2,
        configured_fraction=0.40,
    ) == 4
    assert effective_min_model_support_sources(
        prediction_source_count=2,
        configured_minimum=4,
    ) == 4


@pytest.mark.parametrize(
    ("minimum", "source_count", "accepted"),
    [
        (2, 2, True),
        (3, 2, False),
        (3, 3, True),
        (4, 3, False),
        (4, 4, True),
    ],
)
def test_absolute_model_source_requirement_controls_exact_route(
    minimum: int, source_count: int, accepted: bool
) -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(
        candidate,
        primary,
        exact_support=source_count >= 2,
        exact_model_source_count=source_count,
        mosaic_min_component_source_count=1,
    )
    result, _, _ = _selection(
        primary,
        [candidate],
        [row],
        config=AlternativeSpliceSelectionConfig.resolve(
            min_model_support_sources=minimum
        ),
    )
    assert bool(result.accepted_alternatives) is accepted
    decision = result.candidate_decisions[0]
    assert decision.exact_model_source_requirement_met is accepted


def test_mosaic_source_requirement_is_a_raw_k_qualified_fact() -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(
        candidate,
        primary,
        exact_support=True,
        full_mosaic=True,
        exact_model_source_count=3,
        mosaic_min_component_source_count=3,
    )
    result, _, _ = _selection(
        primary,
        [candidate],
        [row],
        config=AlternativeSpliceSelectionConfig.resolve(
            support_policy="exact-plus-mosaic",
            min_model_support_sources=3,
        ),
    )
    decision = result.candidate_decisions[0]
    assert decision.exact_model_source_requirement_met
    assert decision.mosaic_source_requirement_met
    assert decision.qualifying_support_routes == (
        "exact_multisource_model_support",
    )


def test_raised_k_exact_below_threshold_can_qualify_through_external_routes() -> None:
    primary, candidate = _distinct_pair()
    for protein_support, transcript_support in ((True, False), (False, True)):
        row = _candidate_row(
            candidate,
            primary,
            exact_support=True,
            exact_model_source_count=2,
            mosaic_min_component_source_count=1,
            protein_support=protein_support,
            transcript_support=transcript_support,
        )
        result, _, _ = _selection(
            primary,
            [candidate],
            [row],
            config=AlternativeSpliceSelectionConfig.resolve(
                min_model_support_sources=3
            ),
        )
        decision = result.candidate_decisions[0]
        assert decision.stage2b_decision == "accepted"
        assert not decision.exact_model_source_requirement_met


def test_raised_k_exact_below_threshold_does_not_block_mosaic_route() -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(
        candidate,
        primary,
        exact_support=True,
        full_mosaic=False,
        exact_model_source_count=2,
        mosaic_min_component_source_count=3,
    )
    exact, _, _ = _selection(
        primary,
        [candidate],
        [row],
        config=AlternativeSpliceSelectionConfig.resolve(
            support_policy="exact", min_model_support_sources=3
        ),
    )
    mosaic, _, _ = _selection(
        primary,
        [candidate],
        [row],
        config=AlternativeSpliceSelectionConfig.resolve(
            support_policy="exact-plus-mosaic", min_model_support_sources=3
        ),
    )
    assert exact.candidate_decisions[0].stage2b_decision == (
        "rejected_no_qualifying_support"
    )
    mosaic_decision = mosaic.candidate_decisions[0]
    assert mosaic_decision.stage2b_decision == "accepted"
    assert not mosaic_decision.exact_model_source_requirement_met
    assert mosaic_decision.mosaic_source_requirement_met
    assert mosaic_decision.mosaic_policy_required


def test_absolute_three_and_fraction_half_are_equivalent_with_five_sources() -> None:
    primary, candidate = _distinct_pair()
    extra_models = [
        _model(f"support_{source}", list(candidate.cds_parts), source=source)
        for source in ("B", "C", "D")
    ]
    row = _candidate_row(
        candidate,
        primary,
        exact_support=True,
        exact_model_source_count=3,
        mosaic_min_component_source_count=1,
    )
    results = []
    for config in (
        AlternativeSpliceSelectionConfig.resolve(
            min_model_support_sources=3
        ),
        AlternativeSpliceSelectionConfig.resolve(
            min_model_support_source_fraction=0.50
        ),
    ):
        result, _, _ = _selection(
            primary,
            [candidate],
            [row],
            config=config,
            extra_models=extra_models,
        )
        results.append(result)
    assert all(
        result.config.prediction_source_count == 5
        and result.config.effective_min_model_support_sources == 3
        for result in results
    )
    assert [
        item.candidate_model.output_transcript_id
        for item in results[0].accepted_alternatives
    ] == [
        item.candidate_model.output_transcript_id
        for item in results[1].accepted_alternatives
    ]
    assert [
        item.accepted_public_transcript_id
        for item in results[0].accepted_alternatives
    ] == [
        item.accepted_public_transcript_id
        for item in results[1].accepted_alternatives
    ]


def test_one_weak_mosaic_component_fails_raised_k() -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(
        candidate,
        primary,
        exact_support=False,
        full_mosaic=True,
        exact_model_source_count=1,
        mosaic_min_component_source_count=2,
    )
    result, _, _ = _selection(
        primary,
        [candidate],
        [row],
        config=AlternativeSpliceSelectionConfig.resolve(
            support_policy="exact-plus-mosaic", min_model_support_sources=3
        ),
    )
    assert not result.accepted_alternatives
    assert not result.candidate_decisions[0].mosaic_source_requirement_met


@pytest.mark.parametrize(
    "requirement",
    ("protein-or-transcript", "protein", "transcript", "protein-and-transcript"),
)
def test_raised_k_preserves_non_none_exact_vs_mosaic_accepted_invariant(
    requirement: str,
) -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(
        candidate,
        primary,
        exact_support=True,
        exact_model_source_count=2,
        mosaic_min_component_source_count=3,
        protein_support=True,
        transcript_support=True,
    )
    results = []
    for policy in ("exact", "exact-plus-mosaic"):
        result, _, _ = _selection(
            primary,
            [candidate],
            [row],
            config=AlternativeSpliceSelectionConfig.resolve(
                support_policy=policy,
                evidence_requirement=requirement,
                min_model_support_sources=3,
            ),
        )
        results.append(result)
    assert [
        item.candidate_model.output_transcript_id
        for item in results[0].accepted_alternatives
    ] == [
        item.candidate_model.output_transcript_id
        for item in results[1].accepted_alternatives
    ]


SUPPORT_COMBINATIONS = (
    ("model_only", True, False, False, False, False, False),
    ("protein_only", False, True, False, False, False, False),
    ("transcript_only", False, False, True, False, False, False),
    ("model_protein", True, True, False, False, False, False),
    ("model_transcript", True, False, True, False, False, False),
    ("protein_transcript", False, True, True, False, False, False),
    ("all_exact", True, True, True, False, False, False),
    ("mosaic_only", False, False, False, True, False, False),
    ("mosaic_partial_protein", False, False, False, True, True, False),
    ("mosaic_partial_transcript", False, False, False, True, False, True),
)


@pytest.mark.parametrize("support_policy", ("exact", "exact-plus-mosaic"))
@pytest.mark.parametrize(
    ("requirement", "accepted_names"),
    (
        ("none", None),
        (
            "protein-or-transcript",
            {
                "protein_only", "transcript_only", "model_protein",
                "model_transcript", "protein_transcript", "all_exact",
            },
        ),
        (
            "protein",
            {"protein_only", "model_protein", "protein_transcript", "all_exact"},
        ),
        (
            "transcript",
            {"transcript_only", "model_transcript", "protein_transcript", "all_exact"},
        ),
        ("protein-and-transcript", {"protein_transcript", "all_exact"}),
    ),
)
def test_support_policy_and_external_evidence_requirement_matrix(
    support_policy: str,
    requirement: str,
    accepted_names: set[str] | None,
) -> None:
    primary, candidate = _distinct_pair()
    for (
        name, model_support, protein_support, transcript_support, mosaic,
        partial_protein, partial_transcript,
    ) in SUPPORT_COMBINATIONS:
        row = _candidate_row(
            candidate,
            primary,
            exact_support=model_support,
            protein_support=protein_support,
            transcript_support=transcript_support,
            full_mosaic=mosaic,
            protein_strong_partial=partial_protein,
            transcript_strong_partial=partial_transcript,
        )
        result, _, _ = _selection(
            primary,
            [candidate],
            [row],
            config=AlternativeSpliceSelectionConfig.resolve(
                support_policy=support_policy,
                evidence_requirement=requirement,
            ),
        )
        accepted = bool(result.accepted_alternatives)
        if requirement == "none":
            expected = model_support or protein_support or transcript_support
            if support_policy == "exact-plus-mosaic":
                expected = expected or mosaic
        else:
            assert accepted_names is not None
            expected = name in accepted_names
        assert accepted is expected, (support_policy, requirement, name)
        decision = result.candidate_decisions[0]
        expected_external_routes = tuple(
            route
            for present, route in (
                (protein_support, "exact_protein_support"),
                (transcript_support, "exact_coherent_transcript_support"),
            )
            if present
        )
        assert decision.external_evidence_requirement_routes == expected_external_routes
        if decision.support_policy_eligible and not expected:
            assert decision.stage2b_decision == "rejected_external_evidence_requirement"


@pytest.mark.parametrize(
    "requirement",
    ("protein-or-transcript", "protein", "transcript", "protein-and-transcript"),
)
def test_non_none_external_requirement_makes_exact_and_mosaic_sets_identical(
    requirement: str,
) -> None:
    primary = _model(
        "primary", [(0, 90), (180, 270), (360, 450)], source="P"
    )
    candidates: list[CodingModel] = []
    rows: list[dict[str, object]] = []
    for index, combination in enumerate(SUPPORT_COMBINATIONS, start=1):
        (
            name, model_support, protein_support, transcript_support, mosaic,
            partial_protein, partial_transcript,
        ) = combination
        candidate = _model(
            name,
            [(0, 90), (180, 225 + index), (315 + index, 450)],
            source=f"S{index}",
        )
        candidates.append(candidate)
        rows.append(_candidate_row(
            candidate,
            primary,
            rank=index,
            exact_support=model_support,
            protein_support=protein_support,
            transcript_support=transcript_support,
            full_mosaic=mosaic,
            protein_strong_partial=partial_protein,
            transcript_strong_partial=partial_transcript,
        ))
    accepted_by_policy = {}
    for policy in ("exact", "exact-plus-mosaic"):
        result, _, _ = _selection(
            primary,
            candidates,
            rows,
            config=AlternativeSpliceSelectionConfig.resolve(
                support_policy=policy,
                evidence_requirement=requirement,
            ),
        )
        accepted_by_policy[policy] = tuple(
            item.candidate_model.output_transcript_id
            for item in result.accepted_alternatives
        )
    assert accepted_by_policy["exact"] == accepted_by_policy["exact-plus-mosaic"]


@pytest.mark.parametrize("value", [-0.01, 1.01, float("inf"), float("nan")])
def test_invalid_threshold_ranges(value: float) -> None:
    with pytest.raises(ValueError):
        AlternativeSpliceSelectionConfig.resolve(
            min_primary_score_ratio=value
        )


def test_cds_interval_union_subset_does_not_double_count() -> None:
    assert union_cds_intervals([(0, 10), (0, 10), (5, 20), (20, 25)]) == ((0, 25),)
    assert cds_interval_union_is_subset(
        [(2, 8), (6, 12), (2, 8)],
        [(0, 5), (5, 20), (0, 20)],
    )
    assert not cds_interval_union_is_subset([(2, 21)], [(0, 20)])


def test_cds_union_overlap_metrics_do_not_double_count_duplicate_parts() -> None:
    retained = _model(
        "retained", [(0, 100), (0, 100), (50, 150)], source="P"
    )
    candidate = _model("candidate", [(25, 125), (25, 125)], source="A")
    metrics = cds_union_overlap_metrics(retained, candidate)
    assert metrics.shared_bp == 100
    assert metrics.retained_fraction == pytest.approx(100 / 150)
    assert metrics.candidate_fraction == 1.0
    assert metrics.reciprocal_fraction == pytest.approx(100 / 150)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_broad_cds_base_containment_both_directions(strand: str) -> None:
    primary = _model("primary", [(0, 100), (200, 300)], source="P", strand=strand)
    inner = _model("inner", [(20, 80), (220, 280)], source="A", strand=strand)
    outer = _model("outer", [(-20, 120), (180, 320)], source="B", strand=strand)
    assert cds_base_containment_relation(inner, primary).candidate_in_retained
    assert cds_base_containment_relation(outer, primary).retained_in_candidate


def test_genomic_span_containment_is_not_cds_base_containment() -> None:
    primary = _model("primary", [(0, 100), (200, 300)], source="P")
    intronic = _model("intronic", [(120, 180)], source="A")
    assert primary.cds_start0 < intronic.cds_start0 < intronic.cds_end0 < primary.cds_end0
    assert not cds_base_containment_relation(intronic, primary).any_direction


@pytest.mark.parametrize("strand", ["+", "-"])
def test_relaxed_contiguous_chain_terminal_trim_is_containment(strand: str) -> None:
    full = _model(
        "full",
        [(0, 100), (200, 300), (400, 500), (600, 700)],
        source="P",
        strand=strand,
    )
    trimmed_e2_e3 = _model(
        "trimmed",
        [(220, 300), (400, 480)],
        source="A",
        strand=strand,
    )
    relation = relaxed_contiguous_cds_chain_containment_relation(
        trimmed_e2_e3, full
    )
    assert relation.candidate_in_retained
    assert not relation.retained_in_candidate

    row = _candidate_row(trimmed_e2_e3, full)
    result, _, _ = _selection(full, [trimmed_e2_e3], [row])
    decision = result.candidate_decisions[0]
    assert decision.stage2b_decision == "rejected_candidate_contained_in_primary"
    assert decision.candidate_relaxed_cds_chain_contained_in_primary


def test_relaxed_contiguous_chain_reverse_e2_e3_in_full_chain() -> None:
    primary_e2_e3 = _model(
        "primary", [(200, 300), (400, 500)], source="P"
    )
    candidate_full = _model(
        "candidate",
        [(0, 100), (200, 300), (400, 500), (600, 700)],
        source="A",
    )
    relation = relaxed_contiguous_cds_chain_containment_relation(
        candidate_full, primary_e2_e3
    )
    assert relation.retained_in_candidate
    assert not relation.candidate_in_retained

    row = _candidate_row(candidate_full, primary_e2_e3)
    result, _, _ = _selection(primary_e2_e3, [candidate_full], [row])
    decision = result.candidate_decisions[0]
    assert decision.stage2b_decision == "rejected_primary_contained_in_candidate"
    assert decision.primary_relaxed_cds_chain_contained_in_candidate


@pytest.mark.parametrize(
    ("name", "primary_parts", "candidate_parts"),
    [
        (
            "exon_skipping",
            [(0, 100), (200, 300), (400, 500), (600, 700)],
            [(0, 100), (400, 500), (600, 700)],
        ),
        (
            "intron_retention",
            [(0, 100), (200, 300), (400, 500)],
            [(0, 100), (200, 500)],
        ),
        (
            "additional_internal_intron",
            [(0, 100), (200, 300), (400, 500)],
            [(0, 100), (200, 240), (260, 300), (400, 500)],
        ),
        (
            "alternative_internal_donor_acceptor",
            [(0, 100), (200, 300), (400, 500)],
            [(0, 90), (210, 300), (400, 500)],
        ),
        (
            "same_outer_span_different_internal_splicing",
            [(0, 100), (200, 300), (400, 500)],
            [(0, 120), (220, 280), (380, 500)],
        ),
    ],
)
@pytest.mark.parametrize("strand", ["+", "-"])
def test_splice_distinct_structures_are_not_relaxed_chain_containment(
    name: str,
    primary_parts: list[tuple[int, int]],
    candidate_parts: list[tuple[int, int]],
    strand: str,
) -> None:
    primary = _model("primary", primary_parts, source="P", strand=strand)
    candidate = _model(name, candidate_parts, source="A", strand=strand)
    relation = relaxed_contiguous_cds_chain_containment_relation(candidate, primary)
    assert not relation.any_direction


def test_single_exon_primary_two_exon_candidate_ignores_base_union_containment() -> None:
    primary = _model("primary", [(0, 300)], source="P")
    candidate = _model("candidate", [(0, 120), (180, 300)], source="A")
    assert cds_base_containment_relation(candidate, primary).candidate_in_retained
    assert not relaxed_contiguous_cds_chain_containment_relation(
        candidate, primary
    ).any_direction
    row = _candidate_row(candidate, primary)
    result, _, _ = _selection(
        primary,
        [candidate],
        [row],
        config=AlternativeSpliceSelectionConfig.resolve(
            min_reciprocal_cds_overlap=0.0
        ),
    )
    assert result.candidate_decisions[0].stage2b_decision == "accepted"


def test_structural_rejection_precedence_exact_then_intron_chain_then_containment() -> None:
    primary = _model("primary", [(0, 100), (200, 300)], source="P")
    exact = _model("exact", [(0, 100), (200, 300)], source="A")
    exact_row = _candidate_row(exact, primary)
    exact_row["exact_cds_equivalent_to_primary"] = True
    exact_result, _, _ = _selection(primary, [exact], [exact_row])
    exact_decision = exact_result.candidate_decisions[0]
    assert exact_decision.bidirectional_cds_base_containment is True
    assert exact_decision.bidirectional_relaxed_cds_chain_containment is True
    assert exact_decision.stage2b_decision == "rejected_exact_equivalent"

    terminal = _model("terminal", [(20, 100), (200, 280)], source="B")
    terminal_row = _candidate_row(terminal, primary)
    terminal_row["same_intron_chain_terminal_variant"] = True
    terminal_result, _, _ = _selection(primary, [terminal], [terminal_row])
    terminal_decision = terminal_result.candidate_decisions[0]
    assert terminal_decision.candidate_cds_bases_contained_in_primary is True
    assert terminal_decision.candidate_relaxed_cds_chain_contained_in_primary is True
    assert terminal_decision.stage2b_decision == "rejected_same_intron_chain_terminal_variant"


@pytest.mark.parametrize(
    ("primary_raw", "candidate_raw", "candidate_scaled", "status", "ratio", "decision"),
    [
        (1.0, 0.7, 0.7, "available", 0.7, "accepted"),
        (1.0, 0.0, 0.0, "available", 0.0, "rejected_primary_score_ratio"),
        (1.0, -0.2, 0.0, "available", -0.2, "rejected_primary_score_ratio"),
        (0.0, 0.7, 0.7, "unavailable_zero_primary", None, "rejected_primary_score_ratio_unavailable"),
        (-1.0, -0.5, 0.0, "unavailable_negative_primary", None, "rejected_primary_score_ratio_unavailable"),
        (-1.0, -2.0, 0.0, "unavailable_negative_primary", None, "rejected_primary_score_ratio_unavailable"),
        (-1.0, 0.7, 0.7, "unavailable_negative_primary", None, "rejected_primary_score_ratio_unavailable"),
    ],
)
def test_stage2b_primary_score_ratio_policy(
    primary_raw: float,
    candidate_raw: float,
    candidate_scaled: float,
    status: str,
    ratio: float | None,
    decision: str,
) -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(
        candidate,
        primary,
        primary_raw=primary_raw,
        candidate_raw=candidate_raw,
        candidate_scaled=candidate_scaled,
    )
    config = AlternativeSpliceSelectionConfig.resolve(
        min_model_selection_score_scaled=0.0,
        min_primary_score_ratio=0.70,
    )
    result, _, _ = _selection(primary, [candidate], [row], config=config)
    actual = result.candidate_decisions[0]
    assert actual.stage2b_primary_score_ratio_status == status
    assert actual.stage2b_candidate_primary_model_selection_score_ratio == ratio
    assert actual.stage2b_decision == decision


def test_historical_negative_denominator_ratio_is_preserved_but_never_eligible() -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(
        candidate, primary, primary_raw=-1.0, candidate_raw=-0.5, candidate_scaled=0.0
    )
    assert row["candidate_primary_model_selection_score_ratio"] == 0.5
    config = AlternativeSpliceSelectionConfig.resolve(
        min_model_selection_score_scaled=0.0,
        min_primary_score_ratio=0.0,
    )
    selection, _, _ = _selection(primary, [candidate], [row], config=config)
    historical_report = {
        "candidate_rows": [row],
        "gene_rows": [{"primary_output_transcript_id": primary.output_transcript_id}],
        "run_rows": [{}],
    }
    extended = extend_alternative_splice_report(
        historical_report, selection, alternatives_emitted=True
    )
    extended_row = extended["candidate_rows"][0]
    assert extended_row["candidate_primary_model_selection_score_ratio"] == 0.5
    assert extended_row["alternative_splice_candidate_primary_model_selection_score_ratio"] == "."
    assert extended_row["alternative_splice_primary_score_ratio_status"] == "unavailable_negative_primary"
    assert extended_row["alternative_splice_decision"] == "rejected_primary_score_ratio_unavailable"


def test_ratio_unavailability_summary_distinguishes_candidates_and_assessments() -> None:
    positive_primary = _model(
        "positive_primary", [(0, 90), (180, 270), (360, 450)], source="P"
    )
    negative_primary_a = _model(
        "negative_primary_a", [(0, 90), (180, 260), (350, 450)], source="N1"
    )
    negative_primary_b = _model(
        "negative_primary_b", [(0, 90), (190, 270), (360, 450)], source="N2"
    )
    candidate = _model(
        "candidate", [(0, 90), (180, 240), (330, 450)], source="A"
    )
    primaries = [positive_primary, negative_primary_a, negative_primary_b]
    primary_scores = [1.0, -1.0, -2.0]
    rows = []
    for index, (primary, primary_score) in enumerate(
        zip(primaries, primary_scores), start=1
    ):
        row = _candidate_row(candidate, primary, primary_raw=primary_score)
        row["alternative_candidate_assessment_id"] = f"assessment:{index}"
        row["ownership_status"] = "ambiguous_multiple_owners"
        rows.append(row)

    score_rows = {
        primary.output_transcript_id: _score_row(score)
        for primary, score in zip(primaries, primary_scores)
    }
    score_rows[candidate.output_transcript_id] = _score_row(0.8, 0.8)
    annotation = build_primary_annotation(primaries, score_rows)
    selection = select_alternative_splice_candidates(
        candidate_rows=rows,
        all_scored_models=[*primaries, candidate],
        primary_annotation=annotation,
        config=AlternativeSpliceSelectionConfig.resolve(),
        terminal_cds_3prime_wobble_bp=3,
    )
    extended = extend_alternative_splice_report(
        {"candidate_rows": rows, "gene_rows": [], "run_rows": [{}]},
        selection,
        alternatives_emitted=True,
    )
    summary = extended["run_rows"][0]
    assert summary["alternative_splice_primary_score_ratio_unavailable_negative_count"] == 1
    assert (
        summary[
            "alternative_splice_primary_score_ratio_unavailable_negative_assessment_count"
        ]
        == 2
    )
    assert sum(
        row["alternative_splice_primary_score_ratio_status"] == "unavailable_negative_primary"
        for row in extended["candidate_rows"]
    ) == 2


def test_external_requirement_reporting_keeps_eligibility_populations_distinct() -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(candidate, primary, exact_support=True)
    selection, _, _ = _selection(
        primary,
        [candidate],
        [row],
        config=AlternativeSpliceSelectionConfig.resolve(
            evidence_requirement="protein"
        ),
    )
    extended = extend_alternative_splice_report(
        {
            "candidate_rows": [row],
            "gene_rows": [{"primary_output_transcript_id": primary.output_transcript_id}],
            "run_rows": [{}],
        },
        selection,
        alternatives_emitted=False,
    )
    candidate_row = extended["candidate_rows"][0]
    assert candidate_row[
        "alternative_splice_exact_model_source_requirement_met"
    ] is True
    assert candidate_row[
        "alternative_splice_mosaic_source_requirement_met"
    ] is False
    assert candidate_row["alternative_splice_fixed_rule_eligible"] is True
    assert candidate_row["alternative_splice_support_policy_eligible"] is True
    assert candidate_row["alternative_splice_external_evidence_requirement_met"] is False
    assert candidate_row["alternative_splice_prethreshold_eligible"] is False
    assert candidate_row["alternative_splice_decision"] == (
        "rejected_external_evidence_requirement"
    )
    gene_row = extended["gene_rows"][0]
    assert gene_row["alternative_splice_fixed_rule_eligible_candidates"] == 1
    assert gene_row["alternative_splice_support_policy_eligible_candidates"] == 1
    assert gene_row[
        "alternative_splice_external_evidence_requirement_failure_candidates"
    ] == 1
    assert gene_row["alternative_splice_prethreshold_eligible_candidates"] == 0
    run_row = extended["run_rows"][0]
    assert run_row["alternative_splice_selector_evaluated"] is True
    assert run_row["alternative_splice_alternatives_emitted"] is False
    assert run_row["alternative_splice_fixed_rule_eligible_candidate_count"] == 1
    assert run_row["alternative_splice_support_policy_eligible_candidate_count"] == 1
    assert run_row[
        "alternative_splice_external_evidence_requirement_failure_candidate_count"
    ] == 1
    assert run_row["alternative_splice_prethreshold_eligible_candidate_count"] == 0
    assert run_row["alternative_splice_prediction_source_count"] == 2
    assert run_row["alternative_splice_min_model_support_sources"] == 2
    assert run_row[
        "alternative_splice_min_model_support_source_fraction"
    ] == "."
    assert run_row["alternative_splice_effective_min_model_support_sources"] == 2


def test_threshold_boundaries_are_inclusive_and_epsilon_below_fails() -> None:
    primary, candidate = _distinct_pair()
    base = _candidate_row(candidate, primary, candidate_raw=0.7, candidate_scaled=0.6)
    actual_overlap = cds_union_overlap_metrics(primary, candidate).reciprocal_fraction
    config = AlternativeSpliceSelectionConfig.resolve(
        min_reciprocal_cds_overlap=actual_overlap
    )
    accepted, _, _ = _selection(primary, [candidate], [base], config=config)
    assert accepted.candidate_decisions[0].stage2b_decision == "accepted"

    below_overlap_config = AlternativeSpliceSelectionConfig.resolve(
        min_reciprocal_cds_overlap=actual_overlap + 1e-8
    )
    result, _, _ = _selection(
        primary, [candidate], [base], config=below_overlap_config
    )
    assert result.candidate_decisions[0].stage2b_decision == "rejected_reciprocal_overlap"

    for field, expected, threshold in (
        ("same_frame_overlap_fraction", "rejected_frame_concordance", 0.95),
        ("candidate_model_selection_score_scaled", "rejected_absolute_model_selection_score", 0.60),
    ):
        row = dict(base)
        row[field] = threshold - 1e-8
        result, _, _ = _selection(primary, [candidate], [row], config=config)
        assert result.candidate_decisions[0].stage2b_decision == expected
    row = dict(base)
    row["candidate_model_selection_score"] = 0.7 - 1e-8
    row["candidate_primary_model_selection_score_ratio"] = 0.7 - 1e-8
    result, _, _ = _selection(primary, [candidate], [row], config=config)
    assert result.candidate_decisions[0].stage2b_decision == "rejected_primary_score_ratio"


def test_exact_plus_mosaic_is_independent_of_threshold_preset() -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(
        candidate, primary, exact_support=False, full_mosaic=True
    )
    exact_result, _, _ = _selection(
        primary,
        [candidate],
        [row],
        config=AlternativeSpliceSelectionConfig.resolve(support_policy="exact"),
    )
    exact_decision = exact_result.candidate_decisions[0]
    assert exact_decision.stage2b_decision == "rejected_no_qualifying_support"
    assert exact_decision.qualifying_support_routes == ()
    assert exact_decision.factual_support_routes == ("full_structural_mosaic",)
    mosaic_result, _, _ = _selection(
        primary,
        [candidate],
        [row],
        config=AlternativeSpliceSelectionConfig.resolve(
            support_policy="exact-plus-mosaic"
        ),
    )
    decision = mosaic_result.candidate_decisions[0]
    assert decision.stage2b_decision == "accepted"
    assert decision.mosaic_policy_required is True
    assert decision.qualifying_support_routes == ("full_structural_mosaic",)
    assert decision.factual_support_routes == ("full_structural_mosaic",)
    assert decision.stage2b_frame_concordance_status == "calculated"
    assert decision.stage2b_frame_concordance == pytest.approx(0.95)


def test_greedy_order_conflict_cap_and_public_numbering() -> None:
    primary = _model(
        "primary", [(0, 90), (180, 270), (360, 450)], source="P"
    )
    first = _model("first", [(0, 90), (180, 240), (330, 450)], source="A")
    equivalent = _model(
        "equivalent", [(0, 90), (180, 240), (330, 450)], source="B"
    )
    later = _model("later", [(0, 90), (180, 246), (336, 450)], source="C")
    beyond_cap = _model(
        "beyond_cap", [(0, 90), (180, 252), (342, 450)], source="D"
    )
    candidates = [first, equivalent, later, beyond_cap]
    rows = [
        _candidate_row(model, primary, candidate_raw=score, candidate_scaled=score, rank=rank)
        for model, score, rank in zip(candidates, (0.95, 0.90, 0.85, 0.80), (1, 2, 3, 4))
    ]
    result, annotation, score_rows = _selection(
        primary,
        candidates,
        rows,
        config=AlternativeSpliceSelectionConfig.resolve(max_alternatives_per_gene=2),
    )
    by_id = {
        decision.candidate_internal_transcript_id: decision
        for decision in result.candidate_decisions
    }
    assert by_id["first"].stage2b_decision == "accepted"
    assert by_id["first"].accepted_public_transcript_id == "MoSelecta_000001.2"
    assert by_id["equivalent"].stage2b_decision == "rejected_greedy_exact_equivalent"
    assert by_id["equivalent"].greedy_conflict_internal_transcript_id == "first"
    assert by_id["equivalent"].greedy_conflict_public_transcript_id == "MoSelecta_000001.2"
    assert by_id["later"].stage2b_decision == "accepted"
    assert by_id["later"].accepted_public_transcript_id == "MoSelecta_000001.3"
    assert by_id["beyond_cap"].stage2b_decision == "rejected_max_alternatives_per_gene"

    augmented = augment_annotation(annotation, result, score_rows)
    gene = augmented.genes[0]
    assert gene.gene_id == annotation.genes[0].gene_id
    assert gene.primary_transcript is annotation.genes[0].primary_transcript
    assert [transcript.transcript_id for transcript in gene.transcripts] == [
        "MoSelecta_000001.1",
        "MoSelecta_000001.2",
        "MoSelecta_000001.3",
    ]
    assert [transcript.transcript_role for transcript in gene.transcripts] == [
        "primary", "alternative", "alternative"
    ]
    assert gene.transcripts[1].source_model is first
    assert gene.transcripts[1].source_model.source_name == "A"
    assert gene.transcripts[1].annotation_summary.model_selection_score_scaled == 0.95
    assert gene.start0 == min(model.model_start0 for model in (primary, first, later))
    assert gene.end0 == max(model.model_end0 for model in (primary, first, later))


def test_active_unmodified_preset_matches_corresponding_canonical_comparison() -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(candidate, primary)
    for preset in AlternativeSplicePreset:
        result, _, _ = _selection(
            primary,
            [candidate],
            [row],
            config=AlternativeSpliceSelectionConfig.resolve(preset=preset),
        )
        active = {
            decision.candidate_internal_transcript_id
            for decision in result.candidate_decisions
            if decision.stage2b_decision == "accepted"
        }
        comparison = next(
            item for item in result.canonical_comparisons if item.preset is preset
        )
        assert active == set(comparison.accepted_internal_transcript_ids)
        assert all(
            accepted.accepted_public_transcript_id.endswith(".2")
            for accepted in result.accepted_alternatives
        )
        assert all(
            not hasattr(comparison, "accepted_public_transcript_ids")
            for _ in (0,)
        )


def test_cap_allows_fifty_and_reports_every_later_candidate() -> None:
    primary = _model(
        "primary", [(0, 100), (200, 300), (400, 500)], source="P"
    )
    candidates = [
        _model(
            f"candidate_{index:02d}",
            [
                (0, 100),
                (200, 275 + (index if index < 25 else index + 1)),
                (375 + (index if index < 25 else index + 1), 500),
            ],
            source=f"S{index:02d}",
        )
        for index in range(52)
    ]
    rows = [
        _candidate_row(
            candidate,
            primary,
            candidate_raw=0.9 - index / 1000,
            candidate_scaled=0.9 - index / 1000,
            rank=index + 1,
        )
        for index, candidate in enumerate(candidates)
    ]
    result, _, _ = _selection(primary, candidates, rows)
    decisions = sorted(
        result.candidate_decisions,
        key=lambda item: item.alternative_candidate_order_rank or 0,
    )
    assert sum(item.stage2b_decision == "accepted" for item in decisions) == 50
    assert decisions[49].accepted_alternative_ordinal == 51
    assert decisions[50].stage2b_decision == "rejected_max_alternatives_per_gene"
    assert decisions[51].stage2b_decision == "rejected_max_alternatives_per_gene"


def test_native_greedy_preset_nonmonotonicity_is_reported_not_closed() -> None:
    primary = _model(
        "primary", [(105, 310), (390, 595)], source="P"
    )
    permissive_blocker = _model(
        "blocker",
        [(0, 40), (100, 300), (400, 600), (700, 740)],
        source="A",
    )
    stricter_survivor = _model(
        "survivor", [(105, 300), (400, 595)], source="B"
    )
    blocker_row = _candidate_row(
        permissive_blocker,
        primary,
        candidate_raw=0.90,
        candidate_scaled=0.90,
        rank=1,
    )
    survivor_row = _candidate_row(
        stricter_survivor,
        primary,
        candidate_raw=0.85,
        candidate_scaled=0.85,
        rank=2,
    )
    result, _, _ = _selection(
        primary,
        [permissive_blocker, stricter_survivor],
        [blocker_row, survivor_row],
    )
    comparisons = {
        comparison.preset: set(comparison.accepted_internal_transcript_ids)
        for comparison in result.canonical_comparisons
    }
    assert comparisons[AlternativeSplicePreset.PERMISSIVE] == {"blocker"}
    assert comparisons[AlternativeSplicePreset.INTERMEDIATE] == {"survivor"}
    assert comparisons[AlternativeSplicePreset.STRINGENT] == {"survivor"}
    assert {
        item.candidate_model.output_transcript_id
        for item in result.accepted_alternatives
    } == {"survivor"}
    assert len(result.preset_nonmonotonicities) == 1
    nonmonotonicity = result.preset_nonmonotonicities[0]
    assert nonmonotonicity.candidate_internal_transcript_id == "survivor"
    assert nonmonotonicity.looser_decision == "rejected_greedy_containment"
    assert nonmonotonicity.looser_conflict_internal_transcript_id == "blocker"
    extended = extend_alternative_splice_report(
        {
            "candidate_rows": [blocker_row, survivor_row],
            "gene_rows": [],
            "run_rows": [{}],
        },
        result,
        alternatives_emitted=True,
    )
    run_row = extended["run_rows"][0]
    details = str(
        run_row["alternative_splice_preset_nonmonotonicity_details"]
    ).split(";")
    assert run_row["alternative_splice_preset_nonmonotonicity_count"] == 1
    assert len(details) == len(set(details)) == 1
    assert details[0].startswith("permissive->intermediate:survivor:")


def test_zero_preset_nonmonotonicity_reporting_is_clean() -> None:
    primary, candidate = _distinct_pair()
    row = _candidate_row(candidate, primary, candidate_raw=0.9, candidate_scaled=0.9)
    result, _, _ = _selection(primary, [candidate], [row])
    assert result.preset_nonmonotonicities == ()
    extended = extend_alternative_splice_report(
        {
            "candidate_rows": [row],
            "gene_rows": [],
            "run_rows": [{}],
        },
        result,
        alternatives_emitted=True,
    )
    run_row = extended["run_rows"][0]
    assert run_row["alternative_splice_preset_nonmonotonicity_count"] == 0
    assert run_row["alternative_splice_preset_nonmonotonicity_details"] == "."


def test_intermediate_to_stringent_nonmonotonicity_is_reported_once() -> None:
    primary = _model("primary", [(105, 310), (390, 595)], source="P")
    intermediate_blocker = _model(
        "blocker",
        [(70, 90), (100, 300), (400, 600), (610, 630)],
        source="A",
    )
    stringent_survivor = _model(
        "survivor", [(105, 300), (400, 595)], source="B"
    )
    blocker_row = _candidate_row(
        intermediate_blocker,
        primary,
        candidate_raw=0.75,
        candidate_scaled=0.90,
        rank=1,
    )
    survivor_row = _candidate_row(
        stringent_survivor,
        primary,
        candidate_raw=0.85,
        candidate_scaled=0.90,
        rank=2,
    )
    result, _, _ = _selection(
        primary,
        [intermediate_blocker, stringent_survivor],
        [blocker_row, survivor_row],
    )
    comparisons = {
        comparison.preset: set(comparison.accepted_internal_transcript_ids)
        for comparison in result.canonical_comparisons
    }
    assert comparisons[AlternativeSplicePreset.PERMISSIVE] == {"blocker"}
    assert comparisons[AlternativeSplicePreset.INTERMEDIATE] == {"blocker"}
    assert comparisons[AlternativeSplicePreset.STRINGENT] == {"survivor"}
    assert {
        item.candidate_model.output_transcript_id
        for item in result.accepted_alternatives
    } == {"blocker"}
    assert len(result.preset_nonmonotonicities) == 1
    case = result.preset_nonmonotonicities[0]
    assert (
        case.looser_preset,
        case.stricter_preset,
        case.candidate_internal_transcript_id,
    ) == (
        AlternativeSplicePreset.INTERMEDIATE,
        AlternativeSplicePreset.STRINGENT,
        "survivor",
    )
    extended = extend_alternative_splice_report(
        {
            "candidate_rows": [blocker_row, survivor_row],
            "gene_rows": [],
            "run_rows": [{}],
        },
        result,
        alternatives_emitted=True,
    )
    run_row = extended["run_rows"][0]
    assert run_row["alternative_splice_preset_nonmonotonicity_count"] == 1
    details = str(
        run_row["alternative_splice_preset_nonmonotonicity_details"]
    ).split(";")
    assert len(details) == len(set(details)) == 1
    assert details[0].startswith("intermediate->stringent:survivor:")
