"""Stage 2B noncanonical alternative-splice selection.

The selector consumes the biological facts produced by ``alternative_splice``
after the retained primary population has been frozen. It never mutates a
``CodingModel``, the retained-primary set, or the canonical final annotation;
Stage 2C may only serialize its existing active result.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from decimal import Decimal, ROUND_CEILING
from enum import Enum
from itertools import pairwise
from types import MappingProxyType

from .final_annotation import FinalAnnotation, public_transcript_id
from .models import CodingModel
from .structural import (
    cds_exon_chains_match,
    cds_intron_chains_match,
)


class AlternativeSplicePreset(str, Enum):
    """Named canonical Stage 2B numerical threshold presets."""

    PERMISSIVE = "permissive"
    INTERMEDIATE = "intermediate"
    STRINGENT = "stringent"


class AlternativeSpliceSupportPolicy(str, Enum):
    """Biological support routes allowed to qualify a candidate."""

    EXACT = "exact"
    EXACT_PLUS_MOSAIC = "exact-plus-mosaic"


class AlternativeSpliceEvidenceRequirement(str, Enum):
    """Optional qualifying exact external-evidence requirement."""

    NONE = "none"
    PROTEIN_OR_TRANSCRIPT = "protein-or-transcript"
    PROTEIN = "protein"
    TRANSCRIPT = "transcript"
    PROTEIN_AND_TRANSCRIPT = "protein-and-transcript"


ADJACENT_CANONICAL_PRESET_COMPARISONS = tuple(pairwise(AlternativeSplicePreset))


@dataclass(frozen=True, slots=True)
class AlternativeSpliceThresholds:
    """The four inclusive numerical Stage 2B eligibility thresholds."""

    min_reciprocal_cds_overlap: float
    min_frame_concordance: float
    min_model_selection_score_scaled: float
    min_primary_score_ratio: float

    def __post_init__(self) -> None:
        for name, value in (
            ("min_reciprocal_cds_overlap", self.min_reciprocal_cds_overlap),
            ("min_frame_concordance", self.min_frame_concordance),
            ("min_model_selection_score_scaled", self.min_model_selection_score_scaled),
            ("min_primary_score_ratio", self.min_primary_score_ratio),
        ):
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be a finite value between 0 and 1 inclusive")


CANONICAL_ALTERNATIVE_SPLICE_PRESETS: Mapping[
    AlternativeSplicePreset, AlternativeSpliceThresholds
] = MappingProxyType({
    AlternativeSplicePreset.PERMISSIVE: AlternativeSpliceThresholds(0.80, 0.95, 0.50, 0.70),
    AlternativeSplicePreset.INTERMEDIATE: AlternativeSpliceThresholds(0.85, 0.95, 0.60, 0.70),
    AlternativeSplicePreset.STRINGENT: AlternativeSpliceThresholds(0.85, 0.95, 0.60, 0.80),
})


def effective_min_model_support_sources(
    *,
    prediction_source_count: int,
    configured_minimum: int = 2,
    configured_fraction: float | None = None,
) -> int:
    """Resolve the absolute/fractional prediction-source requirement."""
    if prediction_source_count < 0:
        raise ValueError("prediction_source_count must be >= 0")
    if configured_minimum < 2:
        raise ValueError("min_model_support_sources must be >= 2")
    if configured_fraction is not None and (
        not math.isfinite(configured_fraction)
        or not 0.0 < configured_fraction <= 1.0
    ):
        raise ValueError(
            "min_model_support_source_fraction must be finite and > 0 and <= 1"
        )
    fraction_minimum = 0
    if configured_fraction is not None:
        fraction_minimum = int(
            (Decimal(prediction_source_count) * Decimal(str(configured_fraction)))
            .to_integral_value(rounding=ROUND_CEILING)
        )
    return max(2, configured_minimum, fraction_minimum)


@dataclass(frozen=True, slots=True)
class AlternativeSpliceSelectionConfig:
    """Resolved active Stage 2B settings and their canonical provenance."""

    requested_preset: AlternativeSplicePreset
    canonical_thresholds: AlternativeSpliceThresholds
    effective_thresholds: AlternativeSpliceThresholds
    support_policy: AlternativeSpliceSupportPolicy
    evidence_requirement: AlternativeSpliceEvidenceRequirement
    max_alternatives_per_gene: int = 50
    min_model_support_sources: int = 2
    min_model_support_source_fraction: float | None = None
    prediction_source_count: int | None = None
    effective_min_model_support_sources: int | None = None
    overridden_thresholds: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.max_alternatives_per_gene < 0:
            raise ValueError("max_alternatives_per_gene must be >= 0")
        if self.min_model_support_sources < 2:
            raise ValueError("min_model_support_sources must be >= 2")
        if self.min_model_support_source_fraction is not None and (
            not math.isfinite(self.min_model_support_source_fraction)
            or not 0.0 < self.min_model_support_source_fraction <= 1.0
        ):
            raise ValueError(
                "min_model_support_source_fraction must be finite and > 0 and <= 1"
            )
        if self.prediction_source_count is None:
            if self.effective_min_model_support_sources is not None:
                raise ValueError(
                    "effective model-source requirement needs prediction_source_count"
                )
        else:
            expected = effective_min_model_support_sources(
                prediction_source_count=self.prediction_source_count,
                configured_minimum=self.min_model_support_sources,
                configured_fraction=self.min_model_support_source_fraction,
            )
            if self.effective_min_model_support_sources != expected:
                raise ValueError(
                    "effective_min_model_support_sources does not match configured values"
                )

    @classmethod
    def resolve(
        cls,
        *,
        preset: str | AlternativeSplicePreset = AlternativeSplicePreset.INTERMEDIATE,
        support_policy: str | AlternativeSpliceSupportPolicy = AlternativeSpliceSupportPolicy.EXACT,
        evidence_requirement: str | AlternativeSpliceEvidenceRequirement = AlternativeSpliceEvidenceRequirement.NONE,
        max_alternatives_per_gene: int = 50,
        min_model_support_sources: int = 2,
        min_model_support_source_fraction: float | None = None,
        min_reciprocal_cds_overlap: float | None = None,
        min_frame_concordance: float | None = None,
        min_model_selection_score_scaled: float | None = None,
        min_primary_score_ratio: float | None = None,
    ) -> AlternativeSpliceSelectionConfig:
        """Resolve a named preset plus independent explicit overrides."""
        requested = AlternativeSplicePreset(preset)
        active_support_policy = AlternativeSpliceSupportPolicy(support_policy)
        active_evidence_requirement = AlternativeSpliceEvidenceRequirement(
            evidence_requirement
        )
        canonical = CANONICAL_ALTERNATIVE_SPLICE_PRESETS[requested]
        values = {
            "min_reciprocal_cds_overlap": min_reciprocal_cds_overlap,
            "min_frame_concordance": min_frame_concordance,
            "min_model_selection_score_scaled": min_model_selection_score_scaled,
            "min_primary_score_ratio": min_primary_score_ratio,
        }
        overridden = tuple(name for name, value in values.items() if value is not None)
        effective = AlternativeSpliceThresholds(
            min_reciprocal_cds_overlap=(
                canonical.min_reciprocal_cds_overlap
                if min_reciprocal_cds_overlap is None
                else min_reciprocal_cds_overlap
            ),
            min_frame_concordance=(
                canonical.min_frame_concordance
                if min_frame_concordance is None
                else min_frame_concordance
            ),
            min_model_selection_score_scaled=(
                canonical.min_model_selection_score_scaled
                if min_model_selection_score_scaled is None
                else min_model_selection_score_scaled
            ),
            min_primary_score_ratio=(
                canonical.min_primary_score_ratio
                if min_primary_score_ratio is None
                else min_primary_score_ratio
            ),
        )
        return cls(
            requested_preset=requested,
            canonical_thresholds=canonical,
            effective_thresholds=effective,
            support_policy=active_support_policy,
            evidence_requirement=active_evidence_requirement,
            max_alternatives_per_gene=max_alternatives_per_gene,
            min_model_support_sources=min_model_support_sources,
            min_model_support_source_fraction=min_model_support_source_fraction,
            overridden_thresholds=overridden,
        )

    def for_prediction_source_count(
        self, prediction_source_count: int
    ) -> AlternativeSpliceSelectionConfig:
        """Return this configuration with its dataset-level requirement resolved."""
        effective = effective_min_model_support_sources(
            prediction_source_count=prediction_source_count,
            configured_minimum=self.min_model_support_sources,
            configured_fraction=self.min_model_support_source_fraction,
        )
        return replace(
            self,
            prediction_source_count=prediction_source_count,
            effective_min_model_support_sources=effective,
        )

    @property
    def is_unmodified_canonical_preset(self) -> bool:
        return not self.overridden_thresholds


@dataclass(frozen=True, slots=True)
class CdsBaseContainment:
    """Diagnostic subset relation between two unions of half-open CDS intervals."""

    candidate_in_retained: bool
    retained_in_candidate: bool

    @property
    def bidirectional(self) -> bool:
        return self.candidate_in_retained and self.retained_in_candidate

    @property
    def any_direction(self) -> bool:
        return self.candidate_in_retained or self.retained_in_candidate


@dataclass(frozen=True, slots=True)
class RelaxedCdsChainContainment:
    """Directional relaxed contiguous CDS-exon-chain containment."""

    candidate_in_retained: bool
    retained_in_candidate: bool

    @property
    def bidirectional(self) -> bool:
        return self.candidate_in_retained and self.retained_in_candidate

    @property
    def any_direction(self) -> bool:
        return self.candidate_in_retained or self.retained_in_candidate


@dataclass(frozen=True, slots=True)
class CdsUnionOverlap:
    """Shared-base metrics calculated from disjoint CDS interval unions."""

    shared_bp: int
    retained_fraction: float
    candidate_fraction: float
    reciprocal_fraction: float


@dataclass(frozen=True, slots=True)
class AlternativeCandidateDecision:
    """One immutable active or canonical-preset candidate assessment."""

    assessment_id: str
    candidate_internal_transcript_id: str
    primary_internal_transcript_id: str | None
    alternative_candidate_order_rank: int | None
    stage2b_decision: str
    stage2b_rejection_reason: str
    fixed_eligible: bool
    support_policy_eligible: bool
    external_evidence_requirement_met: bool
    external_evidence_requirement_routes: tuple[str, ...]
    prethreshold_eligible: bool
    preset_threshold_survivor: bool
    stage2b_primary_score_ratio_status: str
    stage2b_candidate_primary_model_selection_score_ratio: float | None
    stage2b_cds_overlap_bp: int | None
    stage2b_primary_cds_overlap_fraction: float | None
    stage2b_candidate_cds_overlap_fraction: float | None
    stage2b_reciprocal_cds_overlap: float | None
    stage2b_frame_concordance_status: str
    stage2b_frame_concordance: float | None
    candidate_cds_bases_contained_in_primary: bool
    primary_cds_bases_contained_in_candidate: bool
    bidirectional_cds_base_containment: bool
    candidate_relaxed_cds_chain_contained_in_primary: bool
    primary_relaxed_cds_chain_contained_in_candidate: bool
    bidirectional_relaxed_cds_chain_containment: bool
    mosaic_components_complete: bool
    mosaic_additional_corroboration: bool
    exact_model_source_requirement_met: bool
    mosaic_source_requirement_met: bool
    mosaic_policy_required: bool
    factual_support_routes: tuple[str, ...]
    qualifying_support_routes: tuple[str, ...]
    greedy_conflict_public_transcript_id: str | None = None
    greedy_conflict_internal_transcript_id: str | None = None
    greedy_conflict_type: str | None = None
    accepted_alternative_ordinal: int | None = None
    accepted_public_transcript_id: str | None = None


@dataclass(frozen=True, slots=True)
class AcceptedAlternative:
    """The actual representative model accepted for one frozen public gene."""

    primary_internal_transcript_id: str
    public_gene_id: str
    candidate_model: CodingModel
    alternative_candidate_order_rank: int
    accepted_alternative_ordinal: int
    accepted_public_transcript_id: str
    factual_support_routes: tuple[str, ...]
    qualifying_support_routes: tuple[str, ...]
    evidence_requirement: str
    mosaic_policy_required: bool
    stage2b_reciprocal_cds_overlap: float
    stage2b_frame_concordance: float
    stage2b_primary_score_ratio: float


@dataclass(frozen=True, slots=True)
class AlternativeGeneDecision:
    """Active accepted representatives for one retained primary."""

    primary_internal_transcript_id: str
    public_gene_id: str
    accepted_alternatives: tuple[AcceptedAlternative, ...]


@dataclass(frozen=True, slots=True)
class AlternativePresetComparison:
    """Canonical comparison result without public alternative numbering."""

    preset: AlternativeSplicePreset
    accepted_internal_transcript_ids: tuple[str, ...]
    decisions_by_assessment_id: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class AlternativePresetNonmonotonicity:
    """A candidate accepted by a stricter preset but not a looser preset."""

    looser_preset: AlternativeSplicePreset
    stricter_preset: AlternativeSplicePreset
    candidate_internal_transcript_id: str
    looser_decision: str
    looser_conflict_internal_transcript_id: str | None
    looser_conflict_type: str | None


@dataclass(frozen=True, slots=True)
class AlternativeSpliceSelectionResult:
    """Active Stage 2B decisions plus independent canonical comparisons."""

    config: AlternativeSpliceSelectionConfig
    candidate_decisions: tuple[AlternativeCandidateDecision, ...]
    gene_decisions: tuple[AlternativeGeneDecision, ...]
    accepted_alternatives: tuple[AcceptedAlternative, ...]
    canonical_comparisons: tuple[AlternativePresetComparison, ...]
    preset_nonmonotonicities: tuple[AlternativePresetNonmonotonicity, ...]

    def decisions_by_assessment_id(self) -> dict[str, AlternativeCandidateDecision]:
        return {decision.assessment_id: decision for decision in self.candidate_decisions}

    def accepted_by_primary_internal_id(self) -> dict[str, tuple[AcceptedAlternative, ...]]:
        return {
            gene.primary_internal_transcript_id: gene.accepted_alternatives
            for gene in self.gene_decisions
        }


def _diagnose_preset_nonmonotonicities(
    comparison_evaluations: Mapping[AlternativeSplicePreset, _Evaluation],
) -> tuple[AlternativePresetNonmonotonicity, ...]:
    """Report each native adjacent-preset nonmonotonicity exactly once."""
    nonmonotonicities: list[AlternativePresetNonmonotonicity] = []
    seen_cases: set[tuple[AlternativeSplicePreset, AlternativeSplicePreset, str]] = set()
    for looser, stricter in ADJACENT_CANONICAL_PRESET_COMPARISONS:
        looser_eval = comparison_evaluations[looser]
        stricter_eval = comparison_evaluations[stricter]
        looser_accepted = {
            decision.candidate_internal_transcript_id
            for decision in looser_eval.decisions
            if decision.stage2b_decision == "accepted"
        }
        stricter_accepted = {
            decision.candidate_internal_transcript_id
            for decision in stricter_eval.decisions
            if decision.stage2b_decision == "accepted"
        }
        looser_by_candidate = {
            decision.candidate_internal_transcript_id: decision
            for decision in looser_eval.decisions
        }
        for candidate_id in sorted(stricter_accepted - looser_accepted):
            case_key = (looser, stricter, candidate_id)
            if case_key in seen_cases:
                continue
            seen_cases.add(case_key)
            looser_decision = looser_by_candidate[candidate_id]
            nonmonotonicities.append(AlternativePresetNonmonotonicity(
                looser_preset=looser,
                stricter_preset=stricter,
                candidate_internal_transcript_id=candidate_id,
                looser_decision=looser_decision.stage2b_decision,
                looser_conflict_internal_transcript_id=(
                    looser_decision.greedy_conflict_internal_transcript_id
                ),
                looser_conflict_type=looser_decision.greedy_conflict_type,
            ))
    return tuple(nonmonotonicities)


STAGE2B_CANDIDATE_COLUMNS = [
    "alternative_splice_preset_requested",
    "alternative_splice_canonical_preset",
    "alternative_splice_support_policy",
    "alternative_splice_evidence_requirement",
    "alternative_splice_max_alternatives_per_gene",
    "alternative_splice_prediction_source_count",
    "alternative_splice_min_model_support_sources",
    "alternative_splice_min_model_support_source_fraction",
    "alternative_splice_effective_min_model_support_sources",
    "canonical_min_reciprocal_cds_overlap",
    "canonical_min_frame_concordance",
    "canonical_min_model_selection_score_scaled",
    "canonical_min_primary_score_ratio",
    "effective_min_reciprocal_cds_overlap",
    "effective_min_frame_concordance",
    "effective_min_model_selection_score_scaled",
    "effective_min_primary_score_ratio",
    "override_min_reciprocal_cds_overlap",
    "override_min_frame_concordance",
    "override_min_model_selection_score_scaled",
    "override_min_primary_score_ratio",
    "alternative_splice_candidate_cds_bases_contained_in_primary",
    "alternative_splice_primary_cds_bases_contained_in_candidate",
    "alternative_splice_bidirectional_cds_base_containment",
    "alternative_splice_candidate_relaxed_cds_chain_contained_in_primary",
    "alternative_splice_primary_relaxed_cds_chain_contained_in_candidate",
    "alternative_splice_bidirectional_relaxed_cds_chain_containment",
    "alternative_splice_cds_overlap_bp",
    "alternative_splice_primary_cds_overlap_fraction",
    "alternative_splice_candidate_cds_overlap_fraction",
    "alternative_splice_reciprocal_cds_overlap",
    "alternative_splice_primary_score_ratio_status",
    "alternative_splice_candidate_primary_model_selection_score_ratio",
    "mosaic_components_complete",
    "mosaic_additional_corroboration",
    "alternative_splice_exact_model_source_requirement_met",
    "alternative_splice_mosaic_source_requirement_met",
    "mosaic_policy_required",
    "alternative_splice_qualifying_support_routes",
    "alternative_splice_fixed_rule_eligible",
    "alternative_splice_support_policy_eligible",
    "alternative_splice_external_evidence_requirement_met",
    "alternative_splice_external_evidence_requirement_routes",
    "alternative_splice_prethreshold_eligible",
    "alternative_splice_numerical_threshold_survivor",
    "alternative_splice_decision",
    "alternative_splice_rejection_reason",
    "greedy_conflict_public_transcript_id",
    "greedy_conflict_internal_transcript_id",
    "greedy_conflict_type",
    "accepted_alternative_ordinal",
    "accepted_public_transcript_id",
    "canonical_permissive_alternative_splice_decision",
    "canonical_intermediate_alternative_splice_decision",
    "canonical_stringent_alternative_splice_decision",
]


STAGE2B_GENE_COLUMNS = [
    "alternative_splice_fixed_rule_eligible_candidates",
    "alternative_splice_support_policy_eligible_candidates",
    "alternative_splice_external_evidence_requirement_failure_candidates",
    "alternative_splice_prethreshold_eligible_candidates",
    "alternative_splice_numerical_threshold_survivors",
    "alternative_splice_greedy_accepted_alternatives",
    "alternative_splice_greedy_redundancy_rejections",
    "alternative_splice_cap_rejections",
    "alternative_splice_max_alternatives_per_gene",
    "accepted_alternative_count",
    "max_candidate_rank_accepted",
    "accepted_internal_transcript_ids",
    "canonical_permissive_accepted_count",
    "canonical_permissive_accepted_internal_transcript_ids",
    "canonical_intermediate_accepted_count",
    "canonical_intermediate_accepted_internal_transcript_ids",
    "canonical_stringent_accepted_count",
    "canonical_stringent_accepted_internal_transcript_ids",
]


STAGE2B_RUN_COLUMNS = [
    "alternative_splice_selector_evaluated",
    "alternative_splice_alternatives_emitted",
    "alternative_splice_preset_requested",
    "alternative_splice_canonical_preset",
    "alternative_splice_support_policy",
    "alternative_splice_evidence_requirement",
    "alternative_splice_max_alternatives_per_gene",
    "alternative_splice_prediction_source_count",
    "alternative_splice_min_model_support_sources",
    "alternative_splice_min_model_support_source_fraction",
    "alternative_splice_effective_min_model_support_sources",
    "canonical_min_reciprocal_cds_overlap",
    "canonical_min_frame_concordance",
    "canonical_min_model_selection_score_scaled",
    "canonical_min_primary_score_ratio",
    "effective_min_reciprocal_cds_overlap",
    "effective_min_frame_concordance",
    "effective_min_model_selection_score_scaled",
    "effective_min_primary_score_ratio",
    "override_min_reciprocal_cds_overlap",
    "override_min_frame_concordance",
    "override_min_model_selection_score_scaled",
    "override_min_primary_score_ratio",
    "alternative_splice_fixed_rule_eligible_candidate_count",
    "alternative_splice_support_policy_eligible_candidate_count",
    "alternative_splice_external_evidence_requirement_failure_candidate_count",
    "alternative_splice_prethreshold_eligible_candidate_count",
    "alternative_splice_numerical_threshold_survivor_count",
    "alternative_splice_genes_with_accepted_alternatives",
    "alternative_splice_total_accepted_alternatives",
    "alternative_splice_max_accepted_alternatives_per_gene",
    "alternative_splice_exact_model_support_accepted_count",
    "alternative_splice_exact_protein_support_accepted_count",
    "alternative_splice_exact_transcript_support_accepted_count",
    "alternative_splice_mosaic_only_accepted_count",
    "alternative_splice_accepted_support_route_composition",
    "alternative_splice_cap_rejected_count",
    "alternative_splice_greedy_exact_conflict_count",
    "alternative_splice_greedy_same_intron_chain_conflict_count",
    "alternative_splice_greedy_containment_conflict_count",
    "alternative_splice_greedy_conflict_type_counts",
    "alternative_splice_primary_score_ratio_unavailable_zero_count",
    "alternative_splice_primary_score_ratio_unavailable_negative_count",
    "alternative_splice_primary_score_ratio_unavailable_zero_assessment_count",
    "alternative_splice_primary_score_ratio_unavailable_negative_assessment_count",
    "canonical_permissive_accepted_count",
    "canonical_permissive_accepted_internal_transcript_ids",
    "canonical_intermediate_accepted_count",
    "canonical_intermediate_accepted_internal_transcript_ids",
    "canonical_stringent_accepted_count",
    "canonical_stringent_accepted_internal_transcript_ids",
    "alternative_splice_preset_nonmonotonicity_count",
    "alternative_splice_preset_nonmonotonicity_details",
]


def _truth(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "passed"}


def _optional_float(value: object) -> float | None:
    if value in {None, ".", ""}:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _optional_int(value: object) -> int | None:
    numeric = _optional_float(value)
    return None if numeric is None else int(numeric)


def _csv(values: Iterable[str]) -> str:
    ordered = tuple(sorted(set(values)))
    return ",".join(ordered) if ordered else "."


def union_cds_intervals(intervals: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Return a sorted disjoint union, robust to overlaps and duplicates."""
    ordered = sorted((int(start0), int(end0)) for start0, end0 in intervals if end0 > start0)
    merged: list[list[int]] = []
    for start0, end0 in ordered:
        if not merged or start0 > merged[-1][1]:
            merged.append([start0, end0])
        else:
            merged[-1][1] = max(merged[-1][1], end0)
    return tuple((start0, end0) for start0, end0 in merged)


def cds_interval_union_is_subset(
    subset_intervals: Iterable[tuple[int, int]],
    superset_intervals: Iterable[tuple[int, int]],
) -> bool:
    """Return whether every base in one CDS union occurs in the other union."""
    subset = union_cds_intervals(subset_intervals)
    superset = union_cds_intervals(superset_intervals)
    if not subset:
        return True
    if not superset:
        return False
    super_index = 0
    for sub_start, sub_end in subset:
        while super_index < len(superset) and superset[super_index][1] <= sub_start:
            super_index += 1
        if super_index >= len(superset):
            return False
        super_start, super_end = superset[super_index]
        if super_start > sub_start or super_end < sub_end:
            return False
    return True


def cds_union_overlap_metrics(
    retained: CodingModel,
    candidate: CodingModel,
) -> CdsUnionOverlap:
    """Calculate overlap without double-counting duplicate/overlapping CDS."""
    retained_union = union_cds_intervals(retained.cds_parts)
    candidate_union = union_cds_intervals(candidate.cds_parts)
    retained_bp = sum(end0 - start0 for start0, end0 in retained_union)
    candidate_bp = sum(end0 - start0 for start0, end0 in candidate_union)
    retained_index = candidate_index = 0
    shared_bp = 0
    while retained_index < len(retained_union) and candidate_index < len(candidate_union):
        retained_start, retained_end = retained_union[retained_index]
        candidate_start, candidate_end = candidate_union[candidate_index]
        shared_bp += max(
            0, min(retained_end, candidate_end) - max(retained_start, candidate_start)
        )
        if retained_end <= candidate_end:
            retained_index += 1
        else:
            candidate_index += 1
    retained_fraction = shared_bp / retained_bp if retained_bp else 0.0
    candidate_fraction = shared_bp / candidate_bp if candidate_bp else 0.0
    return CdsUnionOverlap(
        shared_bp=shared_bp,
        retained_fraction=retained_fraction,
        candidate_fraction=candidate_fraction,
        reciprocal_fraction=min(retained_fraction, candidate_fraction),
    )


def cds_base_containment_relation(
    candidate: CodingModel,
    retained: CodingModel,
) -> CdsBaseContainment:
    """Compare true CDS-base unions, never genomic spans."""
    if candidate.chrom != retained.chrom or candidate.strand != retained.strand:
        return CdsBaseContainment(False, False)
    return CdsBaseContainment(
        candidate_in_retained=cds_interval_union_is_subset(
            candidate.cds_parts, retained.cds_parts
        ),
        retained_in_candidate=cds_interval_union_is_subset(
            retained.cds_parts, candidate.cds_parts
        ),
    )


def _transcript_order_cds_parts(model: CodingModel) -> tuple[tuple[int, int], ...]:
    parts = tuple(model.sorted_cds_parts)
    return parts if model.strand == "+" else tuple(reversed(parts))


def _cds_junction(
    first: tuple[int, int],
    second: tuple[int, int],
) -> tuple[int, int]:
    """Return the genomic intron interval between transcript-adjacent CDS exons."""
    left, right = sorted((first, second))
    return left[1], right[0]


def relaxed_contiguous_cds_chain_is_contained(
    contained_model: CodingModel,
    containing_model: CodingModel,
) -> bool:
    """Return whether one multi-exon CDS chain is a relaxed contiguous subchain.

    Exons map in strand-aware transcript order to a contiguous slice. Every
    mapped splice junction is exact. Internal exons are exact, while only the
    outer boundary of the first and last mapped exon may be trimmed inside the
    corresponding containing exon. A single CDS exon has no junction-defined
    chain and is therefore never contained by this Stage 2B rule.
    """
    if (
        contained_model.chrom != containing_model.chrom
        or contained_model.strand != containing_model.strand
    ):
        return False
    contained = _transcript_order_cds_parts(contained_model)
    containing = _transcript_order_cds_parts(containing_model)
    if len(contained) < 2 or len(contained) > len(containing):
        return False

    width = len(contained)
    contained_junctions = tuple(
        _cds_junction(contained[index], contained[index + 1])
        for index in range(width - 1)
    )
    for offset in range(len(containing) - width + 1):
        mapped = containing[offset : offset + width]
        mapped_junctions = tuple(
            _cds_junction(mapped[index], mapped[index + 1])
            for index in range(width - 1)
        )
        if contained_junctions != mapped_junctions:
            continue
        if any(
            contained[index] != mapped[index]
            for index in range(1, width - 1)
        ):
            continue
        first_start, first_end = contained[0]
        mapped_first_start, mapped_first_end = mapped[0]
        last_start, last_end = contained[-1]
        mapped_last_start, mapped_last_end = mapped[-1]
        if not (
            mapped_first_start <= first_start
            and first_end <= mapped_first_end
            and mapped_last_start <= last_start
            and last_end <= mapped_last_end
        ):
            continue
        return True
    return False


def relaxed_contiguous_cds_chain_containment_relation(
    candidate: CodingModel,
    retained: CodingModel,
) -> RelaxedCdsChainContainment:
    """Compare relaxed contiguous CDS-chain containment bidirectionally."""
    return RelaxedCdsChainContainment(
        candidate_in_retained=relaxed_contiguous_cds_chain_is_contained(
            candidate, retained
        ),
        retained_in_candidate=relaxed_contiguous_cds_chain_is_contained(
            retained, candidate
        ),
    )


def _mosaic_components_complete(row: Mapping[str, object]) -> bool:
    states = [
        str(row.get(field, "incomplete"))
        for field in (
            "mosaic_cds_introns_status",
            "mosaic_cds_splice_sites_status",
            "mosaic_internal_cds_exons_status",
            "mosaic_terminal_cds_exons_status",
            "mosaic_cds_base_consensus_status",
        )
    ]
    return any(state == "complete" for state in states) and all(
        state in {"complete", "not_applicable"} for state in states
    )


def _ratio_state(row: Mapping[str, object]) -> tuple[str, float | None]:
    primary_score = _optional_float(row.get("primary_model_selection_score"))
    candidate_score = _optional_float(row.get("candidate_model_selection_score"))
    if primary_score is None or candidate_score is None:
        return "unavailable_missing_primary", None
    if primary_score == 0.0:
        return "unavailable_zero_primary", None
    if primary_score < 0.0:
        return "unavailable_negative_primary", None
    return "available", candidate_score / primary_score


def _support_state(
    row: Mapping[str, object],
    policy: AlternativeSpliceSupportPolicy,
    effective_model_source_minimum: int,
) -> tuple[
    bool,
    tuple[str, ...],
    tuple[str, ...],
    bool,
    bool,
    bool,
    bool,
    bool,
]:
    exact_routes: list[str] = []
    exact_model_source_count = _optional_int(
        row.get("alternative_splice_exact_model_source_count")
    )
    if exact_model_source_count is None:
        exact_model_source_count = _optional_int(
            row.get("exact_structural_all_source_count")
        )
    if exact_model_source_count is None:
        exact_model_source_count = (
            2 if _truth(row.get("exact_multisource_model_support")) else 0
        )
    mosaic_minimum_source_count = _optional_int(
        row.get("alternative_splice_mosaic_min_component_source_count")
    )
    if mosaic_minimum_source_count is None:
        mosaic_minimum_source_count = (
            2 if _truth(row.get("full_structural_mosaic")) else 0
        )
    exact_model_requirement_met = (
        exact_model_source_count >= effective_model_source_minimum
    )
    mosaic_source_requirement_met = (
        mosaic_minimum_source_count >= effective_model_source_minimum
    )
    mosaic_route_met = (
        not exact_model_requirement_met and mosaic_source_requirement_met
    )
    if exact_model_requirement_met:
        exact_routes.append("exact_multisource_model_support")
    if _truth(row.get("protein_exact_cds_chain_support")):
        exact_routes.append("exact_protein_support")
    if _truth(row.get("exact_coherent_transcript_support")):
        exact_routes.append("exact_coherent_transcript_support")
    components_complete = _mosaic_components_complete(row)
    factual_routes = list(exact_routes)
    if mosaic_route_met:
        factual_routes.append("full_structural_mosaic")
    mosaic_required = (
        not exact_routes
        and policy is AlternativeSpliceSupportPolicy.EXACT_PLUS_MOSAIC
        and mosaic_route_met
    )
    routes = list(exact_routes)
    if mosaic_required:
        routes.append("full_structural_mosaic")
    eligible = bool(routes)
    mosaic_additional = components_complete and bool(exact_routes)
    return (
        eligible,
        tuple(factual_routes),
        tuple(routes),
        components_complete,
        mosaic_additional,
        exact_model_requirement_met,
        mosaic_source_requirement_met,
        mosaic_required,
    )


def _external_evidence_state(
    row: Mapping[str, object],
    requirement: AlternativeSpliceEvidenceRequirement,
) -> tuple[bool, tuple[str, ...]]:
    """Evaluate the configured requirement from exact external routes only."""
    routes: list[str] = []
    protein = _truth(row.get("protein_exact_cds_chain_support"))
    transcript = _truth(row.get("exact_coherent_transcript_support"))
    if protein:
        routes.append("exact_protein_support")
    if transcript:
        routes.append("exact_coherent_transcript_support")
    met = (
        requirement is AlternativeSpliceEvidenceRequirement.NONE
        or (
            requirement is AlternativeSpliceEvidenceRequirement.PROTEIN_OR_TRANSCRIPT
            and (protein or transcript)
        )
        or (
            requirement is AlternativeSpliceEvidenceRequirement.PROTEIN
            and protein
        )
        or (
            requirement is AlternativeSpliceEvidenceRequirement.TRANSCRIPT
            and transcript
        )
        or (
            requirement is AlternativeSpliceEvidenceRequirement.PROTEIN_AND_TRANSCRIPT
            and protein
            and transcript
        )
    )
    return met, tuple(routes)


def _base_decision(
    row: Mapping[str, object],
    candidate: CodingModel | None,
    primary: CodingModel | None,
    thresholds: AlternativeSpliceThresholds,
    support_policy: AlternativeSpliceSupportPolicy,
    evidence_requirement: AlternativeSpliceEvidenceRequirement,
    effective_model_source_minimum: int,
) -> AlternativeCandidateDecision:
    assessment_id = str(row.get("alternative_candidate_assessment_id", "."))
    candidate_id = str(row.get("candidate_output_transcript_id", "."))
    primary_id_text = str(row.get("primary_output_transcript_id", "."))
    primary_id = None if primary_id_text == "." else primary_id_text
    rank = _optional_int(row.get("alternative_candidate_order_rank"))
    ownership = str(row.get("ownership_status", "zero_owner"))
    (
        support_eligible,
        factual_routes,
        routes,
        mosaic_complete,
        mosaic_additional,
        exact_model_requirement_met,
        mosaic_source_requirement_met,
        mosaic_required,
    ) = _support_state(
        row, support_policy, effective_model_source_minimum
    )
    external_evidence_met, external_evidence_routes = _external_evidence_state(
        row, evidence_requirement
    )
    ratio_status, stage2b_ratio = _ratio_state(row)
    base_containment = (
        cds_base_containment_relation(candidate, primary)
        if candidate is not None and primary is not None
        else CdsBaseContainment(False, False)
    )
    chain_containment = (
        relaxed_contiguous_cds_chain_containment_relation(candidate, primary)
        if candidate is not None and primary is not None
        else RelaxedCdsChainContainment(False, False)
    )
    overlap = (
        cds_union_overlap_metrics(primary, candidate)
        if candidate is not None and primary is not None
        else None
    )

    fixed_eligible = False
    support_policy_eligible = False
    prethreshold_eligible = False
    threshold_survivor = False
    if ownership == "zero_owner":
        decision = "rejected_no_owner"
        reason = "no_structurally_associated_retained_primary"
    elif ownership == "ambiguous_multiple_owners":
        decision = "rejected_ambiguous_ownership"
        reason = "multiple_structurally_associated_retained_primaries"
    elif _truth(row.get("exact_cds_equivalent_to_primary")):
        decision = "rejected_exact_equivalent"
        reason = "wobble_aware_exact_cds_equivalence"
    elif _truth(row.get("same_intron_chain_terminal_variant")):
        decision = "rejected_same_intron_chain_terminal_variant"
        reason = "same_cds_intron_chain_terminal_only_difference"
    elif chain_containment.bidirectional:
        decision = "rejected_bidirectional_relaxed_cds_chain_containment"
        reason = "bidirectional_relaxed_contiguous_cds_chain_containment"
    elif chain_containment.candidate_in_retained:
        decision = "rejected_candidate_contained_in_primary"
        reason = "candidate_relaxed_contiguous_cds_chain_contained_in_primary"
    elif chain_containment.retained_in_candidate:
        decision = "rejected_primary_contained_in_candidate"
        reason = "primary_relaxed_contiguous_cds_chain_contained_in_candidate"
    elif str(row.get("frame_concordance_status", ".")) != "calculated":
        decision = "rejected_frame_unavailable"
        reason = str(row.get("frame_concordance_status", "unavailable"))
    elif _truth(row.get("sequence_integrity_default_exclusion")):
        decision = "rejected_sequence_integrity"
        reason = "alternative_specific_sequence_integrity_exclusion"
    elif _truth(row.get("intrinsic_default_exclusion")):
        decision = "rejected_intrinsic_concern"
        reason = "alternative_specific_intrinsic_exclusion"
    else:
        fixed_eligible = True
        support_policy_eligible = support_eligible
        if not support_policy_eligible:
            decision = "rejected_no_qualifying_support"
            reason = f"no_qualifying_support_under_{support_policy.value}"
        elif not external_evidence_met:
            decision = "rejected_external_evidence_requirement"
            reason = (
                "external_evidence_requirement_"
                f"{evidence_requirement.value}_not_met"
            )
        else:
            prethreshold_eligible = True
            reciprocal = overlap.reciprocal_fraction if overlap is not None else None
            frame = _optional_float(row.get("same_frame_overlap_fraction"))
            score = _optional_float(row.get("candidate_model_selection_score_scaled"))
            if reciprocal is None or reciprocal + 1e-15 < thresholds.min_reciprocal_cds_overlap:
                decision = "rejected_reciprocal_overlap"
                reason = "below_effective_min_reciprocal_cds_overlap"
            elif frame is None or frame + 1e-15 < thresholds.min_frame_concordance:
                decision = "rejected_frame_concordance"
                reason = "below_effective_min_frame_concordance"
            elif score is None or score + 1e-15 < thresholds.min_model_selection_score_scaled:
                decision = "rejected_absolute_model_selection_score"
                reason = "below_effective_min_model_selection_score_scaled"
            elif ratio_status != "available":
                decision = "rejected_primary_score_ratio_unavailable"
                reason = ratio_status
            elif stage2b_ratio is None or stage2b_ratio + 1e-15 < thresholds.min_primary_score_ratio:
                decision = "rejected_primary_score_ratio"
                reason = "below_effective_min_primary_score_ratio"
            else:
                threshold_survivor = True
                decision = "pending_greedy_selection"
                reason = "."

    return AlternativeCandidateDecision(
        assessment_id=assessment_id,
        candidate_internal_transcript_id=candidate_id,
        primary_internal_transcript_id=primary_id,
        alternative_candidate_order_rank=rank,
        stage2b_decision=decision,
        stage2b_rejection_reason=reason,
        fixed_eligible=fixed_eligible,
        support_policy_eligible=support_policy_eligible,
        external_evidence_requirement_met=external_evidence_met,
        external_evidence_requirement_routes=external_evidence_routes,
        prethreshold_eligible=prethreshold_eligible,
        preset_threshold_survivor=threshold_survivor,
        stage2b_primary_score_ratio_status=ratio_status,
        stage2b_candidate_primary_model_selection_score_ratio=stage2b_ratio,
        stage2b_cds_overlap_bp=overlap.shared_bp if overlap is not None else None,
        stage2b_primary_cds_overlap_fraction=(
            overlap.retained_fraction if overlap is not None else None
        ),
        stage2b_candidate_cds_overlap_fraction=(
            overlap.candidate_fraction if overlap is not None else None
        ),
        stage2b_reciprocal_cds_overlap=(
            overlap.reciprocal_fraction if overlap is not None else None
        ),
        stage2b_frame_concordance_status=str(
            row.get("frame_concordance_status", ".")
        ),
        stage2b_frame_concordance=_optional_float(
            row.get("same_frame_overlap_fraction")
        ),
        candidate_cds_bases_contained_in_primary=base_containment.candidate_in_retained,
        primary_cds_bases_contained_in_candidate=base_containment.retained_in_candidate,
        bidirectional_cds_base_containment=base_containment.bidirectional,
        candidate_relaxed_cds_chain_contained_in_primary=(
            chain_containment.candidate_in_retained
        ),
        primary_relaxed_cds_chain_contained_in_candidate=(
            chain_containment.retained_in_candidate
        ),
        bidirectional_relaxed_cds_chain_containment=chain_containment.bidirectional,
        mosaic_components_complete=mosaic_complete,
        mosaic_additional_corroboration=mosaic_additional,
        exact_model_source_requirement_met=exact_model_requirement_met,
        mosaic_source_requirement_met=mosaic_source_requirement_met,
        mosaic_policy_required=mosaic_required,
        factual_support_routes=factual_routes,
        qualifying_support_routes=routes,
    )


def _retained_conflict(
    candidate: CodingModel,
    retained: CodingModel,
    wobble_bp: int,
) -> str | None:
    """Return the first conflict in the authoritative structural precedence."""
    if cds_exon_chains_match(candidate, retained, wobble_bp):
        return "exact_equivalent"
    if cds_intron_chains_match(candidate, retained):
        return "same_intron_chain_terminal_variant"
    containment = relaxed_contiguous_cds_chain_containment_relation(
        candidate, retained
    )
    if containment.bidirectional:
        return "bidirectional_relaxed_cds_chain_containment"
    if containment.candidate_in_retained:
        return "candidate_relaxed_cds_chain_contained_in_retained"
    if containment.retained_in_candidate:
        return "retained_relaxed_cds_chain_contained_in_candidate"
    return None


@dataclass(frozen=True, slots=True)
class _Evaluation:
    decisions: tuple[AlternativeCandidateDecision, ...]
    accepted: tuple[AcceptedAlternative, ...]


def _evaluate(
    *,
    candidate_rows: Iterable[Mapping[str, object]],
    model_by_internal_id: Mapping[str, CodingModel],
    primary_public_gene_by_internal_id: Mapping[str, str],
    thresholds: AlternativeSpliceThresholds,
    support_policy: AlternativeSpliceSupportPolicy,
    evidence_requirement: AlternativeSpliceEvidenceRequirement,
    effective_model_source_minimum: int,
    max_alternatives_per_gene: int,
    terminal_cds_3prime_wobble_bp: int,
    assign_public_ids: bool,
) -> _Evaluation:
    rows = [dict(row) for row in candidate_rows]
    decisions: list[AlternativeCandidateDecision] = []
    pending_by_primary: dict[str, list[tuple[Mapping[str, object], AlternativeCandidateDecision]]] = defaultdict(list)
    for row in rows:
        candidate_id = str(row.get("candidate_output_transcript_id", "."))
        primary_id_text = str(row.get("primary_output_transcript_id", "."))
        primary_id = None if primary_id_text == "." else primary_id_text
        decision = _base_decision(
            row,
            model_by_internal_id.get(candidate_id),
            model_by_internal_id.get(primary_id) if primary_id is not None else None,
            thresholds,
            support_policy,
            evidence_requirement,
            effective_model_source_minimum,
        )
        if decision.stage2b_decision == "pending_greedy_selection" and primary_id is not None:
            pending_by_primary[primary_id].append((row, decision))
        else:
            decisions.append(decision)

    accepted: list[AcceptedAlternative] = []
    for primary_id in sorted(
        pending_by_primary,
        key=lambda internal_id: (
            primary_public_gene_by_internal_id.get(internal_id, "~"), internal_id
        ),
    ):
        primary = model_by_internal_id[primary_id]
        public_gene_id = primary_public_gene_by_internal_id.get(primary_id, ".")
        if assign_public_ids and public_gene_id == ".":
            raise ValueError(
                "Alternative-splice primary has no public gene ID: " + primary_id
            )
        retained: list[tuple[CodingModel, str, str]] = [
            (
                primary,
                public_transcript_id(public_gene_id, 1) if assign_public_ids else ".",
                primary_id,
            )
        ]
        accepted_for_gene = 0
        ordered = sorted(
            pending_by_primary[primary_id],
            key=lambda item: (
                item[1].alternative_candidate_order_rank
                if item[1].alternative_candidate_order_rank is not None
                else 2**31,
                item[1].candidate_internal_transcript_id,
            ),
        )
        for _row, decision in ordered:
            candidate = model_by_internal_id[decision.candidate_internal_transcript_id]
            conflict: tuple[str, str, str] | None = None
            for retained_model, public_id, internal_id in retained:
                conflict_type = _retained_conflict(
                    candidate, retained_model, terminal_cds_3prime_wobble_bp
                )
                if conflict_type is not None:
                    conflict = (conflict_type, public_id, internal_id)
                    break
            if conflict is not None:
                conflict_type, conflict_public, conflict_internal = conflict
                conflict_decision = (
                    "rejected_greedy_exact_equivalent"
                    if conflict_type == "exact_equivalent"
                    else "rejected_greedy_same_intron_chain"
                    if conflict_type == "same_intron_chain_terminal_variant"
                    else "rejected_greedy_containment"
                )
                decisions.append(replace(
                    decision,
                    stage2b_decision=conflict_decision,
                    stage2b_rejection_reason=conflict_type,
                    greedy_conflict_public_transcript_id=(
                        conflict_public if assign_public_ids else None
                    ),
                    greedy_conflict_internal_transcript_id=conflict_internal,
                    greedy_conflict_type=conflict_type,
                ))
                continue
            if accepted_for_gene >= max_alternatives_per_gene:
                decisions.append(replace(
                    decision,
                    stage2b_decision="rejected_max_alternatives_per_gene",
                    stage2b_rejection_reason="configured_max_alternatives_per_gene_reached",
                ))
                continue
            accepted_for_gene += 1
            ordinal = accepted_for_gene + 1
            public_id = (
                public_transcript_id(public_gene_id, ordinal) if assign_public_ids else None
            )
            accepted_decision = replace(
                decision,
                stage2b_decision="accepted",
                stage2b_rejection_reason=".",
                accepted_alternative_ordinal=ordinal if assign_public_ids else None,
                accepted_public_transcript_id=public_id,
            )
            decisions.append(accepted_decision)
            if assign_public_ids:
                assert public_id is not None
                assert decision.stage2b_reciprocal_cds_overlap is not None
                assert decision.stage2b_frame_concordance is not None
                assert decision.stage2b_candidate_primary_model_selection_score_ratio is not None
                accepted.append(AcceptedAlternative(
                    primary_internal_transcript_id=primary_id,
                    public_gene_id=public_gene_id,
                    candidate_model=candidate,
                    alternative_candidate_order_rank=(
                        decision.alternative_candidate_order_rank or 0
                    ),
                    accepted_alternative_ordinal=ordinal,
                    accepted_public_transcript_id=public_id,
                    factual_support_routes=decision.factual_support_routes,
                    qualifying_support_routes=decision.qualifying_support_routes,
                    evidence_requirement=evidence_requirement.value,
                    mosaic_policy_required=decision.mosaic_policy_required,
                    stage2b_reciprocal_cds_overlap=(
                        decision.stage2b_reciprocal_cds_overlap
                    ),
                    stage2b_frame_concordance=decision.stage2b_frame_concordance,
                    stage2b_primary_score_ratio=(
                        decision.stage2b_candidate_primary_model_selection_score_ratio
                    ),
                ))
                retained.append((candidate, public_id, candidate.output_transcript_id))
            else:
                retained.append((candidate, ".", candidate.output_transcript_id))

    return _Evaluation(
        decisions=tuple(sorted(decisions, key=lambda decision: decision.assessment_id)),
        accepted=tuple(accepted),
    )


def select_alternative_splice_candidates(
    *,
    candidate_rows: Iterable[Mapping[str, object]],
    all_scored_models: Iterable[CodingModel],
    primary_annotation: FinalAnnotation,
    config: AlternativeSpliceSelectionConfig,
    terminal_cds_3prime_wobble_bp: int,
) -> AlternativeSpliceSelectionResult:
    """Apply active and canonical Stage 2B greedy selection to existing facts."""
    rows = tuple(dict(row) for row in candidate_rows)
    scored_models = tuple(all_scored_models)
    config = config.for_prediction_source_count(
        len({model.source_name for model in scored_models})
    )
    model_by_id = {
        model.output_transcript_id: model for model in scored_models
    }
    public_gene_by_primary = {
        internal_id: gene.gene_id
        for internal_id, gene in primary_annotation.gene_by_primary_internal_transcript_id().items()
    }
    active = _evaluate(
        candidate_rows=rows,
        model_by_internal_id=model_by_id,
        primary_public_gene_by_internal_id=public_gene_by_primary,
        thresholds=config.effective_thresholds,
        support_policy=config.support_policy,
        evidence_requirement=config.evidence_requirement,
        effective_model_source_minimum=(
            config.effective_min_model_support_sources or 2
        ),
        max_alternatives_per_gene=config.max_alternatives_per_gene,
        terminal_cds_3prime_wobble_bp=terminal_cds_3prime_wobble_bp,
        assign_public_ids=True,
    )

    comparisons: list[AlternativePresetComparison] = []
    comparison_evaluations: dict[AlternativeSplicePreset, _Evaluation] = {}
    for preset in AlternativeSplicePreset:
        evaluation = _evaluate(
            candidate_rows=rows,
            model_by_internal_id=model_by_id,
            primary_public_gene_by_internal_id=public_gene_by_primary,
            thresholds=CANONICAL_ALTERNATIVE_SPLICE_PRESETS[preset],
            support_policy=config.support_policy,
            evidence_requirement=config.evidence_requirement,
            effective_model_source_minimum=(
                config.effective_min_model_support_sources or 2
            ),
            max_alternatives_per_gene=config.max_alternatives_per_gene,
            terminal_cds_3prime_wobble_bp=terminal_cds_3prime_wobble_bp,
            assign_public_ids=False,
        )
        comparison_evaluations[preset] = evaluation
        accepted_ids = tuple(
            decision.candidate_internal_transcript_id
            for decision in evaluation.decisions
            if decision.stage2b_decision == "accepted"
        )
        comparisons.append(AlternativePresetComparison(
            preset=preset,
            accepted_internal_transcript_ids=accepted_ids,
            decisions_by_assessment_id=tuple(
                (decision.assessment_id, decision.stage2b_decision)
                for decision in evaluation.decisions
            ),
        ))

    nonmonotonicities = _diagnose_preset_nonmonotonicities(
        comparison_evaluations
    )

    accepted_by_primary: dict[str, list[AcceptedAlternative]] = defaultdict(list)
    for accepted in active.accepted:
        accepted_by_primary[accepted.primary_internal_transcript_id].append(accepted)
    gene_decisions = tuple(
        AlternativeGeneDecision(
            primary_internal_transcript_id=primary_id,
            public_gene_id=public_gene_by_primary[primary_id],
            accepted_alternatives=tuple(accepted_by_primary.get(primary_id, [])),
        )
        for primary_id in sorted(
            public_gene_by_primary, key=lambda item: public_gene_by_primary[item]
        )
    )
    return AlternativeSpliceSelectionResult(
        config=config,
        candidate_decisions=active.decisions,
        gene_decisions=gene_decisions,
        accepted_alternatives=active.accepted,
        canonical_comparisons=tuple(comparisons),
        preset_nonmonotonicities=nonmonotonicities,
    )


def _config_fields(config: AlternativeSpliceSelectionConfig) -> dict[str, object]:
    canonical = config.canonical_thresholds
    effective = config.effective_thresholds
    overridden = set(config.overridden_thresholds)
    return {
        "alternative_splice_preset_requested": config.requested_preset.value,
        "alternative_splice_canonical_preset": config.requested_preset.value,
        "alternative_splice_support_policy": config.support_policy.value,
        "alternative_splice_evidence_requirement": config.evidence_requirement.value,
        "alternative_splice_max_alternatives_per_gene": config.max_alternatives_per_gene,
        "alternative_splice_prediction_source_count": (
            config.prediction_source_count
            if config.prediction_source_count is not None else "."
        ),
        "alternative_splice_min_model_support_sources": (
            config.min_model_support_sources
        ),
        "alternative_splice_min_model_support_source_fraction": (
            config.min_model_support_source_fraction
            if config.min_model_support_source_fraction is not None else "."
        ),
        "alternative_splice_effective_min_model_support_sources": (
            config.effective_min_model_support_sources
            if config.effective_min_model_support_sources is not None else "."
        ),
        "canonical_min_reciprocal_cds_overlap": canonical.min_reciprocal_cds_overlap,
        "canonical_min_frame_concordance": canonical.min_frame_concordance,
        "canonical_min_model_selection_score_scaled": canonical.min_model_selection_score_scaled,
        "canonical_min_primary_score_ratio": canonical.min_primary_score_ratio,
        "effective_min_reciprocal_cds_overlap": effective.min_reciprocal_cds_overlap,
        "effective_min_frame_concordance": effective.min_frame_concordance,
        "effective_min_model_selection_score_scaled": effective.min_model_selection_score_scaled,
        "effective_min_primary_score_ratio": effective.min_primary_score_ratio,
        "override_min_reciprocal_cds_overlap": "min_reciprocal_cds_overlap" in overridden,
        "override_min_frame_concordance": "min_frame_concordance" in overridden,
        "override_min_model_selection_score_scaled": "min_model_selection_score_scaled" in overridden,
        "override_min_primary_score_ratio": "min_primary_score_ratio" in overridden,
    }


def extend_alternative_splice_report(
    report: Mapping[str, list[dict[str, object]]],
    selection: AlternativeSpliceSelectionResult,
    *,
    alternatives_emitted: bool,
) -> dict[str, list[dict[str, object]]]:
    """Return copied 0.2.9 report rows extended with Stage 2B decisions."""
    config_fields = _config_fields(selection.config)
    active_by_assessment = selection.decisions_by_assessment_id()
    canonical_by_preset = {
        comparison.preset: dict(comparison.decisions_by_assessment_id)
        for comparison in selection.canonical_comparisons
    }
    candidate_rows: list[dict[str, object]] = []
    for historical in report["candidate_rows"]:
        row = dict(historical)
        assessment_id = str(row["alternative_candidate_assessment_id"])
        decision = active_by_assessment[assessment_id]
        row.update(config_fields)
        row.update({
            "alternative_splice_candidate_cds_bases_contained_in_primary": decision.candidate_cds_bases_contained_in_primary,
            "alternative_splice_primary_cds_bases_contained_in_candidate": decision.primary_cds_bases_contained_in_candidate,
            "alternative_splice_bidirectional_cds_base_containment": decision.bidirectional_cds_base_containment,
            "alternative_splice_candidate_relaxed_cds_chain_contained_in_primary": decision.candidate_relaxed_cds_chain_contained_in_primary,
            "alternative_splice_primary_relaxed_cds_chain_contained_in_candidate": decision.primary_relaxed_cds_chain_contained_in_candidate,
            "alternative_splice_bidirectional_relaxed_cds_chain_containment": decision.bidirectional_relaxed_cds_chain_containment,
            "alternative_splice_cds_overlap_bp": decision.stage2b_cds_overlap_bp if decision.stage2b_cds_overlap_bp is not None else ".",
            "alternative_splice_primary_cds_overlap_fraction": decision.stage2b_primary_cds_overlap_fraction if decision.stage2b_primary_cds_overlap_fraction is not None else ".",
            "alternative_splice_candidate_cds_overlap_fraction": decision.stage2b_candidate_cds_overlap_fraction if decision.stage2b_candidate_cds_overlap_fraction is not None else ".",
            "alternative_splice_reciprocal_cds_overlap": decision.stage2b_reciprocal_cds_overlap if decision.stage2b_reciprocal_cds_overlap is not None else ".",
            "alternative_splice_primary_score_ratio_status": decision.stage2b_primary_score_ratio_status,
            "alternative_splice_candidate_primary_model_selection_score_ratio": (
                decision.stage2b_candidate_primary_model_selection_score_ratio
                if decision.stage2b_candidate_primary_model_selection_score_ratio is not None
                else "."
            ),
            "mosaic_components_complete": decision.mosaic_components_complete,
            "mosaic_additional_corroboration": decision.mosaic_additional_corroboration,
            "alternative_splice_exact_model_source_requirement_met": (
                decision.exact_model_source_requirement_met
            ),
            "alternative_splice_mosaic_source_requirement_met": (
                decision.mosaic_source_requirement_met
            ),
            "mosaic_policy_required": decision.mosaic_policy_required,
            "alternative_splice_qualifying_support_routes": (
                ",".join(decision.qualifying_support_routes)
                if decision.qualifying_support_routes else "."
            ),
            "alternative_splice_fixed_rule_eligible": decision.fixed_eligible,
            "alternative_splice_support_policy_eligible": decision.support_policy_eligible,
            "alternative_splice_external_evidence_requirement_met": decision.external_evidence_requirement_met,
            "alternative_splice_external_evidence_requirement_routes": (
                ",".join(decision.external_evidence_requirement_routes)
                if decision.external_evidence_requirement_routes else "."
            ),
            "alternative_splice_prethreshold_eligible": decision.prethreshold_eligible,
            "alternative_splice_numerical_threshold_survivor": decision.preset_threshold_survivor,
            "alternative_splice_decision": decision.stage2b_decision,
            "alternative_splice_rejection_reason": decision.stage2b_rejection_reason,
            "greedy_conflict_public_transcript_id": decision.greedy_conflict_public_transcript_id or ".",
            "greedy_conflict_internal_transcript_id": decision.greedy_conflict_internal_transcript_id or ".",
            "greedy_conflict_type": decision.greedy_conflict_type or ".",
            "accepted_alternative_ordinal": decision.accepted_alternative_ordinal or ".",
            "accepted_public_transcript_id": decision.accepted_public_transcript_id or ".",
            "canonical_permissive_alternative_splice_decision": canonical_by_preset[AlternativeSplicePreset.PERMISSIVE][assessment_id],
            "canonical_intermediate_alternative_splice_decision": canonical_by_preset[AlternativeSplicePreset.INTERMEDIATE][assessment_id],
            "canonical_stringent_alternative_splice_decision": canonical_by_preset[AlternativeSplicePreset.STRINGENT][assessment_id],
        })
        candidate_rows.append(row)

    comparison_by_preset = {
        comparison.preset: comparison for comparison in selection.canonical_comparisons
    }
    accepted_by_primary = selection.accepted_by_primary_internal_id()
    decisions_by_primary: dict[str, list[AlternativeCandidateDecision]] = defaultdict(list)
    for decision in selection.candidate_decisions:
        if decision.primary_internal_transcript_id is not None:
            decisions_by_primary[decision.primary_internal_transcript_id].append(decision)
    def comparison_ids_for_primary(
        preset: AlternativeSplicePreset, primary_id: str
    ) -> tuple[str, ...]:
        accepted_set = set(comparison_by_preset[preset].accepted_internal_transcript_ids)
        return tuple(sorted(
            decision.candidate_internal_transcript_id
            for decision in decisions_by_primary.get(primary_id, [])
            if decision.candidate_internal_transcript_id in accepted_set
        ))

    gene_rows: list[dict[str, object]] = []
    for historical in report["gene_rows"]:
        row = dict(historical)
        primary_id = str(row["primary_output_transcript_id"])
        decisions = decisions_by_primary.get(primary_id, [])
        accepted = accepted_by_primary.get(primary_id, ())
        ranks = [
            item.alternative_candidate_order_rank for item in accepted
            if item.alternative_candidate_order_rank is not None
        ]
        per_preset = {
            preset: comparison_ids_for_primary(preset, primary_id)
            for preset in AlternativeSplicePreset
        }
        row.update({
            "alternative_splice_fixed_rule_eligible_candidates": sum(decision.fixed_eligible for decision in decisions),
            "alternative_splice_support_policy_eligible_candidates": sum(decision.support_policy_eligible for decision in decisions),
            "alternative_splice_external_evidence_requirement_failure_candidates": sum(
                decision.stage2b_decision == "rejected_external_evidence_requirement"
                for decision in decisions
            ),
            "alternative_splice_prethreshold_eligible_candidates": sum(decision.prethreshold_eligible for decision in decisions),
            "alternative_splice_numerical_threshold_survivors": sum(decision.preset_threshold_survivor for decision in decisions),
            "alternative_splice_greedy_accepted_alternatives": len(accepted),
            "alternative_splice_greedy_redundancy_rejections": sum(
                decision.stage2b_decision.startswith("rejected_greedy_") for decision in decisions
            ),
            "alternative_splice_cap_rejections": sum(
                decision.stage2b_decision == "rejected_max_alternatives_per_gene" for decision in decisions
            ),
            "alternative_splice_max_alternatives_per_gene": selection.config.max_alternatives_per_gene,
            "accepted_alternative_count": len(accepted),
            "max_candidate_rank_accepted": max(ranks) if ranks else ".",
            "accepted_internal_transcript_ids": _csv(
                item.candidate_model.output_transcript_id for item in accepted
            ),
            "canonical_permissive_accepted_count": len(per_preset[AlternativeSplicePreset.PERMISSIVE]),
            "canonical_permissive_accepted_internal_transcript_ids": _csv(per_preset[AlternativeSplicePreset.PERMISSIVE]),
            "canonical_intermediate_accepted_count": len(per_preset[AlternativeSplicePreset.INTERMEDIATE]),
            "canonical_intermediate_accepted_internal_transcript_ids": _csv(per_preset[AlternativeSplicePreset.INTERMEDIATE]),
            "canonical_stringent_accepted_count": len(per_preset[AlternativeSplicePreset.STRINGENT]),
            "canonical_stringent_accepted_internal_transcript_ids": _csv(per_preset[AlternativeSplicePreset.STRINGENT]),
        })
        gene_rows.append(row)

    unique_decisions: dict[str, AlternativeCandidateDecision] = {}
    decisions_by_candidate: dict[str, list[AlternativeCandidateDecision]] = defaultdict(list)
    for decision in selection.candidate_decisions:
        unique_decisions.setdefault(decision.candidate_internal_transcript_id, decision)
        decisions_by_candidate[decision.candidate_internal_transcript_id].append(decision)
    decision_values = tuple(unique_decisions.values())
    assessment_values = selection.candidate_decisions

    def candidate_count_with_ratio_status(status: str) -> int:
        return sum(
            any(item.stage2b_primary_score_ratio_status == status for item in items)
            for items in decisions_by_candidate.values()
        )
    accepted = selection.accepted_alternatives
    route_counts = Counter(
        route for item in accepted for route in item.qualifying_support_routes
    )
    route_combinations = Counter(
        "+".join(item.qualifying_support_routes) for item in accepted
    )
    conflict_type_counts = Counter(
        item.greedy_conflict_type
        for item in decision_values
        if item.greedy_conflict_type is not None
    )
    gene_accepted_counts = Counter(
        item.primary_internal_transcript_id for item in accepted
    )
    comparisons = comparison_by_preset
    nonmonotonicity_details = ";".join(
        f"{item.looser_preset.value}->{item.stricter_preset.value}:"
        f"{item.candidate_internal_transcript_id}:"
        f"{item.looser_decision}:"
        f"{item.looser_conflict_internal_transcript_id or '.'}:"
        f"{item.looser_conflict_type or '.'}"
        for item in selection.preset_nonmonotonicities
    ) or "."
    run_rows: list[dict[str, object]] = []
    for historical in report["run_rows"]:
        row = dict(historical)
        row.update(config_fields)
        row.update({
            "alternative_splice_selector_evaluated": True,
            "alternative_splice_alternatives_emitted": alternatives_emitted,
            "alternative_splice_fixed_rule_eligible_candidate_count": sum(item.fixed_eligible for item in decision_values),
            "alternative_splice_support_policy_eligible_candidate_count": sum(item.support_policy_eligible for item in decision_values),
            "alternative_splice_external_evidence_requirement_failure_candidate_count": sum(
                item.stage2b_decision == "rejected_external_evidence_requirement"
                for item in decision_values
            ),
            "alternative_splice_prethreshold_eligible_candidate_count": sum(item.prethreshold_eligible for item in decision_values),
            "alternative_splice_numerical_threshold_survivor_count": sum(item.preset_threshold_survivor for item in decision_values),
            "alternative_splice_genes_with_accepted_alternatives": len(gene_accepted_counts),
            "alternative_splice_total_accepted_alternatives": len(accepted),
            "alternative_splice_max_accepted_alternatives_per_gene": max(gene_accepted_counts.values(), default=0),
            "alternative_splice_exact_model_support_accepted_count": route_counts["exact_multisource_model_support"],
            "alternative_splice_exact_protein_support_accepted_count": route_counts["exact_protein_support"],
            "alternative_splice_exact_transcript_support_accepted_count": route_counts["exact_coherent_transcript_support"],
            "alternative_splice_mosaic_only_accepted_count": sum(item.mosaic_policy_required for item in accepted),
            "alternative_splice_accepted_support_route_composition": (
                ",".join(
                    f"{route}:{count}"
                    for route, count in sorted(route_combinations.items())
                ) or "."
            ),
            "alternative_splice_cap_rejected_count": sum(item.stage2b_decision == "rejected_max_alternatives_per_gene" for item in decision_values),
            "alternative_splice_greedy_exact_conflict_count": sum(item.greedy_conflict_type == "exact_equivalent" for item in decision_values),
            "alternative_splice_greedy_same_intron_chain_conflict_count": sum(item.greedy_conflict_type == "same_intron_chain_terminal_variant" for item in decision_values),
            "alternative_splice_greedy_containment_conflict_count": sum(
                item.greedy_conflict_type is not None
                and "contain" in item.greedy_conflict_type
                for item in decision_values
            ),
            "alternative_splice_greedy_conflict_type_counts": (
                ",".join(
                    f"{conflict_type}:{count}"
                    for conflict_type, count in sorted(conflict_type_counts.items())
                ) or "."
            ),
            "alternative_splice_primary_score_ratio_unavailable_zero_count": candidate_count_with_ratio_status("unavailable_zero_primary"),
            "alternative_splice_primary_score_ratio_unavailable_negative_count": candidate_count_with_ratio_status("unavailable_negative_primary"),
            "alternative_splice_primary_score_ratio_unavailable_zero_assessment_count": sum(item.stage2b_primary_score_ratio_status == "unavailable_zero_primary" for item in assessment_values),
            "alternative_splice_primary_score_ratio_unavailable_negative_assessment_count": sum(item.stage2b_primary_score_ratio_status == "unavailable_negative_primary" for item in assessment_values),
            "canonical_permissive_accepted_count": len(comparisons[AlternativeSplicePreset.PERMISSIVE].accepted_internal_transcript_ids),
            "canonical_permissive_accepted_internal_transcript_ids": _csv(comparisons[AlternativeSplicePreset.PERMISSIVE].accepted_internal_transcript_ids),
            "canonical_intermediate_accepted_count": len(comparisons[AlternativeSplicePreset.INTERMEDIATE].accepted_internal_transcript_ids),
            "canonical_intermediate_accepted_internal_transcript_ids": _csv(comparisons[AlternativeSplicePreset.INTERMEDIATE].accepted_internal_transcript_ids),
            "canonical_stringent_accepted_count": len(comparisons[AlternativeSplicePreset.STRINGENT].accepted_internal_transcript_ids),
            "canonical_stringent_accepted_internal_transcript_ids": _csv(comparisons[AlternativeSplicePreset.STRINGENT].accepted_internal_transcript_ids),
            "alternative_splice_preset_nonmonotonicity_count": len(selection.preset_nonmonotonicities),
            "alternative_splice_preset_nonmonotonicity_details": nonmonotonicity_details,
        })
        run_rows.append(row)
    return {
        "candidate_rows": candidate_rows,
        "gene_rows": gene_rows,
        "run_rows": run_rows,
    }
