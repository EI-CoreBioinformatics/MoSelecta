"""Reference-free automatic final-pruning assessment.

The expensive model, protein, and transcript calculations happen before this
module is used.  It only summarizes cached selected-model/pruning results and
chooses one of the fixed medium, high, or very-high scenarios.
"""

from __future__ import annotations

import csv
import math
import statistics
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .models import CodingModel, Sublocus
from .scoring import value_as_float


AUTO_SCENARIO_LEVELS = ("medium", "high", "very-high")
LEVEL_ORDER = {"low": 0, "medium": 1, "high": 2, "very-high": 3}

# Source-consensus classification thresholds.  The only source-count-specific
# strong criterion is deep-consensus incidence: with three inputs, deep
# consensus is unanimity, whereas with four or more inputs it is support from
# at least three sources.
CONSENSUS_MIN_SOURCE_COUNT = 3
CONSENSUS_MIN_MULTISOURCE_SUBLOCI = 100
STRONG_MAX_SINGLETON_INCIDENCE_FRACTION = 0.35
STRONG_MIN_CORROBORATION_FRACTION = 0.45
STRONG_MIN_DEEP_CONSENSUS_INCIDENCE_THREE_SOURCES = 0.55
STRONG_MIN_DEEP_CONSENSUS_INCIDENCE_FOUR_PLUS_SOURCES = 0.60
STRONG_MIN_DEEP_CONSENSUS_SUBLOCUS_FRACTION = 0.25
INTERMEDIATE_MAX_SINGLETON_INCIDENCE_FRACTION = 0.55
INTERMEDIATE_MIN_CORROBORATION_FRACTION = 0.30
INTERMEDIATE_MIN_DEEP_CONSENSUS_INCIDENCE_FRACTION = 0.35
INTERMEDIATE_MIN_DEEP_CONSENSUS_SUBLOCUS_FRACTION = 0.15

STRONG_SOURCE_EXCESS_MIN_SINGLETON_SHARE = 0.60
STRONG_SOURCE_EXCESS_MIN_SINGLETON_RATE = 0.20
STRONG_SOURCE_EXCESS_MIN_RATE_RATIO = 2.0

INPUT_MODEL_COUNT_HIGH_FACTOR = 2.0
INPUT_MODEL_COUNT_LOW_FACTOR = 0.5


def _safe_fraction(numerator: int | float, denominator: int | float) -> float:
    return float(numerator) / float(denominator) if denominator else 0.0


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _model_id(model: CodingModel) -> str:
    return model.output_transcript_id


def _serialize_counter(counter: Mapping[str, int]) -> str:
    if not counter:
        return "."
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter))


def _serialize_float_map(values: Mapping[str, float]) -> str:
    if not values:
        return "."
    return ",".join(f"{key}:{values[key]:.6f}" for key in sorted(values))


@dataclass(frozen=True)
class EvidenceAdequacy:
    evidence_class: str
    state: str
    supplied: bool
    first_value: float | None
    second_value: float | None
    first_threshold: float
    second_threshold: float
    first_margin: float | None
    second_margin: float | None
    limiting_margin: float | None

    @property
    def passed(self) -> bool:
        return self.state == "pass"


def classify_evidence_adequacy(
    evidence_class: str,
    *,
    supplied: bool,
    first_value: float,
    second_value: float,
    first_threshold: float,
    second_threshold: float,
    marginal_band: float = 0.05,
) -> EvidenceAdequacy:
    """Classify a two-metric adequacy gate as absent/pass/marginal/fail."""
    if not supplied:
        return EvidenceAdequacy(
            evidence_class, "absent", False, None, None,
            first_threshold, second_threshold, None, None, None,
        )
    first_margin = first_value - first_threshold
    second_margin = second_value - second_threshold
    limiting = min(first_margin, second_margin)
    state = "pass" if limiting >= -1e-15 else ("marginal" if limiting >= -marginal_band - 1e-15 else "fail")
    return EvidenceAdequacy(
        evidence_class, state, True, first_value, second_value,
        first_threshold, second_threshold, first_margin, second_margin, limiting,
    )


@dataclass
class SourceConsensusMetrics:
    source_count_distribution: dict[int, int]
    total_subloci: int
    total_source_incidences: int
    raw_source_count: int
    effective_source_count: int
    singleton_sublocus_fraction: float
    multisource_sublocus_fraction: float
    deep_consensus_sublocus_fraction: float
    full_consensus_sublocus_fraction: float
    singleton_incidence_fraction: float
    consensus_incidence_fraction: float
    deep_consensus_incidence_fraction: float
    corroboration_fraction: float
    mean_sources_per_sublocus: float
    mean_sources_at_multisource_subloci: float
    multisource_subloci: int
    per_source_singleton_count: dict[str, int]
    per_source_total_subloci: dict[str, int]
    per_source_singleton_share: dict[str, float]
    per_source_singleton_rate: dict[str, float]
    per_source_other_median_rate: dict[str, float]
    per_source_singleton_rate_difference: dict[str, float]
    per_source_singleton_rate_ratio: dict[str, float]
    strong_source_excess_sources: tuple[str, ...]
    backbone_class: str


def strong_deep_consensus_incidence_threshold(source_count: int) -> float:
    """Return the strong threshold appropriate to the original input count."""
    if source_count == 3:
        return STRONG_MIN_DEEP_CONSENSUS_INCIDENCE_THREE_SOURCES
    return STRONG_MIN_DEEP_CONSENSUS_INCIDENCE_FOUR_PLUS_SOURCES


def consensus_classification_criteria(
    metrics: SourceConsensusMetrics,
) -> dict[str, dict[str, object]]:
    """Return auditable strong/intermediate criterion values and margins.

    Margins are positive when a criterion passes, zero at its threshold, and
    negative when it fails.  This helper is diagnostic only; classification
    below consumes the same named constants directly.
    """
    deep_threshold = strong_deep_consensus_incidence_threshold(
        metrics.effective_source_count
    )

    def minimum(value: float | int, threshold: float | int) -> dict[str, object]:
        margin = float(value) - float(threshold)
        return {
            "value": value,
            "operator": ">=",
            "threshold": threshold,
            "margin": margin,
            "passed": margin >= -1e-15,
        }

    def maximum(value: float | int, threshold: float | int) -> dict[str, object]:
        margin = float(threshold) - float(value)
        return {
            "value": value,
            "operator": "<=",
            "threshold": threshold,
            "margin": margin,
            "passed": margin >= -1e-15,
        }

    return {
        "strong.minimum_source_count": minimum(
            metrics.effective_source_count, CONSENSUS_MIN_SOURCE_COUNT
        ),
        "strong.minimum_multisource_subloci": minimum(
            metrics.multisource_subloci, CONSENSUS_MIN_MULTISOURCE_SUBLOCI
        ),
        "strong.maximum_singleton_incidence_fraction": maximum(
            metrics.singleton_incidence_fraction,
            STRONG_MAX_SINGLETON_INCIDENCE_FRACTION,
        ),
        "strong.minimum_corroboration_fraction": minimum(
            metrics.corroboration_fraction, STRONG_MIN_CORROBORATION_FRACTION
        ),
        "strong.minimum_deep_consensus_incidence_fraction": minimum(
            metrics.deep_consensus_incidence_fraction, deep_threshold
        ),
        "strong.minimum_deep_consensus_sublocus_fraction": minimum(
            metrics.deep_consensus_sublocus_fraction,
            STRONG_MIN_DEEP_CONSENSUS_SUBLOCUS_FRACTION,
        ),
        "intermediate.minimum_source_count": minimum(
            metrics.effective_source_count, CONSENSUS_MIN_SOURCE_COUNT
        ),
        "intermediate.minimum_multisource_subloci": minimum(
            metrics.multisource_subloci, CONSENSUS_MIN_MULTISOURCE_SUBLOCI
        ),
        "intermediate.maximum_singleton_incidence_fraction": maximum(
            metrics.singleton_incidence_fraction,
            INTERMEDIATE_MAX_SINGLETON_INCIDENCE_FRACTION,
        ),
        "intermediate.minimum_corroboration_fraction": minimum(
            metrics.corroboration_fraction,
            INTERMEDIATE_MIN_CORROBORATION_FRACTION,
        ),
        "intermediate.minimum_deep_consensus_incidence_fraction": minimum(
            metrics.deep_consensus_incidence_fraction,
            INTERMEDIATE_MIN_DEEP_CONSENSUS_INCIDENCE_FRACTION,
        ),
        "intermediate.minimum_deep_consensus_sublocus_fraction": minimum(
            metrics.deep_consensus_sublocus_fraction,
            INTERMEDIATE_MIN_DEEP_CONSENSUS_SUBLOCUS_FRACTION,
        ),
    }


def classify_consensus_backbone(metrics: SourceConsensusMetrics) -> str:
    """Apply the v0.2.14 source-consensus thresholds."""
    common = (
        metrics.effective_source_count >= CONSENSUS_MIN_SOURCE_COUNT
        and metrics.multisource_subloci >= CONSENSUS_MIN_MULTISOURCE_SUBLOCI
    )
    strong = (
        common
        and metrics.singleton_incidence_fraction
        <= STRONG_MAX_SINGLETON_INCIDENCE_FRACTION + 1e-15
        and metrics.corroboration_fraction
        >= STRONG_MIN_CORROBORATION_FRACTION - 1e-15
        and metrics.deep_consensus_incidence_fraction
        >= strong_deep_consensus_incidence_threshold(
            metrics.effective_source_count
        ) - 1e-15
        and metrics.deep_consensus_sublocus_fraction
        >= STRONG_MIN_DEEP_CONSENSUS_SUBLOCUS_FRACTION - 1e-15
    )
    if strong:
        return "strong"
    intermediate = (
        common
        and metrics.singleton_incidence_fraction
        <= INTERMEDIATE_MAX_SINGLETON_INCIDENCE_FRACTION + 1e-15
        and metrics.corroboration_fraction
        >= INTERMEDIATE_MIN_CORROBORATION_FRACTION - 1e-15
        and metrics.deep_consensus_incidence_fraction
        >= INTERMEDIATE_MIN_DEEP_CONSENSUS_INCIDENCE_FRACTION - 1e-15
        and metrics.deep_consensus_sublocus_fraction
        >= INTERMEDIATE_MIN_DEEP_CONSENSUS_SUBLOCUS_FRACTION - 1e-15
    )
    return "intermediate" if intermediate else "weak_or_discordant"


def calculate_source_consensus_metrics(
    subloci: Sequence[Sublocus],
    *,
    selected_model_ids: set[str] | None = None,
) -> SourceConsensusMetrics:
    """Summarize source agreement for subloci represented in final selection."""
    selected_model_ids = selected_model_ids or {
        _model_id(model) for sub in subloci for model in sub.models
    }
    distribution: Counter[int] = Counter()
    source_total: Counter[str] = Counter()
    source_singleton: Counter[str] = Counter()
    all_sources: set[str] = set()
    multisource_incidence_total = 0
    multisource_subloci = 0
    kept_source_sets: list[set[str]] = []
    for sub in subloci:
        if not any(_model_id(model) in selected_model_ids for model in sub.models):
            continue
        sources = {model.source_name for model in sub.models}
        if not sources:
            continue
        kept_source_sets.append(sources)
        count = len(sources)
        distribution[count] += 1
        all_sources.update(sources)
        for source in sources:
            source_total[source] += 1
        if count == 1:
            source_singleton[next(iter(sources))] += 1
        else:
            multisource_subloci += 1
            multisource_incidence_total += count

    total_subloci = sum(distribution.values())
    incidences = sum(count * number for count, number in distribution.items())
    raw_source_count = len(all_sources)
    singleton_loci = distribution.get(1, 0)
    deep_loci = sum(number for count, number in distribution.items() if count >= 3)
    deep_incidence = sum(count * number for count, number in distribution.items() if count >= 3)
    full_loci = distribution.get(raw_source_count, 0) if raw_source_count else 0
    singleton_share: dict[str, float] = {}
    singleton_rate: dict[str, float] = {}
    other_median: dict[str, float] = {}
    rate_difference: dict[str, float] = {}
    rate_ratio: dict[str, float] = {}
    for source in sorted(all_sources):
        singleton_share[source] = _safe_fraction(source_singleton[source], singleton_loci)
        singleton_rate[source] = _safe_fraction(source_singleton[source], source_total[source])
    for source in sorted(all_sources):
        other_rates = [singleton_rate[other] for other in sorted(all_sources) if other != source]
        median = statistics.median(other_rates) if other_rates else 0.0
        other_median[source] = median
        rate_difference[source] = singleton_rate[source] - median
        rate_ratio[source] = (
            singleton_rate[source] / median
            if median > 0 else (math.inf if singleton_rate[source] > 0 else 1.0)
        )
    excess = tuple(
        source for source in sorted(all_sources)
        if singleton_share[source]
        >= STRONG_SOURCE_EXCESS_MIN_SINGLETON_SHARE - 1e-15
        and singleton_rate[source]
        >= STRONG_SOURCE_EXCESS_MIN_SINGLETON_RATE - 1e-15
        and rate_ratio[source]
        >= STRONG_SOURCE_EXCESS_MIN_RATE_RATIO - 1e-15
    )
    mean_multisource = _safe_fraction(multisource_incidence_total, multisource_subloci)
    placeholder = SourceConsensusMetrics(
        dict(sorted(distribution.items())),
        total_subloci,
        incidences,
        raw_source_count,
        raw_source_count,
        _safe_fraction(singleton_loci, total_subloci),
        _safe_fraction(multisource_subloci, total_subloci),
        _safe_fraction(deep_loci, total_subloci),
        _safe_fraction(full_loci, total_subloci),
        _safe_fraction(singleton_loci, incidences),
        _safe_fraction(incidences - singleton_loci, incidences),
        _safe_fraction(deep_incidence, incidences),
        _safe_fraction(incidences - total_subloci, incidences),
        _safe_fraction(incidences, total_subloci),
        mean_multisource,
        multisource_subloci,
        dict(sorted(source_singleton.items())),
        dict(sorted(source_total.items())),
        singleton_share,
        singleton_rate,
        other_median,
        rate_difference,
        rate_ratio,
        excess,
        ".",
    )
    placeholder.backbone_class = classify_consensus_backbone(placeholder)
    return placeholder


def consensus_metrics_from_distribution(
    source_count_distribution: Mapping[int, int],
    *,
    raw_source_count: int | None = None,
    per_source_singleton_count: Mapping[str, int] | None = None,
    per_source_total_subloci: Mapping[str, int] | None = None,
) -> SourceConsensusMetrics:
    """Build consensus metrics directly from an existing ``L_k`` summary.

    This is also the compact fixture interface used for the multispecies
    calibration distributions.  Per-source maps are optional because backbone
    classification itself depends only on the source-count distribution.
    """
    distribution = Counter({
        int(count): int(number)
        for count, number in source_count_distribution.items()
        if int(count) > 0 and int(number) > 0
    })
    raw_source_count = raw_source_count or (max(distribution) if distribution else 0)
    total_subloci = sum(distribution.values())
    incidences = sum(count * number for count, number in distribution.items())
    singleton_loci = distribution.get(1, 0)
    multisource_subloci = sum(number for count, number in distribution.items() if count >= 2)
    multisource_incidence = sum(count * number for count, number in distribution.items() if count >= 2)
    deep_loci = sum(number for count, number in distribution.items() if count >= 3)
    deep_incidence = sum(count * number for count, number in distribution.items() if count >= 3)
    full_loci = distribution.get(raw_source_count, 0) if raw_source_count else 0
    source_singleton = Counter(per_source_singleton_count or {})
    source_total = Counter(per_source_total_subloci or {})
    sources = sorted(set(source_singleton) | set(source_total))
    singleton_share = {
        source: _safe_fraction(source_singleton[source], singleton_loci)
        for source in sources
    }
    singleton_rate = {
        source: _safe_fraction(source_singleton[source], source_total[source])
        for source in sources
    }
    other_median: dict[str, float] = {}
    difference: dict[str, float] = {}
    ratio: dict[str, float] = {}
    for source in sources:
        other_rates = [singleton_rate[other] for other in sources if other != source]
        median = statistics.median(other_rates) if other_rates else 0.0
        other_median[source] = median
        difference[source] = singleton_rate[source] - median
        ratio[source] = (
            singleton_rate[source] / median
            if median > 0 else (math.inf if singleton_rate[source] > 0 else 1.0)
        )
    excess = tuple(
        source for source in sources
        if singleton_share[source]
        >= STRONG_SOURCE_EXCESS_MIN_SINGLETON_SHARE - 1e-15
        and singleton_rate[source]
        >= STRONG_SOURCE_EXCESS_MIN_SINGLETON_RATE - 1e-15
        and ratio[source]
        >= STRONG_SOURCE_EXCESS_MIN_RATE_RATIO - 1e-15
    )
    metrics = SourceConsensusMetrics(
        dict(sorted(distribution.items())),
        total_subloci,
        incidences,
        raw_source_count,
        raw_source_count,
        _safe_fraction(singleton_loci, total_subloci),
        _safe_fraction(multisource_subloci, total_subloci),
        _safe_fraction(deep_loci, total_subloci),
        _safe_fraction(full_loci, total_subloci),
        _safe_fraction(singleton_loci, incidences),
        _safe_fraction(incidences - singleton_loci, incidences),
        _safe_fraction(deep_incidence, incidences),
        _safe_fraction(incidences - total_subloci, incidences),
        _safe_fraction(incidences, total_subloci),
        _safe_fraction(multisource_incidence, multisource_subloci),
        multisource_subloci,
        dict(sorted(source_singleton.items())),
        dict(sorted(source_total.items())),
        singleton_share,
        singleton_rate,
        other_median,
        difference,
        ratio,
        excess,
        ".",
    )
    metrics.backbone_class = classify_consensus_backbone(metrics)
    return metrics


@dataclass(frozen=True)
class FloodAdjustedConsensusAssessment:
    """Raw and optional flood-adjusted views of one consensus distribution."""

    unadjusted: SourceConsensusMetrics
    adjusted: SourceConsensusMetrics | None
    adjustment_enabled: bool
    qualifying_excess_source: str | None
    discounted_source_singleton_subloci: int

    @property
    def adjustment_applied(self) -> bool:
        return self.adjusted is not None


def assess_flood_adjusted_consensus(
    consensus: SourceConsensusMetrics,
    *,
    enabled: bool = True,
) -> FloodAdjustedConsensusAssessment:
    """Discount one qualifying source's isolated subloci for classification.

    Eligibility is determined once from ``consensus``.  The adjusted metrics
    never feed back into the source-excess detector, and the original source
    count is retained even if all isolated subloci from that source are
    discounted.  No biological sublocus or model collection is modified.
    """
    qualifying = (
        consensus.strong_source_excess_sources[0]
        if len(consensus.strong_source_excess_sources) == 1
        else None
    )
    if qualifying is None or not enabled:
        return FloodAdjustedConsensusAssessment(
            consensus, None, enabled, qualifying, 0
        )

    singleton_total = consensus.source_count_distribution.get(1, 0)
    discount = min(
        max(consensus.per_source_singleton_count.get(qualifying, 0), 0),
        singleton_total,
    )
    if discount <= 0:
        return FloodAdjustedConsensusAssessment(
            consensus, None, enabled, qualifying, 0
        )

    distribution = dict(consensus.source_count_distribution)
    remaining_singletons = singleton_total - discount
    if remaining_singletons:
        distribution[1] = remaining_singletons
    else:
        distribution.pop(1, None)

    adjusted_singletons = dict(consensus.per_source_singleton_count)
    adjusted_singletons[qualifying] = max(
        adjusted_singletons.get(qualifying, 0) - discount, 0
    )
    adjusted_totals = dict(consensus.per_source_total_subloci)
    adjusted_totals[qualifying] = max(
        adjusted_totals.get(qualifying, 0) - discount, 0
    )
    adjusted = consensus_metrics_from_distribution(
        distribution,
        raw_source_count=consensus.raw_source_count,
        per_source_singleton_count=adjusted_singletons,
        per_source_total_subloci=adjusted_totals,
    )
    return FloodAdjustedConsensusAssessment(
        consensus, adjusted, enabled, qualifying, discount
    )


@dataclass(frozen=True)
class SourceInputCountMetric:
    source: str
    input_coding_model_count: int
    run_median_input_count: float
    ratio_to_run_median: float | None
    other_source_median_count: float
    ratio_to_other_source_median: float
    warning: str


def calculate_source_input_count_metrics(
    counts: Mapping[str, int],
) -> dict[str, SourceInputCountMetric]:
    """Calculate reporting-only model-count comparisons for current inputs."""
    normalized = {str(source): max(int(count), 0) for source, count in counts.items()}
    run_median = (
        float(statistics.median(normalized.values())) if normalized else 0.0
    )
    result: dict[str, SourceInputCountMetric] = {}
    for source in sorted(normalized):
        others = [normalized[other] for other in sorted(normalized) if other != source]
        other_median = float(statistics.median(others)) if others else 0.0
        count = normalized[source]
        legacy_ratio = (
            count / other_median
            if other_median > 0
            else (math.inf if count > 0 and others else 1.0)
        )
        ratio_to_run_median = count / run_median if run_median > 0 else None
        warning = "."
        if (
            ratio_to_run_median is not None
            and ratio_to_run_median >= INPUT_MODEL_COUNT_HIGH_FACTOR - 1e-15
        ):
            warning = "INPUT_CODING_MODEL_COUNT_HIGH"
        elif (
            ratio_to_run_median is not None
            and ratio_to_run_median <= INPUT_MODEL_COUNT_LOW_FACTOR + 1e-15
        ):
            warning = "INPUT_CODING_MODEL_COUNT_LOW"
        result[source] = SourceInputCountMetric(
            source=source,
            input_coding_model_count=count,
            run_median_input_count=run_median,
            ratio_to_run_median=ratio_to_run_median,
            other_source_median_count=other_median,
            ratio_to_other_source_median=legacy_ratio,
            warning=warning,
        )
    return result


PRIMARY_REASON_PRECEDENCE = (
    ("opposite_strand_shadow", "opposite_strand_shadow"),
    ("sequence_integrity_nonsingleton", "sequence_integrity_nonsingleton"),
    ("sequence_integrity_singleton", "sequence_integrity_singleton"),
    ("low_score_nonsingleton", "low_score_nonsingleton"),
    ("low_score_singleton", "low_score_singleton"),
    ("singleton_sublocus", "singleton_sublocus"),
)


def primary_reason_for_row(row: Mapping[str, object]) -> str:
    if _truthy(row.get("final_monotonic_inherited_removal")):
        return "inherited_lower_level_removal"
    classes = str(row.get("final_prune_classes", ".")).split(",")
    classes = [item for item in classes if item and item != "."]
    for class_name, primary in PRIMARY_REASON_PRECEDENCE:
        if class_name in classes:
            return primary
    return classes[0] if classes else "."


def protection_flags(row: Mapping[str, object]) -> tuple[bool, bool, bool]:
    protein = (
        _truthy(row.get("final_protein_singleton_rescued"))
        or _truthy(row.get("final_sequence_integrity_nonsingleton_protein_protected"))
        or _truthy(row.get("final_singleton_overlap_rescue_strong_local_protein"))
    )
    transcript = (
        _truthy(row.get("final_transcript_singleton_rescued"))
        or _truthy(row.get("final_sequence_integrity_nonsingleton_transcript_protected"))
        or _truthy(row.get("final_singleton_overlap_rescue_strong_local_transcript"))
    )
    overlap = str(row.get("final_singleton_overlap_rescue_decision", ".")) in {
        "rescued_strong_local_evidence", "rescued_no_adequate_evidence", "rescued_unconditional",
    }
    return protein, transcript, overlap


@dataclass
class ScenarioSummary:
    level: str
    starting_models: int
    candidate_models: int
    protected_models: int
    removed_models: int
    retained_models: int
    removed_fraction: float
    affected_subloci: int
    retained_subloci: int
    removed_ids: set[str] = field(repr=False)
    retained_ids: set[str] = field(repr=False)
    primary_reason_by_id: dict[str, str] = field(repr=False)
    all_triggered_by_id: dict[str, tuple[str, ...]] = field(repr=False)
    native_removed_ids: set[str] = field(default_factory=set, repr=False)
    inherited_removed_ids: set[str] = field(default_factory=set, repr=False)
    source_by_id: dict[str, str] = field(default_factory=dict, repr=False)
    inherited_from_level_by_id: dict[str, str] = field(default_factory=dict, repr=False)
    inherited_original_level_by_id: dict[str, str] = field(default_factory=dict, repr=False)
    inherited_original_reason_by_id: dict[str, str] = field(default_factory=dict, repr=False)
    inherited_original_classes_by_id: dict[str, str] = field(default_factory=dict, repr=False)
    removal_provenance_by_id: dict[str, dict[str, str]] = field(default_factory=dict, repr=False)
    protein_protected_ids: set[str] = field(default_factory=set, repr=False)
    transcript_protected_ids: set[str] = field(default_factory=set, repr=False)
    overlap_protected_ids: set[str] = field(default_factory=set, repr=False)
    primary_route_counts: dict[str, int] = field(default_factory=dict)
    native_primary_route_counts: dict[str, int] = field(default_factory=dict)
    inherited_original_route_counts: dict[str, int] = field(default_factory=dict)
    removal_provenance_route_counts: dict[str, int] = field(default_factory=dict)
    all_triggered_route_counts: dict[str, int] = field(default_factory=dict)
    removal_by_source: dict[str, int] = field(default_factory=dict)
    source_selected_counts: dict[str, int] = field(default_factory=dict)
    source_removed_fractions: dict[str, float] = field(default_factory=dict)
    source_share_all_removals: dict[str, float] = field(default_factory=dict)
    dimension_counts: dict[str, int] = field(default_factory=dict)
    protection_counts: dict[str, int] = field(default_factory=dict)
    incremental_removed_ids: set[str] = field(default_factory=set, repr=False)
    incremental_removals: int = 0
    incremental_fraction_original: float = 0.0
    incremental_fraction_previous_retained: float = 0.0
    incremental_primary_route_counts: dict[str, int] = field(default_factory=dict)
    incremental_triggered_route_counts: dict[str, int] = field(default_factory=dict)
    incremental_by_source: dict[str, int] = field(default_factory=dict)
    source_share_incremental_removals: dict[str, float] = field(default_factory=dict)
    thresholds: str = "."
    low_score_nonsingleton_scope: str = "."
    sequence_integrity_nonsingleton_scope: str = "."


def removal_provenance_from_result(
    result: Mapping[str, object],
    *,
    level: str,
) -> dict[str, dict[str, str]]:
    """Return the original biological/pruning provenance for removed models."""
    provenance: dict[str, dict[str, str]] = {}
    for row in result.get("removed_rows", []):  # type: ignore[assignment]
        model_id = str(row.get("output_transcript_id", "."))
        inherited = _truthy(row.get("final_monotonic_inherited_removal"))
        provenance[model_id] = {
            "original_level": str(
                row.get("final_monotonic_inherited_original_level", level)
                if inherited else level
            ),
            "original_reason": str(
                row.get(
                    "final_monotonic_inherited_original_primary_reason",
                    primary_reason_for_row(row),
                )
                if inherited else primary_reason_for_row(row)
            ),
            "original_prune_classes": str(
                row.get(
                    "final_monotonic_inherited_original_prune_classes",
                    row.get("final_prune_classes", "."),
                )
                if inherited else row.get("final_prune_classes", ".")
            ),
        }
    return provenance


def apply_monotonic_pruning_closure(
    result: dict[str, object],
    *,
    required_removed_ids: set[str],
    lower_level: str,
    current_level: str,
    lower_provenance_by_id: Mapping[str, Mapping[str, str]],
) -> dict[str, object]:
    """Inherit lower-level removals that a stronger native scenario restored.

    Native pruning is still evaluated independently against the immutable
    pre-pruning population.  This closure is applied afterwards so the fixed
    preset contract remains monotonic while preserving the original lower-level
    removal provenance for audit.
    """
    retained_models = list(result.get("retained_models", []))  # type: ignore[arg-type]
    retained_rows = list(result.get("retained_rows", []))  # type: ignore[arg-type]
    removed_rows = [dict(row) for row in result.get("removed_rows", [])]  # type: ignore[arg-type]
    retained_info = dict(result.get("retained_selection_info_by_transcript", {}))  # type: ignore[arg-type]
    row_by_id = {str(row.get("output_transcript_id", ".")): dict(row) for row in retained_rows}
    model_by_id = {_model_id(model): model for model in retained_models}
    native_removed_ids = {
        str(row.get("output_transcript_id", ".")) for row in removed_rows
    }
    inherited_ids = sorted((required_removed_ids - native_removed_ids) & set(model_by_id))

    for row in removed_rows:
        row["final_monotonic_inherited_removal"] = "false"
        row["final_monotonic_inherited_from_level"] = "."
        row["final_monotonic_inherited_original_level"] = "."
        row["final_monotonic_inherited_original_primary_reason"] = "."
        row["final_monotonic_inherited_original_prune_classes"] = "."
        row["final_monotonic_native_status"] = "removed"

    inherited_set = set(inherited_ids)
    retained_models_out = [
        model for model in retained_models if _model_id(model) not in inherited_set
    ]
    retained_rows_out: list[dict[str, object]] = []
    for row in retained_rows:
        model_id = str(row.get("output_transcript_id", "."))
        if model_id in inherited_set:
            continue
        out_row = dict(row)
        out_row["final_monotonic_inherited_removal"] = "false"
        out_row["final_monotonic_inherited_from_level"] = "."
        out_row["final_monotonic_inherited_original_level"] = "."
        out_row["final_monotonic_inherited_original_primary_reason"] = "."
        out_row["final_monotonic_inherited_original_prune_classes"] = "."
        out_row["final_monotonic_native_status"] = "retained"
        retained_rows_out.append(out_row)

    for model_id in inherited_ids:
        row = dict(row_by_id[model_id])
        provenance = lower_provenance_by_id.get(model_id, {})
        row["final_prune_level"] = current_level
        classes = [
            item for item in str(row.get("final_prune_classes", ".")).split(",")
            if item and item != "."
        ]
        if "inherited_lower_level_removal" not in classes:
            classes.append("inherited_lower_level_removal")
        row["final_prune_classes"] = ",".join(classes)
        row.update({
            "final_monotonic_inherited_removal": "true",
            "final_monotonic_inherited_from_level": lower_level,
            "final_monotonic_inherited_original_level": provenance.get(
                "original_level", lower_level
            ),
            "final_monotonic_inherited_original_primary_reason": provenance.get(
                "original_reason", "."
            ),
            "final_monotonic_inherited_original_prune_classes": provenance.get(
                "original_prune_classes", "."
            ),
            "final_monotonic_native_status": "retained",
        })
        removed_rows.append(row)
        retained_info.pop(model_id, None)

    out = dict(result)
    out.update({
        "retained_models": retained_models_out,
        "retained_rows": retained_rows_out,
        "removed_rows": removed_rows,
        "retained_selection_info_by_transcript": retained_info,
        "monotonic_closure_applied": bool(inherited_ids),
        "monotonic_closure_inherited_count": len(inherited_ids),
        "monotonic_closure_inherited_ids": tuple(inherited_ids),
        "monotonic_closure_inherited_from_level": lower_level,
    })
    return out


def _exact_structure_supported(row: Mapping[str, object]) -> bool:
    return (
        value_as_float(row.get("exact_structure_support_score", 0.0), 0.0) > 0
        or value_as_float(row.get("exact_cds_exon_chain_support_source_count", 0), 0.0) > 0
        or value_as_float(row.get("single_exon_cds_exact_support_source_count", 0), 0.0) > 0
    )


def summarize_pruning_scenario(
    level: str,
    result: Mapping[str, object],
    *,
    initial_models: Sequence[CodingModel],
    subloci: Sequence[Sublocus],
    structural_rows_by_id: Mapping[str, Mapping[str, object]],
    protein_adequacy_state: str,
    transcript_adequacy_state: str,
    protein_evidence_supplied: bool,
    transcript_evidence_supplied: bool,
    scenario_configuration: Mapping[str, object],
    thresholds: str,
    low_score_nonsingleton_scope: str,
    sequence_integrity_nonsingleton_scope: str,
) -> ScenarioSummary:
    """Create a compact, stable summary from an existing pruning result."""
    retained_rows = list(result.get("retained_rows", []))  # type: ignore[arg-type]
    removed_rows = list(result.get("removed_rows", []))  # type: ignore[arg-type]
    rows = retained_rows + removed_rows
    row_by_id = {str(row.get("output_transcript_id")): row for row in rows}
    model_by_id = {_model_id(model): model for model in initial_models}
    sub_by_id = {model.output_transcript_id: sub for sub in subloci for model in sub.models}
    initial_ids = set(model_by_id)
    retained_ids = {str(row.get("output_transcript_id")) for row in retained_rows}
    removed_ids = {str(row.get("output_transcript_id")) for row in removed_rows}
    inherited_ids = {
        model_id for model_id in removed_ids
        if _truthy(row_by_id[model_id].get("final_monotonic_inherited_removal"))
    }
    native_removed_ids = removed_ids - inherited_ids
    inherited_from_level_by_id = {
        model_id: str(row_by_id[model_id].get("final_monotonic_inherited_from_level", "."))
        for model_id in inherited_ids
    }
    inherited_original_level_by_id = {
        model_id: str(row_by_id[model_id].get("final_monotonic_inherited_original_level", "."))
        for model_id in inherited_ids
    }
    inherited_original_reason_by_id = {
        model_id: str(
            row_by_id[model_id].get(
                "final_monotonic_inherited_original_primary_reason", "."
            )
        )
        for model_id in inherited_ids
    }
    inherited_original_classes_by_id = {
        model_id: str(
            row_by_id[model_id].get(
                "final_monotonic_inherited_original_prune_classes", "."
            )
        )
        for model_id in inherited_ids
    }
    primary_by_id = {model_id: primary_reason_for_row(row_by_id[model_id]) for model_id in removed_ids}
    triggered_by_id: dict[str, tuple[str, ...]] = {}
    primary_counts: Counter[str] = Counter(primary_by_id.values())
    native_primary_counts: Counter[str] = Counter(
        primary_by_id[model_id] for model_id in native_removed_ids
    )
    inherited_original_counts: Counter[str] = Counter(
        inherited_original_reason_by_id.values()
    )
    provenance_counts: Counter[str] = Counter(native_primary_counts)
    provenance_counts.update(inherited_original_counts)
    triggered_counts: Counter[str] = Counter()
    external_adequacy_passed = (
        protein_adequacy_state == "pass" or transcript_adequacy_state == "pass"
    )
    singleton_thresholds = scenario_configuration["singleton_thresholds"]
    nonsingleton_thresholds = scenario_configuration["nonsingleton_thresholds"]
    for model_id in removed_ids:
        if model_id in inherited_ids:
            classes = ("inherited_lower_level_removal",)
            triggered_by_id[model_id] = classes
            triggered_counts.update(classes)
            continue
        row = row_by_id[model_id]
        model = model_by_id[model_id]
        sub = sub_by_id.get(model_id)
        structural = structural_rows_by_id.get(model_id, row)
        is_singleton = bool(sub and len(sub.models) == 1)
        is_backfill = str(row.get("final_selection_stage", ".")) == "backfill"
        triggers = {
            item for item in str(row.get("final_prune_classes", ".")).split(",")
            if item and item != "."
        }
        nonsingleton_prune_route = str(
            row.get("final_low_score_nonsingleton_prune_route", ".")
        )
        if nonsingleton_prune_route != ".":
            triggers.add(nonsingleton_prune_route)
        singleton_level = str(scenario_configuration["singleton_level"])
        intrinsic = (
            value_as_float(structural.get("intrinsic_penalty_score", 0.0), 0.0) > 0
            and str(structural.get("intrinsic_penalty_components", ".")) != "."
        )
        sequence_concern = _truthy(row.get("sequence_integrity_pruning_concern"))
        severe_sequence_concern = _truthy(row.get("sequence_integrity_severe_pruning_concern"))
        if is_singleton and singleton_level != "none":
            singleton_trigger = (
                singleton_level == "all"
                or (singleton_level == "intrinsic" and intrinsic)
                or (singleton_level == "sequence_integrity" and sequence_concern)
                or (
                    singleton_level == "intrinsic_or_sequence_integrity"
                    and (intrinsic or sequence_concern)
                )
            )
            if singleton_trigger:
                triggers.add("singleton_sublocus")
                if sequence_concern:
                    triggers.add("sequence_integrity_singleton")

        model_score = value_as_float(row.get("model_structural_support_score_scaled", 0.0), 0.0)
        protein_score = value_as_float(row.get("protein_evidence_score_scaled", 0.0), 0.0)
        transcript_score = value_as_float(row.get("transcript_evidence_score_scaled", 0.0), 0.0)
        model_selection_score = value_as_float(row.get("model_selection_score_scaled", 0.0), 0.0)
        context_concern = _truthy(row.get("transcript_cdna_context_pruning_concern"))

        def score_or_context_trigger(threshold_values: object) -> tuple[bool, bool]:
            threshold = threshold_values  # four-item fixed scenario tuple
            model_low = model_score <= threshold[0] + 1e-15  # type: ignore[index]
            protein_low = (not protein_evidence_supplied) or protein_score <= threshold[1] + 1e-15  # type: ignore[index]
            transcript_low = (not transcript_evidence_supplied) or transcript_score <= threshold[2] + 1e-15  # type: ignore[index]
            model_selection_low = model_selection_score <= threshold[3] + 1e-15  # type: ignore[index]
            score_trigger = bool(
                external_adequacy_passed
                and model_low and protein_low and transcript_low and model_selection_low
            )
            context_trigger = bool(
                external_adequacy_passed
                and transcript_evidence_supplied
                and transcript_adequacy_state == "pass"
                and context_concern and model_low and protein_low
            )
            return score_trigger, context_trigger

        if is_singleton and scenario_configuration["low_score_singletons"]:
            score_trigger, context_trigger = score_or_context_trigger(singleton_thresholds)
            if score_trigger:
                triggers.add("low_score_singleton_score_thresholds")
            if context_trigger:
                triggers.add("transcript_cdna_context_singleton")
        if not is_singleton and scenario_configuration["low_score_nonsingletons"]:
            in_scope = (
                scenario_configuration["low_score_nonsingleton_scope"] == "all"
                or is_backfill
            )
            if in_scope:
                score_trigger, context_trigger = score_or_context_trigger(nonsingleton_thresholds)
                if score_trigger:
                    triggers.add("low_score_nonsingleton_score_thresholds")
                if context_trigger:
                    triggers.add("transcript_cdna_context_nonsingleton")
        if not is_singleton and scenario_configuration["sequence_integrity_nonsingletons"]:
            sequence_scope = str(scenario_configuration["sequence_integrity_nonsingleton_scope"])
            if (is_backfill and sequence_concern) or (
                sequence_scope == "all" and severe_sequence_concern
            ):
                triggers.add("sequence_integrity_nonsingleton")
        classes = tuple(sorted(triggers))
        triggered_by_id[model_id] = classes
        triggered_counts.update(classes)

    protein_ids: set[str] = set()
    transcript_ids: set[str] = set()
    overlap_ids: set[str] = set()
    protection_counts: Counter[str] = Counter()
    candidate_ids: set[str] = set(removed_ids)
    for model_id in retained_ids:
        row = row_by_id[model_id]
        protein, transcript, overlap = protection_flags(row)
        if protein:
            protein_ids.add(model_id)
            protection_counts["protein"] += 1
        if transcript:
            transcript_ids.add(model_id)
            protection_counts["transcript"] += 1
        if protein and transcript:
            protection_counts["protein_and_transcript"] += 1
        if overlap:
            overlap_ids.add(model_id)
            protection_counts["monosublocus_overlap"] += 1
        if protein or transcript or overlap:
            candidate_ids.add(model_id)
        if _truthy(row.get("final_protein_singleton_rescued")):
            reason = str(row.get("final_protein_singleton_rescue_applied_reason", "."))
            protection_counts["protein_exact" if "exact" in reason else "protein_score"] += 1
            protection_counts[f"protein_adequacy_{protein_adequacy_state}"] += 1
        if _truthy(row.get("final_transcript_singleton_rescued")):
            reason = str(row.get("final_transcript_singleton_rescue_applied_reason", "."))
            protection_counts["transcript_exact" if "exact" in reason else "transcript_score"] += 1
            protection_counts[f"transcript_adequacy_{transcript_adequacy_state}"] += 1
            context_state = str(row.get("final_transcript_singleton_rescue_cdna_context_fraction_state", "."))
            if context_state in {"not_evaluable", "not_calculated"}:
                protection_counts["transcript_context_not_evaluable"] += 1
            else:
                fraction = value_as_float(
                    row.get("final_transcript_singleton_rescue_cdna_context_fraction", 0.0), 0.0
                )
                protection_counts[
                    "transcript_context_ge_0.30" if fraction >= 0.30 - 1e-15
                    else "transcript_context_lt_0.30"
                ] += 1
        decision = str(row.get("final_singleton_overlap_rescue_decision", "."))
        if decision != ".":
            protection_counts[f"overlap_{decision}"] += 1
    for key, label in (
        ("singleton_overlap_rescue_overlap_bp_geometry_candidate_count", "overlap_bp_geometry_candidates"),
        ("singleton_overlap_rescue_reciprocal_rejected_count", "overlap_reciprocal_rejected"),
        ("singleton_overlap_rescue_models_with_same_source_edge_count", "models_with_same_source_edge"),
        ("singleton_overlap_rescue_geometry_candidate_count", "overlap_geometry_candidates"),
        ("singleton_overlap_rescue_policy_denied_count", "overlap_denied_adequate_evidence"),
        ("singleton_overlap_rescued_by_strong_local_evidence_count", "overlap_rescued_strong_local_evidence"),
        ("singleton_overlap_rescued_without_adequate_evidence_count", "overlap_rescued_no_adequate_evidence"),
        ("singleton_rescued_by_protein_evidence_count", "protein_rescued"),
        ("singleton_rescued_by_transcript_evidence_count", "transcript_rescued"),
    ):
        protection_counts[label] = int(value_as_float(result.get(key, 0), 0.0))

    source_by_id = {model_id: model.source_name for model_id, model in model_by_id.items()}
    source_selected = Counter(model.source_name for model in initial_models)
    source_removed = Counter(model_by_id[model_id].source_name for model_id in removed_ids)
    source_removed_fractions = {
        source: _safe_fraction(source_removed[source], count)
        for source, count in sorted(source_selected.items())
    }
    source_share = {
        source: _safe_fraction(count, len(removed_ids))
        for source, count in sorted(source_removed.items())
    }
    dimensions: Counter[str] = Counter()
    removed_sublocus_ids: set[str] = set()
    retained_sublocus_ids: set[str] = set()
    for model_id in retained_ids:
        sub = sub_by_id.get(model_id)
        if sub:
            retained_sublocus_ids.add(sub.sublocus_id)
    for model_id in removed_ids:
        model = model_by_id[model_id]
        sub = sub_by_id.get(model_id)
        row = structural_rows_by_id.get(model_id, row_by_id[model_id])
        if sub:
            removed_sublocus_ids.add(sub.sublocus_id)
        dimensions["singleton" if sub and len(sub.models) == 1 else "nonsingleton"] += 1
        dimensions["single_exon" if model.is_single_exon else "multi_exon"] += 1
        stage = str(row_by_id[model_id].get("final_selection_stage", "."))
        dimensions["backfill" if stage == "backfill" else "nonbackfill"] += 1
        dimensions[f"selection_stage:{stage}"] += 1
        source_count = len({member.source_name for member in sub.models}) if sub else 1
        dimensions["single_source" if source_count == 1 else "multisource"] += 1
        dimensions["exact_structure_supported" if _exact_structure_supported(row) else "exact_structure_unsupported"] += 1

    return ScenarioSummary(
        level=level,
        starting_models=len(initial_ids),
        candidate_models=len(candidate_ids),
        protected_models=len(protein_ids | transcript_ids | overlap_ids),
        removed_models=len(removed_ids),
        retained_models=len(retained_ids),
        removed_fraction=_safe_fraction(len(removed_ids), len(initial_ids)),
        affected_subloci=len(removed_sublocus_ids),
        retained_subloci=len(retained_sublocus_ids),
        removed_ids=removed_ids,
        retained_ids=retained_ids,
        native_removed_ids=native_removed_ids,
        inherited_removed_ids=inherited_ids,
        primary_reason_by_id=primary_by_id,
        all_triggered_by_id=triggered_by_id,
        source_by_id=source_by_id,
        inherited_from_level_by_id=inherited_from_level_by_id,
        inherited_original_level_by_id=inherited_original_level_by_id,
        inherited_original_reason_by_id=inherited_original_reason_by_id,
        inherited_original_classes_by_id=inherited_original_classes_by_id,
        removal_provenance_by_id=removal_provenance_from_result(result, level=level),
        protein_protected_ids=protein_ids,
        transcript_protected_ids=transcript_ids,
        overlap_protected_ids=overlap_ids,
        primary_route_counts=dict(sorted(primary_counts.items())),
        native_primary_route_counts=dict(sorted(native_primary_counts.items())),
        inherited_original_route_counts=dict(sorted(inherited_original_counts.items())),
        removal_provenance_route_counts=dict(sorted(provenance_counts.items())),
        all_triggered_route_counts=dict(sorted(triggered_counts.items())),
        removal_by_source=dict(sorted(source_removed.items())),
        source_selected_counts=dict(sorted(source_selected.items())),
        source_removed_fractions=source_removed_fractions,
        source_share_all_removals=source_share,
        dimension_counts=dict(sorted(dimensions.items())),
        protection_counts=dict(sorted(protection_counts.items())),
        thresholds=thresholds,
        low_score_nonsingleton_scope=low_score_nonsingleton_scope,
        sequence_integrity_nonsingleton_scope=sequence_integrity_nonsingleton_scope,
    )


def add_incremental_metrics(
    scenarios: Mapping[str, ScenarioSummary],
) -> None:
    previous: ScenarioSummary | None = None
    for level in AUTO_SCENARIO_LEVELS:
        scenario = scenarios[level]
        if previous is None:
            previous = scenario
            continue
        incremental = scenario.removed_ids - previous.removed_ids
        scenario.incremental_removed_ids = incremental
        scenario.incremental_removals = len(incremental)
        scenario.incremental_fraction_original = _safe_fraction(
            len(incremental), scenario.starting_models
        )
        scenario.incremental_fraction_previous_retained = _safe_fraction(
            len(incremental), previous.retained_models
        )
        scenario.incremental_primary_route_counts = dict(sorted(Counter(
            scenario.primary_reason_by_id.get(model_id, ".") for model_id in incremental
        ).items()))
        all_counts: Counter[str] = Counter()
        for model_id in incremental:
            all_counts.update(scenario.all_triggered_by_id.get(model_id, ()))
        scenario.incremental_triggered_route_counts = dict(sorted(all_counts.items()))
        by_source = Counter(
            scenario.source_by_id.get(model_id, ".") for model_id in incremental
        )
        scenario.incremental_by_source = dict(sorted(by_source.items()))
        scenario.source_share_incremental_removals = {
            source: _safe_fraction(count, len(incremental))
            for source, count in sorted(by_source.items())
        }
        previous = scenario


@dataclass(frozen=True)
class ScenarioValidation:
    reconciliation_passed: bool
    medium_high_monotonicity_passed: bool
    high_very_high_monotonicity_passed: bool
    unexpected_restored_medium_to_high: int
    unexpected_restored_high_to_very_high: int
    native_restored_medium_to_high: int = 0
    native_restored_high_to_very_high: int = 0

    @property
    def passed(self) -> bool:
        return (
            self.reconciliation_passed
            and self.medium_high_monotonicity_passed
            and self.high_very_high_monotonicity_passed
        )


def validate_scenarios(scenarios: Mapping[str, ScenarioSummary]) -> ScenarioValidation:
    reconciliation = all(
        scenario.starting_models == scenario.removed_models + scenario.retained_models
        and not (scenario.removed_ids & scenario.retained_ids)
        for scenario in scenarios.values()
    )
    medium = scenarios["medium"]
    high = scenarios["high"]
    very_high = scenarios["very-high"]
    restored_mh = high.retained_ids - medium.retained_ids
    restored_hv = very_high.retained_ids - high.retained_ids
    return ScenarioValidation(
        reconciliation,
        not restored_mh,
        not restored_hv,
        len(restored_mh),
        len(restored_hv),
        len(high.inherited_removed_ids),
        len(very_high.inherited_removed_ids),
    )


@dataclass(frozen=True)
class AutoDecision:
    mode: str
    unconstrained_level: str
    effective_level: str
    confidence: str
    reason_codes: tuple[str, ...]
    warning_codes: tuple[str, ...]
    user_min_level: str
    user_max_level: str
    backbone_class_used: str = "."
    used_flood_adjustment: bool = False


def _clamp_level(level: str, minimum: str, maximum: str) -> str:
    order = LEVEL_ORDER[level]
    if order < LEVEL_ORDER[minimum]:
        return minimum
    if order > LEVEL_ORDER[maximum]:
        return maximum
    return level


def choose_auto_pruning_level(
    mode: str,
    *,
    scenarios: Mapping[str, ScenarioSummary],
    validation: ScenarioValidation,
    consensus: SourceConsensusMetrics,
    adequacy_pass_count: int,
    selected_singleton_model_fraction: float,
    adjusted_consensus: SourceConsensusMetrics | None = None,
    min_level: str = "medium",
    max_level: str = "very-high",
    max_level_explicit: bool = False,
    max_extra_removal_fraction: float = 0.05,
    high_min_singleton_fraction: float = 0.10,
    high_min_prospective_removal_fraction: float = 0.02,
    use_source_excess: bool = True,
) -> AutoDecision:
    """Choose a fixed pruning scenario using the agreed v0.2.14 rules.

    Standard ``auto`` may use an eligible flood-adjusted backbone only on the
    no-adequate-evidence pathway.  ``auto-conservative`` and both modes' fully
    evidence-driven pathways retain the unadjusted backbone semantics.
    """
    reasons: list[str] = []
    warnings: list[str] = []
    decision_consensus = consensus
    used_flood_adjustment = False
    if (
        adequacy_pass_count == 0
        and mode == "auto"
        and use_source_excess
        and adjusted_consensus is not None
    ):
        decision_consensus = adjusted_consensus
        used_flood_adjustment = True
    if not validation.passed:
        unconstrained = "medium"
        confidence = "low"
        reasons.append("SCENARIO_VALIDATION_FAILED_FALLBACK_MEDIUM")
        warnings.append("SCENARIO_RECONCILIATION_OR_MONOTONICITY_FAILURE")
        # Consistency failure is a hard safety fallback; user escalation
        # constraints do not override it.
        effective = "medium"
        return AutoDecision(
            mode, unconstrained, effective, confidence,
            tuple(reasons), tuple(warnings), min_level, max_level,
            decision_consensus.backbone_class, used_flood_adjustment,
        )

    if (
        validation.native_restored_medium_to_high
        or validation.native_restored_high_to_very_high
    ):
        warnings.append("MONOTONIC_CLOSURE_APPLIED")

    high = scenarios["high"]
    very_high = scenarios["very-high"]
    high_increment_meaningful = (
        high.incremental_removals >= 10
        and high.incremental_fraction_original >= 0.0005 - 1e-15
    )
    vh_increment_meaningful = (
        very_high.incremental_removals >= 10
        and very_high.incremental_fraction_original >= 0.0005 - 1e-15
    )
    confidence = "moderate"

    if adequacy_pass_count == 0:
        if decision_consensus.effective_source_count < 3:
            unconstrained = "medium"
            reasons.append("NO_ADEQUATE_EVIDENCE_FEWER_THAN_THREE_MODEL_SOURCES")
        elif decision_consensus.backbone_class == "weak_or_discordant":
            unconstrained = "medium"
            reasons.append("NO_ADEQUATE_EVIDENCE_WEAK_OR_DISCORDANT_BACKBONE")
            warnings.append("BALANCED_OR_WEAK_MODEL_SOURCE_CONSENSUS")
        elif decision_consensus.backbone_class == "strong":
            unconstrained = "high"
            reasons.append("NO_ADEQUATE_EVIDENCE_STRONG_CONSENSUS_BACKBONE")
            if used_flood_adjustment:
                reasons.append("FLOOD_ADJUSTED_CONSENSUS_BACKBONE_USED")
            confidence = "high" if decision_consensus.effective_source_count >= 5 else "moderate"
            if decision_consensus.effective_source_count >= 5:
                warnings.append("SOURCE_INDEPENDENCE_NOT_MODELLED")
        elif mode == "auto-conservative":
            unconstrained = "medium"
            reasons.append("INTERMEDIATE_BACKBONE_CONSERVATIVE_MODE")
        else:
            excess_sources = consensus.strong_source_excess_sources if use_source_excess else ()
            concentration = 0.0
            if len(excess_sources) == 1 and high.incremental_removals:
                concentration = _safe_fraction(
                    high.incremental_by_source.get(excess_sources[0], 0),
                    high.incremental_removals,
                )
            if (
                len(excess_sources) == 1
                and high_increment_meaningful
                and concentration >= 0.70 - 1e-15
            ):
                unconstrained = "high"
                reasons.extend((
                    "INTERMEDIATE_BACKBONE_SINGLE_STRONG_SOURCE_EXCESS",
                    "HIGH_INCREMENT_CONCENTRATED_IN_EXCESS_SOURCE",
                ))
                confidence = "moderate"
            else:
                unconstrained = "medium"
                reasons.append("INTERMEDIATE_BACKBONE_ESCALATION_CRITERIA_NOT_MET")
        if LEVEL_ORDER[unconstrained] > LEVEL_ORDER["high"]:
            unconstrained = "high"
        reasons.append("NO_ADEQUATE_EVIDENCE_AUTOMATIC_MAXIMUM_HIGH")
    else:
        conservative = mode == "auto-conservative"
        singleton_threshold = 0.15 if conservative else high_min_singleton_fraction
        prospective_threshold = 0.03 if conservative else high_min_prospective_removal_fraction
        singleton_trigger = selected_singleton_model_fraction >= singleton_threshold - 1e-15
        prospective_trigger = (
            high.incremental_removals >= 10
            and high.incremental_fraction_original >= prospective_threshold - 1e-15
        )
        if singleton_trigger or prospective_trigger:
            unconstrained = "high"
            if singleton_trigger:
                reasons.append("ADEQUATE_EVIDENCE_SINGLETON_BURDEN_HIGH_TRIGGER")
            if prospective_trigger:
                reasons.append("ADEQUATE_EVIDENCE_PROSPECTIVE_HIGH_INCREMENT_TRIGGER")
        else:
            unconstrained = "medium"
            reasons.append("ADEQUATE_EVIDENCE_HIGH_TRIGGER_NOT_MET")

        may_consider_vh = (
            unconstrained == "high"
            and (not conservative or (max_level_explicit and max_level == "very-high"))
        )
        vh_limit = 0.02 if conservative else max_extra_removal_fraction
        if may_consider_vh and vh_increment_meaningful:
            if very_high.incremental_fraction_previous_retained <= vh_limit + 1e-15:
                unconstrained = "very-high"
                reasons.append("VERY_HIGH_INCREMENT_WITHIN_SAFETY_LIMIT")
            else:
                warnings.append("VERY_HIGH_EXTRA_REMOVAL_FRACTION_EXCEEDS_LIMIT")
                reasons.append("VERY_HIGH_REJECTED_BY_SAFETY_LIMIT")
        elif conservative and unconstrained == "high" and not (
            max_level_explicit and max_level == "very-high"
        ):
            reasons.append("AUTO_CONSERVATIVE_DEFAULT_MAXIMUM_HIGH")
        elif unconstrained == "high" and not vh_increment_meaningful:
            reasons.append("VERY_HIGH_INCREMENT_NOT_MEANINGFUL")
        confidence = "high" if adequacy_pass_count == 2 and len(reasons) > 1 else "moderate"

    effective = _clamp_level(unconstrained, min_level, max_level)
    if effective != unconstrained:
        reasons.append("USER_AUTO_LEVEL_CONSTRAINT_APPLIED")
    return AutoDecision(
        mode, unconstrained, effective, confidence,
        tuple(dict.fromkeys(reasons)), tuple(dict.fromkeys(warnings)),
        min_level, max_level,
        decision_consensus.backbone_class, used_flood_adjustment,
    )


SCENARIO_SUMMARY_COLUMNS = (
    "level", "thresholds", "low_score_nonsingleton_scope",
    "sequence_integrity_nonsingleton_scope", "starting_models",
    "candidate_models", "protected_models", "removed_models", "retained_models",
    "native_removed_models", "inherited_removed_models",
    "removed_fraction", "affected_subloci", "retained_subloci",
    "incremental_removals", "incremental_fraction_original",
    "incremental_fraction_previous_retained", "primary_route_counts",
    "native_primary_route_counts", "inherited_original_route_counts",
    "removal_provenance_route_counts",
    "all_triggered_route_counts", "incremental_primary_route_counts",
    "incremental_all_triggered_route_counts", "removal_by_source",
    "source_selected_counts", "source_removed_fractions",
    "source_share_all_removals", "incremental_by_source",
    "source_share_incremental_removals", "dimension_counts", "protection_counts",
)


def scenario_summary_row(scenario: ScenarioSummary) -> dict[str, object]:
    return {
        "level": scenario.level,
        "thresholds": scenario.thresholds,
        "low_score_nonsingleton_scope": scenario.low_score_nonsingleton_scope,
        "sequence_integrity_nonsingleton_scope": scenario.sequence_integrity_nonsingleton_scope,
        "starting_models": scenario.starting_models,
        "candidate_models": scenario.candidate_models,
        "protected_models": scenario.protected_models,
        "removed_models": scenario.removed_models,
        "retained_models": scenario.retained_models,
        "native_removed_models": len(scenario.native_removed_ids),
        "inherited_removed_models": len(scenario.inherited_removed_ids),
        "removed_fraction": f"{scenario.removed_fraction:.6f}",
        "affected_subloci": scenario.affected_subloci,
        "retained_subloci": scenario.retained_subloci,
        "incremental_removals": scenario.incremental_removals,
        "incremental_fraction_original": f"{scenario.incremental_fraction_original:.6f}",
        "incremental_fraction_previous_retained": f"{scenario.incremental_fraction_previous_retained:.6f}",
        "primary_route_counts": _serialize_counter(scenario.primary_route_counts),
        "native_primary_route_counts": _serialize_counter(scenario.native_primary_route_counts),
        "inherited_original_route_counts": _serialize_counter(scenario.inherited_original_route_counts),
        "removal_provenance_route_counts": _serialize_counter(scenario.removal_provenance_route_counts),
        "all_triggered_route_counts": _serialize_counter(scenario.all_triggered_route_counts),
        "incremental_primary_route_counts": _serialize_counter(scenario.incremental_primary_route_counts),
        "incremental_all_triggered_route_counts": _serialize_counter(scenario.incremental_triggered_route_counts),
        "removal_by_source": _serialize_counter(scenario.removal_by_source),
        "source_selected_counts": _serialize_counter(scenario.source_selected_counts),
        "source_removed_fractions": _serialize_float_map(scenario.source_removed_fractions),
        "source_share_all_removals": _serialize_float_map(scenario.source_share_all_removals),
        "incremental_by_source": _serialize_counter(scenario.incremental_by_source),
        "source_share_incremental_removals": _serialize_float_map(scenario.source_share_incremental_removals),
        "dimension_counts": _serialize_counter(scenario.dimension_counts),
        "protection_counts": _serialize_counter(scenario.protection_counts),
    }


def write_scenario_summary(path: Path, scenarios: Mapping[str, ScenarioSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SCENARIO_SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for level in AUTO_SCENARIO_LEVELS:
            writer.writerow(scenario_summary_row(scenarios[level]))


AUTO_DIAGNOSTIC_COLUMNS = ("section", "metric", "value")


def automatic_diagnostics_rows(
    *,
    protein: EvidenceAdequacy,
    transcript: EvidenceAdequacy,
    consensus: SourceConsensusMetrics,
    selected_singleton_model_fraction: float,
    validation: ScenarioValidation,
    standard_decision: AutoDecision,
    conservative_decision: AutoDecision,
    flood_assessment: FloodAdjustedConsensusAssessment | None = None,
    input_count_metrics: Mapping[str, SourceInputCountMetric] | None = None,
    scenarios: Mapping[str, ScenarioSummary] | None = None,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def add(section: str, metric: str, value: object) -> None:
        if isinstance(value, float):
            value = "inf" if math.isinf(value) else f"{value:.6f}"
        rows.append({"section": section, "metric": metric, "value": value})

    for adequacy in (protein, transcript):
        prefix = adequacy.evidence_class
        add("evidence_adequacy", f"{prefix}_state", adequacy.state)
        add("evidence_adequacy", f"{prefix}_first_value", adequacy.first_value if adequacy.first_value is not None else ".")
        add("evidence_adequacy", f"{prefix}_second_value", adequacy.second_value if adequacy.second_value is not None else ".")
        add("evidence_adequacy", f"{prefix}_first_threshold", adequacy.first_threshold)
        add("evidence_adequacy", f"{prefix}_second_threshold", adequacy.second_threshold)
        add("evidence_adequacy", f"{prefix}_first_margin", adequacy.first_margin if adequacy.first_margin is not None else ".")
        add("evidence_adequacy", f"{prefix}_second_margin", adequacy.second_margin if adequacy.second_margin is not None else ".")
        add("evidence_adequacy", f"{prefix}_limiting_margin", adequacy.limiting_margin if adequacy.limiting_margin is not None else ".")
    add("evidence_adequacy", "adequacy_pass_count", int(protein.passed) + int(transcript.passed))

    add("source_consensus", "source_count_distribution", _serialize_counter({
        str(key): value for key, value in consensus.source_count_distribution.items()
    }))
    for name in (
        "total_subloci", "total_source_incidences", "raw_source_count",
        "effective_source_count", "multisource_subloci",
        "singleton_sublocus_fraction", "multisource_sublocus_fraction",
        "deep_consensus_sublocus_fraction", "full_consensus_sublocus_fraction",
        "singleton_incidence_fraction", "consensus_incidence_fraction",
        "deep_consensus_incidence_fraction", "corroboration_fraction",
        "mean_sources_per_sublocus", "mean_sources_at_multisource_subloci",
        "backbone_class",
    ):
        add("source_consensus", name, getattr(consensus, name))
    add("source_consensus", "strong_source_excess_sources", ",".join(consensus.strong_source_excess_sources) or ".")
    for source in sorted(consensus.per_source_total_subloci):
        add("per_source", f"{source}.singleton_count", consensus.per_source_singleton_count.get(source, 0))
        add("per_source", f"{source}.total_subloci", consensus.per_source_total_subloci[source])
        add("per_source", f"{source}.singleton_share", consensus.per_source_singleton_share[source])
        add("per_source", f"{source}.singleton_rate", consensus.per_source_singleton_rate[source])
        add("per_source", f"{source}.other_source_median_rate", consensus.per_source_other_median_rate[source])
        add("per_source", f"{source}.singleton_rate_difference", consensus.per_source_singleton_rate_difference[source])
        add("per_source", f"{source}.singleton_rate_ratio", consensus.per_source_singleton_rate_ratio[source])

    # Existing source_consensus row names above deliberately retain their raw
    # meanings.  The v0.2.14 rows below are append-only and explicit.
    add("source_consensus_v0.2.14", "unadjusted_backbone_class", consensus.backbone_class)
    add(
        "source_consensus_v0.2.14",
        "strong_deep_consensus_incidence_threshold",
        strong_deep_consensus_incidence_threshold(consensus.effective_source_count),
    )
    for name, criterion in consensus_classification_criteria(consensus).items():
        add("source_consensus_unadjusted_criteria", f"{name}.value", criterion["value"])
        add("source_consensus_unadjusted_criteria", f"{name}.operator", criterion["operator"])
        add("source_consensus_unadjusted_criteria", f"{name}.threshold", criterion["threshold"])
        add("source_consensus_unadjusted_criteria", f"{name}.margin", criterion["margin"])
        add(
            "source_consensus_unadjusted_criteria",
            f"{name}.passed",
            "true" if criterion["passed"] else "false",
        )

    if flood_assessment is not None:
        add(
            "source_consensus_adjustment",
            "source_excess_detection_result",
            "exactly_one_qualifying_source"
            if flood_assessment.qualifying_excess_source is not None
            else "no_unique_qualifying_source",
        )
        add(
            "source_consensus_adjustment",
            "adjustment_enabled",
            "true" if flood_assessment.adjustment_enabled else "false",
        )
        add(
            "source_consensus_adjustment",
            "adjustment_applied",
            "true" if flood_assessment.adjustment_applied else "false",
        )
        add(
            "source_consensus_adjustment",
            "discounted_source",
            flood_assessment.qualifying_excess_source or ".",
        )
        add(
            "source_consensus_adjustment",
            "discounted_source_singleton_subloci",
            flood_assessment.discounted_source_singleton_subloci,
        )
        adjusted = flood_assessment.adjusted
        if adjusted is not None:
            for name in (
                "source_count_distribution", "total_subloci",
                "total_source_incidences", "raw_source_count",
                "effective_source_count", "multisource_subloci",
                "singleton_sublocus_fraction", "multisource_sublocus_fraction",
                "deep_consensus_sublocus_fraction",
                "full_consensus_sublocus_fraction",
                "singleton_incidence_fraction", "consensus_incidence_fraction",
                "deep_consensus_incidence_fraction", "corroboration_fraction",
                "mean_sources_per_sublocus",
                "mean_sources_at_multisource_subloci", "backbone_class",
            ):
                value = getattr(adjusted, name)
                if name == "source_count_distribution":
                    value = _serialize_counter({
                        str(key): number for key, number in value.items()
                    })
                add("source_consensus_adjusted", name, value)
            for name, criterion in consensus_classification_criteria(adjusted).items():
                add("source_consensus_adjusted_criteria", f"{name}.value", criterion["value"])
                add("source_consensus_adjusted_criteria", f"{name}.operator", criterion["operator"])
                add("source_consensus_adjusted_criteria", f"{name}.threshold", criterion["threshold"])
                add("source_consensus_adjusted_criteria", f"{name}.margin", criterion["margin"])
                add(
                    "source_consensus_adjusted_criteria",
                    f"{name}.passed",
                    "true" if criterion["passed"] else "false",
                )

    if input_count_metrics:
        first_input_metric = next(iter(input_count_metrics.values()))
        add(
            "input_count_summary",
            "run_median_input_count",
            first_input_metric.run_median_input_count,
        )
    for source, metric in sorted((input_count_metrics or {}).items()):
        add("per_source_input", f"{source}.input_coding_model_count", metric.input_coding_model_count)
        add(
            "per_source_input",
            f"{source}.ratio_to_run_median",
            metric.ratio_to_run_median if metric.ratio_to_run_median is not None else ".",
        )
        add("per_source_input", f"{source}.other_source_median_count", metric.other_source_median_count)
        add("per_source_input", f"{source}.ratio_to_other_source_median", metric.ratio_to_other_source_median)
        add("per_source_input", f"{source}.warning_code", metric.warning)

    add("selected_population", "selected_singleton_model_fraction", selected_singleton_model_fraction)
    add("validation", "count_reconciliation", "pass" if validation.reconciliation_passed else "fail")
    add("validation", "medium_to_high_monotonicity", "pass" if validation.medium_high_monotonicity_passed else "fail")
    add("validation", "high_to_very_high_monotonicity", "pass" if validation.high_very_high_monotonicity_passed else "fail")
    add("validation", "unexpected_restored_medium_to_high", validation.unexpected_restored_medium_to_high)
    add("validation", "unexpected_restored_high_to_very_high", validation.unexpected_restored_high_to_very_high)
    add("validation", "native_restored_medium_to_high", validation.native_restored_medium_to_high)
    add("validation", "native_restored_high_to_very_high", validation.native_restored_high_to_very_high)
    add(
        "validation",
        "monotonic_closure_applied",
        "true" if (
            validation.native_restored_medium_to_high
            or validation.native_restored_high_to_very_high
        ) else "false",
    )
    for decision in (standard_decision, conservative_decision):
        prefix = decision.mode
        add("decision", f"{prefix}.unconstrained_level", decision.unconstrained_level)
        add("decision", f"{prefix}.effective_level", decision.effective_level)
        add("decision", f"{prefix}.confidence", decision.confidence)
        add("decision", f"{prefix}.min_level", decision.user_min_level)
        add("decision", f"{prefix}.max_level", decision.user_max_level)
        add("decision", f"{prefix}.reason_codes", ",".join(decision.reason_codes) or ".")
        add("decision", f"{prefix}.warning_codes", ",".join(decision.warning_codes) or ".")
        add("decision", f"{prefix}.backbone_class_used", decision.backbone_class_used)
        add(
            "decision",
            f"{prefix}.used_flood_adjustment",
            "true" if decision.used_flood_adjustment else "false",
        )
        add(
            "decision",
            f"{prefix}.explanation",
            explain_auto_pruning_decision(
                decision,
                consensus=consensus,
                flood_assessment=flood_assessment,
            ),
        )

    if scenarios is not None and "medium" in scenarios and "high" in scenarios:
        medium = scenarios["medium"]
        high = scenarios["high"]
        additional = high.removed_ids - medium.removed_ids
        by_source = Counter(
            high.source_by_id.get(model_id, ".") for model_id in additional
        )
        add("alternative_pruning", "high_vs_medium_additional_models_removed", len(additional))
        add(
            "alternative_pruning",
            "high_vs_medium_fraction_starting_selected_models",
            _safe_fraction(len(additional), high.starting_models),
        )
        add(
            "alternative_pruning",
            "high_vs_medium_additional_removals_by_source",
            _serialize_counter(by_source),
        )
    return rows


def _failed_strong_criterion_descriptions(
    metrics: SourceConsensusMetrics,
) -> list[str]:
    descriptions = {
        "strong.minimum_source_count": "source count",
        "strong.minimum_multisource_subloci": "multisource-sublocus count",
        "strong.maximum_singleton_incidence_fraction": "singleton incidence",
        "strong.minimum_corroboration_fraction": "corroboration",
        "strong.minimum_deep_consensus_incidence_fraction": "deep-consensus incidence",
        "strong.minimum_deep_consensus_sublocus_fraction": "deep-consensus sublocus fraction",
    }
    failed: list[str] = []
    for name, item in consensus_classification_criteria(metrics).items():
        if not name.startswith("strong.") or item["passed"]:
            continue
        value = item["value"]
        threshold = item["threshold"]
        operator = item["operator"]
        if isinstance(value, float) or isinstance(threshold, float):
            failed.append(
                f"{descriptions[name]} {float(value):.3f} does not satisfy "
                f"{operator} {float(threshold):.3f}"
            )
        else:
            failed.append(
                f"{descriptions[name]} {value} does not satisfy {operator} {threshold}"
            )
    return failed


def flood_adjusted_consensus_log_lines(
    assessment: FloodAdjustedConsensusAssessment,
    *,
    requested_mode: str,
    standard_decision: AutoDecision | None,
    conservative_decision: AutoDecision | None,
) -> list[str]:
    """Format an applied flood adjustment for the normal run log.

    The values and thresholds come from the same metric and criterion helpers
    used by classification.  This function is reporting-only.
    """
    adjusted = assessment.adjusted
    if adjusted is None:
        return []
    raw = assessment.unadjusted
    criteria = consensus_classification_criteria(raw)
    metric_rows = (
        (
            "singleton incidence",
            "singleton_incidence_fraction",
            "strong.maximum_singleton_incidence_fraction",
        ),
        (
            "corroboration",
            "corroboration_fraction",
            "strong.minimum_corroboration_fraction",
        ),
        (
            "deep-consensus incidence",
            "deep_consensus_incidence_fraction",
            "strong.minimum_deep_consensus_incidence_fraction",
        ),
        (
            "deep-consensus sublocus fraction",
            "deep_consensus_sublocus_fraction",
            "strong.minimum_deep_consensus_sublocus_fraction",
        ),
    )
    lines = [
        "Flood-adjusted consensus detail:",
        (
            f"  excessive_source={assessment.qualifying_excess_source or '.'}; "
            f"discounted_source_singleton_subloci="
            f"{assessment.discounted_source_singleton_subloci:,}; "
            f"source_count={raw.effective_source_count}"
        ),
        "  Metric                              Raw  Adjusted  Strong threshold",
    ]
    for label, attribute, criterion_name in metric_rows:
        criterion = criteria[criterion_name]
        lines.append(
            f"  {label:<31} "
            f"{getattr(raw, attribute):>7.4f}  "
            f"{getattr(adjusted, attribute):>8.4f}  "
            f"{criterion['operator']} {float(criterion['threshold']):.4f}"
        )
    lines.append(
        f"  Backbone class: {raw.backbone_class} -> {adjusted.backbone_class}"
    )
    if (
        requested_mode == "auto"
        and standard_decision is not None
        and standard_decision.used_flood_adjustment
    ):
        lines.append(
            "  Decision basis: standard auto uses the adjusted backbone on the "
            "no-adequate-evidence pathway; "
            f"class_changed={str(raw.backbone_class != adjusted.backbone_class).lower()}; "
            f"recommended={standard_decision.unconstrained_level}; "
            f"effective={standard_decision.effective_level}."
        )
    elif requested_mode == "auto" and standard_decision is not None:
        lines.append(
            "  Decision basis: the adjusted assessment did not determine this "
            "standard-auto recommendation; "
            f"recommended={standard_decision.unconstrained_level}; "
            f"effective={standard_decision.effective_level}."
        )
    elif requested_mode == "auto-conservative" and conservative_decision is not None:
        lines.append(
            "  Decision basis: auto-conservative uses the unadjusted backbone; "
            "the adjusted class cannot promote the recommendation; "
            f"recommended={conservative_decision.unconstrained_level}; "
            f"effective={conservative_decision.effective_level}."
        )
    else:
        lines.append(
            f"  Decision basis: requested level {requested_mode} is fixed; the "
            "adjusted assessment is reporting-only."
        )
    lines.append(
        "  Adjustment scope: consensus-backbone statistics only; no subloci or "
        "prediction models were removed."
    )
    return lines


def _with_effective_level_constraint(
    decision: AutoDecision,
    explanation: str,
) -> str:
    """Append the user-bound effect without changing the decision itself."""
    if decision.effective_level == decision.unconstrained_level:
        return explanation
    if (
        LEVEL_ORDER[decision.effective_level]
        < LEVEL_ORDER[decision.unconstrained_level]
    ):
        return (
            f"{explanation} The user-specified maximum pruning level is "
            f"{decision.user_max_level}, so the effective level is "
            f"{decision.effective_level}."
        )
    return (
        f"{explanation} The user-specified minimum pruning level is "
        f"{decision.user_min_level}, so the effective level is "
        f"{decision.effective_level}."
    )


def explain_auto_pruning_decision(
    decision: AutoDecision,
    *,
    consensus: SourceConsensusMetrics,
    flood_assessment: FloodAdjustedConsensusAssessment | None = None,
) -> str:
    """Build one concise human-readable explanation of an automatic decision."""
    metrics = consensus
    if (
        decision.used_flood_adjustment
        and flood_assessment is not None
        and flood_assessment.adjusted is not None
    ):
        metrics = flood_assessment.adjusted
    if "SCENARIO_VALIDATION_FAILED_FALLBACK_MEDIUM" in decision.reason_codes:
        return _with_effective_level_constraint(
            decision,
            "Automatic pruning recommended medium because scenario validation "
            "failed; the conservative safety fallback was applied.",
        )
    if any(code.startswith("ADEQUATE_EVIDENCE_") for code in decision.reason_codes):
        triggers: list[str] = []
        if "ADEQUATE_EVIDENCE_SINGLETON_BURDEN_HIGH_TRIGGER" in decision.reason_codes:
            triggers.append("selected-singleton burden")
        if "ADEQUATE_EVIDENCE_PROSPECTIVE_HIGH_INCREMENT_TRIGGER" in decision.reason_codes:
            triggers.append("the prospective high-pruning increment")
        if "VERY_HIGH_INCREMENT_WITHIN_SAFETY_LIMIT" in decision.reason_codes:
            triggers.append("the very-high increment remaining within its safety limit")
        trigger_text = ", ".join(triggers) or "the evidence-aware high trigger"
        return _with_effective_level_constraint(
            decision,
            f"Automatic pruning recommended {decision.unconstrained_level} on the "
            f"adequate-evidence pathway because of {trigger_text}; flood-adjusted "
            "backbone classification does not alter this pathway.",
        )
    if "INTERMEDIATE_BACKBONE_SINGLE_STRONG_SOURCE_EXCESS" in decision.reason_codes:
        source = (
            consensus.strong_source_excess_sources[0]
            if len(consensus.strong_source_excess_sources) == 1
            else "."
        )
        return _with_effective_level_constraint(
            decision,
            f"Automatic pruning recommended {decision.unconstrained_level}. The backbone "
            f"is intermediate, source {source} has a strong source-singleton "
            "excess, and the prospective high-versus-medium increment is "
            "concentrated in that source.",
        )
    if (
        decision.mode == "auto-conservative"
        and flood_assessment is not None
        and flood_assessment.adjusted is not None
    ):
        adjusted_class = flood_assessment.adjusted.backbone_class
        return _with_effective_level_constraint(
            decision,
            f"Flood-adjusted consensus is {adjusted_class}, but auto-conservative "
            f"uses the unadjusted {consensus.backbone_class} backbone. "
            f"{decision.unconstrained_level.capitalize()} pruning was recommended "
            "from the unadjusted backbone; the statistical adjustment did not "
            "remove any prediction models.",
        )
    if decision.used_flood_adjustment and metrics.backbone_class == "strong":
        source = flood_assessment.qualifying_excess_source if flood_assessment else "."
        count = flood_assessment.discounted_source_singleton_subloci if flood_assessment else 0
        return _with_effective_level_constraint(
            decision,
            f"Automatic pruning recommended {decision.unconstrained_level}. Source {source} "
            f"qualified as strongly excessive; {count:,} source-singleton subloci "
            "were discounted only for backbone assessment, and the adjusted "
            "consensus backbone is strong. No prediction models were removed by "
            "this statistical adjustment.",
        )
    if metrics.backbone_class == "strong":
        threshold = strong_deep_consensus_incidence_threshold(
            metrics.effective_source_count
        )
        return _with_effective_level_constraint(
            decision,
            f"Automatic pruning recommended {decision.unconstrained_level}. The "
            f"{metrics.effective_source_count}-source consensus backbone is strong; "
            f"deep-consensus incidence is {metrics.deep_consensus_incidence_fraction:.3f} "
            f"against the {threshold:.2f} threshold.",
        )
    failed = _failed_strong_criterion_descriptions(metrics)
    details = "; ".join(failed) if failed else "strong consensus was not established"
    excess_text = (
        " Exactly one source satisfies the strong source-excess detector."
        if len(consensus.strong_source_excess_sources) == 1
        else " No single source satisfies the strong source-excess detector."
    )
    return _with_effective_level_constraint(
        decision,
        f"Automatic pruning recommended {decision.unconstrained_level}. Strong consensus "
        f"was not established: {details}.{excess_text} High pruning can be "
        "evaluated manually.",
    )


def write_auto_diagnostics(path: Path, rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUTO_DIAGNOSTIC_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def annotate_result_with_scenario_audit(
    result: dict[str, object],
    scenarios: Mapping[str, ScenarioSummary],
) -> None:
    """Add compact prospective fields to the chosen retained/removed rows."""
    all_rows = list(result.get("retained_rows", [])) + list(result.get("removed_rows", []))  # type: ignore[arg-type]
    for row in all_rows:
        model_id = str(row.get("output_transcript_id", "."))
        for level, suffix in (("medium", "medium"), ("high", "high"), ("very-high", "very_high")):
            scenario = scenarios[level]
            row[f"removed_by_{suffix}"] = "true" if model_id in scenario.removed_ids else "false"
            row[f"native_removed_by_{suffix}"] = "true" if model_id in scenario.native_removed_ids else "false"
            row[f"inherited_removal_by_{suffix}"] = "true" if model_id in scenario.inherited_removed_ids else "false"
            row[f"inherited_from_level_{suffix}"] = scenario.inherited_from_level_by_id.get(model_id, ".")
            row[f"inherited_original_level_{suffix}"] = scenario.inherited_original_level_by_id.get(model_id, ".")
            row[f"inherited_original_reason_{suffix}"] = scenario.inherited_original_reason_by_id.get(model_id, ".")
            row[f"primary_reason_{suffix}"] = scenario.primary_reason_by_id.get(model_id, ".")
        protein = any(model_id in scenario.protein_protected_ids for scenario in scenarios.values())
        transcript = any(model_id in scenario.transcript_protected_ids for scenario in scenarios.values())
        overlap = any(model_id in scenario.overlap_protected_ids for scenario in scenarios.values())
        row["protected_by_protein_in_any_scenario"] = "true" if protein else "false"
        row["protected_by_transcript_in_any_scenario"] = "true" if transcript else "false"
        row["protected_by_overlap_policy_in_any_scenario"] = "true" if overlap else "false"
