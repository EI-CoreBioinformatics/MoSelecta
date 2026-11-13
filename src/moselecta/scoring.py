"""Generic scoring helpers and aggregate model-selection-score calculations."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Callable, Optional

from .constants import LOW_SOURCE_CONFLICT_AGREEMENT_FACTOR, LOW_SOURCE_CONFLICT_PAIRWISE_THRESHOLD


def saturating_corroboration_score(
    corroborating_other_source_count: int,
    run_possible_other_source_count: int,
    alpha: float,
) -> float:
    """Return the bounded run-depth-normalised corroboration score."""

    k = int(corroborating_other_source_count)
    n = int(run_possible_other_source_count)
    if n < 0 or k < 0 or k > n:
        raise ValueError(f"corroboration counts must satisfy 0 <= k <= N; got k={k}, N={n}")
    if not math.isfinite(alpha) or alpha <= 0:
        raise ValueError(f"saturating alpha must be finite and > 0; got {alpha}")
    if n == 0 or k == 0:
        return 0.0
    return (k / (k + alpha)) / (n / (n + alpha))

def fmt_optional_float(value: Optional[float], digits: int = 6) -> str:
    if value is None:
        return "."
    if math.isnan(value):
        return "."
    return f"{value:.{digits}f}"


def fmt_sources(sources: set[str]) -> str:
    return ",".join(sorted(sources)) if sources else "."


def fmt_model_ids(model_ids: set[str]) -> str:
    return ",".join(sorted(model_ids)) if model_ids else "."


def fraction_or_none(numerator: int, denominator: int) -> Optional[float]:
    if denominator == 0:
        return None
    return numerator / denominator


def source_scaled_fraction(numerator: int, denominator: int) -> float:
    """Scale support by possible other-source count.

    If there are no possible other sources in a source-singleton sublocus, this
    returns 0.0 rather than NA because there is no independent source support.
    """
    if denominator <= 0:
        return 0.0
    return numerator / denominator


def weighted_score(components: list[tuple[str, float, Optional[float]]]) -> tuple[Optional[float], str]:
    available = [(name, weight, value) for name, weight, value in components if value is not None]
    if not available:
        return None, "."
    total_weight = sum(weight for _name, weight, _value in available)
    if total_weight <= 0:
        return None, "."
    score = sum(weight * value for _name, weight, value in available) / total_weight
    component_text = ",".join(f"{name}:{value:.6f}*{weight:.3f}" for name, weight, value in available)
    return score, component_text


def value_as_float(value: object, default: float = 0.0) -> float:
    if value is None or value == ".":
        return default
    try:
        f = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(f):
        return default
    return f


PRIMARY_MODEL_SELECTION_TIE_BREAKERS = (
    ("model_completeness_score", "desc"),
    ("intrinsic_penalty_score", "asc"),
    ("consensus_similarity_score", "desc"),
    ("structural_support_score", "desc"),
    ("source_prior_score", "desc"),
    ("cds_length", "desc"),
    ("output_transcript_id", "asc"),
)


def primary_model_selection_order_key(
    row: dict[str, object],
    *,
    score_field: str = "model_selection_score",
) -> tuple:
    """Return the authoritative deterministic primary-selection total order."""
    return (
        -value_as_float(row.get(score_field)),
        -value_as_float(row.get("model_completeness_score")),
        value_as_float(row.get("intrinsic_penalty_score")),
        -value_as_float(row.get("consensus_similarity_score")),
        -value_as_float(row.get("structural_support_score")),
        -value_as_float(row.get("source_prior_score")),
        -int(value_as_float(row.get("cds_length"))),
        str(row.get("output_transcript_id", "")),
    )


def primary_model_selection_tie_break_value(row: dict[str, object], field: str) -> object:
    """Return a field value normalised exactly as the primary ordering key."""
    if field == "output_transcript_id":
        return str(row.get(field, ""))
    if field == "cds_length":
        return int(value_as_float(row.get(field)))
    return value_as_float(row.get(field))


TRANSCRIPT_CDNA_CONTEXT_MODERATE_PRUNING_FLAGS = frozenset({
    "cds_chain_contained_in_longer_cdna",
    "many_extra_5prime_transcript_exons",
    "many_extra_3prime_transcript_exons",
    "single_exon_cds_in_multiexon_cdna_context",
})


def transcript_cdna_context_concern_fields(
    transcript_row: Optional[dict[str, object]],
    *,
    neutral_fraction: float = 0.30,
    severe_fraction: float = 0.15,
) -> dict[str, object]:
    """Classify CDS/cDNA context while treating missing context as neutral.

    The continuous concern score is used for the modest ranking adjustment. A
    pruning concern is more conservative: a fraction below ``severe_fraction``
    is sufficient, whereas a fraction between the severe and neutral thresholds
    must also have a structural context flag.
    """
    if not (0.0 <= severe_fraction <= neutral_fraction <= 1.0) or neutral_fraction <= 0.0:
        raise ValueError("cDNA-context fractions must satisfy 0 <= severe <= neutral <= 1 and neutral > 0")
    transcript_row = transcript_row or {}
    status = str(transcript_row.get("transcript_cdna_context_status", "."))
    raw_fraction = transcript_row.get("transcript_cdna_context_cds_to_min_cdna_fraction", ".")
    try:
        fraction = float(raw_fraction)
    except (TypeError, ValueError):
        fraction = math.nan
    evaluable = status == "calculated" and math.isfinite(fraction) and 0.0 < fraction <= 1.0 + 1e-12
    if not evaluable:
        return {
            "transcript_cdna_context_evaluable": "false",
            "transcript_cdna_context_concern_score": 0.0,
            "transcript_cdna_context_concern_class": "not_evaluable",
            "transcript_cdna_context_pruning_concern": "false",
            "transcript_cdna_context_pruning_concern_class": "not_evaluable",
            "transcript_cdna_context_pruning_supporting_flags": ".",
            "transcript_cdna_context_neutral_fraction": neutral_fraction,
            "transcript_cdna_context_severe_fraction": severe_fraction,
        }

    concern_score = max(0.0, min(1.0, (neutral_fraction - fraction) / neutral_fraction))
    raw_flags = str(transcript_row.get("transcript_cdna_context_flags", "."))
    flags = {flag for flag in raw_flags.split(",") if flag and flag != "."}
    supporting_flags = sorted(flags & TRANSCRIPT_CDNA_CONTEXT_MODERATE_PRUNING_FLAGS)

    pruning_concern = False
    pruning_class = "evaluated_no_concern"
    if fraction < severe_fraction:
        concern_class = "evaluated_severe_concern"
        pruning_class = "severe_low_fraction"
        pruning_concern = True
    elif fraction < neutral_fraction:
        if supporting_flags:
            concern_class = "evaluated_moderate_concern"
            pruning_class = "moderate_low_fraction_with_structural_context"
            pruning_concern = True
        else:
            concern_class = "evaluated_low_fraction"
            pruning_class = "low_fraction_without_structural_context"
    else:
        concern_class = "evaluated_no_concern"

    return {
        "transcript_cdna_context_evaluable": "true",
        "transcript_cdna_context_concern_score": concern_score,
        "transcript_cdna_context_concern_class": concern_class,
        "transcript_cdna_context_pruning_concern": "true" if pruning_concern else "false",
        "transcript_cdna_context_pruning_concern_class": pruning_class,
        "transcript_cdna_context_pruning_supporting_flags": ",".join(supporting_flags) if supporting_flags else ".",
        "transcript_cdna_context_neutral_fraction": neutral_fraction,
        "transcript_cdna_context_severe_fraction": severe_fraction,
    }


def annotate_transcript_cdna_context_concern_rows(
    rows: list[dict[str, object]],
    *,
    neutral_fraction: float = 0.30,
    severe_fraction: float = 0.15,
) -> None:
    """Add cDNA-context concern fields to transcript evidence rows."""
    for row in rows:
        row.update(
            transcript_cdna_context_concern_fields(
                row,
                neutral_fraction=neutral_fraction,
                severe_fraction=severe_fraction,
            )
        )


def values_equal_for_tie(a: object, b: object, epsilon: float) -> bool:
    try:
        fa = float(a); fb = float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)
    return abs(fa - fb) <= epsilon


def model_selection_score_positive_scale_max(selection_weights: dict[str, float], *, include_protein_evidence: bool = True, include_transcript_evidence: bool = True) -> float:
    """Return the theoretical positive maximum for model_selection_score_scaled.

    The scaled score is intentionally a zero-based positive transform:
    raw scores <= 0 map to 0, and positive raw scores are divided by the
    theoretical maximum possible positive contribution under the configured
    selection weights. Negative penalty bounds are not needed for this mapping
    and raw-score ranking is unchanged.
    """
    return (
        selection_weights["structural_support"]
        + selection_weights["consensus_similarity"]
        + selection_weights["pairwise_similarity"]
        + selection_weights["completeness"]
        + selection_weights["source_prior"]
        + selection_weights["terminal_extension"]
        + selection_weights["exact_structure_support"]
        + selection_weights.get("structural_corroboration", 0.0)
        + (selection_weights.get("protein_evidence", 0.0) if include_protein_evidence else 0.0)
        + (selection_weights.get("transcript_evidence", 0.0) if include_transcript_evidence else 0.0)
    )


def scale_model_selection_score(raw_score: float, selection_weights: dict[str, float], *, include_protein_evidence: bool = True, include_transcript_evidence: bool = True) -> float:
    """Map the raw model selection score onto a review-friendly 0-1 positive scale.

    Negative and zero raw scores are mapped to 0. Positive scores are scaled by
    the theoretical maximum positive contribution under the configured weights
    and clamped at 1. This makes the scaled value easier to interpret in reports
    and source summaries, while leaving raw-score ranking unchanged.
    """
    theoretical_positive_max = model_selection_score_positive_scale_max(selection_weights, include_protein_evidence=include_protein_evidence, include_transcript_evidence=include_transcript_evidence)
    if theoretical_positive_max <= 0 or raw_score <= 0:
        return 0.0
    return min(raw_score / theoretical_positive_max, 1.0)


def _score_model_row_components(
    row: dict[str, object],
    source_prior_score: float,
    selection_weights: dict[str, float],
    show_zero_weight_components: bool,
    low_source_conflict_dampen_agreement_scores: bool,
) -> tuple[float, str, bool, float, str]:
    """Return model-only score/components and low-source conflict annotations."""
    sublocus_source_count = int(value_as_float(row.get("sublocus_source_count")))
    sublocus_model_count = int(value_as_float(row.get("sublocus_model_count")))
    pairwise_other_source_count = int(value_as_float(row.get("pairwise_other_source_count")))
    exact_structure_support_score = value_as_float(row.get("exact_structure_support_score"))
    pairwise_similarity_score = value_as_float(row.get("pairwise_similarity_score"))
    is_low_source_conflict = (
        sublocus_source_count == 2
        and sublocus_model_count == 2
        and pairwise_other_source_count > 0
        and exact_structure_support_score <= 0.0
        and pairwise_similarity_score < LOW_SOURCE_CONFLICT_PAIRWISE_THRESHOLD
    )
    agreement_factor = (
        LOW_SOURCE_CONFLICT_AGREEMENT_FACTOR
        if (low_source_conflict_dampen_agreement_scores and is_low_source_conflict)
        else 1.0
    )
    if is_low_source_conflict:
        reason = (
            f"sources_eq_2,models_eq_2,no_exact_structure_support,"
            f"pairwise_similarity_lt_{LOW_SOURCE_CONFLICT_PAIRWISE_THRESHOLD:.2f}"
        )
        if not low_source_conflict_dampen_agreement_scores:
            reason += ",dampening_disabled"
    else:
        reason = "."

    components = [
        ("structural_support_score", "+", selection_weights["structural_support"], value_as_float(row.get("structural_support_score")), agreement_factor),
        ("consensus_similarity_score", "+", selection_weights["consensus_similarity"], value_as_float(row.get("consensus_similarity_score")), agreement_factor),
        ("pairwise_similarity_score", "+", selection_weights["pairwise_similarity"], pairwise_similarity_score, 1.0),
        ("model_completeness_score", "+", selection_weights["completeness"], value_as_float(row.get("model_completeness_score")), 1.0),
        ("source_prior_score", "+", selection_weights["source_prior"], source_prior_score, 1.0),
        ("terminal_extension_score", "+", selection_weights["terminal_extension"], value_as_float(row.get("terminal_extension_score")), 1.0),
        ("exact_structure_support_score", "+", selection_weights["exact_structure_support"], exact_structure_support_score, 1.0),
        ("structural_corroboration_score", "+", selection_weights.get("structural_corroboration", 0.0), value_as_float(row.get("structural_corroboration_score")), 1.0),
        ("unsupported_short_cds_exon_penalty_score", "-", selection_weights.get("unsupported_short_cds_exon_penalty", 0.0), value_as_float(row.get("unsupported_short_cds_exon_penalty_score")), 1.0),
        ("structural_coherence_raw_score", "-", selection_weights.get("structural_coherence_penalty", 0.0), value_as_float(row.get("structural_coherence_raw_score")), 1.0),
        ("sequence_integrity_penalty_score", "-", selection_weights.get("sequence_integrity_penalty", 0.0), value_as_float(row.get("sequence_integrity_penalty_score")), 1.0),
        ("intrinsic_penalty_score", "-", selection_weights["intrinsic_penalty"], value_as_float(row.get("intrinsic_penalty_score")), 1.0),
        ("partial_structure_penalty", "-", selection_weights["partial_structure_penalty"], value_as_float(row.get("partial_structure_penalty")), 1.0),
        ("internal_extra_fraction", "-", selection_weights["internal_extra_penalty"], value_as_float(row.get("internal_extra_fraction")), 1.0),
        ("terminal_extra_fraction", "-", selection_weights["terminal_extra_penalty"], value_as_float(row.get("terminal_extra_fraction")), 1.0),
    ]
    score = 0.0
    component_text: list[str] = []
    for name, sign, weight, value, factor in components:
        contribution = weight * value * factor
        if sign == "+":
            score += contribution
        else:
            score -= contribution
        if show_zero_weight_components or weight != 0:
            if factor != 1.0:
                component_text.append(f"{sign}{name}:{value:.6f}*{weight:.3f}*low_source_factor:{factor:.3f}")
            else:
                component_text.append(f"{sign}{name}:{value:.6f}*{weight:.3f}")
    return score, ",".join(component_text), is_low_source_conflict, agreement_factor, reason


def _assign_score_ranks(
    rows_by_sublocus: dict[str, list[dict[str, object]]],
    *,
    score_field: str,
    score_rank_field: str,
    order_rank_field: str,
    tie_group_field: str,
    tie_break_reason_field: str,
    tie_epsilon: float,
    order_key: Optional[Callable[[dict[str, object]], tuple]] = None,
    tie_breakers_override: Optional[tuple[tuple[str, str], ...]] = None,
    tie_break_value: Optional[Callable[[dict[str, object], str], object]] = None,
    tie_break_reason_epsilon: Optional[float] = None,
) -> None:
    def legacy_sort_key(row: dict[str, object]) -> tuple:
        return (
            -value_as_float(row.get(score_field)),
            -value_as_float(row.get("model_completeness_score")),
            value_as_float(row.get("intrinsic_penalty_score")),
            value_as_float(row.get("unsupported_short_cds_exon_penalty_score")),
            value_as_float(row.get("sequence_integrity_penalty_contribution")),
            -value_as_float(row.get("consensus_similarity_score")),
            -value_as_float(row.get("structural_support_score")),
            -value_as_float(row.get("exact_structure_support_score")),
            -value_as_float(row.get("source_prior_score")),
            -int(value_as_float(row.get("cds_length"))),
            str(row.get("output_transcript_id", "")),
        )

    legacy_tie_breakers = (
        ("model_completeness_score", "desc"),
        ("intrinsic_penalty_score", "asc"),
        ("unsupported_short_cds_exon_penalty_score", "asc"),
        ("sequence_integrity_penalty_contribution", "asc"),
        ("consensus_similarity_score", "desc"),
        ("structural_support_score", "desc"),
        ("exact_structure_support_score", "desc"),
        ("source_prior_score", "desc"),
        ("cds_length", "desc"),
        ("output_transcript_id", "asc"),
    )
    effective_order_key = order_key or legacy_sort_key
    tie_breakers = tie_breakers_override or legacy_tie_breakers
    reason_epsilon = tie_epsilon if tie_break_reason_epsilon is None else tie_break_reason_epsilon

    def reason_value(row: dict[str, object], field: str) -> object:
        if tie_break_value is not None:
            return tie_break_value(row, field)
        return row.get(field, ".")

    for _sublocus_id, rows in rows_by_sublocus.items():
        ordered = sorted(rows, key=effective_order_key)
        for idx, row in enumerate(ordered, start=1):
            row[order_rank_field] = idx

        rank = 1
        idx = 0
        tie_group_index = 1
        while idx < len(ordered):
            group = [ordered[idx]]
            score = value_as_float(ordered[idx].get(score_field))
            j = idx + 1
            while j < len(ordered) and abs(value_as_float(ordered[j].get(score_field)) - score) <= tie_epsilon:
                group.append(ordered[j])
                j += 1
            for row in group:
                row[score_rank_field] = rank
            if len(group) > 1:
                tie_group = f"tie{tie_group_index:03d}"
                reason = "output_transcript_id"
                for name, _direction in tie_breakers:
                    first = reason_value(group[0], name)
                    if any(
                        not values_equal_for_tie(first, reason_value(other, name), reason_epsilon)
                        for other in group[1:]
                    ):
                        reason = name
                        break
                for row in group:
                    row[tie_group_field] = tie_group
                    row[tie_break_reason_field] = reason
                tie_group_index += 1
            rank += len(group)
            idx = j


def _assign_model_selection_score_ranks(
    rows_by_sublocus: dict[str, list[dict[str, object]]],
    *,
    tie_epsilon: float,
) -> None:
    """Assign score groups and selector-consistent total-order ranks."""
    _assign_score_ranks(
        rows_by_sublocus,
        score_field="model_selection_score",
        score_rank_field="model_selection_score_rank",
        order_rank_field="model_selection_order_rank",
        tie_group_field="model_selection_tie_group",
        tie_break_reason_field="model_selection_tie_break_reason",
        tie_epsilon=tie_epsilon,
        order_key=lambda row: primary_model_selection_order_key(row, score_field="model_selection_score"),
        tie_breakers_override=(("model_selection_score", "desc"),) + PRIMARY_MODEL_SELECTION_TIE_BREAKERS,
        tie_break_value=primary_model_selection_tie_break_value,
        tie_break_reason_epsilon=0.0,
    )


def apply_model_selection_scores(
    model_rows: list[dict[str, object]],
    source_prior_scores: dict[str, float],
    selection_weights: dict[str, float],
    tie_epsilon: float,
    show_zero_weight_components: bool = False,
    low_source_conflict_dampen_agreement_scores: bool = True,
    protein_evidence_rows_by_transcript: Optional[dict[str, dict[str, object]]] = None,
    protein_evidence_active: bool = False,
    transcript_evidence_rows_by_transcript: Optional[dict[str, dict[str, object]]] = None,
    transcript_evidence_active: bool = False,
    transcript_cdna_context_neutral_fraction: float = 0.30,
    transcript_cdna_context_severe_fraction: float = 0.15,
    transcript_cdna_context_penalty_weight: float = 0.15,
    structural_corroboration_alpha: float = 16.0,
) -> None:
    """Add model-only and final model selection score/rank fields.

    ``model_structural_support_score*`` captures the model-only ranking before
    external evidence integration. ``model_selection_score*`` is the final
    score used by selection and includes protein and/or transcript evidence when
    the corresponding evidence class is active with a positive weight.
    """
    protein_evidence_rows_by_transcript = protein_evidence_rows_by_transcript or {}
    transcript_evidence_rows_by_transcript = transcript_evidence_rows_by_transcript or {}
    protein_weight = selection_weights.get("protein_evidence", 0.0) if protein_evidence_active else 0.0
    transcript_weight = selection_weights.get("transcript_evidence", 0.0) if transcript_evidence_active else 0.0
    corroboration_weight = selection_weights.get("structural_corroboration", 0.0)

    for row in model_rows:
        source_prior_score = float(source_prior_scores.get(str(row.get("source_name", "")), 0.0))
        row["source_prior_score"] = source_prior_score
        corroboration_k = int(value_as_float(row.get("corroborating_other_source_count")))
        corroboration_n = int(value_as_float(row.get("run_possible_other_source_count")))
        corroboration_score = saturating_corroboration_score(
            corroboration_k,
            corroboration_n,
            structural_corroboration_alpha,
        )
        row["structural_corroboration_score"] = corroboration_score
        row["structural_corroboration_alpha"] = structural_corroboration_alpha
        row["structural_corroboration_weight"] = corroboration_weight
        row["structural_corroboration_contribution"] = corroboration_weight * corroboration_score
        row["structural_corroboration_in_selection"] = (
            "true" if corroboration_weight > 0 else "false"
        )

        structural_score, structural_components, is_low_source_conflict, agreement_factor, reason = _score_model_row_components(
            row,
            source_prior_score,
            selection_weights,
            show_zero_weight_components,
            low_source_conflict_dampen_agreement_scores,
        )
        row["low_source_conflict_flag"] = "true" if is_low_source_conflict else "false"
        row["low_source_conflict_agreement_factor"] = agreement_factor
        row["low_source_conflict_reason"] = reason

        sequence_penalty_weight = selection_weights.get("sequence_integrity_penalty", 0.0)
        sequence_penalty_score = value_as_float(row.get("sequence_integrity_penalty_score"))
        sequence_penalty_contribution = sequence_penalty_weight * sequence_penalty_score
        score_before_sequence = structural_score + sequence_penalty_contribution
        row["model_structural_support_score_before_sequence_integrity"] = score_before_sequence
        row["model_structural_support_score_before_sequence_integrity_scaled"] = scale_model_selection_score(
            score_before_sequence,
            selection_weights,
            include_protein_evidence=False,
            include_transcript_evidence=False,
        )
        row["sequence_integrity_penalty_weight"] = sequence_penalty_weight
        row["sequence_integrity_penalty_contribution"] = sequence_penalty_contribution
        row["sequence_integrity_penalty_in_selection"] = (
            "true" if sequence_penalty_weight > 0 and sequence_penalty_score > 0 else "false"
        )
        row["model_structural_support_score"] = structural_score
        row["model_structural_support_score_scaled"] = scale_model_selection_score(
            structural_score,
            selection_weights,
            include_protein_evidence=False,
            include_transcript_evidence=False,
        )
        row["model_structural_support_score_components"] = structural_components
        unsupported_short_penalty_weight = selection_weights.get("unsupported_short_cds_exon_penalty", 0.0)
        unsupported_short_penalty_score = value_as_float(row.get("unsupported_short_cds_exon_penalty_score"))
        row["unsupported_short_cds_exon_penalty_weight"] = unsupported_short_penalty_weight
        row["unsupported_short_cds_exon_penalty_in_selection"] = (
            "true" if unsupported_short_penalty_weight > 0 and unsupported_short_penalty_score > 0 else "false"
        )
        structural_coherence_weight = selection_weights.get("structural_coherence_penalty", 0.0)
        structural_coherence_raw = value_as_float(row.get("structural_coherence_raw_score"))
        row["structural_coherence_penalty_weight"] = structural_coherence_weight
        row["structural_coherence_penalty_score"] = structural_coherence_weight * structural_coherence_raw
        row["structural_coherence_penalty_in_selection"] = (
            "true" if structural_coherence_weight > 0 and structural_coherence_raw > 0 else "false"
        )
        row["model_structural_support_score_rank"] = "."
        row["model_structural_support_order_rank"] = "."
        row["model_structural_support_tie_group"] = "."
        row["model_structural_support_tie_break_reason"] = "."

        protein_row = protein_evidence_rows_by_transcript.get(str(row.get("output_transcript_id", "")))
        protein_score = value_as_float(protein_row.get("protein_evidence_score") if protein_row else 0.0)
        protein_score_scaled = value_as_float(protein_row.get("protein_evidence_score_scaled") if protein_row else protein_score)
        row["protein_evidence_score"] = protein_score
        row["protein_evidence_score_scaled"] = protein_score_scaled
        row["protein_evidence_score_components"] = protein_row.get("protein_evidence_score_components", ".") if protein_row else "."
        row["protein_evidence_status"] = protein_row.get("protein_evidence_status", "protein_evidence_not_configured") if protein_row else "protein_evidence_not_configured"
        row["protein_evidence_flags"] = protein_row.get("protein_evidence_flags", ".") if protein_row else "."
        row["protein_evidence_score_weight"] = protein_weight
        row["protein_evidence_score_in_selection"] = "true" if protein_weight > 0 else "false"

        transcript_row = transcript_evidence_rows_by_transcript.get(str(row.get("output_transcript_id", "")))
        transcript_score = value_as_float(transcript_row.get("transcript_evidence_score") if transcript_row else 0.0)
        transcript_score_scaled = value_as_float(transcript_row.get("transcript_evidence_score_scaled") if transcript_row else transcript_score)
        context_fields = transcript_cdna_context_concern_fields(
            transcript_row,
            neutral_fraction=transcript_cdna_context_neutral_fraction,
            severe_fraction=transcript_cdna_context_severe_fraction,
        )
        requested_context_penalty = transcript_cdna_context_penalty_weight * value_as_float(
            context_fields.get("transcript_cdna_context_concern_score"), 0.0
        )
        transcript_score_context_adjusted = max(0.0, transcript_score - requested_context_penalty)
        transcript_score_context_adjusted_scaled = max(0.0, transcript_score_scaled - requested_context_penalty)
        context_penalty_applied = transcript_score - transcript_score_context_adjusted
        row["transcript_evidence_score"] = transcript_score
        row["transcript_evidence_score_scaled"] = transcript_score_scaled
        row["transcript_evidence_score_context_adjusted"] = transcript_score_context_adjusted
        row["transcript_evidence_score_context_adjusted_scaled"] = transcript_score_context_adjusted_scaled
        row["transcript_evidence_score_components"] = transcript_row.get("transcript_evidence_score_components", ".") if transcript_row else "."
        row["transcript_evidence_status"] = transcript_row.get("transcript_evidence_status", "transcript_evidence_not_configured") if transcript_row else "transcript_evidence_not_configured"
        row["transcript_evidence_flags"] = transcript_row.get("transcript_evidence_flags", ".") if transcript_row else "."
        for field in (
            "transcript_cdna_context_status",
            "transcript_cdna_context_cds_to_min_cdna_fraction",
            "transcript_cdna_context_flags",
        ):
            row[field] = transcript_row.get(field, ".") if transcript_row else "."
        row.update(context_fields)
        row["transcript_cdna_context_penalty_weight"] = transcript_cdna_context_penalty_weight
        row["transcript_cdna_context_penalty"] = context_penalty_applied
        row["transcript_cdna_context_penalty_in_selection"] = "true" if transcript_weight > 0 and context_penalty_applied > 0 else "false"
        row["transcript_evidence_score_weight"] = transcript_weight
        row["transcript_evidence_score_in_selection"] = "true" if transcript_weight > 0 else "false"

        score = structural_score
        component_text = [structural_components] if structural_components else []
        if protein_weight > 0 or (show_zero_weight_components and selection_weights.get("protein_evidence", 0.0) != 0):
            score += protein_weight * protein_score
            component_text.append(f"+protein_evidence_score:{protein_score:.6f}*{protein_weight:.3f}")
        if transcript_weight > 0 or (show_zero_weight_components and selection_weights.get("transcript_evidence", 0.0) != 0):
            score += transcript_weight * transcript_score_context_adjusted
            component_text.append(f"+transcript_evidence_score:{transcript_score:.6f}*{transcript_weight:.3f}")
            if context_penalty_applied > 0 or show_zero_weight_components:
                component_text.append(f"-transcript_cdna_context_penalty:{context_penalty_applied:.6f}*{transcript_weight:.3f}")

        row["model_selection_score"] = score
        row["model_selection_score_scaled"] = scale_model_selection_score(
            score,
            selection_weights,
            include_protein_evidence=protein_weight > 0,
            include_transcript_evidence=transcript_weight > 0,
        )
        row["model_selection_score_components"] = ",".join(part for part in component_text if part and part != ".") or "."
        row["model_selection_score_rank"] = "."
        row["model_selection_order_rank"] = "."
        row["model_selection_tie_group"] = "."
        row["model_selection_tie_break_reason"] = "."
        row["score_delta_from_evidence"] = score - structural_score

    rows_by_sublocus: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in model_rows:
        rows_by_sublocus[str(row.get("sublocus_id", "."))].append(row)

    _assign_score_ranks(
        rows_by_sublocus,
        score_field="model_structural_support_score",
        score_rank_field="model_structural_support_score_rank",
        order_rank_field="model_structural_support_order_rank",
        tie_group_field="model_structural_support_tie_group",
        tie_break_reason_field="model_structural_support_tie_break_reason",
        tie_epsilon=tie_epsilon,
    )
    _assign_model_selection_score_ranks(rows_by_sublocus, tie_epsilon=tie_epsilon)

    for _sublocus_id, rows in rows_by_sublocus.items():
        structural_top = {str(row.get("output_transcript_id")) for row in rows if int(value_as_float(row.get("model_structural_support_order_rank"))) == 1}
        model_selection_top = {str(row.get("output_transcript_id")) for row in rows if int(value_as_float(row.get("model_selection_order_rank"))) == 1}
        top_changed = structural_top != model_selection_top
        for row in rows:
            structural_order = int(value_as_float(row.get("model_structural_support_order_rank")))
            model_selection_order = int(value_as_float(row.get("model_selection_order_rank")))
            delta = model_selection_order - structural_order
            row["order_rank_delta_from_evidence"] = delta
            row["rank_changed_by_evidence"] = "true" if delta != 0 else "false"
            row["evidence_changed_top_rank"] = "true" if top_changed and (structural_order == 1 or model_selection_order == 1) else "false"

def weighted_metric_score(components: list[tuple[str, float, Optional[float]]]) -> tuple[Optional[float], str]:
    return weighted_score(components)


def weighted_metric_value(components: list[tuple[str, float, Optional[float]]]) -> Optional[float]:
    available = [(weight, value) for _name, weight, value in components if value is not None]
    if not available:
        return None
    total_weight = sum(weight for weight, _value in available)
    if total_weight <= 0:
        return None
    return sum(weight * value for weight, value in available) / total_weight


def prf_values(prefix: str, prf: dict[str, Optional[float]]) -> dict[str, Optional[float]]:
    return {
        f"{prefix}_precision": prf.get("precision"),
        f"{prefix}_sensitivity": prf.get("sensitivity"),
        f"{prefix}_f1": prf.get("f1"),
    }


def safe_mean(values: list[Optional[float]]) -> Optional[float]:
    usable = [v for v in values if v is not None]
    if not usable:
        return None
    return sum(usable) / len(usable)
