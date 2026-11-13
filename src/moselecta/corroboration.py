"""Run-aware whole-gene structural corroboration.

The comparison definition is shared by the optional scoring component, strict
R1 backfill replacement, and the candidate-wide diagnostic.  Voting is binary
and source-collapsed; the denominator is the number of distinct model-input
sources in the run other than the candidate's source.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Optional, Sequence

from .intervals import interval_list_overlap_bp_sorted
from .models import CodingModel, Sublocus
from .scoring import saturating_corroboration_score
from .tsvio import write_dict_tsv


SINGLE_EXON_THREE_PRIME_WOBBLE_BP = 3
MULTI_EXON_RECIPROCAL_CDS_COVERAGE = 0.85
SATURATING_NORMALIZATION_ALPHAS = (2.0, 4.0, 8.0, 16.0)


@dataclass(frozen=True, slots=True)
class CorroborationComparison:
    """Outcome of one complete-model pairwise corroboration comparison."""

    corroborates: bool
    reason: str
    shared_intron_count: int = 0
    candidate_cds_coverage: Optional[float] = None
    other_cds_coverage: Optional[float] = None


def strand_aware_cds_boundaries(model: CodingModel) -> tuple[int, int]:
    """Return the biological 5-prime and 3-prime CDS outer boundaries.

    Coordinates retain MoSelecta's internal 0-based half-open representation.
    Equality and absolute boundary differences are invariant to whether an end
    coordinate denotes a boundary immediately after the last covered base.
    """

    if model.strand == "+":
        return model.cds_start0, model.cds_end0
    return model.cds_end0, model.cds_start0


def intron_agreement_requirement_met(shared_count: int, intron_count: int) -> bool:
    """Apply the requested exact-intron threshold to one model direction."""

    if intron_count <= 0:
        return False
    if intron_count <= 2:
        return shared_count == intron_count
    if intron_count <= 4:
        return shared_count >= intron_count - 1
    # Integer arithmetic makes the inclusive 80% boundary deterministic.
    return shared_count * 5 >= intron_count * 4


def compare_structural_corroboration(
    candidate: CodingModel,
    other: CodingModel,
    *,
    pairwise_metrics: Optional[Mapping[str, object]] = None,
) -> CorroborationComparison:
    """Return whether ``other`` corroborates ``candidate`` as a complete model.

    ``pairwise_metrics`` may be supplied by the existing structural pairwise
    loop.  Its CDS-base precision and sensitivity are the candidate and other
    reciprocal spliced-CDS coverages, respectively, avoiding another interval
    overlap calculation in pipeline use.
    """

    if candidate.chrom != other.chrom or candidate.strand != other.strand:
        return CorroborationComparison(False, "different_chromosome_or_strand")
    if candidate.is_single_exon != other.is_single_exon:
        return CorroborationComparison(False, "exon_class_mismatch")

    if candidate.is_single_exon:
        candidate_five, candidate_three = strand_aware_cds_boundaries(candidate)
        other_five, other_three = strand_aware_cds_boundaries(other)
        if candidate_five != other_five:
            return CorroborationComparison(False, "single_exon_five_prime_boundary")
        if abs(candidate_three - other_three) > SINGLE_EXON_THREE_PRIME_WOBBLE_BP:
            return CorroborationComparison(False, "single_exon_three_prime_boundary")
        return CorroborationComparison(True, "corroborates")

    shared_introns = len(candidate.cds_intron_set & other.cds_intron_set)
    if not (
        intron_agreement_requirement_met(shared_introns, len(candidate.cds_introns))
        and intron_agreement_requirement_met(shared_introns, len(other.cds_introns))
    ):
        return CorroborationComparison(
            False,
            "intron_agreement",
            shared_intron_count=shared_introns,
        )

    candidate_coverage: Optional[float] = None
    other_coverage: Optional[float] = None
    if pairwise_metrics is not None:
        raw_candidate = pairwise_metrics.get("cds_base_precision")
        raw_other = pairwise_metrics.get("cds_base_sensitivity")
        if isinstance(raw_candidate, (int, float)) and isinstance(raw_other, (int, float)):
            candidate_coverage = float(raw_candidate)
            other_coverage = float(raw_other)
    if candidate_coverage is None or other_coverage is None:
        overlap_bp = interval_list_overlap_bp_sorted(candidate.cds_parts, other.cds_parts)
        candidate_coverage = overlap_bp / candidate.cds_length if candidate.cds_length else 0.0
        other_coverage = overlap_bp / other.cds_length if other.cds_length else 0.0

    if (
        candidate_coverage + 1e-12 < MULTI_EXON_RECIPROCAL_CDS_COVERAGE
        or other_coverage + 1e-12 < MULTI_EXON_RECIPROCAL_CDS_COVERAGE
    ):
        return CorroborationComparison(
            False,
            "reciprocal_cds_coverage",
            shared_intron_count=shared_introns,
            candidate_cds_coverage=candidate_coverage,
            other_cds_coverage=other_coverage,
        )
    return CorroborationComparison(
        True,
        "corroborates",
        shared_intron_count=shared_introns,
        candidate_cds_coverage=candidate_coverage,
        other_cds_coverage=other_coverage,
    )


def corroboration_normalizations(
    corroborating_other_source_count: int,
    run_possible_other_source_count: int,
    *,
    alphas: Sequence[float] = SATURATING_NORMALIZATION_ALPHAS,
) -> dict[str, float]:
    """Return bounded linear, logarithmic, and saturating normalizations.

    When there is no possible other source (``N == 0``), every transformation
    is defined as 0.0: independent corroboration is unavailable rather than
    complete.  ``k`` is required to be a valid source count in ``[0, N]``.
    """

    k = int(corroborating_other_source_count)
    n = int(run_possible_other_source_count)
    if n < 0 or k < 0 or k > n:
        raise ValueError(f"corroboration counts must satisfy 0 <= k <= N; got k={k}, N={n}")
    for alpha in alphas:
        if not math.isfinite(alpha) or alpha <= 0:
            raise ValueError(f"saturating alpha must be finite and > 0; got {alpha}")

    out = {
        "structural_corroboration_linear_score": 0.0,
        "structural_corroboration_log_score": 0.0,
    }
    for alpha in alphas:
        out[f"structural_corroboration_saturating_alpha_{alpha:g}_score"] = 0.0
    if n == 0:
        return out

    out["structural_corroboration_linear_score"] = k / n
    out["structural_corroboration_log_score"] = math.log1p(k) / math.log1p(n)
    for alpha in alphas:
        out[f"structural_corroboration_saturating_alpha_{alpha:g}_score"] = (
            saturating_corroboration_score(k, n, alpha)
        )
    return out


def make_structural_corroboration_row(
    candidate: CodingModel,
    sublocus: Sublocus,
    comparisons: Sequence[tuple[CodingModel, CorroborationComparison]],
    *,
    run_source_count: int,
    exact_cds_exon_chain_support_source_count: int,
) -> dict[str, object]:
    """Collapse accepted model comparisons to binary source votes."""

    corroborating_models = {
        other.output_transcript_id
        for other, result in comparisons
        if result.corroborates
    }
    corroborating_sources = {
        other.source_name
        for other, result in comparisons
        if result.corroborates
    }
    all_other_sources = {other.source_name for other, _result in comparisons}
    reason_counts: dict[str, int] = {}
    for _other, result in comparisons:
        reason_counts[result.reason] = reason_counts.get(result.reason, 0) + 1

    k = len(corroborating_sources)
    n = max(0, int(run_source_count) - 1)
    row: dict[str, object] = {
        "superlocus_id": candidate.superlocus_id or sublocus.superlocus_id,
        "sublocus_id": candidate.sublocus_id or sublocus.sublocus_id,
        "sublocus_type": sublocus.sublocus_type,
        "sublocus_flags": sublocus.flags,
        "candidate_model_id": candidate.output_transcript_id,
        "candidate_source": candidate.source_name,
        "input_gene_id": candidate.input_gene_id,
        "input_transcript_id": candidate.input_transcript_id,
        "chrom": candidate.chrom,
        "strand": candidate.strand,
        "candidate_cds_exon_count": candidate.cds_exon_count,
        "candidate_cds_intron_count": len(candidate.cds_introns),
        "candidate_cds_length_bp": candidate.cds_length,
        "run_model_source_count": int(run_source_count),
        "run_possible_other_source_count": n,
        "sublocus_model_count": len(sublocus.models),
        "sublocus_source_count": len({model.source_name for model in sublocus.models}),
        "sublocus_other_source_count": len(all_other_sources),
        "corroborating_other_source_count": k,
        "corroborating_other_sources": ",".join(sorted(corroborating_sources)) or ".",
        "corroborating_model_count": len(corroborating_models),
        "corroborating_model_ids": ",".join(sorted(corroborating_models)) or ".",
        "noncorroborating_sublocus_other_source_count": len(all_other_sources - corroborating_sources),
        "exact_cds_exon_chain_support_source_count": int(exact_cds_exon_chain_support_source_count),
        "other_source_model_comparison_count": len(comparisons),
        "corroborating_model_comparison_count": reason_counts.get("corroborates", 0),
        "exon_class_mismatch_rejected_model_comparison_count": reason_counts.get("exon_class_mismatch", 0),
        "single_exon_five_prime_rejected_model_comparison_count": reason_counts.get("single_exon_five_prime_boundary", 0),
        "single_exon_three_prime_rejected_model_comparison_count": reason_counts.get("single_exon_three_prime_boundary", 0),
        "intron_agreement_rejected_model_comparison_count": reason_counts.get("intron_agreement", 0),
        "reciprocal_cds_coverage_rejected_model_comparison_count": reason_counts.get("reciprocal_cds_coverage", 0),
        "single_exon_three_prime_wobble_bp": SINGLE_EXON_THREE_PRIME_WOBBLE_BP,
        "multi_exon_reciprocal_cds_coverage_threshold": MULTI_EXON_RECIPROCAL_CDS_COVERAGE,
    }
    row.update(corroboration_normalizations(k, n))
    return row


STRUCTURAL_CORROBORATION_COLUMNS = [
    "superlocus_id",
    "sublocus_id",
    "sublocus_type",
    "sublocus_flags",
    "candidate_model_id",
    "candidate_source",
    "input_gene_id",
    "input_transcript_id",
    "chrom",
    "strand",
    "candidate_cds_exon_count",
    "candidate_cds_intron_count",
    "candidate_cds_length_bp",
    "run_model_source_count",
    "run_possible_other_source_count",
    "sublocus_model_count",
    "sublocus_source_count",
    "sublocus_other_source_count",
    "corroborating_other_source_count",
    "corroborating_other_sources",
    "corroborating_model_count",
    "corroborating_model_ids",
    "noncorroborating_sublocus_other_source_count",
    "exact_cds_exon_chain_support_source_count",
    "structural_support_score",
    "consensus_similarity_score",
    "pairwise_similarity_score",
    "model_completeness_score",
    "exact_structure_support_score",
    "structural_corroboration_score",
    "structural_corroboration_alpha",
    "structural_corroboration_weight",
    "structural_corroboration_contribution",
    "structural_corroboration_in_selection",
    "model_structural_support_score",
    "model_structural_support_score_scaled",
    "model_selection_score",
    "model_selection_score_scaled",
    "model_selection_order_rank",
    "structural_corroboration_linear_score",
    "structural_corroboration_log_score",
    "structural_corroboration_saturating_alpha_2_score",
    "structural_corroboration_saturating_alpha_4_score",
    "structural_corroboration_saturating_alpha_8_score",
    "structural_corroboration_saturating_alpha_16_score",
    "other_source_model_comparison_count",
    "corroborating_model_comparison_count",
    "exon_class_mismatch_rejected_model_comparison_count",
    "single_exon_five_prime_rejected_model_comparison_count",
    "single_exon_three_prime_rejected_model_comparison_count",
    "intron_agreement_rejected_model_comparison_count",
    "reciprocal_cds_coverage_rejected_model_comparison_count",
    "single_exon_three_prime_wobble_bp",
    "multi_exon_reciprocal_cds_coverage_threshold",
]


def write_structural_corroboration_tsv(
    rows: list[dict[str, object]], path: Path
) -> None:
    """Write deterministic candidate-wide diagnostic rows."""

    ordered = sorted(
        rows,
        key=lambda row: (
            str(row.get("sublocus_id", ".")),
            str(row.get("candidate_model_id", ".")),
            str(row.get("candidate_source", ".")),
        ),
    )
    write_dict_tsv(ordered, path, STRUCTURAL_CORROBORATION_COLUMNS)
