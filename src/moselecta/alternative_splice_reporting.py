"""Concise console reporting for alternative-splice and final-annotation results."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from .alternative_splice_selection import AlternativeSpliceSelectionResult
from .final_annotation import FinalAnnotation


@dataclass(frozen=True, slots=True)
class AlternativeSpliceConsoleSummary:
    """Display-only counts copied from immutable selector/report representations."""

    preset: str
    support_policy: str
    evidence_requirement: str
    prediction_source_count: int
    min_model_support_sources: int
    min_model_support_source_fraction: float | None
    effective_min_model_support_sources: int
    max_alternatives_per_gene: int
    alternatives_emitted: bool
    uniquely_owned_candidates: int
    fixed_rule_eligible_candidates: int
    support_policy_eligible_candidates: int
    external_evidence_requirement_failures: int
    prethreshold_eligible_candidates: int
    numerical_threshold_survivors: int
    accepted_alternative_transcripts: int
    genes_with_alternative_transcripts: int
    maximum_accepted_alternatives_in_one_gene: int
    final_annotation_genes: int
    final_annotation_transcripts: int
    final_annotation_primary_transcripts: int
    final_annotation_alternative_transcripts: int


def _count(row: Mapping[str, object], field: str) -> int:
    return int(row[field])


def build_alternative_splice_console_summary(
    selection: AlternativeSpliceSelectionResult,
    run_row: Mapping[str, object],
    final_annotation: FinalAnnotation,
    *,
    alternatives_emitted: bool,
) -> AlternativeSpliceConsoleSummary:
    """Build console state without recalculating any biological classification."""
    config = selection.config
    return AlternativeSpliceConsoleSummary(
        preset=config.requested_preset.value,
        support_policy=config.support_policy.value,
        evidence_requirement=config.evidence_requirement.value,
        prediction_source_count=config.prediction_source_count or 0,
        min_model_support_sources=config.min_model_support_sources,
        min_model_support_source_fraction=(
            config.min_model_support_source_fraction
        ),
        effective_min_model_support_sources=(
            config.effective_min_model_support_sources or 2
        ),
        max_alternatives_per_gene=config.max_alternatives_per_gene,
        alternatives_emitted=alternatives_emitted,
        uniquely_owned_candidates=_count(run_row, "unique_owner_candidate_count"),
        fixed_rule_eligible_candidates=_count(
            run_row, "alternative_splice_fixed_rule_eligible_candidate_count"
        ),
        support_policy_eligible_candidates=_count(
            run_row, "alternative_splice_support_policy_eligible_candidate_count"
        ),
        external_evidence_requirement_failures=_count(
            run_row,
            "alternative_splice_external_evidence_requirement_failure_candidate_count",
        ),
        prethreshold_eligible_candidates=_count(
            run_row, "alternative_splice_prethreshold_eligible_candidate_count"
        ),
        numerical_threshold_survivors=_count(
            run_row, "alternative_splice_numerical_threshold_survivor_count"
        ),
        accepted_alternative_transcripts=_count(
            run_row, "alternative_splice_total_accepted_alternatives"
        ),
        genes_with_alternative_transcripts=_count(
            run_row, "alternative_splice_genes_with_accepted_alternatives"
        ),
        maximum_accepted_alternatives_in_one_gene=_count(
            run_row, "alternative_splice_max_accepted_alternatives_per_gene"
        ),
        final_annotation_genes=final_annotation.gene_count,
        final_annotation_transcripts=final_annotation.transcript_count,
        final_annotation_primary_transcripts=final_annotation.primary_transcript_count,
        final_annotation_alternative_transcripts=(
            final_annotation.alternative_transcript_count
        ),
    )


def alternative_splice_console_lines(
    summary: AlternativeSpliceConsoleSummary,
) -> tuple[str, ...]:
    """Render a concise, deterministic stdout summary."""
    source_requirement = (
        "Prediction-model support requirement: "
        f">={summary.effective_min_model_support_sources:,} of "
        f"{summary.prediction_source_count:,} sources"
    )
    if summary.min_model_support_source_fraction is not None:
        source_requirement += (
            " (configured fraction "
            f">={summary.min_model_support_source_fraction:.3f})"
        )
    return (
        f"Threshold preset: {summary.preset}",
        f"Support policy: {summary.support_policy}",
        f"External evidence requirement: {summary.evidence_requirement}",
        source_requirement,
        f"Maximum alternatives per gene: {summary.max_alternatives_per_gene:,}",
        "",
        "Alternatives emitted into final annotation: "
        f"{str(summary.alternatives_emitted).lower()}",
        f"Uniquely owned candidates: {summary.uniquely_owned_candidates:,}",
        f"Fixed-rule eligible candidates: {summary.fixed_rule_eligible_candidates:,}",
        "Support-policy eligible candidates: "
        f"{summary.support_policy_eligible_candidates:,}",
        "External-evidence requirement failures: "
        f"{summary.external_evidence_requirement_failures:,}",
        f"Prethreshold eligible candidates: {summary.prethreshold_eligible_candidates:,}",
        f"Numerical-threshold survivors: {summary.numerical_threshold_survivors:,}",
        "Accepted alternative transcripts (selector result): "
        f"{summary.accepted_alternative_transcripts:,}",
        "Genes with alternative transcripts (selector result): "
        f"{summary.genes_with_alternative_transcripts:,}",
        "Maximum accepted alternatives in one gene: "
        f"{summary.maximum_accepted_alternatives_in_one_gene:,}",
        "",
        "Final annotation: "
        f"genes={summary.final_annotation_genes:,}; "
        f"transcripts={summary.final_annotation_transcripts:,}",
        f"  primary transcripts: {summary.final_annotation_primary_transcripts:,}",
        "  alternative transcripts: "
        f"{summary.final_annotation_alternative_transcripts:,}",
    )
