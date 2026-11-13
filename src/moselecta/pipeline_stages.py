"""Typed orchestration boundaries used by the MoSelecta command pipeline.

This module coordinates already-established biological and annotation
components.  It does not implement selection rules or serialization.  Keeping
these boundaries outside :mod:`moselecta.pipeline` makes the final pipeline
states directly testable while preserving the lower-level dependency direction.
"""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass

from .alternative_splice import (
    AlternativeSpliceReportPolicy,
    analyze_alternative_splice_report,
)
from .alternative_splice_selection import (
    AlternativeSpliceSelectionConfig,
    AlternativeSpliceSelectionResult,
    extend_alternative_splice_report,
    select_alternative_splice_candidates,
)
from .final_annotation import FinalAnnotation, augment_annotation, build_primary_annotation
from .models import CodingModel
from .protein_evidence import ProteinEvidenceCollapsedAlignment
from .transcript_evidence import TranscriptEvidenceSuperlocus


ALTERNATIVE_SPLICE_TSV_OUTPUT_KEYS = frozenset({
    "alternative_splice_candidates_tsv",
    "alternative_splice_genes_tsv",
    "alternative_splice_run_summary_tsv",
})

ALTERNATIVE_SPLICE_TUNING_OPTIONS = (
    "--alternative-splice-preset",
    "--alternative-splice-support-policy",
    "--alternative-splice-evidence-requirement",
    "--alternative-splice-max-alternatives-per-gene",
    "--alternative-splice-min-model-support-sources",
    "--alternative-splice-min-model-support-source-fraction",
    "--alternative-splice-min-reciprocal-cds-overlap",
    "--alternative-splice-min-frame-concordance",
    "--alternative-splice-min-model-selection-score-scaled",
    "--alternative-splice-min-primary-score-ratio",
)


@dataclass(frozen=True, slots=True)
class AlternativeSpliceExecutionPlan:
    """Resolved run-state decisions for the optional selector and diagnostics."""

    enabled: bool
    audit_requested: bool
    tsv_requested: bool
    diagnostics_explicitly_requested: bool
    analysis_requested: bool
    supplied_inactive_tuning_options: tuple[str, ...]


def resolve_alternative_splice_execution_plan(
    *,
    stop_after: str,
    no_alternative_splice: bool,
    alternative_splice_report: bool,
    explicitly_requested_outputs: Collection[str],
    enabled_outputs: Collection[str],
    raw_args: Sequence[str],
) -> AlternativeSpliceExecutionPlan:
    """Resolve selector evaluation and output-routing decisions without mutation."""
    audit_requested = (
        "all" in explicitly_requested_outputs
        or "alternative_splice_candidate_audit_gff3"
        in explicitly_requested_outputs
    )
    explicitly_requested_tsv = (
        "all" in explicitly_requested_outputs
        or bool(
            ALTERNATIVE_SPLICE_TSV_OUTPUT_KEYS
            & set(explicitly_requested_outputs)
        )
    )
    tsv_requested = (
        alternative_splice_report
        or explicitly_requested_tsv
        or (
            stop_after == "final"
            and not no_alternative_splice
            and any(
                key in enabled_outputs
                for key in ALTERNATIVE_SPLICE_TSV_OUTPUT_KEYS
            )
        )
    )
    diagnostics_explicitly_requested = (
        alternative_splice_report
        or explicitly_requested_tsv
        or audit_requested
    )
    enabled = not no_alternative_splice
    analysis_requested = stop_after == "final" and (
        enabled or tsv_requested or audit_requested
    )
    supplied_inactive_tuning_options: tuple[str, ...] = ()
    if no_alternative_splice and stop_after == "final" and not analysis_requested:
        supplied_inactive_tuning_options = tuple(
            option
            for option in ALTERNATIVE_SPLICE_TUNING_OPTIONS
            if any(
                argument == option or argument.startswith(option + "=")
                for argument in raw_args
            )
        )
    return AlternativeSpliceExecutionPlan(
        enabled=enabled,
        audit_requested=audit_requested,
        tsv_requested=tsv_requested,
        diagnostics_explicitly_requested=diagnostics_explicitly_requested,
        analysis_requested=analysis_requested,
        supplied_inactive_tuning_options=supplied_inactive_tuning_options,
    )


def alternative_splice_evidence_requirement_available(
    requirement: str,
    *,
    protein_configured: bool,
    transcript_configured: bool,
) -> bool:
    """Return whether configured evidence can satisfy a selector requirement."""
    return {
        "none": True,
        "protein-or-transcript": protein_configured or transcript_configured,
        "protein": protein_configured,
        "transcript": transcript_configured,
        "protein-and-transcript": protein_configured and transcript_configured,
    }[requirement]


@dataclass(frozen=True, slots=True)
class FrozenPrimaryAnnotation:
    """The immutable retained-primary population and its public annotation."""

    models: tuple[CodingModel, ...]
    annotation: FinalAnnotation


def build_standardised_primary_annotation(
    *,
    final_selected_models: Sequence[CodingModel],
    final_pruned_models: Sequence[CodingModel] | None,
    model_score_rows_by_transcript: Mapping[str, Mapping[str, object]],
) -> FrozenPrimaryAnnotation:
    """Freeze final primaries and assign public IDs from that population alone."""
    frozen_models = tuple(
        final_pruned_models
        if final_pruned_models is not None
        else final_selected_models
    )
    annotation = build_primary_annotation(
        frozen_models,
        model_score_rows_by_transcript=model_score_rows_by_transcript,
    )
    return FrozenPrimaryAnnotation(models=frozen_models, annotation=annotation)


AlternativeSpliceReport = dict[str, list[dict[str, object]]]


@dataclass(frozen=True, slots=True)
class AlternativeSpliceStageResult:
    """Immutable boundaries for one evaluated alternative-splice selector run."""

    fact_report: AlternativeSpliceReport
    selection: AlternativeSpliceSelectionResult
    report: AlternativeSpliceReport
    final_annotation: FinalAnnotation

    @property
    def config(self) -> AlternativeSpliceSelectionConfig:
        """Return the selector's fully resolved effective configuration."""
        return self.selection.config


class AlternativeSpliceConfigurationError(ValueError):
    """An invalid selector configuration resolved at its historical stage point."""


@dataclass(frozen=True, slots=True)
class AlternativeSpliceSelectionSettings:
    """Unresolved CLI settings carried to the selector's configuration boundary."""

    preset: str
    support_policy: str
    evidence_requirement: str
    max_alternatives_per_gene: int
    min_model_support_sources: int
    min_model_support_source_fraction: float | None
    min_reciprocal_cds_overlap: float | None
    min_frame_concordance: float | None
    min_model_selection_score_scaled: float | None
    min_primary_score_ratio: float | None

    def resolve(self) -> AlternativeSpliceSelectionConfig:
        """Resolve the established selector configuration without changing values."""
        return AlternativeSpliceSelectionConfig.resolve(
            preset=self.preset,
            support_policy=self.support_policy,
            evidence_requirement=self.evidence_requirement,
            max_alternatives_per_gene=self.max_alternatives_per_gene,
            min_model_support_sources=self.min_model_support_sources,
            min_model_support_source_fraction=self.min_model_support_source_fraction,
            min_reciprocal_cds_overlap=self.min_reciprocal_cds_overlap,
            min_frame_concordance=self.min_frame_concordance,
            min_model_selection_score_scaled=self.min_model_selection_score_scaled,
            min_primary_score_ratio=self.min_primary_score_ratio,
        )


def evaluate_alternative_splice_stage(
    *,
    frozen_primary: FrozenPrimaryAnnotation,
    all_scored_models: list[CodingModel],
    model_score_rows_by_transcript: dict[str, dict[str, object]],
    protein_rows_by_transcript: dict[str, dict[str, object]],
    protein_collapsed_alignments: list[ProteinEvidenceCollapsedAlignment],
    transcript_rows_by_transcript: dict[str, dict[str, object]],
    transcript_superloci: list[TranscriptEvidenceSuperlocus],
    report_policy: AlternativeSpliceReportPolicy,
    run_metadata: dict[str, object],
    selection_settings: AlternativeSpliceSelectionSettings,
    terminal_cds_3prime_wobble_bp: int,
    alternatives_emitted: bool,
) -> AlternativeSpliceStageResult:
    """Evaluate frozen candidate facts, selection, and annotation augmentation."""
    fact_report = analyze_alternative_splice_report(
        anchor_models=list(frozen_primary.models),
        all_scored_models=all_scored_models,
        model_score_rows_by_transcript=model_score_rows_by_transcript,
        protein_rows_by_transcript=protein_rows_by_transcript,
        protein_collapsed_alignments=protein_collapsed_alignments,
        transcript_rows_by_transcript=transcript_rows_by_transcript,
        transcript_superloci=transcript_superloci,
        policy=report_policy,
        run_metadata=run_metadata,
    )
    try:
        selection_config = selection_settings.resolve()
    except ValueError as exc:
        raise AlternativeSpliceConfigurationError(str(exc)) from exc
    selection = select_alternative_splice_candidates(
        candidate_rows=fact_report["candidate_rows"],
        all_scored_models=all_scored_models,
        primary_annotation=frozen_primary.annotation,
        config=selection_config,
        terminal_cds_3prime_wobble_bp=terminal_cds_3prime_wobble_bp,
    )
    final_annotation = (
        augment_annotation(
            frozen_primary.annotation,
            selection,
            model_score_rows_by_transcript,
        )
        if alternatives_emitted
        else frozen_primary.annotation
    )
    report = extend_alternative_splice_report(
        fact_report,
        selection,
        alternatives_emitted=alternatives_emitted,
    )
    return AlternativeSpliceStageResult(
        fact_report=fact_report,
        selection=selection,
        report=report,
        final_annotation=final_annotation,
    )
