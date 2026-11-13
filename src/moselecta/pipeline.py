#!/usr/bin/env python3
"""
MoSelecta pipeline

CDS-structure-aware gene-model selection with optional protein and transcript evidence integration.

CDS-structure-aware coding gene-model selection from multiple annotation sources.

MoSelecta parses one or more GFF/GFF3-like annotation files, extracts coding
transcript models, validates CDS coordinates/strands, builds stranded CDS-overlap
superloci, clusters conservative subloci, calculates coordinate-only intrinsic
and sublocus-context model metrics, evaluates structural support and source agreement, and ranks
candidate models using transparent score components.

The default pipeline scores all candidates, selects one or more non-overlapping
coding models per locus through monosublocus de-overlap and superlocus-local
backfill rescue, then applies light-touch final selected-model pruning. Outputs
are controlled by output profiles, with concise selected/final products by
default and fuller diagnostic reports available for review/debugging.

Coordinate conventions:
  * input/output GFF3 coordinates are 1-based inclusive
  * internal interval calculations use 0-based half-open coordinates
  * all structural scoring and de-overlap/pruning overlap calculations are based
    on CDS exon bases, not gene/transcript spans or UTRs
"""

from __future__ import annotations

import sys
import gc
import math
from collections import Counter
from pathlib import Path
from typing import Optional

from .constants import (
    COMPLETENESS_PENALTY_WEIGHTS,
    CONSENSUS_MULTI_SCORE_WEIGHTS,
    CONSENSUS_SINGLE_SCORE_WEIGHTS,
    LOW_SOURCE_CONFLICT_AGREEMENT_FACTOR,
    LOW_SOURCE_CONFLICT_PAIRWISE_THRESHOLD,
    PROGRAM_NAME,
    PROGRAM_VERSION,
    PRUNING_LEVEL_PRESETS,
    PRUNING_SCENARIO_THRESHOLDS,
    STRUCTURAL_MULTI_SCORE_WEIGHTS,
    STRUCTURAL_SINGLE_SCORE_WEIGHTS,
)
from .config import (
    apply_pruning_level_presets,
    build_arg_parser,
    count_pruned_single_model_subloci,
    expand_name_list,
    parse_manual_source_weights,
    parse_source_filter_values,
    validate_source_prior_args,
    with_default,
)
from .gffio import make_unique_source_names, parse_input
from .input_manifests import resolve_input_family
from .loci import (
    build_subloci,
    build_superloci,
    calculate_source_auto_weights,
    superlocus_gene_count,
    superlocus_sources,
)
from .log_utils import ProgressLogger, current_rss_mb, format_duration, peak_rss_mb
from .models import ParsedInput, ValidationFailure
from .final_annotation import FinalAnnotation
from .scoring import (
    annotate_transcript_cdna_context_concern_rows,
    apply_model_selection_scores,
    model_selection_score_positive_scale_max,
    value_as_float,
)
from .structural import (
    calculate_structural_support,
    write_model_selection_score_summary_tsv,
    write_model_selection_scores_tsv,
    write_model_structural_support_tsv,
    write_structural_features_tsv,
)
from .corroboration import write_structural_corroboration_tsv
from .structural_coherence import (
    annotate_structural_coherence,
    write_structural_coherence_regions_tsv,
)
from .validation import (
    apply_model_validation,
    calculate_all_model_intrinsic_metrics,
    collect_superlocus_validation_failures,
)
from .selection import (
    count_retained_internal_stop_evidence_overrides,
    final_prune_low_score_nonsingletons,
    final_prune_low_score_singletons,
    final_prune_opposite_strand_shadows,
    final_prune_sequence_integrity_nonsingletons,
    final_prune_singleton_selected_models,
    make_initial_final_prune_result,
    stage4_selection_from_sublocus_winners,
    stage5_backfill_selection,
    write_backfill_outputs,
    write_r1_replacement_outputs,
    write_final_pruned_outputs,
    write_monosublocus_outputs,
    write_sublocus_winners_tsv,
)
from .outputs import (
    collect_cds_intron_lengths_by_source,
    format_intron_stats,
    resolve_long_cds_intron_cutoff,
    write_combined_gff3,
    write_excluded_tsv,
    write_input_summary_tsv,
    write_members_tsv,
    write_model_intrinsic_metrics_tsv,
    write_model_subset_gff3,
    write_models_tsv,
    write_other_features_gff3,
    write_relationship_edges_tsv,
    write_source_auto_weights_tsv,
    write_subloci_tsv,
    write_sublocus_members_tsv,
    write_superloci_tsv,
    write_validation_failures_tsv,
)
from .run_context import RunContext
from .protein_evidence import (
    build_protein_subloci,
    build_protein_superloci,
    calculate_protein_model_support,
    collapse_protein_alignments,
    parse_protein_evidence,
    protein_source_labels,
    protein_summary_counts,
    protein_summary_log_lines,
    protein_summary_tsv_rows,
    release_protein_raw_alignments,
    split_path_list,
    write_protein_collapsed_alignments_tsv,
    write_protein_model_support_tsv,
    write_protein_raw_to_collapsed_tsv,
    write_protein_subloci_tsv,
    write_protein_sublocus_members_tsv,
    write_protein_summary_tsv,
    write_protein_superloci_tsv,
)
from .transcript_evidence import (
    build_transcript_subloci,
    build_transcript_superloci,
    calculate_transcript_model_support,
    collapse_transcript_evidence,
    parse_transcript_evidence,
    release_transcript_raw_transcripts,
    transcript_model_support_log_lines,
    transcript_model_support_summary_counts,
    transcript_source_labels,
    transcript_summary_counts,
    transcript_summary_log_lines,
    transcript_summary_tsv_rows,
    write_transcript_collapsed_transcripts_tsv,
    write_transcript_model_support_tsv,
    write_transcript_raw_to_collapsed_tsv,
    write_transcript_subloci_tsv,
    write_transcript_sublocus_members_tsv,
    write_transcript_summary_tsv,
    write_transcript_superloci_tsv,
)
from .genome_sequence import (
    IndexedFasta,
    calculate_and_attach_sequence_integrity,
    write_sequence_integrity_model_metrics_from_models,
    write_sequence_integrity_summary_tsv,
)
from .terminal_stop import (
    canonicalise_parsed_inputs,
    write_terminal_stop_actions_tsv,
    write_terminal_stop_audit_tsv,
    write_terminal_stop_summary_tsv,
)
from .five_prime_boundary import (
    canonicalise_parsed_input_five_prime_boundaries,
    write_five_prime_boundary_actions_tsv,
    write_five_prime_boundary_audit_tsv,
    write_five_prime_boundary_summary_tsv,
)
from .phase_recovery import (
    recover_parsed_input_phases,
    write_phase_recovery_actions_tsv,
    write_phase_recovery_audit_tsv,
    write_phase_recovery_summary_tsv,
)
from .definite_internal_stop import (
    DefiniteInternalStopResult,
    empty_definite_internal_stop_result,
    recover_phases_and_exclude_definite_internal_stops,
    write_definite_internal_stop_exclusions_tsv,
    write_definite_internal_stop_summary_tsv,
)
from .extreme_short_cds import (
    ExtremeShortCDSResult,
    exclude_extreme_short_cds,
    write_extreme_short_cds_exclusions_tsv,
    write_extreme_short_cds_summary_tsv,
)
from .cds_refinement import (
    CDSTerminusRefinementSettings,
    empty_cds_refinement_result,
    refine_cds_termini,
    write_cds_refinement_audit_tsv,
    write_cds_refinement_summary_tsv,
)
from .auto_pruning import (
    AUTO_SCENARIO_LEVELS,
    LEVEL_ORDER,
    AutoDecision,
    ScenarioSummary,
    ScenarioValidation,
    add_incremental_metrics,
    annotate_result_with_scenario_audit,
    apply_monotonic_pruning_closure,
    assess_flood_adjusted_consensus,
    automatic_diagnostics_rows,
    calculate_source_consensus_metrics,
    calculate_source_input_count_metrics,
    choose_auto_pruning_level,
    classify_evidence_adequacy,
    explain_auto_pruning_decision,
    flood_adjusted_consensus_log_lines,
    removal_provenance_from_result,
    summarize_pruning_scenario,
    validate_scenarios,
    write_auto_diagnostics,
    write_scenario_summary,
)
from .alternative_splice import (
    AlternativeSpliceReportPolicy,
    write_alternative_splice_reports,
)
from .alternative_splice_outputs import write_alternative_splice_candidate_audit_gff3
from .alternative_splice_reporting import (
    alternative_splice_console_lines,
    build_alternative_splice_console_summary,
)
from .annotation_outputs import (
    FinalAnnotationOutputMetadata,
    inherited_utr_validations,
    write_final_annotation_gff3,
    write_standardised_primary_gff3,
)
from .utr_validation import UTRValidationRecorder, write_utr_validation_summary
from .console_reporting import (
    format_counter,
    log_distribution,
    log_section,
    log_source_combination_breakdown,
    log_source_presence_breakdown,
    log_subsection,
    summarise_score_rows_by_source,
    summarise_source_scaled_scores,
)
from .pipeline_stages import (
    AlternativeSpliceConfigurationError,
    AlternativeSpliceSelectionSettings,
    alternative_splice_evidence_requirement_available,
    build_standardised_primary_annotation,
    evaluate_alternative_splice_stage,
    resolve_alternative_splice_execution_plan,
)
from .utr_augmentation import (
    UTRSettings,
    augment_final_annotation,
    build_prediction_utr_catalogue,
    build_transcript_utr_catalogue,
    learn_utr_limits_for_annotation,
    write_utr_audit_tsv,
    write_utr_augmented_gff3,
    write_utr_summary_tsv,
)
from .run_configuration import (
    build_run_configuration,
    write_run_configuration_json,
)


# Model-vs-consensus similarity is deliberately separate from the feature-support
# score above.  The support score is precision-like over the query model's own
# features; consensus similarity also asks whether the query recovers the
# supported consensus structure in the sublocus.


# Completeness/containment is deliberately separate from strict consensus
# similarity. It asks whether a model is missing source-supported coding
# structure and whether any non-consensus CDS is terminal-like versus internal.

# Model selection score is a transparent ranking helper, not a final
# selection decision. Positive terms are added and penalty terms are subtracted.


# Low-source conflict dampening is applied to the model selection score only.
# It reduces the influence of asymmetric prediction-agreement terms in two-source,
# two-model conflicts where there is no exact model-level support and the two models
# differ substantially. Raw metric columns are left unchanged.


def main(argv: Optional[list[str]] = None) -> int:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    if "--expert-help" in raw_args:
        build_arg_parser(expert_help=True).print_help()
        return 0
    if "--full-help" in raw_args:
        build_arg_parser(full_help=True).print_help()
        return 0
    args = build_arg_parser(full_help=False).parse_args(argv)
    model_inputs = resolve_input_family(
        args.inputs,
        args.source_name,
        family_name="model",
        input_option="--inputs",
        label_option="--source-name",
        comma_separated_inputs=False,
        comma_separated_labels=False,
    )
    protein_inputs = resolve_input_family(
        args.evidence_proteins,
        args.evidence_protein_source_name,
        family_name="protein evidence",
        input_option="--evidence-proteins",
        label_option="--evidence-protein-source-name",
        comma_separated_inputs=True,
        comma_separated_labels=True,
        fallback_prefix="protein",
    )
    transcript_inputs = resolve_input_family(
        args.evidence_transcripts,
        args.evidence_transcript_source_name,
        family_name="transcript evidence",
        input_option="--evidence-transcripts",
        label_option="--evidence-transcript-source-name",
        comma_separated_inputs=True,
        comma_separated_labels=True,
        fallback_prefix="transcript",
    )
    args.inputs = model_inputs.paths
    args.source_name = list(model_inputs.labels) if model_inputs.labels else None
    args.evidence_proteins = [str(path) for path in protein_inputs.paths] or None
    args.evidence_protein_source_name = (
        list(protein_inputs.labels) if protein_inputs.labels else None
    )
    args.evidence_transcripts = [str(path) for path in transcript_inputs.paths] or None
    args.evidence_transcript_source_name = (
        list(transcript_inputs.labels) if transcript_inputs.labels else None
    )
    args.auto_pruning_max_level_explicit = args.auto_pruning_max_level is not None
    if args.debug_progress:
        args.progress_log_level = "debug"
    apply_pruning_level_presets(args)

    args.strict_r1_backfill_replacement_configured = (
        args.strict_r1_backfill_replacement
    )
    args.strict_r1_backfill_replacement_explicit = any(
        option in raw_args
        for option in (
            "--strict-r1-backfill-replacement",
            "--no-strict-r1-backfill-replacement",
        )
    )
    args.strict_r1_allow_multiexon_blockers_configured = (
        args.strict_r1_allow_multiexon_blockers
    )
    args.strict_r1_allow_multiexon_blockers_explicit = any(
        option in raw_args
        for option in (
            "--strict-r1-allow-multiexon-blockers",
            "--no-strict-r1-allow-multiexon-blockers",
        )
    )
    if (
        args.strict_r1_backfill_replacement
        and args.stop_after not in {"backfill", "final"}
    ):
        if "--strict-r1-backfill-replacement" in raw_args:
            raise SystemExit(
                "--strict-r1-backfill-replacement requires --stop-after "
                "backfill or final."
            )
        args.strict_r1_backfill_replacement = False
    if (
        not args.strict_r1_backfill_replacement
        and args.strict_r1_allow_multiexon_blockers
    ):
        if "--strict-r1-allow-multiexon-blockers" in raw_args:
            raise SystemExit(
                "--strict-r1-allow-multiexon-blockers requires "
                "--strict-r1-backfill-replacement."
            )
        args.strict_r1_allow_multiexon_blockers = False

    # The full pipeline is now the default. --stop-after is the user-facing way to halt earlier.
    args.run_monosublocus_selection = args.stop_after in {"monosublocus", "backfill", "final"}
    args.run_backfill_selection = args.stop_after in {"backfill", "final"}

    context = RunContext.from_args(args)
    timer = context.timer
    progress = context.progress
    output_enabled = context.output_enabled
    note_output = context.note_output
    if args.strict_r1_backfill_replacement:
        if not args.run_backfill_selection:
            raise SystemExit(
                "--strict-r1-backfill-replacement requires --stop-after backfill or final."
            )
        context.ensure_output("r1_backfill_replacements_tsv")
        context.ensure_output("r1_backfill_summary_tsv")
    for message in (
        *model_inputs.expansion_messages,
        *protein_inputs.expansion_messages,
        *transcript_inputs.expansion_messages,
    ):
        progress.log(message)

    # Explicit stop stages retain useful terminal products even when those files
    # are otherwise reserved for a fuller reporting profile. Explicit skips are
    # authoritative and are enforced by RunContext.ensure_output().
    if args.stop_after == "scored-candidates":
        for key in (
            "combined_gff3",
            "model_selection_scores_tsv",
            "model_structural_support_tsv",
            "model_selection_score_summary_tsv",
        ):
            context.ensure_output(key)
    elif args.stop_after == "monosublocus":
        for key in ("selected_models_gff3", "selected_models_tsv"):
            context.ensure_output(key)
    elif args.stop_after == "backfill":
        for key in ("final_selected_models_gff3", "final_selected_models_tsv"):
            context.ensure_output(key)
    if args.add_utrs:
        if args.stop_after != "final":
            raise SystemExit("--add-utrs requires --stop-after final.")
        for key in (
            "utr_augmented_final_annotation_gff3",
            "utr_augmentation_summary_tsv",
        ):
            context.ensure_output(key)
    if args.prune_singleton_subloci_from_scored_gff3:
        context.ensure_output("scored_pruned_gff3")
    automatic_assessment_requested = (
        args.pruning_level in {"auto", "auto-conservative"}
        or bool(args.auto_pruning_report_only)
    )
    final_selected_pruning_requested = args.run_backfill_selection and (
        automatic_assessment_requested
        or
        args.final_prune_singleton_subloci_level != "none"
        or bool(args.final_prune_low_score_singletons)
        or bool(args.final_prune_low_score_nonsingletons)
        or bool(args.final_prune_sequence_integrity_nonsingletons)
        or bool(args.final_prune_opposite_strand_shadows)
    )
    if automatic_assessment_requested:
        context.ensure_output("pruning_scenario_summary_tsv")
        context.ensure_output("auto_pruning_diagnostics_tsv")
    explicitly_requested_outputs = expand_name_list(args.write_output)
    alternative_splice_plan = resolve_alternative_splice_execution_plan(
        stop_after=args.stop_after,
        no_alternative_splice=args.no_alternative_splice,
        alternative_splice_report=args.alternative_splice_report,
        explicitly_requested_outputs=explicitly_requested_outputs,
        enabled_outputs=context.enabled_outputs,
        raw_args=raw_args,
    )
    alternative_splice_audit_requested = alternative_splice_plan.audit_requested
    alternative_splice_analysis_requested = alternative_splice_plan.analysis_requested
    alternative_splice_enabled = alternative_splice_plan.enabled
    if (
        alternative_splice_plan.diagnostics_explicitly_requested
        and args.stop_after != "final"
    ):
        raise SystemExit(
            "Alternative-splice diagnostics require --stop-after final."
        )
    if args.alternative_splice_report:
        for key in (
            "alternative_splice_candidates_tsv",
            "alternative_splice_genes_tsv",
            "alternative_splice_run_summary_tsv",
        ):
            context.ensure_output(key)
    protein_evidence_paths = split_path_list(args.evidence_proteins)
    transcript_evidence_paths = split_path_list(args.evidence_transcripts)
    if alternative_splice_analysis_requested:
        evidence_requirement = args.alternative_splice_evidence_requirement
        if not alternative_splice_evidence_requirement_available(
            evidence_requirement,
            protein_configured=bool(protein_evidence_paths),
            transcript_configured=bool(transcript_evidence_paths),
        ):
            raise SystemExit(
                "--alternative-splice-evidence-requirement "
                f"{evidence_requirement} requires compatible configured evidence inputs."
            )
    elif alternative_splice_plan.supplied_inactive_tuning_options:
        print(
            "Warning: alternative-splice selector options have no effect because "
            "--no-alternative-splice was used without a diagnostic report or audit: "
            + ", ".join(
                alternative_splice_plan.supplied_inactive_tuning_options
            ),
            file=sys.stderr,
        )
    genome_fasta_path = args.genome_fasta
    if args.final_prune_rescue_protein_supported_singletons is None:
        args.final_prune_rescue_protein_supported_singletons = bool(protein_evidence_paths)
    else:
        args.final_prune_rescue_protein_supported_singletons = bool(args.final_prune_rescue_protein_supported_singletons) and bool(protein_evidence_paths)
    if args.final_prune_rescue_transcript_supported_singletons is None:
        args.final_prune_rescue_transcript_supported_singletons = bool(transcript_evidence_paths)
    else:
        args.final_prune_rescue_transcript_supported_singletons = bool(args.final_prune_rescue_transcript_supported_singletons) and bool(transcript_evidence_paths)
    # Evidence outputs are controlled by --output-profile/--write-output.
    # Transcript evidence is scored separately and can contribute to selection and
    # evidence-aware pruning when the corresponding options/presets are active.
    progress.log(
        f"Starting {PROGRAM_NAME} v{PROGRAM_VERSION} run: inputs={len(args.inputs):,}; out_prefix={args.out_prefix}; "
        f"stop_after={args.stop_after}; output_profile={args.output_profile}; pruning_level={args.pruning_level}"
    )
    validate_source_prior_args(args)
    if args.terminal_cds_3prime_wobble_bp < 0:
        raise SystemExit("--terminal-cds-3prime-wobble-bp must be >= 0.")
    protein_evidence_wobble_bp = (
        args.protein_evidence_terminal_cds_3prime_wobble_bp
        if args.protein_evidence_terminal_cds_3prime_wobble_bp is not None
        else args.terminal_cds_3prime_wobble_bp
    )
    if protein_evidence_wobble_bp < 0:
        raise SystemExit("--protein-evidence-terminal-cds-3prime-wobble-bp must be >= 0.")
    if not (0 <= args.protein_evidence_support_score_threshold <= 1):
        raise SystemExit("--protein-evidence-support-score-threshold must be between 0 and 1.")
    if not (0 <= args.single_exon_min_reciprocal_overlap <= 1):
        raise SystemExit("--single-exon-min-reciprocal-overlap must be between 0 and 1.")
    if args.single_exon_support_min_reciprocal_overlap is not None and not (0 <= args.single_exon_support_min_reciprocal_overlap <= 1):
        raise SystemExit("--single-exon-support-min-reciprocal-overlap must be between 0 and 1.")
    if not (0 <= args.long_cds_intron_percentile <= 100):
        raise SystemExit("--long-cds-intron-percentile must be between 0 and 100.")
    if not (0 <= args.singleton_rescue_min_transcript_cdna_context_fraction <= 1):
        raise SystemExit(
            "--singleton-rescue-min-transcript-cdna-context-fraction must be between 0 and 1."
        )
    if args.model_selection_tie_epsilon < 0:
        raise SystemExit("--model-selection-tie-epsilon must be >= 0.")
    if args.selection_penalty_transcript_cdna_context < 0:
        raise SystemExit("--selection-penalty-transcript-cdna-context must be >= 0.")
    if args.selection_penalty_structural_coherence < 0:
        raise SystemExit("--selection-penalty-structural-coherence must be >= 0.")
    for option_name in (
        "auto_pruning_max_extra_removal_fraction",
        "auto_pruning_high_min_singleton_fraction",
        "auto_pruning_high_min_prospective_removal_fraction",
    ):
        option_value = float(getattr(args, option_name))
        if not 0.0 <= option_value <= 1.0:
            raise SystemExit(
                "--" + option_name.replace("_", "-") + " must be between 0 and 1."
            )
    if LEVEL_ORDER[args.auto_pruning_min_level] > LEVEL_ORDER[args.auto_pruning_max_level]:
        raise SystemExit("--auto-pruning-min-level cannot be more aggressive than --auto-pruning-max-level.")
    if args.genetic_code != 1:
        raise SystemExit(f"{PROGRAM_NAME} {PROGRAM_VERSION} currently supports --genetic-code 1 only.")
    if args.minimum_coding_length_bp <= 0:
        raise SystemExit("--minimum-coding-length-bp must be > 0.")
    if args.selection_penalty_sequence_integrity < 0:
        raise SystemExit("--selection-penalty-sequence-integrity must be >= 0.")
    if args.sequence_internal_stop_count_cap < 1:
        raise SystemExit("--sequence-internal-stop-count-cap must be >= 1.")
    for option_name, value in (
        ("--sequence-penalty-internal-stop", args.sequence_penalty_internal_stop),
        ("--sequence-penalty-phase-inconsistent", args.sequence_penalty_phase_inconsistent),
        ("--sequence-penalty-missing-stop", args.sequence_penalty_missing_stop),
        ("--sequence-penalty-missing-start", args.sequence_penalty_missing_start),
        ("--sequence-penalty-noncanonical-splice", args.sequence_penalty_noncanonical_splice),
    ):
        if value < 0:
            raise SystemExit(f"{option_name} must be >= 0.")
    for option_name, value in (
        ("--unsupported-short-5prime-terminal-cds-exon-bp", args.unsupported_short_5prime_terminal_cds_exon_bp),
        ("--unsupported-short-3prime-terminal-cds-exon-bp", args.unsupported_short_3prime_terminal_cds_exon_bp),
        ("--unsupported-short-internal-cds-exon-bp", args.unsupported_short_internal_cds_exon_bp),
    ):
        if value < 1:
            raise SystemExit(f"{option_name} must be >= 1.")
    unsupported_short_component_weights = (
        args.unsupported_short_5prime_terminal_weight,
        args.unsupported_short_3prime_terminal_weight,
        args.unsupported_short_internal_weight,
    )
    if any(value < 0 for value in unsupported_short_component_weights):
        raise SystemExit("Unsupported-short CDS exon component weights must be >= 0.")
    if sum(unsupported_short_component_weights) <= 0:
        raise SystemExit("At least one unsupported-short CDS exon component weight must be > 0.")
    if not (0 < args.transcript_cdna_context_neutral_fraction <= 1):
        raise SystemExit("--transcript-cdna-context-neutral-fraction must be > 0 and <= 1.")
    if not (0 <= args.transcript_cdna_context_severe_fraction <= args.transcript_cdna_context_neutral_fraction):
        raise SystemExit(
            "--transcript-cdna-context-severe-fraction must be >= 0 and <= "
            "--transcript-cdna-context-neutral-fraction."
        )
    if not (0 < args.monosublocus_min_cds_exon_overlap_fraction <= 1):
        raise SystemExit("--monosublocus-min-cds-exon-overlap-fraction must be > 0 and <= 1.")
    # Final selected-model pruning is only applicable when final_selected_models are produced.
    # With earlier --stop-after values, pruning-level presets are ignored.
    if args.shadow_prune_min_overlap_bp < 1:
        raise SystemExit("--shadow-prune-min-overlap-bp must be >= 1.")
    if not (0 < args.shadow_prune_min_loser_cds_overlap_fraction <= 1):
        raise SystemExit("--shadow-prune-min-loser-cds-overlap-fraction must be > 0 and <= 1.")
    if args.shadow_prune_min_score_margin_scaled < 0:
        raise SystemExit("--shadow-prune-min-score-margin-scaled must be >= 0.")
    if args.singleton_rescue_min_monosublocus_member_count < 2:
        raise SystemExit("--singleton-rescue-min-monosublocus-member-count must be >= 2.")
    if args.singleton_rescue_min_cds_overlap_bp < 1:
        raise SystemExit("--singleton-rescue-min-cds-overlap-bp must be >= 1.")
    if not (0 < args.singleton_rescue_min_self_overlap_fraction <= 1):
        raise SystemExit("--singleton-rescue-min-self-overlap-fraction must be > 0 and <= 1.")
    if not (0 < args.singleton_rescue_min_reciprocal_overlap_fraction <= 1):
        raise SystemExit("--singleton-rescue-min-reciprocal-overlap-fraction must be > 0 and <= 1.")
    if not (0 <= args.singleton_rescue_min_protein_evidence_score_scaled <= 1):
        raise SystemExit("--singleton-rescue-min-protein-evidence-score-scaled must be between 0 and 1.")
    if args.singleton_rescue_min_protein_exact_chain_source_count < 0:
        raise SystemExit("--singleton-rescue-min-protein-exact-chain-source-count must be >= 0.")
    for option_name, value in (
        ("--low-score-singleton-max-model-structural-support-score-scaled", args.low_score_singleton_max_model_structural_support_score_scaled),
        ("--low-score-singleton-max-protein-evidence-score-scaled", args.low_score_singleton_max_protein_evidence_score_scaled),
        ("--low-score-singleton-max-transcript-evidence-score-scaled", args.low_score_singleton_max_transcript_evidence_score_scaled),
        ("--low-score-singleton-max-model-selection-score-scaled", args.low_score_singleton_max_model_selection_score_scaled),
        ("--low-score-singleton-min-protein-overlap-model-fraction", args.low_score_singleton_min_protein_overlap_model_fraction),
        ("--low-score-singleton-min-protein-supported-model-fraction", args.low_score_singleton_min_protein_supported_model_fraction),
        ("--low-score-singleton-min-transcript-scoring-eligible-model-fraction", args.low_score_singleton_min_transcript_scoring_eligible_model_fraction),
        ("--low-score-singleton-min-transcript-supported-model-fraction", args.low_score_singleton_min_transcript_supported_model_fraction),
        ("--low-score-nonsingleton-max-model-structural-support-score-scaled", args.low_score_nonsingleton_max_model_structural_support_score_scaled),
        ("--low-score-nonsingleton-max-protein-evidence-score-scaled", args.low_score_nonsingleton_max_protein_evidence_score_scaled),
        ("--low-score-nonsingleton-max-transcript-evidence-score-scaled", args.low_score_nonsingleton_max_transcript_evidence_score_scaled),
        ("--low-score-nonsingleton-max-model-selection-score-scaled", args.low_score_nonsingleton_max_model_selection_score_scaled),
        ("--low-score-nonsingleton-min-protein-overlap-model-fraction", args.low_score_nonsingleton_min_protein_overlap_model_fraction),
        ("--low-score-nonsingleton-min-protein-supported-model-fraction", args.low_score_nonsingleton_min_protein_supported_model_fraction),
        ("--low-score-nonsingleton-min-transcript-scoring-eligible-model-fraction", args.low_score_nonsingleton_min_transcript_scoring_eligible_model_fraction),
        ("--low-score-nonsingleton-min-transcript-supported-model-fraction", args.low_score_nonsingleton_min_transcript_supported_model_fraction),
    ):
        if not (0 <= value <= 1):
            raise SystemExit(f"{option_name} must be between 0 and 1.")
    selection_weight_values = [
        args.selection_weight_structural_support,
        args.selection_weight_consensus_similarity,
        args.selection_weight_pairwise_similarity,
        args.selection_weight_completeness,
        args.selection_weight_source_prior,
        args.selection_weight_terminal_extension,
        args.selection_weight_exact_structure_support,
        args.selection_weight_structural_corroboration,
        args.selection_weight_protein_evidence,
        args.selection_weight_transcript_evidence,
        args.selection_penalty_transcript_cdna_context,
        args.selection_penalty_unsupported_short_cds_exon,
        args.selection_penalty_structural_coherence,
        args.unsupported_short_5prime_terminal_weight,
        args.unsupported_short_3prime_terminal_weight,
        args.unsupported_short_internal_weight,
        args.selection_penalty_intrinsic,
        args.selection_penalty_partial_structure,
        args.selection_penalty_internal_extra,
        args.selection_penalty_terminal_extra,
    ]
    if any(value < 0 for value in selection_weight_values):
        raise SystemExit("Model selection weights/penalties must be >= 0.")
    if not math.isfinite(args.structural_corroboration_alpha) or args.structural_corroboration_alpha <= 0:
        raise SystemExit("--structural-corroboration-alpha must be finite and > 0.")

    source_names = make_unique_source_names(args.inputs, args.source_name)
    manual_source_weights = parse_manual_source_weights(args.source_weight, source_names)
    parsed_inputs: list[ParsedInput] = []
    indexed_fasta: Optional[IndexedFasta] = None
    genome_fasta_index_action = "."
    genome_fasta_index_rebuild_reason = "."
    genome_fasta_index_path: Optional[Path] = None
    terminal_stop_canonicalisation_counts: Counter[str] = Counter()
    five_prime_boundary_canonicalisation_counts: Counter[str] = Counter()
    cds_phase_recovery_counts: Counter[str] = Counter()
    terminal_stop_audit_models = []
    five_prime_boundary_audit_models = []
    cds_phase_recovery_audit_models = []
    definite_internal_stop_result: DefiniteInternalStopResult | None = None
    extreme_short_cds_result: ExtremeShortCDSResult | None = None

    for input_index, (path, source_name) in enumerate(zip(args.inputs, source_names), start=1):
        progress.log(f"Parsing input {input_index:,}/{len(args.inputs):,}: source={source_name}; file={path}", level="verbose")
        parsed_inputs.append(parse_input(path, source_name, strict=args.strict))
        progress.log(
            f"Parsed {source_name}: coding_models={len(parsed_inputs[-1].coding_models):,}; "
            f"features={len(parsed_inputs[-1].records):,}",
            level="verbose",
        )
    progress.log(f"Parsed all inputs; total coding models before validation={sum(len(p.coding_models) for p in parsed_inputs):,}")
    input_coding_model_counts = {
        source_name: len(parsed.coding_models)
        for source_name, parsed in zip(source_names, parsed_inputs)
    }
    input_count_metrics = calculate_source_input_count_metrics(
        input_coding_model_counts
    )
    timer.mark("Input parsing")
    definite_internal_stop_result = empty_definite_internal_stop_result(
        sum(len(parsed.coding_models) for parsed in parsed_inputs)
    )

    if genome_fasta_path is not None:
        progress.log(f"Opening indexed genome FASTA for early CDS-boundary canonicalisation: {genome_fasta_path}")
        try:
            indexed_fasta = IndexedFasta(genome_fasta_path)
        except (OSError, ValueError) as exc:
            raise SystemExit(f"Genome FASTA error: {exc}") from exc
        genome_fasta_index_action = indexed_fasta.index_action
        genome_fasta_index_rebuild_reason = indexed_fasta.index_rebuild_reason
        genome_fasta_index_path = indexed_fasta.index_path
        if genome_fasta_index_action == "created":
            progress.log(f"Built missing genome FASTA index automatically: {genome_fasta_index_path}")
        elif genome_fasta_index_action == "rebuilt":
            progress.log(
                "Existing genome FASTA index was stale or invalid and was rebuilt automatically: "
                f"{genome_fasta_index_path}; reason={genome_fasta_index_rebuild_reason}"
            )
        with indexed_fasta as fasta_reader:
            five_prime_boundary_canonicalisation_counts = (
                canonicalise_parsed_input_five_prime_boundaries(
                    parsed_inputs, fasta_reader
                )
            )
            terminal_stop_canonicalisation_counts = canonicalise_parsed_inputs(
                parsed_inputs, fasta_reader
            )
            if args.definite_internal_stop_exclusion:
                (
                    cds_phase_recovery_counts,
                    early_audit_models,
                    definite_internal_stop_result,
                ) = recover_phases_and_exclude_definite_internal_stops(
                    parsed_inputs, fasta_reader
                )
            else:
                cds_phase_recovery_counts = recover_parsed_input_phases(
                    parsed_inputs, fasta_reader
                )
                early_audit_models = [
                    model
                    for parsed in parsed_inputs
                    for model in parsed.coding_models
                ]
                definite_internal_stop_result = empty_definite_internal_stop_result(
                    len(early_audit_models), execution_status="disabled_by_user"
                )
            extreme_short_cds_result = exclude_extreme_short_cds(
                parsed_inputs,
                minimum_coding_length_bp=args.minimum_coding_length_bp,
                fasta=fasta_reader,
            )
        terminal_stop_audit_models = early_audit_models
        five_prime_boundary_audit_models = terminal_stop_audit_models
        cds_phase_recovery_audit_models = terminal_stop_audit_models
        progress.log(
            "Early class-B 5-prime CDS-boundary canonicalisation: "
            + ", ".join(
                f"{key}={value:,}"
                for key, value in sorted(five_prime_boundary_canonicalisation_counts.items())
            )
        )
        timer.mark("5-prime CDS-boundary canonicalisation")
        progress.log(
            "Early terminal-stop canonicalisation: "
            + ", ".join(
                f"{key}={value:,}"
                for key, value in sorted(terminal_stop_canonicalisation_counts.items())
            )
        )
        timer.mark("Terminal-stop canonicalisation")
        progress.log(
            "Early biological CDS phase recovery: "
            + ", ".join(
                f"{key}={value:,}"
                for key, value in sorted(cds_phase_recovery_counts.items())
            )
        )
        timer.mark("CDS phase recovery")
        progress.log(
            "Early definite internal-stop exclusion: "
            f"status={definite_internal_stop_result.execution_status}; "
            f"excluded={definite_internal_stop_result.excluded_models:,}; "
            f"resolved={definite_internal_stop_result.excluded_resolved_frame_internal_stop:,}; "
            f"unresolved_all_frames={definite_internal_stop_result.excluded_unresolved_all_frames_internal_stop:,}; "
            f"invalid_cds_geometry={definite_internal_stop_result.excluded_invalid_cds_geometry:,}"
        )
        timer.mark("Definite internal-stop exclusion")

    if extreme_short_cds_result is None:
        extreme_short_cds_result = exclude_extreme_short_cds(
            parsed_inputs,
            minimum_coding_length_bp=args.minimum_coding_length_bp,
        )
    progress.log(
        "Early extreme-short CDS exclusion: "
        f"threshold={extreme_short_cds_result.minimum_coding_length_bp:,}; "
        f"evaluated={extreme_short_cds_result.candidates_evaluated:,}; "
        f"excluded={extreme_short_cds_result.excluded_models:,}; "
        f"retained={extreme_short_cds_result.retained_models:,}"
    )
    timer.mark("Extreme-short CDS exclusion")

    validation_failures: list[ValidationFailure] = []
    model_failures = apply_model_validation(
        parsed_inputs,
        action=args.model_validation_action,
        max_model_cds_span_bp=args.max_model_cds_span_bp,
        max_model_cds_intron_bp=args.max_model_cds_intron_bp,
    )
    validation_failures.extend(model_failures)

    coding_models = [model for parsed in parsed_inputs for model in parsed.coding_models]
    excluded_models = [item for parsed in parsed_inputs for item in parsed.excluded_models]

    combined_gff3 = with_default(args.combined_gff3, args.out_prefix, ".coding_models_with_superloci.gff3")
    superloci_tsv = with_default(args.superloci_tsv, args.out_prefix, ".superloci.tsv")
    members_tsv = with_default(args.members_tsv, args.out_prefix, ".superlocus_members.tsv")
    models_tsv = with_default(args.models_tsv, args.out_prefix, ".coding_models.tsv")
    excluded_tsv = with_default(args.excluded_models_tsv, args.out_prefix, ".excluded_models.tsv")
    other_gff3 = with_default(args.other_features_gff3, args.out_prefix, ".other_or_no_cds_features.gff3")
    input_summary_tsv = with_default(args.input_summary_tsv, args.out_prefix, ".input_summary.tsv")
    validation_failures_tsv = with_default(args.validation_failures_tsv, args.out_prefix, ".validation_failures.tsv")
    subloci_tsv = with_default(args.subloci_tsv, args.out_prefix, ".subloci.tsv")
    sublocus_members_tsv = with_default(args.sublocus_members_tsv, args.out_prefix, ".sublocus_members.tsv")
    relationship_edges_tsv = with_default(args.relationship_edges_tsv, args.out_prefix, ".relationship_edges.tsv")
    model_metrics_tsv = with_default(args.model_metrics_tsv, args.out_prefix, ".model_intrinsic_metrics.tsv")
    sequence_integrity_summary_tsv = with_default(args.sequence_integrity_summary_tsv, args.out_prefix, ".sequence_integrity_summary.tsv")
    sequence_integrity_model_metrics_tsv = with_default(args.sequence_integrity_model_metrics_tsv, args.out_prefix, ".sequence_integrity_model_metrics.tsv")
    run_configuration_json = with_default(
        args.run_configuration_json, args.out_prefix, ".run_configuration.json"
    )
    terminal_stop_canonicalisation_summary_tsv = with_default(
        args.terminal_stop_canonicalisation_summary_tsv,
        args.out_prefix,
        ".terminal_stop_canonicalisation.summary.tsv",
    )
    terminal_stop_canonicalisation_actions_tsv = with_default(
        args.terminal_stop_canonicalisation_actions_tsv,
        args.out_prefix,
        ".terminal_stop_canonicalisation.actions.tsv",
    )
    terminal_stop_canonicalisation_audit_tsv = with_default(
        args.terminal_stop_canonicalisation_audit_tsv,
        args.out_prefix,
        ".terminal_stop_canonicalisation.audit.tsv",
    )
    five_prime_boundary_canonicalisation_summary_tsv = with_default(
        args.cds_5prime_boundary_canonicalisation_summary_tsv,
        args.out_prefix,
        ".cds_5prime_boundary_canonicalisation.summary.tsv",
    )
    five_prime_boundary_canonicalisation_actions_tsv = with_default(
        args.cds_5prime_boundary_canonicalisation_actions_tsv,
        args.out_prefix,
        ".cds_5prime_boundary_canonicalisation.actions.tsv",
    )
    five_prime_boundary_canonicalisation_audit_tsv = with_default(
        args.cds_5prime_boundary_canonicalisation_audit_tsv,
        args.out_prefix,
        ".cds_5prime_boundary_canonicalisation.audit.tsv",
    )
    cds_phase_recovery_summary_tsv = with_default(
        args.cds_phase_recovery_summary_tsv,
        args.out_prefix,
        ".cds_phase_recovery.summary.tsv",
    )
    cds_phase_recovery_actions_tsv = with_default(
        args.cds_phase_recovery_actions_tsv,
        args.out_prefix,
        ".cds_phase_recovery.actions.tsv",
    )
    cds_phase_recovery_audit_tsv = with_default(
        args.cds_phase_recovery_audit_tsv,
        args.out_prefix,
        ".cds_phase_recovery.audit.tsv",
    )
    definite_internal_stop_exclusion_summary_tsv = with_default(
        args.definite_internal_stop_exclusion_summary_tsv,
        args.out_prefix,
        ".definite_internal_stop_exclusion.summary.tsv",
    )
    definite_internal_stop_exclusions_tsv = with_default(
        args.definite_internal_stop_exclusions_tsv,
        args.out_prefix,
        ".definite_internal_stop_exclusions.tsv",
    )
    extreme_short_cds_exclusion_summary_tsv = with_default(
        args.extreme_short_cds_exclusion_summary_tsv,
        args.out_prefix,
        ".extreme_short_cds_exclusion.summary.tsv",
    )
    extreme_short_cds_exclusions_tsv = with_default(
        args.extreme_short_cds_exclusions_tsv,
        args.out_prefix,
        ".extreme_short_cds_exclusions.tsv",
    )
    source_auto_weights_tsv = with_default(args.source_auto_weights_tsv, args.out_prefix, ".source_auto_weights.tsv")
    model_structural_support_tsv = with_default(args.model_structural_support_tsv, args.out_prefix, ".model_structural_support.tsv")
    structural_corroboration_tsv = with_default(
        args.structural_corroboration_tsv,
        args.out_prefix,
        ".structural_corroboration.tsv",
    )
    structural_coherence_regions_tsv = with_default(args.structural_coherence_regions_tsv, args.out_prefix, ".structural_coherence_regions.tsv")
    structural_features_tsv = with_default(args.structural_features_tsv, args.out_prefix, ".structural_features.tsv")
    model_selection_scores_tsv = with_default(args.model_selection_scores_tsv, args.out_prefix, ".model_selection_scores.tsv")
    model_selection_score_summary_tsv = with_default(args.model_selection_score_summary_tsv, args.out_prefix, ".model_selection_score_summary.tsv")
    scored_pruned_gff3 = with_default(
        args.scored_pruned_gff3,
        args.out_prefix,
        ".coding_models_with_superloci.no_singleton_subloci.gff3",
    )
    sublocus_winners_gff3 = with_default(args.sublocus_winners_gff3, args.out_prefix, ".sublocus_winners.gff3")
    sublocus_winners_tsv = with_default(args.sublocus_winners_tsv, args.out_prefix, ".sublocus_winners.tsv")
    monosubloci_tsv = with_default(args.monosubloci_tsv, args.out_prefix, ".monosubloci.tsv")
    monosublocus_members_tsv = with_default(args.monosublocus_members_tsv, args.out_prefix, ".monosublocus_members.tsv")
    monosublocus_edges_tsv = with_default(args.monosublocus_edges_tsv, args.out_prefix, ".monosublocus_edges.tsv")
    monosublocus_review_gff3 = with_default(args.monosublocus_review_gff3, args.out_prefix, ".monosublocus_review.gff3")
    selected_models_gff3 = with_default(args.selected_models_gff3, args.out_prefix, ".selected_models.gff3")
    selected_models_tsv = with_default(args.selected_models_tsv, args.out_prefix, ".selected_models.tsv")
    suppressed_monosublocus_models_tsv = with_default(args.suppressed_monosublocus_models_tsv, args.out_prefix, ".suppressed_monosublocus_models.tsv")
    backfill_candidates_tsv = with_default(args.backfill_candidates_tsv, args.out_prefix, ".backfill_candidates.tsv")
    backfilled_models_tsv = with_default(args.backfilled_models_tsv, args.out_prefix, ".backfilled_models.tsv")
    backfill_rejected_models_tsv = with_default(args.backfill_rejected_models_tsv, args.out_prefix, ".backfill_rejected_models.tsv")
    r1_backfill_replacements_tsv = with_default(
        args.r1_backfill_replacements_tsv,
        args.out_prefix,
        ".r1_backfill_replacements.tsv",
    )
    r1_backfill_summary_tsv = with_default(
        args.r1_backfill_summary_tsv,
        args.out_prefix,
        ".r1_backfill_summary.tsv",
    )
    final_selected_models_tsv = with_default(args.final_selected_models_tsv, args.out_prefix, ".final_selected_models.tsv")
    final_selected_models_gff3 = with_default(args.final_selected_models_gff3, args.out_prefix, ".final_selected_models.gff3")
    final_pruned_selected_models_tsv = with_default(args.final_pruned_selected_models_tsv, args.out_prefix, ".final_selected_models.pruned.tsv")
    final_pruned_selected_models_gff3 = with_default(args.final_pruned_selected_models_gff3, args.out_prefix, ".final_selected_models.pruned.gff3")
    final_pruned_selected_removed_tsv = with_default(args.final_pruned_selected_removed_tsv, args.out_prefix, ".final_selected_models.pruned.removed.tsv")
    final_shadow_gene_pairs_tsv = with_default(args.final_shadow_gene_pairs_tsv, args.out_prefix, ".final_shadow_gene_pairs.tsv")
    pruning_scenario_summary_tsv = with_default(args.pruning_scenario_summary_tsv, args.out_prefix, ".pruning_scenario_summary.tsv")
    auto_pruning_diagnostics_tsv = with_default(args.auto_pruning_diagnostics_tsv, args.out_prefix, ".auto_pruning_diagnostics.tsv")
    standardised_primary_models_gff3 = with_default(
        args.standardised_primary_models_gff3,
        args.out_prefix,
        ".standardised_primary_models.gff3",
    )
    final_annotation_gff3 = with_default(
        args.final_annotation_gff3,
        args.out_prefix,
        ".final_annotation.gff3",
    )
    utr_augmented_final_annotation_gff3 = with_default(
        args.utr_augmented_final_annotation_gff3,
        args.out_prefix,
        ".final_annotation.with_utrs.gff3",
    )
    utr_augmentation_audit_tsv = with_default(
        args.utr_augmentation_audit_tsv,
        args.out_prefix,
        ".utr_augmentation.audit.tsv",
    )
    utr_junction_validation_audit_tsv = with_default(
        args.utr_junction_validation_audit_tsv,
        args.out_prefix,
        ".utr_junction_validation.audit.tsv",
    )
    utr_junction_validation_summary_tsv = with_default(
        args.utr_junction_validation_summary_tsv,
        args.out_prefix,
        ".utr_junction_validation.summary.tsv",
    )
    utr_augmentation_summary_tsv = with_default(
        args.utr_augmentation_summary_tsv,
        args.out_prefix,
        ".utr_augmentation.summary.tsv",
    )
    cds_terminus_refinement_audit_tsv = with_default(
        args.cds_terminus_refinement_audit_tsv,
        args.out_prefix,
        ".cds_terminus_refinement.audit.tsv",
    )
    cds_terminus_refinement_summary_tsv = with_default(
        args.cds_terminus_refinement_summary_tsv,
        args.out_prefix,
        ".cds_terminus_refinement.summary.tsv",
    )
    alternative_splice_candidates_tsv = with_default(args.alternative_splice_candidates_tsv, args.out_prefix, ".alternative_splice_candidates.tsv")
    alternative_splice_genes_tsv = with_default(args.alternative_splice_genes_tsv, args.out_prefix, ".alternative_splice_genes.tsv")
    alternative_splice_run_summary_tsv = with_default(args.alternative_splice_run_summary_tsv, args.out_prefix, ".alternative_splice_run_summary.tsv")
    alternative_splice_candidate_audit_gff3 = with_default(args.alternative_splice_candidate_audit_gff3, args.out_prefix, ".alternative_splice_candidate_audit.gff3")
    protein_evidence_summary_tsv = with_default(args.protein_evidence_summary_tsv, args.out_prefix, ".protein_evidence_summary.tsv")
    protein_evidence_raw_to_collapsed_tsv = with_default(args.protein_evidence_raw_to_collapsed_tsv, args.out_prefix, ".protein_evidence_raw_to_collapsed.tsv")
    protein_evidence_collapsed_alignments_tsv = with_default(args.protein_evidence_collapsed_alignments_tsv, args.out_prefix, ".protein_evidence_collapsed_alignments.tsv")
    protein_evidence_superloci_tsv = with_default(args.protein_evidence_superloci_tsv, args.out_prefix, ".protein_evidence_superloci.tsv")
    protein_evidence_subloci_tsv = with_default(args.protein_evidence_subloci_tsv, args.out_prefix, ".protein_evidence_subloci.tsv")
    protein_evidence_sublocus_members_tsv = with_default(args.protein_evidence_sublocus_members_tsv, args.out_prefix, ".protein_evidence_sublocus_members.tsv")
    protein_evidence_model_support_tsv = with_default(args.protein_evidence_model_support_tsv, args.out_prefix, ".protein_evidence_model_support.tsv")
    transcript_evidence_summary_tsv = with_default(args.transcript_evidence_summary_tsv, args.out_prefix, ".transcript_evidence_summary.tsv")
    transcript_evidence_raw_to_collapsed_tsv = with_default(args.transcript_evidence_raw_to_collapsed_tsv, args.out_prefix, ".transcript_evidence_raw_to_collapsed.tsv")
    transcript_evidence_collapsed_transcripts_tsv = with_default(args.transcript_evidence_collapsed_transcripts_tsv, args.out_prefix, ".transcript_evidence_collapsed_transcripts.tsv")
    transcript_evidence_superloci_tsv = with_default(args.transcript_evidence_superloci_tsv, args.out_prefix, ".transcript_evidence_superloci.tsv")
    transcript_evidence_subloci_tsv = with_default(args.transcript_evidence_subloci_tsv, args.out_prefix, ".transcript_evidence_subloci.tsv")
    transcript_evidence_sublocus_members_tsv = with_default(args.transcript_evidence_sublocus_members_tsv, args.out_prefix, ".transcript_evidence_sublocus_members.tsv")
    transcript_evidence_model_support_tsv = with_default(args.transcript_evidence_model_support_tsv, args.out_prefix, ".transcript_evidence_model_support.tsv")
    scored_gff3_prune_singleton_sources = parse_source_filter_values(
        args.scored_gff3_prune_singleton_sources,
        source_names,
        "--scored-gff3-prune-singleton-sources",
    )
    final_prune_singleton_subloci_sources = parse_source_filter_values(
        args.final_prune_singleton_subloci_source,
        source_names,
        "--final-prune-singleton-subloci-source",
    )
    shadow_prune_sources = parse_source_filter_values(
        args.shadow_prune_source,
        source_names,
        "--shadow-prune-source",
    )

    if args.model_validation_action == "stop" and model_failures:
        write_models_tsv(coding_models, models_tsv)
        write_excluded_tsv(excluded_models, excluded_tsv)
        write_other_features_gff3(parsed_inputs, other_gff3)
        write_input_summary_tsv(parsed_inputs, input_summary_tsv)
        source_auto_weights = calculate_source_auto_weights(
            coding_models,
            source_names=source_names,
            source_prior_mode=args.source_prior_mode,
            manual_source_weights=manual_source_weights,
            source_prior_scale=args.source_prior_scale,
            auto_source_weight_min=args.auto_source_weight_min,
            auto_source_weight_max=args.auto_source_weight_max,
        )
        write_source_auto_weights_tsv(source_auto_weights, source_auto_weights_tsv)
        write_validation_failures_tsv(validation_failures, validation_failures_tsv)
        print("ERROR: Model-level validation failed before superlocus construction.")
        print(f"Review validation failures TSV: {validation_failures_tsv}")
        print("Rerun with relaxed --max-model-* thresholds, or use --model-validation-action exclude/warn.")
        return 2

    if not coding_models:
        # Early biological filters are intentionally capable of emptying the
        # candidate universe. Preserve their compact reports and enabled
        # ledgers before terminating at the existing no-candidate boundary.
        assert definite_internal_stop_result is not None
        assert extreme_short_cds_result is not None
        if output_enabled("definite_internal_stop_exclusion_summary_tsv"):
            write_definite_internal_stop_summary_tsv(
                definite_internal_stop_result,
                definite_internal_stop_exclusion_summary_tsv,
            )
        if output_enabled("definite_internal_stop_exclusions_tsv"):
            write_definite_internal_stop_exclusions_tsv(
                definite_internal_stop_result,
                definite_internal_stop_exclusions_tsv,
            )
        if output_enabled("extreme_short_cds_exclusion_summary_tsv"):
            write_extreme_short_cds_summary_tsv(
                extreme_short_cds_result,
                extreme_short_cds_exclusion_summary_tsv,
            )
        if output_enabled("extreme_short_cds_exclusions_tsv"):
            write_extreme_short_cds_exclusions_tsv(
                extreme_short_cds_result,
                extreme_short_cds_exclusions_tsv,
            )
        raise SystemExit(
            "No valid coding models remained after early biological exclusion "
            "and model-level validation."
        )
    timer.mark("Model validation")

    progress.log(f"Valid coding models after validation={len(coding_models):,}; excluded_models={len(excluded_models):,}")
    progress.log("Building superloci")
    superloci = build_superloci(coding_models)
    progress.log(f"Built superloci: {len(superloci):,}")
    timer.mark("Superlocus clustering")
    superlocus_failures = collect_superlocus_validation_failures(
        superloci,
        action=args.superlocus_validation_action,
        max_superlocus_span_bp=args.max_superlocus_span_bp,
        max_superlocus_gene_count=args.max_superlocus_gene_count,
        max_genes_per_source_in_superlocus=args.max_genes_per_source_in_superlocus,
    )
    validation_failures.extend(superlocus_failures)

    if args.superlocus_validation_action == "stop" and superlocus_failures:
        # Stop immediately after superlocus construction. Write only the
        # diagnostics that are valid at this stage; do not continue into
        # sublocus construction, evidence processing, scoring, or final outputs.
        write_models_tsv(coding_models, models_tsv)
        write_excluded_tsv(excluded_models, excluded_tsv)
        write_other_features_gff3(parsed_inputs, other_gff3)
        write_input_summary_tsv(parsed_inputs, input_summary_tsv)
        source_auto_weights = calculate_source_auto_weights(
            coding_models,
            source_names=source_names,
            source_prior_mode=args.source_prior_mode,
            manual_source_weights=manual_source_weights,
            source_prior_scale=args.source_prior_scale,
            auto_source_weight_min=args.auto_source_weight_min,
            auto_source_weight_max=args.auto_source_weight_max,
        )
        write_source_auto_weights_tsv(source_auto_weights, source_auto_weights_tsv)
        write_superloci_tsv(superloci, [], superloci_tsv, validation_failures)
        write_members_tsv(superloci, members_tsv)
        write_validation_failures_tsv(validation_failures, validation_failures_tsv)
        print("ERROR: Superlocus validation failed after clustering.")
        print(f"Review validation failures TSV: {validation_failures_tsv}")
        print("Rerun with relaxed --max-superlocus-* thresholds, or use --superlocus-validation-action warn.")
        return 2

    progress.log("Building conservative subloci")
    subloci, relationship_edges = build_subloci(
        superloci,
        single_exon_min_reciprocal_overlap=args.single_exon_min_reciprocal_overlap,
    )
    progress.log(f"Built subloci: {len(subloci):,}; relationship_edges={len(relationship_edges):,}")
    timer.mark("Sublocus clustering")

    protein_evidence_raw_alignments = []
    protein_evidence_superloci = []
    protein_evidence_collapsed_alignments = []
    protein_evidence_subloci = []
    protein_evidence_model_support_rows = []
    protein_evidence_summary_rows = []
    protein_evidence_source_labels = []
    protein_raw_mapping_written_early = False
    transcript_evidence_raw_transcripts = []
    transcript_utr_evidence_catalogue = []
    transcript_evidence_superloci = []
    transcript_evidence_collapsed_transcripts = []
    transcript_evidence_subloci = []
    transcript_evidence_summary_rows = []
    transcript_evidence_model_support_rows = []
    transcript_evidence_source_labels = []
    sequence_integrity_summary: list[dict[str, object]] = []
    sequence_integrity_log_messages: list[str] = []
    transcript_raw_mapping_written_early = False
    protein_summary: dict[str, object] = {}
    transcript_summary: dict[str, object] = {}
    transcript_model_summary: dict[str, object] = {}
    protein_overlap_model_fraction = 0.0
    protein_supported_model_fraction = 0.0
    transcript_scoring_eligible_model_fraction = 0.0
    transcript_supported_model_fraction = 0.0
    if protein_evidence_paths:
        protein_evidence_source_labels = protein_source_labels(protein_evidence_paths, args.evidence_protein_source_name)
        progress.log(
            f"Parsing protein evidence: files={len(protein_evidence_paths):,}; "
            f"sources={','.join(protein_evidence_source_labels)}"
        )
        protein_evidence_raw_alignments = parse_protein_evidence(
            protein_evidence_paths,
            protein_evidence_source_labels,
            strict=args.strict,
        )
        progress.log(f"Parsed protein evidence alignments: raw_alignments={len(protein_evidence_raw_alignments):,}")
        timer.mark("Protein evidence parsing")
        progress.log("Building and collapsing protein evidence loci")
        protein_evidence_superloci = build_protein_superloci(protein_evidence_raw_alignments)
        protein_evidence_collapsed_alignments = collapse_protein_alignments(protein_evidence_superloci)
        protein_evidence_subloci = build_protein_subloci(
            protein_evidence_superloci,
            single_exon_min_reciprocal_overlap=args.single_exon_min_reciprocal_overlap,
        )
        timer.mark("Protein evidence construction")
        progress.log("Associating protein evidence with candidate models")
        protein_evidence_model_support_rows = calculate_protein_model_support(
            coding_models,
            protein_evidence_superloci,
            protein_evidence_source_labels,
            terminal_cds_3prime_wobble_bp=protein_evidence_wobble_bp,
            support_score_threshold=args.protein_evidence_support_score_threshold,
        )
        timer.mark("Protein model association/scoring")
        protein_evidence_summary_rows = protein_summary_tsv_rows(
            protein_paths=protein_evidence_paths,
            protein_source_labels=protein_evidence_source_labels,
            raw_alignments=protein_evidence_raw_alignments,
            superloci=protein_evidence_superloci,
            collapsed_alignments=protein_evidence_collapsed_alignments,
            subloci=protein_evidence_subloci,
            model_support_rows=protein_evidence_model_support_rows,
            raw_to_collapsed_mapping_requested=output_enabled("protein_evidence_raw_to_collapsed_tsv"),
        )
        protein_summary = protein_summary_counts(
            protein_paths=protein_evidence_paths,
            protein_source_labels=protein_evidence_source_labels,
            raw_alignments=protein_evidence_raw_alignments,
            superloci=protein_evidence_superloci,
            collapsed_alignments=protein_evidence_collapsed_alignments,
            subloci=protein_evidence_subloci,
            model_support_rows=protein_evidence_model_support_rows,
            raw_to_collapsed_mapping_requested=output_enabled("protein_evidence_raw_to_collapsed_tsv"),
        )
        protein_model_rows_total = int(protein_summary.get("model_support_row_count", 0) or 0)
        if protein_model_rows_total > 0:
            protein_overlap_model_fraction = (
                int(protein_summary.get("models_with_overlapping_protein_evidence", 0) or 0) / protein_model_rows_total
            )
            protein_supported_model_fraction = (
                (
                    int(protein_summary.get("models_with_strong_protein_support", 0) or 0)
                    + int(protein_summary.get("models_with_partial_protein_support", 0) or 0)
                )
                / protein_model_rows_total
            )
        for line in protein_summary_log_lines(protein_summary):
            progress.log(line)
        if output_enabled("protein_evidence_raw_to_collapsed_tsv"):
            # Write the optional raw mapping before releasing raw alignment
            # objects. The writer is streaming, so it does not create another
            # full in-memory row table.
            write_protein_raw_to_collapsed_tsv(
                protein_evidence_raw_alignments,
                protein_evidence_raw_to_collapsed_tsv,
            )
            protein_raw_mapping_written_early = True
        released_raw = release_protein_raw_alignments(
            protein_evidence_raw_alignments,
            protein_evidence_superloci,
        )
        gc.collect()
        progress.log(
            f"Released raw protein evidence after collapse: records={released_raw:,}",
            level="detailed",
        )
        timer.mark("Protein evidence summaries")
    if transcript_evidence_paths:
        transcript_evidence_source_labels = transcript_source_labels(transcript_evidence_paths, args.evidence_transcript_source_name)
        progress.log(
            f"Parsing transcript evidence: files={len(transcript_evidence_paths):,}; "
            f"sources={','.join(transcript_evidence_source_labels)}"
        )
        transcript_evidence_raw_transcripts = parse_transcript_evidence(
            transcript_evidence_paths,
            transcript_evidence_source_labels,
            strict=args.strict,
        )
        progress.log(f"Parsed transcript evidence records: raw_transcripts={len(transcript_evidence_raw_transcripts):,}")
        timer.mark("Transcript evidence parsing")
        progress.log("Building and collapsing transcript evidence loci")
        transcript_evidence_superloci = build_transcript_superloci(transcript_evidence_raw_transcripts)
        transcript_evidence_collapsed_transcripts = collapse_transcript_evidence(transcript_evidence_superloci)
        transcript_evidence_subloci = build_transcript_subloci(
            transcript_evidence_superloci,
            single_exon_min_reciprocal_overlap=args.single_exon_min_reciprocal_overlap,
        )
        timer.mark("Transcript evidence construction")
        progress.log("Associating transcript evidence with candidate models")
        transcript_evidence_model_support_rows = calculate_transcript_model_support(
            coding_models,
            transcript_evidence_superloci,
            transcript_evidence_source_labels,
        )
        timer.mark("Transcript model association/scoring")
        annotate_transcript_cdna_context_concern_rows(
            transcript_evidence_model_support_rows,
            neutral_fraction=args.transcript_cdna_context_neutral_fraction,
            severe_fraction=args.transcript_cdna_context_severe_fraction,
        )
        timer.mark("Transcript cDNA-context annotation")
        transcript_model_summary = transcript_model_support_summary_counts(transcript_evidence_model_support_rows)
        transcript_model_rows_total = int(transcript_model_summary.get("model_support_row_count", 0) or 0)
        if transcript_model_rows_total > 0:
            transcript_scoring_eligible_model_fraction = (
                int(transcript_model_summary.get("models_with_scoring_eligible_transcript_evidence", 0) or 0)
                / transcript_model_rows_total
            )
            transcript_supported_model_fraction = (
                int(transcript_model_summary.get("models_with_partial_or_strong_transcript_support", 0) or 0)
                / transcript_model_rows_total
            )
        for line in transcript_model_support_log_lines(transcript_model_summary):
            progress.log(line)
        transcript_evidence_summary_rows = transcript_summary_tsv_rows(
            transcript_paths=transcript_evidence_paths,
            transcript_source_labels=transcript_evidence_source_labels,
            raw_transcripts=transcript_evidence_raw_transcripts,
            superloci=transcript_evidence_superloci,
            collapsed_transcripts=transcript_evidence_collapsed_transcripts,
            subloci=transcript_evidence_subloci,
            raw_to_collapsed_mapping_requested=output_enabled("transcript_evidence_raw_to_collapsed_tsv"),
            model_support_summary=transcript_model_summary,
        )
        transcript_summary = transcript_summary_counts(
            transcript_paths=transcript_evidence_paths,
            transcript_source_labels=transcript_evidence_source_labels,
            raw_transcripts=transcript_evidence_raw_transcripts,
            superloci=transcript_evidence_superloci,
            collapsed_transcripts=transcript_evidence_collapsed_transcripts,
            subloci=transcript_evidence_subloci,
            raw_to_collapsed_mapping_requested=output_enabled("transcript_evidence_raw_to_collapsed_tsv"),
        )
        for line in transcript_summary_log_lines(transcript_summary):
            progress.log(line)
        if output_enabled("transcript_evidence_raw_to_collapsed_tsv"):
            write_transcript_raw_to_collapsed_tsv(
                transcript_evidence_raw_transcripts,
                transcript_evidence_raw_to_collapsed_tsv,
            )
            transcript_raw_mapping_written_early = True
        if args.add_utrs:
            progress.log("Building compact endpoint-paired transcript UTR catalogue")
            transcript_utr_evidence_catalogue = build_transcript_utr_catalogue(
                transcript_evidence_raw_transcripts
            )
            progress.log(
                "Built transcript UTR catalogue: "
                f"records={len(transcript_utr_evidence_catalogue):,}"
            )
            timer.mark("Transcript UTR catalogue")
        released_raw = release_transcript_raw_transcripts(
            transcript_evidence_raw_transcripts,
            transcript_evidence_superloci,
        )
        gc.collect()
        progress.log(
            f"Released raw transcript evidence after collapse: records={released_raw:,}",
            level="detailed",
        )
        timer.mark("Transcript evidence summaries")
    if genome_fasta_path is not None:
        assert indexed_fasta is not None
        progress.log("Calculating genome-derived sequence-integrity metrics by superlocus")
        with indexed_fasta as fasta_reader:
            sequence_integrity_summary, sequence_integrity_log_messages = (
                calculate_and_attach_sequence_integrity(
                    superloci,
                    fasta_reader,
                    strict=args.strict,
                    genetic_code=args.genetic_code,
                    internal_stop_penalty=args.sequence_penalty_internal_stop,
                    internal_stop_count_cap=args.sequence_internal_stop_count_cap,
                    phase_inconsistent_penalty=args.sequence_penalty_phase_inconsistent,
                    missing_stop_penalty=args.sequence_penalty_missing_stop,
                    missing_start_penalty=args.sequence_penalty_missing_start,
                    noncanonical_splice_penalty=args.sequence_penalty_noncanonical_splice,
                    penalty_active=(
                        not args.no_sequence_integrity_scoring
                        and args.selection_penalty_sequence_integrity > 0
                    ),
                    progress=progress,
                )
            )
        for line in sequence_integrity_log_messages:
            progress.log(line)
        timer.mark("Genome sequence integrity metrics")

    long_cds_intron_cutoff, long_cds_intron_cutoff_source = resolve_long_cds_intron_cutoff(
        coding_models,
        user_long_cds_intron_bp=args.long_cds_intron_bp,
        percentile_value=args.long_cds_intron_percentile,
    )
    min_terminal_cds_exon_bp = (
        args.min_terminal_cds_exon_bp
        if args.min_terminal_cds_exon_bp is not None
        else args.min_cds_exon_bp
        if args.min_cds_exon_bp is not None
        else 3
    )
    min_internal_cds_exon_bp = (
        args.min_internal_cds_exon_bp
        if args.min_internal_cds_exon_bp is not None
        else args.min_cds_exon_bp
        if args.min_cds_exon_bp is not None
        else 6
    )
    for label, intrinsic_threshold, soft_threshold in (
        ("5-prime terminal", min_terminal_cds_exon_bp, args.unsupported_short_5prime_terminal_cds_exon_bp),
        ("3-prime terminal", min_terminal_cds_exon_bp, args.unsupported_short_3prime_terminal_cds_exon_bp),
        ("internal", min_internal_cds_exon_bp, args.unsupported_short_internal_cds_exon_bp),
    ):
        if soft_threshold <= intrinsic_threshold:
            progress.log(
                f"Unsupported-short {label} CDS exon penalty disabled: soft threshold {soft_threshold} bp is not greater "
                f"than active intrinsic threshold {intrinsic_threshold} bp; shorter exons remain handled by the intrinsic penalty."
            )
    progress.log("Calculating intrinsic model metrics")
    model_metrics = calculate_all_model_intrinsic_metrics(
        coding_models,
        min_cds_length_bp=args.min_cds_length_bp,
        min_terminal_cds_exon_bp=min_terminal_cds_exon_bp,
        min_internal_cds_exon_bp=min_internal_cds_exon_bp,
        min_cds_intron_bp=args.min_cds_intron_bp,
        long_cds_intron_bp=long_cds_intron_cutoff,
    )
    timer.mark("Intrinsic model metrics")
    progress.log("Calculating source prior / auto weights")
    source_auto_weights = calculate_source_auto_weights(
        coding_models,
        source_names=source_names,
        source_prior_mode=args.source_prior_mode,
        manual_source_weights=manual_source_weights,
        source_prior_scale=args.source_prior_scale,
        auto_source_weight_min=args.auto_source_weight_min,
        auto_source_weight_max=args.auto_source_weight_max,
    )
    timer.mark("Source priors")
    single_exon_support_min_reciprocal_overlap = (
        args.single_exon_support_min_reciprocal_overlap
        if args.single_exon_support_min_reciprocal_overlap is not None
        else args.single_exon_min_reciprocal_overlap
    )
    intrinsic_penalty_weights = {
        "short_cds": args.intrinsic_penalty_short_cds,
        "cds_length_not_multiple_of_3": args.intrinsic_penalty_cds_length_not_multiple_of_3,
        "short_terminal_cds_exon": args.intrinsic_penalty_short_terminal_cds_exon,
        "short_internal_cds_exon": args.intrinsic_penalty_short_internal_cds_exon,
        "short_cds_intron": args.intrinsic_penalty_short_cds_intron,
        "long_cds_intron": args.intrinsic_penalty_long_cds_intron,
        "contained_model": args.intrinsic_penalty_contained_model,
        "unsupported_cds_intron_contained_in_other_cds_exon": args.intrinsic_penalty_unsupported_contained_cds_intron,
    }
    progress.log("Calculating structural support, consensus similarity, completeness, exact-structure support")
    structural_corroboration_rows: Optional[list[dict[str, object]]] = (
        [] if output_enabled("structural_corroboration_tsv") else None
    )
    structural_corroboration_calculation_enabled = (
        structural_corroboration_rows is not None
        or args.selection_weight_structural_corroboration > 0
        or args.strict_r1_backfill_replacement
    )
    structural_feature_rows, model_structural_rows = calculate_structural_support(
        subloci,
        terminal_cds_3prime_wobble_bp=args.terminal_cds_3prime_wobble_bp,
        single_exon_support_min_reciprocal_overlap=single_exon_support_min_reciprocal_overlap,
        min_cds_length_bp=args.min_cds_length_bp,
        min_terminal_cds_exon_bp=min_terminal_cds_exon_bp,
        min_internal_cds_exon_bp=min_internal_cds_exon_bp,
        unsupported_short_5prime_terminal_cds_exon_bp=args.unsupported_short_5prime_terminal_cds_exon_bp,
        unsupported_short_3prime_terminal_cds_exon_bp=args.unsupported_short_3prime_terminal_cds_exon_bp,
        unsupported_short_internal_cds_exon_bp=args.unsupported_short_internal_cds_exon_bp,
        unsupported_short_cds_exon_component_weights={
            "5prime_terminal": args.unsupported_short_5prime_terminal_weight,
            "3prime_terminal": args.unsupported_short_3prime_terminal_weight,
            "internal": args.unsupported_short_internal_weight,
        },
        min_cds_intron_bp=args.min_cds_intron_bp,
        long_cds_intron_bp=long_cds_intron_cutoff,
        intrinsic_penalty_weights=intrinsic_penalty_weights,
        run_source_count=len(source_names),
        structural_corroboration_enabled=structural_corroboration_calculation_enabled,
        structural_corroboration_rows=structural_corroboration_rows,
        progress=progress,
    )
    progress.log(f"Calculated structural rows: models={len(model_structural_rows):,}; feature_rows={len(structural_feature_rows):,}")
    progress.log("Calculating candidate structural-coherence regions")
    structural_coherence_region_rows = annotate_structural_coherence(
        superloci,
        model_structural_rows,
        penalty_weight=args.selection_penalty_structural_coherence,
    )
    progress.log(
        f"Calculated structural-coherence regions: {len(structural_coherence_region_rows):,}"
    )
    timer.mark("Model scoring / structural support")
    source_prior_scores = {item.source_name: item.source_prior_score for item in source_auto_weights}
    sequence_integrity_active_for_selection = (
        genome_fasta_path is not None
        and not args.no_sequence_integrity_scoring
        and args.selection_penalty_sequence_integrity > 0
    )
    effective_sequence_integrity_penalty_weight = (
        args.selection_penalty_sequence_integrity if sequence_integrity_active_for_selection else 0.0
    )
    model_selection_weights = {
        "structural_support": args.selection_weight_structural_support,
        "consensus_similarity": args.selection_weight_consensus_similarity,
        "pairwise_similarity": args.selection_weight_pairwise_similarity,
        "completeness": args.selection_weight_completeness,
        "source_prior": args.selection_weight_source_prior,
        "terminal_extension": args.selection_weight_terminal_extension,
        "exact_structure_support": args.selection_weight_exact_structure_support,
        "structural_corroboration": args.selection_weight_structural_corroboration,
        "protein_evidence": args.selection_weight_protein_evidence,
        "transcript_evidence": args.selection_weight_transcript_evidence,
        "transcript_cdna_context_penalty": args.selection_penalty_transcript_cdna_context,
        "unsupported_short_cds_exon_penalty": args.selection_penalty_unsupported_short_cds_exon,
        "structural_coherence_penalty": args.selection_penalty_structural_coherence,
        "sequence_integrity_penalty": effective_sequence_integrity_penalty_weight,
        "intrinsic_penalty": args.selection_penalty_intrinsic,
        "partial_structure_penalty": args.selection_penalty_partial_structure,
        "internal_extra_penalty": args.selection_penalty_internal_extra,
        "terminal_extra_penalty": args.selection_penalty_terminal_extra,
    }
    protein_evidence_rows_by_transcript = {
        str(row.get("output_transcript_id")): row for row in protein_evidence_model_support_rows
    }
    transcript_evidence_rows_by_transcript = {
        str(row.get("output_transcript_id")): row for row in transcript_evidence_model_support_rows
    }
    protein_evidence_active_for_selection = bool(protein_evidence_paths) and args.selection_weight_protein_evidence > 0
    transcript_evidence_active_for_selection = bool(transcript_evidence_paths) and args.selection_weight_transcript_evidence > 0
    evidence_weight_text = []
    if sequence_integrity_active_for_selection:
        evidence_weight_text.append(
            f"sequence_integrity_penalty weight={effective_sequence_integrity_penalty_weight:.3f}"
        )
    if protein_evidence_active_for_selection:
        evidence_weight_text.append(f"protein_evidence_score weight={args.selection_weight_protein_evidence:.3f}")
    if transcript_evidence_active_for_selection:
        evidence_weight_text.append(
            f"transcript_evidence_score weight={args.selection_weight_transcript_evidence:.3f}; "
            f"cDNA-context max_penalty={args.selection_penalty_transcript_cdna_context:.3f}"
        )
    progress.log(
        "Applying model selection scores and ranks"
        + (" with " + "; ".join(evidence_weight_text) if evidence_weight_text else "")
    )
    apply_model_selection_scores(
        model_structural_rows,
        source_prior_scores=source_prior_scores,
        selection_weights=model_selection_weights,
        tie_epsilon=args.model_selection_tie_epsilon,
        show_zero_weight_components=args.show_zero_weight_score_components,
        low_source_conflict_dampen_agreement_scores=args.low_source_conflict_dampen_agreement_scores,
        protein_evidence_rows_by_transcript=protein_evidence_rows_by_transcript,
        protein_evidence_active=protein_evidence_active_for_selection,
        transcript_evidence_rows_by_transcript=transcript_evidence_rows_by_transcript,
        transcript_evidence_active=transcript_evidence_active_for_selection,
        transcript_cdna_context_neutral_fraction=args.transcript_cdna_context_neutral_fraction,
        transcript_cdna_context_severe_fraction=args.transcript_cdna_context_severe_fraction,
        transcript_cdna_context_penalty_weight=args.selection_penalty_transcript_cdna_context,
        structural_corroboration_alpha=args.structural_corroboration_alpha,
    )
    model_structural_rows_by_transcript = {
        str(row.get("output_transcript_id")): row for row in model_structural_rows
    }
    if structural_corroboration_rows is not None:
        for diagnostic_row in structural_corroboration_rows:
            scoring_row = model_structural_rows_by_transcript.get(
                str(diagnostic_row.get("candidate_model_id")), {}
            )
            for field in (
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
            ):
                diagnostic_row[field] = scoring_row.get(field)
    timer.mark("Model selection ranking")

    progress.log("Writing scored candidate outputs")
    if output_enabled("combined_gff3"):
        write_combined_gff3(
            superloci,
            subloci,
            combined_gff3,
            validation_failures,
            model_structural_rows_by_transcript=model_structural_rows_by_transcript,
            include_model_selection_components=args.gff3_include_model_selection_score_components,
        )
        note_output("combined_gff3", "scored candidate combined GFF3", combined_gff3)
    if args.prune_singleton_subloci_from_scored_gff3 and output_enabled("scored_pruned_gff3"):
        write_combined_gff3(
            superloci,
            subloci,
            scored_pruned_gff3,
            validation_failures,
            model_structural_rows_by_transcript=model_structural_rows_by_transcript,
            include_model_selection_components=args.gff3_include_model_selection_score_components,
            prune_single_model_subloci=True,
            scored_gff3_prune_singleton_sources=scored_gff3_prune_singleton_sources,
            keep_empty_superloci=args.scored_gff3_keep_empty_superloci,
        )
        note_output("scored_pruned_gff3", "scored-candidate pruned-singletons GFF3", scored_pruned_gff3)
    if output_enabled("superloci_tsv"):
        write_superloci_tsv(superloci, subloci, superloci_tsv, validation_failures)
        note_output("superloci_tsv", "superloci TSV", superloci_tsv)
    if output_enabled("superlocus_members_tsv"):
        write_members_tsv(superloci, members_tsv)
        note_output("superlocus_members_tsv", "superlocus members TSV", members_tsv)
    if output_enabled("subloci_tsv"):
        write_subloci_tsv(subloci, subloci_tsv)
        note_output("subloci_tsv", "subloci TSV", subloci_tsv)
    if output_enabled("sublocus_members_tsv"):
        write_sublocus_members_tsv(subloci, sublocus_members_tsv)
        note_output("sublocus_members_tsv", "sublocus members TSV", sublocus_members_tsv)
    if output_enabled("relationship_edges_tsv"):
        write_relationship_edges_tsv(relationship_edges, relationship_edges_tsv)
        note_output("relationship_edges_tsv", "relationship edges TSV", relationship_edges_tsv)
    if output_enabled("model_intrinsic_metrics_tsv"):
        write_model_intrinsic_metrics_tsv(model_metrics, model_metrics_tsv)
        note_output("model_intrinsic_metrics_tsv", "model intrinsic metrics TSV", model_metrics_tsv)
    if genome_fasta_path is not None and output_enabled("sequence_integrity_summary_tsv"):
        write_sequence_integrity_summary_tsv(sequence_integrity_summary, sequence_integrity_summary_tsv)
        note_output("sequence_integrity_summary_tsv", "sequence integrity summary TSV", sequence_integrity_summary_tsv)
    if genome_fasta_path is not None and output_enabled("sequence_integrity_model_metrics_tsv"):
        write_sequence_integrity_model_metrics_from_models(
            coding_models, sequence_integrity_model_metrics_tsv
        )
        note_output("sequence_integrity_model_metrics_tsv", "sequence integrity model metrics TSV", sequence_integrity_model_metrics_tsv)
    terminal_report_models = terminal_stop_audit_models or coding_models
    five_prime_report_models = five_prime_boundary_audit_models or coding_models
    phase_report_models = cds_phase_recovery_audit_models or coding_models
    early_execution_status = (
        "performed"
        if genome_fasta_path is not None
        else "prerequisite_unavailable:no_genome_fasta"
    )
    if output_enabled("terminal_stop_canonicalisation_summary_tsv"):
        write_terminal_stop_summary_tsv(
            terminal_report_models,
            terminal_stop_canonicalisation_summary_tsv,
            execution_status=early_execution_status,
        )
        note_output(
            "terminal_stop_canonicalisation_summary_tsv",
            "terminal-stop canonicalisation summary TSV",
            terminal_stop_canonicalisation_summary_tsv,
        )
    if output_enabled("terminal_stop_canonicalisation_actions_tsv"):
        write_terminal_stop_actions_tsv(
            terminal_report_models, terminal_stop_canonicalisation_actions_tsv
        )
        note_output(
            "terminal_stop_canonicalisation_actions_tsv",
            "terminal-stop canonicalisation actions TSV",
            terminal_stop_canonicalisation_actions_tsv,
        )
    if output_enabled("terminal_stop_canonicalisation_audit_tsv"):
        write_terminal_stop_audit_tsv(
            terminal_report_models,
            terminal_stop_canonicalisation_audit_tsv,
        )
        note_output(
            "terminal_stop_canonicalisation_audit_tsv",
            "terminal-stop canonicalisation audit TSV",
            terminal_stop_canonicalisation_audit_tsv,
        )
    if output_enabled("cds_5prime_boundary_canonicalisation_summary_tsv"):
        write_five_prime_boundary_summary_tsv(
            five_prime_report_models,
            five_prime_boundary_canonicalisation_summary_tsv,
            execution_status=early_execution_status,
        )
        note_output(
            "cds_5prime_boundary_canonicalisation_summary_tsv",
            "5-prime CDS-boundary canonicalisation summary TSV",
            five_prime_boundary_canonicalisation_summary_tsv,
        )
    if output_enabled("cds_5prime_boundary_canonicalisation_actions_tsv"):
        write_five_prime_boundary_actions_tsv(
            five_prime_report_models,
            five_prime_boundary_canonicalisation_actions_tsv,
        )
        note_output(
            "cds_5prime_boundary_canonicalisation_actions_tsv",
            "5-prime CDS-boundary canonicalisation actions TSV",
            five_prime_boundary_canonicalisation_actions_tsv,
        )
    if output_enabled("cds_5prime_boundary_canonicalisation_audit_tsv"):
        write_five_prime_boundary_audit_tsv(
            five_prime_report_models,
            five_prime_boundary_canonicalisation_audit_tsv,
        )
        note_output(
            "cds_5prime_boundary_canonicalisation_audit_tsv",
            "5-prime CDS-boundary canonicalisation audit TSV",
            five_prime_boundary_canonicalisation_audit_tsv,
        )
    if output_enabled("cds_phase_recovery_summary_tsv"):
        write_phase_recovery_summary_tsv(
            phase_report_models,
            cds_phase_recovery_summary_tsv,
            execution_status=early_execution_status,
        )
        note_output(
            "cds_phase_recovery_summary_tsv",
            "CDS phase-recovery summary TSV",
            cds_phase_recovery_summary_tsv,
        )
    if output_enabled("cds_phase_recovery_actions_tsv"):
        write_phase_recovery_actions_tsv(
            phase_report_models, cds_phase_recovery_actions_tsv
        )
        note_output(
            "cds_phase_recovery_actions_tsv",
            "CDS phase-recovery actions TSV",
            cds_phase_recovery_actions_tsv,
        )
    if output_enabled("cds_phase_recovery_audit_tsv"):
        write_phase_recovery_audit_tsv(
            phase_report_models,
            cds_phase_recovery_audit_tsv,
        )
        note_output(
            "cds_phase_recovery_audit_tsv",
            "CDS phase-recovery audit TSV",
            cds_phase_recovery_audit_tsv,
        )
    assert definite_internal_stop_result is not None
    if output_enabled("definite_internal_stop_exclusion_summary_tsv"):
        write_definite_internal_stop_summary_tsv(
            definite_internal_stop_result,
            definite_internal_stop_exclusion_summary_tsv,
        )
        note_output(
            "definite_internal_stop_exclusion_summary_tsv",
            "definite internal-stop exclusion summary TSV",
            definite_internal_stop_exclusion_summary_tsv,
        )
    if output_enabled("definite_internal_stop_exclusions_tsv"):
        write_definite_internal_stop_exclusions_tsv(
            definite_internal_stop_result,
            definite_internal_stop_exclusions_tsv,
        )
        note_output(
            "definite_internal_stop_exclusions_tsv",
            "definite internal-stop exclusions TSV",
            definite_internal_stop_exclusions_tsv,
        )
    assert extreme_short_cds_result is not None
    if output_enabled("extreme_short_cds_exclusion_summary_tsv"):
        write_extreme_short_cds_summary_tsv(
            extreme_short_cds_result,
            extreme_short_cds_exclusion_summary_tsv,
        )
        note_output(
            "extreme_short_cds_exclusion_summary_tsv",
            "extreme-short CDS exclusion summary TSV",
            extreme_short_cds_exclusion_summary_tsv,
        )
    if output_enabled("extreme_short_cds_exclusions_tsv"):
        write_extreme_short_cds_exclusions_tsv(
            extreme_short_cds_result,
            extreme_short_cds_exclusions_tsv,
        )
        note_output(
            "extreme_short_cds_exclusions_tsv",
            "extreme-short CDS exclusions TSV",
            extreme_short_cds_exclusions_tsv,
        )
    if output_enabled("source_auto_weights_tsv"):
        write_source_auto_weights_tsv(source_auto_weights, source_auto_weights_tsv)
        note_output("source_auto_weights_tsv", "source auto-weights TSV", source_auto_weights_tsv)
    if output_enabled("model_structural_support_tsv"):
        write_model_structural_support_tsv(model_structural_rows, model_structural_support_tsv)
        note_output("model_structural_support_tsv", "model structural support TSV", model_structural_support_tsv)
    if output_enabled("structural_corroboration_tsv"):
        assert structural_corroboration_rows is not None
        write_structural_corroboration_tsv(
            structural_corroboration_rows,
            structural_corroboration_tsv,
        )
        note_output(
            "structural_corroboration_tsv",
            "structural corroboration diagnostic TSV",
            structural_corroboration_tsv,
        )
    if output_enabled("structural_coherence_regions_tsv"):
        write_structural_coherence_regions_tsv(structural_coherence_region_rows, structural_coherence_regions_tsv)
        note_output("structural_coherence_regions_tsv", "structural coherence regions TSV", structural_coherence_regions_tsv)
    if output_enabled("model_selection_scores_tsv"):
        write_model_selection_scores_tsv(model_structural_rows, model_selection_scores_tsv)
        note_output("model_selection_scores_tsv", "model selection scores TSV", model_selection_scores_tsv)
    if output_enabled("model_selection_score_summary_tsv"):
        write_model_selection_score_summary_tsv(model_structural_rows, model_selection_score_summary_tsv)
        note_output("model_selection_score_summary_tsv", "model selection score summary TSV", model_selection_score_summary_tsv)
    if output_enabled("structural_features_tsv"):
        write_structural_features_tsv(structural_feature_rows, structural_features_tsv)
        note_output("structural_features_tsv", "structural features TSV", structural_features_tsv)
    if protein_evidence_paths and output_enabled("protein_evidence_summary_tsv"):
        write_protein_summary_tsv(protein_evidence_summary_rows, protein_evidence_summary_tsv)
        note_output("protein_evidence_summary_tsv", "protein evidence summary TSV", protein_evidence_summary_tsv)
    if protein_evidence_paths and output_enabled("protein_evidence_raw_to_collapsed_tsv"):
        if not protein_raw_mapping_written_early:
            write_protein_raw_to_collapsed_tsv(protein_evidence_raw_alignments, protein_evidence_raw_to_collapsed_tsv)
        note_output("protein_evidence_raw_to_collapsed_tsv", "protein evidence raw-to-collapsed TSV", protein_evidence_raw_to_collapsed_tsv)
    if protein_evidence_paths and output_enabled("protein_evidence_collapsed_alignments_tsv"):
        write_protein_collapsed_alignments_tsv(protein_evidence_collapsed_alignments, protein_evidence_collapsed_alignments_tsv)
        note_output("protein_evidence_collapsed_alignments_tsv", "protein evidence collapsed alignments TSV", protein_evidence_collapsed_alignments_tsv)
    if protein_evidence_paths and output_enabled("protein_evidence_superloci_tsv"):
        write_protein_superloci_tsv(protein_evidence_superloci, protein_evidence_superloci_tsv)
        note_output("protein_evidence_superloci_tsv", "protein evidence superloci TSV", protein_evidence_superloci_tsv)
    if protein_evidence_paths and output_enabled("protein_evidence_subloci_tsv"):
        write_protein_subloci_tsv(protein_evidence_subloci, protein_evidence_subloci_tsv)
        note_output("protein_evidence_subloci_tsv", "protein evidence subloci TSV", protein_evidence_subloci_tsv)
    if protein_evidence_paths and output_enabled("protein_evidence_sublocus_members_tsv"):
        write_protein_sublocus_members_tsv(protein_evidence_subloci, protein_evidence_sublocus_members_tsv)
        note_output("protein_evidence_sublocus_members_tsv", "protein evidence sublocus members TSV", protein_evidence_sublocus_members_tsv)
    if protein_evidence_paths and output_enabled("protein_evidence_model_support_tsv"):
        write_protein_model_support_tsv(protein_evidence_model_support_rows, protein_evidence_model_support_tsv)
        note_output("protein_evidence_model_support_tsv", "protein evidence model support TSV", protein_evidence_model_support_tsv)
    if transcript_evidence_paths and output_enabled("transcript_evidence_summary_tsv"):
        write_transcript_summary_tsv(transcript_evidence_summary_rows, transcript_evidence_summary_tsv)
        note_output("transcript_evidence_summary_tsv", "transcript evidence summary TSV", transcript_evidence_summary_tsv)
    if transcript_evidence_paths and output_enabled("transcript_evidence_raw_to_collapsed_tsv"):
        if not transcript_raw_mapping_written_early:
            write_transcript_raw_to_collapsed_tsv(transcript_evidence_raw_transcripts, transcript_evidence_raw_to_collapsed_tsv)
        note_output("transcript_evidence_raw_to_collapsed_tsv", "transcript evidence raw-to-collapsed TSV", transcript_evidence_raw_to_collapsed_tsv)
    if transcript_evidence_paths and output_enabled("transcript_evidence_collapsed_transcripts_tsv"):
        write_transcript_collapsed_transcripts_tsv(transcript_evidence_collapsed_transcripts, transcript_evidence_collapsed_transcripts_tsv)
        note_output("transcript_evidence_collapsed_transcripts_tsv", "transcript evidence collapsed transcripts TSV", transcript_evidence_collapsed_transcripts_tsv)
    if transcript_evidence_paths and output_enabled("transcript_evidence_superloci_tsv"):
        write_transcript_superloci_tsv(transcript_evidence_superloci, transcript_evidence_superloci_tsv)
        note_output("transcript_evidence_superloci_tsv", "transcript evidence superloci TSV", transcript_evidence_superloci_tsv)
    if transcript_evidence_paths and output_enabled("transcript_evidence_subloci_tsv"):
        write_transcript_subloci_tsv(transcript_evidence_subloci, transcript_evidence_subloci_tsv)
        note_output("transcript_evidence_subloci_tsv", "transcript evidence subloci TSV", transcript_evidence_subloci_tsv)
    if transcript_evidence_paths and output_enabled("transcript_evidence_sublocus_members_tsv"):
        write_transcript_sublocus_members_tsv(transcript_evidence_subloci, transcript_evidence_sublocus_members_tsv)
        note_output("transcript_evidence_sublocus_members_tsv", "transcript evidence sublocus members TSV", transcript_evidence_sublocus_members_tsv)
    if transcript_evidence_paths and output_enabled("transcript_evidence_model_support_tsv"):
        write_transcript_model_support_tsv(transcript_evidence_model_support_rows, transcript_evidence_model_support_tsv)
        note_output("transcript_evidence_model_support_tsv", "transcript evidence model support TSV", transcript_evidence_model_support_tsv)
    if output_enabled("coding_models_tsv"):
        write_models_tsv(coding_models, models_tsv)
        note_output("coding_models_tsv", "coding models TSV", models_tsv)
    if output_enabled("excluded_models_tsv"):
        write_excluded_tsv(excluded_models, excluded_tsv)
        note_output("excluded_models_tsv", "excluded models TSV", excluded_tsv)
    if output_enabled("other_features_gff3"):
        write_other_features_gff3(parsed_inputs, other_gff3)
        note_output("other_features_gff3", "other/no-CDS features GFF3", other_gff3)
    if output_enabled("input_summary_tsv"):
        write_input_summary_tsv(parsed_inputs, input_summary_tsv)
        note_output("input_summary_tsv", "input summary TSV", input_summary_tsv)
    if output_enabled("validation_failures_tsv"):
        write_validation_failures_tsv(validation_failures, validation_failures_tsv)
        note_output("validation_failures_tsv", "validation failures TSV", validation_failures_tsv)
    progress.log("Finished writing scored candidate outputs")
    timer.mark("Writing scored candidate outputs")

    stage4_result: Optional[dict[str, object]] = None
    stage5_result: Optional[dict[str, object]] = None
    final_prune_result: Optional[dict[str, object]] = None
    if args.run_monosublocus_selection:
        progress.log("Starting Monosublocus de-overlap")
        stage4_result = stage4_selection_from_sublocus_winners(
            subloci,
            model_structural_rows_by_transcript,
            min_cds_exon_overlap_fraction=args.monosublocus_min_cds_exon_overlap_fraction,
            progress=progress,
            progress_interval=args.progress_interval,
        )
        sublocus_winners = stage4_result["sublocus_winners"]  # type: ignore[index]
        selection_info_by_transcript = stage4_result["selection_info_by_transcript"]  # type: ignore[index]
        monosublocus_rows = stage4_result["monosublocus_rows"]  # type: ignore[index]
        selected_models = stage4_result["selected_models"]  # type: ignore[index]
        timer.mark("Monosublocus de-overlap")

        progress.log("Writing Monosublocus de-overlap monosublocus outputs")
        if output_enabled("sublocus_winners_tsv"):
            write_sublocus_winners_tsv(
                sublocus_winners,  # type: ignore[arg-type]
                subloci,
                model_structural_rows_by_transcript,
                sublocus_winners_tsv,
            )
            note_output("sublocus_winners_tsv", "sublocus winners TSV", sublocus_winners_tsv)
        if output_enabled("sublocus_winners_gff3"):
            write_model_subset_gff3(
                sublocus_winners,  # type: ignore[arg-type]
                sublocus_winners_gff3,
                model_structural_rows_by_transcript=model_structural_rows_by_transcript,
                subloci=subloci,
                include_model_selection_components=args.gff3_include_model_selection_score_components,
                selection_label="sublocus_winner",
            )
            note_output("sublocus_winners_gff3", "sublocus winners GFF3", sublocus_winners_gff3)
        write_monosublocus_outputs(
            stage4_result,
            monosubloci_tsv=monosubloci_tsv if output_enabled("monosubloci_tsv") else None,
            monosublocus_members_tsv=monosublocus_members_tsv if output_enabled("monosublocus_members_tsv") else None,
            monosublocus_edges_tsv=monosublocus_edges_tsv if output_enabled("monosublocus_edges_tsv") else None,
            selected_models_tsv=selected_models_tsv if output_enabled("selected_models_tsv") else None,
            suppressed_models_tsv=suppressed_monosublocus_models_tsv if output_enabled("suppressed_monosublocus_models_tsv") else None,
        )
        note_output("monosubloci_tsv", "monosubloci TSV", monosubloci_tsv)
        note_output("monosublocus_members_tsv", "monosublocus members TSV", monosublocus_members_tsv)
        note_output("monosublocus_edges_tsv", "monosublocus edges TSV", monosublocus_edges_tsv)
        note_output("selected_models_tsv", "selected models TSV", selected_models_tsv)
        note_output("suppressed_monosublocus_models_tsv", "suppressed monosublocus models TSV", suppressed_monosublocus_models_tsv)
        if output_enabled("monosublocus_review_gff3"):
            write_model_subset_gff3(
                sublocus_winners,  # type: ignore[arg-type]
                monosublocus_review_gff3,
                model_structural_rows_by_transcript=model_structural_rows_by_transcript,
                subloci=subloci,
                selection_info_by_transcript=selection_info_by_transcript,  # type: ignore[arg-type]
                monosublocus_rows=monosublocus_rows,  # type: ignore[arg-type]
                include_model_selection_components=args.gff3_include_model_selection_score_components,
                include_monosublocus_features=True,
                selection_label="sublocus_winner",
            )
            note_output("monosublocus_review_gff3", "monosublocus review GFF3", monosublocus_review_gff3)
        if output_enabled("selected_models_gff3"):
            write_model_subset_gff3(
                selected_models,  # type: ignore[arg-type]
                selected_models_gff3,
                model_structural_rows_by_transcript=model_structural_rows_by_transcript,
                subloci=subloci,
                selection_info_by_transcript=selection_info_by_transcript,  # type: ignore[arg-type]
                monosublocus_rows=monosublocus_rows,  # type: ignore[arg-type]
                include_model_selection_components=args.gff3_include_model_selection_score_components,
                include_monosublocus_features=True,
                selection_label="final_selected_model",
            )
            note_output("selected_models_gff3", "selected models GFF3", selected_models_gff3)

        progress.log("Finished writing Monosublocus de-overlap monosublocus outputs")
        timer.mark("Writing monosublocus outputs")

        if args.run_backfill_selection:
            progress.log("Starting Backfill rescue")
            stage5_result = stage5_backfill_selection(
                subloci,
                model_structural_rows_by_transcript,
                stage4_result,
                min_cds_exon_overlap_fraction=args.monosublocus_min_cds_exon_overlap_fraction,
                strict_r1_backfill_replacement=args.strict_r1_backfill_replacement,
                strict_r1_allow_multiexon_blockers=args.strict_r1_allow_multiexon_blockers,
                structural_corroboration_weight=args.selection_weight_structural_corroboration,
                structural_corroboration_alpha=args.structural_corroboration_alpha,
                sequence_penalty_missing_start=args.sequence_penalty_missing_start,
                sequence_penalty_missing_stop=args.sequence_penalty_missing_stop,
                progress=progress,
                progress_interval=args.progress_interval,
            )
            timer.mark("Backfill rescue")
            progress.log("Writing Backfill rescue backfill/final selected outputs")
            write_backfill_outputs(
                stage5_result,
                backfill_candidates_tsv=backfill_candidates_tsv if output_enabled("backfill_candidates_tsv") else None,
                backfilled_models_tsv=backfilled_models_tsv if output_enabled("backfilled_models_tsv") else None,
                backfill_rejected_models_tsv=backfill_rejected_models_tsv if output_enabled("backfill_rejected_models_tsv") else None,
                final_selected_models_tsv=final_selected_models_tsv if output_enabled("final_selected_models_tsv") else None,
            )
            note_output("backfill_candidates_tsv", "backfill candidates TSV", backfill_candidates_tsv)
            note_output("backfilled_models_tsv", "backfilled models TSV", backfilled_models_tsv)
            note_output("backfill_rejected_models_tsv", "backfill rejected models TSV", backfill_rejected_models_tsv)
            note_output("final_selected_models_tsv", "final selected models TSV", final_selected_models_tsv)
            if args.strict_r1_backfill_replacement:
                write_r1_replacement_outputs(
                    stage5_result,
                    replacements_tsv=(
                        r1_backfill_replacements_tsv
                        if output_enabled("r1_backfill_replacements_tsv") else None
                    ),
                    summary_tsv=(
                        r1_backfill_summary_tsv
                        if output_enabled("r1_backfill_summary_tsv") else None
                    ),
                )
                note_output(
                    "r1_backfill_replacements_tsv",
                    "strict R1 backfill replacements TSV",
                    r1_backfill_replacements_tsv,
                )
                note_output(
                    "r1_backfill_summary_tsv",
                    "strict R1 backfill summary TSV",
                    r1_backfill_summary_tsv,
                )
            if output_enabled("final_selected_models_gff3"):
                write_model_subset_gff3(
                    stage5_result["final_models"],  # type: ignore[arg-type]
                    final_selected_models_gff3,
                    model_structural_rows_by_transcript=model_structural_rows_by_transcript,
                    subloci=subloci,
                    selection_info_by_transcript=stage5_result["final_selection_info_by_transcript"],  # type: ignore[arg-type]
                    monosublocus_rows=monosublocus_rows,  # type: ignore[arg-type]
                    include_model_selection_components=args.gff3_include_model_selection_score_components,
                    include_monosublocus_features=True,
                    selection_label="final_selected_model",
                    cds_only_model_structure=args.final_gff3_cds_only_structure,
                )
                note_output("final_selected_models_gff3", "final selected models GFF3", final_selected_models_gff3)
            if final_selected_pruning_requested:
                progress.log("Starting final selected-model pruning")
                def pruning_configuration(level: str, *, effective_args: bool = False) -> dict[str, object]:
                    if effective_args:
                        return {
                            "level": level,
                            "singleton_level": args.final_prune_singleton_subloci_level,
                            "low_score_singletons": args.final_prune_low_score_singletons,
                            "low_score_nonsingletons": args.final_prune_low_score_nonsingletons,
                            "low_score_nonsingleton_scope": args.low_score_nonsingleton_scope,
                            "sequence_integrity_nonsingletons": args.final_prune_sequence_integrity_nonsingletons,
                            "sequence_integrity_nonsingleton_scope": args.sequence_integrity_nonsingleton_scope,
                            "shadow_enabled": args.final_prune_opposite_strand_shadows,
                            "shadow_require_loser_single_exon": args.shadow_prune_require_loser_single_exon,
                            "singleton_thresholds": (
                                args.low_score_singleton_max_model_structural_support_score_scaled,
                                args.low_score_singleton_max_protein_evidence_score_scaled,
                                args.low_score_singleton_max_transcript_evidence_score_scaled,
                                args.low_score_singleton_max_model_selection_score_scaled,
                            ),
                            "nonsingleton_thresholds": (
                                args.low_score_nonsingleton_max_model_structural_support_score_scaled,
                                args.low_score_nonsingleton_max_protein_evidence_score_scaled,
                                args.low_score_nonsingleton_max_transcript_evidence_score_scaled,
                                args.low_score_nonsingleton_max_model_selection_score_scaled,
                            ),
                        }
                    preset = PRUNING_LEVEL_PRESETS[level]
                    scenario = PRUNING_SCENARIO_THRESHOLDS[level]
                    singleton_level = str(preset["singleton_level"])
                    if singleton_level == "intrinsic_or_sequence_integrity" and not genome_fasta_path:
                        singleton_level = "intrinsic"
                    external_supplied = bool(protein_evidence_paths or transcript_evidence_paths)
                    return {
                        "level": level,
                        "singleton_level": singleton_level,
                        "low_score_singletons": bool(
                            preset["low_score_singletons_when_protein_evidence"] and external_supplied
                        ),
                        "low_score_nonsingletons": bool(
                            preset["low_score_nonsingletons"]
                        ),
                        "low_score_nonsingleton_scope": scenario["low_score_nonsingleton_scope"],
                        "sequence_integrity_nonsingletons": bool(
                            preset["sequence_integrity_nonsingletons_when_genome"]
                            and genome_fasta_path
                            and not args.no_sequence_integrity_scoring
                        ),
                        "sequence_integrity_nonsingleton_scope": scenario[
                            "sequence_integrity_nonsingleton_scope"
                        ],
                        "shadow_enabled": bool(preset["shadow_enabled"]),
                        "shadow_require_loser_single_exon": bool(
                            preset["shadow_require_loser_single_exon"]
                        ),
                        "singleton_thresholds": scenario["singleton"],
                        "nonsingleton_thresholds": scenario["nonsingleton"],
                    }

                overlap_protein_adequacy_passed = bool(
                    protein_evidence_paths
                    and protein_overlap_model_fraction + 1e-15
                    >= args.low_score_singleton_min_protein_overlap_model_fraction
                    and protein_supported_model_fraction + 1e-15
                    >= args.low_score_singleton_min_protein_supported_model_fraction
                )
                overlap_transcript_adequacy_passed = bool(
                    transcript_evidence_paths
                    and transcript_scoring_eligible_model_fraction + 1e-15
                    >= args.low_score_singleton_min_transcript_scoring_eligible_model_fraction
                    and transcript_supported_model_fraction + 1e-15
                    >= args.low_score_singleton_min_transcript_supported_model_fraction
                )

                def run_pruning_configuration(
                    configuration: dict[str, object],
                    *,
                    scenario_progress: Optional[ProgressLogger],
                ) -> dict[str, object]:
                    level = str(configuration["level"])
                    singleton_thresholds = configuration["singleton_thresholds"]
                    nonsingleton_thresholds = configuration["nonsingleton_thresholds"]
                    result = make_initial_final_prune_result(
                        stage5_result,
                        subloci,
                        model_structural_rows_by_transcript,
                    )
                    if configuration["singleton_level"] != "none":
                        result = final_prune_singleton_selected_models(
                            result,
                            subloci,
                            model_structural_rows_by_transcript,
                            prune_level=str(configuration["singleton_level"]),
                            prune_sources=final_prune_singleton_subloci_sources,
                            pruning_level_name=level,
                            rescue_monosublocus_overlap=args.singleton_prune_rescue_monosublocus_overlap,
                            stage4_result=stage4_result,
                            rescue_min_holder_member_count=args.singleton_rescue_min_monosublocus_member_count,
                            rescue_min_overlap_bp=args.singleton_rescue_min_cds_overlap_bp,
                            rescue_min_self_overlap_fraction=args.singleton_rescue_min_self_overlap_fraction,
                            rescue_min_reciprocal_overlap_fraction=args.singleton_rescue_min_reciprocal_overlap_fraction,
                            monosublocus_overlap_rescue_evidence_policy=args.monosublocus_overlap_rescue_evidence_policy,
                            protein_adequacy_passed=overlap_protein_adequacy_passed,
                            transcript_adequacy_passed=overlap_transcript_adequacy_passed,
                            rescue_protein_supported_singletons=args.final_prune_rescue_protein_supported_singletons,
                            protein_evidence_rows_by_transcript=protein_evidence_rows_by_transcript,
                            rescue_min_protein_evidence_score_scaled=args.singleton_rescue_min_protein_evidence_score_scaled,
                            rescue_min_protein_exact_chain_source_count=args.singleton_rescue_min_protein_exact_chain_source_count,
                            rescue_transcript_supported_singletons=args.final_prune_rescue_transcript_supported_singletons,
                            transcript_evidence_rows_by_transcript=transcript_evidence_rows_by_transcript,
                            rescue_min_transcript_evidence_score_scaled=args.singleton_rescue_min_transcript_evidence_score_scaled,
                            rescue_min_transcript_supporting_source_count=args.singleton_rescue_min_transcript_supporting_source_count,
                            rescue_min_transcript_scoring_eligible_count=args.singleton_rescue_min_transcript_scoring_eligible_count,
                            rescue_min_transcript_cdna_context_fraction=args.singleton_rescue_min_transcript_cdna_context_fraction,
                            transcript_rescue_require_stranded_support=args.transcript_singleton_rescue_require_stranded_support,
                            transcript_rescue_block_on_cdna_context_concern=args.transcript_singleton_rescue_block_on_cdna_context_concern,
                            allow_internal_stop_evidence_rescue=args.allow_internal_stop_evidence_rescue,
                            progress=scenario_progress,
                        )
                    if configuration["low_score_singletons"]:
                        result = final_prune_low_score_singletons(
                            result,
                            subloci,
                            model_structural_rows_by_transcript,
                            enabled=True,
                            protein_evidence_supplied=bool(protein_evidence_paths),
                            transcript_evidence_supplied=bool(transcript_evidence_paths),
                            protein_adequacy_required=args.low_score_singleton_require_protein_evidence_adequacy,
                            transcript_adequacy_required=args.low_score_singleton_require_transcript_evidence_adequacy,
                            protein_overlap_model_fraction=protein_overlap_model_fraction,
                            protein_supported_model_fraction=protein_supported_model_fraction,
                            transcript_scoring_eligible_model_fraction=transcript_scoring_eligible_model_fraction,
                            transcript_supported_model_fraction=transcript_supported_model_fraction,
                            min_protein_overlap_model_fraction=args.low_score_singleton_min_protein_overlap_model_fraction,
                            min_protein_supported_model_fraction=args.low_score_singleton_min_protein_supported_model_fraction,
                            min_transcript_scoring_eligible_model_fraction=args.low_score_singleton_min_transcript_scoring_eligible_model_fraction,
                            min_transcript_supported_model_fraction=args.low_score_singleton_min_transcript_supported_model_fraction,
                            max_model_structural_support_score_scaled=singleton_thresholds[0],  # type: ignore[index]
                            max_protein_evidence_score_scaled=singleton_thresholds[1],  # type: ignore[index]
                            max_transcript_evidence_score_scaled=singleton_thresholds[2],  # type: ignore[index]
                            max_model_selection_score_scaled=singleton_thresholds[3],  # type: ignore[index]
                            pruning_level_name=level,
                            use_transcript_cdna_context=args.low_score_singleton_use_transcript_cdna_context,
                            progress=scenario_progress,
                        )
                    if configuration["low_score_nonsingletons"]:
                        result = final_prune_low_score_nonsingletons(
                            result,
                            subloci,
                            model_structural_rows_by_transcript,
                            enabled=True,
                            scope=str(configuration["low_score_nonsingleton_scope"]),
                            protein_evidence_supplied=bool(protein_evidence_paths),
                            transcript_evidence_supplied=bool(transcript_evidence_paths),
                            protein_adequacy_required=args.low_score_nonsingleton_require_protein_evidence_adequacy,
                            transcript_adequacy_required=args.low_score_nonsingleton_require_transcript_evidence_adequacy,
                            protein_overlap_model_fraction=protein_overlap_model_fraction,
                            protein_supported_model_fraction=protein_supported_model_fraction,
                            transcript_scoring_eligible_model_fraction=transcript_scoring_eligible_model_fraction,
                            transcript_supported_model_fraction=transcript_supported_model_fraction,
                            min_protein_overlap_model_fraction=args.low_score_nonsingleton_min_protein_overlap_model_fraction,
                            min_protein_supported_model_fraction=args.low_score_nonsingleton_min_protein_supported_model_fraction,
                            min_transcript_scoring_eligible_model_fraction=args.low_score_nonsingleton_min_transcript_scoring_eligible_model_fraction,
                            min_transcript_supported_model_fraction=args.low_score_nonsingleton_min_transcript_supported_model_fraction,
                            max_model_structural_support_score_scaled=nonsingleton_thresholds[0],  # type: ignore[index]
                            max_protein_evidence_score_scaled=nonsingleton_thresholds[1],  # type: ignore[index]
                            max_transcript_evidence_score_scaled=nonsingleton_thresholds[2],  # type: ignore[index]
                            max_model_selection_score_scaled=nonsingleton_thresholds[3],  # type: ignore[index]
                            pruning_level_name=level,
                            run_model_source_count=len(source_names),
                            use_transcript_cdna_context=args.low_score_nonsingleton_use_transcript_cdna_context,
                            progress=scenario_progress,
                        )
                    if configuration["sequence_integrity_nonsingletons"]:
                        result = final_prune_sequence_integrity_nonsingletons(
                            result,
                            subloci,
                            model_structural_rows_by_transcript,
                            enabled=True,
                            scope=str(configuration["sequence_integrity_nonsingleton_scope"]),
                            pruning_level_name=level,
                            protect_protein_supported=args.final_prune_rescue_protein_supported_singletons,
                            protein_evidence_rows_by_transcript=protein_evidence_rows_by_transcript,
                            min_protein_evidence_score_scaled=args.singleton_rescue_min_protein_evidence_score_scaled,
                            min_protein_exact_chain_source_count=args.singleton_rescue_min_protein_exact_chain_source_count,
                            protect_transcript_supported=args.final_prune_rescue_transcript_supported_singletons,
                            transcript_evidence_rows_by_transcript=transcript_evidence_rows_by_transcript,
                            min_transcript_evidence_score_scaled=args.singleton_rescue_min_transcript_evidence_score_scaled,
                            min_transcript_supporting_source_count=args.singleton_rescue_min_transcript_supporting_source_count,
                            min_transcript_scoring_eligible_count=args.singleton_rescue_min_transcript_scoring_eligible_count,
                            min_transcript_cdna_context_fraction=args.singleton_rescue_min_transcript_cdna_context_fraction,
                            transcript_require_stranded_support=args.transcript_singleton_rescue_require_stranded_support,
                            transcript_block_on_cdna_context_concern=args.transcript_singleton_rescue_block_on_cdna_context_concern,
                            allow_internal_stop_evidence_rescue=args.allow_internal_stop_evidence_rescue,
                            progress=scenario_progress,
                        )
                    if configuration["shadow_enabled"]:
                        result = final_prune_opposite_strand_shadows(
                            result,
                            model_structural_rows_by_transcript,
                            prune_sources=shadow_prune_sources,
                            require_loser_single_exon=bool(configuration["shadow_require_loser_single_exon"]),
                            min_overlap_bp=args.shadow_prune_min_overlap_bp,
                            min_loser_cds_overlap_fraction=args.shadow_prune_min_loser_cds_overlap_fraction,
                            min_score_margin_scaled=args.shadow_prune_min_score_margin_scaled,
                            pruning_level_name=level,
                            progress=scenario_progress,
                        )
                    return result

                scenarios: dict[str, ScenarioSummary] = {}
                scenario_validation: Optional[ScenarioValidation] = None
                standard_auto_decision: Optional[AutoDecision] = None
                conservative_auto_decision: Optional[AutoDecision] = None
                flood_consensus_assessment = None
                if automatic_assessment_requested:
                    protein_adequacy = classify_evidence_adequacy(
                        "protein",
                        supplied=bool(protein_evidence_paths),
                        first_value=protein_overlap_model_fraction,
                        second_value=protein_supported_model_fraction,
                        first_threshold=args.low_score_singleton_min_protein_overlap_model_fraction,
                        second_threshold=args.low_score_singleton_min_protein_supported_model_fraction,
                    )
                    transcript_adequacy = classify_evidence_adequacy(
                        "transcript",
                        supplied=bool(transcript_evidence_paths),
                        first_value=transcript_scoring_eligible_model_fraction,
                        second_value=transcript_supported_model_fraction,
                        first_threshold=args.low_score_singleton_min_transcript_scoring_eligible_model_fraction,
                        second_threshold=args.low_score_singleton_min_transcript_supported_model_fraction,
                    )
                    initial_models = list(stage5_result["final_models"])  # type: ignore[arg-type,index]
                    # Consensus is measured from the complete model-input
                    # sublocus distribution, matching the existing Sublocus
                    # summary and the calibration fixtures.
                    consensus_metrics = calculate_source_consensus_metrics(subloci)
                    flood_consensus_assessment = assess_flood_adjusted_consensus(
                        consensus_metrics,
                        enabled=args.auto_pruning_use_source_excess,
                    )
                    sublocus_model_count_by_id = {
                        member.output_transcript_id: len(sub.models)
                        for sub in subloci
                        for member in sub.models
                    }
                    singleton_selected = sum(
                        sublocus_model_count_by_id.get(model.output_transcript_id, 0) == 1
                        for model in initial_models
                    )
                    selected_singleton_model_fraction = (
                        singleton_selected / len(initial_models) if initial_models else 0.0
                    )
                    previous_scenario: Optional[ScenarioSummary] = None
                    for scenario_level in AUTO_SCENARIO_LEVELS:
                        configuration = pruning_configuration(scenario_level)
                        scenario_result = run_pruning_configuration(
                            configuration,
                            scenario_progress=None,
                        )
                        if previous_scenario is not None:
                            scenario_result = apply_monotonic_pruning_closure(
                                scenario_result,
                                required_removed_ids=previous_scenario.removed_ids,
                                lower_level=previous_scenario.level,
                                current_level=scenario_level,
                                lower_provenance_by_id=previous_scenario.removal_provenance_by_id,
                            )
                        singleton_thresholds = configuration["singleton_thresholds"]
                        nonsingleton_thresholds = configuration["nonsingleton_thresholds"]
                        thresholds_text = (
                            "singleton="
                            + "/".join(f"{float(value):.2f}" for value in singleton_thresholds)  # type: ignore[arg-type]
                            + ";nonsingleton="
                            + "/".join(f"{float(value):.2f}" for value in nonsingleton_thresholds)  # type: ignore[arg-type]
                            + ";nonsingleton_policy="
                            + str(scenario_result.get("low_score_nonsingleton_policy", "."))
                            + f";run_model_sources={len(source_names)}"
                            + (
                                ";fallback_model=medium:0.20,high:0.25,very-high:0.25;"
                                "fallback_depth_tiers=0.35/0.45"
                                if scenario_result.get("low_score_nonsingleton_policy")
                                == "model_confidence_fallback"
                                else ""
                            )
                        )
                        scenarios[scenario_level] = summarize_pruning_scenario(
                            scenario_level,
                            scenario_result,
                            initial_models=initial_models,
                            subloci=subloci,
                            structural_rows_by_id=model_structural_rows_by_transcript,
                            protein_adequacy_state=protein_adequacy.state,
                            transcript_adequacy_state=transcript_adequacy.state,
                            protein_evidence_supplied=bool(protein_evidence_paths),
                            transcript_evidence_supplied=bool(transcript_evidence_paths),
                            scenario_configuration=configuration,
                            thresholds=thresholds_text,
                            low_score_nonsingleton_scope=str(
                                configuration["low_score_nonsingleton_scope"]
                            ),
                            sequence_integrity_nonsingleton_scope=str(
                                configuration["sequence_integrity_nonsingleton_scope"]
                            ),
                        )
                        previous_scenario = scenarios[scenario_level]
                        del scenario_result
                    add_incremental_metrics(scenarios)
                    scenario_validation = validate_scenarios(scenarios)
                    adequacy_pass_count = int(protein_adequacy.passed) + int(transcript_adequacy.passed)
                    standard_max = (
                        args.auto_pruning_max_level
                        if args.auto_pruning_max_level_explicit
                        else "very-high"
                    )
                    conservative_max = (
                        args.auto_pruning_max_level
                        if args.auto_pruning_max_level_explicit
                        else "high"
                    )
                    standard_auto_decision = choose_auto_pruning_level(
                        "auto",
                        scenarios=scenarios,
                        validation=scenario_validation,
                        consensus=consensus_metrics,
                        adjusted_consensus=flood_consensus_assessment.adjusted,
                        adequacy_pass_count=adequacy_pass_count,
                        selected_singleton_model_fraction=selected_singleton_model_fraction,
                        min_level=args.auto_pruning_min_level,
                        max_level=standard_max,
                        max_level_explicit=args.auto_pruning_max_level_explicit,
                        max_extra_removal_fraction=args.auto_pruning_max_extra_removal_fraction,
                        high_min_singleton_fraction=args.auto_pruning_high_min_singleton_fraction,
                        high_min_prospective_removal_fraction=args.auto_pruning_high_min_prospective_removal_fraction,
                        use_source_excess=args.auto_pruning_use_source_excess,
                    )
                    conservative_auto_decision = choose_auto_pruning_level(
                        "auto-conservative",
                        scenarios=scenarios,
                        validation=scenario_validation,
                        consensus=consensus_metrics,
                        adjusted_consensus=flood_consensus_assessment.adjusted,
                        adequacy_pass_count=adequacy_pass_count,
                        selected_singleton_model_fraction=selected_singleton_model_fraction,
                        min_level=args.auto_pruning_min_level,
                        max_level=conservative_max,
                        max_level_explicit=args.auto_pruning_max_level_explicit,
                        max_extra_removal_fraction=args.auto_pruning_max_extra_removal_fraction,
                        high_min_singleton_fraction=args.auto_pruning_high_min_singleton_fraction,
                        high_min_prospective_removal_fraction=args.auto_pruning_high_min_prospective_removal_fraction,
                        use_source_excess=args.auto_pruning_use_source_excess,
                    )
                    if output_enabled("pruning_scenario_summary_tsv"):
                        write_scenario_summary(pruning_scenario_summary_tsv, scenarios)
                        note_output(
                            "pruning_scenario_summary_tsv",
                            "pruning scenario summary TSV",
                            pruning_scenario_summary_tsv,
                        )
                    if output_enabled("auto_pruning_diagnostics_tsv"):
                        diagnostic_rows = automatic_diagnostics_rows(
                            protein=protein_adequacy,
                            transcript=transcript_adequacy,
                            consensus=consensus_metrics,
                            selected_singleton_model_fraction=selected_singleton_model_fraction,
                            validation=scenario_validation,
                            standard_decision=standard_auto_decision,
                            conservative_decision=conservative_auto_decision,
                            flood_assessment=flood_consensus_assessment,
                            input_count_metrics=input_count_metrics,
                            scenarios=scenarios,
                        )
                        write_auto_diagnostics(auto_pruning_diagnostics_tsv, diagnostic_rows)
                        note_output(
                            "auto_pruning_diagnostics_tsv",
                            "automatic pruning diagnostics TSV",
                            auto_pruning_diagnostics_tsv,
                        )

                if args.pruning_level == "auto":
                    assert standard_auto_decision is not None
                    effective_pruning_level = standard_auto_decision.effective_level
                    effective_configuration = pruning_configuration(effective_pruning_level)
                elif args.pruning_level == "auto-conservative":
                    assert conservative_auto_decision is not None
                    effective_pruning_level = conservative_auto_decision.effective_level
                    effective_configuration = pruning_configuration(effective_pruning_level)
                else:
                    effective_pruning_level = args.pruning_level
                    effective_configuration = pruning_configuration(
                        effective_pruning_level,
                        effective_args=True,
                    )
                final_prune_result = run_pruning_configuration(
                    effective_configuration,
                    scenario_progress=progress,
                )
                if effective_pruning_level in {"high", "very-high"}:
                    lower_level = "medium" if effective_pruning_level == "high" else "high"
                    if scenarios:
                        lower_scenario = scenarios[lower_level]
                        final_prune_result = apply_monotonic_pruning_closure(
                            final_prune_result,
                            required_removed_ids=lower_scenario.removed_ids,
                            lower_level=lower_level,
                            current_level=effective_pruning_level,
                            lower_provenance_by_id=lower_scenario.removal_provenance_by_id,
                        )
                    elif effective_configuration == pruning_configuration(effective_pruning_level):
                        medium_result = run_pruning_configuration(
                            pruning_configuration("medium"),
                            scenario_progress=None,
                        )
                        if effective_pruning_level == "high":
                            lower_result = medium_result
                            lower_level = "medium"
                        else:
                            high_result = run_pruning_configuration(
                                pruning_configuration("high"),
                                scenario_progress=None,
                            )
                            high_result = apply_monotonic_pruning_closure(
                                high_result,
                                required_removed_ids={
                                    str(row.get("output_transcript_id", "."))
                                    for row in medium_result.get("removed_rows", [])
                                },
                                lower_level="medium",
                                current_level="high",
                                lower_provenance_by_id=removal_provenance_from_result(
                                    medium_result, level="medium"
                                ),
                            )
                            lower_result = high_result
                            lower_level = "high"
                        final_prune_result = apply_monotonic_pruning_closure(
                            final_prune_result,
                            required_removed_ids={
                                str(row.get("output_transcript_id", "."))
                                for row in lower_result.get("removed_rows", [])
                            },
                            lower_level=lower_level,
                            current_level=effective_pruning_level,
                            lower_provenance_by_id=removal_provenance_from_result(
                                lower_result, level=lower_level
                            ),
                        )
                final_prune_result["requested_pruning_level"] = args.pruning_level
                final_prune_result["effective_pruning_level"] = effective_pruning_level
                final_prune_result["effective_pruning_configuration"] = dict(effective_configuration)
                final_prune_result["standard_auto_decision"] = standard_auto_decision
                final_prune_result["conservative_auto_decision"] = conservative_auto_decision
                final_prune_result["auto_scenarios"] = scenarios
                final_prune_result["auto_scenario_validation"] = scenario_validation
                final_prune_result["auto_consensus_metrics"] = (
                    consensus_metrics if automatic_assessment_requested else None
                )
                final_prune_result["auto_flood_consensus_assessment"] = (
                    flood_consensus_assessment
                )
                final_prune_result["auto_input_count_metrics"] = input_count_metrics
                if scenarios:
                    annotate_result_with_scenario_audit(final_prune_result, scenarios)
                write_final_pruned_outputs(
                    final_prune_result,
                    final_pruned_selected_models_tsv=final_pruned_selected_models_tsv if output_enabled("final_pruned_selected_models_tsv") else None,
                    final_pruned_selected_removed_tsv=final_pruned_selected_removed_tsv if output_enabled("final_pruned_selected_removed_tsv") else None,
                    final_shadow_gene_pairs_tsv=final_shadow_gene_pairs_tsv if output_enabled("final_shadow_gene_pairs_tsv") else None,
                )
                note_output("final_pruned_selected_models_tsv", "final selected pruned TSV", final_pruned_selected_models_tsv)
                note_output("final_pruned_selected_removed_tsv", "final selected pruned removed TSV", final_pruned_selected_removed_tsv)
                note_output("final_shadow_gene_pairs_tsv", "final shadow gene pairs TSV", final_shadow_gene_pairs_tsv)
                if output_enabled("final_pruned_selected_models_gff3"):
                    write_model_subset_gff3(
                        final_prune_result["retained_models"],  # type: ignore[arg-type]
                        final_pruned_selected_models_gff3,
                        model_structural_rows_by_transcript=model_structural_rows_by_transcript,
                        subloci=subloci,
                        selection_info_by_transcript=final_prune_result["retained_selection_info_by_transcript"],  # type: ignore[arg-type]
                        monosublocus_rows=monosublocus_rows,  # type: ignore[arg-type]
                        include_model_selection_components=args.gff3_include_model_selection_score_components,
                        include_monosublocus_features=True,
                        selection_label="final_selected_model",
                        cds_only_model_structure=args.final_gff3_cds_only_structure,
                    )
                    note_output("final_pruned_selected_models_gff3", "final selected pruned GFF3", final_pruned_selected_models_gff3)
                progress.log("Finished final selected-model pruning outputs")
            progress.log("Finished writing Backfill rescue backfill/final selected outputs")
            timer.mark("Writing backfill/final outputs")

    frozen_primary_annotation = None
    final_annotation = None
    if args.stop_after == "final" and (
        output_enabled("standardised_primary_models_gff3")
        or output_enabled("final_annotation_gff3")
        or alternative_splice_analysis_requested
        or args.add_utrs
        or (args.refine_cds_termini and genome_fasta_path is not None)
    ):
        assert stage5_result is not None
        frozen_primary_annotation = build_standardised_primary_annotation(
            final_selected_models=stage5_result["final_models"],  # type: ignore[arg-type]
            final_pruned_models=(
                final_prune_result["retained_models"]  # type: ignore[arg-type]
                if final_prune_result is not None
                else None
            ),
            model_score_rows_by_transcript=model_structural_rows_by_transcript,
        )
        final_annotation = frozen_primary_annotation.annotation
        if output_enabled("standardised_primary_models_gff3"):
            progress.log("Writing standardised primary annotation")
            write_standardised_primary_gff3(
                final_annotation,
                standardised_primary_models_gff3,
                cds_only_structure=args.final_gff3_cds_only_structure,
            )
            note_output(
                "standardised_primary_models_gff3",
                "standardised primary annotation GFF3",
                standardised_primary_models_gff3,
            )
            progress.log("Finished standardised primary annotation")
            timer.mark("Standardised primary annotation")

    alternative_splice_report_result: Optional[dict[str, list[dict[str, object]]]] = None
    alternative_splice_selection_config = None
    alternative_splice_selection_result = None
    public_final_annotation = final_annotation
    if alternative_splice_analysis_requested:
        assert stage5_result is not None
        assert final_annotation is not None
        assert frozen_primary_annotation is not None
        progress.log("Starting alternative-splice selection")
        alternative_splice_selection_settings = AlternativeSpliceSelectionSettings(
            preset=args.alternative_splice_preset,
            support_policy=args.alternative_splice_support_policy,
            evidence_requirement=args.alternative_splice_evidence_requirement,
            max_alternatives_per_gene=args.alternative_splice_max_alternatives_per_gene,
            min_model_support_sources=args.alternative_splice_min_model_support_sources,
            min_model_support_source_fraction=(
                args.alternative_splice_min_model_support_source_fraction
            ),
            min_reciprocal_cds_overlap=(
                args.alternative_splice_min_reciprocal_cds_overlap
            ),
            min_frame_concordance=args.alternative_splice_min_frame_concordance,
            min_model_selection_score_scaled=(
                args.alternative_splice_min_model_selection_score_scaled
            ),
            min_primary_score_ratio=args.alternative_splice_min_primary_score_ratio,
        )
        try:
            alternative_splice_stage_result = evaluate_alternative_splice_stage(
                frozen_primary=frozen_primary_annotation,
                all_scored_models=coding_models,
                model_score_rows_by_transcript=model_structural_rows_by_transcript,
                protein_rows_by_transcript=protein_evidence_rows_by_transcript,
                protein_collapsed_alignments=protein_evidence_collapsed_alignments,
                transcript_rows_by_transcript=transcript_evidence_rows_by_transcript,
                transcript_superloci=transcript_evidence_superloci,
                report_policy=AlternativeSpliceReportPolicy(
                    terminal_cds_3prime_wobble_bp=args.terminal_cds_3prime_wobble_bp,
                    protein_terminal_cds_3prime_wobble_bp=protein_evidence_wobble_bp,
                    transcript_min_cdna_context_fraction=args.singleton_rescue_min_transcript_cdna_context_fraction,
                    transcript_require_stranded_support=args.transcript_singleton_rescue_require_stranded_support,
                    transcript_block_on_cdna_context_concern=args.transcript_singleton_rescue_block_on_cdna_context_concern,
                    transcript_cdna_context_neutral_fraction=args.transcript_cdna_context_neutral_fraction,
                    transcript_cdna_context_severe_fraction=args.transcript_cdna_context_severe_fraction,
                    proposed_max_alternatives_per_gene=args.alternative_splice_max_alternatives_per_gene,
                ),
                run_metadata={
                    "requested_pruning_level": args.pruning_level,
                    "effective_pruning_level": (
                        final_prune_result.get("effective_pruning_level", args.pruning_level)
                        if final_prune_result is not None
                        else args.pruning_level
                    ),
                    "protein_evidence_source_count": len(protein_evidence_source_labels),
                    "transcript_evidence_source_count": len(transcript_evidence_source_labels),
                    "genome_sequence_configured": genome_fasta_path is not None,
                },
                selection_settings=alternative_splice_selection_settings,
                terminal_cds_3prime_wobble_bp=args.terminal_cds_3prime_wobble_bp,
                alternatives_emitted=alternative_splice_enabled,
            )
        except AlternativeSpliceConfigurationError as exc:
            raise SystemExit(f"Invalid alternative-splice configuration: {exc}") from exc
        alternative_splice_fact_result = alternative_splice_stage_result.fact_report
        alternative_splice_selection_config = alternative_splice_stage_result.config
        alternative_splice_selection_result = alternative_splice_stage_result.selection
        alternative_splice_report_result = alternative_splice_stage_result.report
        public_final_annotation = alternative_splice_stage_result.final_annotation
        write_alternative_splice_reports(
            alternative_splice_report_result,
            candidate_tsv=(
                alternative_splice_candidates_tsv
                if output_enabled("alternative_splice_candidates_tsv") else None
            ),
            gene_tsv=(
                alternative_splice_genes_tsv
                if output_enabled("alternative_splice_genes_tsv") else None
            ),
            run_tsv=(
                alternative_splice_run_summary_tsv
                if output_enabled("alternative_splice_run_summary_tsv") else None
            ),
        )
        note_output(
            "alternative_splice_candidates_tsv",
            "alternative-splice candidate diagnostics TSV",
            alternative_splice_candidates_tsv,
        )
        note_output(
            "alternative_splice_genes_tsv",
            "alternative-splice retained-gene summary TSV",
            alternative_splice_genes_tsv,
        )
        note_output(
            "alternative_splice_run_summary_tsv",
            "alternative-splice run summary TSV",
            alternative_splice_run_summary_tsv,
        )
        if (
            alternative_splice_audit_requested
            and output_enabled("alternative_splice_candidate_audit_gff3")
        ):
            write_alternative_splice_candidate_audit_gff3(
                candidate_rows=alternative_splice_fact_result["candidate_rows"],
                selection_result=alternative_splice_selection_result,
                all_scored_models=coding_models,
                primary_annotation=final_annotation,
                path=alternative_splice_candidate_audit_gff3,
            )
            note_output(
                "alternative_splice_candidate_audit_gff3",
                "alternative-splice candidate audit GFF3",
                alternative_splice_candidate_audit_gff3,
            )
        genes_with_alternatives = sum(
            bool(gene.accepted_alternatives)
            for gene in alternative_splice_selection_result.gene_decisions
        )
        progress.log(
            "Finished alternative-splice selection: "
            f"accepted={len(alternative_splice_selection_result.accepted_alternatives):,}; "
            f"genes_with_alternatives={genes_with_alternatives:,}"
        )
        timer.mark("Alternative-splice selection")

    prediction_utr_evidence_catalogue = []
    utr_settings = None
    learned_utr_limits = None
    cds_refinement_result = None
    has_incomplete_final_cds = bool(
        public_final_annotation is not None
        and any(
            transcript.source_model.sequence_integrity is not None
            and (
                transcript.source_model.sequence_integrity.start_codon_status == "missing"
                or transcript.source_model.sequence_integrity.stop_codon_status == "missing"
            )
            for gene in public_final_annotation.genes
            for transcript in gene.transcripts
        )
    )
    refinement_execution_requested = bool(
        args.refine_cds_termini
        and genome_fasta_path is not None
        and has_incomplete_final_cds
    )
    if (
        refinement_execution_requested
        and not args.add_utrs
        and transcript_evidence_paths
    ):
        progress.log(
            "Re-reading transcript evidence for deferred CDS-completion catalogue"
        )
        completion_raw_transcripts = parse_transcript_evidence(
            transcript_evidence_paths,
            transcript_evidence_source_labels,
            strict=args.strict,
        )
        transcript_utr_evidence_catalogue = build_transcript_utr_catalogue(
            completion_raw_transcripts
        )
        released_completion_raw = len(completion_raw_transcripts)
        completion_raw_transcripts.clear()
        gc.collect()
        progress.log(
            "Built deferred CDS-completion transcript catalogue: "
            f"records={len(transcript_utr_evidence_catalogue):,}; "
            f"released_reparsed_raw={released_completion_raw:,}",
            level="detailed",
        )
        timer.mark("Deferred CDS-completion transcript catalogue")
    postselection_structure_requested = bool(
        public_final_annotation is not None
        and (
            args.add_utrs
            or refinement_execution_requested
        )
    )
    if postselection_structure_requested:
        assert public_final_annotation is not None
        selected_internal_ids = {
            transcript.source_model.output_transcript_id
            for gene in public_final_annotation.genes
            for transcript in gene.transcripts
        }
        prediction_utr_evidence_catalogue = build_prediction_utr_catalogue(
            coding_models,
            selected_internal_ids=selected_internal_ids,
        )
        utr_settings = UTRSettings(
            min_complete_exon_bp=args.utr_min_complete_exon_bp,
            min_intron_bp=(
                args.min_cds_intron_bp
                if args.utr_min_intron_bp is None
                else args.utr_min_intron_bp or None
            ),
            max_intron_bp=(
                long_cds_intron_cutoff
                if args.utr_max_intron_bp is None
                else args.utr_max_intron_bp or None
            ),
            max_five_prime_exons=args.utr_max_five_prime_exons or None,
            max_three_prime_exons=args.utr_max_three_prime_exons or None,
            exonic_failsafe_bp=args.utr_exonic_failsafe_bp or None,
            five_prime_exonic_override_bp=args.utr_five_prime_exonic_override_bp,
            three_prime_exonic_override_bp=args.utr_three_prime_exonic_override_bp,
            five_prime_genomic_override_bp=args.utr_five_prime_genomic_override_bp,
            three_prime_genomic_override_bp=args.utr_three_prime_genomic_override_bp,
            member_depth_cap=args.utr_member_depth_cap,
            endpoint_window_min_bp=args.utr_endpoint_window_min_bp,
            endpoint_window_max_bp=args.utr_endpoint_window_max_bp,
            local_consensus_min_fraction=args.utr_local_consensus_min_fraction,
            limit_min_observations=args.utr_limit_min_observations,
        )

    if args.add_utrs:
        assert public_final_annotation is not None and utr_settings is not None
        if indexed_fasta is None:
            learned_utr_limits = learn_utr_limits_for_annotation(
                public_final_annotation,
                transcript_utr_evidence_catalogue,
                prediction_utr_evidence_catalogue,
                utr_settings,
            )
        else:
            with indexed_fasta as fasta_reader:
                learned_utr_limits = learn_utr_limits_for_annotation(
                    public_final_annotation,
                    transcript_utr_evidence_catalogue,
                    prediction_utr_evidence_catalogue,
                    utr_settings,
                    fasta=fasta_reader,
                )

    if (
        refinement_execution_requested
        and public_final_annotation is not None
    ):
        assert utr_settings is not None
        progress.log("Starting post-selection CDS terminus refinement")
        refinement_settings = CDSTerminusRefinementSettings(
            downstream_start_search_bp=args.cds_refinement_downstream_start_search_bp,
            downstream_stop_search_bp=args.cds_refinement_downstream_stop_search_bp,
            max_relative_extension=args.cds_refinement_max_relative_extension,
            max_relative_truncation=args.cds_refinement_max_relative_truncation,
            min_cds_length_bp=args.min_cds_length_bp,
            single_exon_min_reciprocal_overlap=args.single_exon_min_reciprocal_overlap,
            genetic_code=args.genetic_code,
            internal_stop_penalty=args.sequence_penalty_internal_stop,
            internal_stop_count_cap=args.sequence_internal_stop_count_cap,
            phase_inconsistent_penalty=args.sequence_penalty_phase_inconsistent,
            missing_stop_penalty=args.sequence_penalty_missing_stop,
            missing_start_penalty=args.sequence_penalty_missing_start,
            noncanonical_splice_penalty=args.sequence_penalty_noncanonical_splice,
        )
        try:
            with indexed_fasta as fasta_reader:
                cds_refinement_result = refine_cds_termini(
                    public_final_annotation,
                    transcript_utr_evidence_catalogue,
                    prediction_utr_evidence_catalogue,
                    fasta_reader,
                    refinement_settings,
                    utr_settings=utr_settings,
                    learned_utr_limits=learned_utr_limits,
                )
        except ValueError as exc:
            raise SystemExit(f"Invalid CDS terminus refinement: {exc}") from exc
        public_final_annotation = cds_refinement_result.annotation
        refinement_summary = cds_refinement_result.summary
        progress.log(
            "Finished CDS terminus refinement: "
            f"audit_sides={refinement_summary.audit_sides:,}; "
            f"missing_start_sides={refinement_summary.incomplete_five_prime_sides:,}; "
            f"upstream_refined_sides={refinement_summary.refined_sides_upstream_start:,}; "
            f"downstream_fallback_refined_sides="
            f"{refinement_summary.refined_sides_downstream_start_fallback:,}; "
            f"missing_stop_sides={refinement_summary.incomplete_three_prime_sides:,}; "
            f"stop_refined_sides={refinement_summary.refined_sides_downstream_stop:,}; "
            f"final_refined_sides={refinement_summary.final_refined_sides:,}; "
            f"final_failed_sides={refinement_summary.final_failed_sides:,}; "
            f"refined_transcripts={refinement_summary.refined_transcripts:,}; "
            f"unchanged_transcripts={refinement_summary.unchanged_transcripts:,}; "
            f"encountered_rejection_events="
            f"{sum(dict(refinement_summary.encountered_rejection_event_counts).values()):,}"
        )
        timer.mark("CDS terminus refinement")

    if cds_refinement_result is None:
        summary_annotation = public_final_annotation or FinalAnnotation(())
        if not args.refine_cds_termini:
            refinement_status = "disabled"
        elif genome_fasta_path is None:
            refinement_status = "prerequisite_unavailable:no_genome_fasta"
        elif public_final_annotation is None:
            refinement_status = "not_applicable:no_final_annotation"
        elif not has_incomplete_final_cds:
            refinement_status = "zero_incomplete_models"
        else:
            refinement_status = "not_run"
        cds_refinement_result = empty_cds_refinement_result(
            summary_annotation,
            execution_status=refinement_status,
        )
    if output_enabled("cds_terminus_refinement_audit_tsv"):
        write_cds_refinement_audit_tsv(
            cds_refinement_result, cds_terminus_refinement_audit_tsv
        )
        note_output(
            "cds_terminus_refinement_audit_tsv",
            "CDS terminus refinement audit TSV",
            cds_terminus_refinement_audit_tsv,
        )
    if output_enabled("cds_terminus_refinement_summary_tsv"):
        write_cds_refinement_summary_tsv(
            cds_refinement_result, cds_terminus_refinement_summary_tsv
        )
        note_output(
            "cds_terminus_refinement_summary_tsv",
            "CDS terminus refinement summary TSV",
            cds_terminus_refinement_summary_tsv,
        )

    inherited_utr_validation_results = {}
    utr_validation_recorder = UTRValidationRecorder(
        utr_junction_validation_audit_tsv
        if output_enabled("utr_junction_validation_audit_tsv")
        else None
    )
    if public_final_annotation is not None:
        effective_utr_min_intron_bp = (
            args.min_cds_intron_bp
            if args.utr_min_intron_bp is None
            else args.utr_min_intron_bp or 0
        )
        if indexed_fasta is None:
            inherited_utr_validation_results = inherited_utr_validations(
                public_final_annotation,
                minimum_intron_bp=effective_utr_min_intron_bp,
            )
        else:
            with indexed_fasta as fasta_reader:
                inherited_utr_validation_results = inherited_utr_validations(
                    public_final_annotation,
                    minimum_intron_bp=effective_utr_min_intron_bp,
                    fasta=fasta_reader,
                )
        for by_side in inherited_utr_validation_results.values():
            for validation in by_side.values():
                utr_validation_recorder.record(
                    validation, rejected_action="suppress_source_side"
                )

    if (
        args.stop_after == "final"
        and output_enabled("final_annotation_gff3")
        and public_final_annotation is not None
    ):
        progress.log("Writing final annotation")
        write_final_annotation_gff3(
            public_final_annotation,
            final_annotation_gff3,
            FinalAnnotationOutputMetadata(
                alternative_splice_enabled=alternative_splice_enabled,
                selector_evaluated=alternative_splice_selection_config is not None,
                alternatives_emitted=(
                    alternative_splice_enabled
                    and alternative_splice_selection_config is not None
                ),
                selection_config=alternative_splice_selection_config,
            ),
            cds_only_structure=args.final_gff3_cds_only_structure,
            utr_validations=inherited_utr_validation_results,
        )
        note_output(
            "final_annotation_gff3",
            "final annotation GFF3",
            final_annotation_gff3,
        )
        progress.log(
            "Finished final annotation: "
            f"genes={public_final_annotation.gene_count:,}; "
            f"transcripts={public_final_annotation.transcript_count:,}"
        )
        timer.mark("Final annotation")

    utr_augmentation_result = None
    if args.add_utrs:
        assert public_final_annotation is not None and utr_settings is not None
        progress.log("Starting post-selection UTR augmentation")
        try:
            augmentation_kwargs = {
                "learned_limits": learned_utr_limits,
                "excluded_training_sides": (
                    set(cds_refinement_result.originally_incomplete_sides)
                    if cds_refinement_result is not None
                    else None
                ),
                "validation_recorder": lambda validation: (
                    utr_validation_recorder.record(
                        validation, rejected_action="exclude_candidate_side"
                    )
                ),
            }
            if indexed_fasta is None:
                utr_augmentation_result = augment_final_annotation(
                    public_final_annotation,
                    transcript_utr_evidence_catalogue,
                    prediction_utr_evidence_catalogue,
                    utr_settings,
                    **augmentation_kwargs,
                )
            else:
                with indexed_fasta as fasta_reader:
                    utr_augmentation_result = augment_final_annotation(
                        public_final_annotation,
                        transcript_utr_evidence_catalogue,
                        prediction_utr_evidence_catalogue,
                        utr_settings,
                        fasta=fasta_reader,
                        **augmentation_kwargs,
                    )
        except ValueError as exc:
            raise SystemExit(f"Invalid UTR augmentation configuration: {exc}") from exc
        if output_enabled("utr_augmented_final_annotation_gff3"):
            write_utr_augmented_gff3(
                utr_augmentation_result,
                utr_augmented_final_annotation_gff3,
            )
            note_output(
                "utr_augmented_final_annotation_gff3",
                "UTR-augmented final annotation GFF3",
                utr_augmented_final_annotation_gff3,
            )
        if output_enabled("utr_augmentation_audit_tsv"):
            write_utr_audit_tsv(utr_augmentation_result, utr_augmentation_audit_tsv)
            note_output(
                "utr_augmentation_audit_tsv",
                "UTR augmentation audit TSV",
                utr_augmentation_audit_tsv,
            )
        if output_enabled("utr_augmentation_summary_tsv"):
            write_utr_summary_tsv(
                utr_augmentation_result, utr_augmentation_summary_tsv
            )
            note_output(
                "utr_augmentation_summary_tsv",
                "UTR augmentation summary TSV",
                utr_augmentation_summary_tsv,
            )
        for warning in utr_augmentation_result.limits.warnings:
            progress.log(f"UTR augmentation warning: {warning}")
        utr_summary = utr_augmentation_result.summary
        progress.log(
            "Finished UTR augmentation: "
            f"transcripts={utr_summary.transcript_count:,}; "
            f"both={utr_summary.both_assigned:,}; "
            f"five_prime_only={utr_summary.five_prime_only:,}; "
            f"three_prime_only={utr_summary.three_prime_only:,}; "
            f"neither={utr_summary.neither_assigned:,}; "
            f"runtime={utr_summary.runtime_seconds:.1f}s"
        )
        timer.mark("UTR augmentation")

    utr_validation_recorder.close()
    if output_enabled("utr_junction_validation_audit_tsv"):
        note_output(
            "utr_junction_validation_audit_tsv",
            "UTR junction validation audit TSV",
            utr_junction_validation_audit_tsv,
        )
    if output_enabled("utr_junction_validation_summary_tsv"):
        write_utr_validation_summary(
            utr_validation_recorder,
            utr_junction_validation_summary_tsv,
        )
        note_output(
            "utr_junction_validation_summary_tsv",
            "UTR junction validation summary TSV",
            utr_junction_validation_summary_tsv,
        )

    if final_prune_result is not None:
        resolved_pruning_level = str(
            final_prune_result.get("effective_pruning_level", args.pruning_level)
        )
    elif args.pruning_level in {"auto", "auto-conservative"}:
        resolved_pruning_level = None
    else:
        resolved_pruning_level = str(args.pruning_level)
    if output_enabled("run_configuration_json"):
        write_run_configuration_json(
            build_run_configuration(
                args,
                raw_arguments=raw_args,
                enabled_output_keys=context.enabled_outputs,
                effective_pruning_level=resolved_pruning_level,
            ),
            run_configuration_json,
        )
        note_output(
            "run_configuration_json",
            "resolved run configuration JSON",
            run_configuration_json,
        )

    all_sources = sorted(source_names)

    log_section("Input summary")
    print(f"Input files: {len(args.inputs):,}")
    print("source\tfeatures\tcoding_models\texcluded_models\tunused_features\tfile")
    for parsed in parsed_inputs:
        print(
            f"{parsed.source_name}\t{len(parsed.records):,}\t"
            f"{len(parsed.coding_models):,}\t{len(parsed.excluded_models):,}\t"
            f"{len(parsed.records) - len(parsed.used_record_indices):,}\t{parsed.input_file}"
        )
    print(f"Valid coding models after model-level validation: {len(coding_models):,}")
    print(f"Excluded CDS models: {len(excluded_models):,}")
    if excluded_models:
        reason_counts = Counter(item.reason for item in excluded_models)
        print(f"Excluded CDS model reasons: {format_counter(reason_counts)}")

    log_section("Validation settings")
    print(
        "Model validation: "
        f"action={args.model_validation_action}; "
        f"max_model_cds_span_bp={args.max_model_cds_span_bp:,}; "
        f"max_model_cds_intron_bp={args.max_model_cds_intron_bp:,}"
    )
    if model_failures:
        model_failure_counts = Counter(f.failure_type for f in model_failures)
        print(f"Model validation failure counts: {format_counter(model_failure_counts)}")
    else:
        print("Model validation failure counts: none")
    print(
        "Superlocus validation: "
        f"action={args.superlocus_validation_action}; "
        f"max_superlocus_span_bp={args.max_superlocus_span_bp:,}; "
        f"max_superlocus_gene_count={args.max_superlocus_gene_count:,}; "
        f"max_genes_per_source_in_superlocus={args.max_genes_per_source_in_superlocus:,}"
    )
    if superlocus_failures:
        superlocus_failure_counts = Counter(f.failure_type for f in superlocus_failures)
        print(f"Superlocus validation failure counts: {format_counter(superlocus_failure_counts)}")
    else:
        print("Superlocus validation failure counts: none")

    log_section("Superlocus summary")
    print(f"Superloci: {len(superloci):,}")
    model_counts = Counter(len(sl.models) for sl in superloci)
    source_counts = Counter(len(superlocus_sources(sl)) for sl in superloci)
    gene_counts = Counter(superlocus_gene_count(sl) for sl in superloci)
    log_distribution("Model-count distribution", model_counts, "models", "superloci")
    log_distribution("Source-count distribution", source_counts, "sources", "superloci")
    log_distribution("Gene-count distribution", gene_counts, "genes", "superloci")
    log_subsection("Source-count breakdown by source presence")
    log_source_presence_breakdown(
        "Superloci containing each source within each source-count class",
        "superloci",
        superloci,
        superlocus_sources,
        all_sources,
    )
    log_subsection("Source-count breakdown by source combination")
    log_source_combination_breakdown(
        "Superlocus source combinations",
        superloci,
        superlocus_sources,
    )

    log_section("Sublocus summary")
    print(f"Subloci: {len(subloci):,}")
    sublocus_type_counts = Counter(sub.sublocus_type for sub in subloci)
    sublocus_model_counts = Counter(len(sub.models) for sub in subloci)
    sublocus_source_counts = Counter(len({m.source_name for m in sub.models}) for sub in subloci)
    print(f"Sublocus type counts: {format_counter(sublocus_type_counts)}")
    log_distribution("Model-count distribution", sublocus_model_counts, "models", "subloci")
    log_distribution("Source-count distribution", sublocus_source_counts, "sources", "subloci")
    log_subsection("Source-count breakdown by source presence")
    log_source_presence_breakdown(
        "Subloci containing each source within each source-count class",
        "subloci",
        subloci,
        lambda sub: {m.source_name for m in sub.models},
        all_sources,
    )
    log_subsection("Source-count breakdown by source combination")
    log_source_combination_breakdown(
        "Sublocus source combinations",
        subloci,
        lambda sub: {m.source_name for m in sub.models},
    )

    sublocus_flag_counts = Counter()
    for sub in subloci:
        if sub.flags != ".":
            for flag in sub.flags.split(","):
                sublocus_flag_counts[flag] += 1
    log_subsection("Sublocus and relationship diagnostics")
    print(f"Sublocus flag counts: {format_counter(sublocus_flag_counts)}")
    edge_type_counts = Counter(edge.edge_type for edge in relationship_edges)
    print(f"Relationship edge counts: {format_counter(edge_type_counts)}")

    if genome_fasta_path is not None:
        log_section("Genome sequence integrity")
        print(f"Genome FASTA: {genome_fasta_path}")
        print(
            f"Genome FASTA index: {genome_fasta_index_path}; action={genome_fasta_index_action}; "
            f"rebuild_reason={genome_fasta_index_rebuild_reason}"
        )
        print(
            "Canonical splice motifs: GT-AG, GC-AG, AT-AC; all other motifs, including ambiguous motifs, are noncanonical."
        )
        for line in sequence_integrity_log_messages:
            print(line)

    log_section("Intrinsic model metrics")
    metric_flag_counts = Counter()
    for metric in model_metrics:
        if metric.flags != ".":
            for flag in metric.flags.split(","):
                metric_flag_counts[flag] += 1
    print(f"Model intrinsic metric flag counts: {format_counter(metric_flag_counts)}")
    if long_cds_intron_cutoff is None:
        print(f"Long CDS intron metric cutoff: {long_cds_intron_cutoff_source}")
    else:
        print(
            f"Long CDS intron metric cutoff: > {long_cds_intron_cutoff:.3f} bp "
            f"({long_cds_intron_cutoff_source})"
        )
    print(
        "CDS exon metric cutoffs: "
        f"terminal < {min_terminal_cds_exon_bp} bp -> short_terminal_cds_exon; "
        f"internal < {min_internal_cds_exon_bp} bp -> short_internal_cds_exon"
    )
    print(
        "CDS intron metric cutoffs: "
        f"short < {args.min_cds_intron_bp} bp -> short_cds_intron; "
        f"long > {long_cds_intron_cutoff:.3f} bp -> long_cds_intron"
        if long_cds_intron_cutoff is not None
        else f"short < {args.min_cds_intron_bp} bp -> short_cds_intron; long intron flag not applied"
    )
    log_subsection("CDS intron length statistics by source")
    print("source\tn\tmin\tmedian\tp99\tp99.99\tp99.999\tmax")
    introns_by_source = collect_cds_intron_lengths_by_source(coding_models)
    all_introns_for_log: list[int] = []
    for source in all_sources:
        source_introns = introns_by_source.get(source, [])
        all_introns_for_log.extend(source_introns)
        stats_text = format_intron_stats(source_introns)
        stats_map = {}
        if stats_text != "n=0":
            for part in stats_text.split("; "):
                key, value = part.split("=", 1)
                stats_map[key] = value
            print(
                f"{source}\t{stats_map.get('n', '.')}\t{stats_map.get('min', '.')}\t"
                f"{stats_map.get('median', '.')}\t{stats_map.get('p99', '.')}\t"
                f"{stats_map.get('p99.99', '.')}\t{stats_map.get('p99.999', '.')}\t"
                f"{stats_map.get('max', '.')}"
            )
        else:
            print(f"{source}\t0\t.\t.\t.\t.\t.\t.")
    stats_text = format_intron_stats(all_introns_for_log)
    stats_map = {}
    if stats_text != "n=0":
        for part in stats_text.split("; "):
            key, value = part.split("=", 1)
            stats_map[key] = value
        print(
            f"ALL\t{stats_map.get('n', '.')}\t{stats_map.get('min', '.')}\t"
            f"{stats_map.get('median', '.')}\t{stats_map.get('p99', '.')}\t"
            f"{stats_map.get('p99.99', '.')}\t{stats_map.get('p99.999', '.')}\t"
            f"{stats_map.get('max', '.')}"
        )
    else:
        print("ALL\t0\t.\t.\t.\t.\t.\t.")


    log_section("Structural support metrics")
    print(
        "Structural support is calculated within first-pass subloci using other-source support fractions. "
        "Raw counts are retained in the model structural support TSV."
    )
    print(
        f"Terminal 3' CDS exon wobble: {args.terminal_cds_3prime_wobble_bp} bp; "
        f"single-exon support reciprocal overlap threshold: {single_exon_support_min_reciprocal_overlap:.3f}"
    )
    print(
        "Multi-exon preliminary score weights: "
        + ", ".join(f"{name}={weight:.2f}" for name, weight in STRUCTURAL_MULTI_SCORE_WEIGHTS.items())
    )
    print(
        "Single-exon preliminary score weights: "
        + ", ".join(f"{name}={weight:.2f}" for name, weight in STRUCTURAL_SINGLE_SCORE_WEIGHTS.items())
    )
    print(
        "Multi-exon consensus similarity weights: "
        + ", ".join(f"{name}={weight:.2f}" for name, weight in CONSENSUS_MULTI_SCORE_WEIGHTS.items())
    )
    print(
        "Single-exon consensus similarity weights: "
        + ", ".join(f"{name}={weight:.2f}" for name, weight in CONSENSUS_SINGLE_SCORE_WEIGHTS.items())
    )
    print(
        "Pairwise similarity is source-collapsed: for each other source, the best matching model is selected by weighted pairwise F1, then best hits are averaged across sources."
    )
    print(
        "Completeness penalty weights: "
        + ", ".join(f"{name}={weight:.2f}" for name, weight in COMPLETENESS_PENALTY_WEIGHTS.items())
    )
    print(
        "Completeness diagnostics classify missing supported CDS bases relative to the query model, and classify extra non-consensus CDS bases as terminal versus internal relative to the consensus CDS span."
    )
    print(
        "Intrinsic penalty weights: "
        + ", ".join(f"{name}={weight:.2f}" for name, weight in intrinsic_penalty_weights.items())
    )
    unsupported_short_penalty_models = sum(
        1 for row in model_structural_rows
        if value_as_float(row.get("unsupported_short_cds_exon_penalty_score")) > 0
    )
    unsupported_short_5prime_models = sum(
        1 for row in model_structural_rows
        if value_as_float(row.get("unsupported_short_5prime_terminal_cds_exon_count")) > 0
    )
    unsupported_short_3prime_models = sum(
        1 for row in model_structural_rows
        if value_as_float(row.get("unsupported_short_3prime_terminal_cds_exon_count")) > 0
    )
    unsupported_short_internal_models = sum(
        1 for row in model_structural_rows
        if value_as_float(row.get("unsupported_moderately_short_internal_cds_exon_count")) > 0
    )
    print(
        "Unsupported moderately short CDS exon selection penalty: "
        f"5prime interval=[{min_terminal_cds_exon_bp},{args.unsupported_short_5prime_terminal_cds_exon_bp}) bp weight={args.unsupported_short_5prime_terminal_weight:.2f}; "
        f"3prime interval=[{min_terminal_cds_exon_bp},{args.unsupported_short_3prime_terminal_cds_exon_bp}) bp weight={args.unsupported_short_3prime_terminal_weight:.2f}; "
        f"internal interval=[{min_internal_cds_exon_bp},{args.unsupported_short_internal_cds_exon_bp}) bp weight={args.unsupported_short_internal_weight:.2f}; "
        f"selection_penalty={args.selection_penalty_unsupported_short_cds_exon:.2f}; "
        f"affected_models={unsupported_short_penalty_models:,} "
        f"(5prime={unsupported_short_5prime_models:,},3prime={unsupported_short_3prime_models:,},internal={unsupported_short_internal_models:,})."
    )
    coherence_eligible_regions = sum(
        str(row.get("penalty_gate_passed", "false")) == "true"
        for row in structural_coherence_region_rows
    )
    coherence_penalized_models = sum(
        value_as_float(row.get("structural_coherence_penalty_score")) > 0
        for row in model_structural_rows
    )
    coherence_max_raw = max(
        (value_as_float(row.get("structural_coherence_raw_score")) for row in model_structural_rows),
        default=0.0,
    )
    print(
        "Candidate structural-coherence penalty: "
        "gate=S>=2,S>C,left_anchor_independent_sources>=2,right_anchor_independent_sources>=2; "
        "raw=max((S-C)/(S+C)*min(1.5,1+0.10*(S-2))); "
        f"selection_penalty={args.selection_penalty_structural_coherence:.2f}; "
        f"eligible_regions={coherence_eligible_regions:,}; affected_models={coherence_penalized_models:,}; "
        f"max_raw={coherence_max_raw:.3f}."
    )
    if genome_fasta_path is not None:
        print(
            "Genome sequence-integrity model-only penalty: "
            f"effective_weight={effective_sequence_integrity_penalty_weight:.2f}; "
            f"internal_stop={args.sequence_penalty_internal_stop:.2f} per stop cap={args.sequence_internal_stop_count_cap}; "
            f"phase_inconsistent={args.sequence_penalty_phase_inconsistent:.2f}; "
            f"missing_stop={args.sequence_penalty_missing_stop:.2f}; missing_start={args.sequence_penalty_missing_start:.2f}; "
            f"noncanonical_splice_fraction={args.sequence_penalty_noncanonical_splice:.2f}; "
            f"diagnostic_only={'true' if args.no_sequence_integrity_scoring else 'false'}."
        )
    print(
        "Model selection weights: "
        + ", ".join(f"{name}={weight:.2f}" for name, weight in model_selection_weights.items())
        + f"; tie_epsilon={args.model_selection_tie_epsilon:g}"
    )
    print(
        "Model selection ranking: score_rank groups sublocus scores within tie_epsilon using competition ranks; "
        "order_rank reports the exact shared primary-selector total order before later suppression, backfill, and pruning."
    )
    print(
        "Low-source conflict dampening: "
        + ("enabled" if args.low_source_conflict_dampen_agreement_scores else "disabled")
        + f"; criteria=sublocus_source_count==2, sublocus_model_count==2, exact_structure_support_score==0, pairwise_similarity_score<{LOW_SOURCE_CONFLICT_PAIRWISE_THRESHOLD:.2f}; "
        + f"agreement_factor={LOW_SOURCE_CONFLICT_AGREEMENT_FACTOR:.2f} applied only to structural_support_score and consensus_similarity_score."
    )
    protein_evidence_active_for_selection = bool(protein_evidence_paths) and args.selection_weight_protein_evidence > 0
    transcript_evidence_active_for_selection = bool(transcript_evidence_paths) and args.selection_weight_transcript_evidence > 0
    model_only_scale_max = model_selection_score_positive_scale_max(model_selection_weights, include_protein_evidence=False, include_transcript_evidence=False)
    final_scale_max = model_selection_score_positive_scale_max(
        model_selection_weights,
        include_protein_evidence=protein_evidence_active_for_selection,
        include_transcript_evidence=transcript_evidence_active_for_selection,
    )
    print(
        "Model selection scaled score: model_structural_support_score_scaled is scaled by the model-only theoretical positive maximum "
        f"{model_only_scale_max:.3f}; model_selection_score_scaled is scaled by the final theoretical positive maximum {final_scale_max:.3f}. "
        "Raw scores <= 0 map to 0 and values above the maximum are clamped to 1."
    )
    structural_score_type_counts = Counter(str(row.get("structural_support_score_type", ".")) for row in model_structural_rows)
    print(f"Structural support rows: {len(model_structural_rows):,}; feature rows: {len(structural_feature_rows):,}")
    print(f"Structural score type counts: {format_counter(structural_score_type_counts)}")
    feature_type_counts = Counter(str(row.get("feature_type", ".")) for row in structural_feature_rows)
    print(f"Structural feature type counts: {format_counter(feature_type_counts)}")
    tied_model_count = sum(1 for row in model_structural_rows if str(row.get("model_selection_tie_group", ".")) != ".")
    top_rank_count = sum(1 for row in model_structural_rows if str(row.get("model_selection_score_rank", ".")) == "1")
    low_source_conflict_count = sum(1 for row in model_structural_rows if str(row.get("low_source_conflict_flag", ".")) == "true")
    low_source_dampened_count = sum(1 for row in model_structural_rows if value_as_float(row.get("low_source_conflict_agreement_factor"), 1.0) < 1.0)
    print(f"Model selection ranking: top_score_rank_models={top_rank_count:,}; tied_models={tied_model_count:,}")
    print(f"Low-source conflict rows: {low_source_conflict_count:,}; dampened rows: {low_source_dampened_count:,}")
    evidence_rank_changed_rows = sum(1 for row in model_structural_rows if str(row.get("rank_changed_by_evidence", ".")) == "true")
    evidence_top_changed_rows = sum(1 for row in model_structural_rows if str(row.get("evidence_changed_top_rank", ".")) == "true")
    evidence_top_changed_subloci = len({str(row.get("sublocus_id", ".")) for row in model_structural_rows if str(row.get("evidence_changed_top_rank", ".")) == "true"})
    print(
        "Evidence rank-change diagnostics: "
        f"rows_rank_changed={evidence_rank_changed_rows:,}; "
        f"rows_evidence_changed_top_rank={evidence_top_changed_rows:,}; "
        f"subloci_evidence_changed_top_rank={evidence_top_changed_subloci:,}"
    )

    if protein_evidence_paths and protein_summary:
        log_section("Protein evidence diagnostics")
        print(
            "Protein evidence selection integration: "
            + (
                f"active with selection_weight_protein_evidence={args.selection_weight_protein_evidence:.3f}; "
                "protein_evidence_score contributes to model_selection_score and final ranks."
                if protein_evidence_active_for_selection
                else "inactive; protein evidence is parsed, collapsed, grouped, and scored for diagnostics but does not alter selection."
            )
        )
        print(
            "Protein evidence inputs: "
            f"files={int(protein_summary['file_count']):,}; "
            f"sources={int(protein_summary['source_count']):,}; "
            f"source_labels={protein_summary['source_labels']}"
        )
        print(
            "Protein evidence parsing/collapse: "
            f"raw_alignments={int(protein_summary['raw_alignment_count']):,}; "
            f"collapsed_alignments={int(protein_summary['collapsed_alignment_count']):,}; "
            f"collapse_ratio={float(protein_summary['collapse_ratio']):.3f}; "
            f"raw_to_collapsed_mapping_requested={str(protein_summary['raw_to_collapsed_mapping_requested']).lower()}"
        )
        print(
            "Protein evidence settings: "
            f"terminal_cds_3prime_wobble_bp={protein_evidence_wobble_bp:,}; "
            f"support_score_threshold={args.protein_evidence_support_score_threshold:.3f}; "
            f"single_exon_reciprocal_overlap={args.single_exon_min_reciprocal_overlap:.3f}"
        )
        print(
            "Protein evidence loci: "
            f"superloci={int(protein_summary['superlocus_count']):,}; "
            f"subloci={int(protein_summary['sublocus_count']):,}; "
            f"single_exon_subloci={int(protein_summary['single_exon_sublocus_count']):,}; "
            f"multi_exon_subloci={int(protein_summary['multi_exon_sublocus_count']):,}; "
            f"mixed_subloci={int(protein_summary['mixed_sublocus_count']):,}"
        )
        print(
            "Protein evidence model associations: "
            f"models={int(protein_summary['model_support_row_count']):,}; "
            f"with_overlap={int(protein_summary['models_with_overlapping_protein_evidence']):,}; "
            f"exact_chain={int(protein_summary['models_with_exact_chain_support']):,}; "
            f"supporting_source={int(protein_summary['models_with_supporting_protein_source']):,}; "
            f"multi_source_support={int(protein_summary['models_with_multi_source_support']):,}; "
            f"strong={int(protein_summary['models_with_strong_protein_support']):,}; "
            f"partial={int(protein_summary['models_with_partial_protein_support']):,}; "
            f"weak={int(protein_summary['models_with_weak_protein_support']):,}; "
            f"no_overlap={int(protein_summary['models_with_no_overlapping_protein_evidence']):,}"
        )
        print(f"Protein evidence model status counts: {protein_summary['model_status_counts']}")
        print(f"Protein evidence model flag counts: {protein_summary['model_flag_counts']}")
        print(
            "Protein evidence per-source summaries: "
            f"raw_alignments({protein_summary['per_source_raw_alignment_summary']}); "
            f"collapsed_alignments({protein_summary['per_source_collapsed_alignment_summary']}); "
            f"raw_members_in_collapsed_alignments({protein_summary['per_source_member_count_summary']})"
        )

    if transcript_evidence_paths and transcript_summary:
        log_section("Transcript evidence diagnostics")
        print(
            "Transcript evidence selection integration: "
            + (
                f"active with selection_weight_transcript_evidence={args.selection_weight_transcript_evidence:.3f}; "
                "transcript_evidence_score contributes to model_selection_score and final ranks."
                if transcript_evidence_active_for_selection
                else "inactive; transcript evidence is parsed, collapsed, grouped, associated, and scored for diagnostics but does not alter selection."
            )
        )
        print(
            "Transcript evidence inputs: "
            f"files={int(transcript_summary['file_count']):,}; "
            f"sources={int(transcript_summary['source_count']):,}; "
            f"source_labels={transcript_summary['source_labels']}"
        )
        print(
            "Transcript evidence parsing/collapse: "
            f"raw_transcripts={int(transcript_summary['raw_transcript_count']):,}; "
            f"collapsed_transcripts={int(transcript_summary['collapsed_transcript_count']):,}; "
            f"collapse_ratio={float(transcript_summary['collapse_ratio']):.3f}; "
            f"raw_to_collapsed_mapping_requested={str(transcript_summary['raw_to_collapsed_mapping_requested']).lower()}"
        )
        print(
            "Transcript evidence loci: "
            f"superloci={int(transcript_summary['superlocus_count']):,}; "
            f"subloci={int(transcript_summary['sublocus_count']):,}; "
            f"single_exon_subloci={int(transcript_summary['single_exon_sublocus_count']):,}; "
            f"multi_exon_subloci={int(transcript_summary['multi_exon_sublocus_count']):,}; "
            f"mixed_subloci={int(transcript_summary['mixed_sublocus_count']):,}"
        )
        print(
            "Transcript evidence collapsed records: "
            f"single_exon={int(transcript_summary['single_exon_collapsed_transcript_count']):,}; "
            f"multi_exon={int(transcript_summary['multi_exon_collapsed_transcript_count']):,}; "
            f"collapse_types={transcript_summary['collapse_type_counts']}; "
            f"flags={transcript_summary['collapsed_flag_counts']}"
        )
        print(
            "Transcript evidence unstranded records: "
            f"raw={int(transcript_summary['unstranded_raw_transcript_count']):,}; "
            f"raw_single_exon={int(transcript_summary['unstranded_single_exon_raw_transcript_count']):,}; "
            f"raw_multi_exon={int(transcript_summary['unstranded_multi_exon_raw_transcript_count']):,}; "
            f"collapsed={int(transcript_summary['unstranded_collapsed_transcript_count']):,}; "
            f"collapsed_single_exon={int(transcript_summary['unstranded_single_exon_collapsed_transcript_count']):,}; "
            f"collapsed_multi_exon={int(transcript_summary['unstranded_multi_exon_collapsed_transcript_count']):,}"
        )
        print(
            "Transcript evidence model associations: "
            f"models={int(transcript_model_summary.get('model_support_row_count', 0)):,}; "
            f"broad_association={int(transcript_model_summary.get('models_with_broad_transcript_association', 0)):,}; "
            f"scoring_eligible={int(transcript_model_summary.get('models_with_scoring_eligible_transcript_evidence', 0)):,}; "
            f"shared_intron={int(transcript_model_summary.get('models_with_shared_intron_transcript_evidence', 0)):,}; "
            f"shared_splice_site={int(transcript_model_summary.get('models_with_shared_splice_site_transcript_evidence', 0)):,}; "
            f"broad_only={int(transcript_model_summary.get('models_with_broad_only_transcript_association', 0)):,}; "
            f"score_gt_0={int(transcript_model_summary.get('models_with_transcript_evidence_score_gt_0', 0)):,}; "
            f"score_ge_0.75={int(transcript_model_summary.get('models_with_transcript_evidence_score_scaled_ge_0.75', 0)):,}"
        )
        print(f"Transcript evidence model status counts: {transcript_model_summary.get('model_status_counts', '.')}")
        print(f"Transcript evidence model flag counts: {transcript_model_summary.get('model_flag_counts', '.')}")
        print(
            "Transcript evidence association classes: "
            f"all={transcript_model_summary.get('association_class_counts', '.')}; "
            f"scoring={transcript_model_summary.get('scoring_association_class_counts', '.')}"
        )
        print(
            "Transcript evidence unstranded model associations: "
            f"models_broad={int(transcript_model_summary.get('models_with_unstranded_broad_transcript_association', 0)):,}; "
            f"models_scoring={int(transcript_model_summary.get('models_with_unstranded_scoring_eligible_transcript_evidence', 0)):,}; "
            f"models_only_unstranded_scoring={int(transcript_model_summary.get('models_with_only_unstranded_scoring_eligible_transcript_evidence', 0)):,}; "
            f"unstranded_broad_records={int(transcript_model_summary.get('unstranded_broad_association_total', 0)):,}; "
            f"unstranded_scoring_records={int(transcript_model_summary.get('unstranded_scoring_eligible_total', 0)):,}; "
            f"stranded_scoring_records={int(transcript_model_summary.get('stranded_scoring_eligible_total', 0)):,}"
        )
        print(
            "Transcript cDNA-context diagnostics: "
            f"calculated={int(transcript_model_summary.get('models_with_cdna_context_calculated', 0)):,}; "
            f"low_cds_fraction={int(transcript_model_summary.get('models_with_cdna_context_low_cds_fraction', 0)):,}; "
            f"very_low_cds_fraction={int(transcript_model_summary.get('models_with_cdna_context_very_low_cds_fraction', 0)):,}; "
            f"pruning_concern={int(transcript_model_summary.get('models_with_cdna_context_pruning_concern', 0)):,}; "
            f"severe_concern={int(transcript_model_summary.get('models_with_cdna_context_severe_concern', 0)):,}; "
            f"moderate_structural_concern={int(transcript_model_summary.get('models_with_cdna_context_moderate_structural_concern', 0)):,}; "
            f"context_best_not_pairwise_best={int(transcript_model_summary.get('models_with_cdna_context_not_best_pairwise', 0)):,}"
        )
        print(
            "Transcript evidence per-source summaries: "
            f"raw_transcripts({transcript_summary['per_source_raw_transcript_summary']}); "
            f"collapsed_transcripts({transcript_summary['per_source_collapsed_transcript_summary']}); "
            f"raw_members_in_collapsed_transcripts({transcript_summary['per_source_member_count_summary']})"
        )

    if args.prune_singleton_subloci_from_scored_gff3:
        pruned_sublocus_count, pruned_model_count, pruned_source_counts = count_pruned_single_model_subloci(
            subloci,
            scored_gff3_prune_singleton_sources,
        )
        log_section("Scored-candidate GFF3 pruning")
        source_scope = ",".join(sorted(scored_gff3_prune_singleton_sources)) if scored_gff3_prune_singleton_sources else "ALL"
        print(f"Scored-candidate pruned GFF3 enabled: {scored_pruned_gff3}")
        print(f"Singleton sublocus source scope: {source_scope}")
        print(f"Pruned singleton subloci: {pruned_sublocus_count:,}")
        print(f"Pruned candidate models: {pruned_model_count:,}")
        print(f"Pruned singleton subloci by source: {format_counter(pruned_source_counts)}")
        print("Only the additional scored-candidate pruned GFF3 is filtered; TSV scoring reports, selected-model outputs, and the main combined GFF3 remain unfiltered.")

    if args.run_monosublocus_selection and stage4_result is not None:
        log_section("Monosublocus de-overlap")
        mono_rows = stage4_result["monosublocus_rows"]  # type: ignore[index]
        edge_rows = stage4_result["monosublocus_edge_rows"]  # type: ignore[index]
        selected_rows = stage4_result["selected_rows"]  # type: ignore[index]
        suppressed_rows = stage4_result["suppressed_rows"]  # type: ignore[index]
        print(
            "Sublocus winners were clustered into monosublocus holders using non-reciprocal CDS-exon overlap; "
            "overlap is measured over CDS exons, not CDS span."
        )
        print(f"Monosublocus CDS-exon overlap threshold: {args.monosublocus_min_cds_exon_overlap_fraction:.3f}")
        print(f"Sublocus winners: {len(stage4_result['sublocus_winners']):,}")
        print(f"Monosublocus holders: {len(mono_rows):,}; overlap edges: {len(edge_rows):,}")
        print(f"Final selected models: {len(selected_rows):,}; suppressed winners: {len(suppressed_rows):,}")
        chain_like = sum(1 for row in mono_rows if str(row.get("not_all_pairwise_overlapping", "false")) == "true")
        print(f"Non-all-pairwise / chain-like monosubloci: {chain_like:,}")
    else:
        log_section("Monosublocus de-overlap")
        print("Not run because --stop-after scored-candidates was used.")

    if args.run_backfill_selection and stage5_result is not None:
        log_section("Backfill rescue")
        candidate_rows = stage5_result["backfill_candidate_rows"]  # type: ignore[index]
        backfilled_rows = stage5_result["backfilled_rows"]  # type: ignore[index]
        rejected_rows = stage5_result["backfill_rejected_rows"]  # type: ignore[index]
        final_rows = stage5_result["final_rows"]  # type: ignore[index]
        reason_counts = Counter(str(row.get("backfill_rejection_reason", ".")) for row in rejected_rows)
        print(
            "Backfill scanned non-selected, non-suppressed scored candidate models within each superlocus. "
            "Candidates were added greedily by model-selection rank only if they did not share CDS splice sites "
            "with selected multi-exon models and did not exceed the CDS-exon overlap threshold against any selected model."
        )
        print(f"Backfill CDS-exon overlap threshold: {args.monosublocus_min_cds_exon_overlap_fraction:.3f}")
        print(f"Backfill candidates considered: {len(candidate_rows):,}")
        print(f"Backfilled models accepted: {len(backfilled_rows):,}; rejected: {len(rejected_rows):,}")
        print(f"Final selected models after backfill: {len(final_rows):,}")
        print(f"Backfill rejection reasons: {format_counter(reason_counts)}")
        if args.strict_r1_backfill_replacement:
            r1_rows = stage5_result.get("r1_replacement_rows", [])
            r1_summary = {
                str(row.get("metric")): row.get("value")
                for row in stage5_result.get("r1_summary_rows", [])  # type: ignore[union-attr]
                if str(row.get("metric")) != "rejections"
            }
            print(
                "Strict R1 backfill replacement: enabled; "
                f"allow_multiexon_blockers={'true' if args.strict_r1_allow_multiexon_blockers else 'false'}; "
                f"candidates_evaluated={int(value_as_float(r1_summary.get('candidates_evaluated'))):,}; "
                f"accepted={len(r1_rows):,}; "
                f"one_blocker={int(value_as_float(r1_summary.get('one_blocker_replacements'))):,}; "
                f"two_blocker={int(value_as_float(r1_summary.get('two_blocker_replacements'))):,}; "
                f"corroboration_weight={args.selection_weight_structural_corroboration:.3f}; "
                f"alpha={args.structural_corroboration_alpha:g}."
            )
    elif args.run_monosublocus_selection:
        log_section("Backfill rescue")
        print("Not run because --stop-after monosublocus was used.")

    if automatic_assessment_requested and final_prune_result is not None:
        log_section("Automatic pruning assessment")
        scenario_map = final_prune_result.get("auto_scenarios", {})
        validation = final_prune_result.get("auto_scenario_validation")
        standard_decision = final_prune_result.get("standard_auto_decision")
        conservative_decision = final_prune_result.get("conservative_auto_decision")
        consensus_metrics = final_prune_result.get("auto_consensus_metrics")
        flood_assessment = final_prune_result.get(
            "auto_flood_consensus_assessment"
        )
        input_metrics = final_prune_result.get("auto_input_count_metrics", {})
        print(
            f"Requested mode: {args.pruning_level}; "
            f"effective output level: {final_prune_result.get('effective_pruning_level', args.pruning_level)}"
        )
        if scenario_map:
            print("Prospective scenarios:")
            for scenario_level in AUTO_SCENARIO_LEVELS:
                scenario = scenario_map[scenario_level]
                print(
                    f"  {scenario_level}: remove {scenario.removed_models:,}; "
                    f"retain {scenario.retained_models:,}; removed_fraction={scenario.removed_fraction:.4f}; "
                    f"native_removed={len(scenario.native_removed_ids):,}; "
                    f"inherited_removed={len(scenario.inherited_removed_ids):,}; "
                    f"additional={scenario.incremental_removals:,}; "
                    f"additional_fraction_previous_retained="
                    f"{scenario.incremental_fraction_previous_retained:.4f}; "
                    f"thresholds={scenario.thresholds}"
                )
            high_scenario = scenario_map["high"]
            print(
                "High versus medium: "
                f"additional_models_removed={high_scenario.incremental_removals:,}; "
                f"fraction_starting_selected_models="
                f"{high_scenario.incremental_fraction_original:.4f}; "
                f"additional_removals_by_source="
                f"{format_counter(Counter(high_scenario.incremental_by_source))}"
            )
        if consensus_metrics is not None:
            print(
                "Consensus backbone: "
                f"sources={consensus_metrics.effective_source_count}; "
                f"unadjusted={consensus_metrics.backbone_class}; "
                f"deep_consensus_incidence="
                f"{consensus_metrics.deep_consensus_incidence_fraction:.4f}"
            )
        if flood_assessment is not None:
            adjusted = flood_assessment.adjusted
            if adjusted is not None:
                for line in flood_adjusted_consensus_log_lines(
                    flood_assessment,
                    requested_mode=args.pruning_level,
                    standard_decision=(
                        standard_decision
                        if isinstance(standard_decision, AutoDecision)
                        else None
                    ),
                    conservative_decision=(
                        conservative_decision
                        if isinstance(conservative_decision, AutoDecision)
                        else None
                    ),
                ):
                    print(line)
            else:
                print(
                    "Flood-aware backbone assessment: "
                    f"enabled={str(flood_assessment.adjustment_enabled).lower()}; "
                    f"qualifying_source="
                    f"{flood_assessment.qualifying_excess_source or '.'}; "
                    "adjustment_applied=false"
                )
        if input_metrics:
            run_median = next(iter(input_metrics.values())).run_median_input_count
            print(
                "Input coding-model count comparisons (reporting only): "
                f"run_median={run_median:.1f}"
            )
            for source, metric in sorted(input_metrics.items()):
                run_ratio_text = (
                    "."
                    if metric.ratio_to_run_median is None
                    else f"{metric.ratio_to_run_median:.2f}"
                )
                legacy_ratio_text = (
                    "inf"
                    if math.isinf(metric.ratio_to_other_source_median)
                    else f"{metric.ratio_to_other_source_median:.2f}"
                )
                print(
                    f"  {source}: models={metric.input_coding_model_count:,}; "
                    f"ratio_to_run_median={run_ratio_text}; "
                    f"legacy_other_source_median={metric.other_source_median_count:.1f}; "
                    f"legacy_ratio={legacy_ratio_text}; warning={metric.warning}"
                )
                if metric.warning == "INPUT_CODING_MODEL_COUNT_HIGH":
                    print(
                        f"  Warning: Source {source} has {run_ratio_text}x the "
                        "run's median input coding-model count. This may indicate "
                        "over-prediction or differences among input annotations. "
                        "Review the input annotation and pruning diagnostics."
                    )
                elif metric.warning == "INPUT_CODING_MODEL_COUNT_LOW":
                    print(
                        f"  Warning: Source {source} has {run_ratio_text}x the "
                        "run's median input coding-model count. Check whether the "
                        "prediction file is complete and whether the upstream "
                        "prediction job finished successfully."
                    )
        if isinstance(validation, ScenarioValidation):
            print(
                "Scenario checks: "
                f"reconciliation={'pass' if validation.reconciliation_passed else 'fail'}; "
                f"medium_to_high_monotonicity="
                f"{'pass' if validation.medium_high_monotonicity_passed else 'fail'}; "
                f"high_to_very_high_monotonicity="
                f"{'pass' if validation.high_very_high_monotonicity_passed else 'fail'}"
            )
            if (
                validation.native_restored_medium_to_high
                or validation.native_restored_high_to_very_high
            ):
                print(
                    "Monotonic closure: "
                    f"native_medium_to_high_restorations="
                    f"{validation.native_restored_medium_to_high:,}; "
                    f"native_high_to_very_high_restorations="
                    f"{validation.native_restored_high_to_very_high:,}; "
                    "lower-level removal provenance inherited."
                )
        for label, decision in (
            ("auto", standard_decision),
            ("auto-conservative", conservative_decision),
        ):
            if not isinstance(decision, AutoDecision):
                continue
            print(
                f"{label}: unconstrained={decision.unconstrained_level}; "
                f"effective={decision.effective_level}; confidence={decision.confidence}; "
                f"constraints={decision.user_min_level}..{decision.user_max_level}"
            )
            print(f"  reasons: {','.join(decision.reason_codes) or '.'}")
            if consensus_metrics is not None:
                print(
                    "  explanation: "
                    + explain_auto_pruning_decision(
                        decision,
                        consensus=consensus_metrics,
                        flood_assessment=flood_assessment,
                    )
                )
            if decision.warning_codes:
                print(f"  warnings: {','.join(decision.warning_codes)}")

    if final_selected_pruning_requested and final_prune_result is not None:
        log_section("Final selected-model pruning")
        retained_rows = final_prune_result["retained_rows"]  # type: ignore[index]
        removed_rows = final_prune_result["removed_rows"]  # type: ignore[index]
        shadow_pair_rows = final_prune_result.get("shadow_pair_rows", [])  # type: ignore[assignment]
        class_counts: Counter[str] = Counter()
        for row in removed_rows:
            for cls in str(row.get("final_prune_classes", ".")).split(","):
                if cls and cls != ".":
                    class_counts[cls] += 1
        singleton_reason_counts = Counter(str(row.get("final_singleton_prune_reason", ".")) for row in removed_rows if row.get("final_singleton_pruned") == "true")
        shadow_reason_counts = Counter(str(row.get("final_shadow_prune_reason", ".")) for row in removed_rows if row.get("final_shadow_pruned") == "true")
        source_counts = Counter(str(row.get("source_name", ".")) for row in removed_rows)
        singleton_source_scope = ",".join(sorted(final_prune_singleton_subloci_sources)) if final_prune_singleton_subloci_sources else "ALL"
        shadow_source_scope = ",".join(sorted(shadow_prune_sources)) if shadow_prune_sources else "ALL"
        effective_level = str(
            final_prune_result.get("effective_pruning_level", args.pruning_level)
        )
        effective_config = dict(
            final_prune_result.get("effective_pruning_configuration", {})
        )
        singleton_thresholds = effective_config.get(
            "singleton_thresholds",
            (
                args.low_score_singleton_max_model_structural_support_score_scaled,
                args.low_score_singleton_max_protein_evidence_score_scaled,
                args.low_score_singleton_max_transcript_evidence_score_scaled,
                args.low_score_singleton_max_model_selection_score_scaled,
            ),
        )
        nonsingleton_thresholds = effective_config.get(
            "nonsingleton_thresholds",
            (
                args.low_score_nonsingleton_max_model_structural_support_score_scaled,
                args.low_score_nonsingleton_max_protein_evidence_score_scaled,
                args.low_score_nonsingleton_max_transcript_evidence_score_scaled,
                args.low_score_nonsingleton_max_model_selection_score_scaled,
            ),
        )
        print(
            "Final selected pruning is an output filter on final_selected_models after monosublocus de-overlap and backfill rescue. "
            "It does not alter the unpruned final_selected_models files."
        )
        print(
            f"Pruning level: requested={args.pruning_level}; "
            f"effective={effective_level}"
        )
        effective_low_score_singletons = bool(
            effective_config.get("low_score_singletons", args.final_prune_low_score_singletons)
        )
        effective_low_score_nonsingletons = bool(
            effective_config.get("low_score_nonsingletons", args.final_prune_low_score_nonsingletons)
        )
        effective_low_score_nonsingleton_scope = str(
            effective_config.get("low_score_nonsingleton_scope", args.low_score_nonsingleton_scope)
        )
        effective_sequence_nonsingletons = bool(
            effective_config.get(
                "sequence_integrity_nonsingletons",
                args.final_prune_sequence_integrity_nonsingletons,
            )
        )
        effective_sequence_nonsingleton_scope = str(
            effective_config.get(
                "sequence_integrity_nonsingleton_scope",
                args.sequence_integrity_nonsingleton_scope,
            )
        )
        effective_singleton_level = str(
            effective_config.get(
                "singleton_level", args.final_prune_singleton_subloci_level
            )
        )
        effective_shadow_enabled = bool(
            effective_config.get("shadow_enabled", args.final_prune_opposite_strand_shadows)
        )
        effective_shadow_require_single_exon = bool(
            effective_config.get(
                "shadow_require_loser_single_exon",
                args.shadow_prune_require_loser_single_exon,
            )
        )
        low_score_singletons_text = "enabled" if effective_low_score_singletons else "disabled"
        if effective_low_score_nonsingletons:
            low_score_nonsingletons_text = f"enabled(scope={effective_low_score_nonsingleton_scope})"
        else:
            low_score_nonsingletons_text = "disabled"
        if effective_sequence_nonsingletons:
            sequence_nonsingletons_text = (
                f"enabled(scope={effective_sequence_nonsingleton_scope})"
            )
        else:
            sequence_nonsingletons_text = "disabled"
        print(
            "Effective pruning preset behaviour: "
            f"pruning_level={effective_level}; "
            f"genome_fasta={str(genome_fasta_path is not None).lower()}; "
            f"protein_evidence={str(bool(protein_evidence_paths)).lower()}; "
            f"transcript_evidence={str(bool(transcript_evidence_paths)).lower()}; "
            f"low_score_singletons={low_score_singletons_text}; "
            f"low_score_nonsingletons={low_score_nonsingletons_text}; "
            f"sequence_integrity_nonsingletons={sequence_nonsingletons_text}"
        )
        print(f"Singleton pruning level: {effective_singleton_level}; source scope: {singleton_source_scope}")
        print(
            "Singleton monosublocus-overlap rescue: "
            f"enabled={str(args.singleton_prune_rescue_monosublocus_overlap).lower()}; "
            f"evidence_policy={args.monosublocus_overlap_rescue_evidence_policy}; "
            f"min_holder_member_count={args.singleton_rescue_min_monosublocus_member_count:,}; "
            f"min_overlap_bp={args.singleton_rescue_min_cds_overlap_bp:,}; "
            f"min_self_overlap_fraction={args.singleton_rescue_min_self_overlap_fraction:.3f}; "
            f"min_reciprocal_overlap_fraction={args.singleton_rescue_min_reciprocal_overlap_fraction:.3f}; "
            "partner_source=different_required"
        )
        print(
            "Singleton protein-evidence rescue: "
            f"enabled={str(args.final_prune_rescue_protein_supported_singletons).lower()}; "
            f"min_protein_evidence_score_scaled={args.singleton_rescue_min_protein_evidence_score_scaled:.3f}; "
            f"min_exact_chain_source_count={args.singleton_rescue_min_protein_exact_chain_source_count:,}"
        )
        print(
            "Singleton transcript-evidence rescue: "
            f"enabled={str(args.final_prune_rescue_transcript_supported_singletons).lower()}; "
            f"min_transcript_evidence_score_scaled={args.singleton_rescue_min_transcript_evidence_score_scaled:.3f}; "
            f"min_supporting_source_count={args.singleton_rescue_min_transcript_supporting_source_count:,}; "
            f"min_scoring_eligible_count={args.singleton_rescue_min_transcript_scoring_eligible_count:,}; "
            f"require_stranded_support={str(args.transcript_singleton_rescue_require_stranded_support).lower()}; "
            f"block_on_cdna_context_concern={str(args.transcript_singleton_rescue_block_on_cdna_context_concern).lower()}"
        )
        print(
            "Transcript singleton-rescue cDNA context: "
            f"min_cds_to_min_cdna_fraction={args.singleton_rescue_min_transcript_cdna_context_fraction:.3f}; "
            "not_evaluable_context=neutral"
        )
        low_score_adequacy_passed = bool(final_prune_result.get("low_score_singleton_external_adequacy_passed", False))
        low_score_protein_adequacy_passed = bool(final_prune_result.get("low_score_singleton_protein_adequacy_passed", False))
        low_score_transcript_adequacy_passed = bool(final_prune_result.get("low_score_singleton_transcript_adequacy_passed", False))
        print(
            "Evidence-aware low-score singleton pruning: "
            f"enabled={str(effective_low_score_singletons).lower()}; "
            f"use_transcript_cdna_context={str(args.low_score_singleton_use_transcript_cdna_context).lower()}; "
            f"removed_score_route={int(final_prune_result.get('low_score_singleton_removed_score_route', 0)):,}; "
            f"removed_context_route={int(final_prune_result.get('low_score_singleton_removed_context_route', 0)):,}; "
            f"external_adequacy_passed={str(low_score_adequacy_passed).lower()}; "
            f"protein_supplied={str(bool(protein_evidence_paths)).lower()}; "
            f"protein_adequacy_required={str(args.low_score_singleton_require_protein_evidence_adequacy).lower()}; "
            f"protein_adequacy_passed={str(low_score_protein_adequacy_passed).lower()}; "
            f"protein_overlap_model_fraction={protein_overlap_model_fraction:.3f}; "
            f"protein_supported_model_fraction={protein_supported_model_fraction:.3f}; "
            f"min_protein_overlap_fraction={args.low_score_singleton_min_protein_overlap_model_fraction:.3f}; "
            f"min_protein_supported_fraction={args.low_score_singleton_min_protein_supported_model_fraction:.3f}; "
            f"transcript_supplied={str(bool(transcript_evidence_paths)).lower()}; "
            f"transcript_adequacy_required={str(args.low_score_singleton_require_transcript_evidence_adequacy).lower()}; "
            f"transcript_adequacy_passed={str(low_score_transcript_adequacy_passed).lower()}; "
            f"transcript_scoring_eligible_model_fraction={transcript_scoring_eligible_model_fraction:.3f}; "
            f"transcript_supported_model_fraction={transcript_supported_model_fraction:.3f}; "
            f"min_transcript_scoring_eligible_fraction={args.low_score_singleton_min_transcript_scoring_eligible_model_fraction:.3f}; "
            f"min_transcript_supported_fraction={args.low_score_singleton_min_transcript_supported_model_fraction:.3f}; "
            f"max_model_score_scaled={float(singleton_thresholds[0]):.3f}; "
            f"max_protein_score_scaled={float(singleton_thresholds[1]):.3f}; "
            f"max_transcript_score_scaled={float(singleton_thresholds[2]):.3f}; "
            f"max_model_selection_score_scaled={float(singleton_thresholds[3]):.3f}"
        )
        low_score_nonsingleton_adequacy_passed = bool(final_prune_result.get("low_score_nonsingleton_external_adequacy_passed", False))
        low_score_nonsingleton_protein_adequacy_passed = bool(final_prune_result.get("low_score_nonsingleton_protein_adequacy_passed", False))
        low_score_nonsingleton_transcript_adequacy_passed = bool(final_prune_result.get("low_score_nonsingleton_transcript_adequacy_passed", False))
        print(
            "Low-score nonsingleton pruning: "
            f"enabled={str(effective_low_score_nonsingletons).lower()}; "
            f"scope={effective_low_score_nonsingleton_scope}; "
            f"policy={final_prune_result.get('low_score_nonsingleton_policy', '.')}; "
            f"run_model_source_count={len(source_names):,}; "
            f"use_transcript_cdna_context={str(args.low_score_nonsingleton_use_transcript_cdna_context).lower()}; "
            f"removed_score_route={int(final_prune_result.get('low_score_nonsingleton_removed_score_route', 0)):,}; "
            f"removed_context_route={int(final_prune_result.get('low_score_nonsingleton_removed_context_route', 0)):,}; "
            f"fallback_tiers={format_counter(Counter(final_prune_result.get('low_score_nonsingleton_model_confidence_tier_counts', {})))}; "
            f"external_adequacy_passed={str(low_score_nonsingleton_adequacy_passed).lower()}; "
            f"protein_supplied={str(bool(protein_evidence_paths)).lower()}; "
            f"protein_adequacy_required={str(args.low_score_nonsingleton_require_protein_evidence_adequacy).lower()}; "
            f"protein_adequacy_passed={str(low_score_nonsingleton_protein_adequacy_passed).lower()}; "
            f"protein_overlap_model_fraction={protein_overlap_model_fraction:.3f}; "
            f"protein_supported_model_fraction={protein_supported_model_fraction:.3f}; "
            f"min_protein_overlap_fraction={args.low_score_nonsingleton_min_protein_overlap_model_fraction:.3f}; "
            f"min_protein_supported_fraction={args.low_score_nonsingleton_min_protein_supported_model_fraction:.3f}; "
            f"transcript_supplied={str(bool(transcript_evidence_paths)).lower()}; "
            f"transcript_adequacy_required={str(args.low_score_nonsingleton_require_transcript_evidence_adequacy).lower()}; "
            f"transcript_adequacy_passed={str(low_score_nonsingleton_transcript_adequacy_passed).lower()}; "
            f"transcript_scoring_eligible_model_fraction={transcript_scoring_eligible_model_fraction:.3f}; "
            f"transcript_supported_model_fraction={transcript_supported_model_fraction:.3f}; "
            f"min_transcript_scoring_eligible_fraction={args.low_score_nonsingleton_min_transcript_scoring_eligible_model_fraction:.3f}; "
            f"min_transcript_supported_fraction={args.low_score_nonsingleton_min_transcript_supported_model_fraction:.3f}; "
            f"max_model_score_scaled={float(nonsingleton_thresholds[0]):.3f}; "
            f"max_protein_score_scaled={float(nonsingleton_thresholds[1]):.3f}; "
            f"max_transcript_score_scaled={float(nonsingleton_thresholds[2]):.3f}; "
            f"max_model_selection_score_scaled={float(nonsingleton_thresholds[3]):.3f}"
        )
        print(
            "Sequence-integrity nonsingleton pruning: "
            f"enabled={str(effective_sequence_nonsingletons).lower()}; "
            f"scope={effective_sequence_nonsingleton_scope}; "
            "backfill_scope_concerns=internal_stop,missing_start,missing_stop,all_splices_noncanonical; "
            "all_scope_severe_concerns=internal_stop,missing_start+missing_stop,all_splices_noncanonical; "
            f"protein_protection={str(args.final_prune_rescue_protein_supported_singletons).lower()}; "
            f"transcript_protection={str(args.final_prune_rescue_transcript_supported_singletons).lower()}"
        )
        print(
            "Internal-stop evidence rescue policy: "
            f"allow_internal_stop_evidence_rescue={str(args.allow_internal_stop_evidence_rescue).lower()}; "
            "resolved_action="
            + (
                "allow_protein_transcript_and_exact_chain_override"
                if args.allow_internal_stop_evidence_rescue else
                "block_protein_transcript_and_exact_chain_override"
            )
        )
        print(
            "Opposite-strand shadow pruning: "
            f"enabled={str(effective_shadow_enabled).lower()}; source scope={shadow_source_scope}; "
            f"min_overlap_bp={args.shadow_prune_min_overlap_bp:,}; "
            f"min_loser_cds_overlap_fraction={args.shadow_prune_min_loser_cds_overlap_fraction:.3f}; "
            f"min_score_margin_scaled={args.shadow_prune_min_score_margin_scaled:.3f}; "
            f"require_loser_single_exon={str(effective_shadow_require_single_exon).lower()}"
        )
        print(f"Final selected retained after pruning: {len(retained_rows):,}; removed: {len(removed_rows):,}")
        inherited_removed_count = sum(
            str(row.get("final_monotonic_inherited_removal", "false")) == "true"
            for row in removed_rows
        )
        if inherited_removed_count:
            inherited_reason_counts = Counter(
                str(row.get("final_monotonic_inherited_original_primary_reason", "."))
                for row in removed_rows
                if str(row.get("final_monotonic_inherited_removal", "false")) == "true"
            )
            print(
                "Monotonic inherited removals: "
                f"count={inherited_removed_count:,}; "
                f"original_reasons={format_counter(inherited_reason_counts)}"
            )
        print(f"Final pruning removed by class: {format_counter(class_counts)}")
        print(f"Final singleton removed by reason: {format_counter(singleton_reason_counts)}")
        print(
            "Final sequence-integrity singleton pruning summary: "
            f"candidates={int(final_prune_result.get('sequence_integrity_singleton_candidate_count', 0)):,}; "
            f"removed={int(final_prune_result.get('sequence_integrity_singleton_pruned_count', 0)):,}; "
            f"rescued={int(final_prune_result.get('sequence_integrity_singleton_rescued_count', 0)):,}"
        )
        singleton_rescue_counts = Counter(
            str(row.get("final_singleton_rescue_reason", "."))
            for row in retained_rows
            if str(row.get("final_singleton_rescued", "false")) == "true"
        )
        print(f"Final singleton rescued by reason: {format_counter(singleton_rescue_counts)}")
        print(
            "Final monosublocus-overlap rescue evidence-policy summary: "
            f"policy={args.monosublocus_overlap_rescue_evidence_policy}; "
            f"overlap_bp_geometry_candidates={int(final_prune_result.get('singleton_overlap_rescue_overlap_bp_geometry_candidate_count', 0)):,}; "
            f"rejected_reciprocal_overlap={int(final_prune_result.get('singleton_overlap_rescue_reciprocal_rejected_count', 0)):,}; "
            f"models_with_same_source_edge={int(final_prune_result.get('singleton_overlap_rescue_models_with_same_source_edge_count', 0)):,}; "
            f"geometry_candidates={int(final_prune_result.get('singleton_overlap_rescue_geometry_candidate_count', 0)):,}; "
            f"rescued_strong_local_evidence={int(final_prune_result.get('singleton_overlap_rescued_by_strong_local_evidence_count', 0)):,}; "
            f"rescued_no_adequate_evidence={int(final_prune_result.get('singleton_overlap_rescued_without_adequate_evidence_count', 0)):,}; "
            f"denied_adequate_evidence_without_strong_local_support={int(final_prune_result.get('singleton_overlap_rescue_policy_denied_count', 0)):,}"
        )
        actual_protein_singleton_rescue_counts = Counter(
            str(row.get("final_protein_singleton_rescue_applied_reason", "."))
            for row in retained_rows
            if str(row.get("final_protein_singleton_rescue_applied", "false")) == "true"
        )
        eligible_protein_singleton_rescue_counts = Counter(
            str(row.get("final_protein_singleton_rescue_criteria_reason", "."))
            for row in retained_rows
            if str(row.get("final_protein_singleton_rescue_criteria_met", "false")) == "true"
        )
        print(f"Final singleton protein-evidence actual rescue criteria: {format_counter(actual_protein_singleton_rescue_counts)}")
        print(f"Retained models meeting protein singleton-rescue criteria: {format_counter(eligible_protein_singleton_rescue_counts)}")
        actual_transcript_singleton_rescue_counts = Counter(
            str(row.get("final_transcript_singleton_rescue_applied_reason", "."))
            for row in retained_rows
            if str(row.get("final_transcript_singleton_rescue_applied", "false")) == "true"
        )
        eligible_transcript_singleton_rescue_counts = Counter(
            str(row.get("final_transcript_singleton_rescue_criteria_reason", "."))
            for row in retained_rows
            if str(row.get("final_transcript_singleton_rescue_criteria_met", "false")) == "true"
        )
        transcript_rescue_support_type_counts = Counter(
            str(row.get("final_transcript_singleton_rescue_support_type", "."))
            for row in retained_rows
            if str(row.get("final_transcript_singleton_rescue_applied", "false")) == "true"
        )
        print(f"Final singleton transcript-evidence actual rescue criteria: {format_counter(actual_transcript_singleton_rescue_counts)}")
        print(f"Retained models meeting transcript singleton-rescue criteria: {format_counter(eligible_transcript_singleton_rescue_counts)}")
        print(f"Final singleton transcript-evidence rescue support type: {format_counter(transcript_rescue_support_type_counts)}")
        low_score_reason_counts = Counter(
            str(row.get("final_low_score_singleton_prune_reason", "."))
            for row in removed_rows
            if str(row.get("final_low_score_singleton_pruned", "false")) == "true"
        )
        print(f"Final low-score singleton removed by reason: {format_counter(low_score_reason_counts)}")
        if args.final_prune_low_score_singletons:
            print(
                "Final low-score singleton pruning summary: "
                f"candidates={int(final_prune_result.get('low_score_singleton_candidates', 0)):,}; "
                f"removed={int(final_prune_result.get('low_score_singleton_removed_count', 0)):,}; "
                f"skipped_protein_rescue={int(final_prune_result.get('low_score_singleton_skipped_protein_rescue', 0)):,}; "
                f"skipped_transcript_rescue={int(final_prune_result.get('low_score_singleton_skipped_transcript_rescue', 0)):,}; "
                f"skipped_adequacy={int(final_prune_result.get('low_score_singleton_skipped_adequacy', 0)):,}; "
                f"skipped_thresholds={int(final_prune_result.get('low_score_singleton_skipped_thresholds', 0)):,}"
            )
        low_score_nonsingleton_reason_counts = Counter(
            str(row.get("final_low_score_nonsingleton_prune_reason", "."))
            for row in removed_rows
            if str(row.get("final_low_score_nonsingleton_pruned", "false")) == "true"
        )
        print(f"Final low-score nonsingleton removed by reason: {format_counter(low_score_nonsingleton_reason_counts)}")
        if args.final_prune_low_score_nonsingletons:
            print(
                "Final low-score nonsingleton pruning summary: "
                f"scope={final_prune_result.get('low_score_nonsingleton_scope', args.low_score_nonsingleton_scope)}; "
                f"candidates={int(final_prune_result.get('low_score_nonsingleton_candidates', 0)):,}; "
                f"removed={int(final_prune_result.get('low_score_nonsingleton_removed_count', 0)):,}; "
                f"policy={final_prune_result.get('low_score_nonsingleton_policy', '.')}; "
                f"skipped_scope={int(final_prune_result.get('low_score_nonsingleton_skipped_scope', 0)):,}; "
                f"skipped_model_source_guard={int(final_prune_result.get('low_score_nonsingleton_skipped_model_source_guard', 0)):,}; "
                f"skipped_thresholds={int(final_prune_result.get('low_score_nonsingleton_skipped_thresholds', 0)):,}; "
                f"removed_from_backfill={int(final_prune_result.get('low_score_nonsingleton_removed_from_backfill', 0)):,}; "
                f"removed_from_nonbackfill={int(final_prune_result.get('low_score_nonsingleton_removed_from_nonbackfill', 0)):,}"
            )
        sequence_nonsingleton_reason_counts = Counter(
            str(row.get("final_sequence_integrity_nonsingleton_prune_reason", "."))
            for row in removed_rows
            if str(row.get("final_sequence_integrity_nonsingleton_pruned", "false")) == "true"
        )
        print(
            "Final sequence-integrity nonsingleton removed by reason: "
            f"{format_counter(sequence_nonsingleton_reason_counts)}"
        )
        if args.final_prune_sequence_integrity_nonsingletons:
            print(
                "Final sequence-integrity nonsingleton pruning summary: "
                f"scope={final_prune_result.get('sequence_integrity_nonsingleton_scope', args.sequence_integrity_nonsingleton_scope)}; "
                f"candidates={int(final_prune_result.get('sequence_integrity_nonsingleton_candidates', 0)):,}; "
                f"removed={int(final_prune_result.get('sequence_integrity_nonsingleton_removed_count', 0)):,}; "
                f"protected_by_protein={int(final_prune_result.get('sequence_integrity_nonsingleton_protected_by_protein', 0)):,}; "
                f"protected_by_transcript={int(final_prune_result.get('sequence_integrity_nonsingleton_protected_by_transcript', 0)):,}; "
                f"skipped_scope={int(final_prune_result.get('sequence_integrity_nonsingleton_skipped_scope', 0)):,}; "
                f"skipped_no_concern={int(final_prune_result.get('sequence_integrity_nonsingleton_skipped_no_concern', 0)):,}; "
                f"backfill_standard_candidates={int(final_prune_result.get('sequence_integrity_nonsingleton_backfill_standard_candidates', 0)):,}; "
                f"all_severe_candidates={int(final_prune_result.get('sequence_integrity_nonsingleton_all_severe_candidates', 0)):,}; "
                f"removed_backfill_standard_route={int(final_prune_result.get('sequence_integrity_nonsingleton_removed_backfill_standard_route', 0)):,}; "
                f"removed_all_severe_route={int(final_prune_result.get('sequence_integrity_nonsingleton_removed_all_severe_route', 0)):,}; "
                f"removed_from_backfill={int(final_prune_result.get('sequence_integrity_nonsingleton_removed_from_backfill', 0)):,}; "
                f"removed_from_nonbackfill={int(final_prune_result.get('sequence_integrity_nonsingleton_removed_from_nonbackfill', 0)):,}"
            )
        print(f"Final shadow removed by reason: {format_counter(shadow_reason_counts)}")
        print(f"Final pruning removed by source: {format_counter(source_counts)}")
        if args.final_prune_opposite_strand_shadows:
            pruned_pair_count = sum(1 for row in shadow_pair_rows if str(row.get("shadow_pruned")) == "true")
            before_singleton_pairs = int(final_prune_result.get("shadow_candidate_pairs_before_singleton", len(shadow_pair_rows)))
            after_singleton_pairs = int(final_prune_result.get("shadow_candidate_pairs_after_singleton", len(shadow_pair_rows)))
            skipped_after_singleton = int(final_prune_result.get(
                "shadow_candidate_pairs_skipped_after_singleton",
                max(0, before_singleton_pairs - after_singleton_pairs),
            ))
            log_subsection("Opposite-strand shadow pruning")
            print(f"Candidate pairs before singleton pruning: {before_singleton_pairs:,}")
            print(f"Candidate pairs after singleton pruning: {after_singleton_pairs:,}")
            print(f"Pairs skipped because one member was already removed: {skipped_after_singleton:,}")
            print(f"Pruned pairs: {pruned_pair_count:,}")
    else:
        log_section("Final selected-model pruning")
        print("Not run. Use --pruning-level low/medium/high/very-high/auto, or advanced final pruning options in --full-help, with --stop-after backfill/final.")

    if args.stop_after == "final":
        log_section("Alternative-splice selection")
        if (
            alternative_splice_selection_result is not None
            and alternative_splice_report_result is not None
            and public_final_annotation is not None
        ):
            console_summary = build_alternative_splice_console_summary(
                alternative_splice_selection_result,
                alternative_splice_report_result["run_rows"][0],
                public_final_annotation,
                alternatives_emitted=alternative_splice_enabled,
            )
            for line in alternative_splice_console_lines(console_summary):
                print(line)
        else:
            print("Alternative-splice augmentation: disabled")

        log_section("Post-selection UTR augmentation")
        if utr_augmentation_result is None:
            print("UTR augmentation: disabled")
        else:
            utr_summary = utr_augmentation_result.summary
            utr_limits = utr_augmentation_result.limits
            print(
                "Assignments: "
                f"both={utr_summary.both_assigned:,}; "
                f"5-prime-only={utr_summary.five_prime_only:,}; "
                f"3-prime-only={utr_summary.three_prime_only:,}; "
                f"neither={utr_summary.neither_assigned:,}"
            )
            print(
                "Evidence contribution: "
                f"transcript-only={utr_summary.transcript_only:,}; "
                f"prediction-only={utr_summary.prediction_only:,}; "
                f"both={utr_summary.both_evidence_kinds:,}"
            )
            print(
                "Learned/effective limits (bp): "
                f"5-prime exonic={utr_limits.five_prime_exonic_bp or 'disabled'}; "
                f"3-prime exonic={utr_limits.three_prime_exonic_bp or 'disabled'}; "
                f"5-prime genomic={utr_limits.five_prime_genomic_bp or 'unavailable'}; "
                f"3-prime genomic={utr_limits.three_prime_genomic_bp or 'unavailable'}"
            )
            print(
                "Endpoint decisions: "
                f"local-consensus={utr_summary.local_consensus_sides:,}; "
                f"diffuse={utr_summary.diffuse_sides:,}; "
                f"strong-alternative-structure={utr_summary.strong_competing_structure_sides:,}; "
                f"strong-multimodality={utr_summary.strongly_multimodal_sides:,}"
            )
            print(
                "Evidence catalogue: "
                f"transcript={utr_summary.transcript_evidence_records:,}; "
                f"prediction={utr_summary.prediction_evidence_records:,}; "
                f"estimated-bytes={utr_summary.catalogue_estimated_bytes:,}; "
                f"multi-gene-excluded={utr_summary.multi_gene_evidence_excluded:,}"
            )
            print(
                "Cross-final-gene structural safeguard: "
                f"authoritative-candidate-sides-rejected="
                f"{utr_summary.authoritative_candidate_sides_rejected:,}; "
                f"initial-UTR-vs-UTR-conflict-sides="
                f"{utr_summary.provisional_selected_sides_in_utr_conflicts:,}; "
                f"fallback-sides={utr_summary.conflict_fallback_sides:,}; "
                f"no-alternative-sides={utr_summary.conflict_no_alternative_sides:,}; "
                f"affected-genes={utr_summary.conflict_affected_genes:,}; "
                f"affected-transcripts={utr_summary.conflict_affected_transcripts:,}; "
                f"resolution-rounds={utr_summary.conflict_resolution_rounds:,}"
            )
            print(
                "Quality rejections: "
                + format_counter(Counter(dict(utr_summary.rejection_counts)))
            )

    log_section("Source prior diagnostics")
    print(
        "Source prior mode: "
        f"{args.source_prior_mode}; source_prior_scale={args.source_prior_scale:.3f}; "
        f"auto_weight_range={args.auto_source_weight_min:.3f}-{args.auto_source_weight_max:.3f}"
    )
    if len(source_names) < 3:
        print(
            "Automatic source weights are neutral at 1.0 because fewer than three sources were supplied; "
            "with --source-prior-mode auto this makes the source-prior component neutral for all sources."
        )
    print(
        "Automatic agreement uses exact, unique, source-level coordinate features: "
        "0.50*supported CDS splice-site fraction + 0.30*supported CDS intron fraction + "
        "0.20*supported internal CDS exon fraction, re-normalised when a class is absent."
    )
    log_subsection("Source agreement and effective weights")
    print(
        "source\tcoding_models\tagreement\tauto_weight\tauto_status\tmanual_weight\teffective_weight\t"
        "prior_score\tunsupported_splice_sites\tunsupported_introns\tunsupported_internal_exons\t"
        "consensus_splice_capture\tconsensus_intron_capture"
    )
    for item in source_auto_weights:
        print(
            f"{item.source_name}\t{item.coding_models:,}\t{item.auto_agreement_score:.3f}\t"
            f"{item.auto_source_weight:.3f}\t{item.auto_source_weight_status}\t{item.manual_source_weight:.3f}\t"
            f"{item.effective_source_weight:.3f}\t{item.source_prior_score:.3f}\t"
            f"{item.unsupported_cds_splice_site_fraction:.3f}\t"
            f"{item.unsupported_cds_intron_fraction:.3f}\t"
            f"{item.unsupported_internal_cds_exon_fraction:.3f}\t"
            f"{item.consensus_cds_splice_site_capture:.3f}\t"
            f"{item.consensus_cds_intron_capture:.3f}"
        )

    log_section("Runtime summary")
    print("stage	duration	cumulative	current_rss_mib	peak_rss_mib")
    for stage_name, duration, cumulative, current_rss, peak_rss in timer.rows:
        print(
            f"{stage_name}	{format_duration(duration)}	{format_duration(cumulative)}	"
            f"{current_rss:.1f}	{peak_rss:.1f}"
        )
    total_elapsed = timer.total_elapsed()
    print(
        f"Total	{format_duration(total_elapsed)}	{format_duration(total_elapsed)}	"
        f"{current_rss_mb():.1f}	{peak_rss_mb():.1f}"
    )

    log_section("Output summary")
    print(f"Output profile: {args.output_profile}; stop_after: {args.stop_after}")
    print(f"Candidate coding models: {len(coding_models):,}")
    print(f"Superloci: {len(superloci):,}")
    print(f"Subloci: {len(subloci):,}")
    print(f"Model structural rows: {len(model_structural_rows):,}")
    if stage4_result is not None:
        print(f"Sublocus winners: {len(stage4_result['sublocus_winners']):,}")
        print(f"Monosublocus selected models: {len(stage4_result['selected_rows']):,}")
        print(f"Monosublocus suppressed winners: {len(stage4_result['suppressed_rows']):,}")
    if stage5_result is not None:
        print(f"Backfilled models: {len(stage5_result['backfilled_rows']):,}")
        print(f"Final selected models: {len(stage5_result['final_rows']):,}")
    if final_prune_result is not None:
        print(
            "Final-pruned retained primary models: "
            f"{len(final_prune_result['retained_rows']):,}"
        )
        print(f"Final-pruned removed models: {len(final_prune_result['removed_rows']):,}")
    if public_final_annotation is not None:
        print(f"Final annotation genes: {public_final_annotation.gene_count:,}")
        print(
            "Final annotation transcripts: "
            f"{public_final_annotation.transcript_count:,}"
        )
        print(
            "  primary transcripts: "
            f"{public_final_annotation.primary_transcript_count:,}"
        )
        print(
            "  alternative transcripts: "
            f"{public_final_annotation.alternative_transcript_count:,}"
        )

    selected_pre_pruning_models = list(stage5_result.get("final_models", [])) if stage5_result else []
    selected_pre_pruning_internal_stops = sum(
        bool(model.sequence_integrity and model.sequence_integrity.internal_stop_count > 0)
        for model in selected_pre_pruning_models
    )
    if public_final_annotation is not None:
        final_internal_stop_models = [
            transcript.source_model
            for gene in public_final_annotation.genes
            for transcript in gene.transcripts
            if transcript.source_model.sequence_integrity
            and transcript.source_model.sequence_integrity.internal_stop_count > 0
        ]
        final_internal_stop_denominator = public_final_annotation.transcript_count
    elif final_prune_result is not None:
        final_models = list(final_prune_result.get("retained_models", []))
        final_internal_stop_models = [
            model for model in final_models
            if model.sequence_integrity and model.sequence_integrity.internal_stop_count > 0
        ]
        final_internal_stop_denominator = len(final_models)
    else:
        final_internal_stop_models = [
            model for model in selected_pre_pruning_models
            if model.sequence_integrity and model.sequence_integrity.internal_stop_count > 0
        ]
        final_internal_stop_denominator = len(selected_pre_pruning_models)
    override_retained_internal_stops = (
        count_retained_internal_stop_evidence_overrides(
            list(final_prune_result.get("retained_rows", []))
        )
        if final_prune_result is not None else 0
    )
    if stage5_result is not None:
        print(
            "Internal-stop selected/final visibility: "
            f"selected_pre_pruning={selected_pre_pruning_internal_stops:,}/{len(selected_pre_pruning_models):,}; "
            f"final_retained={len(final_internal_stop_models):,}/{final_internal_stop_denominator:,}; "
            f"retained_via_evidence_override={override_retained_internal_stops:,}; "
            f"override_enabled={str(args.allow_internal_stop_evidence_rescue).lower()}"
        )
    if final_internal_stop_models:
        print(
            "WARNING: Final annotation contains retained internal-stop CDS models "
            f"({len(final_internal_stop_models):,}). Check the genetic code, source annotations/projections, "
            "genome assembly/sequence, and --allow-internal-stop-evidence-rescue settings.",
            file=sys.stderr,
        )

    scored_source_stats = summarise_source_scaled_scores(model_structural_rows)
    selected_rows_for_summary: list[dict] = []
    selected_label = "."
    removed_rows_for_summary: list[dict] = []
    if final_prune_result is not None:
        selected_rows_for_summary = final_prune_result["retained_rows"]
        selected_label = "final_pruned_retained"
        removed_rows_for_summary = final_prune_result["removed_rows"]
    elif stage5_result is not None:
        selected_rows_for_summary = stage5_result["final_rows"]
        selected_label = "final_selected"
    elif stage4_result is not None:
        selected_rows_for_summary = stage4_result["selected_rows"]
        selected_label = "monosublocus_selected"

    selected_source_stats = summarise_score_rows_by_source(selected_rows_for_summary)
    removed_source_counts = Counter(str(row.get("source_name", ".")) for row in removed_rows_for_summary)

    log_section("Source scaled-score summary")
    print(
        "Mean model_selection_score_scaled is calculated across all scored candidate models per source. "
        "The non-singleton mean excludes rows whose canonical sublocus_model_count is 1. "
        "Selected columns refer to the last selected model set produced by this run."
    )
    print(
        "source\tscored_models\tmean_scaled_score_all\tmean_scaled_score_non_singletons\t"
        "n_scored_non_singletons\tn_scored_singletons\tsingleton_fraction\tmean_raw_score_all\t"
        f"{selected_label}_models\t{selected_label}_mean_scaled_score\tfinal_pruned_removed_models"
    )
    all_sources_for_summary = sorted(
        scored_source_stats,
        key=lambda src: (
            -(scored_source_stats[src]["scaled_sum"] / scored_source_stats[src]["count"] if scored_source_stats[src]["count"] else 0.0),
            src,
        ),
    )
    for source in all_sources_for_summary:
        scored = scored_source_stats[source]
        scored_count = int(scored["count"])
        mean_scaled = scored["scaled_sum"] / scored["count"] if scored["count"] else 0.0
        non_singleton_count = int(scored["non_singleton_count"])
        singleton_count = int(scored["singleton_count"])
        mean_non_singletons = (
            scored["non_singleton_scaled_sum"] / scored["non_singleton_count"]
            if scored["non_singleton_count"] else None
        )
        mean_non_singletons_text = (
            f"{mean_non_singletons:.6f}" if mean_non_singletons is not None else "."
        )
        singleton_fraction = singleton_count / scored_count if scored_count else 0.0
        mean_raw = scored["raw_sum"] / scored["count"] if scored["count"] else 0.0
        selected = selected_source_stats.get(source, {"count": 0.0, "scaled_sum": 0.0})
        selected_count = int(selected["count"])
        selected_mean_scaled = selected["scaled_sum"] / selected["count"] if selected["count"] else None
        selected_count_text = str(selected_count) if selected_label != "." else "."
        selected_mean_text = f"{selected_mean_scaled:.6f}" if selected_mean_scaled is not None and selected_label != "." else "."
        removed_text = str(removed_source_counts.get(source, 0)) if final_prune_result is not None else "."
        print(
            f"{source}\t{scored_count}\t{mean_scaled:.6f}\t{mean_non_singletons_text}\t"
            f"{non_singleton_count}\t{singleton_count}\t{singleton_fraction:.6f}\t{mean_raw:.6f}\t"
            f"{selected_count_text}\t{selected_mean_text}\t{removed_text}"
        )

    log_section("Output files")
    if not context.written_output_items:
        print("  No output files were written under the selected --output-profile/--skip-output settings.")
    else:
        for label, path in context.written_output_items:
            print(f"  {label}: {path}")
    if args.output_profile not in {"full", "debug"}:
        print("  Note: use --output-profile full to write all intermediate diagnostic files, or --write-output NAME to add specific files.")

    progress.log("Run complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
