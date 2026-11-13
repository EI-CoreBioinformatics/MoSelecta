"""Constants and default configuration for MoSelecta."""

from __future__ import annotations

PROGRAM_NAME = "MoSelecta"

PROGRAM_VERSION = "0.3.0"

GFF3_SOURCE = "MoSelecta"

TRANSCRIPT_FEATURES = {
    "mRNA",
    "transcript",
    "lnc_RNA",
    "ncRNA",
    "rRNA",
    "tRNA",
    "snRNA",
    "snoRNA",
    "miRNA",
    "pseudogenic_transcript",
}

CHILD_FEATURE_ORDER = {
    "exon": 0,
    "five_prime_UTR": 1,
    "three_prime_UTR": 2,
    "CDS": 3,
}

PROGRESS_LEVELS = {"none": 0, "basic": 1, "verbose": 2, "debug": 3}

OUTPUT_PROFILES = {"minimal", "standard", "full", "debug"}

STOP_AFTER_CHOICES = ("scored-candidates", "monosublocus", "backfill", "final")

PRUNING_LEVEL_CHOICES = (
    "none", "low", "medium", "high", "very-high",
    "auto", "auto-conservative",
)

PRUNING_LEVEL_PRESETS = {
    "none": {
        "singleton_level": "none",
        "shadow_enabled": False,
        "shadow_require_loser_single_exon": False,
        "singleton_rescue_monosublocus_overlap": False,
        "low_score_singletons_when_protein_evidence": False,
        "low_score_nonsingletons": False,
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingletons_when_genome": False,
        "sequence_integrity_nonsingleton_scope": "all",
    },
    "low": {
        "singleton_level": "intrinsic",
        "shadow_enabled": True,
        "shadow_require_loser_single_exon": True,
        "singleton_rescue_monosublocus_overlap": True,
        "low_score_singletons_when_protein_evidence": False,
        "low_score_nonsingletons": False,
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingletons_when_genome": False,
        "sequence_integrity_nonsingleton_scope": "all",
    },
    "medium": {
        "singleton_level": "intrinsic_or_sequence_integrity",
        "shadow_enabled": True,
        "shadow_require_loser_single_exon": False,
        "singleton_rescue_monosublocus_overlap": True,
        "low_score_singletons_when_protein_evidence": True,
        "low_score_nonsingletons": True,
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingletons_when_genome": True,
        "sequence_integrity_nonsingleton_scope": "backfill",
    },
    "high": {
        "singleton_level": "all",
        "shadow_enabled": True,
        "shadow_require_loser_single_exon": False,
        "singleton_rescue_monosublocus_overlap": True,
        "low_score_singletons_when_protein_evidence": True,
        "low_score_nonsingletons": True,
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingletons_when_genome": True,
        "sequence_integrity_nonsingleton_scope": "all",
    },
    "very-high": {
        "singleton_level": "all",
        "shadow_enabled": True,
        "shadow_require_loser_single_exon": False,
        "singleton_rescue_monosublocus_overlap": True,
        "low_score_singletons_when_protein_evidence": True,
        "low_score_nonsingletons": True,
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingletons_when_genome": True,
        "sequence_integrity_nonsingleton_scope": "all",
    },
    # Automatic modes simulate the three fixed scenarios.  Medium is used as
    # the provisional configuration until the post-selection decision is made.
    "auto": {
        "singleton_level": "intrinsic_or_sequence_integrity",
        "shadow_enabled": True,
        "shadow_require_loser_single_exon": False,
        "singleton_rescue_monosublocus_overlap": True,
        "low_score_singletons_when_protein_evidence": True,
        "low_score_nonsingletons": True,
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingletons_when_genome": True,
        "sequence_integrity_nonsingleton_scope": "backfill",
    },
    "auto-conservative": {
        "singleton_level": "intrinsic_or_sequence_integrity",
        "shadow_enabled": True,
        "shadow_require_loser_single_exon": False,
        "singleton_rescue_monosublocus_overlap": True,
        "low_score_singletons_when_protein_evidence": True,
        "low_score_nonsingletons": True,
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingletons_when_genome": True,
        "sequence_integrity_nonsingleton_scope": "backfill",
    },
}

PRUNING_SCENARIO_THRESHOLDS = {
    "medium": {
        "singleton": (0.25, 0.25, 0.40, 0.25),
        "nonsingleton": (0.20, 0.20, 0.40, 0.20),
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingleton_scope": "backfill",
    },
    "high": {
        "singleton": (0.35, 0.35, 0.50, 0.35),
        "nonsingleton": (0.35, 0.35, 0.50, 0.35),
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingleton_scope": "all",
    },
    "very-high": {
        "singleton": (0.40, 0.40, 0.60, 0.40),
        "nonsingleton": (0.40, 0.40, 0.60, 0.40),
        "low_score_nonsingleton_scope": "all",
        "sequence_integrity_nonsingleton_scope": "all",
    },
}

# Canonical structural-confidence thresholds used only when no supplied
# external evidence class passes the run-wide adequacy gate. These deliberately
# remain separate from the configurable evidence-aware thresholds above.
MODEL_CONFIDENCE_NONSINGLETON_THRESHOLDS = {
    "medium": 0.20,
    "high": 0.25,
    "very-high": 0.25,
}
MODEL_CONFIDENCE_LOW_DEPTH_5PLUS_THRESHOLD = 0.35
MODEL_CONFIDENCE_VERY_LOW_DEPTH_7PLUS_THRESHOLD = 0.45

STRUCTURAL_MULTI_SCORE_WEIGHTS = {
    # Balance feature-level agreement with exact model-level support. Exact
    # CDS-exon-chain support now has the same influence as splice-site and
    # internal CDS-exon support, while intron support remains the largest single
    # component. These defaults sum to 1.0 when all multi-exon components exist.
    "cds_intron_support_fraction": 0.25,
    "cds_splice_site_support_fraction": 0.20,
    "internal_cds_exon_support_fraction": 0.20,
    "terminal_cds_exon_support_fraction": 0.15,
    "exact_intron_chain_support_fraction": 0.20,
}

STRUCTURAL_SINGLE_SCORE_WEIGHTS = {
    "single_exon_overlap_support_fraction": 0.60,
    "single_exon_exact_support_fraction": 0.40,
}

CONSENSUS_MULTI_SCORE_WEIGHTS = {
    "cds_intron_consensus_f1": 0.30,
    "cds_splice_site_consensus_f1": 0.25,
    "internal_cds_exon_consensus_f1": 0.20,
    "cds_base_consensus_f1": 0.15,
    "terminal_cds_exon_consensus_f1": 0.10,
}

CONSENSUS_SINGLE_SCORE_WEIGHTS = {
    "cds_base_consensus_f1": 0.75,
    "single_exon_overlap_support_fraction": 0.25,
}

PAIRWISE_MULTI_SCORE_WEIGHTS = {
    "cds_intron_f1": 0.30,
    "cds_splice_site_f1": 0.25,
    "internal_cds_exon_f1": 0.20,
    "cds_base_f1": 0.15,
    "terminal_cds_exon_f1": 0.10,
}

PAIRWISE_SINGLE_SCORE_WEIGHTS = {
    "cds_base_f1": 1.00,
}

COMPLETENESS_PENALTY_WEIGHTS = {
    "missing_supported_terminal_structure_fraction": 0.50,
    "missing_supported_internal_structure_fraction": 0.35,
    "contained_model_source_fraction": 0.15,
}

MODEL_SELECTION_DEFAULT_WEIGHTS = {
    "structural_support": 1.00,
    "consensus_similarity": 0.75,
    "pairwise_similarity": 1.00,
    "completeness": 0.50,
    "source_prior": 0.25,
    "terminal_extension": 0.10,
    "exact_structure_support": 0.75,
    # Accepted v0.2.15 configuration G: run-aware, source-collapsed whole-gene
    # corroboration with saturating alpha 16.
    "structural_corroboration": 0.50,
    "protein_evidence": 3.00,
    "transcript_evidence": 2.00,
    "unsupported_short_cds_exon_penalty": 0.50,
    "structural_coherence_penalty": 0.30,
    "sequence_integrity_penalty": 1.00,
    "intrinsic_penalty": 1.00,
    "partial_structure_penalty": 0.50,
    "internal_extra_penalty": 0.50,
    "terminal_extra_penalty": 0.00,
}

# P1 is the accepted v0.2.15 default. P0 remains available through the same
# terminal-codon penalty controls for historical scientific reproduction.
SEQUENCE_INTEGRITY_P0_MISSING_START_PENALTY = 0.30
SEQUENCE_INTEGRITY_P0_MISSING_STOP_PENALTY = 0.50
SEQUENCE_INTEGRITY_P1_MISSING_START_PENALTY = 0.50
SEQUENCE_INTEGRITY_P1_MISSING_STOP_PENALTY = 0.75
SEQUENCE_INTEGRITY_DEFAULT_MISSING_START_PENALTY = (
    SEQUENCE_INTEGRITY_P1_MISSING_START_PENALTY
)
SEQUENCE_INTEGRITY_DEFAULT_MISSING_STOP_PENALTY = (
    SEQUENCE_INTEGRITY_P1_MISSING_STOP_PENALTY
)

LOW_SOURCE_CONFLICT_AGREEMENT_FACTOR = 0.25

LOW_SOURCE_CONFLICT_PAIRWISE_THRESHOLD = 0.80

INTRINSIC_GFF_FLAG_COUNT_KEYS = [
    ("short_cds", "intrinsic_penalty_short_cds_count"),
    ("cds_length_not_multiple_of_3", "intrinsic_penalty_cds_length_not_multiple_of_3_count"),
    ("short_terminal_cds_exon", "intrinsic_penalty_short_terminal_cds_exon_count"),
    ("short_internal_cds_exon", "intrinsic_penalty_short_internal_cds_exon_count"),
    ("short_cds_intron", "intrinsic_penalty_short_cds_intron_count"),
    ("long_cds_intron", "intrinsic_penalty_long_cds_intron_count"),
    ("contained_model", "intrinsic_penalty_contained_model_count"),
    ("unsupported_cds_intron_contained_in_other_cds_exon", "intrinsic_penalty_unsupported_cds_intron_contained_in_other_cds_exon_count"),
]

FINAL_PRUNE_EXTRA_COLUMNS = [
    "final_prune_level",
    "final_prune_classes",
    "final_singleton_prune_level",
    "final_singleton_prune_source_scope",
    "final_singleton_pruned",
    "final_singleton_prune_reason",
    "final_singleton_prune_sublocus_model_count",
    "final_singleton_prune_sublocus_source_count",
    "final_singleton_rescue_enabled",
    "final_singleton_rescued",
    "final_singleton_rescue_reason",
    "final_singleton_rescue_monosublocus_id",
    "final_singleton_rescue_partner_model_id",
    "final_singleton_rescue_partner_source",
    "final_singleton_rescue_overlap_bp",
    "final_singleton_rescue_self_overlap_fraction",
    "final_singleton_rescue_partner_overlap_fraction",
    "final_singleton_rescue_reciprocal_overlap_fraction",
    "final_singleton_rescue_min_reciprocal_overlap_fraction",
    "final_singleton_rescue_holder_member_count",
    "final_singleton_overlap_rescue_evidence_policy",
    "final_singleton_overlap_rescue_geometry_passed",
    "final_singleton_overlap_rescue_geometry_reason",
    "final_singleton_overlap_rescue_overlap_bp_geometry_candidate",
    "final_singleton_overlap_rescue_reciprocal_rejected",
    "final_singleton_overlap_rescue_same_source_edge",
    "final_singleton_overlap_rescue_any_adequacy_passed",
    "final_singleton_overlap_rescue_protein_adequacy_passed",
    "final_singleton_overlap_rescue_transcript_adequacy_passed",
    "final_singleton_overlap_rescue_strong_local_protein",
    "final_singleton_overlap_rescue_strong_local_transcript",
    "final_singleton_overlap_rescue_decision",
    "final_singleton_overlap_rescue_denial_reason",
    "final_protein_singleton_rescue_enabled",
    "final_protein_singleton_rescued",
    "final_protein_singleton_rescue_reason",
    "final_protein_singleton_rescue_criteria_met",
    "final_protein_singleton_rescue_criteria_reason",
    "final_protein_singleton_rescue_applied",
    "final_protein_singleton_rescue_applied_reason",
    "final_protein_singleton_rescue_min_score_scaled",
    "final_protein_singleton_rescue_min_exact_chain_source_count",
    "final_protein_singleton_rescue_score_scaled",
    "final_protein_singleton_rescue_exact_chain_source_count",
    "final_protein_singleton_rescue_supporting_source_count",
    "final_protein_singleton_rescue_status",
    "final_protein_singleton_rescue_flags",
    "final_transcript_singleton_rescue_enabled",
    "final_transcript_singleton_rescued",
    "final_transcript_singleton_rescue_reason",
    "final_transcript_singleton_rescue_criteria_met",
    "final_transcript_singleton_rescue_criteria_reason",
    "final_transcript_singleton_rescue_applied",
    "final_transcript_singleton_rescue_applied_reason",
    "final_transcript_singleton_rescue_min_score_scaled",
    "final_transcript_singleton_rescue_min_supporting_source_count",
    "final_transcript_singleton_rescue_min_scoring_eligible_count",
    "final_transcript_singleton_rescue_min_cdna_context_fraction",
    "final_transcript_singleton_rescue_require_stranded_support",
    "final_transcript_singleton_rescue_block_on_cdna_context_concern",
    "final_transcript_singleton_rescue_cdna_context_eligible",
    "final_transcript_singleton_rescue_cdna_context_reason",
    "final_transcript_singleton_rescue_cdna_context_fraction_state",
    "final_transcript_singleton_rescue_cdna_context_fraction",
    "final_transcript_singleton_rescue_cdna_context_concern_class",
    "final_transcript_singleton_rescue_score_scaled",
    "final_transcript_singleton_rescue_supporting_source_count",
    "final_transcript_singleton_rescue_scoring_eligible_count",
    "final_transcript_singleton_rescue_stranded_scoring_eligible_count",
    "final_transcript_singleton_rescue_unstranded_scoring_eligible_count",
    "final_transcript_singleton_rescue_support_type",
    "final_transcript_singleton_rescue_status",
    "final_transcript_singleton_rescue_flags",
    "final_low_score_singleton_prune_enabled",
    "final_low_score_singleton_pruned",
    "final_low_score_singleton_prune_reason",
    "final_low_score_singleton_external_evidence_adequacy_passed",
    "final_low_score_singleton_protein_evidence_supplied",
    "final_low_score_singleton_protein_adequacy_required",
    "final_low_score_singleton_protein_adequacy_passed",
    "final_low_score_singleton_protein_overlap_model_fraction",
    "final_low_score_singleton_protein_supported_model_fraction",
    "final_low_score_singleton_min_protein_overlap_model_fraction",
    "final_low_score_singleton_min_protein_supported_model_fraction",
    "final_low_score_singleton_transcript_evidence_supplied",
    "final_low_score_singleton_transcript_adequacy_required",
    "final_low_score_singleton_transcript_adequacy_passed",
    "final_low_score_singleton_transcript_scoring_eligible_model_fraction",
    "final_low_score_singleton_transcript_supported_model_fraction",
    "final_low_score_singleton_min_transcript_scoring_eligible_model_fraction",
    "final_low_score_singleton_min_transcript_supported_model_fraction",
    "final_low_score_singleton_max_model_structural_support_score_scaled",
    "final_low_score_singleton_max_protein_evidence_score_scaled",
    "final_low_score_singleton_max_transcript_evidence_score_scaled",
    "final_low_score_singleton_max_model_selection_score_scaled",
    "final_low_score_singleton_model_structural_support_score_scaled",
    "final_low_score_singleton_protein_evidence_score_scaled",
    "final_low_score_singleton_transcript_evidence_score_scaled",
    "final_low_score_singleton_model_selection_score_scaled",
    "final_low_score_singleton_prune_route",
    "final_low_score_singleton_use_transcript_cdna_context",
    "final_low_score_singleton_cdna_context_pruning_concern",
    "final_low_score_singleton_cdna_context_pruning_concern_class",
    "final_low_score_singleton_cdna_context_pruning_supporting_flags",
    "final_low_score_nonsingleton_prune_enabled",
    "final_low_score_nonsingleton_prune_scope",
    "final_low_score_nonsingleton_pruned",
    "final_low_score_nonsingleton_prune_reason",
    "final_low_score_nonsingleton_policy",
    "final_low_score_nonsingleton_run_model_source_count",
    "final_low_score_nonsingleton_sublocus_source_count",
    "final_low_score_nonsingleton_source_depth_fraction",
    "final_low_score_nonsingleton_effective_model_confidence_threshold",
    "final_low_score_nonsingleton_model_confidence_tier",
    "final_low_score_nonsingleton_external_evidence_adequacy_passed",
    "final_low_score_nonsingleton_protein_evidence_supplied",
    "final_low_score_nonsingleton_protein_adequacy_required",
    "final_low_score_nonsingleton_protein_adequacy_passed",
    "final_low_score_nonsingleton_protein_overlap_model_fraction",
    "final_low_score_nonsingleton_protein_supported_model_fraction",
    "final_low_score_nonsingleton_min_protein_overlap_model_fraction",
    "final_low_score_nonsingleton_min_protein_supported_model_fraction",
    "final_low_score_nonsingleton_transcript_evidence_supplied",
    "final_low_score_nonsingleton_transcript_adequacy_required",
    "final_low_score_nonsingleton_transcript_adequacy_passed",
    "final_low_score_nonsingleton_transcript_scoring_eligible_model_fraction",
    "final_low_score_nonsingleton_transcript_supported_model_fraction",
    "final_low_score_nonsingleton_min_transcript_scoring_eligible_model_fraction",
    "final_low_score_nonsingleton_min_transcript_supported_model_fraction",
    "final_low_score_nonsingleton_max_model_structural_support_score_scaled",
    "final_low_score_nonsingleton_max_protein_evidence_score_scaled",
    "final_low_score_nonsingleton_max_transcript_evidence_score_scaled",
    "final_low_score_nonsingleton_max_model_selection_score_scaled",
    "final_low_score_nonsingleton_model_structural_support_score_scaled",
    "final_low_score_nonsingleton_protein_evidence_score_scaled",
    "final_low_score_nonsingleton_transcript_evidence_score_scaled",
    "final_low_score_nonsingleton_model_selection_score_scaled",
    "final_low_score_nonsingleton_prune_route",
    "final_low_score_nonsingleton_use_transcript_cdna_context",
    "final_low_score_nonsingleton_cdna_context_pruning_concern",
    "final_low_score_nonsingleton_cdna_context_pruning_concern_class",
    "final_low_score_nonsingleton_cdna_context_pruning_supporting_flags",
    "sequence_integrity_pruning_concern",
    "sequence_integrity_pruning_concern_class",
    "sequence_integrity_pruning_reasons",
    "sequence_integrity_severe_pruning_concern",
    "final_model_has_internal_stop_codon",
    "final_allow_internal_stop_evidence_rescue",
    "final_internal_stop_evidence_criteria_met",
    "final_internal_stop_evidence_criteria_routes",
    "final_internal_stop_evidence_rescue_blocked",
    "final_internal_stop_evidence_rescue_block_reason",
    "final_internal_stop_evidence_rescue_applied",
    "final_internal_stop_evidence_rescue_applied_route",
    "final_sequence_integrity_singleton_candidate",
    "final_sequence_integrity_singleton_rescued",
    "final_sequence_integrity_singleton_pruned",
    "final_sequence_integrity_singleton_prune_reason",
    "final_sequence_integrity_nonsingleton_prune_enabled",
    "final_sequence_integrity_nonsingleton_prune_scope",
    "final_sequence_integrity_nonsingleton_candidate",
    "final_sequence_integrity_nonsingleton_pruned",
    "final_sequence_integrity_nonsingleton_prune_route",
    "final_sequence_integrity_nonsingleton_prune_reason",
    "final_sequence_integrity_nonsingleton_protein_protected",
    "final_sequence_integrity_nonsingleton_transcript_protected",
    "final_shadow_prune_enabled",
    "final_shadow_pruned",
    "final_shadow_prune_reason",
    "final_shadow_prune_source_scope",
    "final_shadow_prune_require_loser_single_exon",
    "final_shadow_prune_retained_model_id",
    "final_shadow_prune_retained_model_source",
    "final_shadow_prune_overlap_bp",
    "final_shadow_prune_retained_model_overlap_fraction",
    "final_shadow_prune_removed_model_overlap_fraction",
    "final_shadow_prune_score_margin_scaled",
    "final_shadow_prune_removed_model_count",
    "final_shadow_prune_removed_model_ids",
    "final_monotonic_inherited_removal",
    "final_monotonic_inherited_from_level",
    "final_monotonic_inherited_original_level",
    "final_monotonic_inherited_original_primary_reason",
    "final_monotonic_inherited_original_prune_classes",
    "final_monotonic_native_status",
    "removed_by_medium",
    "removed_by_high",
    "removed_by_very_high",
    "native_removed_by_medium",
    "native_removed_by_high",
    "native_removed_by_very_high",
    "inherited_removal_by_medium",
    "inherited_removal_by_high",
    "inherited_removal_by_very_high",
    "inherited_from_level_medium",
    "inherited_from_level_high",
    "inherited_from_level_very_high",
    "inherited_original_level_medium",
    "inherited_original_level_high",
    "inherited_original_level_very_high",
    "inherited_original_reason_medium",
    "inherited_original_reason_high",
    "inherited_original_reason_very_high",
    "primary_reason_medium",
    "primary_reason_high",
    "primary_reason_very_high",
    "protected_by_protein_in_any_scenario",
    "protected_by_transcript_in_any_scenario",
    "protected_by_overlap_policy_in_any_scenario",
]

SHADOW_PAIR_COLUMNS = [
    "shadow_pair_id", "shadow_pruned", "prune_reason", "chrom",
    "winner_model_id", "winner_source", "winner_strand", "winner_score_scaled", "winner_score_raw",
    "winner_cds_length", "winner_cds_exon_count",
    "loser_model_id", "loser_source", "loser_strand", "loser_score_scaled", "loser_score_raw",
    "loser_cds_length", "loser_cds_exon_count", "overlap_bp",
    "winner_overlap_fraction", "loser_overlap_fraction", "score_margin_scaled",
    "min_overlap_bp", "min_loser_cds_overlap_fraction", "min_score_margin_scaled",
    "require_loser_single_exon", "source_scope",
]

ALL_OUTPUT_KEYS = {
    "combined_gff3", "scored_pruned_gff3", "superloci_tsv", "superlocus_members_tsv",
    "coding_models_tsv", "excluded_models_tsv", "other_features_gff3", "input_summary_tsv",
    "validation_failures_tsv", "subloci_tsv", "sublocus_members_tsv", "relationship_edges_tsv",
    "model_intrinsic_metrics_tsv", "sequence_integrity_summary_tsv", "sequence_integrity_model_metrics_tsv",
    "source_auto_weights_tsv", "model_structural_support_tsv", "structural_corroboration_tsv", "structural_coherence_regions_tsv",
    "structural_features_tsv", "sublocus_winners_gff3", "sublocus_winners_tsv",
    "monosubloci_tsv", "monosublocus_members_tsv", "monosublocus_edges_tsv",
    "monosublocus_review_gff3", "selected_models_gff3", "selected_models_tsv",
    "suppressed_monosublocus_models_tsv", "backfill_candidates_tsv", "backfilled_models_tsv",
    "backfill_rejected_models_tsv", "final_selected_models_tsv", "final_selected_models_gff3",
    "r1_backfill_replacements_tsv", "r1_backfill_summary_tsv",
    "final_pruned_selected_models_tsv", "final_pruned_selected_models_gff3",
    "final_pruned_selected_removed_tsv", "final_shadow_gene_pairs_tsv",
    "standardised_primary_models_gff3", "final_annotation_gff3",
    "run_configuration_json",
    "terminal_stop_canonicalisation_summary_tsv",
    "terminal_stop_canonicalisation_actions_tsv",
    "terminal_stop_canonicalisation_audit_tsv",
    "cds_5prime_boundary_canonicalisation_summary_tsv",
    "cds_5prime_boundary_canonicalisation_actions_tsv",
    "cds_5prime_boundary_canonicalisation_audit_tsv",
    "cds_phase_recovery_summary_tsv", "cds_phase_recovery_actions_tsv",
    "cds_phase_recovery_audit_tsv",
    "definite_internal_stop_exclusion_summary_tsv",
    "definite_internal_stop_exclusions_tsv",
    "extreme_short_cds_exclusion_summary_tsv",
    "extreme_short_cds_exclusions_tsv",
    "cds_terminus_refinement_audit_tsv", "cds_terminus_refinement_summary_tsv",
    "utr_augmented_final_annotation_gff3", "utr_augmentation_audit_tsv",
    "utr_junction_validation_audit_tsv", "utr_junction_validation_summary_tsv",
    "utr_augmentation_summary_tsv",
    "pruning_scenario_summary_tsv", "auto_pruning_diagnostics_tsv",
    "alternative_splice_candidates_tsv", "alternative_splice_genes_tsv",
    "alternative_splice_run_summary_tsv",
    "alternative_splice_candidate_audit_gff3",
    "protein_evidence_summary_tsv", "protein_evidence_raw_to_collapsed_tsv", "protein_evidence_collapsed_alignments_tsv",
    "protein_evidence_superloci_tsv", "protein_evidence_subloci_tsv",
    "protein_evidence_sublocus_members_tsv", "protein_evidence_model_support_tsv",
    "transcript_evidence_summary_tsv", "transcript_evidence_raw_to_collapsed_tsv",
    "transcript_evidence_collapsed_transcripts_tsv", "transcript_evidence_superloci_tsv",
    "transcript_evidence_subloci_tsv", "transcript_evidence_sublocus_members_tsv",
    "transcript_evidence_model_support_tsv",
    "model_selection_scores_tsv", "model_selection_score_summary_tsv",
}
