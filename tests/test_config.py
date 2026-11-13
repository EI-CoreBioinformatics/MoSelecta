from __future__ import annotations

import argparse

import pytest

from moselecta.config import (
    apply_pruning_level_presets,
    build_arg_parser,
    parse_manual_source_weights,
    resolve_enabled_outputs,
)

LOW_SCORE_THRESHOLD_OPTIONS = (
    (
        "--low-score-singleton-max-model-structural-support-score-scaled",
        "low_score_singleton_max_model_structural_support_score_scaled",
    ),
    (
        "--low-score-singleton-max-protein-evidence-score-scaled",
        "low_score_singleton_max_protein_evidence_score_scaled",
    ),
    (
        "--low-score-singleton-max-transcript-evidence-score-scaled",
        "low_score_singleton_max_transcript_evidence_score_scaled",
    ),
    (
        "--low-score-singleton-max-model-selection-score-scaled",
        "low_score_singleton_max_model_selection_score_scaled",
    ),
    (
        "--low-score-nonsingleton-max-model-structural-support-score-scaled",
        "low_score_nonsingleton_max_model_structural_support_score_scaled",
    ),
    (
        "--low-score-nonsingleton-max-protein-evidence-score-scaled",
        "low_score_nonsingleton_max_protein_evidence_score_scaled",
    ),
    (
        "--low-score-nonsingleton-max-transcript-evidence-score-scaled",
        "low_score_nonsingleton_max_transcript_evidence_score_scaled",
    ),
    (
        "--low-score-nonsingleton-max-model-selection-score-scaled",
        "low_score_nonsingleton_max_model_selection_score_scaled",
    ),
)


def test_definite_internal_stop_exclusion_defaults_on_and_has_explicit_opt_out() -> None:
    parser = build_arg_parser()
    required = ["--inputs", "input.gff3", "--out-prefix", "out"]
    assert parser.parse_args(required).definite_internal_stop_exclusion is True
    assert parser.parse_args(
        [*required, "--no-definite-internal-stop-exclusion"]
    ).definite_internal_stop_exclusion is False


def test_extreme_short_cds_exclusion_defaults_on_at_30_coding_bp() -> None:
    parser = build_arg_parser(full_help=True)
    required = ["--inputs", "input.gff3", "--out-prefix", "out"]
    args = parser.parse_args(required)
    assert args.minimum_coding_length_bp == 30
    overridden = parser.parse_args([
        *required, "--minimum-coding-length-bp", "33",
    ])
    assert overridden.minimum_coding_length_bp == 33
    help_text = " ".join(parser.format_help().split())
    assert "--minimum-coding-length-bp" in help_text
    assert "Exactly 30 coding bp is eligible" in help_text
    assert "distinct from --min-cds-length-bp" in help_text


def test_extreme_short_cds_ledger_supports_explicit_output_override() -> None:
    key = "extreme_short_cds_exclusions_tsv"
    assert key not in resolve_enabled_outputs("minimal", None, None)
    assert key in resolve_enabled_outputs("minimal", [key], None)
    assert key not in resolve_enabled_outputs("full", None, [key])


def _preset_namespace(
    level: str,
    rescue=None,
    *,
    evidence_proteins=None,
    evidence_transcripts=None,
    genome_fasta=None,
    low_score_singletons=None,
    low_score_nonsingletons=None,
    low_score_nonsingleton_scope=None,
    sequence_nonsingletons=None,
    sequence_nonsingleton_scope=None,
) -> argparse.Namespace:
    return argparse.Namespace(
        pruning_level=level,
        final_prune_singleton_subloci_level=None,
        final_prune_opposite_strand_shadows=None,
        shadow_prune_require_loser_single_exon=None,
        singleton_prune_rescue_monosublocus_overlap=rescue,
        final_prune_low_score_singletons=low_score_singletons,
        final_prune_low_score_nonsingletons=low_score_nonsingletons,
        low_score_nonsingleton_scope=low_score_nonsingleton_scope,
        final_prune_sequence_integrity_nonsingletons=sequence_nonsingletons,
        sequence_integrity_nonsingleton_scope=sequence_nonsingleton_scope,
        evidence_proteins=evidence_proteins,
        evidence_transcripts=evidence_transcripts,
        genome_fasta=genome_fasta,
    )


def test_pruning_presets_enable_singleton_rescue_for_fixed_and_auto_levels() -> None:
    expected = {
        "low": ("intrinsic", True, True, True),
        "medium": ("intrinsic", True, False, True),
        "high": ("all", True, False, True),
        "very-high": ("all", True, False, True),
        "auto": ("intrinsic", True, False, True),
        "auto-conservative": ("intrinsic", True, False, True),
    }
    for level, values in expected.items():
        args = _preset_namespace(level)
        apply_pruning_level_presets(args)
        assert (
            args.final_prune_singleton_subloci_level,
            args.final_prune_opposite_strand_shadows,
            args.shadow_prune_require_loser_single_exon,
            args.singleton_prune_rescue_monosublocus_overlap,
        ) == values


def test_pruning_none_disables_singleton_rescue() -> None:
    args = _preset_namespace("none")
    apply_pruning_level_presets(args)
    assert args.final_prune_singleton_subloci_level == "none"
    assert args.final_prune_opposite_strand_shadows is False
    assert args.singleton_prune_rescue_monosublocus_overlap is False


def test_pruning_cli_choices_and_default() -> None:
    parser = build_arg_parser()
    action = next(action for action in parser._actions if action.dest == "pruning_level")
    assert action.default == "auto"
    assert tuple(action.choices) == (
        "none", "low", "medium", "high", "very-high", "auto", "auto-conservative",
    )


def test_internal_stop_evidence_rescue_override_defaults_off_and_can_be_enabled() -> None:
    parser = build_arg_parser(expert_help=True)
    required = ["--inputs", "input.gff3", "--out-prefix", "out"]
    assert parser.parse_args(required).allow_internal_stop_evidence_rescue is False
    assert parser.parse_args([
        *required, "--allow-internal-stop-evidence-rescue",
    ]).allow_internal_stop_evidence_rescue is True
    help_text = " ".join(parser.format_help().split())
    assert "--allow-internal-stop-evidence-rescue" in help_text
    assert "including exact-chain protein support" in help_text


def test_r1me_defaults_on_and_explicit_contradiction_is_rejected() -> None:
    from moselecta.pipeline import main

    parser = build_arg_parser(expert_help=True)
    required = ["--inputs", "input.gff3", "--out-prefix", "out"]
    defaults = parser.parse_args(required)
    assert defaults.strict_r1_backfill_replacement is True
    assert defaults.strict_r1_allow_multiexon_blockers is True

    enabled = parser.parse_args([
        *required,
        "--strict-r1-backfill-replacement",
        "--strict-r1-allow-multiexon-blockers",
    ])
    assert enabled.strict_r1_backfill_replacement is True
    assert enabled.strict_r1_allow_multiexon_blockers is True

    with pytest.raises(
        SystemExit,
        match="--strict-r1-allow-multiexon-blockers requires --strict-r1-backfill-replacement",
    ):
        main([
            *required,
            "--no-strict-r1-backfill-replacement",
            "--strict-r1-allow-multiexon-blockers",
        ])


def test_release_default_and_reproducibility_configurations_parse() -> None:
    parser = build_arg_parser(expert_help=True)
    required = ["--inputs", "input.gff3", "--out-prefix", "out"]

    def values(extra: list[str]) -> tuple[float, bool, bool, float, float]:
        args = parser.parse_args([*required, *extra])
        return (
            args.selection_weight_structural_corroboration,
            args.strict_r1_backfill_replacement,
            args.strict_r1_allow_multiexon_blockers,
            args.sequence_penalty_missing_start,
            args.sequence_penalty_missing_stop,
        )

    assert values([]) == (0.5, True, True, 0.50, 0.75)  # G
    assert values([
        "--no-strict-r1-allow-multiexon-blockers",
        "--sequence-penalty-missing-start", "0.30",
        "--sequence-penalty-missing-stop", "0.50",
    ]) == (0.5, True, False, 0.30, 0.50)  # D
    assert values([
        "--sequence-penalty-missing-start", "0.30",
        "--sequence-penalty-missing-stop", "0.50",
    ]) == (0.5, True, True, 0.30, 0.50)  # E
    assert values([
        "--no-strict-r1-allow-multiexon-blockers",
    ]) == (0.5, True, False, 0.50, 0.75)  # F
    assert values([
        "--selection-weight-structural-corroboration", "0",
        "--no-strict-r1-backfill-replacement",
        "--no-strict-r1-allow-multiexon-blockers",
        "--sequence-penalty-missing-start", "0.30",
        "--sequence-penalty-missing-stop", "0.50",
    ]) == (0.0, False, False, 0.30, 0.50)  # frozen-equivalent controls


def test_obsolete_upstream_start_distance_option_is_completely_removed() -> None:
    parser = build_arg_parser(full_help=True)
    help_text = parser.format_help()
    assert "--cds-refinement-upstream-start-search-bp" not in help_text
    with pytest.raises(SystemExit):
        parser.parse_args([
            "--inputs",
            "input.gff3",
            "--out-prefix",
            "out",
            "--cds-refinement-upstream-start-search-bp",
            "100",
        ])


def test_alternative_splice_report_cli_defaults_and_explicit_paths() -> None:
    parser = build_arg_parser(full_help=True)
    defaults = parser.parse_args(["--inputs", "input.gff3", "--out-prefix", "out"])
    assert defaults.alternative_splice_report is False
    assert defaults.no_alternative_splice is False
    assert defaults.alternative_splice_candidates_tsv is None
    assert defaults.alternative_splice_genes_tsv is None
    assert defaults.alternative_splice_run_summary_tsv is None
    assert defaults.alternative_splice_candidate_audit_gff3 is None
    assert defaults.final_annotation_gff3 is None
    assert defaults.alternative_splice_preset == "intermediate"
    assert defaults.alternative_splice_support_policy == "exact"
    assert defaults.alternative_splice_evidence_requirement == "none"
    assert defaults.alternative_splice_max_alternatives_per_gene == 50
    assert defaults.alternative_splice_min_model_support_sources == 2
    assert defaults.alternative_splice_min_model_support_source_fraction is None
    assert defaults.alternative_splice_min_reciprocal_cds_overlap is None
    assert defaults.alternative_splice_min_frame_concordance is None
    assert defaults.alternative_splice_min_model_selection_score_scaled is None
    assert defaults.alternative_splice_min_primary_score_ratio is None

    explicit = parser.parse_args([
        "--inputs", "input.gff3", "--out-prefix", "out",
        "--alternative-splice-report",
        "--no-alternative-splice",
        "--alternative-splice-candidates-tsv", "candidates.tsv",
        "--alternative-splice-genes-tsv", "genes.tsv",
        "--alternative-splice-run-summary-tsv", "run.tsv",
        "--final-annotation-gff3", "final.gff3",
        "--alternative-splice-candidate-audit-gff3", "audit.gff3",
        "--alternative-splice-preset", "stringent",
        "--alternative-splice-support-policy", "exact-plus-mosaic",
        "--alternative-splice-evidence-requirement", "protein-and-transcript",
        "--alternative-splice-max-alternatives-per-gene", "12",
        "--alternative-splice-min-model-support-sources", "3",
        "--alternative-splice-min-model-support-source-fraction", "0.50",
        "--alternative-splice-min-reciprocal-cds-overlap", "0.81",
        "--alternative-splice-min-frame-concordance", "0.92",
        "--alternative-splice-min-model-selection-score-scaled", "0.51",
        "--alternative-splice-min-primary-score-ratio", "0.62",
    ])
    assert explicit.alternative_splice_report is True
    assert explicit.no_alternative_splice is True
    assert str(explicit.alternative_splice_candidates_tsv) == "candidates.tsv"
    assert str(explicit.alternative_splice_genes_tsv) == "genes.tsv"
    assert str(explicit.alternative_splice_run_summary_tsv) == "run.tsv"
    assert str(explicit.final_annotation_gff3) == "final.gff3"
    assert str(explicit.alternative_splice_candidate_audit_gff3) == "audit.gff3"
    assert explicit.alternative_splice_preset == "stringent"
    assert explicit.alternative_splice_support_policy == "exact-plus-mosaic"
    assert explicit.alternative_splice_evidence_requirement == "protein-and-transcript"
    assert explicit.alternative_splice_max_alternatives_per_gene == 12
    assert explicit.alternative_splice_min_model_support_sources == 3
    assert explicit.alternative_splice_min_model_support_source_fraction == 0.50
    assert explicit.alternative_splice_min_reciprocal_cds_overlap == 0.81
    assert explicit.alternative_splice_min_frame_concordance == 0.92
    assert explicit.alternative_splice_min_model_selection_score_scaled == 0.51
    assert explicit.alternative_splice_min_primary_score_ratio == 0.62


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--alternative-splice-min-reciprocal-cds-overlap", "-0.1"),
        ("--alternative-splice-min-frame-concordance", "1.1"),
        ("--alternative-splice-min-model-selection-score-scaled", "nan"),
        ("--alternative-splice-min-primary-score-ratio", "2"),
        ("--alternative-splice-max-alternatives-per-gene", "-1"),
        ("--alternative-splice-min-model-support-sources", "1"),
        ("--alternative-splice-min-model-support-source-fraction", "0"),
        ("--alternative-splice-min-model-support-source-fraction", "1.01"),
        ("--alternative-splice-min-model-support-source-fraction", "nan"),
    ],
)
def test_alternative_splice_cli_rejects_invalid_thresholds_and_cap(
    option: str, value: str
) -> None:
    parser = build_arg_parser(full_help=True)
    with pytest.raises(SystemExit):
        parser.parse_args([
            "--inputs", "input.gff3", "--out-prefix", "out", option, value
        ])


def test_alternative_splice_controls_follow_help_hierarchy() -> None:
    everyday_help = build_arg_parser().format_help()
    full_help = build_arg_parser(full_help=True).format_help()
    expert_help = build_arg_parser(expert_help=True).format_help()
    assert "--alternative-splice-report" not in everyday_help
    assert "--no-alternative-splice" in everyday_help
    assert "--alternative-splice-report" in full_help
    assert "--alternative-splice-evidence-requirement" in full_help
    assert "--alternative-splice-min-model-support-sources" in full_help
    assert "--alternative-splice-min-model-support-source-fraction" in full_help
    assert "--alternative-splice-preview" not in full_help
    for obsolete_term in (
        "stage2", "stage 2", "contract-v1", "pending_calibration", "gate3",
    ):
        assert obsolete_term not in full_help.lower()
    assert "--alternative-splice-candidates-tsv" not in everyday_help
    assert "--final-annotation-gff3" not in everyday_help
    assert "--alternative-splice-candidates-tsv" not in full_help
    assert "--final-annotation-gff3" not in full_help
    assert "--alternative-splice-candidate-audit-gff3" not in full_help
    assert "--alternative-splice-candidates-tsv" in expert_help
    assert "--final-annotation-gff3" in expert_help
    assert "--alternative-splice-candidate-audit-gff3" in expert_help


def test_final_and_audit_output_profile_defaults() -> None:
    for profile in ("minimal", "standard"):
        enabled = resolve_enabled_outputs(profile, None, None)
        assert "final_annotation_gff3" in enabled
        assert "alternative_splice_candidate_audit_gff3" not in enabled
    enabled = resolve_enabled_outputs(
        "minimal",
        ["alternative_splice_candidate_audit_gff3"],
        ["alternative_splice_candidate_audit_gff3"],
    )
    assert "final_annotation_gff3" in enabled
    assert "alternative_splice_candidate_audit_gff3" not in enabled


def test_standardised_primary_output_defaults_path_and_help() -> None:
    parser = build_arg_parser(expert_help=True)
    defaults = parser.parse_args(["--inputs", "input.gff3", "--out-prefix", "out"])
    explicit = parser.parse_args([
        "--inputs", "input.gff3", "--out-prefix", "out",
        "--standardised-primary-models-gff3", "public.gff3",
    ])
    assert defaults.standardised_primary_models_gff3 is None
    assert str(explicit.standardised_primary_models_gff3) == "public.gff3"
    assert "standardised_primary_models_gff3" not in resolve_enabled_outputs(
        "minimal", None, None
    )
    assert "standardised_primary_models_gff3" in resolve_enabled_outputs(
        "standard", None, None
    )
    assert "--standardised-primary-models-gff3" not in build_arg_parser().format_help()
    assert "--standardised-primary-models-gff3" in parser.format_help()


@pytest.mark.parametrize("level", ("none", "low", "medium", "high", "very-high"))
@pytest.mark.parametrize(("option", "destination"), LOW_SCORE_THRESHOLD_OPTIONS)
def test_explicit_low_score_threshold_overrides_fixed_presets(
    level: str,
    option: str,
    destination: str,
) -> None:
    parser = build_arg_parser(full_help=True)
    args = parser.parse_args(
        ["--inputs", "input.gff3", "--out-prefix", "out", "--pruning-level", level, option, "0.731"]
    )
    apply_pruning_level_presets(args)
    assert getattr(args, destination) == pytest.approx(0.731)


@pytest.mark.parametrize("level", ("auto", "auto-conservative"))
@pytest.mark.parametrize(("option", "destination"), LOW_SCORE_THRESHOLD_OPTIONS)
def test_automatic_modes_reject_per_threshold_overrides(
    level: str,
    option: str,
    destination: str,
) -> None:
    del destination
    parser = build_arg_parser(full_help=True)
    args = parser.parse_args(
        ["--inputs", "input.gff3", "--out-prefix", "out", "--pruning-level", level, option, "0.731"]
    )
    with pytest.raises(SystemExit, match="automatic modes simulate canonical") as exc_info:
        apply_pruning_level_presets(args)
    assert option in str(exc_info.value)


@pytest.mark.parametrize(
    ("level", "singleton", "nonsingleton", "low_scope", "sequence_scope"),
    (
        ("medium", (0.25, 0.25, 0.40, 0.25), (0.20, 0.20, 0.40, 0.20), "all", "backfill"),
        ("high", (0.35, 0.35, 0.50, 0.35), (0.35, 0.35, 0.50, 0.35), "all", "all"),
        ("very-high", (0.40, 0.40, 0.60, 0.40), (0.40, 0.40, 0.60, 0.40), "all", "all"),
    ),
)
def test_revised_fixed_preset_thresholds_and_scopes(
    level: str,
    singleton: tuple[float, ...],
    nonsingleton: tuple[float, ...],
    low_scope: str,
    sequence_scope: str,
) -> None:
    args = _preset_namespace(
        level,
        evidence_proteins=["protein.gff3"],
        genome_fasta="genome.fa",
    )
    apply_pruning_level_presets(args)
    assert (
        args.low_score_singleton_max_model_structural_support_score_scaled,
        args.low_score_singleton_max_protein_evidence_score_scaled,
        args.low_score_singleton_max_transcript_evidence_score_scaled,
        args.low_score_singleton_max_model_selection_score_scaled,
    ) == singleton
    assert (
        args.low_score_nonsingleton_max_model_structural_support_score_scaled,
        args.low_score_nonsingleton_max_protein_evidence_score_scaled,
        args.low_score_nonsingleton_max_transcript_evidence_score_scaled,
        args.low_score_nonsingleton_max_model_selection_score_scaled,
    ) == nonsingleton
    assert args.low_score_nonsingleton_scope == low_scope
    assert args.sequence_integrity_nonsingleton_scope == sequence_scope


def test_pruning_explicit_rescue_override_is_preserved() -> None:
    args = _preset_namespace("high", rescue=False)
    apply_pruning_level_presets(args)
    assert args.singleton_prune_rescue_monosublocus_overlap is False


def test_low_score_singleton_pruning_defaults_on_for_medium_high_very_high_with_external_evidence() -> None:
    for level in ("medium", "high", "very-high", "auto", "auto-conservative"):
        args = _preset_namespace(level, evidence_proteins=["protein.gff3"])
        apply_pruning_level_presets(args)
        assert args.final_prune_low_score_singletons is True

        args = _preset_namespace(level, evidence_transcripts=["transcripts.gtf"])
        apply_pruning_level_presets(args)
        assert args.final_prune_low_score_singletons is True


def test_low_score_singleton_pruning_defaults_off_for_low_or_no_external_evidence() -> None:
    args = _preset_namespace("low", evidence_proteins=["protein.gff3"])
    apply_pruning_level_presets(args)
    assert args.final_prune_low_score_singletons is False

    args = _preset_namespace("medium", evidence_proteins=None, evidence_transcripts=None)
    apply_pruning_level_presets(args)
    assert args.final_prune_low_score_singletons is False


def test_low_score_singleton_explicit_override_is_preserved() -> None:
    args = _preset_namespace("high", evidence_proteins=["protein.gff3"], low_score_singletons=False)
    apply_pruning_level_presets(args)
    assert args.final_prune_low_score_singletons is False

    args = _preset_namespace("low", evidence_proteins=None, low_score_singletons=True)
    apply_pruning_level_presets(args)
    assert args.final_prune_low_score_singletons is True



def test_low_score_nonsingleton_pruning_defaults_on_for_revised_levels_with_external_evidence_all_scope() -> None:
    for level in ("medium", "high", "very-high", "auto", "auto-conservative"):
        args = _preset_namespace(level, evidence_proteins=["protein.gff3"])
        apply_pruning_level_presets(args)
        assert args.final_prune_low_score_nonsingletons is True
        assert args.low_score_nonsingleton_scope == "all"

        args = _preset_namespace(level, evidence_transcripts=["transcripts.gtf"])
        apply_pruning_level_presets(args)
        assert args.final_prune_low_score_nonsingletons is True
        assert args.low_score_nonsingleton_scope == "all"


def test_very_high_low_score_nonsingleton_pruning_uses_all_scope() -> None:
    for evidence_kwargs in (
        {"evidence_proteins": ["protein.gff3"]},
        {"evidence_transcripts": ["transcripts.gtf"]},
    ):
        args = _preset_namespace("very-high", **evidence_kwargs)
        apply_pruning_level_presets(args)
        assert args.final_prune_low_score_nonsingletons is True
        assert args.low_score_nonsingleton_scope == "all"


def test_low_score_nonsingleton_pruning_defaults_off_for_low_and_on_without_external_evidence() -> None:
    args = _preset_namespace("low", evidence_proteins=["protein.gff3"])
    apply_pruning_level_presets(args)
    assert args.final_prune_low_score_nonsingletons is False
    assert args.low_score_nonsingleton_scope == "all"


def test_revised_levels_enable_nonsingleton_pruning_without_external_evidence() -> None:
    for level in ("medium", "high", "very-high", "auto", "auto-conservative"):
        args = _preset_namespace(level, evidence_proteins=None, evidence_transcripts=None)
        apply_pruning_level_presets(args)
        assert args.final_prune_low_score_nonsingletons is True
        assert args.low_score_nonsingleton_scope == "all"


def test_explicit_no_nonsingleton_pruning_is_preserved_without_external_evidence() -> None:
    args = _preset_namespace(
        "high",
        evidence_proteins=None,
        evidence_transcripts=None,
        low_score_nonsingletons=False,
    )
    apply_pruning_level_presets(args)
    assert args.final_prune_low_score_nonsingletons is False


def test_low_score_nonsingleton_explicit_override_defaults_to_all_scope_and_preserves_scope() -> None:
    args = _preset_namespace("low", evidence_proteins=None, low_score_nonsingletons=True)
    apply_pruning_level_presets(args)
    assert args.final_prune_low_score_nonsingletons is True
    assert args.low_score_nonsingleton_scope == "all"

    args = _preset_namespace("medium", evidence_proteins=["protein.gff3"], low_score_nonsingletons=True, low_score_nonsingleton_scope="backfill")
    apply_pruning_level_presets(args)
    assert args.final_prune_low_score_nonsingletons is True
    assert args.low_score_nonsingleton_scope == "backfill"

    args = _preset_namespace("high", evidence_proteins=["protein.gff3"], low_score_nonsingletons=False)
    apply_pruning_level_presets(args)
    assert args.final_prune_low_score_nonsingletons is False
    assert args.low_score_nonsingleton_scope == "all"

def test_medium_singleton_preset_adds_sequence_integrity_only_when_genome_supplied() -> None:
    args = _preset_namespace("medium", genome_fasta="genome.fa")
    apply_pruning_level_presets(args)
    assert args.final_prune_singleton_subloci_level == "intrinsic_or_sequence_integrity"

    args = _preset_namespace("medium", genome_fasta=None)
    apply_pruning_level_presets(args)
    assert args.final_prune_singleton_subloci_level == "intrinsic"


def test_sequence_integrity_nonsingleton_pruning_presets_require_genome_and_use_revised_scope() -> None:
    for level, expected_scope in (("medium", "backfill"), ("high", "all"), ("very-high", "all")):
        args = _preset_namespace(level, genome_fasta="genome.fa")
        apply_pruning_level_presets(args)
        assert args.final_prune_sequence_integrity_nonsingletons is True
        assert args.sequence_integrity_nonsingleton_scope == expected_scope

    args = _preset_namespace("medium", genome_fasta=None)
    apply_pruning_level_presets(args)
    assert args.final_prune_sequence_integrity_nonsingletons is False
    assert args.sequence_integrity_nonsingleton_scope == "backfill"


def test_very_high_sequence_integrity_nonsingleton_pruning_uses_all_scope() -> None:
    args = _preset_namespace("very-high", genome_fasta="genome.fa")
    apply_pruning_level_presets(args)
    assert args.final_prune_sequence_integrity_nonsingletons is True
    assert args.sequence_integrity_nonsingleton_scope == "all"


def test_sequence_integrity_diagnostic_only_mode_disables_preset_pruning_but_allows_explicit_overrides() -> None:
    args = _preset_namespace(
        "medium",
        genome_fasta="genome.fa",
    )
    args.no_sequence_integrity_scoring = True
    apply_pruning_level_presets(args)
    assert args.final_prune_singleton_subloci_level == "intrinsic"
    assert args.final_prune_sequence_integrity_nonsingletons is False

    args = _preset_namespace(
        "medium",
        genome_fasta="genome.fa",
        sequence_nonsingletons=True,
    )
    args.no_sequence_integrity_scoring = True
    args.final_prune_singleton_subloci_level = "sequence_integrity"
    apply_pruning_level_presets(args)
    assert args.final_prune_singleton_subloci_level == "sequence_integrity"
    assert args.final_prune_sequence_integrity_nonsingletons is True
    assert args.sequence_integrity_nonsingleton_scope == "all"


def test_sequence_integrity_nonsingleton_explicit_enable_defaults_all_and_preserves_scope() -> None:
    args = _preset_namespace("low", genome_fasta="genome.fa", sequence_nonsingletons=True)
    apply_pruning_level_presets(args)
    assert args.final_prune_sequence_integrity_nonsingletons is True
    assert args.sequence_integrity_nonsingleton_scope == "all"

    args = _preset_namespace(
        "medium",
        genome_fasta="genome.fa",
        sequence_nonsingletons=True,
        sequence_nonsingleton_scope="backfill",
    )
    apply_pruning_level_presets(args)
    assert args.sequence_integrity_nonsingleton_scope == "backfill"


def test_resolve_enabled_outputs_rejects_unknown_name() -> None:
    with pytest.raises(SystemExit):
        resolve_enabled_outputs("minimal", ["not_an_output"], None)


def test_parse_manual_source_weights_supports_comma_and_repeat_specs() -> None:
    weights = parse_manual_source_weights(["A=2.5,B=0", "C=1.25"], ["A", "B", "C"])
    assert weights == {"A": 2.5, "B": 0.0, "C": 1.25}


def test_transcript_evidence_outputs_are_valid_write_output_names() -> None:
    enabled = resolve_enabled_outputs(
        "standard",
        [
            "transcript_evidence_collapsed_transcripts_tsv",
            "transcript_evidence_superloci_tsv",
            "transcript_evidence_subloci_tsv",
            "transcript_evidence_sublocus_members_tsv",
            "transcript_evidence_raw_to_collapsed_tsv",
        ],
        None,
    )
    assert "transcript_evidence_collapsed_transcripts_tsv" in enabled
    assert "transcript_evidence_raw_to_collapsed_tsv" in enabled


def test_patch15_selection_and_context_defaults() -> None:
    from moselecta.config import build_arg_parser

    parser = build_arg_parser(full_help=True)
    args = parser.parse_args(["--inputs", "input.gff3", "--out-prefix", "out"])
    apply_pruning_level_presets(args)
    assert args.selection_weight_transcript_evidence == 2.0
    assert args.selection_weight_protein_evidence == 3.0
    assert args.selection_weight_structural_corroboration == 0.5
    assert args.structural_corroboration_alpha == 16.0
    assert args.selection_penalty_unsupported_short_cds_exon == 0.50
    assert args.unsupported_short_5prime_terminal_cds_exon_bp == 15
    assert args.unsupported_short_3prime_terminal_cds_exon_bp == 15
    assert args.unsupported_short_internal_cds_exon_bp == 15
    assert args.unsupported_short_5prime_terminal_weight == 0.25
    assert args.unsupported_short_3prime_terminal_weight == 0.25
    assert args.unsupported_short_internal_weight == 0.50
    assert args.low_score_singleton_max_transcript_evidence_score_scaled == 0.40
    assert args.low_score_nonsingleton_max_transcript_evidence_score_scaled == 0.40
    assert args.selection_penalty_transcript_cdna_context == 0.15
    assert args.transcript_cdna_context_neutral_fraction == 0.30
    assert args.transcript_cdna_context_severe_fraction == 0.15
    assert args.transcript_singleton_rescue_block_on_cdna_context_concern is True
    assert args.low_score_singleton_use_transcript_cdna_context is True
    assert args.low_score_nonsingleton_use_transcript_cdna_context is True
    assert args.genome_fasta is None
    assert args.genetic_code == 1
    assert args.selection_penalty_sequence_integrity == 1.0
    assert args.no_sequence_integrity_scoring is False
    assert args.sequence_penalty_internal_stop == 1.0
    assert args.sequence_internal_stop_count_cap == 2
    assert args.sequence_penalty_phase_inconsistent == 0.75
    assert args.sequence_penalty_missing_stop == 0.75
    assert args.sequence_penalty_missing_start == 0.50
    assert args.sequence_penalty_noncanonical_splice == 1.00


def test_v0211_mod1_selection_defaults_and_explicit_equivalence() -> None:
    from moselecta.config import build_arg_parser

    parser = build_arg_parser(full_help=True)
    required = ["--inputs", "input.gff3", "--out-prefix", "out"]
    default_args = parser.parse_args(required)
    explicit_args = parser.parse_args([
        *required,
        "--selection-weight-structural-support", "1.0",
        "--selection-weight-pairwise-similarity", "1.0",
        "--selection-weight-exact-structure-support", "0.75",
        "--selection-weight-consensus-similarity", "0.75",
        "--selection-weight-completeness", "0.5",
        "--selection-penalty-partial-structure", "0.5",
        "--selection-penalty-structural-coherence", "0.30",
    ])
    fields = (
        "selection_weight_structural_support",
        "selection_weight_pairwise_similarity",
        "selection_weight_exact_structure_support",
        "selection_weight_consensus_similarity",
        "selection_weight_completeness",
        "selection_penalty_partial_structure",
        "selection_penalty_structural_coherence",
    )
    expected = (1.0, 1.0, 0.75, 0.75, 0.5, 0.5, 0.30)
    assert tuple(getattr(default_args, field) for field in fields) == expected
    assert tuple(getattr(default_args, field) for field in fields) == tuple(
        getattr(explicit_args, field) for field in fields
    )


def test_v0211_explicit_selection_overrides_take_precedence() -> None:
    from moselecta.config import build_arg_parser

    args = build_arg_parser(full_help=True).parse_args([
        "--inputs", "input.gff3", "--out-prefix", "out",
        "--selection-weight-structural-support", "1.25",
        "--selection-weight-pairwise-similarity", "0.60",
        "--selection-weight-exact-structure-support", "0.55",
        "--selection-weight-consensus-similarity", "0.45",
        "--selection-weight-completeness", "0.35",
        "--selection-penalty-partial-structure", "0.25",
        "--selection-penalty-structural-coherence", "0.15",
        "--selection-weight-structural-corroboration", "0.5",
        "--structural-corroboration-alpha", "8",
    ])
    assert args.selection_weight_structural_support == 1.25
    assert args.selection_weight_pairwise_similarity == 0.60
    assert args.selection_weight_exact_structure_support == 0.55
    assert args.selection_weight_consensus_similarity == 0.45
    assert args.selection_weight_completeness == 0.35
    assert args.selection_penalty_partial_structure == 0.25
    assert args.selection_penalty_structural_coherence == 0.15
    assert args.selection_weight_structural_corroboration == 0.5
    assert args.structural_corroboration_alpha == 8.0


def test_output_profile_redesign_exact_membership() -> None:
    from moselecta.config import output_keys_for_profile
    from moselecta.constants import ALL_OUTPUT_KEYS

    minimal = output_keys_for_profile("minimal")
    standard = output_keys_for_profile("standard")
    full = output_keys_for_profile("full")
    debug = output_keys_for_profile("debug")

    assert minimal == {
        "alternative_splice_run_summary_tsv",
        "auto_pruning_diagnostics_tsv",
        "cds_5prime_boundary_canonicalisation_summary_tsv",
        "cds_phase_recovery_summary_tsv",
        "definite_internal_stop_exclusion_summary_tsv",
        "extreme_short_cds_exclusion_summary_tsv",
        "cds_terminus_refinement_summary_tsv",
        "excluded_models_tsv",
        "final_annotation_gff3",
        "input_summary_tsv",
        "model_selection_score_summary_tsv",
        "protein_evidence_summary_tsv",
        "pruning_scenario_summary_tsv",
        "run_configuration_json",
        "sequence_integrity_summary_tsv",
        "source_auto_weights_tsv",
        "terminal_stop_canonicalisation_summary_tsv",
        "transcript_evidence_summary_tsv",
        "utr_augmentation_summary_tsv",
        "utr_augmented_final_annotation_gff3",
        "utr_junction_validation_summary_tsv",
        "validation_failures_tsv",
    }
    assert standard == minimal | {
        "backfilled_models_tsv",
        "cds_5prime_boundary_canonicalisation_actions_tsv",
        "cds_phase_recovery_actions_tsv",
        "definite_internal_stop_exclusions_tsv",
        "extreme_short_cds_exclusions_tsv",
        "cds_terminus_refinement_audit_tsv",
        "final_pruned_selected_models_tsv",
        "final_pruned_selected_removed_tsv",
        "final_shadow_gene_pairs_tsv",
        "selected_models_tsv",
        "standardised_primary_models_gff3",
        "suppressed_monosublocus_models_tsv",
        "terminal_stop_canonicalisation_actions_tsv",
        "utr_augmentation_audit_tsv",
    }
    assert len(minimal) == 22
    assert len(standard) == 36
    assert minimal < standard < full
    assert full == debug == set(ALL_OUTPUT_KEYS) - {
        "r1_backfill_replacements_tsv", "r1_backfill_summary_tsv"
    }
    assert len(full) == 83
    assert "structural_corroboration_tsv" in full
    assert "structural_corroboration_tsv" not in standard

    for summary_key in {
        "terminal_stop_canonicalisation_summary_tsv",
        "cds_5prime_boundary_canonicalisation_summary_tsv",
        "cds_phase_recovery_summary_tsv",
        "definite_internal_stop_exclusion_summary_tsv",
        "extreme_short_cds_exclusion_summary_tsv",
    }:
        assert summary_key in minimal
    for action_key in {
        "terminal_stop_canonicalisation_actions_tsv",
        "cds_5prime_boundary_canonicalisation_actions_tsv",
        "cds_phase_recovery_actions_tsv",
        "definite_internal_stop_exclusions_tsv",
        "extreme_short_cds_exclusions_tsv",
    }:
        assert action_key not in minimal
        assert action_key in standard
    for audit_key in {
        "terminal_stop_canonicalisation_audit_tsv",
        "cds_5prime_boundary_canonicalisation_audit_tsv",
        "cds_phase_recovery_audit_tsv",
    }:
        assert audit_key not in minimal | standard
        assert audit_key in full



def test_patch17_transcript_output_keys_support_write_and_skip_overrides() -> None:
    construction_keys = {
        "transcript_evidence_collapsed_transcripts_tsv",
        "transcript_evidence_superloci_tsv",
        "transcript_evidence_subloci_tsv",
        "transcript_evidence_sublocus_members_tsv",
        "transcript_evidence_raw_to_collapsed_tsv",
    }
    for key in construction_keys:
        assert key in resolve_enabled_outputs("standard", [key], None)
        assert key not in resolve_enabled_outputs("full", None, [key])

    assert "transcript_evidence_summary_tsv" not in resolve_enabled_outputs(
        "minimal", None, ["transcript_evidence_summary_tsv"]
    )
    assert "transcript_evidence_model_support_tsv" not in resolve_enabled_outputs(
        "standard", None, ["transcript_evidence_model_support_tsv"]
    )

def test_three_level_help_separates_routine_advanced_and_expert_controls() -> None:
    from moselecta.config import build_arg_parser

    main_help = build_arg_parser(full_help=False).format_help()
    assert "--full-help" in main_help
    assert "--expert-help" in main_help
    assert "--genome-fasta" in main_help
    assert "--pruning-level" in main_help
    assert "--output-profile" in main_help
    assert "--refine-cds-termini" not in main_help
    assert "--alternative-splice-evidence-requirement" in main_help
    assert "--final-gff3-cds-only-structure" in main_help
    assert "--selection-weight-protein-evidence" not in main_help
    assert "--final-prune-low-score-nonsingletons" not in main_help

    # Scientifically meaningful overrides are in --full-help.
    full_help = build_arg_parser(full_help=True).format_help()
    assert "--selection-weight-protein-evidence" in full_help
    assert "--selection-weight-transcript-evidence" in full_help
    assert "--sequence-penalty-missing-start" in full_help
    assert "--sequence-penalty-missing-stop" in full_help
    assert "--refine-cds-termini" in full_help
    assert "--alternative-splice-evidence-requirement" in full_help
    assert "--final-gff3-cds-only-structure" in full_help
    assert "--strict-r1-backfill-replacement" not in full_help

    # Low-level calibration is reserved for complete expert help.
    expert_help = build_arg_parser(expert_help=True).format_help()
    assert "--low-score-nonsingleton-scope" not in main_help
    assert "--sequence-integrity-nonsingleton-scope" not in main_help
    assert "--selection-penalty-unsupported-short-cds-exon" not in main_help
    assert "--unsupported-short-internal-cds-exon-bp" not in main_help
    assert "--selection-penalty-sequence-integrity" not in main_help
    assert "--sequence-penalty-internal-stop" not in main_help

    assert "--low-score-nonsingleton-scope" not in full_help
    assert "--selection-penalty-unsupported-short-cds-exon" not in full_help
    assert "--low-score-nonsingleton-scope" in expert_help
    assert "--sequence-integrity-nonsingleton-scope" in expert_help
    assert "--selection-penalty-unsupported-short-cds-exon" in expert_help
    assert "--unsupported-short-internal-cds-exon-bp" in expert_help
    assert "--selection-penalty-sequence-integrity" in expert_help
    assert "--sequence-penalty-internal-stop" in expert_help
    assert "--strict-r1-backfill-replacement" in expert_help
    assert "--refine-cds-termini" in expert_help
    assert "--alternative-splice-evidence-requirement" in expert_help
    assert "--final-gff3-cds-only-structure" in expert_help


def test_nonsingleton_fallback_help_distinguishes_context_and_threshold_roles() -> None:
    from moselecta.config import build_arg_parser

    expert_help = " ".join(build_arg_parser(expert_help=True).format_help().split())
    assert "In evidence-aware mode, transcript adequacy must pass." in expert_help
    assert "without requiring global transcript adequacy" in expert_help
    assert "does not change the canonical model-confidence fallback thresholds" in expert_help
    assert "local protein veto threshold in model-confidence fallback" in expert_help
    assert "local transcript veto threshold for the conventional model-confidence fallback route" in expert_help
    assert "confidence fallback does not use model selection score" in expert_help
