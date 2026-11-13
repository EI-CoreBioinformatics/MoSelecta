"""Command-line argument parsing and runtime configuration helpers."""

from __future__ import annotations

import argparse
import math
from collections import Counter
from pathlib import Path
from typing import Optional

from .constants import (
    ALL_OUTPUT_KEYS,
    OUTPUT_PROFILES,
    PROGRAM_NAME,
    PROGRAM_VERSION,
    MODEL_SELECTION_DEFAULT_WEIGHTS,
    PRUNING_LEVEL_CHOICES,
    PRUNING_LEVEL_PRESETS,
    PRUNING_SCENARIO_THRESHOLDS,
    SEQUENCE_INTEGRITY_DEFAULT_MISSING_START_PENALTY,
    SEQUENCE_INTEGRITY_DEFAULT_MISSING_STOP_PENALTY,
    SEQUENCE_INTEGRITY_P0_MISSING_START_PENALTY,
    SEQUENCE_INTEGRITY_P0_MISSING_STOP_PENALTY,
    STOP_AFTER_CHOICES,
)
from .models import Sublocus
from .scoring import value_as_float


def value_as_positive_int(value: object) -> int:
    """Return a non-negative integer count from TSV-style numeric values."""
    count = int(round(value_as_float(value, 0.0)))
    return count if count > 0 else 0


def unit_interval_float(value: str) -> float:
    """Argparse type for inclusive biological fraction thresholds."""
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number between 0 and 1") from exc
    if not 0.0 <= parsed <= 1.0:
        raise argparse.ArgumentTypeError("must be between 0 and 1 inclusive")
    return parsed


def nonnegative_int(value: str) -> int:
    """Argparse type for counts that may intentionally be zero."""
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer >= 0") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return parsed


def nonnegative_float(value: str) -> float:
    """Argparse type for non-negative real-valued thresholds."""
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number >= 0") from exc
    if not math.isfinite(parsed) or parsed < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return parsed


def minimum_two_int(value: str) -> int:
    """Argparse type for prediction-model source-consensus counts."""
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer >= 2") from exc
    if parsed < 2:
        raise argparse.ArgumentTypeError("must be >= 2")
    return parsed


def positive_unit_interval_float(value: str) -> float:
    """Argparse type for fractions greater than zero and at most one."""
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number > 0 and <= 1") from exc
    if not 0.0 < parsed <= 1.0:
        raise argparse.ArgumentTypeError("must be > 0 and <= 1")
    return parsed


class MoSelectaHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Keep defaults in option help while preserving full-help epilog tables."""


BASIC_HELP_DESTINATIONS = frozenset(
    {
        "help", "full_help", "expert_help", "version",
        "inputs", "source_name", "genome_fasta", "genetic_code",
        "evidence_proteins", "evidence_protein_source_name",
        "evidence_transcripts", "evidence_transcript_source_name",
        "progress_log_level", "add_utrs",
        "no_alternative_splice", "alternative_splice_evidence_requirement",
        "out_prefix", "output_profile", "final_gff3_cds_only_structure",
        "pruning_level",
    }
)

EXPERT_ONLY_GROUPS = frozenset(
    {
        "scored-candidate GFF3 pruning filters",
        "first-pass sublocus clustering",
        "structural support metrics",
        "intrinsic penalty weights",
    }
)

FULL_HELP_GROUP_ALLOWLISTS = {
    "outputs": frozenset(
        {
            "out_prefix", "output_profile", "write_output", "skip_output",
            "final_gff3_cds_only_structure",
        }
    ),
    "final selected-model pruning": frozenset(
        {"pruning_level", "auto_pruning_min_level", "auto_pruning_max_level"}
    ),
    "model selection ranking": frozenset(
        {
            "selection_weight_protein_evidence",
            "selection_weight_transcript_evidence",
            "no_sequence_integrity_scoring",
            "sequence_penalty_missing_start",
            "sequence_penalty_missing_stop",
        }
    ),
    "coordinate-only intrinsic model metrics": frozenset(
        {
            "min_cds_length_bp", "min_terminal_cds_exon_bp",
            "min_internal_cds_exon_bp", "min_cds_exon_bp",
            "min_cds_intron_bp", "long_cds_intron_bp",
            "long_cds_intron_percentile",
        }
    ),
}

ALWAYS_EXPERT_DESTINATIONS = frozenset(
    {
        "monosublocus_min_cds_exon_overlap_fraction",
        "strict_r1_backfill_replacement",
        "strict_r1_allow_multiexon_blockers",
    }
)


def expert_only_destinations(parser: argparse.ArgumentParser) -> set[str]:
    """Return actions reserved for complete expert help."""
    expert = set(ALWAYS_EXPERT_DESTINATIONS)
    for group in parser._action_groups:
        for action in group._group_actions:
            if group.title in EXPERT_ONLY_GROUPS:
                expert.add(action.dest)
            elif group.title in FULL_HELP_GROUP_ALLOWLISTS:
                if action.dest not in FULL_HELP_GROUP_ALLOWLISTS[group.title]:
                    expert.add(action.dest)
    return expert


def full_help_epilog() -> str:
    """Return additional detailed help text shown only by --full-help."""
    return """
Input manifests
---------------
--inputs, --evidence-transcripts and --evidence-proteins accept @FILE entries.
Each manifest row is PATH or PATH<TAB>LABEL; blank lines and lines beginning
with # are ignored. Relative paths are resolved from the manifest directory.
Direct paths and manifests may be mixed and retain their supplied order.
No shell or home-directory expansion is performed inside manifest rows.
Nested manifests and duplicate resolved paths are rejected.

An explicit source-label option applies positionally to the complete expanded
input family, must contain exactly one label per expanded file, and overrides
all labels in that family's manifests. Without an explicit label list,
manifest labels override the normal automatic fallback for their own rows.

Output profile contents and named output keys
-------------------------------------------
Use --write-output NAME to add a file to the selected profile, or --skip-output NAME to suppress it.
Comma-separated values and repeated arguments are accepted. Use all with --write-output or --skip-output to enable or suppress all known outputs.

Profiles:
  minimal
    Final annotation products, deterministic run configuration, and essential compact summaries.

  standard
    The normal production set. Adds the standardised primary annotation,
    selected/suppressed/backfilled/final-pruning ledgers, action-only early
    biological reports, and final-model refinement/UTR audits.

  full/debug
    Equivalent in 0.2.13. Write exhaustive candidate, evidence, scoring,
    selection and audit diagnostics, including large intermediate tables.

Explicit --skip-output values are applied last and always win. An explicit
--stop-after scored-candidates run writes the combined scored-candidate GFF3,
model_selection_scores.tsv, model_structural_support.tsv, and the score summary
unless individually skipped.

Genome sequence-integrity outputs:
  sequence_integrity_summary_tsv
    Compact run-level start/stop, internal-stop, phase, splice-motif, boundary, ambiguity, and penalty summary. Included in minimal/standard/full/debug when --genome-fasta is supplied.

  sequence_integrity_model_metrics_tsv
    Model-level genome-derived sequence-integrity metrics and penalty components. Included in full/debug when --genome-fasta is supplied, or when explicitly requested.

Protein/model-selection outputs:
  protein_evidence_summary_tsv
    Compact run-level protein evidence coverage and source summary. Included in minimal/standard/full/debug when protein evidence is supplied.

  protein_evidence_model_support_tsv
    Model-level protein evidence support metrics and status/flag columns. Included in full/debug when protein evidence is supplied, or when explicitly requested.

  model_selection_score_summary_tsv
    Compact model-only versus integrated score summary. Included in minimal/standard/full/debug.

  model_selection_scores_tsv
    Model-level selection-score audit table. Included in full/debug, when stopping after scored candidates, or when explicitly requested.

  structural_coherence_regions_tsv
    Region-level source-collapsed structural-coherence evidence, gate decisions, and penalty values. Included in full/debug or when explicitly requested.

Protein evidence construction diagnostics, included in full/debug or by explicit --write-output:
  protein_evidence_raw_to_collapsed_tsv
  protein_evidence_collapsed_alignments_tsv
  protein_evidence_superloci_tsv
  protein_evidence_subloci_tsv
  protein_evidence_sublocus_members_tsv

Transcript evidence diagnostics:
  transcript_evidence_summary_tsv
    Compact run-level transcript parsing, collapse, locus, model-association, support, and cDNA-context summary. Included in minimal/standard/full/debug when transcript evidence is supplied.

  transcript_evidence_model_support_tsv
    Model-level broad/tight transcript association classes, transcript scores, cDNA-context diagnostics, and support metrics. Included in full/debug when transcript evidence is supplied, or when explicitly requested. Scoring-eligible transcript evidence can contribute to ranking, singleton rescue, and evidence-aware low-score pruning according to the active options.

Transcript evidence construction diagnostics, included in full/debug or by explicit --write-output:
  transcript_evidence_collapsed_transcripts_tsv
  transcript_evidence_superloci_tsv
  transcript_evidence_subloci_tsv
  transcript_evidence_sublocus_members_tsv

  transcript_evidence_raw_to_collapsed_tsv
    Potentially large mapping from raw transcript IDs to collapsed transcript evidence records. Included in full/debug or by explicit --write-output.

Optional alternative-splice diagnostics, written when --alternative-splice-report is enabled:
  alternative_splice_candidates_tsv
    One row per candidate/retained-primary assessment, including historical facts, active preset decisions, ratio availability, greedy conflicts, and canonical preset comparisons.

  alternative_splice_genes_tsv
    One row per frozen retained primary with candidate, threshold, conflict, cap, and accepted internal-model counts.

  alternative_splice_run_summary_tsv
    One run-level row with effective settings, support-route composition, accepted counts, unavailable-ratio counts, and preset nonmonotonicity diagnostics.

Alternative-splice browser diagnostics, written only when explicitly requested
with --write-output:
  alternative_splice_candidate_audit_gff3
    Diagnostic browser track for all uniquely owned structural candidates; not a retained gene annotation.

Common final/audit outputs:
  run_configuration_json
    Deterministic JSON containing resolved inputs, settings, profile, pruning
    resolution, and final enabled-output keys. Included in every profile.
  final_selected_models_tsv
  final_selected_models_gff3
  final_pruned_selected_models_tsv
  final_pruned_selected_models_gff3
  final_pruned_selected_removed_tsv
  final_shadow_gene_pairs_tsv
  standardised_primary_models_gff3
    Standardised primary annotation with deterministic MoSelecta gene/transcript IDs and representative-model provenance.
  final_annotation_gff3
    Final annotation with stable primaries and accepted alternative transcripts.
  terminal_stop_canonicalisation_audit_tsv
    Exhaustive early model-level terminal-stop assessment; full/debug only.
  terminal_stop_canonicalisation_actions_tsv
    Canonicalised models only; standard/full/debug.
  terminal_stop_canonicalisation_summary_tsv
    Compact aggregate counts; every profile.
  cds_5prime_boundary_canonicalisation_audit_tsv
    Exhaustive early class-B assessment; full/debug only.
  cds_5prime_boundary_canonicalisation_actions_tsv
    Canonicalised models only; standard/full/debug.
  cds_5prime_boundary_canonicalisation_summary_tsv
    Compact aggregate counts; every profile.
  cds_phase_recovery_audit_tsv
    Exhaustive early phase-recovery assessment; full/debug only.
  cds_phase_recovery_actions_tsv
    Recovered missing-source-phase models only; standard/full/debug.
  cds_phase_recovery_summary_tsv
    Compact aggregate counts; every profile.
  definite_internal_stop_exclusion_summary_tsv
    Compact early hard-exclusion counts; every profile.
  definite_internal_stop_exclusions_tsv
    Ledger of candidates excluded for definite internal stops; standard/full/debug.
  extreme_short_cds_exclusion_summary_tsv
    Compact early coding-length-floor counts; every profile.
  extreme_short_cds_exclusions_tsv
    Ledger of candidates excluded below the coding-length floor; standard/full/debug.
  cds_terminus_refinement_audit_tsv
    One row per explicitly incomplete CDS side considered for post-selection refinement.
  cds_terminus_refinement_summary_tsv
    Compact run-level CDS terminus refinement counts and failure reasons.
  utr_augmented_final_annotation_gff3
    Separate post-selection annotation reconstructed from authoritative CDSs plus selected UTR evidence.
  utr_augmentation_audit_tsv
    One compact decision/audit row per final public transcript when --add-utrs is enabled.
  utr_junction_validation_audit_tsv
    Exception-only inherited and augmentation-candidate UTR validation decisions.
  utr_junction_validation_summary_tsv
    Compact aggregate UTR validation counts by evidence origin.
  utr_augmentation_summary_tsv
    One run-level UTR catalogue, assignment, filtering, learned-limit, and runtime summary row.
  pruning_scenario_summary_tsv
  auto_pruning_diagnostics_tsv
  selected_models_tsv
  selected_models_gff3
  backfilled_models_tsv
  backfill_rejected_models_tsv
  r1_backfill_replacements_tsv
    Successful strict-R1 replacements with candidate, blocker, score, corroboration and conflict provenance; automatically enabled only when R1 is active.
  r1_backfill_summary_tsv
    Compact strict-R1 evaluation, acceptance and rejection-reason counts; automatically enabled only when R1 is active.
  model_structural_support_tsv
  structural_corroboration_tsv
    Candidate-wide diagnostic of binary source-collapsed complete-structure support and run-level normalization alternatives; full/debug or explicit request only. The metric can be scored independently of writing this large table.

Use --write-output NAME to add a construction diagnostic to a standard run, for example:
  --write-output protein_evidence_collapsed_alignments_tsv
  --write-output transcript_evidence_collapsed_transcripts_tsv
"""

def build_arg_parser(
    full_help: bool = False,
    *,
    expert_help: bool = False,
) -> argparse.ArgumentParser:
    detailed_help = full_help or expert_help
    parser = argparse.ArgumentParser(
        prog="moselecta.py",
        description=(
            f"{PROGRAM_NAME} v{PROGRAM_VERSION}: CDS-structure-aware coding gene-model selection from multiple "
            "GFF/GFF3 annotation sources. The default pipeline scores all candidates, performs "
            "monosublocus de-overlap, runs superlocus-local backfill rescue, applies the selected "
            "final-pruning preset, and writes final selected models."
        ),
        epilog=full_help_epilog() if detailed_help else None,
        formatter_class=MoSelectaHelpFormatter,
    )
    parser.add_argument(
        "--full-help",
        action="store_true",
        help="Show advanced scientific and operational options, then exit.",
    )
    parser.add_argument(
        "--expert-help",
        action="store_true",
        help=(
            "Show every supported option, including low-level algorithm, "
            "individual output-path, and development controls, then exit."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"{PROGRAM_NAME} {PROGRAM_VERSION}",
        help="Show the MoSelecta version and exit.",
    )
    required = parser.add_argument_group("required inputs")
    required.add_argument(
        "--inputs",
        nargs="+",
        type=Path,
        required=True,
        help=(
            "Input GFF/GFF3 annotation files containing candidate gene models. "
            "Paths may include @FILE manifests; see --full-help."
        ),
    )

    input_group = parser.add_argument_group("input labelling")
    input_group.add_argument(
        "--source-name",
        nargs="+",
        default=None,
        help=(
            "Optional source/method labels, one per input file, e.g. tiberius annevo helixer. "
            "If omitted, file stems are used."
        ),
    )

    genome_group = parser.add_argument_group("genome sequence input")
    genome_group.add_argument(
        "--genome-fasta",
        type=Path,
        default=None,
        help=(
            "Optional uncompressed genomic FASTA used for start/stop codon, internal-stop, CDS-phase, "
            "canonical splice-motif, ambiguity, and exact contig-boundary metrics. A missing FILE.fai "
            "index is built automatically. When omitted, sequence-integrity scoring is inactive and "
            "current coordinate-only scoring is unchanged."
        ),
    )
    genome_group.add_argument(
        "--genetic-code",
        type=int,
        default=1,
        help=f"NCBI genetic code used for start/stop assessment. {PROGRAM_NAME} {PROGRAM_VERSION} currently supports code 1.",
    )
    genome_group.add_argument(
        "--definite-internal-stop-exclusion",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "With --genome-fasta, exclude resolved-frame models carrying a literal internal stop and unresolved models carrying a literal internal stop in all three frames before locus construction. Enabled by default; use --no-definite-internal-stop-exclusion for historical reproduction or diagnostics. Without a genome FASTA the stage is inactive."
        ),
    )
    eligibility_group = parser.add_argument_group("early biological eligibility")
    eligibility_group.add_argument(
        "--minimum-coding-length-bp",
        type=int,
        default=30,
        metavar="INT",
        help=(
            "Hard pre-locus eligibility floor on coding nucleotides after the accepted "
            "early CDS corrections. A canonical terminal stop established inside the CDS "
            "contributes 3 bp to raw CDS length but not coding length. Exactly 30 coding bp "
            "is eligible. Enabled for all sources by default, with or without genome FASTA. "
            "This is distinct from --min-cds-length-bp, the 150-bp scoring threshold."
        ),
    )


    evidence_group = parser.add_argument_group("protein evidence inputs")
    evidence_proteins_help = (
        "Protein evidence GFF/GTF file(s) containing CDS features. Comma-separated values, repeated arguments and @FILE manifests are accepted. "
        "When supplied, protein evidence contributes to model selection, can rescue supported singleton models, and enables evidence-aware low-score pruning under medium/high/very-high presets. Use --full-help for details."
    )
    if detailed_help:
        evidence_proteins_help = (
            "Protein evidence GFF/GTF file(s) containing CDS features. Comma-separated values, repeated arguments and @FILE manifests are accepted. "
            "Each file is treated as an independent protein evidence source. Protein evidence is parsed as CDS alignment evidence, collapsed by exact CDS-exon chain within same-strand protein evidence superloci, and compared to candidate models. "
            "When supplied, protein evidence contributes to final selection scoring by default, can rescue protein-supported singleton models from final pruning, and can enable evidence-aware low-score singleton/nonsingleton pruning under medium/high/very-high pruning presets. "
            "Use --selection-weight-protein-evidence 0, --no-final-prune-rescue-protein-supported-singletons, --no-final-prune-low-score-singletons, and/or --no-final-prune-low-score-nonsingletons to disable specific behaviours."
        )
    evidence_group.add_argument(
        "--evidence-proteins",
        action="append",
        default=None,
        metavar="GFF[,GFF...]",
        help=evidence_proteins_help,
    )
    evidence_group.add_argument(
        "--evidence-protein-source-name",
        action="append",
        default=None,
        metavar="LABEL[,LABEL...]",
        help=(
            "Optional labels for --evidence-proteins files, one per file after comma/repeated expansion. "
            "If omitted, order-based labels protein_001, protein_002, ... are used."
        ),
    )

    transcript_evidence_group = parser.add_argument_group("transcript evidence inputs")
    transcript_evidence_help = (
        "Transcript evidence GFF/GTF file(s) containing exon features. Comma-separated values, repeated arguments and @FILE manifests are accepted. "
        "Transcript evidence contributes to ranking by default with weight 2.0, and can participate in rescue/pruning; "
        "set --selection-weight-transcript-evidence 0 for diagnostic-only transcript scoring."
    )
    if detailed_help:
        transcript_evidence_help = (
            "Transcript evidence GFF/GTF file(s) containing exon features. Comma-separated values, repeated arguments and @FILE manifests are accepted. "
            "Each file is treated as an independent transcript evidence source. Only exon features are used; CDS, UTR, gene, and transcript container features are ignored for evidence construction. "
            "Raw transcripts are grouped from exon Parent/transcript_id attributes, same-strand transcript superloci are built from exon-overlap connected components, multi-exon transcripts are collapsed by exact intron chain, and single-exon transcripts are collapsed by exact exon interval. "
            "Transcript evidence is associated to candidate models using broad diagnostic and tighter scoring-eligible relationship classes. It contributes to ranking by default with weight 2.0; set --selection-weight-transcript-evidence 0 for diagnostic-only scoring. It can rescue strongly transcript-supported singleton models and participate in evidence-aware low-score pruning when transcript adequacy passes."
        )
    transcript_evidence_group.add_argument(
        "--evidence-transcripts",
        action="append",
        default=None,
        metavar="GFF[,GFF...]",
        help=transcript_evidence_help,
    )
    transcript_evidence_group.add_argument(
        "--evidence-transcript-source-name",
        action="append",
        default=None,
        metavar="LABEL[,LABEL...]",
        help=(
            "Optional labels for --evidence-transcripts files, one per file after comma/repeated expansion. "
            "If omitted, order-based labels transcript_001, transcript_002, ... are used."
        ),
    )
    evidence_group.add_argument(
        "--protein-evidence-terminal-cds-3prime-wobble-bp",
        type=int,
        default=None,
        help=(
            "Strand-aware outer-boundary wobble for 3' terminal CDS comparisons to protein evidence. "
            "If omitted, uses --terminal-cds-3prime-wobble-bp. Default effective value is 3 bp."
        ),
    )
    evidence_group.add_argument(
        "--protein-evidence-support-score-threshold",
        type=float,
        default=0.50,
        help="Minimum source-best protein evidence score used to count a protein source as supporting a model.",
    )

    progress_group = parser.add_argument_group("progress logging")
    progress_group.add_argument(
        "--progress-log-level",
        choices=("none", "basic", "verbose", "debug"),
        default="basic",
        help=(
            "Progress messages written to stderr during long runs. "
            "basic reports stage boundaries and periodic monosublocus de-overlap/backfill progress; "
            "verbose adds per-input/per-group messages; debug adds finer-grained diagnostics."
        ),
    )
    progress_group.add_argument(
        "--progress-interval",
        type=int,
        default=10000,
        help="Approximate item interval for periodic progress messages in long loops.",
    )
    progress_group.add_argument(
        "--debug-progress",
        action="store_true",
        help="Shortcut for --progress-log-level debug.",
    )

    pipeline_group = parser.add_argument_group("pipeline control")
    pipeline_group.add_argument(
        "--stop-after",
        choices=STOP_AFTER_CHOICES,
        default="final",
        help=(
            "Run the pipeline through this named step, then stop. "
            "scored-candidates writes scored candidate outputs only; monosublocus adds sublocus-winner de-overlap; "
            "backfill adds superlocus-local backfill rescue; final is the complete default pipeline."
        ),
    )
    pipeline_group.add_argument(
        "--add-utrs",
        action="store_true",
        help=(
            "After final public transcript IDs are frozen, select supported 5-prime and/or "
            "3-prime UTR structure and write a separate UTR-augmented final annotation. "
            "This switch does not control independent post-selection CDS terminus refinement."
        ),
    )
    pipeline_group.add_argument(
        "--refine-cds-termini",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Automatically attempt bounded post-selection repair of explicitly missing "
            "CDS starts/stops when --genome-fasta is available. This is independent of "
            "--add-utrs; use --no-refine-cds-termini to disable."
        ),
    )
    pipeline_group.add_argument(
        "--monosublocus-min-cds-exon-overlap-fraction",
        type=float,
        default=0.10,
        help=(
            "Non-reciprocal CDS-exon overlap threshold for monosublocus de-overlap and backfill conflict checks. "
            "An edge/conflict is formed when max(overlap_bp / CDS_bp_A, overlap_bp / CDS_bp_B) is at least this value. "
            "Overlap is measured over CDS exons, not CDS span."
        ),
    )
    pipeline_group.add_argument(
        "--strict-r1-backfill-replacement",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "During backfill, allow an eligible corroborated multi-exon "
            "candidate to atomically replace one or two lower-scoring, uncorroborated single-exon "
            "models selected at the monosublocus stage. Enabled by default."
        ),
    )
    pipeline_group.add_argument(
        "--strict-r1-allow-multiexon-blockers",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Allow strict R1 to replace eligible "
            "uncorroborated multi-exon monosublocus holder winners as well as "
            "single-exon winners. Requires --strict-r1-backfill-replacement. "
            "Enabled by default when strict R1 is enabled."
        ),
    )
    pipeline_group.add_argument(
        "--alternative-splice-report",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Write detailed alternative-splice selector TSV diagnostics. Selection is already "
            "performed by default when the pipeline reaches the final annotation stage."
        ),
    )
    pipeline_group.add_argument(
        "--no-alternative-splice",
        action="store_true",
        default=False,
        help=(
            "Write the final annotation with frozen primary .1 transcripts only. "
            "Public gene and .1 identities remain unchanged."
        ),
    )
    pipeline_group.add_argument(
        "--alternative-splice-preset",
        choices=("permissive", "intermediate", "stringent"),
        default="intermediate",
        help=(
            "Alternative-splice numerical threshold preset used for the final annotation "
            "or optional diagnostics."
        ),
    )
    pipeline_group.add_argument(
        "--alternative-splice-support-policy",
        choices=("exact", "exact-plus-mosaic"),
        default="exact",
        help=(
            "Alternative-splice biological support policy. exact requires an exact model, "
            "protein, or coherent transcript route; exact-plus-mosaic additionally permits "
            "a full structural-mosaic model route meeting the active prediction-source "
            "requirement."
        ),
    )
    pipeline_group.add_argument(
        "--alternative-splice-evidence-requirement",
        choices=(
            "none", "protein-or-transcript", "protein", "transcript",
            "protein-and-transcript",
        ),
        default="none",
        help=(
            "Optional qualifying exact external-evidence requirement, applied independently "
            "of the support policy. Strong-partial and aggregate evidence do not satisfy it."
        ),
    )
    pipeline_group.add_argument(
        "--alternative-splice-max-alternatives-per-gene",
        type=nonnegative_int,
        default=50,
        help=(
            "Maximum additional transcripts accepted per frozen primary gene by the "
            "alternative-splice selector. The primary .1 is not counted."
        ),
    )
    pipeline_group.add_argument(
        "--alternative-splice-min-model-support-sources",
        type=minimum_two_int,
        default=2,
        metavar="INT",
        help=(
            "Minimum distinct prediction-model sources required for exact or "
            "full structural-mosaic model support. The candidate source counts."
        ),
    )
    pipeline_group.add_argument(
        "--alternative-splice-min-model-support-source-fraction",
        type=positive_unit_interval_float,
        default=None,
        metavar="FLOAT",
        help=(
            "Optional fraction of all prediction-model sources required for model-based "
            "alternative-splice support. The effective count is the stricter of this "
            "global fraction, the absolute minimum, and the floor of two sources."
        ),
    )
    pipeline_group.add_argument(
        "--alternative-splice-min-reciprocal-cds-overlap",
        type=unit_interval_float,
        default=None,
        metavar="FLOAT",
        help="Override only the active preset's inclusive reciprocal CDS-overlap threshold.",
    )
    pipeline_group.add_argument(
        "--alternative-splice-min-frame-concordance",
        type=unit_interval_float,
        default=None,
        metavar="FLOAT",
        help="Override only the active preset's inclusive same-frame-overlap threshold.",
    )
    pipeline_group.add_argument(
        "--alternative-splice-min-model-selection-score-scaled",
        type=unit_interval_float,
        default=None,
        metavar="FLOAT",
        help="Override only the active preset's inclusive scaled model-selection-score threshold.",
    )
    pipeline_group.add_argument(
        "--alternative-splice-min-primary-score-ratio",
        type=unit_interval_float,
        default=None,
        metavar="FLOAT",
        help=(
            "Override only the active preset's inclusive raw candidate/positive-primary score "
            "ratio threshold. Non-positive primary scores are unavailable and rejected."
        ),
    )

    output_group = parser.add_argument_group("outputs")
    output_group.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("moselecta"),
        help="Prefix for output files. Default writes files beginning with moselecta.",
    )
    output_group.add_argument(
        "--output-profile",
        choices=sorted(OUTPUT_PROFILES),
        default="standard",
        help=(
            "Controls reporting volume. minimal writes final annotation products and essential compact summaries; "
            "standard adds primary annotation, selection/action ledgers, and action-focused biological reports; "
            "full writes exhaustive candidate, evidence, selection, and audit diagnostics; debug is currently an "
            "alias of full. Individual products can be controlled with --write-output and --skip-output."
        ),
    )
    output_group.add_argument(
        "--write-output",
        action="append",
        default=None,
        metavar="NAME[,NAME...]",
        help=(
            "Force named output(s) to be written in addition to the selected --output-profile. "
            "Use --full-help for the output profile/key guide."
        ),
    )
    output_group.add_argument(
        "--skip-output",
        action="append",
        default=None,
        metavar="NAME[,NAME...]",
        help=(
            "Suppress named output(s) even if enabled by the profile, --write-output, a feature, or a stop stage. "
            "Skip is applied last and always wins; use all to suppress every file. Repeatable and comma-separated."
        ),
    )
    output_group.add_argument(
        "--final-gff3-cds-only-structure",
        action="store_true",
        help=(
            "For ordinary final/public GFF3 outputs, including final_selected_models.gff3, "
            "final_selected_models.pruned.gff3, standardised_primary_models.gff3, and "
            "final_annotation.gff3, write gene/mRNA/exon spans from CDS coordinates, "
            "generate exon features from CDS parts, and omit UTR/non-CDS child features. "
            "Default preserves the original model child features and any UTR-extended spans."
        ),
    )
    output_group.add_argument(
        "--combined-gff3",
        type=Path,
        default=None,
        help="Combined GFF3 with superlocus features and candidate coding models.",
    )
    output_group.add_argument(
        "--scored-pruned-gff3",
        type=Path,
        default=None,
        help=(
            "Optional second scored-candidate GFF3 after applying singleton-sublocus pruning. "
            "If omitted and pruning is enabled, uses <out-prefix>.coding_models_with_superloci.no_singleton_subloci.gff3."
        ),
    )
    output_group.add_argument(
        "--superloci-tsv",
        type=Path,
        default=None,
        help="TSV summary of superloci.",
    )
    output_group.add_argument(
        "--members-tsv",
        type=Path,
        default=None,
        help="TSV mapping superloci to input gene/transcript models.",
    )
    output_group.add_argument(
        "--models-tsv",
        type=Path,
        default=None,
        help="TSV summary of all parsed coding models.",
    )
    output_group.add_argument(
        "--excluded-models-tsv",
        type=Path,
        default=None,
        help="TSV report of CDS models excluded due to mixed chromosomes/strands.",
    )
    output_group.add_argument(
        "--other-features-gff3",
        type=Path,
        default=None,
        help="GFF3 containing input features not used as coding candidate models.",
    )
    output_group.add_argument(
        "--input-summary-tsv",
        type=Path,
        default=None,
        help="TSV summary of each input file.",
    )
    output_group.add_argument(
        "--validation-failures-tsv",
        type=Path,
        default=None,
        help="TSV report of model-level and superlocus-level validation failures.",
    )
    output_group.add_argument(
        "--subloci-tsv",
        type=Path,
        default=None,
        help="TSV summary of first-pass subloci.",
    )
    output_group.add_argument(
        "--sublocus-members-tsv",
        type=Path,
        default=None,
        help="TSV mapping subloci to input gene/transcript models.",
    )
    output_group.add_argument(
        "--relationship-edges-tsv",
        type=Path,
        default=None,
        help="TSV report of model-model relationships used to build first-pass subloci.",
    )
    output_group.add_argument(
        "--model-metrics-tsv",
        type=Path,
        default=None,
        help="TSV report of coordinate-only intrinsic metrics for each coding model.",
    )
    output_group.add_argument(
        "--sequence-integrity-summary-tsv",
        type=Path,
        default=None,
        help="Compact TSV summary of genome-derived sequence-integrity metrics when --genome-fasta is supplied.",
    )
    output_group.add_argument(
        "--sequence-integrity-model-metrics-tsv",
        type=Path,
        default=None,
        help="Model-level start/stop, internal-stop, phase, splice-motif, ambiguity, boundary, and sequence-penalty metrics.",
    )
    output_group.add_argument(
        "--source-auto-weights-tsv",
        type=Path,
        default=None,
        help=(
            "TSV report of source-level coordinate agreement metrics and optional "
            "manual/automatic source-prior weights."
        ),
    )

    output_group.add_argument(
        "--model-structural-support-tsv",
        type=Path,
        default=None,
        help="TSV report of per-model structural support metrics within subloci, plus consensus/pairwise similarity and completeness/containment diagnostics.",
    )
    output_group.add_argument(
        "--structural-corroboration-tsv",
        type=Path,
        default=None,
        help=(
            "Optional candidate-wide diagnostic of binary, source-collapsed whole-gene structural "
            "corroboration and run-level normalization alternatives. Enable in a standard run with "
            "--write-output structural_corroboration_tsv. Writing this diagnostic is independent of "
            "the optional corroboration scoring component."
        ),
    )
    output_group.add_argument(
        "--structural-coherence-regions-tsv",
        type=Path,
        default=None,
        help=(
            "Region-level audit of candidate structural-coherence evidence, source-collapsed "
            "anchor states, scoring gates, and penalty contributions."
        ),
    )
    output_group.add_argument(
        "--model-selection-scores-tsv",
        type=Path,
        default=None,
        help=(
            "TSV report of final model selection scores, including model-only structural-support score, "
            "optional evidence scores, final model_selection_score, ranks, and evidence-driven rank-change diagnostics."
        ),
    )
    output_group.add_argument(
        "--model-selection-score-summary-tsv",
        type=Path,
        default=None,
        help=(
            "Compact TSV summary of model-only versus final model selection scores, "
            "protein-evidence integration, and evidence-driven rank changes."
        ),
    )
    output_group.add_argument(
        "--structural-features-tsv",
        type=Path,
        default=None,
        help="TSV report of sublocus-level structural feature support used by structural support metrics.",
    )
    output_group.add_argument(
        "--sublocus-winners-gff3",
        type=Path,
        default=None,
        help="Optional GFF3 of one provisional winner per current sublocus; written when --stop-after is monosublocus, backfill, or final.",
    )
    output_group.add_argument(
        "--sublocus-winners-tsv",
        type=Path,
        default=None,
        help="Optional TSV of one provisional winner per current sublocus; written when --stop-after is monosublocus, backfill, or final.",
    )
    output_group.add_argument(
        "--monosubloci-tsv",
        type=Path,
        default=None,
        help="Optional TSV summary of monosublocus overlap holders among sublocus winners.",
    )
    output_group.add_argument(
        "--monosublocus-members-tsv",
        type=Path,
        default=None,
        help="Optional TSV of sublocus winners assigned to monosublocus holders and greedy selected/suppressed state.",
    )
    output_group.add_argument(
        "--monosublocus-edges-tsv",
        type=Path,
        default=None,
        help="Optional TSV of CDS-exon-overlap edges used to build monosublocus holders.",
    )
    output_group.add_argument(
        "--monosublocus-review-gff3",
        type=Path,
        default=None,
        help="Optional GFF3 containing all sublocus winners with monosublocus selected/suppressed attributes.",
    )
    output_group.add_argument(
        "--selected-models-gff3",
        type=Path,
        default=None,
        help="Optional GFF3 containing the final greedy non-overlapping selected models from monosublocus holders.",
    )
    output_group.add_argument(
        "--selected-models-tsv",
        type=Path,
        default=None,
        help="Optional TSV containing the final greedy non-overlapping selected models from monosublocus holders.",
    )
    output_group.add_argument(
        "--suppressed-monosublocus-models-tsv",
        type=Path,
        default=None,
        help="Optional TSV of sublocus winners suppressed by higher-ranked overlapping selected models.",
    )
    output_group.add_argument(
        "--backfill-candidates-tsv",
        type=Path,
        default=None,
        help="Optional TSV of all non-selected, non-suppressed candidates considered during backfill selection.",
    )
    output_group.add_argument(
        "--backfilled-models-tsv",
        type=Path,
        default=None,
        help="Optional TSV of candidates accepted during backfill selection.",
    )
    output_group.add_argument(
        "--backfill-rejected-models-tsv",
        type=Path,
        default=None,
        help="Optional TSV of backfill candidates rejected because they conflict with an already selected model in the same superlocus.",
    )
    output_group.add_argument(
        "--r1-backfill-replacements-tsv",
        type=Path,
        default=None,
        help="Append-only TSV provenance for successful strict R1 backfill replacements.",
    )
    output_group.add_argument(
        "--r1-backfill-summary-tsv",
        type=Path,
        default=None,
        help="Compact strict R1 counts and rejection-reason summary.",
    )
    output_group.add_argument(
        "--final-selected-models-tsv",
        type=Path,
        default=None,
        help="Optional TSV of the final selected model set after monosublocus selection plus accepted backfill models.",
    )
    output_group.add_argument(
        "--final-selected-models-gff3",
        type=Path,
        default=None,
        help="Optional GFF3 of the final selected model set after monosublocus selection plus accepted backfill models.",
    )
    output_group.add_argument(
        "--final-pruned-selected-models-tsv",
        type=Path,
        default=None,
        help="Optional TSV path for final selected models retained after final selected-model pruning.",
    )
    output_group.add_argument(
        "--final-pruned-selected-models-gff3",
        type=Path,
        default=None,
        help="Optional GFF3 path for final selected models retained after final selected-model pruning.",
    )
    output_group.add_argument(
        "--final-pruned-selected-removed-tsv",
        type=Path,
        default=None,
        help="Optional TSV path for final selected models removed by final selected-model pruning.",
    )
    output_group.add_argument(
        "--final-shadow-gene-pairs-tsv",
        type=Path,
        default=None,
        help="Optional TSV path for opposite-strand shadow-gene candidate pairs considered during final selected-model pruning.",
    )
    output_group.add_argument(
        "--standardised-primary-models-gff3",
        type=Path,
        default=None,
        help=(
            "Additional standardised primary GFF3 with deterministic MoSelecta IDs and "
            "provenance for the actual retained representative models. The legacy final GFF3 "
            "outputs are unchanged."
        ),
    )
    output_group.add_argument(
        "--final-annotation-gff3",
        type=Path,
        default=None,
        help=(
            "Final annotation GFF3 containing stable primary .1 transcripts and, unless "
            "disabled, accepted alternative transcripts."
        ),
    )
    output_group.add_argument(
        "--utr-augmented-final-annotation-gff3",
        type=Path,
        default=None,
        help="Separate final annotation reconstructed with selected post-selection UTRs.",
    )
    output_group.add_argument(
        "--run-configuration-json",
        type=Path,
        default=None,
        help="Deterministic JSON record of resolved inputs, settings, profile, and enabled output keys.",
    )
    output_group.add_argument(
        "--terminal-stop-canonicalisation-summary-tsv",
        type=Path,
        default=None,
        help="Compact run-level terminal-stop canonicalisation summary.",
    )
    output_group.add_argument(
        "--terminal-stop-canonicalisation-actions-tsv",
        type=Path,
        default=None,
        help="Terminal-stop canonicalisation rows for models actually modified.",
    )
    output_group.add_argument(
        "--cds-5prime-boundary-canonicalisation-audit-tsv",
        type=Path,
        default=None,
        help=(
            "Model-level audit of bounded class-B 5-prime CDS-boundary "
            "canonicalisation, including exon context, trim, phase, and sequence checks."
        ),
    )
    output_group.add_argument(
        "--cds-5prime-boundary-canonicalisation-summary-tsv",
        type=Path,
        default=None,
        help="Compact run-level 5-prime CDS-boundary canonicalisation summary.",
    )
    output_group.add_argument(
        "--cds-5prime-boundary-canonicalisation-actions-tsv",
        type=Path,
        default=None,
        help="5-prime CDS-boundary canonicalisation rows for models actually modified.",
    )
    output_group.add_argument(
        "--terminal-stop-canonicalisation-audit-tsv",
        type=Path,
        default=None,
        help=(
            "Model-level audit of early local terminal-stop boundary canonicalisation, "
            "including offsets, frame basis, skips, and newly introduced candidate overlaps."
        ),
    )
    output_group.add_argument(
        "--cds-phase-recovery-audit-tsv",
        type=Path,
        default=None,
        help=(
            "Model-level audit of pre-integrity biological CDS phase recovery "
            "for source first-phase '.' models."
        ),
    )
    output_group.add_argument(
        "--cds-phase-recovery-summary-tsv",
        type=Path,
        default=None,
        help="Compact run-level biological CDS phase-recovery summary.",
    )
    output_group.add_argument(
        "--cds-phase-recovery-actions-tsv",
        type=Path,
        default=None,
        help="CDS phase-recovery rows for models whose missing source phase was recovered.",
    )
    output_group.add_argument(
        "--definite-internal-stop-exclusion-summary-tsv",
        type=Path,
        default=None,
        help="Compact run-level summary of early definite internal-stop exclusion.",
    )
    output_group.add_argument(
        "--definite-internal-stop-exclusions-tsv",
        type=Path,
        default=None,
        help="Candidate ledger for early definite internal-stop exclusions.",
    )
    output_group.add_argument(
        "--extreme-short-cds-exclusion-summary-tsv",
        type=Path,
        default=None,
        help="Compact run-level summary of the early hard coding-length floor.",
    )
    output_group.add_argument(
        "--extreme-short-cds-exclusions-tsv",
        type=Path,
        default=None,
        help="Candidate ledger for early extreme-short CDS exclusions.",
    )
    output_group.add_argument(
        "--cds-terminus-refinement-audit-tsv",
        type=Path,
        default=None,
        help="One-row-per-incomplete-side post-selection CDS refinement audit.",
    )
    output_group.add_argument(
        "--cds-terminus-refinement-summary-tsv",
        type=Path,
        default=None,
        help="Compact run-level CDS terminus refinement summary.",
    )
    output_group.add_argument(
        "--utr-augmentation-audit-tsv",
        type=Path,
        default=None,
        help="One-row-per-final-transcript UTR selection audit.",
    )
    output_group.add_argument(
        "--utr-junction-validation-audit-tsv",
        type=Path,
        default=None,
        help="Exception-only inherited and candidate UTR validation audit.",
    )
    output_group.add_argument(
        "--utr-junction-validation-summary-tsv",
        type=Path,
        default=None,
        help="Compact aggregate inherited and candidate UTR validation summary.",
    )
    output_group.add_argument(
        "--utr-augmentation-summary-tsv",
        type=Path,
        default=None,
        help="Compact run-level UTR augmentation summary.",
    )
    output_group.add_argument(
        "--pruning-scenario-summary-tsv",
        type=Path,
        default=None,
        help="Optional stable one-row-per-level summary of prospective medium, high, and very-high pruning scenarios.",
    )
    output_group.add_argument(
        "--auto-pruning-diagnostics-tsv",
        type=Path,
        default=None,
        help="Optional long-format automatic-pruning evidence, consensus, validation, and decision diagnostics.",
    )
    output_group.add_argument(
        "--alternative-splice-candidates-tsv",
        type=Path,
        default=None,
        help="Candidate/retained-gene alternative-splice selector diagnostic TSV.",
    )
    output_group.add_argument(
        "--alternative-splice-genes-tsv",
        type=Path,
        default=None,
        help="Per-retained-primary alternative-splice selector summary TSV.",
    )
    output_group.add_argument(
        "--alternative-splice-run-summary-tsv",
        type=Path,
        default=None,
        help="Run-level alternative-splice selector diagnostic summary TSV.",
    )
    output_group.add_argument(
        "--alternative-splice-candidate-audit-gff3",
        type=Path,
        default=None,
        help=(
            "Diagnostic browser GFF3 for uniquely owned structural candidates and "
            "active selector decisions; not a retained gene annotation."
        ),
    )

    output_group.add_argument(
        "--protein-evidence-summary-tsv",
        type=Path,
        default=None,
        help="Optional compact TSV summary of protein evidence inputs, grouped loci, and model-support status counts.",
    )
    output_group.add_argument(
        "--protein-evidence-raw-to-collapsed-tsv",
        type=Path,
        default=None,
        help="Optional TSV mapping raw protein evidence alignments to exact-chain collapsed protein evidence alignments.",
    )
    output_group.add_argument(
        "--protein-evidence-collapsed-alignments-tsv",
        type=Path,
        default=None,
        help="Optional TSV of exact-CDS-chain collapsed protein evidence alignments.",
    )
    output_group.add_argument(
        "--protein-evidence-superloci-tsv",
        type=Path,
        default=None,
        help="Optional TSV summary of protein evidence superloci.",
    )
    output_group.add_argument(
        "--protein-evidence-subloci-tsv",
        type=Path,
        default=None,
        help="Optional TSV summary of protein evidence subloci built from collapsed alignments.",
    )
    output_group.add_argument(
        "--protein-evidence-sublocus-members-tsv",
        type=Path,
        default=None,
        help="Optional TSV mapping protein evidence subloci to collapsed alignments.",
    )
    output_group.add_argument(
        "--protein-evidence-model-support-tsv",
        type=Path,
        default=None,
        help="Optional TSV of diagnostic protein evidence support metrics for candidate models.",
    )

    output_group.add_argument(
        "--transcript-evidence-summary-tsv",
        type=Path,
        default=None,
        help="Optional compact TSV summary of transcript evidence inputs, exact-chain collapse, and transcript evidence loci.",
    )
    output_group.add_argument(
        "--transcript-evidence-raw-to-collapsed-tsv",
        type=Path,
        default=None,
        help="Optional TSV mapping raw transcript evidence records to collapsed transcript evidence records. Potentially large; normally written only with --output-profile full/debug or explicit --write-output.",
    )
    output_group.add_argument(
        "--transcript-evidence-collapsed-transcripts-tsv",
        type=Path,
        default=None,
        help="Optional TSV of exact-intron-chain / exact-single-exon collapsed transcript evidence records.",
    )
    output_group.add_argument(
        "--transcript-evidence-superloci-tsv",
        type=Path,
        default=None,
        help="Optional TSV summary of transcript evidence superloci built from raw same-strand exon-overlap components.",
    )
    output_group.add_argument(
        "--transcript-evidence-subloci-tsv",
        type=Path,
        default=None,
        help="Optional TSV summary of transcript evidence subloci built from collapsed transcript evidence records.",
    )
    output_group.add_argument(
        "--transcript-evidence-sublocus-members-tsv",
        type=Path,
        default=None,
        help="Optional TSV mapping transcript evidence subloci to collapsed transcript evidence records.",
    )
    output_group.add_argument(
        "--transcript-evidence-model-support-tsv",
        type=Path,
        default=None,
        help="Optional TSV of diagnostic transcript evidence broad/tight association classes for candidate models.",
    )

    gff3_filter_group = parser.add_argument_group("scored-candidate GFF3 pruning filters")
    gff3_filter_group.add_argument(
        "--prune-singleton-subloci-from-scored-gff3",
        action="store_true",
        help=(
            "Write an additional all-candidate scored GFF3 with singleton subloci removed. This affects only the scored-candidate GFF3; TSV scoring reports and selected-model outputs are unchanged."
        ),
    )
    gff3_filter_group.add_argument(
        "--scored-gff3-prune-singleton-sources",
        action="append",
        default=None,
        metavar="SOURCE[,SOURCE...]",
        help=(
            "Restrict scored-candidate GFF3 singleton pruning to subloci whose only model comes from the named source(s), e.g. Helixer or Helixer,Tiberius. Repeatable. If omitted, singleton subloci from all sources are pruned."
        ),
    )
    gff3_filter_group.add_argument(
        "--scored-gff3-keep-empty-superloci",
        action="store_true",
        help="In the pruned scored-candidate GFF3, keep superlocus features even if all child subloci were removed. Default: omit empty superloci.",
    )

    final_prune_group = parser.add_argument_group("final selected-model pruning")
    pruning_level_help = (
        "Final selected-model pruning level. Automatic modes prospectively evaluate medium, high, and "
        "very-high against the same selected population without a reference; auto is the default. "
        "Use --full-help for the fixed preset thresholds, safeguards, and advanced overrides."
    )
    if detailed_help:
        pruning_level_help = (
            "Shortcut preset for final selected-model pruning. none disables preset pruning and singleton rescue; "
            "low applies singleton-sublocus intrinsic pruning, enables singleton monosublocus-overlap rescue, "
            "and applies conservative opposite-strand shadow pruning requiring the loser to be single-exon; "
            "medium applies singleton intrinsic-or-sequence-integrity pruning when genome metrics are available, enables singleton monosublocus-overlap rescue, "
            "and applies shadow pruning without the single-exon loser requirement; "
            "high applies singleton all pruning, enables singleton monosublocus-overlap rescue, and applies shadow pruning. "
            "medium low-score thresholds are singleton 0.25/0.25/0.40/0.25 and nonsingleton "
            "0.20/0.20/0.40/0.20; high uses 0.35/0.35/0.50/0.35 for both; very-high uses "
            "0.40/0.40/0.60/0.40 for both (model/protein/transcript/provisional). All three use "
            "low-score nonsingleton scope=all. When no supplied evidence class is adequate, nonsingletons instead use canonical model-confidence thresholds: "
            "medium=0.20 and high/very-high=0.25, with 0.35 and 0.45 high/very-high low-depth tiers; this fallback requires at least three run model sources. "
            "Medium retains sequence-integrity scope=backfill, while "
            "high/very-high use all. Evidence adequacy gating and strong local rescue remain active. "
            "auto and auto-conservative simulate all three fixed scenarios and differ only in selection policy."
        )
    final_prune_group.add_argument(
        "--pruning-level",
        choices=PRUNING_LEVEL_CHOICES,
        default="auto",
        help=pruning_level_help,
    )
    final_prune_group.add_argument(
        "--auto-pruning-min-level",
        choices=("low", "medium", "high"),
        default="medium",
        help="Minimum fixed level that automatic pruning may select. Default: medium.",
    )
    final_prune_group.add_argument(
        "--auto-pruning-max-level",
        choices=("medium", "high", "very-high"),
        default=None,
        help=(
            "Maximum fixed level that automatic pruning may select. Defaults to very-high for auto "
            "and high for auto-conservative. An explicit very-high value is required before "
            "auto-conservative may consider very-high."
        ),
    )
    final_prune_group.add_argument(
        "--auto-pruning-max-extra-removal-fraction",
        type=float,
        default=0.05,
        help="Maximum very-high additional removals divided by high-retained models for auto. Default: 0.05.",
    )
    final_prune_group.add_argument(
        "--auto-pruning-high-min-singleton-fraction",
        type=float,
        default=0.10,
        help="Minimum selected singleton-model fraction used by auto when evidence adequacy passes. Default: 0.10.",
    )
    final_prune_group.add_argument(
        "--auto-pruning-high-min-prospective-removal-fraction",
        type=float,
        default=0.02,
        help="Minimum high-versus-medium additional-removal fraction used by auto. Default: 0.02.",
    )
    final_prune_group.add_argument(
        "--auto-pruning-use-source-excess",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use source-specific singleton excess when resolving an intermediate no-evidence consensus backbone.",
    )
    final_prune_group.add_argument(
        "--auto-pruning-report-only",
        action="store_true",
        help=(
            "With a fixed pruning level, preserve that output while simulating medium/high/very-high "
            "and reporting both automatic recommendations."
        ),
    )
    final_prune_group.add_argument(
        "--final-prune-singleton-subloci-level",
        choices=("none", "intrinsic", "sequence_integrity", "intrinsic_or_sequence_integrity", "all"),
        default=None,
        help=(
            "Advanced override for singleton-sublocus final pruning. If omitted, this is set by --pruning-level. "
            "none disables this filter; intrinsic removes only singleton-sublocus final models with intrinsic penalties; "
            "sequence_integrity removes singletons with internal stops, missing starts, missing stops, or all splice junctions noncanonical; "
            "intrinsic_or_sequence_integrity applies either candidate class; all removes all singleton-sublocus final models. "
            "Protein/transcript singleton rescue remains available for sequence-integrity candidates."
        ),
    )
    final_prune_group.add_argument(
        "--final-prune-singleton-subloci-source",
        action="append",
        default=None,
        metavar="SOURCE[,SOURCE...]",
        help=(
            "Restrict final selected singleton-sublocus pruning to models from the named source(s), e.g. Helixer or Helixer,Tiberius. "
            "Repeatable. If omitted, singleton-sublocus final models from all sources are eligible for pruning."
        ),
    )
    final_prune_group.add_argument(
        "--final-prune-opposite-strand-shadows",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Advanced override for opposite-strand shadow-gene pruning. If omitted, this is set by --pruning-level. "
            "Shadow candidates are selected-model pairs on opposite strands with high CDS-exon overlap of the lower-scoring model."
        ),
    )
    final_prune_group.add_argument(
        "--shadow-prune-source",
        action="append",
        default=None,
        metavar="SOURCE[,SOURCE...]",
        help=(
            "Restrict opposite-strand shadow pruning to lower-scoring candidate losers from the named source(s). "
            "If omitted, models from all sources are eligible for removal."
        ),
    )
    final_prune_group.add_argument(
        "--shadow-prune-require-loser-single-exon",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Require the lower-scoring model in an opposite-strand shadow pair to be single-exon before it can be removed. "
            "If omitted, this is set by --pruning-level."
        ),
    )
    final_prune_group.add_argument(
        "--shadow-prune-min-overlap-bp",
        type=int,
        default=150,
        help="Minimum CDS-exon overlap bp required to consider an opposite-strand shadow pair.",
    )
    final_prune_group.add_argument(
        "--shadow-prune-min-loser-cds-overlap-fraction",
        type=float,
        default=0.80,
        help="Minimum fraction of the lower-scoring model's CDS that must be overlapped by the higher-scoring opposite-strand model.",
    )
    final_prune_group.add_argument(
        "--shadow-prune-min-score-margin-scaled",
        type=float,
        default=0.05,
        help="Minimum model_selection_score_scaled margin required before removing the lower-scoring opposite-strand shadow candidate.",
    )
    final_prune_group.add_argument(
        "--singleton-prune-rescue-monosublocus-overlap",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Advanced singleton-pruning rescue. If omitted, this is set by --pruning-level: disabled for none and enabled for low/medium/high/very-high/automatic modes. "
            "If a singleton-sublocus model would otherwise be removed, retain it when a multi-member monosublocus contains a sufficiently overlapping prediction from a different model source covering substantially the same CDS territory. Structural agreement is not required."
        ),
    )
    final_prune_group.add_argument(
        "--monosublocus-overlap-rescue-evidence-policy",
        choices=("unconditional", "adequacy-conditioned"),
        default="adequacy-conditioned",
        help=(
            "Control evidence requirements after monosublocus-overlap rescue geometry passes. "
            "adequacy-conditioned (default): strong local protein or transcript support always rescues; "
            "otherwise overlap-only rescue is denied when either evidence class passes global adequacy, "
            "and retained when no evidence class passes adequacy. unconditional bypasses only the "
            "evidence-adequacy condition; all configured geometry and different-source requirements still apply."
        ),
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-monosublocus-member-count",
        type=int,
        default=2,
        help="Minimum monosublocus holder member count required to rescue an otherwise pruned singleton-sublocus model. Default: 2.",
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-cds-overlap-bp",
        type=int,
        default=150,
        help="Minimum CDS-exon overlap bp to rescue an otherwise pruned singleton-sublocus model. Default: 150.",
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-self-overlap-fraction",
        type=float,
        default=0.80,
        help="Minimum fraction of the singleton model's CDS that must overlap another monosublocus holder member for rescue. Default: 0.80.",
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-reciprocal-overlap-fraction",
        type=float,
        default=0.80,
        help=(
            "Minimum of the singleton-side and partner-side CDS-exon overlap fractions required for monosublocus-overlap rescue. "
            "This requires a prediction from a different model source to cover substantially the same CDS territory; exon/intron structural agreement is not required. Default: 0.80."
        ),
    )
    final_prune_group.add_argument(
        "--final-prune-rescue-protein-supported-singletons",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Rescue otherwise pruned singleton-sublocus final models when they have strong protein evidence support. "
            "This evidence rescue is independent of monosublocus-overlap geometry. If omitted, this is enabled automatically when --evidence-proteins is supplied and disabled otherwise."
        ),
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-protein-evidence-score-scaled",
        type=float,
        default=0.75,
        help=(
            "Minimum protein_evidence_score_scaled required to rescue an otherwise pruned singleton-sublocus model. "
            "Default: 0.75."
        ),
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-protein-exact-chain-source-count",
        type=int,
        default=1,
        help=(
            "Minimum protein evidence exact CDS-chain supporting source count required to rescue an otherwise pruned singleton-sublocus model. "
            "Set to 0 to disable exact-chain based singleton rescue. Default: 1."
        ),
    )
    final_prune_group.add_argument(
        "--final-prune-rescue-transcript-supported-singletons",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Rescue otherwise pruned singleton-sublocus final models when they have strong scoring-eligible transcript evidence support. "
            "This evidence rescue is independent of monosublocus-overlap geometry. If omitted, this is enabled automatically when --evidence-transcripts is supplied and disabled otherwise. "
            "By default, an evaluable CDS/min-cDNA fraction below the configured minimum or a calculated "
            "transcript cDNA-context pruning concern blocks transcript-only rescue. Missing cDNA context is neutral."
        ),
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-transcript-evidence-score-scaled",
        type=float,
        default=0.75,
        help=(
            "Minimum transcript_evidence_score_scaled required to rescue an otherwise pruned singleton-sublocus model. "
            "Only scoring-eligible transcript associations contribute to this score. Default: 0.75."
        ),
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-transcript-supporting-source-count",
        type=int,
        default=1,
        help=(
            "Minimum transcript supporting-source count required to rescue an otherwise pruned singleton-sublocus model. "
            "Default: 1."
        ),
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-transcript-scoring-eligible-count",
        type=int,
        default=1,
        help=(
            "Minimum number of scoring-eligible collapsed transcript evidence records required to rescue an otherwise pruned singleton-sublocus model. "
            "Default: 1."
        ),
    )
    final_prune_group.add_argument(
        "--singleton-rescue-min-transcript-cdna-context-fraction",
        type=float,
        default=0.30,
        help=(
            "Minimum transcript_cdna_context_cds_to_min_cdna_fraction required for transcript-supported "
            "singleton rescue when cDNA context is evaluable. Context that cannot be evaluated is neutral, "
            "so the other transcript rescue criteria still apply. Default: 0.30."
        ),
    )
    final_prune_group.add_argument(
        "--transcript-singleton-rescue-require-stranded-support",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Require at least one strand-resolved scoring-eligible transcript association for transcript-supported singleton rescue. "
            "Default: false, so unstranded transcript evidence can contribute to rescue when other transcript criteria are met."
        ),
    )
    final_prune_group.add_argument(
        "--transcript-singleton-rescue-block-on-cdna-context-concern",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Block transcript-only singleton rescue when calculated cDNA context has a severe low CDS/min-cDNA fraction, "
            "or a moderate low fraction accompanied by a structural context flag. Uncalculable context is neutral and "
            "falls back to the remaining transcript rescue criteria. Default: true."
        ),
    )
    final_prune_group.add_argument(
        "--allow-internal-stop-evidence-rescue",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Advanced override. By default, a model carrying internal_stop_codon cannot be rescued or protected from sequence-integrity removal by protein or transcript evidence, including exact-chain protein support. "
            "Enable this option to restore that later evidence-based rescue/protection for retained internal-stop models. It cannot restore a candidate removed by the earlier definite internal-stop exclusion stage. Missing-start and missing-stop rescue is unaffected. Default: false."
        ),
    )
    final_prune_group.add_argument(
        "--final-prune-low-score-singletons",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Optional evidence-aware low-score pruning for selected singleton-sublocus models. "
            "If omitted, this is enabled by --pruning-level medium/high/very-high/automatic modes when protein or transcript evidence is supplied, "
            "and disabled otherwise. It is applied after monosublocus/protein/transcript singleton rescue and before opposite-strand shadow pruning. "
            "Evidence rescue is a hard exemption; monosublocus-overlap rescue is not."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-max-model-structural-support-score-scaled",
        type=float,
        default=argparse.SUPPRESS,
        help=(
            "Maximum model_structural_support_score_scaled for optional low-score singleton pruning. "
            "Explicit values override fixed pruning presets and are not accepted with auto/auto-conservative. Default: 0.25."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-max-protein-evidence-score-scaled",
        type=float,
        default=argparse.SUPPRESS,
        help=(
            "Maximum protein_evidence_score_scaled for optional low-score singleton pruning when protein evidence is supplied. "
            "Explicit values override fixed pruning presets and are not accepted with auto/auto-conservative. Default: 0.25."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-max-transcript-evidence-score-scaled",
        type=float,
        default=argparse.SUPPRESS,
        help=(
            "Maximum transcript_evidence_score_scaled for optional low-score singleton pruning when transcript evidence is supplied. "
            "Explicit values override fixed pruning presets and are not accepted with auto/auto-conservative. Default: 0.40."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-max-model-selection-score-scaled",
        type=float,
        default=argparse.SUPPRESS,
        help=(
            "Maximum model_selection_score_scaled for optional low-score singleton pruning. "
            "Explicit values override fixed pruning presets and are not accepted with auto/auto-conservative. Default: 0.25."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-use-transcript-cdna-context",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Allow calculated transcript cDNA-context concern to provide an alternative singleton pruning route when model support is low "
            "and protein support is also low when supplied. Transcript evidence adequacy must pass. Default: true."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-require-protein-evidence-adequacy",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Require run-wide protein evidence adequacy before optional low-score singleton pruning is allowed. "
            "Default: true."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-min-protein-overlap-model-fraction",
        type=float,
        default=0.75,
        help=(
            "Minimum fraction of scored models with overlapping protein evidence required by the low-score singleton pruning adequacy gate. "
            "Default: 0.75."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-min-protein-supported-model-fraction",
        type=float,
        default=0.50,
        help=(
            "Minimum fraction of scored models with partial or strong protein support required by the low-score singleton pruning adequacy gate. "
            "Default: 0.50."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-require-transcript-evidence-adequacy",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Require run-wide transcript evidence adequacy before transcript evidence can satisfy the low-score singleton pruning adequacy gate. "
            "Default: true."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-min-transcript-scoring-eligible-model-fraction",
        type=float,
        default=0.75,
        help=(
            "Minimum fraction of scored models with scoring-eligible transcript evidence required by the low-score singleton pruning adequacy gate. "
            "Default: 0.75."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-singleton-min-transcript-supported-model-fraction",
        type=float,
        default=0.50,
        help=(
            "Minimum fraction of scored models with partial or strong transcript support required by the low-score singleton pruning adequacy gate. "
            "Default: 0.50."
        ),
    )
    final_prune_group.add_argument(
        "--final-prune-low-score-nonsingletons",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Low-score pruning for selected non-singleton-sublocus models. It uses the established evidence-aware policy when supplied evidence is adequate, "
            "or a canonical model-confidence fallback when no supplied evidence class is adequate and at least three model sources were run. "
            "If omitted, this is enabled by --pruning-level medium/high/very-high/automatic modes, including runs without external evidence. "
            "The revised medium/high/very-high presets all use all-sublocus scope. If explicitly enabled, the default scope is all unless overridden. "
            "It is applied after singleton low-score pruning and before opposite-strand shadow pruning."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-scope",
        choices=("backfill", "all"),
        default=None,
        help=(
            "Candidate scope for low-score nonsingleton pruning. backfill considers only models introduced by backfill; "
            "all considers every selected nonsingleton model. The revised medium/high/very-high presets use all; "
            "explicit --final-prune-low-score-nonsingletons also defaults to all."
        ),
    )
    final_prune_group.add_argument(
        "--final-prune-sequence-integrity-nonsingletons",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Prune selected nonsingleton models with genome-derived sequence-integrity concerns. "
            "Preset-driven medium uses scope=backfill when --genome-fasta is supplied; high and very-high use scope=all. "
            "Explicit enablement defaults to scope=all. Backfill scope uses internal stop, missing start, missing stop, "
            "or all-splices-noncanonical concerns. All scope retains that backfill rule and additionally applies severe concerns to every selected nonsingleton. Strong configured protein/transcript rescue criteria protect models."
        ),
    )
    final_prune_group.add_argument(
        "--sequence-integrity-nonsingleton-scope",
        choices=("backfill", "all"),
        default=None,
        help=(
            "Candidate scope for sequence-integrity nonsingleton pruning. backfill applies the standard concern class only to models introduced by backfill; "
            "all retains that rule and additionally applies the severe concern class to every selected nonsingleton model. "
            "Medium uses backfill; high, very-high, and explicit enablement use all."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-max-model-structural-support-score-scaled",
        type=float,
        default=argparse.SUPPRESS,
        help=(
            "Maximum model_structural_support_score_scaled for evidence-aware low-score nonsingleton pruning. "
            "This does not change the canonical model-confidence fallback thresholds (0.20/0.25/0.35/0.45), which are fixed by pruning level and source depth. "
            "Explicit values override fixed evidence-aware presets and are not accepted with auto/auto-conservative. Default: 0.20."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-max-protein-evidence-score-scaled",
        type=float,
        default=argparse.SUPPRESS,
        help=(
            "Maximum protein_evidence_score_scaled for evidence-aware nonsingleton pruning when protein evidence is supplied; it also sets the local protein veto threshold in model-confidence fallback when supplied evidence is inadequate. "
            "Explicit values override fixed evidence thresholds and are not accepted with auto/auto-conservative. Default: 0.20."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-max-transcript-evidence-score-scaled",
        type=float,
        default=argparse.SUPPRESS,
        help=(
            "Maximum transcript_evidence_score_scaled for the conventional evidence-aware nonsingleton route when transcript evidence is supplied; it also sets the local transcript veto threshold for the conventional model-confidence fallback route when supplied evidence is inadequate. "
            "Explicit values override fixed evidence thresholds and are not accepted with auto/auto-conservative. Default: 0.40."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-max-model-selection-score-scaled",
        type=float,
        default=argparse.SUPPRESS,
        help=(
            "Maximum model_selection_score_scaled for evidence-aware low-score nonsingleton pruning. Model-confidence fallback does not use model selection score, so this option does not change its canonical thresholds or decisions. "
            "Explicit values override fixed evidence-aware presets and are not accepted with auto/auto-conservative. Default: 0.20."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-use-transcript-cdna-context",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Allow calculated transcript cDNA-context concern to provide an alternative nonsingleton pruning route within the configured scope. "
            "In evidence-aware mode, transcript adequacy must pass. In model-confidence fallback, an evaluated local concern can provide this route when no evidence class passes adequacy, without requiring global transcript adequacy. Default: true."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-require-protein-evidence-adequacy",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Require run-wide protein evidence adequacy before optional low-score nonsingleton pruning is allowed. "
            "Default: true."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-min-protein-overlap-model-fraction",
        type=float,
        default=0.75,
        help=(
            "Minimum fraction of scored models with overlapping protein evidence required by the low-score nonsingleton pruning adequacy gate. "
            "Default: 0.75."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-min-protein-supported-model-fraction",
        type=float,
        default=0.50,
        help=(
            "Minimum fraction of scored models with partial or strong protein support required by the low-score nonsingleton pruning adequacy gate. "
            "Default: 0.50."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-require-transcript-evidence-adequacy",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Require run-wide transcript evidence adequacy before transcript evidence can satisfy the low-score nonsingleton pruning adequacy gate. "
            "Default: true."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-min-transcript-scoring-eligible-model-fraction",
        type=float,
        default=0.75,
        help=(
            "Minimum fraction of scored models with scoring-eligible transcript evidence required by the low-score nonsingleton pruning adequacy gate. "
            "Default: 0.75."
        ),
    )
    final_prune_group.add_argument(
        "--low-score-nonsingleton-min-transcript-supported-model-fraction",
        type=float,
        default=0.50,
        help=(
            "Minimum fraction of scored models with partial or strong transcript support required by the low-score nonsingleton pruning adequacy gate. "
            "Default: 0.50."
        ),
    )

    # Monosublocus de-overlap and backfill are controlled by --stop-after.

    source_prior_group = parser.add_argument_group("source-prior weights for selection scoring")
    source_prior_group.add_argument(
        "--source-prior-mode",
        choices=["manual", "auto", "none"],
        default="auto",
        help=(
            "How to define source_prior_score, which is used as one component of provisional model ranking. "
            "'auto' derives effective source weights from exact coordinate agreement and scales them into the "
            "--auto-source-weight-min/--auto-source-weight-max range; with fewer than three input sources, auto weights are neutral at 1.0. "
            "'manual' uses --source-weight values with unspecified sources at 1.0; 'none' disables the source-prior component."
        ),
    )
    source_prior_group.add_argument(
        "--source-weight",
        action="append",
        default=None,
        metavar="SOURCE=VALUE",
        help=(
            "Manual source prior weight. Repeat for multiple sources, e.g. "
            "--source-weight Tiberius=1.0 --source-weight Augustus=0.7. "
            "Comma-separated entries are also accepted. Unspecified sources default to 1.0."
        ),
    )
    source_prior_group.add_argument(
        "--source-prior-scale",
        type=float,
        default=1.0,
        help=(
            "Multiplier used to calculate source_prior_score = effective_source_weight * source_prior_scale. "
            "source_prior_score is then included in model selection ranking through --selection-weight-source-prior."
        ),
    )
    source_prior_group.add_argument(
        "--auto-source-weight-min",
        type=float,
        default=0.5,
        help="Minimum automatic source weight when --source-prior-mode auto is used.",
    )
    source_prior_group.add_argument(
        "--auto-source-weight-max",
        type=float,
        default=1.0,
        help="Maximum automatic source weight when --source-prior-mode auto is used.",
    )

    model_validation = parser.add_argument_group("model-level validation before superlocus construction")
    model_validation.add_argument(
        "--model-validation-action",
        choices=["exclude", "stop", "warn"],
        default="exclude",
        help=(
            "Action for individual models exceeding model-level size thresholds. "
            "'exclude' removes them before superlocus construction and reports them; "
            "'stop' exits before superlocus construction; 'warn' keeps them."
        ),
    )
    model_validation.add_argument(
        "--max-model-cds-span-bp",
        type=int,
        default=2_000_000,
        help="Maximum allowed genomic span from first to last CDS base for one model.",
    )
    model_validation.add_argument(
        "--max-model-cds-intron-bp",
        type=int,
        default=1_000_000,
        help="Maximum allowed CDS intron length for one model.",
    )

    superlocus_validation = parser.add_argument_group("superlocus-level validation after clustering")
    superlocus_validation.add_argument(
        "--superlocus-validation-action",
        choices=["stop", "warn"],
        default="stop",
        help="Action if superloci exceed superlocus-level thresholds.",
    )
    superlocus_validation.add_argument(
        "--max-superlocus-span-bp",
        type=int,
        default=5_000_000,
        help="Maximum allowed genomic span for one superlocus.",
    )
    superlocus_validation.add_argument(
        "--max-superlocus-gene-count",
        type=int,
        default=250,
        help=(
            "Maximum allowed number of distinct input gene IDs in one superlocus. "
            "This is gene-level, so many transcript isoforms from one input gene do not inflate it."
        ),
    )
    superlocus_validation.add_argument(
        "--max-genes-per-source-in-superlocus",
        type=int,
        default=100,
        help="Maximum allowed number of distinct input gene IDs from any one source in one superlocus.",
    )

    sublocus_group = parser.add_argument_group("first-pass sublocus clustering")
    sublocus_group.add_argument(
        "--single-exon-min-reciprocal-overlap",
        type=float,
        default=0.5,
        help=(
            "Minimum reciprocal CDS-span overlap required to connect two single-exon models. "
            "Multi-exon models are connected by shared CDS splice-site coordinates. "
            "Single-exon and multi-exon models are not connected in this first pass."
        ),
    )


    structural_group = parser.add_argument_group("structural support metrics")
    structural_group.add_argument(
        "--terminal-cds-3prime-wobble-bp",
        type=int,
        default=3,
        help=(
            "Allow this many bp of strand-aware outer-boundary wobble when matching "
            "3' terminal CDS exons for structural support. Default: 3; set 0 for exact only."
        ),
    )
    structural_group.add_argument(
        "--single-exon-support-min-reciprocal-overlap",
        type=float,
        default=None,
        help=(
            "Minimum reciprocal CDS overlap used for single-exon structural support. "
            "If omitted, uses --single-exon-min-reciprocal-overlap."
        ),
    )

    penalty_group = parser.add_argument_group("intrinsic penalty weights")
    penalty_group.add_argument("--intrinsic-penalty-short-cds", type=float, default=1.0, help="Penalty added for a CDS shorter than --min-cds-length-bp. Default: 1.0.")
    penalty_group.add_argument("--intrinsic-penalty-cds-length-not-multiple-of-3", type=float, default=1.0, help="Penalty added when CDS length is not a multiple of 3. Default: 1.0.")
    penalty_group.add_argument("--intrinsic-penalty-short-terminal-cds-exon", type=float, default=0.25, help="Penalty added per short terminal CDS exon. Default: 0.25.")
    penalty_group.add_argument("--intrinsic-penalty-short-internal-cds-exon", type=float, default=0.75, help="Penalty added per short internal CDS exon. Default: 0.75.")
    penalty_group.add_argument("--intrinsic-penalty-short-cds-intron", type=float, default=1.0, help="Penalty added per short CDS intron. Default: 1.0.")
    penalty_group.add_argument("--intrinsic-penalty-long-cds-intron", type=float, default=0.75, help="Penalty added per long CDS intron. Default: 0.75.")
    penalty_group.add_argument("--intrinsic-penalty-contained-model", type=float, default=0.5, help="Penalty multiplier for likely contained models. The contribution is now fractional: contained_model_source_fraction * this value, and containment only applies when the model is missing supported terminal CDS structure. Default: 0.5.")
    penalty_group.add_argument("--intrinsic-penalty-unsupported-contained-cds-intron", type=float, default=0.25, help="Penalty multiplier for unsupported CDS introns fully contained within another model's CDS exon in the same sublocus. Contribution is fraction of query CDS introns affected times this value. Default: 0.25.")

    selection_group = parser.add_argument_group("model selection ranking")
    selection_group.add_argument("--selection-weight-structural-support", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["structural_support"], help="Positive weight for structural_support_score in model selection ranking.")
    selection_group.add_argument("--selection-weight-consensus-similarity", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["consensus_similarity"], help="Positive weight for consensus_similarity_score in model selection ranking.")
    selection_group.add_argument("--selection-weight-pairwise-similarity", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["pairwise_similarity"], help="Positive weight for pairwise_similarity_score in model selection ranking.")
    selection_group.add_argument("--selection-weight-completeness", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["completeness"], help="Positive weight for model_completeness_score in model selection ranking.")
    selection_group.add_argument("--selection-weight-source-prior", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["source_prior"], help="Positive weight for source_prior_score in model selection ranking.")
    selection_group.add_argument("--selection-weight-terminal-extension", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["terminal_extension"], help="Positive weight for terminal_extension_score in model selection ranking.")
    selection_group.add_argument("--selection-weight-exact-structure-support", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["exact_structure_support"], help="Positive weight for exact_structure_support_score in model selection ranking. This rewards exact CDS-exon-chain support scaled by the run-wide number of possible other sources.")
    selection_group.add_argument(
        "--selection-weight-structural-corroboration",
        type=float,
        default=MODEL_SELECTION_DEFAULT_WEIGHTS["structural_corroboration"],
        help=(
            "Positive weight for run-aware, source-collapsed whole-gene structural "
            "corroboration. The accepted release default is 0.5."
        ),
    )
    selection_group.add_argument(
        "--structural-corroboration-alpha",
        type=float,
        default=16.0,
        help=(
            "Positive fixed alpha for saturating run-depth normalization of structural corroboration. "
            "Default: 16; alpha 8 is supported as a comparator."
        ),
    )
    selection_group.add_argument(
        "--selection-weight-protein-evidence",
        type=float,
        default=MODEL_SELECTION_DEFAULT_WEIGHTS["protein_evidence"],
        help=(
            "Positive weight for protein_evidence_score in final model selection ranking. "
            "Default: 3.0 when --evidence-proteins is supplied; if no protein evidence files are supplied, "
            "the effective protein-evidence contribution is 0 and other active components are unchanged."
        ),
    )
    selection_group.add_argument(
        "--selection-weight-transcript-evidence",
        type=float,
        default=MODEL_SELECTION_DEFAULT_WEIGHTS["transcript_evidence"],
        help=(
            "Positive weight for the context-adjusted transcript evidence score in final model selection ranking. "
            "Default: 2.0 when --evidence-transcripts is supplied; if no transcript evidence files are supplied, "
            "the effective transcript-evidence contribution is 0 and other active components are unchanged."
        ),
    )
    selection_group.add_argument(
        "--selection-penalty-unsupported-short-cds-exon",
        type=float,
        default=MODEL_SELECTION_DEFAULT_WEIGHTS["unsupported_short_cds_exon_penalty"],
        help=(
            "Maximum model-only selection penalty multiplier for moderately short CDS exons that lack matching support "
            "from another annotation source. The normalized penalty combines 5-prime terminal, 3-prime terminal, and "
            "internal components. It does not contribute to intrinsic singleton pruning. Default: 0.50."
        ),
    )
    selection_group.add_argument(
        "--selection-penalty-structural-coherence",
        type=float,
        default=MODEL_SELECTION_DEFAULT_WEIGHTS["structural_coherence_penalty"],
        help=(
            "Model-only penalty multiplier for source-collapsed candidate structural-coherence "
            "breaks. The default maximum for a normal 2:0 region is 0.30; broader support can "
            "produce a bounded raw score up to 1.50. Default: 0.30."
        ),
    )
    selection_group.add_argument(
        "--unsupported-short-5prime-terminal-weight",
        type=float,
        default=0.25,
        help="Component weight for an unsupported moderately short biological 5-prime terminal CDS exon. Default: 0.25.",
    )
    selection_group.add_argument(
        "--unsupported-short-3prime-terminal-weight",
        type=float,
        default=0.25,
        help="Component weight for an unsupported moderately short biological 3-prime terminal CDS exon. Default: 0.25.",
    )
    selection_group.add_argument(
        "--unsupported-short-internal-weight",
        type=float,
        default=0.50,
        help=(
            "Component weight for unsupported moderately short internal CDS exons. Per-exon severities accumulate and "
            "are capped at 1 before this weight is applied. Default: 0.50."
        ),
    )
    selection_group.add_argument(
        "--selection-penalty-transcript-cdna-context",
        type=float,
        default=0.15,
        help=(
            "Maximum subtractive adjustment applied to transcript_evidence_score when a calculated CDS/min-cDNA fraction is below "
            "--transcript-cdna-context-neutral-fraction. Missing context is neutral. Default: 0.15."
        ),
    )
    selection_group.add_argument(
        "--transcript-cdna-context-neutral-fraction",
        type=float,
        default=0.30,
        help=(
            "CDS/min-cDNA fraction at or above which cDNA context has no scoring concern. "
            "Below this value a continuous concern score is calculated. Default: 0.30."
        ),
    )
    selection_group.add_argument(
        "--transcript-cdna-context-severe-fraction",
        type=float,
        default=0.15,
        help=(
            "CDS/min-cDNA fraction below which calculated context is a severe pruning concern without requiring an additional structural flag. Default: 0.15."
        ),
    )
    selection_group.add_argument(
        "--selection-penalty-sequence-integrity",
        type=float,
        default=MODEL_SELECTION_DEFAULT_WEIGHTS["sequence_integrity_penalty"],
        help=(
            "Model-only selection penalty multiplier for genome-derived sequence_integrity_penalty_score. "
            "Effective only when --genome-fasta is supplied. Default: 1.0."
        ),
    )
    selection_group.add_argument(
        "--no-sequence-integrity-scoring",
        action="store_true",
        help="Calculate and report genome sequence-integrity metrics but apply zero sequence penalty in model selection.",
    )
    selection_group.add_argument(
        "--sequence-penalty-internal-stop",
        type=float,
        default=1.0,
        help="Penalty per internal in-frame stop codon before the count cap. Default: 1.0.",
    )
    selection_group.add_argument(
        "--sequence-internal-stop-count-cap",
        type=int,
        default=2,
        help="Maximum internal-stop count contributing to the sequence-integrity penalty. Default: 2.",
    )
    selection_group.add_argument(
        "--sequence-penalty-phase-inconsistent",
        type=float,
        default=0.75,
        help="Penalty for an unambiguous inconsistency among declared CDS phases. Default: 0.75.",
    )
    selection_group.add_argument(
        "--sequence-penalty-missing-stop",
        type=float,
        default=SEQUENCE_INTEGRITY_DEFAULT_MISSING_STOP_PENALTY,
        help=(
            "Penalty when neither an included nor immediately adjacent in-frame terminal stop is present. "
            f"Default: {SEQUENCE_INTEGRITY_DEFAULT_MISSING_STOP_PENALTY:.2f}; "
            f"historical P0: {SEQUENCE_INTEGRITY_P0_MISSING_STOP_PENALTY:.2f}."
        ),
    )
    selection_group.add_argument(
        "--sequence-penalty-missing-start",
        type=float,
        default=SEQUENCE_INTEGRITY_DEFAULT_MISSING_START_PENALTY,
        help=(
            "Penalty when the translated CDS does not begin with ATG. Exact contig-boundary and ambiguous "
            f"cases are neutral. Default: {SEQUENCE_INTEGRITY_DEFAULT_MISSING_START_PENALTY:.2f}; "
            f"historical P0: {SEQUENCE_INTEGRITY_P0_MISSING_START_PENALTY:.2f}."
        ),
    )
    selection_group.add_argument(
        "--sequence-penalty-noncanonical-splice",
        type=float,
        default=1.00,
        help=(
            "Penalty multiplier for the fraction of CDS introns whose transcript-oriented motif is not "
            "GT-AG, GC-AG, or AT-AC. Motifs containing ambiguous bases are noncanonical. Default: 1.00."
        ),
    )
    selection_group.add_argument("--selection-penalty-intrinsic", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["intrinsic_penalty"], help="Penalty weight for intrinsic_penalty_score in model selection ranking.")
    selection_group.add_argument("--selection-penalty-partial-structure", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["partial_structure_penalty"], help="Penalty weight for partial_structure_penalty in model selection ranking.")
    selection_group.add_argument("--selection-penalty-internal-extra", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["internal_extra_penalty"], help="Penalty weight for internal_extra_fraction in model selection ranking.")
    selection_group.add_argument("--selection-penalty-terminal-extra", type=float, default=MODEL_SELECTION_DEFAULT_WEIGHTS["terminal_extra_penalty"], help="Penalty weight for terminal_extra_fraction in model selection ranking. Default is 0 because terminal extra CDS can represent useful terminal extension.")
    selection_group.add_argument("--model-selection-tie-epsilon", type=float, default=1e-6, help="Scores within this absolute difference are treated as ties for model_selection_score_rank.")
    selection_group.add_argument("--show-zero-weight-score-components", action="store_true", help="Include zero-weight terms in model_selection_score_components. By default, disabled terms are omitted to avoid misleading score breakdowns.")
    selection_group.add_argument("--gff3-include-model-selection-score-components", action="store_true", help="Include the full model_selection_score_components attribute on gene/mRNA features in the combined GFF3. This can make the GFF3 much larger.")
    selection_group.add_argument("--no-low-source-conflict-dampen-agreement-scores", action="store_false", dest="low_source_conflict_dampen_agreement_scores", default=True, help=("Disable the default low-source conflict dampener. By default, two-source/two-model conflicts with no exact structure support and pairwise_similarity_score < 0.80 have structural_support_score and consensus_similarity_score multiplied by 0.25 in the model selection score only."))

    metrics_group = parser.add_argument_group("coordinate-only intrinsic model metrics")
    metrics_group.add_argument(
        "--min-cds-length-bp",
        type=int,
        default=150,
        help="Flag models with CDS length below this value as short_cds.",
    )
    metrics_group.add_argument(
        "--min-terminal-cds-exon-bp",
        type=int,
        default=None,
        help=(
            "Flag models with a terminal CDS exon shorter than this value as "
            "short_terminal_cds_exon. Terminal CDS exons are the first/last CDS "
            "exons in genomic order; a single-exon CDS is terminal. Default: 3."
        ),
    )
    metrics_group.add_argument(
        "--min-internal-cds-exon-bp",
        type=int,
        default=None,
        help=(
            "Flag models with an internal CDS exon shorter than this value as "
            "short_internal_cds_exon. One- and two-exon CDS models have no internal "
            "CDS exons. Default: 6."
        ),
    )
    metrics_group.add_argument(
        "--unsupported-short-5prime-terminal-cds-exon-bp",
        type=int,
        default=15,
        help=(
            "Soft upper threshold for the unsupported-short-exon selection penalty on biological 5-prime terminal CDS "
            "exons. An exon is eligible only when its length is at least the active intrinsic terminal threshold, shorter "
            "than this value, and unsupported by another annotation source. Default: 15."
        ),
    )
    metrics_group.add_argument(
        "--unsupported-short-3prime-terminal-cds-exon-bp",
        type=int,
        default=15,
        help=(
            "Soft upper threshold for the unsupported-short-exon selection penalty on biological 3-prime terminal CDS "
            "exons. Existing strand-aware 3-prime terminal wobble is used when assessing support. Default: 15."
        ),
    )
    metrics_group.add_argument(
        "--unsupported-short-internal-cds-exon-bp",
        type=int,
        default=15,
        help=(
            "Soft upper threshold for the unsupported-short-exon selection penalty on internal CDS exons. An exon is "
            "eligible only when its length is at least the active intrinsic internal threshold, shorter than this value, "
            "and unsupported by another annotation source. Default: 15."
        ),
    )
    metrics_group.add_argument(
        "--min-cds-exon-bp",
        type=int,
        default=None,
        help=(
            "Deprecated convenience option: set both terminal and internal CDS exon "
            "minimums when the specific options are not supplied."
        ),
    )
    metrics_group.add_argument(
        "--min-cds-intron-bp",
        type=int,
        default=20,
        help="Flag models with any CDS intron shorter than this value as short_cds_intron.",
    )
    metrics_group.add_argument(
        "--long-cds-intron-bp",
        type=int,
        default=None,
        help=(
            "Flag models with any CDS intron longer than this value as long_cds_intron. "
            "If omitted, the cutoff is calculated from all valid input model introns using "
            "--long-cds-intron-percentile."
        ),
    )
    metrics_group.add_argument(
        "--long-cds-intron-percentile",
        type=float,
        default=99.999,
        help=(
            "Percentile of valid input model CDS intron lengths used as the long_cds_intron "
            "flag cutoff when --long-cds-intron-bp is not supplied. Default: 99.999."
        ),
    )

    refinement_group = parser.add_argument_group(
        "advanced post-selection CDS terminus refinement"
    )
    refinement_group.add_argument(
        "--cds-refinement-downstream-start-search-bp",
        type=nonnegative_int,
        default=30,
        help="Downstream missing-start fallback search in CDS bases. Default: 30.",
    )
    refinement_group.add_argument(
        "--cds-refinement-downstream-stop-search-bp",
        type=nonnegative_int,
        default=100,
        help="Downstream missing-stop search in spliced transcript bases. Default: 100.",
    )
    refinement_group.add_argument(
        "--cds-refinement-max-relative-extension",
        type=nonnegative_float,
        default=1.0,
        help="Maximum coding bases added divided by original selected CDS length. Default: 1.0.",
    )
    refinement_group.add_argument(
        "--cds-refinement-max-relative-truncation",
        type=unit_interval_float,
        default=0.33,
        help="Maximum coding bases removed divided by original selected CDS length. Default: 0.33.",
    )

    utr_group = parser.add_argument_group("advanced post-selection UTR augmentation")
    utr_group.add_argument(
        "--utr-min-complete-exon-bp", type=nonnegative_int, default=5,
        help="Minimum complete source/evidence exon length; 0 disables. Default: 5.",
    )
    utr_group.add_argument(
        "--utr-min-intron-bp", type=nonnegative_int, default=None, metavar="INT",
        help="Minimum new UTR intron; omitted inherits --min-cds-intron-bp, 0 disables.",
    )
    utr_group.add_argument(
        "--utr-max-intron-bp", type=nonnegative_int, default=None, metavar="INT",
        help="Maximum individual UTR intron; omitted inherits the resolved long-CDS-intron cutoff, 0 disables.",
    )
    utr_group.add_argument(
        "--utr-max-five-prime-exons", type=nonnegative_int, default=4,
        help="Maximum UTR-contributing exons on the biological 5-prime side; 0 disables. Default: 4.",
    )
    utr_group.add_argument(
        "--utr-max-three-prime-exons", type=nonnegative_int, default=2,
        help="Maximum UTR-contributing exons on the biological 3-prime side; 0 disables. Default: 2.",
    )
    utr_group.add_argument(
        "--utr-exonic-failsafe-bp", type=nonnegative_int, default=20_000,
        help="Absolute exonic-extension fail-safe; 0 disables. Default: 20000.",
    )
    for option, help_text in (
        ("five-prime-exonic", "5-prime total exonic-extension limit"),
        ("three-prime-exonic", "3-prime total exonic-extension limit"),
        ("five-prime-genomic", "5-prime genomic-extension limit"),
        ("three-prime-genomic", "3-prime genomic-extension limit"),
    ):
        utr_group.add_argument(
            f"--utr-{option}-override-bp", type=nonnegative_int, default=None,
            metavar="INT", help=f"Explicit {help_text}; omitted learns automatically, 0 disables.",
        )
    utr_group.add_argument(
        "--utr-member-depth-cap", type=nonnegative_int, default=3,
        help="Cap applied to log2 transcript member-depth support. Default: 3.",
    )
    utr_group.add_argument(
        "--utr-endpoint-window-min-bp", type=nonnegative_int, default=20,
        help="Minimum adaptive endpoint-window width. Default: 20.",
    )
    utr_group.add_argument(
        "--utr-endpoint-window-max-bp", type=nonnegative_int, default=250,
        help="Maximum adaptive endpoint-window width. Default: 250.",
    )
    utr_group.add_argument(
        "--utr-local-consensus-min-fraction", type=unit_interval_float, default=0.5,
        help="Local endpoint support must be strictly greater than this fraction. Default: 0.5.",
    )
    utr_group.add_argument(
        "--utr-limit-min-observations", type=nonnegative_int, default=1000,
        help="Minimum observations for side-specific extension-limit learning. Default: 1000.",
    )

    runtime = parser.add_argument_group("runtime")
    runtime.add_argument(
        "--strict",
        action="store_true",
        help="Raise errors on malformed records instead of warning and skipping.",
    )

    if not detailed_help:
        # Keep everyday --help focused. These options remain valid, but are shown only by --full-help.
        advanced_dests = {
            "combined_gff3", "scored_pruned_gff3", "superloci_tsv", "members_tsv", "models_tsv",
            "excluded_models_tsv", "other_features_gff3", "input_summary_tsv", "validation_failures_tsv",
            "subloci_tsv", "sublocus_members_tsv", "relationship_edges_tsv", "model_metrics_tsv",
            "sequence_integrity_summary_tsv", "sequence_integrity_model_metrics_tsv",
            "source_auto_weights_tsv", "model_structural_support_tsv", "structural_coherence_regions_tsv", "model_selection_scores_tsv", "model_selection_score_summary_tsv", "structural_features_tsv",
            "sublocus_winners_gff3", "sublocus_winners_tsv", "monosubloci_tsv",
            "monosublocus_members_tsv", "monosublocus_edges_tsv", "monosublocus_review_gff3",
            "selected_models_gff3", "selected_models_tsv", "suppressed_monosublocus_models_tsv",
            "backfill_candidates_tsv", "backfilled_models_tsv", "backfill_rejected_models_tsv",
            "r1_backfill_replacements_tsv", "r1_backfill_summary_tsv",
            "final_selected_models_tsv", "final_selected_models_gff3",
            "final_pruned_selected_models_tsv", "final_pruned_selected_models_gff3",
            "final_pruned_selected_removed_tsv", "final_shadow_gene_pairs_tsv",
            "standardised_primary_models_gff3", "final_annotation_gff3",
            "run_configuration_json",
            "terminal_stop_canonicalisation_summary_tsv", "terminal_stop_canonicalisation_actions_tsv",
            "terminal_stop_canonicalisation_audit_tsv",
            "cds_5prime_boundary_canonicalisation_summary_tsv", "cds_5prime_boundary_canonicalisation_actions_tsv",
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
            "utr_min_complete_exon_bp", "utr_min_intron_bp", "utr_max_intron_bp",
            "utr_max_five_prime_exons", "utr_max_three_prime_exons",
            "utr_exonic_failsafe_bp", "utr_five_prime_exonic_override_bp",
            "utr_three_prime_exonic_override_bp", "utr_five_prime_genomic_override_bp",
            "utr_three_prime_genomic_override_bp", "utr_member_depth_cap",
            "utr_endpoint_window_min_bp", "utr_endpoint_window_max_bp",
            "utr_local_consensus_min_fraction", "utr_limit_min_observations",
            "cds_refinement_downstream_start_search_bp",
            "cds_refinement_downstream_stop_search_bp",
            "cds_refinement_max_relative_extension",
            "cds_refinement_max_relative_truncation",
            "final_prune_singleton_subloci_level", "final_prune_singleton_subloci_source",
            "final_prune_opposite_strand_shadows", "shadow_prune_source",
            "shadow_prune_require_loser_single_exon", "shadow_prune_min_overlap_bp",
            "shadow_prune_min_loser_cds_overlap_fraction", "shadow_prune_min_score_margin_scaled",
            "singleton_prune_rescue_monosublocus_overlap", "singleton_rescue_min_monosublocus_member_count",
            "monosublocus_overlap_rescue_evidence_policy",
            "singleton_rescue_min_cds_overlap_bp", "singleton_rescue_min_self_overlap_fraction",
            "singleton_rescue_min_reciprocal_overlap_fraction",
            "final_prune_rescue_protein_supported_singletons", "singleton_rescue_min_protein_evidence_score_scaled",
            "singleton_rescue_min_protein_exact_chain_source_count",
            "final_prune_rescue_transcript_supported_singletons", "singleton_rescue_min_transcript_evidence_score_scaled",
            "singleton_rescue_min_transcript_supporting_source_count", "singleton_rescue_min_transcript_scoring_eligible_count",
            "singleton_rescue_min_transcript_cdna_context_fraction",
            "transcript_singleton_rescue_require_stranded_support", "transcript_singleton_rescue_block_on_cdna_context_concern",
            "allow_internal_stop_evidence_rescue",
            "final_prune_low_score_singletons", "low_score_singleton_max_model_structural_support_score_scaled",
            "low_score_singleton_max_protein_evidence_score_scaled", "low_score_singleton_max_transcript_evidence_score_scaled", "low_score_singleton_max_model_selection_score_scaled", "low_score_singleton_use_transcript_cdna_context",
            "low_score_singleton_require_protein_evidence_adequacy", "low_score_singleton_min_protein_overlap_model_fraction",
            "low_score_singleton_min_protein_supported_model_fraction", "low_score_singleton_require_transcript_evidence_adequacy",
            "low_score_singleton_min_transcript_scoring_eligible_model_fraction", "low_score_singleton_min_transcript_supported_model_fraction",
            "low_score_nonsingleton_scope", "sequence_integrity_nonsingleton_scope", "low_score_nonsingleton_max_model_structural_support_score_scaled",
            "low_score_nonsingleton_max_protein_evidence_score_scaled", "low_score_nonsingleton_max_transcript_evidence_score_scaled", "low_score_nonsingleton_max_model_selection_score_scaled", "low_score_nonsingleton_use_transcript_cdna_context",
            "low_score_nonsingleton_require_protein_evidence_adequacy", "low_score_nonsingleton_min_protein_overlap_model_fraction",
            "low_score_nonsingleton_min_protein_supported_model_fraction", "low_score_nonsingleton_require_transcript_evidence_adequacy",
            "low_score_nonsingleton_min_transcript_scoring_eligible_model_fraction", "low_score_nonsingleton_min_transcript_supported_model_fraction",
            "model_validation_action", "max_model_cds_span_bp", "max_model_cds_intron_bp",
            "superlocus_validation_action", "max_superlocus_span_bp", "max_superlocus_gene_count",
            "max_genes_per_source_in_superlocus", "single_exon_min_reciprocal_overlap",
            "terminal_cds_3prime_wobble_bp", "single_exon_support_min_reciprocal_overlap",
            "intrinsic_penalty_short_cds", "intrinsic_penalty_cds_length_not_multiple_of_3",
            "intrinsic_penalty_short_terminal_cds_exon", "intrinsic_penalty_short_internal_cds_exon",
            "intrinsic_penalty_short_cds_intron", "intrinsic_penalty_long_cds_intron",
            "intrinsic_penalty_contained_model", "intrinsic_penalty_unsupported_contained_cds_intron",
            "selection_weight_structural_support", "selection_weight_consensus_similarity",
            "selection_weight_pairwise_similarity", "selection_weight_completeness",
            "selection_weight_source_prior", "selection_weight_terminal_extension",
            "selection_weight_exact_structure_support", "selection_weight_structural_corroboration",
            "structural_corroboration_alpha",
            "selection_penalty_unsupported_short_cds_exon", "selection_penalty_structural_coherence", "unsupported_short_5prime_terminal_weight", "unsupported_short_3prime_terminal_weight", "unsupported_short_internal_weight",
            "unsupported_short_5prime_terminal_cds_exon_bp", "unsupported_short_3prime_terminal_cds_exon_bp", "unsupported_short_internal_cds_exon_bp",
            "selection_penalty_transcript_cdna_context", "transcript_cdna_context_neutral_fraction", "transcript_cdna_context_severe_fraction", "selection_penalty_intrinsic",
            "selection_penalty_sequence_integrity",
            "sequence_penalty_internal_stop", "sequence_internal_stop_count_cap",
            "sequence_penalty_phase_inconsistent", "sequence_penalty_missing_stop",
            "sequence_penalty_missing_start", "sequence_penalty_noncanonical_splice",
            "selection_penalty_partial_structure", "selection_penalty_internal_extra",
            "selection_penalty_terminal_extra", "model_selection_tie_epsilon",
            "show_zero_weight_score_components", "gff3_include_model_selection_score_components",
            "low_source_conflict_dampen_agreement_scores",
            "strict_r1_backfill_replacement", "strict_r1_allow_multiexon_blockers",
            "structural_corroboration_tsv",
            "protein_evidence_summary_tsv", "protein_evidence_raw_to_collapsed_tsv",
            "protein_evidence_collapsed_alignments_tsv", "protein_evidence_superloci_tsv",
            "protein_evidence_subloci_tsv", "protein_evidence_sublocus_members_tsv",
            "protein_evidence_model_support_tsv",
            "transcript_evidence_summary_tsv", "transcript_evidence_raw_to_collapsed_tsv",
            "transcript_evidence_collapsed_transcripts_tsv", "transcript_evidence_superloci_tsv",
            "transcript_evidence_subloci_tsv", "transcript_evidence_sublocus_members_tsv", "transcript_evidence_model_support_tsv",
        }
        for action in parser._actions:
            if action.dest in advanced_dests:
                action.help = argparse.SUPPRESS
        for action in parser._actions:
            if action.dest not in BASIC_HELP_DESTINATIONS:
                action.help = argparse.SUPPRESS
    elif not expert_help:
        expert_dests = expert_only_destinations(parser)
        for action in parser._actions:
            if action.dest in expert_dests:
                action.help = argparse.SUPPRESS
    return parser


def expand_name_list(values: Optional[list[str]]) -> set[str]:
    names: set[str] = set()
    if values:
        for item in values:
            for part in str(item).split(","):
                part = part.strip()
                if part:
                    names.add(part)
    return names


def output_keys_for_profile(profile: str) -> set[str]:
    minimal = {
        "run_configuration_json",
        "input_summary_tsv",
        "validation_failures_tsv",
        "source_auto_weights_tsv",
        "excluded_models_tsv",
        "model_selection_score_summary_tsv",
        "protein_evidence_summary_tsv",
        "transcript_evidence_summary_tsv",
        "sequence_integrity_summary_tsv",
        "final_annotation_gff3",
        "terminal_stop_canonicalisation_summary_tsv",
        "cds_5prime_boundary_canonicalisation_summary_tsv",
        "cds_phase_recovery_summary_tsv",
        "definite_internal_stop_exclusion_summary_tsv",
        "extreme_short_cds_exclusion_summary_tsv",
        "cds_terminus_refinement_summary_tsv",
        "utr_augmented_final_annotation_gff3",
        "utr_junction_validation_summary_tsv",
        "utr_augmentation_summary_tsv",
        "pruning_scenario_summary_tsv",
        "auto_pruning_diagnostics_tsv",
        "alternative_splice_run_summary_tsv",
    }
    standard = minimal | {
        "standardised_primary_models_gff3",
        "selected_models_tsv",
        "suppressed_monosublocus_models_tsv",
        "backfilled_models_tsv",
        "final_pruned_selected_models_tsv",
        "final_pruned_selected_removed_tsv",
        "final_shadow_gene_pairs_tsv",
        "terminal_stop_canonicalisation_actions_tsv",
        "cds_5prime_boundary_canonicalisation_actions_tsv",
        "cds_phase_recovery_actions_tsv",
        "definite_internal_stop_exclusions_tsv",
        "extreme_short_cds_exclusions_tsv",
        "cds_terminus_refinement_audit_tsv",
        "utr_augmentation_audit_tsv",
    }
    if profile == "minimal":
        return set(minimal)
    if profile == "standard":
        return set(standard)
    if profile in {"full", "debug"}:
        # Feature-gated R1 reports are enabled automatically only when R1 is on,
        # preserving disabled-mode output membership.
        return set(ALL_OUTPUT_KEYS) - {
            "r1_backfill_replacements_tsv", "r1_backfill_summary_tsv"
        }
    return set(standard)


def resolve_enabled_outputs(profile: str, write_values: Optional[list[str]], skip_values: Optional[list[str]]) -> set[str]:
    enabled = output_keys_for_profile(profile)
    forced = expand_name_list(write_values)
    skipped = expand_name_list(skip_values)
    unknown = sorted((forced | skipped) - ALL_OUTPUT_KEYS - {"all"})
    if unknown:
        raise SystemExit(
            "Unknown output name(s): "
            + ", ".join(unknown)
            + ". Known names: "
            + ", ".join(sorted(ALL_OUTPUT_KEYS))
        )
    if "all" in forced:
        enabled = set(ALL_OUTPUT_KEYS)
    else:
        enabled |= forced
    if "all" in skipped:
        enabled = set()
    else:
        enabled -= skipped
    return enabled


def with_default(path: Optional[Path], prefix: Path, suffix: str) -> Path:
    return path if path is not None else Path(str(prefix) + suffix)


def parse_manual_source_weights(weight_specs: Optional[list[str]], source_names: list[str]) -> dict[str, float]:
    """Parse repeated SOURCE=VALUE manual source-weight specifications."""
    weights = {source: 1.0 for source in source_names}
    if not weight_specs:
        return weights
    known_sources = set(source_names)
    expanded_specs: list[str] = []
    for spec in weight_specs:
        expanded_specs.extend(part.strip() for part in spec.split(",") if part.strip())
    for spec in expanded_specs:
        if "=" not in spec:
            raise SystemExit(f"Invalid --source-weight value '{spec}'. Expected SOURCE=VALUE.")
        source, value_text = spec.split("=", 1)
        source = source.strip()
        value_text = value_text.strip()
        if source not in known_sources:
            raise SystemExit(
                f"Unknown source in --source-weight: '{source}'. Known sources: {', '.join(source_names)}"
            )
        try:
            value = float(value_text)
        except ValueError as exc:
            raise SystemExit(f"Invalid numeric weight in --source-weight '{spec}'.") from exc
        if value < 0:
            raise SystemExit(f"Invalid --source-weight '{spec}': weights must be >= 0.")
        weights[source] = value
    return weights


def parse_source_filter_values(values: Optional[list[str]], known_sources: list[str], option_name: str) -> set[str]:
    sources: set[str] = set()
    if values:
        for item in values:
            for part in str(item).split(","):
                part = part.strip()
                if part:
                    sources.add(part)
    unknown = sorted(sources - set(known_sources))
    if unknown:
        raise SystemExit(
            f"Unknown {option_name} value(s): "
            + ", ".join(unknown)
            + ". Known sources: "
            + ", ".join(sorted(known_sources))
        )
    return sources


def count_pruned_single_model_subloci(subloci: list[Sublocus], prune_sources: set[str]) -> tuple[int, int, Counter[str]]:
    sublocus_count = 0
    model_count = 0
    source_counts: Counter[str] = Counter()
    for sub in subloci:
        if len(sub.models) != 1:
            continue
        source = sub.models[0].source_name
        if prune_sources and source not in prune_sources:
            continue
        sublocus_count += 1
        model_count += 1
        source_counts[source] += 1
    return sublocus_count, model_count, source_counts


def validate_source_prior_args(args: argparse.Namespace) -> None:
    if args.auto_source_weight_min < 0 or args.auto_source_weight_max < 0:
        raise SystemExit("--auto-source-weight-min and --auto-source-weight-max must be >= 0.")
    if args.auto_source_weight_min > args.auto_source_weight_max:
        raise SystemExit("--auto-source-weight-min cannot be greater than --auto-source-weight-max.")
    if args.source_prior_scale < 0:
        raise SystemExit("--source-prior-scale must be >= 0.")


def apply_pruning_level_presets(args: argparse.Namespace) -> None:
    """Resolve final-pruning preset shortcuts into effective per-filter options.

    Explicit advanced options override fixed presets. We use None or suppressed
    defaults on advanced options controlled by --pruning-level so argparse can
    distinguish "user did not say" from an explicit override.
    """
    threshold_option_names = (
        "low_score_singleton_max_model_structural_support_score_scaled",
        "low_score_singleton_max_protein_evidence_score_scaled",
        "low_score_singleton_max_transcript_evidence_score_scaled",
        "low_score_singleton_max_model_selection_score_scaled",
        "low_score_nonsingleton_max_model_structural_support_score_scaled",
        "low_score_nonsingleton_max_protein_evidence_score_scaled",
        "low_score_nonsingleton_max_transcript_evidence_score_scaled",
        "low_score_nonsingleton_max_model_selection_score_scaled",
    )
    explicit_threshold_overrides = tuple(
        name for name in threshold_option_names if getattr(args, name, None) is not None
    )
    if args.pruning_level in {"auto", "auto-conservative"} and explicit_threshold_overrides:
        cli_names = ", ".join("--" + name.replace("_", "-") for name in explicit_threshold_overrides)
        raise SystemExit(
            "Per-threshold low-score overrides cannot be combined with "
            f"--pruning-level {args.pruning_level}, because automatic modes simulate canonical "
            "medium/high/very-high scenarios. Choose a fixed --pruning-level to use: "
            + cli_names
        )

    preset = PRUNING_LEVEL_PRESETS.get(args.pruning_level, PRUNING_LEVEL_PRESETS["low"])
    fixed_level = "medium" if args.pruning_level in {"auto", "auto-conservative"} else args.pruning_level
    genome_fasta_supplied = bool(getattr(args, "genome_fasta", None))
    automatic_sequence_pruning_active = bool(
        genome_fasta_supplied and not getattr(args, "no_sequence_integrity_scoring", False)
    )
    if args.final_prune_singleton_subloci_level is None:
        preset_singleton_level = str(preset["singleton_level"])
        if preset_singleton_level == "intrinsic_or_sequence_integrity" and not automatic_sequence_pruning_active:
            preset_singleton_level = "intrinsic"
        args.final_prune_singleton_subloci_level = preset_singleton_level
    if args.final_prune_opposite_strand_shadows is None:
        args.final_prune_opposite_strand_shadows = bool(preset["shadow_enabled"])
    if args.shadow_prune_require_loser_single_exon is None:
        args.shadow_prune_require_loser_single_exon = bool(preset["shadow_require_loser_single_exon"])
    if args.singleton_prune_rescue_monosublocus_overlap is None:
        args.singleton_prune_rescue_monosublocus_overlap = bool(preset["singleton_rescue_monosublocus_overlap"])
    protein_evidence_supplied = bool(getattr(args, "evidence_proteins", None))
    transcript_evidence_supplied = bool(getattr(args, "evidence_transcripts", None))
    external_evidence_supplied = protein_evidence_supplied or transcript_evidence_supplied
    if args.final_prune_low_score_singletons is None:
        args.final_prune_low_score_singletons = bool(
            preset.get("low_score_singletons_when_protein_evidence", False) and external_evidence_supplied
        )

    explicit_low_score_nonsingletons = getattr(args, "final_prune_low_score_nonsingletons", None) is not None
    if getattr(args, "final_prune_low_score_nonsingletons", None) is None:
        args.final_prune_low_score_nonsingletons = bool(
            preset.get("low_score_nonsingletons", False)
        )
    if getattr(args, "low_score_nonsingleton_scope", None) is None:
        if args.final_prune_low_score_nonsingletons and not explicit_low_score_nonsingletons:
            args.low_score_nonsingleton_scope = str(
                preset.get("low_score_nonsingleton_scope", "backfill")
            )
        else:
            args.low_score_nonsingleton_scope = "all"
    explicit_sequence_integrity_nonsingletons = (
        getattr(args, "final_prune_sequence_integrity_nonsingletons", None) is not None
    )
    if getattr(args, "final_prune_sequence_integrity_nonsingletons", None) is None:
        args.final_prune_sequence_integrity_nonsingletons = bool(
            preset.get("sequence_integrity_nonsingletons_when_genome", False)
            and automatic_sequence_pruning_active
        )
    if getattr(args, "sequence_integrity_nonsingleton_scope", None) is None:
        if (
            args.final_prune_sequence_integrity_nonsingletons
            and not explicit_sequence_integrity_nonsingletons
        ):
            args.sequence_integrity_nonsingleton_scope = str(
                preset.get("sequence_integrity_nonsingleton_scope", "backfill")
            )
        else:
            args.sequence_integrity_nonsingleton_scope = "all"

    scenario = PRUNING_SCENARIO_THRESHOLDS.get(fixed_level)
    singleton = scenario["singleton"] if scenario is not None else (0.25, 0.25, 0.40, 0.25)
    nonsingleton = scenario["nonsingleton"] if scenario is not None else (0.20, 0.20, 0.40, 0.20)
    threshold_defaults = dict(zip(threshold_option_names, (*singleton, *nonsingleton)))
    for name, value in threshold_defaults.items():
        if getattr(args, name, None) is None:
            setattr(args, name, value)
    if scenario is not None:
        if not explicit_low_score_nonsingletons:
            args.low_score_nonsingleton_scope = scenario["low_score_nonsingleton_scope"]
        if not explicit_sequence_integrity_nonsingletons:
            args.sequence_integrity_nonsingleton_scope = scenario[
                "sequence_integrity_nonsingleton_scope"
            ]
    if not hasattr(args, "auto_pruning_max_level_explicit"):
        args.auto_pruning_max_level_explicit = (
            getattr(args, "auto_pruning_max_level", None) is not None
        )
    if getattr(args, "auto_pruning_max_level", None) is None:
        args.auto_pruning_max_level = (
            "high" if args.pruning_level == "auto-conservative" else "very-high"
        )
