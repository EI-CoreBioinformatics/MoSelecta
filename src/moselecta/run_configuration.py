"""Deterministic serialization of a resolved MoSelecta command configuration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

from .constants import PROGRAM_NAME, PROGRAM_VERSION


def _json_value(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {
            str(key): _json_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return [_json_value(item) for item in sorted(value, key=str)]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _option_subset(
    resolved: dict[str, Any],
    *,
    prefixes: Iterable[str] = (),
    names: Iterable[str] = (),
) -> dict[str, Any]:
    prefix_tuple = tuple(prefixes)
    name_set = set(names)
    return {
        key: _json_value(value)
        for key, value in sorted(resolved.items())
        if key in name_set or key.startswith(prefix_tuple)
    }


def build_run_configuration(
    args: argparse.Namespace,
    *,
    raw_arguments: list[str],
    enabled_output_keys: set[str],
    effective_pruning_level: str | None,
) -> dict[str, Any]:
    """Return a stable, JSON-native record of requested and resolved settings."""
    resolved = vars(args)
    requested_pruning_level = str(resolved.get("pruning_level", "none"))
    return {
        "alternative_splicing": _option_subset(
            resolved,
            prefixes=("alternative_splice_",),
            names=("no_alternative_splice",),
        ),
        "canonicalisation": _option_subset(
            resolved,
            names=("terminal_cds_3prime_wobble_bp",),
        ),
        "evidence": _option_subset(
            resolved,
            prefixes=("protein_evidence_", "transcript_evidence_", "singleton_rescue_"),
            names=(
                "evidence_proteins", "evidence_protein_source_name",
                "evidence_transcripts", "evidence_transcript_source_name",
            ),
        ),
        "definite_internal_stop_exclusion": {
            "requested_enabled": bool(
                resolved.get("definite_internal_stop_exclusion", True)
            ),
            "effective_enabled": bool(
                resolved.get("definite_internal_stop_exclusion", True)
                and resolved.get("genome_fasta") is not None
            ),
            "execution_status": (
                "performed"
                if (
                    resolved.get("definite_internal_stop_exclusion", True)
                    and resolved.get("genome_fasta") is not None
                )
                else "disabled_by_user"
                if not resolved.get("definite_internal_stop_exclusion", True)
                else "prerequisite_unavailable:no_genome_fasta"
            ),
            "resolved_frame_policy": "exclude_if_active_frame_has_literal_internal_stop",
            "unresolved_frame_policy": "exclude_if_all_three_frames_have_literal_internal_stop",
            "evidence_rescue": False,
        },
        "extreme_short_cds_exclusion": {
            "default_enabled": True,
            "effective_enabled": True,
            "minimum_coding_length_bp": int(
                resolved.get("minimum_coding_length_bp", 30)
            ),
            "included_terminal_stop_subtraction_bp": 3,
            "included_terminal_stop_detection": (
                "genome_sequence_available"
                if resolved.get("genome_fasta") is not None
                else "unavailable_without_genome_fasta"
            ),
            "no_genome_policy": "subtract_nothing_unless_included_stop_established",
            "evidence_rescue": False,
        },
        "inputs": {
            "genome_fasta": _json_value(resolved.get("genome_fasta")),
            "model_files": _json_value(resolved.get("inputs")),
            "model_source_names": _json_value(resolved.get("source_name")),
            "protein_evidence_files": _json_value(resolved.get("evidence_proteins")),
            "protein_evidence_source_names": _json_value(resolved.get("evidence_protein_source_name")),
            "transcript_evidence_files": _json_value(resolved.get("evidence_transcripts")),
            "transcript_evidence_source_names": _json_value(resolved.get("evidence_transcript_source_name")),
            "requested_arguments": list(raw_arguments),
        },
        "moselecta": {
            "name": PROGRAM_NAME,
            "version": PROGRAM_VERSION,
            "run_configuration_schema": 1,
        },
        "output": {
            "enabled_output_keys": sorted(enabled_output_keys),
            "out_prefix": _json_value(resolved.get("out_prefix")),
            "profile": resolved.get("output_profile"),
            "skip_output": _json_value(resolved.get("skip_output")),
            "stop_after": resolved.get("stop_after"),
            "write_output": _json_value(resolved.get("write_output")),
        },
        "phase": {
            "genetic_code": resolved.get("genetic_code"),
            "recovery_enabled_when_genome_supplied": resolved.get("genome_fasta") is not None,
        },
        "pruning": {
            "effective_level": effective_pruning_level,
            "requested_level": requested_pruning_level,
            "settings": _option_subset(
                resolved,
                prefixes=(
                    "auto_pruning_", "final_prune_", "low_score_",
                    "sequence_integrity_nonsingleton_", "shadow_prune_",
                    "singleton_prune_", "singleton_rescue_",
                ),
                names=("allow_internal_stop_evidence_rescue",),
            ),
        },
        "refinement": _option_subset(
            resolved,
            prefixes=("cds_refinement_",),
            names=("refine_cds_termini",),
        ),
        "resolved_options": {
            key: _json_value(value) for key, value in sorted(resolved.items())
        },
        "scoring": _option_subset(
            resolved,
            prefixes=(
                "selection_", "sequence_penalty_", "intrinsic_penalty_",
                "unsupported_short_", "transcript_cdna_context_",
            ),
            names=(
                "min_cds_length_bp", "min_terminal_cds_exon_bp",
                "min_internal_cds_exon_bp", "min_cds_intron_bp",
                "long_cds_intron_bp", "long_cds_intron_percentile",
                "source_prior_mode", "source_prior_scale", "source_weight",
            ),
        ),
        "selection": _option_subset(
            resolved,
            prefixes=("monosublocus_",),
            names=(
                "run_monosublocus_selection", "run_backfill_selection",
                "single_exon_min_reciprocal_overlap",
                "strict_r1_backfill_replacement",
                "strict_r1_backfill_replacement_configured",
                "strict_r1_backfill_replacement_explicit",
                "strict_r1_allow_multiexon_blockers",
                "strict_r1_allow_multiexon_blockers_configured",
                "strict_r1_allow_multiexon_blockers_explicit",
            ),
        ),
        "utr": _option_subset(
            resolved,
            prefixes=("utr_",),
            names=("add_utrs",),
        ),
        "validation": _option_subset(
            resolved,
            prefixes=("max_model_", "max_superlocus_", "max_genes_per_source_"),
            names=("model_validation_action", "superlocus_validation_action", "strict"),
        ),
    }


def write_run_configuration_json(configuration: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(configuration, indent=2, sort_keys=True) + "\n")
