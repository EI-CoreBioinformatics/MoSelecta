"""Report-only alternative-splice biological fact analysis.

The analysis runs after the retained primary set is frozen.  It classifies the
valid/scored models that were not retained, but never mutates model objects,
selection results, pruning results, identifiers, or established output files.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from .alternative_splice_selection import (
    STAGE2B_CANDIDATE_COLUMNS,
    STAGE2B_GENE_COLUMNS,
    STAGE2B_RUN_COLUMNS,
)
from .intervals import IntervalOverlapIndex, merge_intervals
from .models import CodingModel
from .protein_evidence import (
    ProteinEvidenceCollapsedAlignment,
    chains_match_with_3prime_wobble,
)
from .scoring import value_as_float
from .selection import (
    cds_exon_overlap_metrics,
    model_selection_row_sort_key,
    shared_cds_splice_site_count,
)
from .structural import (
    cds_exon_chain_containment_relation,
    cds_exon_chains_match,
    cds_frame_concordance,
    cds_intron_chains_match,
)
from .transcript_evidence import (
    TranscriptEvidenceSuperlocus,
    associated_transcripts_for_model,
    evaluate_coherent_exact_transcript_support,
)
from .tsvio import write_dict_tsv


@dataclass(frozen=True, slots=True)
class AlternativeSpliceReportPolicy:
    """Fixed historical settings and existing context semantics."""

    terminal_cds_3prime_wobble_bp: int
    protein_terminal_cds_3prime_wobble_bp: int
    transcript_min_cdna_context_fraction: float
    transcript_require_stranded_support: bool
    transcript_block_on_cdna_context_concern: bool
    transcript_cdna_context_neutral_fraction: float
    transcript_cdna_context_severe_fraction: float
    proposed_max_alternatives_per_gene: int = 50


ALTERNATIVE_SPLICE_CANDIDATE_COLUMNS = [
    "alternative_candidate_assessment_id",
    "candidate_source_name", "candidate_input_file", "candidate_input_gene_id",
    "candidate_input_transcript_id", "candidate_output_gene_id", "candidate_output_transcript_id",
    "chrom", "strand", "candidate_superlocus_id", "candidate_sublocus_id",
    "candidate_cds_start", "candidate_cds_end", "candidate_cds_length",
    "candidate_cds_exon_count", "candidate_cds_intron_count",
    "primary_source_name", "primary_input_gene_id", "primary_input_transcript_id",
    "primary_output_gene_id", "primary_output_transcript_id", "primary_superlocus_id",
    "primary_sublocus_id", "primary_cds_start", "primary_cds_end", "primary_cds_length",
    "primary_cds_exon_count", "primary_cds_intron_count",
    "structurally_associated_retained_gene_count", "structurally_associated_primary_ids",
    "structurally_associated_primary_sources", "ownership_status", "structural_association_route",
    "structural_candidate", "shared_cds_splice_site_count", "cds_overlap_bp",
    "primary_cds_overlap_fraction", "candidate_cds_overlap_fraction", "reciprocal_cds_overlap",
    "frame_concordance_status", "primary_phase_status", "candidate_phase_status",
    "primary_phase_mismatch_count", "candidate_phase_mismatch_count", "frame_overlap_bp",
    "frame_comparable_overlap_bp", "same_frame_overlap_bp", "same_frame_comparable_fraction",
    "same_frame_overlap_fraction", "same_frame_primary_fraction", "same_frame_candidate_fraction",
    "exact_cds_equivalent_to_primary", "utr_only_difference_from_primary",
    "same_cds_intron_chain_as_primary", "same_intron_chain_terminal_variant",
    "candidate_cds_chain_contained_in_primary", "primary_cds_chain_contained_in_candidate",
    "bidirectional_cds_chain_containment", "genuinely_splice_structure_distinct",
    "exact_equivalent_candidate_ids", "same_intron_chain_candidate_ids",
    "candidate_contained_in_other_candidate_ids", "other_candidate_contained_in_candidate_ids",
    "alternative_candidate_order_rank", "historical_model_selection_order_rank",
    "exact_model_other_source_count", "exact_model_other_sources", "exact_model_supporting_model_count",
    "exact_model_supporting_model_ids", "exact_structural_all_source_count",
    "exact_structural_all_sources", "exact_multisource_model_support",
    "alternative_splice_exact_model_source_count",
    "alternative_splice_exact_model_source_fraction",
    "mosaic_cds_introns_status", "mosaic_cds_splice_sites_status",
    "mosaic_internal_cds_exons_status", "mosaic_terminal_cds_exons_status",
    "mosaic_cds_base_consensus_status", "full_structural_mosaic",
    "alternative_splice_mosaic_min_component_source_count",
    "alternative_splice_mosaic_min_component_source_fraction",
    "protein_exact_cds_chain_support", "protein_exact_supporting_source_count",
    "protein_exact_supporting_sources", "protein_exact_supporting_group_count",
    "protein_exact_supporting_group_ids", "protein_effective_terminal_wobble_bp",
    "protein_strong_partial", "protein_evidence_score", "protein_evidence_score_scaled",
    "protein_best_alignment_score", "protein_mean_best_source_score",
    "protein_supported_cds_base_fraction", "protein_cds_intron_f1",
    "protein_cds_splice_site_f1", "protein_internal_cds_exon_f1",
    "protein_terminal_cds_exon_f1", "protein_evidence_status", "protein_evidence_flags",
    "exact_coherent_transcript_support", "exact_coherent_transcript_support_status",
    "exact_transcript_chain_group_count", "exact_transcript_chain_source_count",
    "exact_transcript_chain_sources", "qualifying_exact_transcript_group_count",
    "qualifying_exact_transcript_source_count", "qualifying_exact_transcript_sources",
    "qualifying_exact_transcript_group_ids", "exact_transcript_scoring_eligibility_status",
    "exact_transcript_stranded_requirement", "exact_transcript_stranded_status",
    "exact_transcript_cdna_context_status", "transcript_strong_partial",
    "transcript_evidence_score", "transcript_evidence_score_scaled",
    "transcript_best_pairwise_score", "transcript_mean_best_source_pairwise_score",
    "transcript_best_intron_chain_f1", "transcript_supported_cds_base_fraction",
    "transcript_scoring_supported_cds_base_fraction", "transcript_evidence_status",
    "transcript_evidence_flags",
    "sequence_integrity_status", "sequence_internal_stop_concern", "sequence_missing_start_concern",
    "sequence_missing_stop_concern", "sequence_phase_inconsistency_concern",
    "sequence_noncanonical_splice_concern", "sequence_neutral_boundary_or_ambiguity_state",
    "sequence_integrity_default_exclusion", "sequence_start_codon_status",
    "sequence_stop_codon_status", "sequence_phase_status", "sequence_internal_stop_count",
    "sequence_noncanonical_splice_intron_count", "sequence_integrity_flags",
    "intrinsic_cds_length_not_multiple_of_3", "intrinsic_short_cds",
    "intrinsic_short_cds_intron", "intrinsic_short_internal_cds_exon",
    "intrinsic_short_terminal_cds_exon", "intrinsic_long_cds_intron",
    "intrinsic_single_exon_model", "intrinsic_default_exclusion",
    "candidate_model_selection_score", "candidate_model_selection_score_scaled",
    "primary_model_selection_score", "primary_model_selection_score_scaled",
    "candidate_primary_model_selection_score_difference", "candidate_primary_model_selection_score_ratio",
    "candidate_primary_scaled_score_difference", "candidate_primary_scaled_score_ratio",
    "candidate_model_structural_support_score", "candidate_model_structural_support_score_scaled",
    "candidate_structural_support_score", "candidate_consensus_similarity_score",
    "candidate_pairwise_similarity_score", "candidate_model_completeness_score",
    "candidate_intrinsic_penalty_score", "candidate_protein_evidence_score_scaled",
    "candidate_transcript_evidence_score_scaled", "fixed_support_route_count",
    "fixed_support_routes", "fixed_rule_potential_candidate",
] + STAGE2B_CANDIDATE_COLUMNS


ALTERNATIVE_SPLICE_GENE_COLUMNS = [
    "primary_output_gene_id", "primary_output_transcript_id", "primary_source_name",
    "primary_input_gene_id", "primary_input_transcript_id", "chrom", "strand",
    "superlocus_id", "sublocus_id", "primary_cds_exon_count", "primary_cds_length",
    "structural_candidates_discovered", "uniquely_owned_candidates",
    "ambiguous_ownership_candidates", "exact_equivalent_candidates",
    "same_intron_chain_terminal_variant_candidates", "contained_candidates",
    "genuinely_splice_structure_distinct_candidates", "exact_multisource_supported_candidates",
    "protein_exact_supported_candidates", "transcript_exact_supported_candidates",
    "full_structural_mosaic_candidates", "protein_strong_partial_candidates",
    "transcript_strong_partial_candidates", "sequence_integrity_concern_candidates",
    "intrinsic_concern_candidates", "frame_status_exclusion_candidates",
    "fixed_rule_potential_candidates", "raw_candidate_population_above_50",
    "fixed_rule_candidate_population_above_50",
] + STAGE2B_GENE_COLUMNS


ALTERNATIVE_SPLICE_RUN_COLUMNS = [
    "requested_pruning_level", "effective_pruning_level",
    "prediction_source_count", "protein_evidence_source_count", "transcript_evidence_source_count",
    "genome_sequence_configured", "retained_primary_count", "valid_scored_model_count",
    "unselected_scored_model_count", "single_exon_unselected_candidate_count",
    "multi_exon_unselected_candidate_count", "zero_owner_candidate_count",
    "unique_owner_candidate_count", "ambiguous_owner_candidate_count",
    "candidate_gene_assessment_count", "exact_equivalent_candidate_count",
    "same_intron_chain_terminal_variant_candidate_count", "contained_candidate_count",
    "genuinely_splice_structure_distinct_candidate_count", "exact_multisource_model_support_count",
    "protein_exact_support_count", "transcript_exact_coherent_support_count",
    "full_structural_mosaic_count", "protein_strong_partial_count",
    "transcript_strong_partial_count", "sequence_integrity_default_exclusion_count",
    "intrinsic_default_exclusion_count", "frame_status_default_exclusion_count",
    "fixed_rule_potential_candidate_count",
    "genes_with_structural_candidates", "genes_with_candidate_population_above_50",
] + STAGE2B_RUN_COLUMNS


def _csv(values: set[str] | list[str]) -> str:
    ordered = sorted(set(values))
    return ",".join(ordered) if ordered else "."


def _count(row: dict[str, object], field: str) -> int:
    return int(value_as_float(row.get(field), 0.0))


def _feature_key(model: CodingModel) -> tuple[tuple[str, int, int, str], ...]:
    """Return non-CDS child geometry for the UTR-only diagnostic."""
    return tuple(sorted(
        (record.feature_type, record.start0, record.end0, record.strand)
        for record in model.child_records
        if record.feature_type != "CDS"
    ))


def _pair_overlap(primary: CodingModel, candidate: CodingModel) -> dict[str, object]:
    metrics = cds_exon_overlap_metrics(primary, candidate)
    primary_fraction = float(metrics["overlap_fraction_of_a"])
    candidate_fraction = float(metrics["overlap_fraction_of_b"])
    return {
        "cds_overlap_bp": int(metrics["cds_exon_overlap_bp"]),
        "primary_cds_overlap_fraction": primary_fraction,
        "candidate_cds_overlap_fraction": candidate_fraction,
        "reciprocal_cds_overlap": min(primary_fraction, candidate_fraction),
    }


def _owner_route(primary: CodingModel, candidate: CodingModel) -> str | None:
    if not candidate.is_multi_exon or primary.superlocus_id != candidate.superlocus_id:
        return None
    if primary.is_multi_exon:
        if primary.sublocus_id == candidate.sublocus_id and shared_cds_splice_site_count(primary, candidate) >= 1:
            return "multi_exon_same_sublocus_shared_cds_splice_site"
        return None
    if candidate.cds_exon_count != 2:
        return None
    overlap = _pair_overlap(primary, candidate)
    if (
        float(overlap["primary_cds_overlap_fraction"]) + 1e-15 >= 0.50
        and float(overlap["candidate_cds_overlap_fraction"]) + 1e-15 >= 0.50
    ):
        return "single_exon_primary_two_exon_candidate_reciprocal_overlap_ge_0.50"
    return None


def _pair_structure(primary: CodingModel, candidate: CodingModel, wobble_bp: int) -> dict[str, object]:
    exact = cds_exon_chains_match(primary, candidate, wobble_bp)
    exact_without_wobble = tuple(primary.sorted_cds_parts) == tuple(candidate.sorted_cds_parts)
    same_intron = cds_intron_chains_match(primary, candidate)
    containment = cds_exon_chain_containment_relation(candidate, primary)
    candidate_in_primary = containment["a_contained_in_b"]
    primary_in_candidate = containment["b_contained_in_a"]
    contained = candidate_in_primary or primary_in_candidate
    same_intron_terminal = same_intron and not exact
    return {
        "exact_cds_equivalent_to_primary": exact,
        "utr_only_difference_from_primary": (
            exact_without_wobble and _feature_key(primary) != _feature_key(candidate)
        ),
        "same_cds_intron_chain_as_primary": same_intron,
        "same_intron_chain_terminal_variant": same_intron_terminal,
        "candidate_cds_chain_contained_in_primary": candidate_in_primary,
        "primary_cds_chain_contained_in_candidate": primary_in_candidate,
        "bidirectional_cds_chain_containment": contained,
        "genuinely_splice_structure_distinct": not (exact or same_intron_terminal or contained),
    }


def _candidate_relationships(
    candidate: CodingModel,
    peers: list[CodingModel],
    wobble_bp: int,
) -> dict[str, object]:
    exact_ids: set[str] = set()
    same_intron_ids: set[str] = set()
    contained_in_ids: set[str] = set()
    contains_ids: set[str] = set()
    for other in peers:
        if other.output_transcript_id == candidate.output_transcript_id:
            continue
        if cds_exon_chains_match(candidate, other, wobble_bp):
            exact_ids.add(other.output_transcript_id)
        if cds_intron_chains_match(candidate, other):
            same_intron_ids.add(other.output_transcript_id)
        relation = cds_exon_chain_containment_relation(candidate, other)
        if relation["a_contained_in_b"]:
            contained_in_ids.add(other.output_transcript_id)
        if relation["b_contained_in_a"]:
            contains_ids.add(other.output_transcript_id)
    return {
        "exact_equivalent_candidate_ids": _csv(exact_ids),
        "same_intron_chain_candidate_ids": _csv(same_intron_ids),
        "candidate_contained_in_other_candidate_ids": _csv(contained_in_ids),
        "other_candidate_contained_in_candidate_ids": _csv(contains_ids),
    }


def _exact_model_support(
    candidate: CodingModel,
    possible_supporters: list[CodingModel],
    wobble_bp: int,
) -> dict[str, object]:
    supporters = [
        model for model in possible_supporters
        if model.output_transcript_id != candidate.output_transcript_id
        and cds_exon_chains_match(candidate, model, wobble_bp)
    ]
    other_sources = {model.source_name for model in supporters if model.source_name != candidate.source_name}
    all_sources = {candidate.source_name} | {model.source_name for model in supporters}
    return {
        "exact_model_other_source_count": len(other_sources),
        "exact_model_other_sources": _csv(other_sources),
        "exact_model_supporting_model_count": len(supporters),
        "exact_model_supporting_model_ids": _csv({model.output_transcript_id for model in supporters}),
        "exact_structural_all_source_count": len(all_sources),
        "exact_structural_all_sources": _csv(all_sources),
        "exact_multisource_model_support": bool(other_sources),
    }


def _exact_chain_invariant_key(item: object) -> tuple[object, ...]:
    """Return the coordinate key shared by every possible direct wobble match.

    The biological 3-prime outer boundary is omitted from the key and is still
    checked pairwise by the authoritative match helper.  This is an index, not
    a transitive equivalence-class construction.
    """
    parts = [list(part) for part in item.sorted_cds_parts]
    strand = str(item.strand)
    if parts:
        if strand == "+":
            parts[-1][1] = -1
        elif strand == "-":
            parts[0][0] = -1
    return (
        str(item.chrom),
        strand,
        tuple((start0, end0) for start0, end0 in parts),
    )


def _component_state(total: int, supported: int) -> str:
    if total <= 0:
        return "not_applicable"
    return "complete" if supported >= total else "incomplete"


def _full_structural_mosaic(
    score_row: dict[str, object],
    *,
    has_exact_other_source: bool,
) -> dict[str, object]:
    states = {
        "mosaic_cds_introns_status": _component_state(
            _count(score_row, "cds_introns_total"),
            _count(score_row, "cds_introns_supported_by_other_source"),
        ),
        "mosaic_cds_splice_sites_status": _component_state(
            _count(score_row, "cds_splice_sites_total"),
            _count(score_row, "cds_splice_sites_supported_by_other_source"),
        ),
        "mosaic_internal_cds_exons_status": _component_state(
            _count(score_row, "internal_cds_exons_total"),
            _count(score_row, "internal_cds_exons_supported_by_other_source"),
        ),
        "mosaic_terminal_cds_exons_status": _component_state(
            _count(score_row, "terminal_cds_exons_total"),
            _count(score_row, "terminal_cds_exons_supported_by_other_source"),
        ),
    }
    base_precision = score_row.get("cds_base_consensus_precision")
    base_complete = base_precision not in {None, "."} and value_as_float(base_precision, -1.0) + 1e-15 >= 1.0
    states["mosaic_cds_base_consensus_status"] = "complete" if base_complete else "incomplete"
    applicable_states = [state for state in states.values() if state != "not_applicable"]
    states["full_structural_mosaic"] = bool(
        not has_exact_other_source
        and applicable_states
        and all(state == "complete" for state in applicable_states)
    )
    return states


def _model_support_group_key(model: CodingModel) -> tuple[str, str, str, str]:
    """Return the structural-scoring context used for model-source consensus."""
    return (
        model.chrom,
        model.strand,
        model.superlocus_id or ".",
        model.sublocus_id or ".",
    )


def _source_depth_segments(models: list[CodingModel]) -> tuple[tuple[int, int, int], ...]:
    """Return disjoint genomic segments labelled by distinct model-source depth."""
    intervals_by_source: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for model in models:
        intervals_by_source[model.source_name].extend(model.sorted_cds_parts)
    events: Counter[int] = Counter()
    for intervals in intervals_by_source.values():
        for start0, end0 in merge_intervals(intervals):
            events[start0] += 1
            events[end0] -= 1
    depth = 0
    previous: int | None = None
    segments: list[tuple[int, int, int]] = []
    for position in sorted(events):
        if previous is not None and previous < position and depth > 0:
            if segments and segments[-1][1] == previous and segments[-1][2] == depth:
                segments[-1] = (segments[-1][0], position, depth)
            else:
                segments.append((previous, position, depth))
        depth += events[position]
        previous = position
    return tuple(segments)


def _minimum_cds_base_source_count(
    candidate: CodingModel,
    depth_segments: tuple[tuple[int, int, int], ...],
) -> int:
    """Return minimum distinct prediction-source depth over candidate CDS bases."""
    minimum: int | None = None
    for part_start, part_end in candidate.sorted_cds_parts:
        cursor = part_start
        for segment_start, segment_end, depth in depth_segments:
            if segment_end <= part_start:
                continue
            if segment_start >= part_end:
                break
            overlap_start = max(part_start, segment_start)
            overlap_end = min(part_end, segment_end)
            if overlap_end <= overlap_start:
                continue
            if overlap_start > cursor:
                return 0
            minimum = depth if minimum is None else min(minimum, depth)
            cursor = max(cursor, overlap_end)
            if cursor >= part_end:
                break
        if cursor < part_end:
            return 0
    return minimum or 0


def _model_support_source_depth(
    candidate: CodingModel,
    score_row: dict[str, object],
    exact_model: dict[str, object],
    depth_segments: tuple[tuple[int, int, int], ...],
    prediction_source_count: int,
) -> dict[str, object]:
    """Return raw exact-model and minimum mosaic component source depths."""
    exact_count = int(exact_model["exact_structural_all_source_count"])
    component_counts: list[int] = []
    for total_field, minimum_field in (
        ("cds_introns_total", "cds_intron_support_source_min"),
        ("cds_splice_sites_total", "cds_splice_site_support_source_min"),
        ("internal_cds_exons_total", "internal_cds_exon_support_source_min"),
        ("terminal_cds_exons_total", "terminal_cds_exon_support_source_min"),
    ):
        if _count(score_row, total_field) > 0:
            component_counts.append(_count(score_row, minimum_field))
    component_counts.append(
        _minimum_cds_base_source_count(candidate, depth_segments)
    )
    mosaic_minimum = min(component_counts)
    return {
        "alternative_splice_exact_model_source_count": exact_count,
        "alternative_splice_exact_model_source_fraction": (
            exact_count / prediction_source_count if prediction_source_count else 0.0
        ),
        "alternative_splice_mosaic_min_component_source_count": mosaic_minimum,
        "alternative_splice_mosaic_min_component_source_fraction": (
            mosaic_minimum / prediction_source_count
            if prediction_source_count else 0.0
        ),
    }


def _protein_support(
    candidate: CodingModel,
    protein_row: dict[str, object],
    alignments: list[ProteinEvidenceCollapsedAlignment],
    wobble_bp: int,
) -> dict[str, object]:
    exact = [
        alignment for alignment in alignments
        if chains_match_with_3prime_wobble(candidate, alignment, wobble_bp)
    ]
    sources = {source for alignment in exact for source in alignment.sources}
    flags = str(protein_row.get("protein_evidence_flags", "."))
    exact_supported = bool(exact)
    return {
        "protein_exact_cds_chain_support": exact_supported,
        "protein_exact_supporting_source_count": len(sources),
        "protein_exact_supporting_sources": _csv(sources),
        "protein_exact_supporting_group_count": len(exact),
        "protein_exact_supporting_group_ids": _csv({alignment.collapsed_id for alignment in exact}),
        "protein_effective_terminal_wobble_bp": wobble_bp,
        "protein_strong_partial": (not exact_supported and "strong_protein_support" in flags.split(",")),
        "protein_evidence_score": protein_row.get("protein_evidence_score", 0.0),
        "protein_evidence_score_scaled": protein_row.get("protein_evidence_score_scaled", 0.0),
        "protein_best_alignment_score": protein_row.get("protein_best_alignment_score", 0.0),
        "protein_mean_best_source_score": protein_row.get("protein_mean_best_source_score", 0.0),
        "protein_supported_cds_base_fraction": protein_row.get("protein_supported_cds_base_fraction", 0.0),
        "protein_cds_intron_f1": protein_row.get("protein_cds_intron_f1", 0.0),
        "protein_cds_splice_site_f1": protein_row.get("protein_cds_splice_site_f1", 0.0),
        "protein_internal_cds_exon_f1": protein_row.get("protein_internal_cds_exon_f1", 0.0),
        "protein_terminal_cds_exon_f1": protein_row.get("protein_terminal_cds_exon_f1", 0.0),
        "protein_evidence_status": protein_row.get("protein_evidence_status", "protein_evidence_not_configured"),
        "protein_evidence_flags": flags,
    }


def _transcript_support(
    candidate: CodingModel,
    transcript_row: dict[str, object],
    transcript_index: IntervalOverlapIndex[str, TranscriptEvidenceSuperlocus] | None,
    policy: AlternativeSpliceReportPolicy,
) -> dict[str, object]:
    associations = associated_transcripts_for_model(candidate, transcript_index) if transcript_index is not None else []
    coherent = evaluate_coherent_exact_transcript_support(
        candidate,
        associations,
        min_cdna_context_fraction=policy.transcript_min_cdna_context_fraction,
        require_stranded_support=policy.transcript_require_stranded_support,
        block_on_cdna_context_concern=policy.transcript_block_on_cdna_context_concern,
        cdna_context_neutral_fraction=policy.transcript_cdna_context_neutral_fraction,
        cdna_context_severe_fraction=policy.transcript_cdna_context_severe_fraction,
    )
    exact_supported = bool(coherent["exact_coherent_transcript_support"])
    transcript_status = str(transcript_row.get("transcript_evidence_status", "transcript_evidence_not_configured"))
    coherent.update({
        "transcript_strong_partial": (not exact_supported and transcript_status == "strong_transcript_support"),
        "transcript_evidence_score": transcript_row.get("transcript_evidence_score", 0.0),
        "transcript_evidence_score_scaled": transcript_row.get("transcript_evidence_score_scaled", 0.0),
        "transcript_best_pairwise_score": transcript_row.get("transcript_best_pairwise_score", 0.0),
        "transcript_mean_best_source_pairwise_score": transcript_row.get("transcript_mean_best_source_pairwise_score", 0.0),
        "transcript_best_intron_chain_f1": transcript_row.get("transcript_best_intron_chain_f1", 0.0),
        "transcript_supported_cds_base_fraction": transcript_row.get("transcript_supported_cds_base_fraction", 0.0),
        "transcript_scoring_supported_cds_base_fraction": transcript_row.get("transcript_scoring_supported_cds_base_fraction", 0.0),
        "transcript_evidence_status": transcript_status,
        "transcript_evidence_flags": transcript_row.get("transcript_evidence_flags", "."),
    })
    return coherent


def _sequence_integrity(candidate: CodingModel) -> dict[str, object]:
    metrics = candidate.sequence_integrity
    if metrics is None:
        return {
            "sequence_integrity_status": "unavailable",
            "sequence_internal_stop_concern": False,
            "sequence_missing_start_concern": False,
            "sequence_missing_stop_concern": False,
            "sequence_phase_inconsistency_concern": False,
            "sequence_noncanonical_splice_concern": False,
            "sequence_neutral_boundary_or_ambiguity_state": "unavailable_genome_sequence_neutral",
            "sequence_integrity_default_exclusion": False,
            "sequence_start_codon_status": ".", "sequence_stop_codon_status": ".",
            "sequence_phase_status": ".", "sequence_internal_stop_count": ".",
            "sequence_noncanonical_splice_intron_count": ".", "sequence_integrity_flags": ".",
        }
    internal = metrics.internal_stop_count > 0
    missing_start = metrics.start_codon_status == "missing"
    missing_stop = metrics.stop_codon_status == "missing"
    phase = metrics.phase_status == "inconsistent"
    splice = metrics.noncanonical_splice_intron_count > 0
    flags = {flag for flag in metrics.flags.split(",") if flag and flag != "."}
    neutral = sorted(flags & {
        "ambiguous_sequence", "five_prime_contig_boundary", "three_prime_contig_boundary",
        "terminal_stop_adjacent_to_cds",
    })
    return {
        "sequence_integrity_status": metrics.status,
        "sequence_internal_stop_concern": internal,
        "sequence_missing_start_concern": missing_start,
        "sequence_missing_stop_concern": missing_stop,
        "sequence_phase_inconsistency_concern": phase,
        "sequence_noncanonical_splice_concern": splice,
        "sequence_neutral_boundary_or_ambiguity_state": _csv(neutral),
        "sequence_integrity_default_exclusion": any((internal, missing_start, missing_stop, phase, splice)),
        "sequence_start_codon_status": metrics.start_codon_status,
        "sequence_stop_codon_status": metrics.stop_codon_status,
        "sequence_phase_status": metrics.phase_status,
        "sequence_internal_stop_count": metrics.internal_stop_count,
        "sequence_noncanonical_splice_intron_count": metrics.noncanonical_splice_intron_count,
        "sequence_integrity_flags": metrics.flags,
    }


def _intrinsic(score_row: dict[str, object], candidate: CodingModel) -> dict[str, object]:
    fields = {
        "intrinsic_cds_length_not_multiple_of_3": _count(score_row, "intrinsic_penalty_cds_length_not_multiple_of_3_count"),
        "intrinsic_short_cds": _count(score_row, "intrinsic_penalty_short_cds_count"),
        "intrinsic_short_cds_intron": _count(score_row, "intrinsic_penalty_short_cds_intron_count"),
        "intrinsic_short_internal_cds_exon": _count(score_row, "intrinsic_penalty_short_internal_cds_exon_count"),
        "intrinsic_short_terminal_cds_exon": _count(score_row, "intrinsic_penalty_short_terminal_cds_exon_count"),
        "intrinsic_long_cds_intron": _count(score_row, "intrinsic_penalty_long_cds_intron_count"),
        "intrinsic_single_exon_model": candidate.is_single_exon,
    }
    fields["intrinsic_default_exclusion"] = any(
        int(fields[name]) > 0
        for name in (
            "intrinsic_cds_length_not_multiple_of_3", "intrinsic_short_cds",
            "intrinsic_short_cds_intron", "intrinsic_short_internal_cds_exon",
            "intrinsic_short_terminal_cds_exon", "intrinsic_long_cds_intron",
        )
    )
    return fields


def _scores(candidate_row: dict[str, object], primary_row: dict[str, object] | None) -> dict[str, object]:
    candidate_raw = value_as_float(candidate_row.get("model_selection_score"), 0.0)
    candidate_scaled = value_as_float(candidate_row.get("model_selection_score_scaled"), 0.0)
    if primary_row is None:
        primary_raw: object = "."
        primary_scaled: object = "."
        raw_diff: object = "."
        raw_ratio: object = "."
        scaled_diff: object = "."
        scaled_ratio: object = "."
    else:
        primary_raw = value_as_float(primary_row.get("model_selection_score"), 0.0)
        primary_scaled = value_as_float(primary_row.get("model_selection_score_scaled"), 0.0)
        raw_diff = candidate_raw - float(primary_raw)
        raw_ratio = candidate_raw / float(primary_raw) if float(primary_raw) != 0 else "."
        scaled_diff = candidate_scaled - float(primary_scaled)
        scaled_ratio = candidate_scaled / float(primary_scaled) if float(primary_scaled) != 0 else "."
    return {
        "candidate_model_selection_score": candidate_raw,
        "candidate_model_selection_score_scaled": candidate_scaled,
        "primary_model_selection_score": primary_raw,
        "primary_model_selection_score_scaled": primary_scaled,
        "candidate_primary_model_selection_score_difference": raw_diff,
        "candidate_primary_model_selection_score_ratio": raw_ratio,
        "candidate_primary_scaled_score_difference": scaled_diff,
        "candidate_primary_scaled_score_ratio": scaled_ratio,
        "candidate_model_structural_support_score": candidate_row.get("model_structural_support_score", "."),
        "candidate_model_structural_support_score_scaled": candidate_row.get("model_structural_support_score_scaled", "."),
        "candidate_structural_support_score": candidate_row.get("structural_support_score", "."),
        "candidate_consensus_similarity_score": candidate_row.get("consensus_similarity_score", "."),
        "candidate_pairwise_similarity_score": candidate_row.get("pairwise_similarity_score", "."),
        "candidate_model_completeness_score": candidate_row.get("model_completeness_score", "."),
        "candidate_intrinsic_penalty_score": candidate_row.get("intrinsic_penalty_score", "."),
        "candidate_protein_evidence_score_scaled": candidate_row.get("protein_evidence_score_scaled", 0.0),
        "candidate_transcript_evidence_score_scaled": candidate_row.get("transcript_evidence_score_scaled", 0.0),
    }


def _primary_fields(primary: CodingModel | None) -> dict[str, object]:
    if primary is None:
        return {field: "." for field in (
            "primary_source_name", "primary_input_gene_id", "primary_input_transcript_id",
            "primary_output_gene_id", "primary_output_transcript_id", "primary_superlocus_id",
            "primary_sublocus_id", "primary_cds_start", "primary_cds_end", "primary_cds_length",
            "primary_cds_exon_count", "primary_cds_intron_count",
        )}
    return {
        "primary_source_name": primary.source_name,
        "primary_input_gene_id": primary.input_gene_id,
        "primary_input_transcript_id": primary.input_transcript_id,
        "primary_output_gene_id": primary.output_gene_id,
        "primary_output_transcript_id": primary.output_transcript_id,
        "primary_superlocus_id": primary.superlocus_id or ".",
        "primary_sublocus_id": primary.sublocus_id or ".",
        "primary_cds_start": primary.cds_start0 + 1,
        "primary_cds_end": primary.cds_end0,
        "primary_cds_length": primary.cds_length,
        "primary_cds_exon_count": primary.cds_exon_count,
        "primary_cds_intron_count": len(primary.cds_introns),
    }


def _final_states(row: dict[str, object]) -> dict[str, object]:
    """Return the historical fixed-rule-potential diagnostic only."""
    ownership = str(row["ownership_status"])
    if row["candidate_cds_exon_count"] == 1:
        overall = "rejected_single_exon_candidate"
    elif ownership == "zero_owner":
        overall = "rejected_no_retained_gene_owner"
    elif ownership == "ambiguous_multiple_owners":
        overall = "rejected_ambiguous_ownership"
    elif bool(row["exact_cds_equivalent_to_primary"]):
        overall = "rejected_exact_cds_equivalent"
    elif bool(row["same_intron_chain_terminal_variant"]):
        overall = "rejected_same_intron_chain_terminal_variant"
    elif bool(row["bidirectional_cds_chain_containment"]):
        overall = "rejected_cds_chain_containment"
    elif str(row["frame_concordance_status"]) != "calculated":
        overall = "rejected_default_frame_status"
    elif bool(row["sequence_integrity_default_exclusion"]):
        overall = "rejected_sequence_integrity_concern"
    elif bool(row["intrinsic_default_exclusion"]):
        overall = "rejected_intrinsic_concern"
    elif int(row["fixed_support_route_count"]) <= 0:
        overall = "rejected_no_qualifying_exact_support_route"
    else:
        overall = "pending_tunable_thresholds"
    pending = overall == "pending_tunable_thresholds"
    return {"fixed_rule_potential_candidate": pending}


def analyze_alternative_splice_report(
    *,
    anchor_models: list[CodingModel],
    all_scored_models: list[CodingModel],
    model_score_rows_by_transcript: dict[str, dict[str, object]],
    protein_rows_by_transcript: dict[str, dict[str, object]] | None = None,
    protein_collapsed_alignments: list[ProteinEvidenceCollapsedAlignment] | None = None,
    transcript_rows_by_transcript: dict[str, dict[str, object]] | None = None,
    transcript_superloci: list[TranscriptEvidenceSuperlocus] | None = None,
    policy: AlternativeSpliceReportPolicy,
    run_metadata: dict[str, object] | None = None,
) -> dict[str, list[dict[str, object]]]:
    """Return candidate, retained-gene and run report rows without mutation."""
    protein_rows_by_transcript = protein_rows_by_transcript or {}
    transcript_rows_by_transcript = transcript_rows_by_transcript or {}
    protein_collapsed_alignments = protein_collapsed_alignments or []
    anchor_ids = {model.output_transcript_id for model in anchor_models}
    candidates = [
        model for model in all_scored_models
        if model.output_transcript_id not in anchor_ids
        and model.output_transcript_id in model_score_rows_by_transcript
    ]
    anchors_by_superlocus: dict[str, list[CodingModel]] = defaultdict(list)
    for anchor in anchor_models:
        anchors_by_superlocus[anchor.superlocus_id or "."].append(anchor)
    potential_owners: dict[str, list[tuple[CodingModel, str]]] = {}
    for candidate in candidates:
        owners: list[tuple[CodingModel, str]] = []
        for primary in anchors_by_superlocus.get(candidate.superlocus_id or ".", []):
            route = _owner_route(primary, candidate)
            if route is not None:
                owners.append((primary, route))
        potential_owners[candidate.output_transcript_id] = sorted(
            owners, key=lambda item: item[0].output_transcript_id
        )

    unique_candidates_by_owner: dict[str, list[CodingModel]] = defaultdict(list)
    for candidate in candidates:
        owners = potential_owners[candidate.output_transcript_id]
        if len(owners) == 1:
            unique_candidates_by_owner[owners[0][0].output_transcript_id].append(candidate)
    alternative_rank: dict[tuple[str, str], int] = {}
    for primary_id, items in unique_candidates_by_owner.items():
        ordered = sorted(
            items,
            key=lambda model: model_selection_row_sort_key(
                model_score_rows_by_transcript[model.output_transcript_id]
            ),
        )
        for rank, candidate in enumerate(ordered, start=1):
            alternative_rank[(primary_id, candidate.output_transcript_id)] = rank

    transcript_index: IntervalOverlapIndex[str, TranscriptEvidenceSuperlocus] | None = None
    if transcript_superloci:
        transcript_index = IntervalOverlapIndex(
            transcript_superloci,
            key=lambda locus: locus.chrom,
            start=lambda locus: locus.start0,
            end=lambda locus: locus.end0,
        )
    protein_by_exact_invariant: dict[tuple[object, ...], list[ProteinEvidenceCollapsedAlignment]] = defaultdict(list)
    for alignment in protein_collapsed_alignments:
        protein_by_exact_invariant[_exact_chain_invariant_key(alignment)].append(alignment)
    models_by_exact_invariant: dict[tuple[object, ...], list[CodingModel]] = defaultdict(list)
    models_by_support_group: dict[
        tuple[str, str, str, str], list[CodingModel]
    ] = defaultdict(list)
    for model in all_scored_models:
        models_by_exact_invariant[_exact_chain_invariant_key(model)].append(model)
        models_by_support_group[_model_support_group_key(model)].append(model)
    source_depth_by_support_group: dict[
        tuple[str, str, str, str], tuple[tuple[int, int, int], ...]
    ] = {}
    prediction_source_count = len(
        {model.source_name for model in all_scored_models}
    )

    candidate_rows: list[dict[str, object]] = []
    for candidate in sorted(
        candidates,
        key=lambda model: (model.chrom, model.cds_start0, model.cds_end0, model.output_transcript_id),
    ):
        candidate_id = candidate.output_transcript_id
        candidate_score_row = model_score_rows_by_transcript[candidate_id]
        owners = potential_owners[candidate_id]
        ownership_status = (
            "zero_owner" if not owners else
            "unique_owner" if len(owners) == 1 else
            "ambiguous_multiple_owners"
        )
        associated_ids = {primary.output_transcript_id for primary, _route in owners}
        associated_sources = {primary.source_name for primary, _route in owners}
        exact_model = _exact_model_support(
            candidate,
            models_by_exact_invariant[_exact_chain_invariant_key(candidate)],
            policy.terminal_cds_3prime_wobble_bp,
        )
        support_group_key = _model_support_group_key(candidate)
        if support_group_key not in source_depth_by_support_group:
            source_depth_by_support_group[support_group_key] = _source_depth_segments(
                models_by_support_group[support_group_key]
            )
        model_support_source_depth = _model_support_source_depth(
            candidate,
            candidate_score_row,
            exact_model,
            source_depth_by_support_group[support_group_key],
            prediction_source_count,
        )
        mosaic = _full_structural_mosaic(
            candidate_score_row,
            has_exact_other_source=bool(exact_model["exact_multisource_model_support"]),
        )
        protein = _protein_support(
            candidate,
            protein_rows_by_transcript.get(candidate_id, {}),
            protein_by_exact_invariant.get(_exact_chain_invariant_key(candidate), []),
            policy.protein_terminal_cds_3prime_wobble_bp,
        )
        transcript = _transcript_support(
            candidate,
            transcript_rows_by_transcript.get(candidate_id, {}),
            transcript_index,
            policy,
        )
        sequence = _sequence_integrity(candidate)
        intrinsic = _intrinsic(candidate_score_row, candidate)
        assessment_owners: list[tuple[CodingModel | None, str]] = (
            [(primary, route) for primary, route in owners] if owners else [(None, ".")]
        )
        for primary, route in assessment_owners:
            if primary is None:
                overlap = {
                    "cds_overlap_bp": ".", "primary_cds_overlap_fraction": ".",
                    "candidate_cds_overlap_fraction": ".", "reciprocal_cds_overlap": ".",
                }
                frame = {field: "." for field in (
                    "frame_concordance_status", "primary_phase_status", "candidate_phase_status",
                    "primary_phase_mismatch_count", "candidate_phase_mismatch_count", "frame_overlap_bp",
                    "frame_comparable_overlap_bp", "same_frame_overlap_bp", "same_frame_comparable_fraction",
                    "same_frame_overlap_fraction", "same_frame_primary_fraction", "same_frame_candidate_fraction",
                )}
                structure = {
                    "exact_cds_equivalent_to_primary": False, "utr_only_difference_from_primary": False,
                    "same_cds_intron_chain_as_primary": False, "same_intron_chain_terminal_variant": False,
                    "candidate_cds_chain_contained_in_primary": False,
                    "primary_cds_chain_contained_in_candidate": False,
                    "bidirectional_cds_chain_containment": False,
                    "genuinely_splice_structure_distinct": False,
                }
                relationships = {
                    "exact_equivalent_candidate_ids": ".", "same_intron_chain_candidate_ids": ".",
                    "candidate_contained_in_other_candidate_ids": ".",
                    "other_candidate_contained_in_candidate_ids": ".",
                }
                rank: object = "."
                primary_score_row = None
                shared_splice: object = "."
                assessment_id = f"alt_assessment:{candidate_id}:no_owner"
            else:
                overlap = _pair_overlap(primary, candidate)
                frame = cds_frame_concordance(primary, candidate)
                structure = _pair_structure(primary, candidate, policy.terminal_cds_3prime_wobble_bp)
                relationships = _candidate_relationships(
                    candidate,
                    unique_candidates_by_owner[primary.output_transcript_id]
                    if ownership_status == "unique_owner" else [],
                    policy.terminal_cds_3prime_wobble_bp,
                )
                rank = alternative_rank.get((primary.output_transcript_id, candidate_id), ".")
                primary_score_row = model_score_rows_by_transcript.get(primary.output_transcript_id)
                shared_splice = shared_cds_splice_site_count(primary, candidate)
                assessment_id = f"alt_assessment:{candidate_id}:{primary.output_transcript_id}"
            support_routes = []
            if bool(exact_model["exact_multisource_model_support"]):
                support_routes.append("exact_multisource_model_support")
            if bool(protein["protein_exact_cds_chain_support"]):
                support_routes.append("exact_protein_support")
            if bool(transcript["exact_coherent_transcript_support"]):
                support_routes.append("exact_coherent_transcript_support")
            row: dict[str, object] = {
                "alternative_candidate_assessment_id": assessment_id,
                "candidate_source_name": candidate.source_name,
                "candidate_input_file": str(candidate.input_file),
                "candidate_input_gene_id": candidate.input_gene_id,
                "candidate_input_transcript_id": candidate.input_transcript_id,
                "candidate_output_gene_id": candidate.output_gene_id,
                "candidate_output_transcript_id": candidate_id,
                "chrom": candidate.chrom, "strand": candidate.strand,
                "candidate_superlocus_id": candidate.superlocus_id or ".",
                "candidate_sublocus_id": candidate.sublocus_id or ".",
                "candidate_cds_start": candidate.cds_start0 + 1,
                "candidate_cds_end": candidate.cds_end0,
                "candidate_cds_length": candidate.cds_length,
                "candidate_cds_exon_count": candidate.cds_exon_count,
                "candidate_cds_intron_count": len(candidate.cds_introns),
                **_primary_fields(primary),
                "structurally_associated_retained_gene_count": len(owners),
                "structurally_associated_primary_ids": _csv(associated_ids),
                "structurally_associated_primary_sources": _csv(associated_sources),
                "ownership_status": ownership_status,
                "structural_association_route": route,
                "structural_candidate": primary is not None,
                "shared_cds_splice_site_count": shared_splice,
                **overlap, **frame, **structure, **relationships,
                "alternative_candidate_order_rank": rank,
                "historical_model_selection_order_rank": candidate_score_row.get("model_selection_order_rank", "."),
                **exact_model, **model_support_source_depth, **mosaic,
                **protein, **transcript, **sequence, **intrinsic,
                **_scores(candidate_score_row, primary_score_row),
                "fixed_support_route_count": len(support_routes),
                "fixed_support_routes": ",".join(support_routes) if support_routes else ".",
            }
            row.update(_final_states(row))
            candidate_rows.append(row)

    gene_rows: list[dict[str, object]] = []
    rows_by_primary: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in candidate_rows:
        primary_id = str(row.get("primary_output_transcript_id", "."))
        if primary_id != ".":
            rows_by_primary[primary_id].append(row)
    for primary in sorted(anchor_models, key=lambda model: (model.chrom, model.cds_start0, model.output_transcript_id)):
        rows = rows_by_primary.get(primary.output_transcript_id, [])
        candidate_ids = {str(row["candidate_output_transcript_id"]) for row in rows}

        def ids_where(field: str, assessment_rows: list[dict[str, object]] = rows) -> set[str]:
            return {
                str(row["candidate_output_transcript_id"])
                for row in assessment_rows if bool(row.get(field))
            }

        fixed_ids = ids_where("fixed_rule_potential_candidate")
        gene_rows.append({
            "primary_output_gene_id": primary.output_gene_id,
            "primary_output_transcript_id": primary.output_transcript_id,
            "primary_source_name": primary.source_name,
            "primary_input_gene_id": primary.input_gene_id,
            "primary_input_transcript_id": primary.input_transcript_id,
            "chrom": primary.chrom, "strand": primary.strand,
            "superlocus_id": primary.superlocus_id or ".", "sublocus_id": primary.sublocus_id or ".",
            "primary_cds_exon_count": primary.cds_exon_count, "primary_cds_length": primary.cds_length,
            "structural_candidates_discovered": len(candidate_ids),
            "uniquely_owned_candidates": len({str(row["candidate_output_transcript_id"]) for row in rows if row["ownership_status"] == "unique_owner"}),
            "ambiguous_ownership_candidates": len({str(row["candidate_output_transcript_id"]) for row in rows if row["ownership_status"] == "ambiguous_multiple_owners"}),
            "exact_equivalent_candidates": len(ids_where("exact_cds_equivalent_to_primary")),
            "same_intron_chain_terminal_variant_candidates": len(ids_where("same_intron_chain_terminal_variant")),
            "contained_candidates": len(ids_where("bidirectional_cds_chain_containment")),
            "genuinely_splice_structure_distinct_candidates": len(ids_where("genuinely_splice_structure_distinct")),
            "exact_multisource_supported_candidates": len(ids_where("exact_multisource_model_support")),
            "protein_exact_supported_candidates": len(ids_where("protein_exact_cds_chain_support")),
            "transcript_exact_supported_candidates": len(ids_where("exact_coherent_transcript_support")),
            "full_structural_mosaic_candidates": len(ids_where("full_structural_mosaic")),
            "protein_strong_partial_candidates": len(ids_where("protein_strong_partial")),
            "transcript_strong_partial_candidates": len(ids_where("transcript_strong_partial")),
            "sequence_integrity_concern_candidates": len(ids_where("sequence_integrity_default_exclusion")),
            "intrinsic_concern_candidates": len(ids_where("intrinsic_default_exclusion")),
            "frame_status_exclusion_candidates": len({str(row["candidate_output_transcript_id"]) for row in rows if row["ownership_status"] == "unique_owner" and row["frame_concordance_status"] != "calculated"}),
            "fixed_rule_potential_candidates": len(fixed_ids),
            "raw_candidate_population_above_50": len(candidate_ids) > policy.proposed_max_alternatives_per_gene,
            "fixed_rule_candidate_population_above_50": len(fixed_ids) > policy.proposed_max_alternatives_per_gene,
        })

    rows_by_candidate: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in candidate_rows:
        rows_by_candidate[str(row["candidate_output_transcript_id"])].append(row)

    def candidate_count(field: str) -> int:
        return sum(any(bool(row.get(field)) for row in rows) for rows in rows_by_candidate.values())

    ownership_counts = Counter(
        str(rows[0]["ownership_status"]) for rows in rows_by_candidate.values()
    )
    metadata = run_metadata or {}
    run_row = {
        "requested_pruning_level": metadata.get("requested_pruning_level", "."),
        "effective_pruning_level": metadata.get("effective_pruning_level", "."),
        "prediction_source_count": prediction_source_count,
        "protein_evidence_source_count": metadata.get("protein_evidence_source_count", 0),
        "transcript_evidence_source_count": metadata.get("transcript_evidence_source_count", 0),
        "genome_sequence_configured": bool(metadata.get("genome_sequence_configured", False)),
        "retained_primary_count": len(anchor_models),
        "valid_scored_model_count": len(all_scored_models),
        "unselected_scored_model_count": len(candidates),
        "single_exon_unselected_candidate_count": sum(model.is_single_exon for model in candidates),
        "multi_exon_unselected_candidate_count": sum(model.is_multi_exon for model in candidates),
        "zero_owner_candidate_count": ownership_counts["zero_owner"],
        "unique_owner_candidate_count": ownership_counts["unique_owner"],
        "ambiguous_owner_candidate_count": ownership_counts["ambiguous_multiple_owners"],
        "candidate_gene_assessment_count": len(candidate_rows),
        "exact_equivalent_candidate_count": candidate_count("exact_cds_equivalent_to_primary"),
        "same_intron_chain_terminal_variant_candidate_count": candidate_count("same_intron_chain_terminal_variant"),
        "contained_candidate_count": candidate_count("bidirectional_cds_chain_containment"),
        "genuinely_splice_structure_distinct_candidate_count": candidate_count("genuinely_splice_structure_distinct"),
        "exact_multisource_model_support_count": candidate_count("exact_multisource_model_support"),
        "protein_exact_support_count": candidate_count("protein_exact_cds_chain_support"),
        "transcript_exact_coherent_support_count": candidate_count("exact_coherent_transcript_support"),
        "full_structural_mosaic_count": candidate_count("full_structural_mosaic"),
        "protein_strong_partial_count": candidate_count("protein_strong_partial"),
        "transcript_strong_partial_count": candidate_count("transcript_strong_partial"),
        "sequence_integrity_default_exclusion_count": candidate_count("sequence_integrity_default_exclusion"),
        "intrinsic_default_exclusion_count": candidate_count("intrinsic_default_exclusion"),
        "frame_status_default_exclusion_count": sum(
            rows[0]["ownership_status"] == "unique_owner"
            and str(rows[0]["frame_concordance_status"]) != "calculated"
            for rows in rows_by_candidate.values()
        ),
        "fixed_rule_potential_candidate_count": candidate_count("fixed_rule_potential_candidate"),
        "genes_with_structural_candidates": sum(int(row["structural_candidates_discovered"]) > 0 for row in gene_rows),
        "genes_with_candidate_population_above_50": sum(bool(row["raw_candidate_population_above_50"]) for row in gene_rows),
    }
    # The selector appends active decision values in a separate pure step.
    for row in candidate_rows:
        for field in STAGE2B_CANDIDATE_COLUMNS:
            row.setdefault(field, ".")
    for row in gene_rows:
        for field in STAGE2B_GENE_COLUMNS:
            row.setdefault(field, ".")
    for field in STAGE2B_RUN_COLUMNS:
        run_row.setdefault(field, ".")
    return {"candidate_rows": candidate_rows, "gene_rows": gene_rows, "run_rows": [run_row]}


def write_alternative_splice_reports(
    result: dict[str, list[dict[str, object]]],
    *,
    candidate_tsv: Path | None,
    gene_tsv: Path | None,
    run_tsv: Path | None,
) -> None:
    """Write alternative-splice facts and selector decisions."""
    if candidate_tsv is not None:
        write_dict_tsv(result["candidate_rows"], candidate_tsv, ALTERNATIVE_SPLICE_CANDIDATE_COLUMNS)
    if gene_tsv is not None:
        write_dict_tsv(result["gene_rows"], gene_tsv, ALTERNATIVE_SPLICE_GENE_COLUMNS)
    if run_tsv is not None:
        write_dict_tsv(result["run_rows"], run_tsv, ALTERNATIVE_SPLICE_RUN_COLUMNS)
