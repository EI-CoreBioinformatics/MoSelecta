"""Structural-support metrics and related report writers for MoSelecta."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from .constants import (
    COMPLETENESS_PENALTY_WEIGHTS,
    CONSENSUS_MULTI_SCORE_WEIGHTS,
    CONSENSUS_SINGLE_SCORE_WEIGHTS,
    PAIRWISE_MULTI_SCORE_WEIGHTS,
    PAIRWISE_SINGLE_SCORE_WEIGHTS,
    STRUCTURAL_MULTI_SCORE_WEIGHTS,
    STRUCTURAL_SINGLE_SCORE_WEIGHTS,
)
from .console_reporting import summarise_source_scaled_scores
from .corroboration import (
    compare_structural_corroboration,
    make_structural_corroboration_row,
)

if TYPE_CHECKING:
    from .log_utils import ProgressLogger
from .intervals import (
    classify_intervals_relative_to_span,
    consensus_intervals_from_source_intervals,
    exact_feature_prf,
    five_three_from_before_after,
    fraction_or_none,
    interval_list_overlap_bp,
    interval_prf,
    interval_total_bp,
    prf_from_counts,
    subtract_intervals,
)
from .loci import model_key, reciprocal_overlap
from .models import CodingModel, Sublocus
from .scoring import (
    fmt_model_ids,
    fmt_optional_float,
    fmt_sources,
    prf_values,
    safe_mean,
    source_scaled_fraction,
    value_as_float,
    weighted_metric_score,
    weighted_metric_value,
    weighted_score,
)
from .tsvio import write_dict_tsv


def compact_sequence_integrity_fields(model: CodingModel) -> dict[str, object]:
    """Return only sequence fields required by scoring and compact audit rows.

    Detailed start/stop, phase, splice, and ambiguity metrics remain in the
    dedicated sequence-integrity model table.  Keeping the main model row
    compact avoids multiplying dictionary entries across million-model runs.
    """
    metrics = model.sequence_integrity
    if metrics is None:
        return {
            "sequence_integrity_status": "genome_not_configured",
            "sequence_integrity_penalty_score": 0.0,
            "sequence_integrity_flags": ".",
        }
    return {
        "sequence_integrity_status": metrics.status,
        "sequence_integrity_penalty_score": metrics.penalty_score,
        "sequence_integrity_flags": metrics.flags,
    }


def terminal_feature_prf(
    query_features: list[dict[str, object]],
    target_features: list[dict[str, object]],
    wobble_bp: int,
) -> dict[str, Optional[float]]:
    query_matched = 0
    for q in query_features:
        if any(terminal_cds_exons_match(q, t, wobble_bp) for t in target_features):
            query_matched += 1
    target_matched = 0
    for t in target_features:
        if any(terminal_cds_exons_match(q, t, wobble_bp) for q in query_features):
            target_matched += 1
    return prf_from_counts(query_matched, len(query_features), target_matched, len(target_features))


def pairwise_model_metrics(query: CodingModel, target: CodingModel, *, terminal_cds_3prime_wobble_bp: int) -> dict[str, Optional[float] | str]:
    q_introns = [exact_interval_feature_key(query, s, e) for s, e in query.cds_introns]
    t_introns = [exact_interval_feature_key(target, s, e) for s, e in target.cds_introns]
    q_splice = [splice_site_feature_key(query, site) for site in sorted(query.cds_splice_sites)]
    t_splice = [splice_site_feature_key(target, site) for site in sorted(target.cds_splice_sites)]
    q_internal = [exact_interval_feature_key(query, s, e) for s, e in model_internal_cds_exons(query)]
    t_internal = [exact_interval_feature_key(target, s, e) for s, e in model_internal_cds_exons(target)]

    base = interval_prf(query.cds_parts, target.cds_parts)
    intron = exact_feature_prf(q_introns, t_introns)
    splice = exact_feature_prf(q_splice, t_splice)
    internal = exact_feature_prf(q_internal, t_internal)
    terminal = terminal_feature_prf(
        model_terminal_cds_exons_with_side(query),
        model_terminal_cds_exons_with_side(target),
        terminal_cds_3prime_wobble_bp,
    )

    if query.is_single_exon:
        score, components = weighted_metric_score([
            ("cds_base_f1", PAIRWISE_SINGLE_SCORE_WEIGHTS["cds_base_f1"], base["f1"]),
        ])
        score_type = "single_exon"
    else:
        score, components = weighted_metric_score([
            ("cds_intron_f1", PAIRWISE_MULTI_SCORE_WEIGHTS["cds_intron_f1"], intron["f1"]),
            ("cds_splice_site_f1", PAIRWISE_MULTI_SCORE_WEIGHTS["cds_splice_site_f1"], splice["f1"]),
            ("internal_cds_exon_f1", PAIRWISE_MULTI_SCORE_WEIGHTS["internal_cds_exon_f1"], internal["f1"]),
            ("cds_base_f1", PAIRWISE_MULTI_SCORE_WEIGHTS["cds_base_f1"], base["f1"]),
            ("terminal_cds_exon_f1", PAIRWISE_MULTI_SCORE_WEIGHTS["terminal_cds_exon_f1"], terminal["f1"]),
        ])
        score_type = "multi_exon"

    out: dict[str, Optional[float] | str] = {
        "score": score,
        "score_type": score_type,
        "score_components": components,
    }
    out.update(prf_values("cds_base", base))
    out.update(prf_values("cds_intron", intron))
    out.update(prf_values("cds_splice_site", splice))
    out.update(prf_values("internal_cds_exon", internal))
    out.update(prf_values("terminal_cds_exon", terminal))
    return out


def intrinsic_penalty_for_model(
    model: CodingModel,
    *,
    min_cds_length_bp: int,
    min_terminal_cds_exon_bp: int,
    min_internal_cds_exon_bp: int,
    min_cds_intron_bp: int,
    long_cds_intron_bp: Optional[float],
    query_chain_contained_in_other_source_count: int,
    other_chain_contained_in_query_source_count: int,
    contained_model_source_fraction: float,
    unsupported_contained_cds_intron_count: int = 0,
    unsupported_contained_cds_intron_fraction: Optional[float] = None,
    weights: dict[str, float],
) -> tuple[float, str, dict[str, int]]:
    short_cds_count = 1 if model.cds_length < min_cds_length_bp else 0
    mod3_count = 1 if model.cds_length % 3 != 0 else 0
    short_terminal_count = sum(1 for length in model.terminal_cds_exon_lengths if length < min_terminal_cds_exon_bp)
    short_internal_count = sum(1 for length in model.internal_cds_exon_lengths if length < min_internal_cds_exon_bp)
    short_intron_count = sum(1 for length in model.cds_intron_lengths if length < min_cds_intron_bp)
    long_intron_count = (
        sum(1 for length in model.cds_intron_lengths if length > long_cds_intron_bp)
        if long_cds_intron_bp is not None
        else 0
    )
    # Directional contained-model indicator: keep the count as a binary flag for
    # reporting, but make the penalty contribution proportional to the fraction
    # of other sources that contain the query.  The source fraction is now
    # calculated only for models with missing supported terminal CDS structure,
    # which avoids penalising cases where another model simply splits a terminal
    # exon by adding an extra intron.
    contained_model_count = 1 if contained_model_source_fraction > 0 else 0
    counts = {
        "short_cds": short_cds_count,
        "cds_length_not_multiple_of_3": mod3_count,
        "short_terminal_cds_exon": short_terminal_count,
        "short_internal_cds_exon": short_internal_count,
        "short_cds_intron": short_intron_count,
        "long_cds_intron": long_intron_count,
        "contained_model": contained_model_count,
        "unsupported_cds_intron_contained_in_other_cds_exon": int(unsupported_contained_cds_intron_count),
    }
    score = sum(
        weights.get(name, 0.0) * count
        for name, count in counts.items()
        if name not in {"contained_model", "unsupported_cds_intron_contained_in_other_cds_exon"}
    )
    if contained_model_source_fraction > 0:
        score += weights.get("contained_model", 0.0) * contained_model_source_fraction
    if unsupported_contained_cds_intron_count > 0 and unsupported_contained_cds_intron_fraction is not None:
        score += weights.get("unsupported_cds_intron_contained_in_other_cds_exon", 0.0) * unsupported_contained_cds_intron_fraction
    components_list = [
        f"{name}:{count}*{weights.get(name, 0.0):.3f}"
        for name, count in counts.items()
        if count and name not in {"contained_model", "unsupported_cds_intron_contained_in_other_cds_exon"}
    ]
    if unsupported_contained_cds_intron_count > 0 and unsupported_contained_cds_intron_fraction is not None:
        components_list.append(
            f"unsupported_cds_intron_contained_in_other_cds_exon_fraction:{unsupported_contained_cds_intron_fraction:.6f}*{weights.get('unsupported_cds_intron_contained_in_other_cds_exon', 0.0):.3f}"
        )
    if contained_model_source_fraction > 0:
        components_list.append(
            f"contained_model_source_fraction:{contained_model_source_fraction:.6f}*{weights.get('contained_model', 0.0):.3f}"
        )
    components = ",".join(components_list) or "."
    return score, components, counts


def is_contiguous_subchain(short_chain: tuple[tuple[int, int], ...], long_chain: tuple[tuple[int, int], ...]) -> bool:
    if not short_chain or len(short_chain) >= len(long_chain):
        return False
    n = len(short_chain)
    return any(tuple(long_chain[i:i + n]) == short_chain for i in range(0, len(long_chain) - n + 1))


def model_internal_cds_exons(model: CodingModel) -> list[tuple[int, int]]:
    parts = sorted(model.cds_parts)
    if len(parts) <= 2:
        return []
    return parts[1:-1]


def model_terminal_cds_exons_with_side(model: CodingModel) -> list[dict[str, object]]:
    """Return terminal CDS exon features for multi-exon models.

    terminal_side is biological 5_prime or 3_prime. For the 3_prime terminal CDS
    exon, support matching can allow a strand-aware outer-boundary wobble to
    accommodate stop-codon inclusion/exclusion differences.
    """
    parts = sorted(model.cds_parts)
    if len(parts) <= 1:
        return []
    if model.strand == "+":
        items = [(parts[0], "5_prime"), (parts[-1], "3_prime")]
    else:
        items = [(parts[-1], "5_prime"), (parts[0], "3_prime")]
    rows: list[dict[str, object]] = []
    for (start0, end0), side in items:
        rows.append(
            {
                "model": model,
                "model_id": model_key(model),
                "source": model.source_name,
                "chrom": model.chrom,
                "strand": model.strand,
                "start0": start0,
                "end0": end0,
                "terminal_side": side,
            }
        )
    return rows


def terminal_cds_exons_match(a: dict[str, object], b: dict[str, object], wobble_bp: int) -> bool:
    if a["chrom"] != b["chrom"] or a["strand"] != b["strand"]:
        return False
    if a["terminal_side"] != b["terminal_side"]:
        return False
    a_start = int(a["start0"]); a_end = int(a["end0"])
    b_start = int(b["start0"]); b_end = int(b["end0"])
    if a["terminal_side"] != "3_prime" or wobble_bp <= 0:
        return a_start == b_start and a_end == b_end
    strand = str(a["strand"])
    if strand == "+":
        # The internal splice boundary is the low-coordinate boundary; allow the
        # coding 3' outer end to differ by up to wobble_bp.
        return a_start == b_start and abs(a_end - b_end) <= wobble_bp
    if strand == "-":
        # The internal splice boundary is the high-coordinate boundary; allow the
        # coding 3' outer start to differ by up to wobble_bp.
        return a_end == b_end and abs(a_start - b_start) <= wobble_bp
    return a_start == b_start and a_end == b_end


def cds_exon_chains_match(a: CodingModel, b: CodingModel, wobble_bp: int) -> bool:
    """Return True if ordered CDS exon chains match, allowing 3' wobble.

    This is stricter than intron-chain identity because it includes terminal CDS
    exon boundaries.  The configured wobble only applies to the biological 3'
    terminal CDS outer boundary; all other CDS exon boundaries must match
    exactly.
    """
    if a.chrom != b.chrom or a.strand != b.strand:
        return False
    a_parts = sorted(a.cds_parts)
    b_parts = sorted(b.cds_parts)
    if len(a_parts) != len(b_parts):
        return False
    if not a_parts:
        return False
    three_prime_index = len(a_parts) - 1 if a.strand == "+" else 0
    for idx, ((a_start, a_end), (b_start, b_end)) in enumerate(zip(a_parts, b_parts)):
        if idx == three_prime_index and wobble_bp > 0:
            if a.strand == "+":
                # Internal/splice-side boundary exact; outer 3' end may wobble.
                if a_start != b_start or abs(a_end - b_end) > wobble_bp:
                    return False
            elif a.strand == "-":
                # Internal/splice-side boundary exact; outer 3' start may wobble.
                if a_end != b_end or abs(a_start - b_start) > wobble_bp:
                    return False
            else:
                if a_start != b_start or a_end != b_end:
                    return False
        else:
            if a_start != b_start or a_end != b_end:
                return False
    return True


def cds_intron_chains_match(a: CodingModel, b: CodingModel) -> bool:
    """Return true when two multi-exon models have the same CDS intron chain.

    Terminal CDS boundaries are intentionally ignored.  The helper is kept
    separate from :func:`cds_exon_chains_match` because alternative-splice
    analysis treats same-intron-chain terminal variants as redundant even when
    their terminal difference exceeds the configured 3-prime wobble.
    """
    return bool(
        a.chrom == b.chrom
        and a.strand == b.strand
        and a.is_multi_exon
        and b.is_multi_exon
        and tuple(a.cds_introns) == tuple(b.cds_introns)
    )


def cds_exon_chain_is_contiguous_subset(
    subset_model: CodingModel,
    container_model: CodingModel,
) -> bool:
    """Return true for a strict, unchanged contiguous CDS-exon subchain.

    This is deliberately based on complete CDS exon intervals rather than
    genomic-span containment or intron chains alone.  Boundaries must match
    exactly and the subset must contain fewer CDS exons than the container.
    """
    if subset_model.chrom != container_model.chrom or subset_model.strand != container_model.strand:
        return False
    subset = tuple(subset_model.sorted_cds_parts)
    container = tuple(container_model.sorted_cds_parts)
    if not subset or len(subset) >= len(container):
        return False
    width = len(subset)
    return any(
        tuple(container[index : index + width]) == subset
        for index in range(len(container) - width + 1)
    )


def cds_exon_chain_containment_relation(a: CodingModel, b: CodingModel) -> dict[str, bool]:
    """Return both directional strict CDS-chain containment relationships."""
    return {
        "a_contained_in_b": cds_exon_chain_is_contiguous_subset(a, b),
        "b_contained_in_a": cds_exon_chain_is_contiguous_subset(b, a),
    }


def parse_cds_phase(value: object) -> Optional[int]:
    """Parse one GFF CDS phase using MoSelecta's existing 0/1/2 semantics."""
    try:
        phase = int(str(value))
    except (TypeError, ValueError):
        return None
    return phase if phase in {0, 1, 2} else None


def model_cds_phases(model: CodingModel) -> dict[tuple[int, int], Optional[int]]:
    """Return the preferred annotation phase for each normalized CDS part.

    This mirrors the sequence-integrity subsystem: the first valid phase for a
    coordinate replaces an unavailable value, while unavailable phase remains
    explicit rather than being interpreted as zero.
    """
    phases: dict[tuple[int, int], Optional[int]] = {}
    for record in model.child_records:
        if record.feature_type != "CDS":
            continue
        key = (record.start0, record.end0)
        parsed = parse_cds_phase(record.phase)
        if key not in phases or (phases[key] is None and parsed is not None):
            phases[key] = parsed
    return phases


def _model_phase_state(
    model: CodingModel,
) -> tuple[str, dict[tuple[int, int], Optional[int]], int]:
    """Return phase availability/consistency state and mismatch count."""
    phases_by_part = model_cds_phases(model)
    transcript_parts = (
        list(reversed(model.sorted_cds_parts))
        if model.strand == "-"
        else list(model.sorted_cds_parts)
    )
    phases = [phases_by_part.get(part) for part in transcript_parts]
    present_count = sum(phase is not None for phase in phases)
    # When the first phase is absent, retain every possible initial phase rather
    # than silently choosing zero.  Known downstream phases may narrow those
    # possibilities; the chain is inconsistent only if none remains valid.
    possible_first_phases = (
        {phases[0]} if phases and phases[0] is not None else {0, 1, 2}
    )
    valid_first_phases: set[int] = set()
    mismatch_counts: list[int] = []
    for first_phase in possible_first_phases:
        cumulative_length = (
            transcript_parts[0][1] - transcript_parts[0][0] - first_phase
            if transcript_parts
            else 0
        )
        possible_mismatches = 0
        for phase, part in zip(phases[1:], transcript_parts[1:]):
            expected = (3 - (cumulative_length % 3)) % 3
            if phase is not None and phase != expected:
                possible_mismatches += 1
            cumulative_length += part[1] - part[0]
        mismatch_counts.append(possible_mismatches)
        if possible_mismatches == 0:
            valid_first_phases.add(first_phase)
    mismatch_count = min(mismatch_counts, default=0)
    if phases and not valid_first_phases:
        status = "phase_inconsistent"
    elif present_count == 0:
        status = "unavailable_missing_phase"
    elif present_count < len(transcript_parts):
        status = "partially_available"
    else:
        status = "calculated"
    return status, phases_by_part, mismatch_count


def _frame_at_genomic_start(
    part: tuple[int, int],
    phase: int,
    strand: str,
) -> int:
    """Return the annotated codon position at the part's low coordinate."""
    start0, end0 = part
    if strand == "-":
        transcript_offset_at_start = end0 - start0 - 1
        return (transcript_offset_at_start - phase) % 3
    return (-phase) % 3


def cds_frame_concordance(a: CodingModel, b: CodingModel) -> dict[str, object]:
    """Compare annotation-derived coding frames across overlapping CDS bases.

    Required fractions are emitted only when both complete phase chains are
    available and internally consistent.  For partial phase, comparable-base
    diagnostics are still returned without pretending that the full overlap was
    evaluable.  No missing phase is silently replaced by zero.
    """
    overlap_bp = interval_list_overlap_bp(a.sorted_cds_parts, b.sorted_cds_parts)
    status_a, phases_a, mismatch_a = _model_phase_state(a)
    status_b, phases_b, mismatch_b = _model_phase_state(b)
    if "phase_inconsistent" in {status_a, status_b}:
        status = "phase_inconsistent"
    elif "unavailable_missing_phase" in {status_a, status_b}:
        status = "unavailable_missing_phase"
    elif "partially_available" in {status_a, status_b}:
        status = "partially_available"
    else:
        status = "calculated"

    comparable_bp = 0
    same_frame_bp = 0
    if a.chrom == b.chrom and a.strand == b.strand:
        for part_a in a.sorted_cds_parts:
            phase_a = phases_a.get(part_a)
            if phase_a is None:
                continue
            frame_a_start = _frame_at_genomic_start(part_a, phase_a, a.strand)
            for part_b in b.sorted_cds_parts:
                if part_b[1] <= part_a[0]:
                    continue
                if part_b[0] >= part_a[1]:
                    break
                phase_b = phases_b.get(part_b)
                if phase_b is None:
                    continue
                start0 = max(part_a[0], part_b[0])
                end0 = min(part_a[1], part_b[1])
                if end0 <= start0:
                    continue
                length = end0 - start0
                comparable_bp += length
                frame_b_start = _frame_at_genomic_start(part_b, phase_b, b.strand)
                direction = -1 if a.strand == "-" else 1
                frame_a = (frame_a_start + direction * (start0 - part_a[0])) % 3
                frame_b = (frame_b_start + direction * (start0 - part_b[0])) % 3
                if frame_a == frame_b:
                    same_frame_bp += length

    calculated = status == "calculated"
    same_frame_overlap_fraction = (
        same_frame_bp / overlap_bp if calculated and overlap_bp else None
    )
    return {
        "frame_concordance_status": status,
        "primary_phase_status": status_a,
        "candidate_phase_status": status_b,
        "primary_phase_mismatch_count": mismatch_a,
        "candidate_phase_mismatch_count": mismatch_b,
        "frame_overlap_bp": overlap_bp,
        "frame_comparable_overlap_bp": comparable_bp,
        "same_frame_overlap_bp": same_frame_bp,
        "same_frame_comparable_fraction": (
            same_frame_bp / comparable_bp if comparable_bp else None
        ),
        "same_frame_overlap_fraction": same_frame_overlap_fraction,
        "same_frame_primary_fraction": (
            same_frame_bp / a.cds_length if calculated and a.cds_length else None
        ),
        "same_frame_candidate_fraction": (
            same_frame_bp / b.cds_length if calculated and b.cds_length else None
        ),
    }


def exact_interval_feature_key(model: CodingModel, start0: int, end0: int) -> tuple[str, str, int, int]:
    return (model.chrom, model.strand, start0, end0)


def splice_site_feature_key(model: CodingModel, site0: int) -> tuple[str, str, int]:
    return (model.chrom, model.strand, site0)


def make_exact_support_map(
    models: list[CodingModel],
    feature_getter,
) -> dict[tuple, dict[str, set[str]]]:
    support: dict[tuple, dict[str, set[str]]] = defaultdict(lambda: {"models": set(), "sources": set()})
    for model in models:
        for feature in feature_getter(model):
            support[feature]["models"].add(model_key(model))
            support[feature]["sources"].add(model.source_name)
    return support


def summarize_feature_support(
    features: list[object],
    *,
    model: CodingModel,
    support_getter,
) -> dict[str, object]:
    total = len(features)
    supported_by_other_model = 0
    supported_by_other_source = 0
    supported_by_2plus_sources = 0
    support_model_sum = 0
    support_source_sum = 0
    support_model_counts: list[int] = []
    support_source_counts: list[int] = []
    for feature in features:
        model_ids, sources = support_getter(feature)
        model_ids = set(model_ids)
        sources = set(sources)
        support_model_counts.append(len(model_ids))
        support_source_counts.append(len(sources))
        support_model_sum += len(model_ids)
        support_source_sum += len(sources)
        if model_ids - {model_key(model)}:
            supported_by_other_model += 1
        if sources - {model.source_name}:
            supported_by_other_source += 1
        if len(sources) >= 2:
            supported_by_2plus_sources += 1
    unsupported = total - supported_by_other_source
    return {
        "total": total,
        "supported_by_other_model": supported_by_other_model,
        "supported_by_other_source": supported_by_other_source,
        "supported_by_2plus_sources": supported_by_2plus_sources,
        "support_model_sum": support_model_sum,
        "support_model_mean": fraction_or_none(support_model_sum, total),
        "support_model_min": min(support_model_counts) if support_model_counts else None,
        "support_source_sum": support_source_sum,
        "support_source_mean": fraction_or_none(support_source_sum, total),
        "support_source_min": min(support_source_counts) if support_source_counts else None,
        "unsupported_count": unsupported,
        "unsupported_fraction": fraction_or_none(unsupported, total),
        "support_fraction": fraction_or_none(supported_by_other_source, total),
    }


def unsupported_short_cds_exon_penalty_for_model(
    model: CodingModel,
    *,
    internal_exon_support: dict[tuple, dict[str, set[str]]],
    terminal_support_getter,
    min_terminal_cds_exon_bp: int,
    min_internal_cds_exon_bp: int,
    soft_5prime_terminal_cds_exon_bp: int,
    soft_3prime_terminal_cds_exon_bp: int,
    soft_internal_cds_exon_bp: int,
    component_weights: dict[str, float],
) -> dict[str, object]:
    """Score unsupported moderately short CDS exons without double counting intrinsic events.

    The active intrinsic thresholds define the lower edge of each soft interval.
    Exons shorter than those intrinsic thresholds are handled only by the existing
    intrinsic penalty.  The soft penalty applies to exons whose length is at least
    the intrinsic threshold and below the corresponding soft threshold, but only
    when no other annotation source supports the same exon structure.
    """

    def severity(length_bp: int, intrinsic_threshold_bp: int, soft_threshold_bp: int) -> float:
        if soft_threshold_bp <= intrinsic_threshold_bp:
            return 0.0
        if length_bp < intrinsic_threshold_bp or length_bp >= soft_threshold_bp:
            return 0.0
        return max(
            0.0,
            min(
                1.0,
                (soft_threshold_bp - length_bp) / (soft_threshold_bp - intrinsic_threshold_bp),
            ),
        )

    out: dict[str, object] = {
        "unsupported_short_5prime_terminal_cds_exon_count": 0,
        "unsupported_short_5prime_terminal_cds_exon_length_bp": ".",
        "unsupported_short_5prime_terminal_cds_exon_supported_by_other_source": ".",
        "unsupported_short_5prime_terminal_cds_exon_severity": 0.0,
        "unsupported_short_3prime_terminal_cds_exon_count": 0,
        "unsupported_short_3prime_terminal_cds_exon_length_bp": ".",
        "unsupported_short_3prime_terminal_cds_exon_supported_by_other_source": ".",
        "unsupported_short_3prime_terminal_cds_exon_severity": 0.0,
        "moderately_short_internal_cds_exon_count": 0,
        "supported_moderately_short_internal_cds_exon_count": 0,
        "unsupported_moderately_short_internal_cds_exon_count": 0,
        "unsupported_moderately_short_internal_cds_exons": ".",
        "unsupported_short_internal_cds_exon_severity": 0.0,
        "unsupported_short_cds_exon_penalty_score": 0.0,
        "unsupported_short_cds_exon_penalty_components": ".",
        "unsupported_short_terminal_intrinsic_threshold_bp": min_terminal_cds_exon_bp,
        "unsupported_short_internal_intrinsic_threshold_bp": min_internal_cds_exon_bp,
        "unsupported_short_5prime_terminal_soft_threshold_bp": soft_5prime_terminal_cds_exon_bp,
        "unsupported_short_3prime_terminal_soft_threshold_bp": soft_3prime_terminal_cds_exon_bp,
        "unsupported_short_internal_soft_threshold_bp": soft_internal_cds_exon_bp,
        "unsupported_short_5prime_terminal_component_weight": component_weights.get("5prime_terminal", 0.0),
        "unsupported_short_3prime_terminal_component_weight": component_weights.get("3prime_terminal", 0.0),
        "unsupported_short_internal_component_weight": component_weights.get("internal", 0.0),
    }

    # Single-exon CDS models are already covered by the short-CDS and intrinsic
    # terminal-exon metrics.  The new comparative penalty is intentionally
    # limited to terminal exons of multi-exon chains and internal exons of 3+
    # exon chains.
    if not model.is_multi_exon:
        return out

    terminal_severities = {"5_prime": 0.0, "3_prime": 0.0}
    terminal_thresholds = {
        "5_prime": soft_5prime_terminal_cds_exon_bp,
        "3_prime": soft_3prime_terminal_cds_exon_bp,
    }
    terminal_prefixes = {
        "5_prime": "unsupported_short_5prime_terminal_cds_exon",
        "3_prime": "unsupported_short_3prime_terminal_cds_exon",
    }
    for feature in model_terminal_cds_exons_with_side(model):
        side = str(feature["terminal_side"])
        prefix = terminal_prefixes[side]
        length_bp = int(feature["end0"]) - int(feature["start0"])
        exon_severity = severity(
            length_bp,
            min_terminal_cds_exon_bp,
            terminal_thresholds[side],
        )
        if exon_severity <= 0:
            continue
        _model_ids, sources = terminal_support_getter(feature)
        supported_by_other_source = bool(set(sources) - {model.source_name})
        out[f"{prefix}_length_bp"] = length_bp
        out[f"{prefix}_supported_by_other_source"] = "true" if supported_by_other_source else "false"
        if not supported_by_other_source:
            out[f"{prefix}_count"] = 1
            out[f"{prefix}_severity"] = exon_severity
            terminal_severities[side] = exon_severity

    internal_severity_sum = 0.0
    moderately_short_internal_count = 0
    supported_internal_count = 0
    unsupported_internal_count = 0
    unsupported_internal_details: list[str] = []
    for start0, end0 in model_internal_cds_exons(model):
        length_bp = end0 - start0
        exon_severity = severity(
            length_bp,
            min_internal_cds_exon_bp,
            soft_internal_cds_exon_bp,
        )
        if exon_severity <= 0:
            continue
        moderately_short_internal_count += 1
        support_entry = internal_exon_support.get(
            exact_interval_feature_key(model, start0, end0),
            {"models": set(), "sources": set()},
        )
        supported_by_other_source = bool(set(support_entry["sources"]) - {model.source_name})
        if supported_by_other_source:
            supported_internal_count += 1
            continue
        unsupported_internal_count += 1
        internal_severity_sum += exon_severity
        unsupported_internal_details.append(
            f"{model.chrom}:{start0 + 1}-{end0}({length_bp}bp;severity={exon_severity:.6f})"
        )

    internal_severity = min(1.0, internal_severity_sum)
    out["moderately_short_internal_cds_exon_count"] = moderately_short_internal_count
    out["supported_moderately_short_internal_cds_exon_count"] = supported_internal_count
    out["unsupported_moderately_short_internal_cds_exon_count"] = unsupported_internal_count
    out["unsupported_moderately_short_internal_cds_exons"] = (
        ",".join(unsupported_internal_details) if unsupported_internal_details else "."
    )
    out["unsupported_short_internal_cds_exon_severity"] = internal_severity

    component_values = [
        (
            "unsupported_5prime_terminal_severity",
            terminal_severities["5_prime"],
            component_weights.get("5prime_terminal", 0.0),
        ),
        (
            "unsupported_3prime_terminal_severity",
            terminal_severities["3_prime"],
            component_weights.get("3prime_terminal", 0.0),
        ),
        (
            "unsupported_internal_severity",
            internal_severity,
            component_weights.get("internal", 0.0),
        ),
    ]
    penalty_score = min(1.0, sum(value * weight for _name, value, weight in component_values))
    out["unsupported_short_cds_exon_penalty_score"] = penalty_score
    out["unsupported_short_cds_exon_penalty_components"] = (
        ",".join(
            f"{name}:{value:.6f}*{weight:.3f}"
            for name, value, weight in component_values
            if value > 0 and weight > 0
        )
        or "."
    )
    return out


def unsupported_cds_introns_contained_in_other_cds_exons(
    model: CodingModel,
    models: list[CodingModel],
    intron_support: dict[tuple, dict[str, set[str]]],
) -> dict[str, object]:
    """Summarise unsupported query CDS introns contained in another model's CDS exon.

    The context is intentionally limited to the current conservative sublocus.
    A query intron is counted only if it is not exactly supported by any other
    source and the full intron interval lies within at least one CDS exon from
    another model in the same sublocus.  The count is per query intron, not per
    containing model/exon.
    """
    total_introns = len(model.cds_introns)
    contained_count = 0
    containing_model_ids: set[str] = set()
    containing_sources: set[str] = set()
    contained_intron_coords: list[str] = []
    mk = model_key(model)

    if total_introns == 0 or len(models) <= 1:
        return {
            "count": 0,
            "fraction": fraction_or_none(0, total_introns),
            "other_model_count": 0,
            "other_source_count": 0,
            "other_sources": ".",
            "other_model_ids": ".",
            "introns": ".",
        }

    other_exons: list[tuple[int, int, CodingModel]] = []
    for other in models:
        if model_key(other) == mk:
            continue
        for exon_start0, exon_end0 in other.cds_parts:
            other_exons.append((exon_start0, exon_end0, other))

    for intron_start0, intron_end0 in model.cds_introns:
        key = exact_interval_feature_key(model, intron_start0, intron_end0)
        support_entry = intron_support.get(key, {"models": set(), "sources": set()})
        # Do not penalise an intron that has exact support from another source.
        if set(support_entry["sources"]) - {model.source_name}:
            continue

        hits: list[CodingModel] = []
        for exon_start0, exon_end0, other in other_exons:
            if exon_start0 <= intron_start0 and intron_end0 <= exon_end0:
                hits.append(other)

        if hits:
            contained_count += 1
            contained_intron_coords.append(f"{model.chrom}:{intron_start0 + 1}-{intron_end0}")
            for other in hits:
                containing_model_ids.add(model_key(other))
                containing_sources.add(other.source_name)

    return {
        "count": contained_count,
        "fraction": fraction_or_none(contained_count, total_introns),
        "other_model_count": len(containing_model_ids),
        "other_source_count": len(containing_sources),
        "other_sources": fmt_sources(containing_sources),
        "other_model_ids": fmt_model_ids(containing_model_ids),
        "introns": ",".join(contained_intron_coords) if contained_intron_coords else ".",
    }


def gff_like_start_end_for_feature(feature_type: str, start0: int, end0: Optional[int]) -> tuple[str, str]:
    if feature_type == "cds_splice_site" or end0 is None:
        # Splice sites are boundaries rather than intervals. Keep a boundary-like
        # coordinate while also writing start0/end0 columns for unambiguous use.
        return str(start0), str(start0)
    return str(start0 + 1), str(end0)


def calculate_structural_support(
    subloci: list[Sublocus],
    *,
    terminal_cds_3prime_wobble_bp: int,
    single_exon_support_min_reciprocal_overlap: float,
    min_cds_length_bp: int,
    min_terminal_cds_exon_bp: int,
    min_internal_cds_exon_bp: int,
    unsupported_short_5prime_terminal_cds_exon_bp: int,
    unsupported_short_3prime_terminal_cds_exon_bp: int,
    unsupported_short_internal_cds_exon_bp: int,
    unsupported_short_cds_exon_component_weights: dict[str, float],
    min_cds_intron_bp: int,
    long_cds_intron_bp: Optional[float],
    intrinsic_penalty_weights: dict[str, float],
    run_source_count: int,
    structural_corroboration_enabled: bool = False,
    structural_corroboration_rows: Optional[list[dict[str, object]]] = None,
    progress: Optional[ProgressLogger] = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Calculate structural support, similarity, and penalty metrics.

    This includes the original feature-support metrics plus model-vs-consensus
    SN/PR/F1, source-collapsed pairwise SN/PR/F1 summaries, partial/extra
    structure penalties, and an intrinsic model-characteristic penalty score.
    """
    feature_rows: list[dict[str, object]] = []
    model_rows: list[dict[str, object]] = []
    run_possible_other_source_count = max(0, run_source_count - 1)
    sorted_subloci = sorted(subloci, key=lambda s: (s.chrom, s.start0, s.end0, s.sublocus_id))
    total_subloci = len(sorted_subloci)
    total_models = sum(len(sub.models) for sub in sorted_subloci)
    processed_models = 0
    if progress is not None:
        progress.log(
            f"Model scoring structural support: starting {total_subloci:,} subloci / {total_models:,} models",
            level="basic",
        )

    for sub_index, sub in enumerate(sorted_subloci, start=1):
        if progress is not None:
            progress.every(
                sub_index,
                total_subloci,
                "Model scoring structural support subloci",
                level="basic",
            )
        models = sorted(sub.models, key=lambda m: (m.cds_start0, m.cds_end0, m.source_name, m.input_transcript_id))
        sub_sources = {m.source_name for m in models}

        intron_support = make_exact_support_map(
            models,
            lambda m: [exact_interval_feature_key(m, s, e) for s, e in m.cds_introns],
        )
        splice_support = make_exact_support_map(
            models,
            lambda m: [splice_site_feature_key(m, site) for site in m.cds_splice_sites],
        )
        internal_exon_support = make_exact_support_map(
            models,
            lambda m: [exact_interval_feature_key(m, s, e) for s, e in model_internal_cds_exons(m)],
        )
        single_exon_exact_support = make_exact_support_map(
            models,
            lambda m: [exact_interval_feature_key(m, m.cds_start0, m.cds_end0)] if m.is_single_exon else [],
        )

        terminal_features: list[dict[str, object]] = []
        single_features: list[dict[str, object]] = []
        model_terminal_features: dict[str, list[dict[str, object]]] = defaultdict(list)
        model_single_features: dict[str, list[dict[str, object]]] = defaultdict(list)
        source_cds_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for model in models:
            source_cds_intervals[model.source_name].extend(model.cds_parts)
            for item in model_terminal_cds_exons_with_side(model):
                terminal_features.append(item)
                model_terminal_features[model_key(model)].append(item)
            if model.is_single_exon:
                item = {
                    "model": model,
                    "model_id": model_key(model),
                    "source": model.source_name,
                    "chrom": model.chrom,
                    "strand": model.strand,
                    "start0": model.cds_start0,
                    "end0": model.cds_end0,
                }
                single_features.append(item)
                model_single_features[model_key(model)].append(item)

        def terminal_support_sets(query_feature: dict[str, object]) -> tuple[set[str], set[str]]:
            model_ids: set[str] = set()
            sources: set[str] = set()
            for other in terminal_features:
                if terminal_cds_exons_match(query_feature, other, terminal_cds_3prime_wobble_bp):
                    model_ids.add(str(other["model_id"]))
                    sources.add(str(other["source"]))
            return model_ids, sources

        def single_overlap_support_sets(query_feature: dict[str, object]) -> tuple[set[str], set[str]]:
            model_ids: set[str] = set()
            sources: set[str] = set()
            q_start = int(query_feature["start0"]); q_end = int(query_feature["end0"])
            for other in single_features:
                if query_feature["chrom"] != other["chrom"] or query_feature["strand"] != other["strand"]:
                    continue
                overlap, ro_q, ro_o = reciprocal_overlap(q_start, q_end, int(other["start0"]), int(other["end0"]))
                if overlap > 0 and min(ro_q, ro_o) >= single_exon_support_min_reciprocal_overlap:
                    model_ids.add(str(other["model_id"]))
                    sources.add(str(other["source"]))
            return model_ids, sources

        # Feature-level support rows. Exact features are collapsed by exact key;
        # terminal 3' rows can include support from wobble-equivalent features.
        for feature_type, support_map in [
            ("cds_intron", intron_support),
            ("cds_splice_site", splice_support),
            ("internal_cds_exon", internal_exon_support),
            ("single_exon_cds_exact", single_exon_exact_support),
        ]:
            for key, support in sorted(support_map.items(), key=lambda kv: kv[0]):
                if feature_type == "cds_splice_site":
                    chrom, strand, site0 = key
                    start0 = site0; end0 = None
                else:
                    chrom, strand, start0, end0 = key
                start_text, end_text = gff_like_start_end_for_feature(feature_type, int(start0), None if end0 is None else int(end0))
                feature_rows.append(
                    {
                        "superlocus_id": sub.superlocus_id,
                        "sublocus_id": sub.sublocus_id,
                        "sublocus_type": sub.sublocus_type,
                        "feature_type": feature_type,
                        "matching_mode": "exact",
                        "terminal_side": ".",
                        "chrom": chrom,
                        "start": start_text,
                        "end": end_text,
                        "strand": strand,
                        "start0": start0,
                        "end0": "." if end0 is None else end0,
                        "support_model_count": len(support["models"]),
                        "support_source_count": len(support["sources"]),
                        "support_sources": fmt_sources(support["sources"]),
                        "support_model_ids": fmt_model_ids(support["models"]),
                    }
                )

        seen_terminal_feature_keys: set[tuple] = set()
        terminal_consensus_features: list[dict[str, object]] = []
        for item in sorted(terminal_features, key=lambda x: (str(x["chrom"]), int(x["start0"]), int(x["end0"]), str(x["strand"]), str(x["terminal_side"]))):
            key = (item["chrom"], item["strand"], item["start0"], item["end0"], item["terminal_side"])
            if key in seen_terminal_feature_keys:
                continue
            seen_terminal_feature_keys.add(key)
            model_ids, sources = terminal_support_sets(item)
            if len(sources) >= 2:
                terminal_consensus_features.append(item)
            start_text, end_text = gff_like_start_end_for_feature("terminal_cds_exon", int(item["start0"]), int(item["end0"]))
            matching_mode = "3prime_wobble" if item["terminal_side"] == "3_prime" and terminal_cds_3prime_wobble_bp > 0 else "exact"
            feature_rows.append(
                {
                    "superlocus_id": sub.superlocus_id,
                    "sublocus_id": sub.sublocus_id,
                    "sublocus_type": sub.sublocus_type,
                    "feature_type": "terminal_cds_exon",
                    "matching_mode": matching_mode,
                    "terminal_side": item["terminal_side"],
                    "chrom": item["chrom"],
                    "start": start_text,
                    "end": end_text,
                    "strand": item["strand"],
                    "start0": item["start0"],
                    "end0": item["end0"],
                    "support_model_count": len(model_ids),
                    "support_source_count": len(sources),
                    "support_sources": fmt_sources(sources),
                    "support_model_ids": fmt_model_ids(model_ids),
                }
            )

        seen_single_overlap_keys: set[tuple] = set()
        for item in sorted(single_features, key=lambda x: (str(x["chrom"]), int(x["start0"]), int(x["end0"]), str(x["strand"]))):
            key = (item["chrom"], item["strand"], item["start0"], item["end0"])
            if key in seen_single_overlap_keys:
                continue
            seen_single_overlap_keys.add(key)
            model_ids, sources = single_overlap_support_sets(item)
            start_text, end_text = gff_like_start_end_for_feature("single_exon_cds_overlap", int(item["start0"]), int(item["end0"]))
            feature_rows.append(
                {
                    "superlocus_id": sub.superlocus_id,
                    "sublocus_id": sub.sublocus_id,
                    "sublocus_type": sub.sublocus_type,
                    "feature_type": "single_exon_cds_overlap",
                    "matching_mode": f"reciprocal_overlap>={single_exon_support_min_reciprocal_overlap:.3f}",
                    "terminal_side": ".",
                    "chrom": item["chrom"],
                    "start": start_text,
                    "end": end_text,
                    "strand": item["strand"],
                    "start0": item["start0"],
                    "end0": item["end0"],
                    "support_model_count": len(model_ids),
                    "support_source_count": len(sources),
                    "support_sources": fmt_sources(sources),
                    "support_model_ids": fmt_model_ids(model_ids),
                }
            )

        # Sublicus consensus feature sets: source-supported by >=2 distinct sources.
        consensus_cds_intervals = consensus_intervals_from_source_intervals(source_cds_intervals, min_source_support=2)
        consensus_introns = {key for key, support in intron_support.items() if len(support["sources"]) >= 2}
        consensus_splice_sites = {key for key, support in splice_support.items() if len(support["sources"]) >= 2}
        consensus_internal_exons = {key for key, support in internal_exon_support.items() if len(support["sources"]) >= 2}

        chain_to_models: dict[tuple[tuple[int, int], ...], list[CodingModel]] = defaultdict(list)
        for model in models:
            chain = tuple(model.cds_introns)
            if chain:
                chain_to_models[chain].append(model)

        for model in models:
            mk = model_key(model)
            possible_other_source_count = max(0, len(sub_sources - {model.source_name}))

            intron_features = [exact_interval_feature_key(model, s, e) for s, e in model.cds_introns]
            splice_features = [splice_site_feature_key(model, site) for site in sorted(model.cds_splice_sites)]
            internal_exon_features = [exact_interval_feature_key(model, s, e) for s, e in model_internal_cds_exons(model)]
            terminal_model_features = model_terminal_features.get(mk, [])
            single_model_features = model_single_features.get(mk, [])

            intron_summary = summarize_feature_support(
                intron_features,
                model=model,
                support_getter=lambda feature: (intron_support[feature]["models"], intron_support[feature]["sources"]),
            )
            splice_summary = summarize_feature_support(
                splice_features,
                model=model,
                support_getter=lambda feature: (splice_support[feature]["models"], splice_support[feature]["sources"]),
            )
            internal_exon_summary = summarize_feature_support(
                internal_exon_features,
                model=model,
                support_getter=lambda feature: (internal_exon_support[feature]["models"], internal_exon_support[feature]["sources"]),
            )
            terminal_summary = summarize_feature_support(
                terminal_model_features,
                model=model,
                support_getter=terminal_support_sets,
            )
            unsupported_short_exon_summary = unsupported_short_cds_exon_penalty_for_model(
                model,
                internal_exon_support=internal_exon_support,
                terminal_support_getter=terminal_support_sets,
                min_terminal_cds_exon_bp=min_terminal_cds_exon_bp,
                min_internal_cds_exon_bp=min_internal_cds_exon_bp,
                soft_5prime_terminal_cds_exon_bp=unsupported_short_5prime_terminal_cds_exon_bp,
                soft_3prime_terminal_cds_exon_bp=unsupported_short_3prime_terminal_cds_exon_bp,
                soft_internal_cds_exon_bp=unsupported_short_internal_cds_exon_bp,
                component_weights=unsupported_short_cds_exon_component_weights,
            )

            best_single_overlap: Optional[float] = None
            best_single_overlap_model = "."
            best_single_overlap_source = "."
            if model.is_single_exon:
                for other in single_features:
                    if str(other["model_id"]) == mk:
                        continue
                    if str(other["source"]) == model.source_name:
                        continue
                    overlap, ro_q, ro_o = reciprocal_overlap(model.cds_start0, model.cds_end0, int(other["start0"]), int(other["end0"]))
                    if overlap <= 0:
                        continue
                    score = min(ro_q, ro_o)
                    if best_single_overlap is None or score > best_single_overlap:
                        best_single_overlap = score
                        best_single_overlap_model = str(other["model_id"])
                        best_single_overlap_source = str(other["source"])

            chain = tuple(model.cds_introns)
            same_models = [m for m in chain_to_models.get(chain, []) if model_key(m) != mk] if chain else []
            same_other_sources = {m.source_name for m in same_models if m.source_name != model.source_name}
            same_other_model_ids = {model_key(m) for m in same_models}
            exact_chain_support_fraction = None if not chain else source_scaled_fraction(len(same_other_sources), possible_other_source_count)

            exact_cds_exon_chain_support_models = [
                other for other in models
                if model_key(other) != mk and cds_exon_chains_match(model, other, terminal_cds_3prime_wobble_bp)
            ]
            exact_cds_exon_chain_support_model_ids = {model_key(m) for m in exact_cds_exon_chain_support_models}
            exact_cds_exon_chain_support_sources = {
                m.source_name for m in exact_cds_exon_chain_support_models
                if m.source_name != model.source_name
            }
            exact_cds_exon_chain_support_source_count = len(exact_cds_exon_chain_support_sources)
            exact_cds_exon_chain_support_fraction_local = source_scaled_fraction(
                exact_cds_exon_chain_support_source_count,
                possible_other_source_count,
            )
            exact_structure_support_score = source_scaled_fraction(
                exact_cds_exon_chain_support_source_count,
                run_possible_other_source_count,
            )

            query_contained_models: list[CodingModel] = []
            other_contained_models: list[CodingModel] = []
            if chain:
                for other in models:
                    if model_key(other) == mk or not other.cds_introns:
                        continue
                    other_chain = tuple(other.cds_introns)
                    if is_contiguous_subchain(chain, other_chain):
                        query_contained_models.append(other)
                    if is_contiguous_subchain(other_chain, chain):
                        other_contained_models.append(other)
            query_contained_sources = {m.source_name for m in query_contained_models if m.source_name != model.source_name}
            other_contained_sources = {m.source_name for m in other_contained_models if m.source_name != model.source_name}

            if model.is_multi_exon:
                structural_score, component_text = weighted_score(
                    [
                        ("intron", STRUCTURAL_MULTI_SCORE_WEIGHTS["cds_intron_support_fraction"], intron_summary["support_fraction"]),
                        ("splice_site", STRUCTURAL_MULTI_SCORE_WEIGHTS["cds_splice_site_support_fraction"], splice_summary["support_fraction"]),
                        ("internal_exon", STRUCTURAL_MULTI_SCORE_WEIGHTS["internal_cds_exon_support_fraction"], internal_exon_summary["support_fraction"]),
                        ("terminal_exon", STRUCTURAL_MULTI_SCORE_WEIGHTS["terminal_cds_exon_support_fraction"], terminal_summary["support_fraction"]),
                        ("exact_chain", STRUCTURAL_MULTI_SCORE_WEIGHTS["exact_intron_chain_support_fraction"], exact_chain_support_fraction),
                    ]
                )
                score_type = "multi_exon"
            else:
                single_exact_support_fraction = source_scaled_fraction(
                    len(set(single_exon_exact_support.get(exact_interval_feature_key(model, model.cds_start0, model.cds_end0), {"sources": set()})["sources"]) - {model.source_name}),
                    possible_other_source_count,
                )
                overlap_sources = set()
                for feature in single_model_features:
                    _model_ids, sources = single_overlap_support_sets(feature)
                    overlap_sources.update(set(sources) - {model.source_name})
                single_overlap_support_fraction = source_scaled_fraction(len(overlap_sources), possible_other_source_count)
                structural_score, component_text = weighted_score(
                    [
                        ("single_overlap", STRUCTURAL_SINGLE_SCORE_WEIGHTS["single_exon_overlap_support_fraction"], single_overlap_support_fraction),
                        ("single_exact", STRUCTURAL_SINGLE_SCORE_WEIGHTS["single_exon_exact_support_fraction"], single_exact_support_fraction),
                    ]
                )
                score_type = "single_exon"

            # Source-scaled single-exon fractions are useful even though raw
            # summaries above are also reported.
            if model.is_single_exon:
                single_exact_key = exact_interval_feature_key(model, model.cds_start0, model.cds_end0)
                single_exact_support_entry = single_exon_exact_support.get(single_exact_key, {"models": set(), "sources": set()})
                exact_model_ids = set(single_exact_support_entry["models"]) - {mk}
                exact_sources = set(single_exact_support_entry["sources"]) - {model.source_name}
                overlap_sources = set()
                overlap_model_ids = set()
                for feature in single_model_features:
                    model_ids, sources = single_overlap_support_sets(feature)
                    overlap_model_ids.update(set(model_ids) - {mk})
                    overlap_sources.update(set(sources) - {model.source_name})
                single_exact_supported_feature_count = 1 if exact_model_ids else 0
                single_exact_support_model_count = len(exact_model_ids)
                single_exact_support_source_count = len(exact_sources)
                single_overlap_support_source_count = len(overlap_sources)
                single_exact_support_fraction = source_scaled_fraction(single_exact_support_source_count, possible_other_source_count)
                single_overlap_support_fraction = source_scaled_fraction(single_overlap_support_source_count, possible_other_source_count)
                single_exact_support_sources = fmt_sources(exact_sources)
                single_overlap_support_sources = fmt_sources(overlap_sources)
                single_overlap_support_model_count = len(overlap_model_ids)
            else:
                single_exact_supported_feature_count = "."
                single_exact_support_model_count = "."
                single_exact_support_source_count = "."
                single_overlap_support_source_count = "."
                single_exact_support_fraction = None
                single_overlap_support_fraction = None
                single_exact_support_sources = "."
                single_overlap_support_sources = "."
                single_overlap_support_model_count = "."

            # Model-vs-consensus SN/PR/F1.
            consensus_base = interval_prf(model.cds_parts, consensus_cds_intervals)
            consensus_intron = exact_feature_prf(intron_features, consensus_introns)
            consensus_splice = exact_feature_prf(splice_features, consensus_splice_sites)
            consensus_internal = exact_feature_prf(internal_exon_features, consensus_internal_exons)
            consensus_terminal = terminal_feature_prf(terminal_model_features, terminal_consensus_features, terminal_cds_3prime_wobble_bp)

            if model.is_single_exon:
                consensus_similarity_score, consensus_components = weighted_metric_score([
                    ("cds_base_consensus_f1", CONSENSUS_SINGLE_SCORE_WEIGHTS["cds_base_consensus_f1"], consensus_base["f1"]),
                    ("single_exon_overlap_support_fraction", CONSENSUS_SINGLE_SCORE_WEIGHTS["single_exon_overlap_support_fraction"], single_overlap_support_fraction),
                ])
                consensus_similarity_score_type = "single_exon"
                consensus_precision = weighted_metric_value([
                    ("cds_base_consensus_precision", CONSENSUS_SINGLE_SCORE_WEIGHTS["cds_base_consensus_f1"], consensus_base["precision"]),
                    ("single_exon_overlap_support_fraction", CONSENSUS_SINGLE_SCORE_WEIGHTS["single_exon_overlap_support_fraction"], single_overlap_support_fraction),
                ])
                consensus_sensitivity = weighted_metric_value([
                    ("cds_base_consensus_sensitivity", CONSENSUS_SINGLE_SCORE_WEIGHTS["cds_base_consensus_f1"], consensus_base["sensitivity"]),
                    ("single_exon_overlap_support_fraction", CONSENSUS_SINGLE_SCORE_WEIGHTS["single_exon_overlap_support_fraction"], single_overlap_support_fraction),
                ])
            else:
                consensus_similarity_score, consensus_components = weighted_metric_score([
                    ("cds_intron_consensus_f1", CONSENSUS_MULTI_SCORE_WEIGHTS["cds_intron_consensus_f1"], consensus_intron["f1"]),
                    ("cds_splice_site_consensus_f1", CONSENSUS_MULTI_SCORE_WEIGHTS["cds_splice_site_consensus_f1"], consensus_splice["f1"]),
                    ("internal_cds_exon_consensus_f1", CONSENSUS_MULTI_SCORE_WEIGHTS["internal_cds_exon_consensus_f1"], consensus_internal["f1"]),
                    ("cds_base_consensus_f1", CONSENSUS_MULTI_SCORE_WEIGHTS["cds_base_consensus_f1"], consensus_base["f1"]),
                    ("terminal_cds_exon_consensus_f1", CONSENSUS_MULTI_SCORE_WEIGHTS["terminal_cds_exon_consensus_f1"], consensus_terminal["f1"]),
                ])
                consensus_similarity_score_type = "multi_exon"
                consensus_precision = weighted_metric_value([
                    ("cds_intron_consensus_precision", CONSENSUS_MULTI_SCORE_WEIGHTS["cds_intron_consensus_f1"], consensus_intron["precision"]),
                    ("cds_splice_site_consensus_precision", CONSENSUS_MULTI_SCORE_WEIGHTS["cds_splice_site_consensus_f1"], consensus_splice["precision"]),
                    ("internal_cds_exon_consensus_precision", CONSENSUS_MULTI_SCORE_WEIGHTS["internal_cds_exon_consensus_f1"], consensus_internal["precision"]),
                    ("cds_base_consensus_precision", CONSENSUS_MULTI_SCORE_WEIGHTS["cds_base_consensus_f1"], consensus_base["precision"]),
                    ("terminal_cds_exon_consensus_precision", CONSENSUS_MULTI_SCORE_WEIGHTS["terminal_cds_exon_consensus_f1"], consensus_terminal["precision"]),
                ])
                consensus_sensitivity = weighted_metric_value([
                    ("cds_intron_consensus_sensitivity", CONSENSUS_MULTI_SCORE_WEIGHTS["cds_intron_consensus_f1"], consensus_intron["sensitivity"]),
                    ("cds_splice_site_consensus_sensitivity", CONSENSUS_MULTI_SCORE_WEIGHTS["cds_splice_site_consensus_f1"], consensus_splice["sensitivity"]),
                    ("internal_cds_exon_consensus_sensitivity", CONSENSUS_MULTI_SCORE_WEIGHTS["internal_cds_exon_consensus_f1"], consensus_internal["sensitivity"]),
                    ("cds_base_consensus_sensitivity", CONSENSUS_MULTI_SCORE_WEIGHTS["cds_base_consensus_f1"], consensus_base["sensitivity"]),
                    ("terminal_cds_exon_consensus_sensitivity", CONSENSUS_MULTI_SCORE_WEIGHTS["terminal_cds_exon_consensus_f1"], consensus_terminal["sensitivity"]),
                ])

            if consensus_precision is None or consensus_sensitivity is None:
                partial_structure_penalty = None
                extra_structure_penalty = None
            else:
                partial_structure_penalty = max(0.0, consensus_precision - consensus_sensitivity)
                extra_structure_penalty = max(0.0, consensus_sensitivity - consensus_precision)

            # Completeness/containment diagnostics. These are base-level and
            # intentionally separate from exact-feature consensus scores. They
            # distinguish models that recover all supported CDS bases but have
            # terminal extensions from models that are precise internal cores but
            # miss supported terminal structure.
            consensus_cds_bp = interval_total_bp(consensus_cds_intervals)
            query_cds_bp = interval_total_bp(model.cds_parts)
            consensus_cds_bp_recovered = interval_list_overlap_bp(model.cds_parts, consensus_cds_intervals)
            missing_supported_intervals = subtract_intervals(consensus_cds_intervals, model.cds_parts)
            extra_query_intervals = subtract_intervals(model.cds_parts, consensus_cds_intervals)
            missing_supported_cds_bp = interval_total_bp(missing_supported_intervals)
            extra_cds_bp = interval_total_bp(extra_query_intervals)
            if consensus_cds_intervals:
                consensus_span_start = min(s for s, _e in consensus_cds_intervals)
                consensus_span_end = max(e for _s, e in consensus_cds_intervals)
                extra_before_bp, internal_extra_bp, extra_after_bp, _extra_total_check = classify_intervals_relative_to_span(
                    extra_query_intervals, consensus_span_start, consensus_span_end
                )
                terminal_extra_5prime_bp, terminal_extra_3prime_bp = five_three_from_before_after(model.strand, extra_before_bp, extra_after_bp)
                terminal_extra_bp = terminal_extra_5prime_bp + terminal_extra_3prime_bp
            else:
                consensus_span_start = None
                consensus_span_end = None
                terminal_extra_bp = 0
                terminal_extra_5prime_bp = 0
                terminal_extra_3prime_bp = 0
                internal_extra_bp = 0

            # Missing supported bases are classified relative to the query CDS
            # span: outside the query span is terminal-like missing structure,
            # whereas missing supported bases within the query span indicate an
            # internal gap/structural disagreement.
            missing_before_query_bp, missing_internal_bp, missing_after_query_bp, _missing_total_check = classify_intervals_relative_to_span(
                missing_supported_intervals, model.cds_start0, model.cds_end0
            )
            missing_supported_terminal_5prime_bp, missing_supported_terminal_3prime_bp = five_three_from_before_after(
                model.strand, missing_before_query_bp, missing_after_query_bp
            )
            missing_supported_terminal_structure_bp = missing_supported_terminal_5prime_bp + missing_supported_terminal_3prime_bp
            missing_supported_internal_structure_bp = missing_internal_bp

            missing_supported_cds_fraction = fraction_or_none(missing_supported_cds_bp, consensus_cds_bp)
            missing_supported_terminal_structure_fraction = fraction_or_none(missing_supported_terminal_structure_bp, consensus_cds_bp)
            missing_supported_internal_structure_fraction = fraction_or_none(missing_supported_internal_structure_bp, consensus_cds_bp)
            extra_cds_fraction = fraction_or_none(extra_cds_bp, query_cds_bp)
            terminal_extra_fraction = fraction_or_none(terminal_extra_bp, query_cds_bp)
            internal_extra_fraction = fraction_or_none(internal_extra_bp, query_cds_bp)
            terminal_extension_score = fraction_or_none(terminal_extra_bp, extra_cds_bp)
            contained_model_source_fraction = (
                source_scaled_fraction(len(query_contained_sources), possible_other_source_count)
                if query_contained_sources and not other_contained_sources and missing_supported_terminal_structure_bp > 0
                else 0.0
            )

            if consensus_cds_bp <= 0:
                model_completeness_penalty = None
                model_completeness_score = None
                model_completeness_score_components = "."
            else:
                completeness_components = [
                    (
                        "missing_supported_terminal_structure_fraction",
                        COMPLETENESS_PENALTY_WEIGHTS["missing_supported_terminal_structure_fraction"],
                        missing_supported_terminal_structure_fraction,
                    ),
                    (
                        "missing_supported_internal_structure_fraction",
                        COMPLETENESS_PENALTY_WEIGHTS["missing_supported_internal_structure_fraction"],
                        missing_supported_internal_structure_fraction,
                    ),
                    (
                        "contained_model_source_fraction",
                        COMPLETENESS_PENALTY_WEIGHTS["contained_model_source_fraction"],
                        contained_model_source_fraction,
                    ),
                ]
                available_completeness = [(name, weight, value) for name, weight, value in completeness_components if value is not None]
                model_completeness_penalty = min(1.0, sum(weight * value for _name, weight, value in available_completeness))
                model_completeness_score = max(0.0, 1.0 - model_completeness_penalty)
                model_completeness_score_components = ",".join(
                    f"{name}:{value:.6f}*{weight:.3f}"
                    for name, weight, value in available_completeness
                ) or "."

            # Source-collapsed pairwise summaries. For each other source, retain
            # the best matching model by combined pairwise similarity score, then
            # average those best hits across sources.
            best_by_source: dict[str, tuple[CodingModel, dict[str, Optional[float] | str]]] = {}
            corroboration_comparisons = []
            calculate_corroboration = (
                structural_corroboration_enabled
                or structural_corroboration_rows is not None
            )
            best_overall_model = "."
            best_overall_source = "."
            best_overall_score: Optional[float] = None
            for other in models:
                if model_key(other) == mk or other.source_name == model.source_name:
                    continue
                metrics = pairwise_model_metrics(model, other, terminal_cds_3prime_wobble_bp=terminal_cds_3prime_wobble_bp)
                if calculate_corroboration:
                    corroboration_comparisons.append(
                        (
                            other,
                            compare_structural_corroboration(
                                model,
                                other,
                                pairwise_metrics=metrics,
                            ),
                        )
                    )
                metric_score = metrics.get("score")
                score_value = float(metric_score) if metric_score is not None else -1.0
                current = best_by_source.get(other.source_name)
                current_score = float(current[1].get("score")) if current and current[1].get("score") is not None else -1.0
                if current is None or score_value > current_score:
                    best_by_source[other.source_name] = (other, metrics)
                if best_overall_score is None or score_value > best_overall_score:
                    best_overall_score = score_value if metric_score is not None else None
                    best_overall_model = model_key(other)
                    best_overall_source = other.source_name

            best_metrics = [metrics for _other, metrics in best_by_source.values()]
            pairwise_source_count = len(best_metrics)
            pairwise_similarity_score = safe_mean([m.get("score") if isinstance(m.get("score"), float) else None for m in best_metrics])
            pairwise_components = ";".join(
                f"{source}:{metrics.get('score'):.6f}" for source, (_other, metrics) in sorted(best_by_source.items()) if isinstance(metrics.get("score"), float)
            ) or "."

            def mean_best(field: str) -> Optional[float]:
                vals: list[Optional[float]] = []
                for metrics in best_metrics:
                    value = metrics.get(field)
                    vals.append(value if isinstance(value, float) else None)
                return safe_mean(vals)

            best_pairwise_metrics: dict[str, object] = {}
            if best_overall_source != "." and best_overall_source in best_by_source:
                best_pairwise_metrics = best_by_source[best_overall_source][1]
            else:
                best_pairwise_metrics = {}

            unsupported_contained_intron_summary = unsupported_cds_introns_contained_in_other_cds_exons(
                model,
                models,
                intron_support,
            )

            intrinsic_penalty_score, intrinsic_penalty_components, intrinsic_counts = intrinsic_penalty_for_model(
                model,
                min_cds_length_bp=min_cds_length_bp,
                min_terminal_cds_exon_bp=min_terminal_cds_exon_bp,
                min_internal_cds_exon_bp=min_internal_cds_exon_bp,
                min_cds_intron_bp=min_cds_intron_bp,
                long_cds_intron_bp=long_cds_intron_bp,
                query_chain_contained_in_other_source_count=len(query_contained_sources),
                other_chain_contained_in_query_source_count=len(other_contained_sources),
                contained_model_source_fraction=contained_model_source_fraction,
                unsupported_contained_cds_intron_count=int(unsupported_contained_intron_summary["count"]),
                unsupported_contained_cds_intron_fraction=unsupported_contained_intron_summary["fraction"] if isinstance(unsupported_contained_intron_summary["fraction"], float) else None,
                weights=intrinsic_penalty_weights,
            )

            row = {
                "superlocus_id": model.superlocus_id or ".",
                "sublocus_id": model.sublocus_id or ".",
                "sublocus_type": sub.sublocus_type,
                "source_name": model.source_name,
                "input_file": str(model.input_file),
                "input_gene_id": model.input_gene_id,
                "input_transcript_id": model.input_transcript_id,
                "output_gene_id": model.output_gene_id,
                "output_transcript_id": model.output_transcript_id,
                "chrom": model.chrom,
                "cds_start": model.cds_start0 + 1,
                "cds_end": model.cds_end0,
                "cds_location": f"{model.chrom}:{model.cds_start0 + 1}..{model.cds_end0}",
                "strand": model.strand,
                "cds_length": model.cds_length,
                "cds_exon_count": model.cds_exon_count,
                "cds_intron_count": len(model.cds_introns),
                "sublocus_model_count": len(models),
                "sublocus_source_count": len(sub_sources),
                "max_possible_other_source_count": possible_other_source_count,
                "same_intron_chain_model_count": len(same_other_model_ids),
                "same_intron_chain_source_count": len(same_other_sources),
                "same_intron_chain_sources": fmt_sources(same_other_sources),
                "same_intron_chain_model_ids": fmt_model_ids(same_other_model_ids),
                "exact_intron_chain_support_fraction": exact_chain_support_fraction,
                "query_chain_contained_in_other_model_count": len({model_key(m) for m in query_contained_models}),
                "query_chain_contained_in_other_source_count": len(query_contained_sources),
                "query_chain_contained_in_other_sources": fmt_sources(query_contained_sources),
                "other_chain_contained_in_query_model_count": len({model_key(m) for m in other_contained_models}),
                "other_chain_contained_in_query_source_count": len(other_contained_sources),
                "other_chain_contained_in_query_sources": fmt_sources(other_contained_sources),
                "cds_introns_total": intron_summary["total"],
                "cds_introns_supported_by_other_model": intron_summary["supported_by_other_model"],
                "cds_introns_supported_by_other_source": intron_summary["supported_by_other_source"],
                "cds_introns_supported_by_2plus_sources": intron_summary["supported_by_2plus_sources"],
                "cds_intron_support_source_sum": intron_summary["support_source_sum"],
                "cds_intron_support_source_mean": intron_summary["support_source_mean"],
                "cds_intron_support_source_min": intron_summary["support_source_min"],
                "unsupported_cds_intron_count": intron_summary["unsupported_count"],
                "unsupported_cds_intron_fraction": intron_summary["unsupported_fraction"],
                "unsupported_cds_intron_contained_in_other_cds_exon_count": unsupported_contained_intron_summary["count"],
                "unsupported_cds_intron_contained_in_other_cds_exon_fraction": unsupported_contained_intron_summary["fraction"],
                "unsupported_cds_intron_contained_in_other_cds_exon_other_model_count": unsupported_contained_intron_summary["other_model_count"],
                "unsupported_cds_intron_contained_in_other_cds_exon_other_source_count": unsupported_contained_intron_summary["other_source_count"],
                "unsupported_cds_intron_contained_in_other_cds_exon_other_sources": unsupported_contained_intron_summary["other_sources"],
                "unsupported_cds_intron_contained_in_other_cds_exon_other_model_ids": unsupported_contained_intron_summary["other_model_ids"],
                "unsupported_cds_intron_contained_in_other_cds_exon_introns": unsupported_contained_intron_summary["introns"],
                "cds_intron_support_fraction": intron_summary["support_fraction"],
                "cds_splice_sites_total": splice_summary["total"],
                "cds_splice_sites_supported_by_other_model": splice_summary["supported_by_other_model"],
                "cds_splice_sites_supported_by_other_source": splice_summary["supported_by_other_source"],
                "cds_splice_sites_supported_by_2plus_sources": splice_summary["supported_by_2plus_sources"],
                "cds_splice_site_support_source_sum": splice_summary["support_source_sum"],
                "cds_splice_site_support_source_mean": splice_summary["support_source_mean"],
                "cds_splice_site_support_source_min": splice_summary["support_source_min"],
                "unsupported_cds_splice_site_count": splice_summary["unsupported_count"],
                "unsupported_cds_splice_site_fraction": splice_summary["unsupported_fraction"],
                "cds_splice_site_support_fraction": splice_summary["support_fraction"],
                "internal_cds_exons_total": internal_exon_summary["total"],
                "internal_cds_exons_supported_by_other_model": internal_exon_summary["supported_by_other_model"],
                "internal_cds_exons_supported_by_other_source": internal_exon_summary["supported_by_other_source"],
                "internal_cds_exons_supported_by_2plus_sources": internal_exon_summary["supported_by_2plus_sources"],
                "internal_cds_exon_support_source_sum": internal_exon_summary["support_source_sum"],
                "internal_cds_exon_support_source_mean": internal_exon_summary["support_source_mean"],
                "internal_cds_exon_support_source_min": internal_exon_summary["support_source_min"],
                "unsupported_internal_cds_exon_count": internal_exon_summary["unsupported_count"],
                "unsupported_internal_cds_exon_fraction": internal_exon_summary["unsupported_fraction"],
                "internal_cds_exon_support_fraction": internal_exon_summary["support_fraction"],
                "terminal_cds_exons_total": terminal_summary["total"],
                "terminal_cds_exons_supported_by_other_model": terminal_summary["supported_by_other_model"],
                "terminal_cds_exons_supported_by_other_source": terminal_summary["supported_by_other_source"],
                "terminal_cds_exons_supported_by_2plus_sources": terminal_summary["supported_by_2plus_sources"],
                "terminal_cds_exon_support_source_sum": terminal_summary["support_source_sum"],
                "terminal_cds_exon_support_source_mean": terminal_summary["support_source_mean"],
                "terminal_cds_exon_support_source_min": terminal_summary["support_source_min"],
                "unsupported_terminal_cds_exon_count": terminal_summary["unsupported_count"],
                "unsupported_terminal_cds_exon_fraction": terminal_summary["unsupported_fraction"],
                "terminal_cds_exon_support_fraction": terminal_summary["support_fraction"],
                "terminal_cds_3prime_wobble_bp": terminal_cds_3prime_wobble_bp,
                "single_exon_cds_exact_supported_feature_count": single_exact_supported_feature_count,
                "single_exon_cds_exact_support_model_count": single_exact_support_model_count,
                "single_exon_cds_exact_support_source_count": single_exact_support_source_count,
                "single_exon_cds_exact_support_sources": single_exact_support_sources,
                "single_exon_exact_support_fraction": single_exact_support_fraction,
                "single_exon_cds_overlap_support_model_count": single_overlap_support_model_count,
                "single_exon_cds_overlap_support_source_count": single_overlap_support_source_count,
                "single_exon_cds_overlap_support_sources": single_overlap_support_sources,
                "single_exon_overlap_support_fraction": single_overlap_support_fraction,
                "single_exon_best_reciprocal_overlap": best_single_overlap,
                "single_exon_best_reciprocal_overlap_model": best_single_overlap_model,
                "single_exon_best_reciprocal_overlap_source": best_single_overlap_source,
                "exact_cds_exon_chain_support_model_count": len(exact_cds_exon_chain_support_model_ids),
                "exact_cds_exon_chain_support_source_count": exact_cds_exon_chain_support_source_count,
                "exact_cds_exon_chain_support_sources": fmt_sources(exact_cds_exon_chain_support_sources),
                "exact_cds_exon_chain_support_model_ids": fmt_model_ids(exact_cds_exon_chain_support_model_ids),
                "exact_cds_exon_chain_support_fraction_local": exact_cds_exon_chain_support_fraction_local,
                "exact_structure_support_score": exact_structure_support_score,
                "structural_support_score": structural_score,
                "structural_support_score_type": score_type,
                "structural_support_score_components": component_text,
                "consensus_similarity_score": consensus_similarity_score,
                "consensus_similarity_score_type": consensus_similarity_score_type,
                "consensus_similarity_score_components": consensus_components,
                "consensus_similarity_precision": consensus_precision,
                "consensus_similarity_sensitivity": consensus_sensitivity,
                "partial_structure_penalty": partial_structure_penalty,
                "extra_structure_penalty": extra_structure_penalty,
                "consensus_cds_bp": consensus_cds_bp,
                "query_cds_bp": query_cds_bp,
                "consensus_cds_bp_recovered": consensus_cds_bp_recovered,
                "missing_supported_cds_bp": missing_supported_cds_bp,
                "missing_supported_cds_fraction": missing_supported_cds_fraction,
                "missing_supported_terminal_structure_bp": missing_supported_terminal_structure_bp,
                "missing_supported_terminal_structure_fraction": missing_supported_terminal_structure_fraction,
                "missing_supported_terminal_5prime_bp": missing_supported_terminal_5prime_bp,
                "missing_supported_terminal_3prime_bp": missing_supported_terminal_3prime_bp,
                "missing_supported_internal_structure_bp": missing_supported_internal_structure_bp,
                "missing_supported_internal_structure_fraction": missing_supported_internal_structure_fraction,
                "extra_cds_bp": extra_cds_bp,
                "extra_cds_fraction": extra_cds_fraction,
                "terminal_extra_bp": terminal_extra_bp,
                "terminal_extra_fraction": terminal_extra_fraction,
                "terminal_extra_5prime_bp": terminal_extra_5prime_bp,
                "terminal_extra_3prime_bp": terminal_extra_3prime_bp,
                "internal_extra_bp": internal_extra_bp,
                "internal_extra_fraction": internal_extra_fraction,
                "terminal_extension_score": terminal_extension_score,
                "contained_model_source_fraction": contained_model_source_fraction,
                "model_completeness_score": model_completeness_score,
                "model_completeness_penalty": model_completeness_penalty,
                "model_completeness_score_components": model_completeness_score_components,
                "pairwise_similarity_score": pairwise_similarity_score,
                "pairwise_similarity_score_type": "single_exon" if model.is_single_exon else "multi_exon",
                "pairwise_similarity_score_components": pairwise_components,
                "pairwise_other_source_count": pairwise_source_count,
                "best_pairwise_model_id": best_overall_model,
                "best_pairwise_source": best_overall_source,
                "best_pairwise_similarity_score": best_pairwise_metrics.get("score"),
                "intrinsic_penalty_score": intrinsic_penalty_score,
                "intrinsic_penalty_components": intrinsic_penalty_components,
            }
            row.update(unsupported_short_exon_summary)
            row.update(compact_sequence_integrity_fields(model))
            row.update(prf_values("cds_base_consensus", consensus_base))
            row.update(prf_values("cds_intron_consensus", consensus_intron))
            row.update(prf_values("cds_splice_site_consensus", consensus_splice))
            row.update(prf_values("internal_cds_exon_consensus", consensus_internal))
            row.update(prf_values("terminal_cds_exon_consensus", consensus_terminal))
            for field in [
                "cds_base", "cds_intron", "cds_splice_site", "internal_cds_exon", "terminal_cds_exon"
            ]:
                for metric in ["precision", "sensitivity", "f1"]:
                    row[f"mean_best_per_source_{field}_{metric}"] = mean_best(f"{field}_{metric}")
                    row[f"best_pairwise_{field}_{metric}"] = best_pairwise_metrics.get(f"{field}_{metric}")
            for name, count in intrinsic_counts.items():
                row[f"intrinsic_penalty_{name}_count"] = count
            if calculate_corroboration:
                corroboration_row = make_structural_corroboration_row(
                    model,
                    sub,
                    corroboration_comparisons,
                    run_source_count=run_source_count,
                    exact_cds_exon_chain_support_source_count=exact_cds_exon_chain_support_source_count,
                )
                for field in (
                    "run_model_source_count",
                    "run_possible_other_source_count",
                    "corroborating_other_source_count",
                    "corroborating_other_sources",
                    "corroborating_model_count",
                    "corroborating_model_ids",
                ):
                    row[field] = corroboration_row[field]
                if structural_corroboration_rows is not None:
                    structural_corroboration_rows.append(corroboration_row)
            model_rows.append(row)
            processed_models += 1
            if progress is not None:
                progress.every(
                    processed_models,
                    total_models,
                    "Model scoring structural support models",
                    level="verbose",
                )

    if progress is not None:
        progress.log(
            f"Model scoring structural support: finished {len(model_rows):,} model rows / {len(feature_rows):,} feature rows",
            level="basic",
        )

    return feature_rows, model_rows


def write_structural_features_tsv(feature_rows: list[dict[str, object]], path: Path) -> None:
    columns = [
        "superlocus_id", "sublocus_id", "sublocus_type", "feature_type", "matching_mode", "terminal_side",
        "chrom", "start", "end", "strand", "start0", "end0",
        "support_model_count", "support_source_count", "support_sources", "support_model_ids",
    ]
    with path.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for row in sorted(feature_rows, key=lambda r: (str(r["sublocus_id"]), str(r["feature_type"]), str(r["chrom"]), str(r["strand"]), str(r["start0"]), str(r["end0"]), str(r["terminal_side"]))):
            out.write("\t".join(str(row.get(col, ".")) for col in columns) + "\n")


def write_model_structural_support_tsv(model_rows: list[dict[str, object]], path: Path) -> None:
    base_columns = [
        "superlocus_id", "sublocus_id", "sublocus_type", "source_name", "input_file", "input_gene_id", "input_transcript_id",
        "output_gene_id", "output_transcript_id", "chrom", "cds_start", "cds_end", "cds_location", "strand", "cds_length", "cds_exon_count", "cds_intron_count",
        "sublocus_model_count", "sublocus_source_count", "max_possible_other_source_count",
        "same_intron_chain_model_count", "same_intron_chain_source_count", "same_intron_chain_sources", "same_intron_chain_model_ids", "exact_intron_chain_support_fraction",
        "query_chain_contained_in_other_model_count", "query_chain_contained_in_other_source_count", "query_chain_contained_in_other_sources",
        "other_chain_contained_in_query_model_count", "other_chain_contained_in_query_source_count", "other_chain_contained_in_query_sources",
        "cds_introns_total", "cds_introns_supported_by_other_model", "cds_introns_supported_by_other_source", "cds_introns_supported_by_2plus_sources",
        "cds_intron_support_source_sum", "cds_intron_support_source_mean", "cds_intron_support_source_min", "unsupported_cds_intron_count", "unsupported_cds_intron_fraction",
        "unsupported_cds_intron_contained_in_other_cds_exon_count", "unsupported_cds_intron_contained_in_other_cds_exon_fraction",
        "unsupported_cds_intron_contained_in_other_cds_exon_other_model_count", "unsupported_cds_intron_contained_in_other_cds_exon_other_source_count",
        "unsupported_cds_intron_contained_in_other_cds_exon_other_sources", "unsupported_cds_intron_contained_in_other_cds_exon_other_model_ids",
        "unsupported_cds_intron_contained_in_other_cds_exon_introns", "cds_intron_support_fraction",
        "cds_splice_sites_total", "cds_splice_sites_supported_by_other_model", "cds_splice_sites_supported_by_other_source", "cds_splice_sites_supported_by_2plus_sources",
        "cds_splice_site_support_source_sum", "cds_splice_site_support_source_mean", "cds_splice_site_support_source_min", "unsupported_cds_splice_site_count", "unsupported_cds_splice_site_fraction", "cds_splice_site_support_fraction",
        "internal_cds_exons_total", "internal_cds_exons_supported_by_other_model", "internal_cds_exons_supported_by_other_source", "internal_cds_exons_supported_by_2plus_sources",
        "internal_cds_exon_support_source_sum", "internal_cds_exon_support_source_mean", "internal_cds_exon_support_source_min", "unsupported_internal_cds_exon_count", "unsupported_internal_cds_exon_fraction", "internal_cds_exon_support_fraction",
        "terminal_cds_exons_total", "terminal_cds_exons_supported_by_other_model", "terminal_cds_exons_supported_by_other_source", "terminal_cds_exons_supported_by_2plus_sources",
        "terminal_cds_exon_support_source_sum", "terminal_cds_exon_support_source_mean", "terminal_cds_exon_support_source_min", "unsupported_terminal_cds_exon_count", "unsupported_terminal_cds_exon_fraction", "terminal_cds_exon_support_fraction", "terminal_cds_3prime_wobble_bp",
        "unsupported_short_5prime_terminal_cds_exon_count", "unsupported_short_5prime_terminal_cds_exon_length_bp", "unsupported_short_5prime_terminal_cds_exon_supported_by_other_source", "unsupported_short_5prime_terminal_cds_exon_severity",
        "unsupported_short_3prime_terminal_cds_exon_count", "unsupported_short_3prime_terminal_cds_exon_length_bp", "unsupported_short_3prime_terminal_cds_exon_supported_by_other_source", "unsupported_short_3prime_terminal_cds_exon_severity",
        "moderately_short_internal_cds_exon_count", "supported_moderately_short_internal_cds_exon_count", "unsupported_moderately_short_internal_cds_exon_count", "unsupported_moderately_short_internal_cds_exons", "unsupported_short_internal_cds_exon_severity",
        "unsupported_short_terminal_intrinsic_threshold_bp", "unsupported_short_internal_intrinsic_threshold_bp",
        "unsupported_short_5prime_terminal_soft_threshold_bp", "unsupported_short_3prime_terminal_soft_threshold_bp", "unsupported_short_internal_soft_threshold_bp",
        "unsupported_short_5prime_terminal_component_weight", "unsupported_short_3prime_terminal_component_weight", "unsupported_short_internal_component_weight",
        "unsupported_short_cds_exon_penalty_score", "unsupported_short_cds_exon_penalty_components",
        "single_exon_cds_exact_supported_feature_count", "single_exon_cds_exact_support_model_count", "single_exon_cds_exact_support_source_count", "single_exon_cds_exact_support_sources", "single_exon_exact_support_fraction",
        "single_exon_cds_overlap_support_model_count", "single_exon_cds_overlap_support_source_count", "single_exon_cds_overlap_support_sources", "single_exon_overlap_support_fraction",
        "single_exon_best_reciprocal_overlap", "single_exon_best_reciprocal_overlap_model", "single_exon_best_reciprocal_overlap_source",
        "exact_cds_exon_chain_support_model_count", "exact_cds_exon_chain_support_source_count", "exact_cds_exon_chain_support_sources", "exact_cds_exon_chain_support_model_ids",
        "exact_cds_exon_chain_support_fraction_local", "exact_structure_support_score",
        "structural_support_score", "structural_support_score_type", "structural_support_score_components",
    ]
    consensus_columns = []
    for prefix in [
        "cds_base_consensus", "cds_intron_consensus", "cds_splice_site_consensus", "internal_cds_exon_consensus", "terminal_cds_exon_consensus"
    ]:
        consensus_columns.extend([f"{prefix}_precision", f"{prefix}_sensitivity", f"{prefix}_f1"])
    consensus_columns.extend([
        "consensus_similarity_score", "consensus_similarity_score_type", "consensus_similarity_score_components",
        "consensus_similarity_precision", "consensus_similarity_sensitivity",
        "partial_structure_penalty", "extra_structure_penalty",
    ])
    completeness_columns = [
        "consensus_cds_bp", "query_cds_bp", "consensus_cds_bp_recovered",
        "missing_supported_cds_bp", "missing_supported_cds_fraction",
        "missing_supported_terminal_structure_bp", "missing_supported_terminal_structure_fraction",
        "missing_supported_terminal_5prime_bp", "missing_supported_terminal_3prime_bp",
        "missing_supported_internal_structure_bp", "missing_supported_internal_structure_fraction",
        "extra_cds_bp", "extra_cds_fraction",
        "terminal_extra_bp", "terminal_extra_fraction", "terminal_extra_5prime_bp", "terminal_extra_3prime_bp",
        "internal_extra_bp", "internal_extra_fraction",
        "terminal_extension_score", "contained_model_source_fraction",
        "model_completeness_score", "model_completeness_penalty", "model_completeness_score_components",
    ]
    pairwise_columns = [
        "pairwise_similarity_score", "pairwise_similarity_score_type", "pairwise_similarity_score_components",
        "pairwise_other_source_count", "best_pairwise_model_id", "best_pairwise_source", "best_pairwise_similarity_score",
    ]
    for prefix in ["cds_base", "cds_intron", "cds_splice_site", "internal_cds_exon", "terminal_cds_exon"]:
        for metric in ["precision", "sensitivity", "f1"]:
            pairwise_columns.append(f"mean_best_per_source_{prefix}_{metric}")
        for metric in ["precision", "sensitivity", "f1"]:
            pairwise_columns.append(f"best_pairwise_{prefix}_{metric}")
    intrinsic_columns = [
        "intrinsic_penalty_score", "intrinsic_penalty_components",
        "intrinsic_penalty_short_cds_count", "intrinsic_penalty_cds_length_not_multiple_of_3_count",
        "intrinsic_penalty_short_terminal_cds_exon_count", "intrinsic_penalty_short_internal_cds_exon_count",
        "intrinsic_penalty_short_cds_intron_count", "intrinsic_penalty_long_cds_intron_count",
        "intrinsic_penalty_contained_model_count",
        "intrinsic_penalty_unsupported_cds_intron_contained_in_other_cds_exon_count",
    ]
    sequence_columns = [
        "sequence_integrity_status", "sequence_integrity_penalty_score", "sequence_integrity_flags",
    ]
    model_selection_columns = [
        "source_prior_score",
        "low_source_conflict_flag", "low_source_conflict_agreement_factor", "low_source_conflict_reason",
        "structural_coherence_penalty_score", "structural_coherence_raw_score",
        "structural_coherence_penalty_weight", "structural_coherence_penalty_in_selection",
        "structural_coherence_region_count", "structural_coherence_max_region",
        "model_structural_support_score_before_sequence_integrity", "model_structural_support_score_before_sequence_integrity_scaled",
        "sequence_integrity_penalty_weight", "sequence_integrity_penalty_contribution", "sequence_integrity_penalty_in_selection",
        "model_structural_support_score", "model_structural_support_score_scaled",
        "model_structural_support_score_rank", "model_structural_support_order_rank",
        "model_structural_support_tie_group", "model_structural_support_tie_break_reason",
        "model_structural_support_score_components",
        "protein_evidence_score", "protein_evidence_score_weight", "protein_evidence_score_in_selection",
        "model_selection_score", "model_selection_score_scaled", "model_selection_score_rank", "model_selection_order_rank",
        "exact_cds_exon_chain_support_model_count", "exact_cds_exon_chain_support_source_count",
        "exact_cds_exon_chain_support_sources", "exact_cds_exon_chain_support_model_ids",
        "exact_cds_exon_chain_support_fraction_local", "exact_structure_support_score",
        "model_selection_tie_group", "model_selection_tie_break_reason",
        "model_selection_score_components",
        "rank_changed_by_evidence", "order_rank_delta_from_evidence", "score_delta_from_evidence", "evidence_changed_top_rank",
    ]
    corroboration_columns = []
    if any(value_as_float(row.get("structural_corroboration_weight")) > 0 for row in model_rows):
        corroboration_columns = [
            "corroborating_other_source_count", "run_possible_other_source_count",
            "structural_corroboration_score", "structural_corroboration_alpha",
            "structural_corroboration_weight", "structural_corroboration_contribution",
            "structural_corroboration_in_selection",
        ]
    columns = (
        base_columns + consensus_columns + completeness_columns + pairwise_columns
        + intrinsic_columns + sequence_columns + corroboration_columns + model_selection_columns
    )

    def format_value(value: object) -> str:
        if value is None:
            return "."
        if isinstance(value, float):
            return fmt_optional_float(value, 6)
        return str(value)

    with path.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for row in sorted(model_rows, key=lambda r: (str(r["chrom"]), int(r["cds_start"]), int(r["cds_end"]), str(r["source_name"]), str(r["input_transcript_id"]))):
            out.write("\t".join(format_value(row.get(col, ".")) for col in columns) + "\n")


MODEL_SELECTION_SCORE_COLUMNS = [
    "superlocus_id", "sublocus_id", "sublocus_type", "source_name", "input_file", "input_gene_id", "input_transcript_id",
    "output_gene_id", "output_transcript_id", "chrom", "cds_start", "cds_end", "cds_location", "strand", "cds_length", "cds_exon_count", "cds_intron_count",
    "sublocus_model_count", "sublocus_source_count", "source_prior_score",
    "structural_support_score", "consensus_similarity_score", "pairwise_similarity_score", "model_completeness_score",
    "terminal_extension_score", "exact_structure_support_score", "intrinsic_penalty_score", "partial_structure_penalty",
    "internal_extra_fraction", "terminal_extra_fraction",
    "unsupported_short_5prime_terminal_cds_exon_count", "unsupported_short_5prime_terminal_cds_exon_length_bp", "unsupported_short_5prime_terminal_cds_exon_supported_by_other_source", "unsupported_short_5prime_terminal_cds_exon_severity",
    "unsupported_short_3prime_terminal_cds_exon_count", "unsupported_short_3prime_terminal_cds_exon_length_bp", "unsupported_short_3prime_terminal_cds_exon_supported_by_other_source", "unsupported_short_3prime_terminal_cds_exon_severity",
    "moderately_short_internal_cds_exon_count", "supported_moderately_short_internal_cds_exon_count", "unsupported_moderately_short_internal_cds_exon_count", "unsupported_moderately_short_internal_cds_exons", "unsupported_short_internal_cds_exon_severity",
    "unsupported_short_cds_exon_penalty_score", "unsupported_short_cds_exon_penalty_components", "unsupported_short_cds_exon_penalty_weight", "unsupported_short_cds_exon_penalty_in_selection",
    "structural_coherence_penalty_score", "structural_coherence_raw_score",
    "structural_coherence_penalty_weight", "structural_coherence_penalty_in_selection",
    "structural_coherence_region_count", "structural_coherence_max_region",
    "sequence_integrity_status", "sequence_integrity_penalty_score", "sequence_integrity_flags",
    "model_structural_support_score_before_sequence_integrity", "model_structural_support_score_before_sequence_integrity_scaled",
    "sequence_integrity_penalty_weight", "sequence_integrity_penalty_contribution", "sequence_integrity_penalty_in_selection",
    "model_structural_support_score", "model_structural_support_score_scaled",
    "model_structural_support_score_rank", "model_structural_support_order_rank",
    "model_structural_support_tie_group", "model_structural_support_tie_break_reason",
    "model_structural_support_score_components",
    "protein_evidence_score", "protein_evidence_score_scaled", "protein_evidence_score_weight", "protein_evidence_score_in_selection",
    "protein_evidence_status", "protein_evidence_flags", "protein_evidence_score_components",
    "transcript_evidence_score", "transcript_evidence_score_scaled", "transcript_evidence_score_context_adjusted", "transcript_evidence_score_context_adjusted_scaled",
    "transcript_evidence_score_weight", "transcript_evidence_score_in_selection",
    "transcript_evidence_status", "transcript_evidence_flags", "transcript_evidence_score_components",
    "transcript_cdna_context_status", "transcript_cdna_context_cds_to_min_cdna_fraction", "transcript_cdna_context_flags",
    "transcript_cdna_context_evaluable", "transcript_cdna_context_concern_score", "transcript_cdna_context_concern_class",
    "transcript_cdna_context_pruning_concern", "transcript_cdna_context_pruning_concern_class", "transcript_cdna_context_pruning_supporting_flags",
    "transcript_cdna_context_neutral_fraction", "transcript_cdna_context_severe_fraction",
    "transcript_cdna_context_penalty_weight", "transcript_cdna_context_penalty", "transcript_cdna_context_penalty_in_selection",
    "model_selection_score", "model_selection_score_scaled", "model_selection_score_rank", "model_selection_order_rank",
    "model_selection_tie_group", "model_selection_tie_break_reason", "model_selection_score_components",
    "rank_changed_by_evidence", "order_rank_delta_from_evidence", "score_delta_from_evidence", "evidence_changed_top_rank",
    "low_source_conflict_flag", "low_source_conflict_agreement_factor", "low_source_conflict_reason",
]


def write_model_selection_scores_tsv(model_rows: list[dict[str, object]], path: Path) -> None:
    """Write final model-selection score audit table.

    This table separates the model-only structural-support aggregate from the
    final model selection score, which may include optional evidence terms.
    """
    def format_value(value: object) -> str:
        if value is None:
            return "."
        if isinstance(value, float):
            return fmt_optional_float(value, 6)
        return str(value)

    columns = list(MODEL_SELECTION_SCORE_COLUMNS)
    if any(value_as_float(row.get("structural_corroboration_weight")) > 0 for row in model_rows):
        insertion = columns.index("model_structural_support_score")
        columns[insertion:insertion] = [
            "corroborating_other_source_count", "run_possible_other_source_count",
            "structural_corroboration_score", "structural_corroboration_alpha",
            "structural_corroboration_weight", "structural_corroboration_contribution",
            "structural_corroboration_in_selection",
        ]

    with path.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for row in sorted(model_rows, key=lambda r: (str(r["chrom"]), int(r["cds_start"]), int(r["cds_end"]), str(r["source_name"]), str(r["input_transcript_id"]))):
            out.write("\t".join(format_value(row.get(col, ".")) for col in columns) + "\n")


MODEL_SELECTION_SCORE_SUMMARY_COLUMNS = [
    "section", "group", "metric", "value", "details",
    "mean_scaled_score_all", "mean_scaled_score_non_singletons",
    "n_scored_all", "n_scored_non_singletons", "n_scored_singletons",
    "singleton_fraction",
]


def write_model_selection_score_summary_tsv(model_rows: list[dict[str, object]], path: Path) -> None:
    """Write a compact summary of model-only versus final selection scoring.

    The summary remains row-oriented and keeps its legacy five columns first;
    source diagnostics use explicitly appended columns. It complements
    model_selection_scores.tsv, which remains the per-model audit table.
    """
    rows: list[dict[str, object]] = []

    def add(section: str, group: str, metric: str, value: object, details: object = ".") -> None:
        rows.append({"section": section, "group": group, "metric": metric, "value": value, "details": details})

    total_models = len(model_rows)
    subloci = {str(row.get("sublocus_id", ".")) for row in model_rows}
    add("overall", ".", "models", total_models)
    add("overall", ".", "subloci", len(subloci))

    source_score_stats = summarise_source_scaled_scores(model_rows)
    source_order = sorted(
        source_score_stats,
        key=lambda source: (
            -(
                source_score_stats[source]["scaled_sum"]
                / source_score_stats[source]["count"]
            ),
            source,
        ),
    )
    for source in source_order:
        stats = source_score_stats[source]
        n_all = int(stats["count"])
        n_non_singletons = int(stats["non_singleton_count"])
        n_singletons = int(stats["singleton_count"])
        mean_all = stats["scaled_sum"] / stats["count"]
        mean_non_singletons = (
            stats["non_singleton_scaled_sum"] / stats["non_singleton_count"]
            if stats["non_singleton_count"]
            else None
        )
        rows.append(
            {
                "section": "source_scaled_score",
                "group": source,
                "metric": "source_scaled_score_summary",
                "value": mean_all,
                "details": "singleton_sublocus=sublocus_model_count==1",
                "mean_scaled_score_all": mean_all,
                "mean_scaled_score_non_singletons": mean_non_singletons,
                "n_scored_all": n_all,
                "n_scored_non_singletons": n_non_singletons,
                "n_scored_singletons": n_singletons,
                "singleton_fraction": n_singletons / n_all,
            }
        )

    protein_weights = [value_as_float(row.get("protein_evidence_score_weight")) for row in model_rows]
    effective_protein_weight = max(protein_weights) if protein_weights else 0.0
    protein_active = any(str(row.get("protein_evidence_score_in_selection", ".")).lower() == "true" for row in model_rows)
    transcript_weights = [value_as_float(row.get("transcript_evidence_score_weight")) for row in model_rows]
    effective_transcript_weight = max(transcript_weights) if transcript_weights else 0.0
    transcript_active = any(str(row.get("transcript_evidence_score_in_selection", ".")).lower() == "true" for row in model_rows)
    add("selection_weights", ".", "effective_protein_evidence_weight", effective_protein_weight)
    add("selection_weights", ".", "protein_evidence_in_selection", "true" if protein_active else "false")
    add("selection_weights", ".", "effective_transcript_evidence_weight", effective_transcript_weight)
    add("selection_weights", ".", "transcript_evidence_in_selection", "true" if transcript_active else "false")
    sequence_weights = [value_as_float(row.get("sequence_integrity_penalty_weight")) for row in model_rows]
    effective_sequence_weight = max(sequence_weights) if sequence_weights else 0.0
    sequence_active = any(str(row.get("sequence_integrity_penalty_in_selection", ".")).lower() == "true" for row in model_rows)
    add("selection_weights", ".", "effective_sequence_integrity_penalty_weight", effective_sequence_weight)
    add("selection_weights", ".", "sequence_integrity_penalty_in_selection", "true" if sequence_active else "false")
    coherence_weights = [value_as_float(row.get("structural_coherence_penalty_weight")) for row in model_rows]
    effective_coherence_weight = max(coherence_weights) if coherence_weights else 0.0
    coherence_active = any(str(row.get("structural_coherence_penalty_in_selection", ".")).lower() == "true" for row in model_rows)
    add("selection_weights", ".", "effective_structural_coherence_penalty_weight", effective_coherence_weight)
    add("selection_weights", ".", "structural_coherence_penalty_in_selection", "true" if coherence_active else "false")

    corroboration_weights = [value_as_float(row.get("structural_corroboration_weight")) for row in model_rows]
    effective_corroboration_weight = max(corroboration_weights) if corroboration_weights else 0.0
    corroboration_active = effective_corroboration_weight > 0
    if corroboration_active:
        corroboration_alphas = [value_as_float(row.get("structural_corroboration_alpha")) for row in model_rows]
        add("selection_weights", ".", "effective_structural_corroboration_weight", effective_corroboration_weight)
        add("selection_weights", ".", "effective_structural_corroboration_alpha", max(corroboration_alphas) if corroboration_alphas else 0.0)
        add("selection_weights", ".", "structural_corroboration_in_selection", "true")

    score_distribution_fields = [
        "model_structural_support_score",
        "model_structural_support_score_scaled",
        "protein_evidence_score",
        "protein_evidence_score_scaled",
        "transcript_evidence_score",
        "transcript_evidence_score_scaled",
        "transcript_evidence_score_context_adjusted",
        "transcript_evidence_score_context_adjusted_scaled",
        "transcript_cdna_context_concern_score",
        "transcript_cdna_context_penalty",
        "unsupported_short_cds_exon_penalty_score",
        "structural_coherence_raw_score",
        "structural_coherence_penalty_score",
        "sequence_integrity_penalty_score",
        "sequence_integrity_penalty_contribution",
        "model_structural_support_score_before_sequence_integrity",
        "model_structural_support_score_before_sequence_integrity_scaled",
        "model_selection_score",
        "model_selection_score_scaled",
        "score_delta_from_evidence",
    ]
    if corroboration_active:
        score_distribution_fields[0:0] = [
            "corroborating_other_source_count",
            "structural_corroboration_score",
            "structural_corroboration_contribution",
        ]
    for field in score_distribution_fields:
        values = [value_as_float(row.get(field)) for row in model_rows if row.get(field, ".") != "."]
        if not values:
            continue
        values_sorted = sorted(values)
        add("score_distribution", field, "min", values_sorted[0])
        add("score_distribution", field, "median", values_sorted[len(values_sorted)//2])
        add("score_distribution", field, "mean", sum(values_sorted) / len(values_sorted))
        add("score_distribution", field, "max", values_sorted[-1])

    add(
        "unsupported_short_cds_exon_penalty",
        ".",
        "models_with_penalty_gt_0",
        sum(1 for row in model_rows if value_as_float(row.get("unsupported_short_cds_exon_penalty_score")) > 0),
    )
    add(
        "structural_coherence_penalty",
        ".",
        "models_with_penalty_gt_0",
        sum(1 for row in model_rows if value_as_float(row.get("structural_coherence_penalty_score")) > 0),
    )
    add(
        "structural_coherence_penalty",
        ".",
        "eligible_region_count",
        sum(int(value_as_float(row.get("structural_coherence_region_count"))) for row in model_rows),
    )
    add(
        "unsupported_short_cds_exon_penalty",
        ".",
        "models_with_5prime_terminal_component",
        sum(1 for row in model_rows if value_as_float(row.get("unsupported_short_5prime_terminal_cds_exon_count")) > 0),
    )
    add(
        "unsupported_short_cds_exon_penalty",
        ".",
        "models_with_3prime_terminal_component",
        sum(1 for row in model_rows if value_as_float(row.get("unsupported_short_3prime_terminal_cds_exon_count")) > 0),
    )
    add(
        "unsupported_short_cds_exon_penalty",
        ".",
        "models_with_internal_component",
        sum(1 for row in model_rows if value_as_float(row.get("unsupported_moderately_short_internal_cds_exon_count")) > 0),
    )

    add(
        "sequence_integrity",
        ".",
        "models_evaluated",
        sum(1 for row in model_rows if str(row.get("sequence_integrity_status", "")).startswith("evaluated_")),
    )
    add(
        "sequence_integrity",
        ".",
        "models_with_penalty_gt_0",
        sum(1 for row in model_rows if value_as_float(row.get("sequence_integrity_penalty_score")) > 0),
    )
    add(
        "sequence_integrity",
        ".",
        "models_with_internal_stop",
        sum(
            1
            for row in model_rows
            if "internal_stop_codon" in str(row.get("sequence_integrity_flags", ".")).split(",")
        ),
    )
    add(
        "sequence_integrity",
        ".",
        "models_with_noncanonical_splice",
        sum(
            1
            for row in model_rows
            if "noncanonical_splice" in str(row.get("sequence_integrity_flags", ".")).split(",")
        ),
    )

    add("rank_changes", ".", "rows_rank_changed_by_evidence", sum(1 for row in model_rows if str(row.get("rank_changed_by_evidence", ".")).lower() == "true"))
    add("rank_changes", ".", "rows_evidence_changed_top_rank", sum(1 for row in model_rows if str(row.get("evidence_changed_top_rank", ".")).lower() == "true"))
    changed_top_subloci = {
        str(row.get("sublocus_id", "."))
        for row in model_rows
        if str(row.get("evidence_changed_top_rank", ".")).lower() == "true"
    }
    add("rank_changes", ".", "subloci_with_evidence_changed_top_rank", len(changed_top_subloci))

    protein_scores_positive = sum(1 for row in model_rows if value_as_float(row.get("protein_evidence_score")) > 0)
    protein_scores_strong = sum(1 for row in model_rows if value_as_float(row.get("protein_evidence_score_scaled")) >= 0.75)
    add("protein_evidence", ".", "models_with_protein_evidence_score_gt_0", protein_scores_positive)
    add("protein_evidence", ".", "models_with_protein_evidence_score_scaled_ge_0.75", protein_scores_strong)

    status_counts = Counter(str(row.get("protein_evidence_status", ".")) for row in model_rows)
    for status, count in sorted(status_counts.items()):
        add("protein_evidence_status", status, "models", count)

    flag_counts: Counter[str] = Counter()
    for row in model_rows:
        flags = str(row.get("protein_evidence_flags", "."))
        if flags and flags != ".":
            for flag in flags.split(","):
                if flag:
                    flag_counts[flag] += 1
    for flag, count in sorted(flag_counts.items()):
        add("protein_evidence_flags", flag, "models", count)

    transcript_scores_positive = sum(1 for row in model_rows if value_as_float(row.get("transcript_evidence_score")) > 0)
    transcript_scores_strong = sum(1 for row in model_rows if value_as_float(row.get("transcript_evidence_score_scaled")) >= 0.75)
    add("transcript_evidence", ".", "models_with_transcript_evidence_score_gt_0", transcript_scores_positive)
    add("transcript_evidence", ".", "models_with_transcript_evidence_score_scaled_ge_0.75", transcript_scores_strong)

    add(
        "transcript_evidence",
        ".",
        "models_with_transcript_cdna_context_penalty_in_selection",
        sum(1 for row in model_rows if str(row.get("transcript_cdna_context_penalty_in_selection", ".")).lower() == "true"),
    )
    add(
        "transcript_evidence",
        ".",
        "models_with_transcript_cdna_context_pruning_concern",
        sum(1 for row in model_rows if str(row.get("transcript_cdna_context_pruning_concern", ".")).lower() == "true"),
    )
    add(
        "transcript_evidence",
        ".",
        "models_with_transcript_cdna_context_calculated",
        sum(1 for row in model_rows if str(row.get("transcript_cdna_context_status", ".")) == "calculated"),
    )

    transcript_status_counts = Counter(str(row.get("transcript_evidence_status", ".")) for row in model_rows)
    for status, count in sorted(transcript_status_counts.items()):
        add("transcript_evidence_status", status, "models", count)

    transcript_flag_counts: Counter[str] = Counter()
    for row in model_rows:
        flags = str(row.get("transcript_evidence_flags", "."))
        if flags and flags != ".":
            for flag in flags.split(","):
                if flag:
                    transcript_flag_counts[flag] += 1
    for flag, count in sorted(transcript_flag_counts.items()):
        add("transcript_evidence_flags", flag, "models", count)

    write_dict_tsv(rows, path, MODEL_SELECTION_SCORE_SUMMARY_COLUMNS)
