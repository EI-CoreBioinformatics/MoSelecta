"""Selection, backfill, and final-pruning stages for MoSelecta."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

from .constants import (
    FINAL_PRUNE_EXTRA_COLUMNS,
    MODEL_CONFIDENCE_LOW_DEPTH_5PLUS_THRESHOLD,
    MODEL_CONFIDENCE_NONSINGLETON_THRESHOLDS,
    MODEL_CONFIDENCE_VERY_LOW_DEPTH_7PLUS_THRESHOLD,
    SHADOW_PAIR_COLUMNS,
)
from .gffio import sanitize_id_part
from .log_utils import ProgressLogger
from .models import CodingModel, Sublocus
from .scoring import primary_model_selection_order_key, value_as_float
from .tsvio import write_dict_tsv


def interval_set_overlap_bp(intervals_a: list[tuple[int, int]], intervals_b: list[tuple[int, int]]) -> int:
    """Return total overlap between two sets of half-open intervals."""
    a = sorted(intervals_a)
    b = sorted(intervals_b)
    i = j = 0
    total = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if e > s:
            total += e - s
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return total


def cds_exon_overlap_metrics(model_a: CodingModel, model_b: CodingModel) -> dict[str, float | int]:
    overlap_bp = interval_set_overlap_bp(model_a.cds_parts, model_b.cds_parts)
    frac_a = overlap_bp / model_a.cds_length if model_a.cds_length else 0.0
    frac_b = overlap_bp / model_b.cds_length if model_b.cds_length else 0.0
    return {
        "cds_exon_overlap_bp": overlap_bp,
        "overlap_fraction_of_a": frac_a,
        "overlap_fraction_of_b": frac_b,
        "max_overlap_fraction": max(frac_a, frac_b),
        "min_overlap_fraction": min(frac_a, frac_b),
    }


def model_selection_row_sort_key(row: dict[str, object]) -> tuple:
    """Deterministic order used for sublocus and monosublocus winner ranking."""
    return primary_model_selection_order_key(row)


def select_sublocus_winners(
    subloci: list[Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
) -> list[CodingModel]:
    """Return one provisional winner per current sublocus."""
    winners: list[CodingModel] = []
    for sub in sorted(subloci, key=lambda s: (s.chrom, s.start0, s.end0, s.sublocus_id)):
        scored_models = [m for m in sub.models if m.output_transcript_id in model_structural_rows_by_transcript]
        if not scored_models:
            continue
        ordered = sorted(scored_models, key=lambda m: model_selection_row_sort_key(model_structural_rows_by_transcript[m.output_transcript_id]))
        winners.append(ordered[0])
    return winners


def stage4_selection_from_sublocus_winners(
    subloci: list[Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    *,
    min_cds_exon_overlap_fraction: float,
    progress: Optional[ProgressLogger] = None,
    progress_interval: int = 10000,
) -> dict[str, object]:
    """Select a non-overlapping set from provisional sublocus winners.

    The overlap graph is built on CDS-exon overlap among the provisional winners.
    Components of that graph are monosublocus holders. Within each holder, models
    are greedily retained in model selection score order; a lower-ranked model is
    suppressed only if it overlaps an already-retained model above the configured
    non-reciprocal threshold.
    """
    if progress:
        progress.log(f"Monosublocus de-overlap: selecting one provisional winner from {len(subloci):,} subloci")
    sub_by_model_id: dict[str, Sublocus] = {}
    for sub in subloci:
        for model in sub.models:
            sub_by_model_id[model.output_transcript_id] = sub

    winners = select_sublocus_winners(subloci, model_structural_rows_by_transcript)
    if progress:
        progress.log(f"Monosublocus de-overlap: selected {len(winners):,} sublocus winners; building CDS-exon overlap graph")
    model_by_id = {m.output_transcript_id: m for m in winners}
    adjacency: dict[str, set[str]] = {m.output_transcript_id: set() for m in winners}
    edge_by_pair: dict[tuple[str, str], dict[str, object]] = {}
    edge_rows: list[dict[str, object]] = []

    grouped: dict[tuple[str, str], list[CodingModel]] = defaultdict(list)
    for model in winners:
        grouped[(model.chrom, model.strand)].append(model)

    group_items = sorted(grouped.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    for group_number, ((_chrom, _strand), items) in enumerate(group_items, start=1):
        ordered = sorted(items, key=lambda m: (m.cds_start0, m.cds_end0, m.output_transcript_id))
        if progress:
            progress.log(
                f"Monosublocus de-overlap: scanning winner overlaps for {_chrom}:{_strand} "
                f"({group_number:,}/{len(group_items):,}; winners={len(ordered):,})",
                level="verbose",
            )
        group_edges_before = len(edge_rows)
        for i, a in enumerate(ordered, start=1):
            if progress:
                progress.every(i, len(ordered), f"Monosublocus de-overlap overlap scan {_chrom}:{_strand}", level="basic", interval=progress_interval)
            for b in ordered[i:]:
                if b.cds_start0 >= a.cds_end0:
                    break
                metrics = cds_exon_overlap_metrics(a, b)
                if int(metrics["cds_exon_overlap_bp"]) <= 0:
                    continue
                if float(metrics["max_overlap_fraction"]) + 1e-15 < min_cds_exon_overlap_fraction:
                    continue
                pair = tuple(sorted((a.output_transcript_id, b.output_transcript_id)))
                sub_a = sub_by_model_id.get(a.output_transcript_id)
                sub_b = sub_by_model_id.get(b.output_transcript_id)
                row = {
                    "model_a_id": a.output_transcript_id,
                    "model_b_id": b.output_transcript_id,
                    "source_a": a.source_name,
                    "source_b": b.source_name,
                    "sublocus_a": sub_a.sublocus_id if sub_a else ".",
                    "sublocus_b": sub_b.sublocus_id if sub_b else ".",
                    "chrom": a.chrom,
                    "strand": a.strand,
                    "model_a_cds_length": a.cds_length,
                    "model_b_cds_length": b.cds_length,
                    **metrics,
                }
                edge_by_pair[pair] = row
                edge_rows.append(row)
                adjacency[a.output_transcript_id].add(b.output_transcript_id)
                adjacency[b.output_transcript_id].add(a.output_transcript_id)

        if progress:
            progress.log(
                f"Monosublocus de-overlap: finished {_chrom}:{_strand}; added {len(edge_rows) - group_edges_before:,} overlap edges "
                f"(total edges={len(edge_rows):,})",
                level="verbose",
            )

    if progress:
        progress.log(f"Monosublocus de-overlap: overlap graph built with {len(edge_rows):,} edges; finding connected components")

    # Connected components over all winners, including singleton components.
    components: list[list[CodingModel]] = []
    visited: set[str] = set()
    sorted_winners_for_components = sorted(winners, key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.output_transcript_id))
    for winner_index, model in enumerate(sorted_winners_for_components, start=1):
        if progress:
            progress.every(winner_index, len(sorted_winners_for_components), "Monosublocus de-overlap connected-component seed scan", level="basic", interval=progress_interval)
        mid = model.output_transcript_id
        if mid in visited:
            continue
        stack = [mid]
        visited.add(mid)
        comp_ids: list[str] = []
        while stack:
            cur = stack.pop()
            comp_ids.append(cur)
            for nxt in sorted(adjacency[cur]):
                if nxt not in visited:
                    visited.add(nxt)
                    stack.append(nxt)
        components.append([model_by_id[x] for x in comp_ids])

    components.sort(key=lambda comp: (
        comp[0].chrom if comp else "",
        min(m.cds_start0 for m in comp),
        max(m.cds_end0 for m in comp),
        comp[0].strand if comp else ".",
        min(m.output_transcript_id for m in comp) if comp else "",
    ))

    monosublocus_rows: list[dict[str, object]] = []
    member_rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []
    suppressed_rows: list[dict[str, object]] = []
    selection_info_by_transcript: dict[str, dict[str, object]] = {}
    selected_models: list[CodingModel] = []

    if progress:
        progress.log(f"Monosublocus de-overlap: built {len(components):,} monosublocus holders; greedy de-overlap selection")
    component_index = 1
    for comp_number, comp in enumerate(components, start=1):
        if progress:
            progress.every(comp_number, len(components), "Monosublocus de-overlap greedy monosublocus selection", level="basic", interval=progress_interval)
        if not comp:
            continue
        chrom = comp[0].chrom
        strand = comp[0].strand
        strand_label = {"+": "plus", "-": "minus"}.get(strand, "unknown")
        start0 = min(m.cds_start0 for m in comp)
        end0 = max(m.cds_end0 for m in comp)
        monosublocus_id = f"monosublocus_{sanitize_id_part(chrom)}_{start0 + 1}_{end0}_{strand_label}.{component_index:06d}"
        component_index += 1

        comp_ids = {m.output_transcript_id for m in comp}
        # Retrieve only edges incident to members of this component.  The previous
        # implementation scanned the full edge_rows list for every component,
        # which is very expensive for large genomes with hundreds of thousands of
        # mostly-singleton monosublocus holders.  Using the adjacency graph keeps
        # this step proportional to the number of edges in the component.
        comp_edges_by_pair: dict[tuple[str, str], dict[str, object]] = {}
        for mid in comp_ids:
            for nbr in adjacency.get(mid, set()):
                if nbr not in comp_ids:
                    continue
                pair = tuple(sorted((mid, nbr)))
                edge = edge_by_pair.get(pair)
                if edge is not None:
                    comp_edges_by_pair[pair] = edge
        comp_edges = [comp_edges_by_pair[pair] for pair in sorted(comp_edges_by_pair)]
        for edge in comp_edges:
            edge["monosublocus_id"] = monosublocus_id
        ordered_members = sorted(comp, key=lambda m: model_selection_row_sort_key(model_structural_rows_by_transcript[m.output_transcript_id]))
        selected_in_comp: list[CodingModel] = []

        # Pre-fill rows after greedy selection so every member has consistent counts.
        temp_member_info: dict[str, dict[str, object]] = {}
        for member_order, model in enumerate(ordered_members, start=1):
            suppressors: list[tuple[float, CodingModel, dict[str, object]]] = []
            for selected_model in selected_in_comp:
                pair = tuple(sorted((model.output_transcript_id, selected_model.output_transcript_id)))
                edge = edge_by_pair.get(pair)
                if not edge:
                    continue
                suppressors.append((float(edge.get("max_overlap_fraction", 0.0)), selected_model, edge))
            if suppressors:
                suppressors.sort(key=lambda x: (-x[0], model_selection_row_sort_key(model_structural_rows_by_transcript[x[1].output_transcript_id])))
                _max_frac, suppressor, edge = suppressors[0]
                if edge["model_a_id"] == model.output_transcript_id:
                    suppressed_fraction = edge["overlap_fraction_of_a"]
                    selected_fraction = edge["overlap_fraction_of_b"]
                else:
                    suppressed_fraction = edge["overlap_fraction_of_b"]
                    selected_fraction = edge["overlap_fraction_of_a"]
                temp_member_info[model.output_transcript_id] = {
                    "monosublocus_id": monosublocus_id,
                    "monosublocus_member_order_rank": member_order,
                    "monosublocus_selected": "false",
                    "monosublocus_selection_order": ".",
                    "monosublocus_suppressed_by": suppressor.output_transcript_id,
                    "monosublocus_suppression_reason": "overlaps_higher_ranked_selected_model",
                    "monosublocus_suppression_overlap_bp": edge["cds_exon_overlap_bp"],
                    "monosublocus_suppressed_model_overlap_fraction": suppressed_fraction,
                    "monosublocus_selected_model_overlap_fraction": selected_fraction,
                    "monosublocus_suppression_max_overlap_fraction": edge["max_overlap_fraction"],
                }
            else:
                selected_in_comp.append(model)
                temp_member_info[model.output_transcript_id] = {
                    "monosublocus_id": monosublocus_id,
                    "monosublocus_member_order_rank": member_order,
                    "monosublocus_selected": "true",
                    "monosublocus_selection_order": len(selected_in_comp),
                    "monosublocus_suppressed_by": ".",
                    "monosublocus_suppression_reason": ".",
                    "monosublocus_suppression_overlap_bp": ".",
                    "monosublocus_suppressed_model_overlap_fraction": ".",
                    "monosublocus_selected_model_overlap_fraction": ".",
                    "monosublocus_suppression_max_overlap_fraction": ".",
                }

        not_all_pairwise = len(comp) > 2 and len(comp_edges) < (len(comp) * (len(comp) - 1) // 2)
        monosublocus_row = {
            "monosublocus_id": monosublocus_id,
            "chrom": chrom,
            "start": start0 + 1,
            "end": end0,
            "strand": strand,
            "member_count": len(comp),
            "selected_count": len(selected_in_comp),
            "suppressed_count": len(comp) - len(selected_in_comp),
            "edge_count": len(comp_edges),
            "min_cds_exon_overlap_fraction": min_cds_exon_overlap_fraction,
            "not_all_pairwise_overlapping": "true" if not_all_pairwise else "false",
            "sources": ",".join(sorted({m.source_name for m in comp})),
            "subloci": ",".join(sorted({sub_by_model_id[m.output_transcript_id].sublocus_id for m in comp if m.output_transcript_id in sub_by_model_id})),
            "selected_model_ids": ",".join(m.output_transcript_id for m in selected_in_comp),
            "suppressed_model_ids": ",".join(m.output_transcript_id for m in comp if temp_member_info[m.output_transcript_id]["monosublocus_selected"] == "false") or ".",
        }
        monosublocus_rows.append(monosublocus_row)

        for model in ordered_members:
            sub = sub_by_model_id.get(model.output_transcript_id)
            row = model_structural_rows_by_transcript[model.output_transcript_id]
            info = temp_member_info[model.output_transcript_id]
            full = {
                "monosublocus_id": monosublocus_id,
                "superlocus_id": sub.superlocus_id if sub else ".",
                "sublocus_id": sub.sublocus_id if sub else ".",
                "source_name": model.source_name,
                "input_gene_id": model.input_gene_id,
                "input_transcript_id": model.input_transcript_id,
                "output_gene_id": model.output_gene_id,
                "output_transcript_id": model.output_transcript_id,
                "chrom": model.chrom,
                "cds_start": model.cds_start0 + 1,
                "cds_end": model.cds_end0,
                "strand": model.strand,
                "cds_length": model.cds_length,
                "model_selection_score": row.get("model_selection_score", "."),
                "model_selection_score_scaled": row.get("model_selection_score_scaled", "."),
                "model_selection_score_rank": row.get("model_selection_score_rank", "."),
                "model_selection_order_rank": row.get("model_selection_order_rank", "."),
                "exact_cds_exon_chain_support_model_count": row.get("exact_cds_exon_chain_support_model_count", "."),
                "exact_cds_exon_chain_support_source_count": row.get("exact_cds_exon_chain_support_source_count", "."),
                "exact_cds_exon_chain_support_sources": row.get("exact_cds_exon_chain_support_sources", "."),
                "exact_cds_exon_chain_support_model_ids": row.get("exact_cds_exon_chain_support_model_ids", "."),
                "exact_cds_exon_chain_support_fraction_local": row.get("exact_cds_exon_chain_support_fraction_local", "."),
                "exact_structure_support_score": row.get("exact_structure_support_score", "."),
                "model_completeness_score": row.get("model_completeness_score", "."),
                "intrinsic_penalty_score": row.get("intrinsic_penalty_score", "."),
                "consensus_similarity_score": row.get("consensus_similarity_score", "."),
                "structural_support_score": row.get("structural_support_score", "."),
                "pairwise_similarity_score": row.get("pairwise_similarity_score", "."),
                "source_prior_score": row.get("source_prior_score", "."),
                "low_source_conflict_flag": row.get("low_source_conflict_flag", "."),
                "low_source_conflict_agreement_factor": row.get("low_source_conflict_agreement_factor", "."),
                "low_source_conflict_reason": row.get("low_source_conflict_reason", "."),
                **info,
            }
            member_rows.append(full)
            selection_info_by_transcript[model.output_transcript_id] = {
                **full,
                "monosublocus_member_count": len(comp),
                "monosublocus_selected_count": len(selected_in_comp),
                "monosublocus_suppressed_count": len(comp) - len(selected_in_comp),
            }
            if info["monosublocus_selected"] == "true":
                selected_rows.append(full)
                selected_models.append(model)
            else:
                suppressed_rows.append(full)

    # Ensure every overlap edge carries its monosublocus ID.  Edges are assigned
    # during component processing above; this fallback is kept for robustness if a
    # future code path appends edges before component IDs are available.
    mono_by_member: dict[str, str] = {row["output_transcript_id"]: row["monosublocus_id"] for row in member_rows}
    for row in edge_rows:
        if row.get("monosublocus_id", ".") == ".":
            row["monosublocus_id"] = mono_by_member.get(str(row["model_a_id"]), ".")

    if progress:
        progress.log(
            f"Monosublocus de-overlap complete: winners={len(winners):,}; monosubloci={len(monosublocus_rows):,}; "
            f"selected={len(selected_rows):,}; suppressed={len(suppressed_rows):,}"
        )

    return {
        "sublocus_winners": winners,
        "monosublocus_rows": monosublocus_rows,
        "monosublocus_member_rows": member_rows,
        "monosublocus_edge_rows": edge_rows,
        "selected_rows": selected_rows,
        "suppressed_rows": suppressed_rows,
        "selected_models": selected_models,
        "selection_info_by_transcript": selection_info_by_transcript,
    }


def shared_cds_splice_site_count(model_a: CodingModel, model_b: CodingModel) -> int:
    """Return count of shared CDS splice-site coordinates between two models."""
    if not (model_a.is_multi_exon and model_b.is_multi_exon):
        return 0
    return len(model_a.cds_splice_sites & model_b.cds_splice_sites)


def selected_model_conflicts_for_backfill(
    candidate: CodingModel,
    selected_models: list[CodingModel],
    *,
    min_cds_exon_overlap_fraction: float,
) -> list[dict[str, object]]:
    """Return every selected-model conflict for a backfill candidate.

    Backfill is deliberately superlocus-local; callers pass only currently selected
    models from the same superlocus. A candidate is rejected if either:
      * it is multi-exon and shares any CDS splice site with an already-selected
        multi-exon model, or
      * it overlaps any already-selected model by the non-reciprocal CDS-exon
        overlap threshold used for monosublocus selection.
    """
    conflicts: list[dict[str, object]] = []
    for selected in selected_models:
        splice_count = shared_cds_splice_site_count(candidate, selected)
        metrics = cds_exon_overlap_metrics(candidate, selected)
        overlap_bp = int(metrics["cds_exon_overlap_bp"])
        max_overlap = float(metrics["max_overlap_fraction"])
        overlap_conflict = overlap_bp > 0 and max_overlap + 1e-15 >= min_cds_exon_overlap_fraction
        if splice_count > 0 or overlap_conflict:
            if splice_count > 0 and overlap_conflict:
                reason = "shared_cds_splice_site_and_cds_exon_overlap_with_selected_model"
                priority = 0
            elif splice_count > 0:
                reason = "shared_cds_splice_site_with_selected_model"
                priority = 1
            else:
                reason = "cds_exon_overlap_with_selected_model"
                priority = 2
            conflicts.append({
                "priority": priority,
                "backfill_rejection_reason": reason,
                "backfill_rejected_by": selected.output_transcript_id,
                "backfill_rejected_by_source": selected.source_name,
                "backfill_conflict_shared_splice_site_count": splice_count,
                "backfill_conflict_overlap_bp": overlap_bp,
                "backfill_candidate_overlap_fraction": metrics["overlap_fraction_of_a"],
                "backfill_selected_model_overlap_fraction": metrics["overlap_fraction_of_b"],
                "backfill_conflict_max_overlap_fraction": metrics["max_overlap_fraction"],
            })
    conflicts.sort(key=lambda row: (
        int(row["priority"]),
        -int(row.get("backfill_conflict_shared_splice_site_count", 0)),
        -float(row.get("backfill_conflict_max_overlap_fraction", 0.0)),
        str(row.get("backfill_rejected_by", "")),
    ))
    for conflict in conflicts:
        conflict.pop("priority", None)
    return conflicts


def selected_model_conflict_for_backfill(
    candidate: CodingModel,
    selected_models: list[CodingModel],
    *,
    min_cds_exon_overlap_fraction: float,
) -> Optional[dict[str, object]]:
    """Return the strongest selected-model conflict for legacy backfill logic."""

    conflicts = selected_model_conflicts_for_backfill(
        candidate,
        selected_models,
        min_cds_exon_overlap_fraction=min_cds_exon_overlap_fraction,
    )
    return conflicts[0] if conflicts else None


def strict_r1_replacement_eligibility(
    candidate: CodingModel,
    conflicts: list[dict[str, object]],
    *,
    model_by_id: dict[str, CodingModel],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    monosublocus_selected_ids: set[str],
    allow_multiexon_blockers: bool = False,
) -> tuple[bool, str, list[CodingModel]]:
    """Evaluate the exact strict-R1 replacement rule without changing state."""

    if not candidate.is_multi_exon:
        return False, "candidate_not_multi_exon", []
    candidate_row = model_structural_rows_by_transcript.get(candidate.output_transcript_id, {})
    if int(value_as_float(candidate_row.get("corroborating_other_source_count"))) < 1:
        return False, "candidate_corroboration_count_lt_1", []
    if len(conflicts) not in {1, 2}:
        reason = "blocking_model_count_gt_2" if len(conflicts) > 2 else "blocking_model_count_not_1_or_2"
        return False, reason, []

    blockers: list[CodingModel] = []
    for conflict in conflicts:
        blocker_id = str(conflict.get("backfill_rejected_by", "."))
        blocker = model_by_id.get(blocker_id)
        if blocker is None:
            return False, "blocking_model_not_found", []
        blockers.append(blocker)
    if not allow_multiexon_blockers and any(not blocker.is_single_exon for blocker in blockers):
        return False, "blocking_model_not_single_exon", blockers
    if any(blocker.output_transcript_id not in monosublocus_selected_ids for blocker in blockers):
        return False, "blocking_model_not_monosublocus_selected", blockers
    if any(
        int(value_as_float(model_structural_rows_by_transcript.get(blocker.output_transcript_id, {}).get("corroborating_other_source_count"))) >= 1
        for blocker in blockers
    ):
        return False, "blocking_model_corroboration_count_ge_1", blockers

    candidate_score = value_as_float(candidate_row.get("model_selection_score"))
    blocker_scores = [
        value_as_float(model_structural_rows_by_transcript.get(blocker.output_transcript_id, {}).get("model_selection_score"))
        for blocker in blockers
    ]
    if any(candidate_score <= blocker_score for blocker_score in blocker_scores):
        return False, "candidate_score_not_strictly_greater_than_every_blocker", blockers
    return True, "strict_r1_eligible", blockers


def stage_model_base_row(
    model: CodingModel,
    sub_by_model_id: dict[str, Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Common row fields for selection-stage TSV outputs."""
    sub = sub_by_model_id.get(model.output_transcript_id)
    row = model_structural_rows_by_transcript.get(model.output_transcript_id, {})
    return {
        "superlocus_id": sub.superlocus_id if sub else ".",
        "sublocus_id": sub.sublocus_id if sub else ".",
        "sublocus_type": sub.sublocus_type if sub else ".",
        "source_name": model.source_name,
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
        "low_source_conflict_flag": row.get("low_source_conflict_flag", "."),
        "low_source_conflict_agreement_factor": row.get("low_source_conflict_agreement_factor", "."),
        "low_source_conflict_reason": row.get("low_source_conflict_reason", "."),
        "model_structural_support_score": row.get("model_structural_support_score", "."),
        "model_structural_support_score_scaled": row.get("model_structural_support_score_scaled", "."),
        "corroborating_other_source_count": row.get("corroborating_other_source_count", "."),
        "run_possible_other_source_count": row.get("run_possible_other_source_count", "."),
        "structural_corroboration_score": row.get("structural_corroboration_score", "."),
        "structural_corroboration_alpha": row.get("structural_corroboration_alpha", "."),
        "structural_corroboration_weight": row.get("structural_corroboration_weight", "."),
        "structural_corroboration_contribution": row.get("structural_corroboration_contribution", "."),
        "structural_corroboration_in_selection": row.get("structural_corroboration_in_selection", "."),
        "sequence_integrity_status": row.get("sequence_integrity_status", "."),
        "sequence_integrity_penalty_score": row.get("sequence_integrity_penalty_score", "."),
        "sequence_integrity_penalty_weight": row.get("sequence_integrity_penalty_weight", "."),
        "sequence_integrity_penalty_in_selection": row.get("sequence_integrity_penalty_in_selection", "."),
        "sequence_integrity_flags": row.get("sequence_integrity_flags", "."),
        "protein_evidence_score": row.get("protein_evidence_score", "."),
        "protein_evidence_score_scaled": row.get("protein_evidence_score_scaled", "."),
        "protein_evidence_score_weight": row.get("protein_evidence_score_weight", "."),
        "protein_evidence_score_in_selection": row.get("protein_evidence_score_in_selection", "."),
        "transcript_evidence_score": row.get("transcript_evidence_score", "."),
        "transcript_evidence_score_scaled": row.get("transcript_evidence_score_scaled", "."),
        "transcript_evidence_score_context_adjusted": row.get("transcript_evidence_score_context_adjusted", "."),
        "transcript_evidence_score_context_adjusted_scaled": row.get("transcript_evidence_score_context_adjusted_scaled", "."),
        "transcript_cdna_context_status": row.get("transcript_cdna_context_status", "."),
        "transcript_cdna_context_cds_to_min_cdna_fraction": row.get("transcript_cdna_context_cds_to_min_cdna_fraction", "."),
        "transcript_cdna_context_flags": row.get("transcript_cdna_context_flags", "."),
        "transcript_cdna_context_evaluable": row.get("transcript_cdna_context_evaluable", "."),
        "transcript_cdna_context_concern_score": row.get("transcript_cdna_context_concern_score", "."),
        "transcript_cdna_context_concern_class": row.get("transcript_cdna_context_concern_class", "."),
        "transcript_cdna_context_pruning_concern": row.get("transcript_cdna_context_pruning_concern", "."),
        "transcript_cdna_context_pruning_concern_class": row.get("transcript_cdna_context_pruning_concern_class", "."),
        "transcript_cdna_context_pruning_supporting_flags": row.get("transcript_cdna_context_pruning_supporting_flags", "."),
        "transcript_cdna_context_penalty": row.get("transcript_cdna_context_penalty", "."),
        "transcript_evidence_score_weight": row.get("transcript_evidence_score_weight", "."),
        "transcript_evidence_score_in_selection": row.get("transcript_evidence_score_in_selection", "."),
        "model_selection_score": row.get("model_selection_score", "."),
        "model_selection_score_scaled": row.get("model_selection_score_scaled", "."),
        "model_selection_score_rank": row.get("model_selection_score_rank", "."),
        "model_selection_order_rank": row.get("model_selection_order_rank", "."),
        "exact_cds_exon_chain_support_model_count": row.get("exact_cds_exon_chain_support_model_count", "."),
        "exact_cds_exon_chain_support_source_count": row.get("exact_cds_exon_chain_support_source_count", "."),
        "exact_cds_exon_chain_support_sources": row.get("exact_cds_exon_chain_support_sources", "."),
        "exact_cds_exon_chain_support_model_ids": row.get("exact_cds_exon_chain_support_model_ids", "."),
        "exact_cds_exon_chain_support_fraction_local": row.get("exact_cds_exon_chain_support_fraction_local", "."),
        "exact_structure_support_score": row.get("exact_structure_support_score", "."),
        "model_completeness_score": row.get("model_completeness_score", "."),
        "intrinsic_penalty_score": row.get("intrinsic_penalty_score", "."),
        "consensus_similarity_score": row.get("consensus_similarity_score", "."),
        "structural_support_score": row.get("structural_support_score", "."),
        "pairwise_similarity_score": row.get("pairwise_similarity_score", "."),
        "source_prior_score": row.get("source_prior_score", "."),
    }


def stage5_backfill_selection(
    subloci: list[Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    stage4: dict[str, object],
    *,
    min_cds_exon_overlap_fraction: float,
    strict_r1_backfill_replacement: bool = False,
    strict_r1_allow_multiexon_blockers: bool = False,
    structural_corroboration_weight: float = 0.0,
    structural_corroboration_alpha: float = 16.0,
    sequence_penalty_missing_start: float = 0.30,
    sequence_penalty_missing_stop: float = 0.50,
    progress: Optional[ProgressLogger] = None,
    progress_interval: int = 10000,
) -> dict[str, object]:
    """Backfill non-conflicting candidate models after monosublocus selection.

    This stage is superlocus-local. It starts with the models retained by Monosublocus de-overlap
    monosublocus selection, permanently removes Monosublocus de-overlap-suppressed winners from
    consideration, and greedily scans all remaining scored models in each
    superlocus by model-selection rank. A candidate is added only if it
    does not conflict with any currently selected model in that same superlocus.
    """
    if strict_r1_allow_multiexon_blockers and not strict_r1_backfill_replacement:
        raise ValueError(
            "strict_r1_allow_multiexon_blockers requires "
            "strict_r1_backfill_replacement"
        )
    if progress:
        progress.log(f"Backfill rescue: preparing backfill over {len(subloci):,} subloci")
    sub_by_model_id: dict[str, Sublocus] = {m.output_transcript_id: sub for sub in subloci for m in sub.models}
    all_model_by_id: dict[str, CodingModel] = {m.output_transcript_id: m for sub in subloci for m in sub.models}
    all_models_by_superlocus: dict[str, list[CodingModel]] = defaultdict(list)
    for sub in subloci:
        for model in sub.models:
            if model.output_transcript_id in model_structural_rows_by_transcript:
                all_models_by_superlocus[sub.superlocus_id].append(model)

    stage4_selected_models: list[CodingModel] = list(stage4.get("selected_models", []))  # type: ignore[arg-type]
    selected_ids: set[str] = {m.output_transcript_id for m in stage4_selected_models}
    monosublocus_selected_ids = set(selected_ids)
    suppressed_rows: list[dict[str, object]] = list(stage4.get("suppressed_rows", []))  # type: ignore[arg-type]
    suppressed_ids: set[str] = {str(row.get("output_transcript_id")) for row in suppressed_rows}
    stage4_selection_info: dict[str, dict[str, object]] = dict(stage4.get("selection_info_by_transcript", {}))  # type: ignore[arg-type]

    selected_by_superlocus: dict[str, list[CodingModel]] = defaultdict(list)
    for model in stage4_selected_models:
        sub = sub_by_model_id.get(model.output_transcript_id)
        if sub:
            selected_by_superlocus[sub.superlocus_id].append(model)

    candidate_rows: list[dict[str, object]] = []
    backfilled_rows: list[dict[str, object]] = []
    rejected_rows: list[dict[str, object]] = []
    final_rows: list[dict[str, object]] = []
    final_models: list[CodingModel] = []
    final_selection_info_by_transcript: dict[str, dict[str, object]] = {}
    r1_replacement_rows: list[dict[str, object]] = []
    r1_rejection_counts: Counter[str] = Counter()
    r1_candidates_evaluated = 0
    r1_replaced_ids: set[str] = set()

    # Seed final set with Monosublocus de-overlap selected models.
    for model in sorted(stage4_selected_models, key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.output_transcript_id)):
        base = stage_model_base_row(model, sub_by_model_id, model_structural_rows_by_transcript)
        info = dict(stage4_selection_info.get(model.output_transcript_id, {}))
        final_info = {
            **info,
            "final_selected": "true",
            "final_selection_stage": "monosublocus",
            "backfill_candidate_order_rank": ".",
            "backfill_selected": ".",
            "backfill_selection_order": ".",
            "backfill_rejection_reason": ".",
            "backfill_rejected_by": ".",
            "backfill_rejected_by_source": ".",
            "backfill_conflict_shared_splice_site_count": ".",
            "backfill_conflict_overlap_bp": ".",
            "backfill_candidate_overlap_fraction": ".",
            "backfill_selected_model_overlap_fraction": ".",
            "backfill_conflict_max_overlap_fraction": ".",
        }
        full = {**base, **final_info}
        final_rows.append(full)
        final_models.append(model)
        final_selection_info_by_transcript[model.output_transcript_id] = {**full, **final_info}

    superlocus_ids_for_backfill = sorted(all_models_by_superlocus)
    if progress:
        progress.log(f"Backfill rescue: scanning {len(superlocus_ids_for_backfill):,} superloci for backfill candidates")
    global_backfill_selection_order = 1
    for superlocus_number, superlocus_id in enumerate(superlocus_ids_for_backfill, start=1):
        if progress:
            progress.every(superlocus_number, len(superlocus_ids_for_backfill), "Backfill rescue backfill superlocus scan", level="basic", interval=progress_interval)
        current_selected = selected_by_superlocus.get(superlocus_id, [])
        candidates = [
            model for model in all_models_by_superlocus[superlocus_id]
            if model.output_transcript_id not in selected_ids
            and model.output_transcript_id not in suppressed_ids
            and model.output_transcript_id in model_structural_rows_by_transcript
        ]
        candidates.sort(key=lambda m: model_selection_row_sort_key(model_structural_rows_by_transcript[m.output_transcript_id]))
        if progress and candidates:
            progress.log(
                f"Backfill rescue: {superlocus_id} candidates={len(candidates):,}; starting selected={len(current_selected):,}",
                level="debug",
            )
        for candidate_order, model in enumerate(candidates, start=1):
            base = stage_model_base_row(model, sub_by_model_id, model_structural_rows_by_transcript)
            conflicts = selected_model_conflicts_for_backfill(
                model,
                current_selected,
                min_cds_exon_overlap_fraction=min_cds_exon_overlap_fraction,
            )
            conflict = conflicts[0] if conflicts else None
            common_info = {
                "final_selected": "false",
                "final_selection_stage": ".",
                "backfill_candidate_order_rank": candidate_order,
                "backfill_selected": "false",
                "backfill_selection_order": ".",
            }
            if conflict is None:
                accept_info = {
                    **common_info,
                    "final_selected": "true",
                    "final_selection_stage": "backfill",
                    "backfill_selected": "true",
                    "backfill_selection_order": global_backfill_selection_order,
                    "backfill_rejection_reason": ".",
                    "backfill_rejected_by": ".",
                    "backfill_rejected_by_source": ".",
                    "backfill_conflict_shared_splice_site_count": ".",
                    "backfill_conflict_overlap_bp": ".",
                    "backfill_candidate_overlap_fraction": ".",
                    "backfill_selected_model_overlap_fraction": ".",
                    "backfill_conflict_max_overlap_fraction": ".",
                }
                full = {**base, **accept_info}
                candidate_rows.append(full)
                backfilled_rows.append(full)
                final_rows.append(full)
                final_models.append(model)
                final_selection_info_by_transcript[model.output_transcript_id] = full
                current_selected.append(model)
                selected_ids.add(model.output_transcript_id)
                global_backfill_selection_order += 1
            else:
                if strict_r1_backfill_replacement:
                    r1_candidates_evaluated += 1
                    eligible, r1_reason, blockers = strict_r1_replacement_eligibility(
                        model,
                        conflicts,
                        model_by_id=all_model_by_id,
                        model_structural_rows_by_transcript=model_structural_rows_by_transcript,
                        monosublocus_selected_ids=monosublocus_selected_ids,
                        allow_multiexon_blockers=strict_r1_allow_multiexon_blockers,
                    )
                    if eligible:
                        blocker_ids = {blocker.output_transcript_id for blocker in blockers}
                        current_selected[:] = [
                            selected for selected in current_selected
                            if selected.output_transcript_id not in blocker_ids
                        ]
                        r1_replaced_ids.update(blocker_ids)
                        for blocker_id in blocker_ids:
                            final_selection_info_by_transcript.pop(blocker_id, None)

                        uses_multiexon_blocker = any(blocker.is_multi_exon for blocker in blockers)
                        eligibility_class = (
                            "extended_multiexon_blocker_eligibility"
                            if uses_multiexon_blocker
                            else "original_single_exon_blocker_eligibility"
                        )
                        replacement_reason = (
                            "strict_r1me_supported_multiexon_replaces_uncorroborated_holder_winner"
                            if uses_multiexon_blocker
                            else "strict_r1_supported_multiexon_replaces_uncorroborated_single_exon_winner"
                        )
                        accept_info = {
                            **common_info,
                            "final_selected": "true",
                            "final_selection_stage": "backfill_r1_replacement",
                            "backfill_selected": "true",
                            "backfill_selection_order": global_backfill_selection_order,
                            "backfill_rejection_reason": ".",
                            "backfill_rejected_by": ".",
                            "backfill_rejected_by_source": ".",
                            "backfill_conflict_shared_splice_site_count": ".",
                            "backfill_conflict_overlap_bp": ".",
                            "backfill_candidate_overlap_fraction": ".",
                            "backfill_selected_model_overlap_fraction": ".",
                            "backfill_conflict_max_overlap_fraction": ".",
                            "r1_replacement_attempted": "true",
                            "r1_replacement_accepted": "true",
                            "r1_replacement_reason": replacement_reason,
                            "r1_replaced_model_ids": ",".join(sorted(blocker_ids)),
                        }
                        full = {**base, **accept_info}
                        candidate_rows.append(full)
                        backfilled_rows.append(full)
                        final_rows.append(full)
                        final_models.append(model)
                        final_selection_info_by_transcript[model.output_transcript_id] = full
                        current_selected.append(model)
                        selected_ids.add(model.output_transcript_id)
                        global_backfill_selection_order += 1

                        candidate_row = model_structural_rows_by_transcript.get(model.output_transcript_id, {})
                        blocker_rows = [
                            model_structural_rows_by_transcript.get(blocker.output_transcript_id, {})
                            for blocker in blockers
                        ]
                        conflict_by_blocker = {
                            str(item.get("backfill_rejected_by", ".")): item for item in conflicts
                        }
                        r1_replacement_rows.append({
                            "superlocus_id": superlocus_id,
                            "candidate_model_id": model.output_transcript_id,
                            "candidate_source": model.source_name,
                            "candidate_sublocus_id": base.get("sublocus_id", "."),
                            "candidate_cds_exon_count": model.cds_exon_count,
                            "candidate_model_selection_score": candidate_row.get("model_selection_score", "."),
                            "candidate_corroborating_other_source_count": candidate_row.get("corroborating_other_source_count", "."),
                            "replaced_model_count": len(blockers),
                            "replaced_model_ids": ",".join(blocker.output_transcript_id for blocker in blockers),
                            "replaced_model_sources": ",".join(blocker.source_name for blocker in blockers),
                            "replaced_model_cds_exon_counts": ",".join(
                                str(blocker.cds_exon_count) for blocker in blockers
                            ),
                            "replaced_model_exon_classes": ",".join(
                                "single_exon" if blocker.is_single_exon else "multi_exon"
                                for blocker in blockers
                            ),
                            "replaced_model_sublocus_ids": ",".join(
                                str(sub_by_model_id.get(blocker.output_transcript_id).sublocus_id)
                                if sub_by_model_id.get(blocker.output_transcript_id) else "."
                                for blocker in blockers
                            ),
                            "replaced_model_selection_scores": ",".join(
                                str(row.get("model_selection_score", ".")) for row in blocker_rows
                            ),
                            "replaced_model_corroborating_other_source_counts": ",".join(
                                str(row.get("corroborating_other_source_count", ".")) for row in blocker_rows
                            ),
                            "conflict_reasons": ";".join(
                                str(conflict_by_blocker[blocker.output_transcript_id].get("backfill_rejection_reason", "."))
                                for blocker in blockers
                            ),
                            "conflict_overlap_bp": ",".join(
                                str(conflict_by_blocker[blocker.output_transcript_id].get("backfill_conflict_overlap_bp", "."))
                                for blocker in blockers
                            ),
                            "conflict_candidate_overlap_fractions": ",".join(
                                str(conflict_by_blocker[blocker.output_transcript_id].get("backfill_candidate_overlap_fraction", "."))
                                for blocker in blockers
                            ),
                            "conflict_replaced_model_overlap_fractions": ",".join(
                                str(conflict_by_blocker[blocker.output_transcript_id].get("backfill_selected_model_overlap_fraction", "."))
                                for blocker in blockers
                            ),
                            "conflict_shared_splice_site_counts": ",".join(
                                str(conflict_by_blocker[blocker.output_transcript_id].get("backfill_conflict_shared_splice_site_count", "."))
                                for blocker in blockers
                            ),
                            "structural_corroboration_weight": structural_corroboration_weight,
                            "structural_corroboration_alpha": structural_corroboration_alpha,
                            "strict_r1_backfill_replacement": "true",
                            "strict_r1_allow_multiexon_blockers": (
                                "true" if strict_r1_allow_multiexon_blockers else "false"
                            ),
                            "sequence_penalty_missing_start": sequence_penalty_missing_start,
                            "sequence_penalty_missing_stop": sequence_penalty_missing_stop,
                            "r1_eligibility_class": eligibility_class,
                            "r1_replacement_reason": accept_info["r1_replacement_reason"],
                        })
                        continue
                    r1_rejection_counts[r1_reason] += 1
                reject_info = {
                    **common_info,
                    **conflict,
                }
                if strict_r1_backfill_replacement:
                    reject_info.update({
                        "r1_replacement_attempted": "true",
                        "r1_replacement_accepted": "false",
                        "r1_replacement_reason": r1_reason,
                        "r1_replaced_model_ids": ".",
                    })
                full = {**base, **reject_info}
                candidate_rows.append(full)
                rejected_rows.append(full)

    if r1_replaced_ids:
        final_rows = [
            row for row in final_rows
            if str(row.get("output_transcript_id")) not in r1_replaced_ids
        ]
        final_models = [
            model for model in final_models
            if model.output_transcript_id not in r1_replaced_ids
        ]

    final_rows.sort(key=lambda row: (
        str(row.get("chrom", "")),
        int(value_as_float(row.get("cds_start"))),
        int(value_as_float(row.get("cds_end"))),
        str(row.get("strand", ".")),
        str(row.get("output_transcript_id", "")),
    ))
    final_models.sort(key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.strand, m.output_transcript_id))

    if progress:
        progress.log(
            f"Backfill rescue complete: candidates={len(candidate_rows):,}; backfilled={len(backfilled_rows):,}; "
            f"rejected={len(rejected_rows):,}; final={len(final_rows):,}"
        )

    one_blocker_replacements = sum(
        1 for row in r1_replacement_rows if int(value_as_float(row.get("replaced_model_count"))) == 1
    )
    two_blocker_replacements = sum(
        1 for row in r1_replacement_rows if int(value_as_float(row.get("replaced_model_count"))) == 2
    )
    r1_summary_rows = [
        {"metric": "r1_enabled", "value": "true" if strict_r1_backfill_replacement else "false", "details": "."},
        {"metric": "candidates_evaluated", "value": r1_candidates_evaluated, "details": "."},
        {"metric": "replacements_accepted", "value": len(r1_replacement_rows), "details": "."},
        {"metric": "one_blocker_replacements", "value": one_blocker_replacements, "details": "."},
        {"metric": "two_blocker_replacements", "value": two_blocker_replacements, "details": "."},
    ]
    if strict_r1_allow_multiexon_blockers:
        r1_summary_rows.extend([
            {"metric": "r1_allow_multiexon_blockers", "value": "true", "details": "."},
            {
                "metric": "original_eligibility_replacements",
                "value": sum(
                    1 for row in r1_replacement_rows
                    if row.get("r1_eligibility_class") == "original_single_exon_blocker_eligibility"
                ),
                "details": ".",
            },
            {
                "metric": "extended_multiexon_eligibility_replacements",
                "value": sum(
                    1 for row in r1_replacement_rows
                    if row.get("r1_eligibility_class") == "extended_multiexon_blocker_eligibility"
                ),
                "details": ".",
            },
            {
                "metric": "single_exon_blockers_replaced",
                "value": sum(
                    str(row.get("replaced_model_exon_classes", "")).split(",").count("single_exon")
                    for row in r1_replacement_rows
                ),
                "details": ".",
            },
            {
                "metric": "multi_exon_blockers_replaced",
                "value": sum(
                    str(row.get("replaced_model_exon_classes", "")).split(",").count("multi_exon")
                    for row in r1_replacement_rows
                ),
                "details": ".",
            },
        ])
    r1_summary_rows.extend(
        {"metric": "rejections", "value": count, "details": reason}
        for reason, count in sorted(r1_rejection_counts.items())
    )

    return {
        "backfill_candidate_rows": candidate_rows,
        "backfilled_rows": backfilled_rows,
        "backfill_rejected_rows": rejected_rows,
        "final_rows": final_rows,
        "final_models": final_models,
        "final_selection_info_by_transcript": final_selection_info_by_transcript,
        "r1_replacement_rows": r1_replacement_rows,
        "r1_summary_rows": r1_summary_rows,
        "r1_replaced_model_ids": sorted(r1_replaced_ids),
        "r1_allow_multiexon_blockers": strict_r1_allow_multiexon_blockers,
    }


def write_backfill_outputs(
    stage5: dict[str, object],
    *,
    backfill_candidates_tsv: Optional[Path],
    backfilled_models_tsv: Optional[Path],
    backfill_rejected_models_tsv: Optional[Path],
    final_selected_models_tsv: Optional[Path],
) -> None:
    stage_rows = list(stage5.get("backfill_candidate_rows", []))  # type: ignore[arg-type]
    include_corroboration = any(
        value_as_float(row.get("structural_corroboration_weight")) > 0 for row in stage_rows
    )
    include_r1 = any(str(row.get("r1_replacement_attempted", ".")) == "true" for row in stage_rows)
    columns = selection_stage_output_columns(
        include_structural_corroboration=include_corroboration,
        include_r1=include_r1,
    )
    if backfill_candidates_tsv is not None:
        write_dict_tsv(stage5["backfill_candidate_rows"], backfill_candidates_tsv, columns)  # type: ignore[arg-type]
    if backfilled_models_tsv is not None:
        write_dict_tsv(stage5["backfilled_rows"], backfilled_models_tsv, columns)  # type: ignore[arg-type]
    if backfill_rejected_models_tsv is not None:
        write_dict_tsv(stage5["backfill_rejected_rows"], backfill_rejected_models_tsv, columns)  # type: ignore[arg-type]
    if final_selected_models_tsv is not None:
        write_dict_tsv(stage5["final_rows"], final_selected_models_tsv, columns)  # type: ignore[arg-type]


R1_REPLACEMENT_COLUMNS = [
    "superlocus_id", "candidate_model_id", "candidate_source", "candidate_sublocus_id",
    "candidate_model_selection_score", "candidate_corroborating_other_source_count",
    "replaced_model_count", "replaced_model_ids", "replaced_model_sources",
    "replaced_model_sublocus_ids", "replaced_model_selection_scores",
    "replaced_model_corroborating_other_source_counts", "conflict_reasons",
    "conflict_overlap_bp", "conflict_candidate_overlap_fractions",
    "conflict_replaced_model_overlap_fractions", "conflict_shared_splice_site_counts",
    "structural_corroboration_weight", "structural_corroboration_alpha",
    "r1_replacement_reason",
]

R1_ME_REPLACEMENT_COLUMNS = [
    "superlocus_id", "candidate_model_id", "candidate_source", "candidate_sublocus_id",
    "candidate_cds_exon_count", "candidate_model_selection_score",
    "candidate_corroborating_other_source_count", "replaced_model_count",
    "replaced_model_ids", "replaced_model_sources", "replaced_model_sublocus_ids",
    "replaced_model_cds_exon_counts", "replaced_model_exon_classes",
    "replaced_model_selection_scores",
    "replaced_model_corroborating_other_source_counts", "conflict_reasons",
    "conflict_overlap_bp", "conflict_candidate_overlap_fractions",
    "conflict_replaced_model_overlap_fractions", "conflict_shared_splice_site_counts",
    "structural_corroboration_weight", "structural_corroboration_alpha",
    "strict_r1_backfill_replacement", "strict_r1_allow_multiexon_blockers",
    "sequence_penalty_missing_start", "sequence_penalty_missing_stop",
    "r1_eligibility_class", "r1_replacement_reason",
]


def write_r1_replacement_outputs(
    stage5: dict[str, object],
    *,
    replacements_tsv: Optional[Path],
    summary_tsv: Optional[Path],
) -> None:
    """Write append-only strict-R1 decision provenance and aggregate counts."""

    if replacements_tsv is not None:
        write_dict_tsv(
            stage5.get("r1_replacement_rows", []),  # type: ignore[arg-type]
            replacements_tsv,
            (
                R1_ME_REPLACEMENT_COLUMNS
                if stage5.get("r1_allow_multiexon_blockers")
                else R1_REPLACEMENT_COLUMNS
            ),
        )
    if summary_tsv is not None:
        write_dict_tsv(
            stage5.get("r1_summary_rows", []),  # type: ignore[arg-type]
            summary_tsv,
            ["metric", "value", "details"],
        )


def write_sublocus_winners_tsv(
    winners: list[CodingModel],
    subloci: list[Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    path: Path,
) -> None:
    sub_by_model_id = {m.output_transcript_id: sub for sub in subloci for m in sub.models}
    rows: list[dict[str, object]] = []
    for model in sorted(winners, key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.output_transcript_id)):
        sub = sub_by_model_id.get(model.output_transcript_id)
        row = model_structural_rows_by_transcript.get(model.output_transcript_id, {})
        rows.append({
            "superlocus_id": sub.superlocus_id if sub else ".",
            "sublocus_id": sub.sublocus_id if sub else ".",
            "sublocus_type": sub.sublocus_type if sub else ".",
            "source_name": model.source_name,
            "low_source_conflict_flag": row.get("low_source_conflict_flag", "."),
            "low_source_conflict_agreement_factor": row.get("low_source_conflict_agreement_factor", "."),
            "low_source_conflict_reason": row.get("low_source_conflict_reason", "."),
            "input_gene_id": model.input_gene_id,
            "input_transcript_id": model.input_transcript_id,
            "output_gene_id": model.output_gene_id,
            "output_transcript_id": model.output_transcript_id,
            "chrom": model.chrom,
            "cds_start": model.cds_start0 + 1,
            "cds_end": model.cds_end0,
            "strand": model.strand,
            "cds_length": model.cds_length,
            "model_selection_score": row.get("model_selection_score", "."),
            "model_selection_score_scaled": row.get("model_selection_score_scaled", "."),
            "model_selection_score_rank": row.get("model_selection_score_rank", "."),
            "model_selection_order_rank": row.get("model_selection_order_rank", "."),
            "exact_cds_exon_chain_support_model_count": row.get("exact_cds_exon_chain_support_model_count", "."),
            "exact_cds_exon_chain_support_source_count": row.get("exact_cds_exon_chain_support_source_count", "."),
            "exact_cds_exon_chain_support_sources": row.get("exact_cds_exon_chain_support_sources", "."),
            "exact_cds_exon_chain_support_model_ids": row.get("exact_cds_exon_chain_support_model_ids", "."),
            "exact_cds_exon_chain_support_fraction_local": row.get("exact_cds_exon_chain_support_fraction_local", "."),
            "exact_structure_support_score": row.get("exact_structure_support_score", "."),
            "model_selection_tie_group": row.get("model_selection_tie_group", "."),
            "model_selection_tie_break_reason": row.get("model_selection_tie_break_reason", "."),
        })
    write_dict_tsv(rows, path, [
        "superlocus_id", "sublocus_id", "sublocus_type", "source_name",
        "low_source_conflict_flag", "low_source_conflict_agreement_factor", "low_source_conflict_reason",
        "input_gene_id", "input_transcript_id", "output_gene_id", "output_transcript_id",
        "chrom", "cds_start", "cds_end", "strand", "cds_length",
        "model_selection_score", "model_selection_score_scaled", "model_selection_score_rank", "model_selection_order_rank",
        "exact_cds_exon_chain_support_model_count", "exact_cds_exon_chain_support_source_count",
        "exact_cds_exon_chain_support_sources", "exact_cds_exon_chain_support_model_ids",
        "exact_cds_exon_chain_support_fraction_local", "exact_structure_support_score",
        "model_selection_tie_group", "model_selection_tie_break_reason",
    ])


def write_monosublocus_outputs(
    stage4: dict[str, object],
    *,
    monosubloci_tsv: Optional[Path],
    monosublocus_members_tsv: Optional[Path],
    monosublocus_edges_tsv: Optional[Path],
    selected_models_tsv: Optional[Path],
    suppressed_models_tsv: Optional[Path],
) -> None:
    monosublocus_columns = [
        "monosublocus_id", "chrom", "start", "end", "strand", "member_count", "selected_count",
        "suppressed_count", "edge_count", "min_cds_exon_overlap_fraction", "not_all_pairwise_overlapping",
        "sources", "subloci", "selected_model_ids", "suppressed_model_ids",
    ]
    member_columns = [
        "monosublocus_id", "superlocus_id", "sublocus_id", "source_name", "input_gene_id",
        "input_transcript_id", "output_gene_id", "output_transcript_id", "chrom", "cds_start", "cds_end",
        "strand", "cds_length", "low_source_conflict_flag", "low_source_conflict_agreement_factor", "low_source_conflict_reason",
        "sequence_integrity_status",
        "sequence_integrity_penalty_score", "sequence_integrity_penalty_weight", "sequence_integrity_penalty_in_selection", "sequence_integrity_flags",
        "structural_coherence_penalty_score", "structural_coherence_raw_score", "structural_coherence_penalty_weight", "structural_coherence_penalty_in_selection",
        "structural_coherence_region_count", "structural_coherence_max_region",
        "model_selection_score", "model_selection_score_scaled", "model_selection_score_rank",
        "model_selection_order_rank",
        "exact_cds_exon_chain_support_model_count", "exact_cds_exon_chain_support_source_count",
        "exact_cds_exon_chain_support_sources", "exact_cds_exon_chain_support_model_ids",
        "exact_cds_exon_chain_support_fraction_local", "exact_structure_support_score",
        "model_completeness_score", "intrinsic_penalty_score",
        "consensus_similarity_score", "structural_support_score", "pairwise_similarity_score", "source_prior_score",
        "monosublocus_member_order_rank", "monosublocus_selected", "monosublocus_selection_order",
        "monosublocus_suppressed_by", "monosublocus_suppression_reason",
        "monosublocus_suppression_overlap_bp", "monosublocus_suppressed_model_overlap_fraction",
        "monosublocus_selected_model_overlap_fraction", "monosublocus_suppression_max_overlap_fraction",
    ]
    edge_columns = [
        "monosublocus_id", "model_a_id", "model_b_id", "source_a", "source_b", "sublocus_a", "sublocus_b",
        "chrom", "strand", "model_a_cds_length", "model_b_cds_length", "cds_exon_overlap_bp",
        "overlap_fraction_of_a", "overlap_fraction_of_b", "max_overlap_fraction", "min_overlap_fraction",
    ]
    if monosubloci_tsv is not None:
        write_dict_tsv(stage4["monosublocus_rows"], monosubloci_tsv, monosublocus_columns)  # type: ignore[arg-type]
    if monosublocus_members_tsv is not None:
        write_dict_tsv(stage4["monosublocus_member_rows"], monosublocus_members_tsv, member_columns)  # type: ignore[arg-type]
    if monosublocus_edges_tsv is not None:
        write_dict_tsv(stage4["monosublocus_edge_rows"], monosublocus_edges_tsv, edge_columns)  # type: ignore[arg-type]
    if selected_models_tsv is not None:
        write_dict_tsv(stage4["selected_rows"], selected_models_tsv, member_columns)  # type: ignore[arg-type]
    if suppressed_models_tsv is not None:
        write_dict_tsv(stage4["suppressed_rows"], suppressed_models_tsv, member_columns)  # type: ignore[arg-type]


def selection_stage_output_columns(
    *,
    include_structural_corroboration: bool = False,
    include_r1: bool = False,
) -> list[str]:
    columns = [
        "superlocus_id", "sublocus_id", "sublocus_type", "source_name", "input_gene_id",
        "input_transcript_id", "output_gene_id", "output_transcript_id", "chrom", "cds_start", "cds_end", "cds_location",
        "strand", "cds_length", "cds_exon_count", "cds_intron_count",
        "low_source_conflict_flag", "low_source_conflict_agreement_factor", "low_source_conflict_reason",
        "model_structural_support_score", "model_structural_support_score_scaled",
        "structural_coherence_penalty_score", "structural_coherence_raw_score", "structural_coherence_penalty_weight", "structural_coherence_penalty_in_selection",
        "structural_coherence_region_count", "structural_coherence_max_region",
        "sequence_integrity_status",
        "sequence_integrity_penalty_score", "sequence_integrity_penalty_weight", "sequence_integrity_penalty_in_selection", "sequence_integrity_flags",
        "protein_evidence_score", "protein_evidence_score_scaled", "protein_evidence_score_weight", "protein_evidence_score_in_selection",
        "transcript_evidence_score", "transcript_evidence_score_scaled", "transcript_evidence_score_weight", "transcript_evidence_score_in_selection",
        "model_selection_score", "model_selection_score_scaled", "model_selection_score_rank",
        "model_selection_order_rank", "exact_cds_exon_chain_support_model_count",
        "exact_cds_exon_chain_support_source_count", "exact_cds_exon_chain_support_sources",
        "exact_cds_exon_chain_support_model_ids", "exact_cds_exon_chain_support_fraction_local",
        "exact_structure_support_score", "model_completeness_score", "intrinsic_penalty_score",
        "consensus_similarity_score", "structural_support_score", "pairwise_similarity_score", "source_prior_score",
        "final_selected", "final_selection_stage", "backfill_candidate_order_rank", "backfill_selected",
        "backfill_selection_order", "backfill_rejection_reason", "backfill_rejected_by", "backfill_rejected_by_source",
        "backfill_conflict_shared_splice_site_count", "backfill_conflict_overlap_bp",
        "backfill_candidate_overlap_fraction", "backfill_selected_model_overlap_fraction",
        "backfill_conflict_max_overlap_fraction",
    ]
    if include_structural_corroboration:
        insertion = columns.index("model_structural_support_score")
        columns[insertion:insertion] = [
            "corroborating_other_source_count", "run_possible_other_source_count",
            "structural_corroboration_score", "structural_corroboration_alpha",
            "structural_corroboration_weight", "structural_corroboration_contribution",
            "structural_corroboration_in_selection",
        ]
    if include_r1:
        columns.extend([
            "r1_replacement_attempted", "r1_replacement_accepted",
            "r1_replacement_reason", "r1_replaced_model_ids",
        ])
    return columns


def final_prune_blank_fields() -> dict[str, object]:
    return {col: "." for col in FINAL_PRUNE_EXTRA_COLUMNS}


def make_initial_final_prune_result(
    stage5: dict[str, object],
    subloci: list[Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Create a generic final-pruning result from the unpruned final selection."""
    sub_by_model_id: dict[str, Sublocus] = {m.output_transcript_id: sub for sub in subloci for m in sub.models}
    final_models: list[CodingModel] = list(stage5.get("final_models", []))  # type: ignore[arg-type]
    final_rows: list[dict[str, object]] = list(stage5.get("final_rows", []))  # type: ignore[arg-type]
    final_rows_by_id = {str(row.get("output_transcript_id")): row for row in final_rows}
    final_info_by_id: dict[str, dict[str, object]] = dict(stage5.get("final_selection_info_by_transcript", {}))  # type: ignore[arg-type]

    retained_rows: list[dict[str, object]] = []
    retained_info_by_id: dict[str, dict[str, object]] = {}
    for model in final_models:
        model_id = model.output_transcript_id
        row = dict(final_rows_by_id.get(model_id, stage_model_base_row(model, sub_by_model_id, model_structural_rows_by_transcript)))
        row.update(final_prune_blank_fields())
        row["final_prune_level"] = "."
        row["final_prune_classes"] = "."
        retained_rows.append(row)
        info = dict(final_info_by_id.get(model_id, row))
        info.update({k: v for k, v in row.items() if k in FINAL_PRUNE_EXTRA_COLUMNS})
        retained_info_by_id[model_id] = info
    return {
        "initial_retained_models": list(final_models),
        "initial_retained_rows": [dict(row) for row in retained_rows],
        "retained_models": final_models,
        "retained_rows": retained_rows,
        "removed_rows": [],
        "retained_selection_info_by_transcript": retained_info_by_id,
        "shadow_pair_rows": [],
        "shadow_candidate_pairs_before_singleton": 0,
        "shadow_candidate_pairs_after_singleton": 0,
        "shadow_candidate_pairs_skipped_after_singleton": 0,
    }


def append_prune_class(row: dict[str, object], prune_class: str) -> None:
    current = str(row.get("final_prune_classes", "."))
    if current in {"", "."}:
        row["final_prune_classes"] = prune_class
    else:
        classes = current.split(",")
        if prune_class not in classes:
            row["final_prune_classes"] = current + "," + prune_class


def build_monosublocus_rescue_context(
    stage4_result: Optional[dict[str, object]],
) -> tuple[dict[str, int], dict[str, str], dict[str, list[dict[str, object]]]]:
    """Return monosublocus holder/member/edge indexes for singleton-pruning rescue."""
    if not stage4_result:
        return {}, {}, {}
    holder_member_count_by_id: dict[str, int] = {}
    for row in stage4_result.get("monosublocus_rows", []):  # type: ignore[union-attr]
        if not isinstance(row, dict):
            continue
        mono_id = str(row.get("monosublocus_id", "."))
        if mono_id == ".":
            continue
        holder_member_count_by_id[mono_id] = int(value_as_float(row.get("member_count", 0), 0.0))

    monosublocus_by_model_id: dict[str, str] = {}
    for row in stage4_result.get("monosublocus_member_rows", []):  # type: ignore[union-attr]
        if not isinstance(row, dict):
            continue
        model_id = str(row.get("output_transcript_id", "."))
        mono_id = str(row.get("monosublocus_id", "."))
        if model_id != "." and mono_id != ".":
            monosublocus_by_model_id[model_id] = mono_id

    edges_by_model_id: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in stage4_result.get("monosublocus_edge_rows", []):  # type: ignore[union-attr]
        if not isinstance(row, dict):
            continue
        model_a = str(row.get("model_a_id", "."))
        model_b = str(row.get("model_b_id", "."))
        if model_a != ".":
            edges_by_model_id[model_a].append(row)
        if model_b != ".":
            edges_by_model_id[model_b].append(row)
    return holder_member_count_by_id, monosublocus_by_model_id, edges_by_model_id


def best_singleton_monosublocus_rescue(
    model: CodingModel,
    *,
    holder_member_count_by_id: dict[str, int],
    monosublocus_by_model_id: dict[str, str],
    edges_by_model_id: dict[str, list[dict[str, object]]],
    min_holder_member_count: int,
    min_overlap_bp: int,
    min_self_overlap_fraction: float,
    min_reciprocal_overlap_fraction: float = 0.80,
) -> tuple[Optional[dict[str, object]], dict[str, object]]:
    """Return the strongest qualifying rescue edge and an interpretable geometry audit."""
    model_id = model.output_transcript_id
    mono_id = monosublocus_by_model_id.get(model_id, ".")
    audit: dict[str, object] = {
        "final_singleton_rescue_monosublocus_id": mono_id,
        "final_singleton_rescue_partner_model_id": ".",
        "final_singleton_rescue_partner_source": ".",
        "final_singleton_rescue_overlap_bp": ".",
        "final_singleton_rescue_self_overlap_fraction": ".",
        "final_singleton_rescue_partner_overlap_fraction": ".",
        "final_singleton_rescue_reciprocal_overlap_fraction": ".",
        "final_singleton_rescue_min_reciprocal_overlap_fraction": min_reciprocal_overlap_fraction,
        "final_singleton_rescue_holder_member_count": ".",
        "final_singleton_overlap_rescue_geometry_reason": "no_monosublocus_holder",
        "final_singleton_overlap_rescue_overlap_bp_geometry_candidate": "false",
        "final_singleton_overlap_rescue_reciprocal_rejected": "false",
        "final_singleton_overlap_rescue_same_source_edge": "false",
    }
    if mono_id == ".":
        return None, audit
    holder_member_count = holder_member_count_by_id.get(mono_id, 1)
    audit["final_singleton_rescue_holder_member_count"] = holder_member_count
    if holder_member_count < min_holder_member_count:
        audit["final_singleton_overlap_rescue_geometry_reason"] = "holder_member_count_below_threshold"
        return None, audit
    candidates: list[dict[str, object]] = []
    diagnostic_edges: list[dict[str, object]] = []
    reciprocal_rejected = False
    same_source_edge_seen = False
    for edge in edges_by_model_id.get(model_id, []):
        edge_mono_id = str(edge.get("monosublocus_id", mono_id))
        if edge_mono_id != mono_id:
            continue
        overlap_bp = int(value_as_float(edge.get("cds_exon_overlap_bp", 0), 0.0))
        if str(edge.get("model_a_id")) == model_id:
            partner_id = str(edge.get("model_b_id", "."))
            partner_source = str(edge.get("source_b", "."))
            self_fraction = value_as_float(edge.get("overlap_fraction_of_a", 0.0), 0.0)
            partner_fraction = value_as_float(edge.get("overlap_fraction_of_b", 0.0), 0.0)
        elif str(edge.get("model_b_id")) == model_id:
            partner_id = str(edge.get("model_a_id", "."))
            partner_source = str(edge.get("source_a", "."))
            self_fraction = value_as_float(edge.get("overlap_fraction_of_b", 0.0), 0.0)
            partner_fraction = value_as_float(edge.get("overlap_fraction_of_a", 0.0), 0.0)
        else:
            continue
        reciprocal_fraction = min(self_fraction, partner_fraction)
        edge_fields = {
            "final_singleton_rescue_enabled": "true",
            "final_singleton_rescued": "true",
            "final_singleton_rescue_reason": "monosublocus_overlap",
            "final_singleton_rescue_monosublocus_id": mono_id,
            "final_singleton_rescue_partner_model_id": partner_id,
            "final_singleton_rescue_partner_source": partner_source,
            "final_singleton_rescue_overlap_bp": overlap_bp,
            "final_singleton_rescue_self_overlap_fraction": self_fraction,
            "final_singleton_rescue_partner_overlap_fraction": partner_fraction,
            "final_singleton_rescue_reciprocal_overlap_fraction": reciprocal_fraction,
            "final_singleton_rescue_min_reciprocal_overlap_fraction": min_reciprocal_overlap_fraction,
            "final_singleton_rescue_holder_member_count": holder_member_count,
        }
        same_source_edge = partner_source == model.source_name
        if same_source_edge:
            same_source_edge_seen = True
        if overlap_bp < min_overlap_bp:
            diagnostic_edges.append({
                **edge_fields,
                "_geometry_failure_stage": 1,
                "_geometry_failure_reason": "overlap_bp_below_threshold",
            })
            continue
        audit["final_singleton_overlap_rescue_overlap_bp_geometry_candidate"] = "true"
        if self_fraction + 1e-15 < min_self_overlap_fraction:
            diagnostic_edges.append({
                **edge_fields,
                "_geometry_failure_stage": 2,
                "_geometry_failure_reason": "self_overlap_below_threshold",
            })
            continue
        if reciprocal_fraction + 1e-15 < min_reciprocal_overlap_fraction:
            reciprocal_rejected = True
            diagnostic_edges.append({
                **edge_fields,
                "_geometry_failure_stage": 3,
                "_geometry_failure_reason": "reciprocal_overlap_below_threshold",
            })
            continue
        if same_source_edge:
            diagnostic_edges.append({
                **edge_fields,
                "_geometry_failure_stage": 4,
                "_geometry_failure_reason": "no_qualifying_different_source_partner",
            })
            continue
        candidates.append(edge_fields)

    def edge_sort_key(row: dict[str, object]) -> tuple[float, int, float, str]:
        return (
            -value_as_float(row.get("final_singleton_rescue_self_overlap_fraction"), 0.0),
            -int(value_as_float(row.get("final_singleton_rescue_overlap_bp"), 0.0)),
            -value_as_float(row.get("final_singleton_rescue_partner_overlap_fraction"), 0.0),
            str(row.get("final_singleton_rescue_partner_model_id", ".")),
        )

    if not candidates:
        if diagnostic_edges:
            diagnostic_edges.sort(key=lambda row: (
                -int(value_as_float(row.get("_geometry_failure_stage"), 0.0)),
                *edge_sort_key(row),
            ))
            selected_diagnostic = diagnostic_edges[0]
            audit.update({
                key: value
                for key, value in selected_diagnostic.items()
                if not key.startswith("_geometry_failure_")
            })
            audit["final_singleton_rescued"] = "false"
            audit["final_singleton_rescue_reason"] = "."
            audit["final_singleton_overlap_rescue_geometry_reason"] = str(
                selected_diagnostic["_geometry_failure_reason"]
            )
        audit["final_singleton_overlap_rescue_reciprocal_rejected"] = (
            "true" if reciprocal_rejected else "false"
        )
        audit["final_singleton_overlap_rescue_same_source_edge"] = (
            "true" if same_source_edge_seen else "false"
        )
        if not diagnostic_edges:
            audit["final_singleton_overlap_rescue_geometry_reason"] = "no_monosublocus_overlap_edge"
        return None, audit
    candidates.sort(key=edge_sort_key)
    selected = candidates[0]
    audit.update(selected)
    audit.update({
        "final_singleton_overlap_rescue_geometry_reason": "passed",
        "final_singleton_overlap_rescue_reciprocal_rejected": "false",
        "final_singleton_overlap_rescue_same_source_edge": (
            "true" if same_source_edge_seen else "false"
        ),
    })
    return selected, audit


def singleton_protein_rescue_fields(
    model_id: str,
    protein_evidence_rows_by_transcript: dict[str, dict[str, object]],
    *,
    rescue_protein_supported_singletons: bool,
    min_score_scaled: float,
    min_exact_chain_source_count: int,
) -> dict[str, object]:
    """Return protein-evidence singleton rescue fields for an otherwise pruned model."""
    protein_row = protein_evidence_rows_by_transcript.get(model_id, {})
    score_scaled = value_as_float(protein_row.get("protein_evidence_score_scaled", 0.0), 0.0)
    exact_source_count = int(value_as_float(protein_row.get("protein_exact_cds_chain_source_count", 0), 0.0))
    supporting_source_count = int(value_as_float(protein_row.get("protein_supporting_source_count", 0), 0.0))
    status = str(protein_row.get("protein_evidence_status", "protein_evidence_not_configured")) if protein_row else "protein_evidence_not_configured"
    flags = str(protein_row.get("protein_evidence_flags", ".")) if protein_row else "."

    reasons: list[str] = []
    if rescue_protein_supported_singletons:
        if score_scaled + 1e-15 >= min_score_scaled:
            reasons.append("protein_evidence_score_scaled")
        if min_exact_chain_source_count > 0 and exact_source_count >= min_exact_chain_source_count:
            reasons.append("protein_exact_cds_chain_support")

    criteria_met = bool(reasons)
    criteria_reason = "+".join(reasons) if reasons else "."
    return {
        "final_protein_singleton_rescue_enabled": "true" if rescue_protein_supported_singletons else "false",
        # Actual rescue is set only if singleton pruning would remove the model and protein evidence spares it.
        "final_protein_singleton_rescued": "false",
        "final_protein_singleton_rescue_reason": criteria_reason,
        "final_protein_singleton_rescue_criteria_met": "true" if criteria_met else "false",
        "final_protein_singleton_rescue_criteria_reason": criteria_reason,
        "final_protein_singleton_rescue_applied": "false",
        "final_protein_singleton_rescue_applied_reason": ".",
        "final_protein_singleton_rescue_min_score_scaled": min_score_scaled,
        "final_protein_singleton_rescue_min_exact_chain_source_count": min_exact_chain_source_count,
        "final_protein_singleton_rescue_score_scaled": score_scaled if protein_row else ".",
        "final_protein_singleton_rescue_exact_chain_source_count": exact_source_count if protein_row else ".",
        "final_protein_singleton_rescue_supporting_source_count": supporting_source_count if protein_row else ".",
        "final_protein_singleton_rescue_status": status,
        "final_protein_singleton_rescue_flags": flags,
    }


def _transcript_rescue_support_type(stranded_count: int, unstranded_count: int) -> str:
    if stranded_count > 0 and unstranded_count > 0:
        return "mixed"
    if stranded_count > 0:
        return "stranded"
    if unstranded_count > 0:
        return "unstranded_only"
    return "."


def singleton_transcript_rescue_fields(
    model_id: str,
    transcript_evidence_rows_by_transcript: dict[str, dict[str, object]],
    *,
    rescue_transcript_supported_singletons: bool,
    min_score_scaled: float,
    min_supporting_source_count: int,
    min_scoring_eligible_count: int,
    min_cdna_context_fraction: float = 0.30,
    require_stranded_support: bool,
    block_on_cdna_context_concern: bool = True,
) -> dict[str, object]:
    """Return transcript-evidence singleton rescue fields for an otherwise pruned model.

    Transcript rescue is deliberately stricter than protein exact-chain rescue: all
    configured transcript criteria must be met, and only scoring-eligible transcript
    associations contribute to the transcript evidence score and counts.
    """
    transcript_row = transcript_evidence_rows_by_transcript.get(model_id, {})
    score_scaled = value_as_float(transcript_row.get("transcript_evidence_score_scaled", 0.0), 0.0)
    supporting_source_count = int(value_as_float(transcript_row.get("transcript_supporting_source_count", 0), 0.0))
    scoring_eligible_count = int(value_as_float(transcript_row.get("transcript_scoring_eligible_count", 0), 0.0))
    stranded_scoring_count = int(value_as_float(transcript_row.get("transcript_stranded_scoring_eligible_count", 0), 0.0))
    unstranded_scoring_count = int(value_as_float(transcript_row.get("transcript_unstranded_scoring_eligible_count", 0), 0.0))
    support_type = _transcript_rescue_support_type(stranded_scoring_count, unstranded_scoring_count)
    status = str(transcript_row.get("transcript_evidence_status", "transcript_evidence_not_configured")) if transcript_row else "transcript_evidence_not_configured"
    flags = str(transcript_row.get("transcript_evidence_flags", ".")) if transcript_row else "."
    context_concern = str(transcript_row.get("transcript_cdna_context_pruning_concern", "false")) == "true"
    context_class = str(transcript_row.get("transcript_cdna_context_pruning_concern_class", "not_evaluable"))
    context_fraction = transcript_row.get("transcript_cdna_context_cds_to_min_cdna_fraction", ".") if transcript_row else "."
    context_evaluable = bool(
        transcript_row
        and str(transcript_row.get("transcript_cdna_context_evaluable", "false")) == "true"
    )
    context_fraction_value = value_as_float(context_fraction, -1.0)
    context_fraction_passed = bool(
        context_evaluable
        and context_fraction_value + 1e-15 >= min_cdna_context_fraction
    )
    context_fraction_failed = bool(
        context_evaluable
        and not context_fraction_passed
    )
    context_fraction_requirement_met = bool(
        context_fraction_passed or not context_evaluable
    )
    context_eligible = bool(
        context_fraction_requirement_met
        and not (block_on_cdna_context_concern and context_concern)
    )
    if not context_evaluable:
        context_reason = "context_fraction_not_evaluable_neutral"
        context_fraction_state = "not_evaluable_neutral"
    elif context_concern and not block_on_cdna_context_concern:
        context_reason = f"{context_class}_not_blocking"
        context_fraction_state = (
            "passed" if context_fraction_passed else "failed_below_minimum"
        )
    elif context_concern and block_on_cdna_context_concern:
        context_reason = f"blocked_by_{context_class}"
        context_fraction_state = (
            "passed" if context_fraction_passed else "failed_below_minimum"
        )
    elif context_fraction_failed:
        context_reason = "context_fraction_failed_below_minimum"
        context_fraction_state = "failed_below_minimum"
    elif context_eligible:
        context_reason = "context_fraction_passed"
        context_fraction_state = "passed"
    else:
        context_reason = f"blocked_by_{context_class}"
        context_fraction_state = (
            "passed" if context_fraction_passed else "failed_below_minimum"
        )

    criteria_parts: list[str] = []
    criteria_met = False
    if rescue_transcript_supported_singletons and transcript_row:
        score_ok = score_scaled + 1e-15 >= min_score_scaled
        source_ok = supporting_source_count >= min_supporting_source_count
        eligible_ok = scoring_eligible_count >= min_scoring_eligible_count
        stranded_ok = (not require_stranded_support) or stranded_scoring_count > 0
        if score_ok:
            criteria_parts.append("transcript_evidence_score_scaled")
        if source_ok:
            criteria_parts.append("transcript_supporting_source_count")
        if eligible_ok:
            criteria_parts.append("transcript_scoring_eligible_count")
        if require_stranded_support and stranded_ok:
            criteria_parts.append("transcript_stranded_scoring_eligible")
        if context_fraction_passed:
            criteria_parts.append("transcript_cdna_context_cds_to_min_cdna_fraction")
        elif not context_evaluable:
            criteria_parts.append("transcript_cdna_context_not_evaluable_neutral")
        criteria_met = bool(score_ok and source_ok and eligible_ok and stranded_ok and context_eligible)

    criteria_reason = "+".join(criteria_parts) if criteria_met else "."
    return {
        "final_transcript_singleton_rescue_enabled": "true" if rescue_transcript_supported_singletons else "false",
        # Actual rescue is set only if singleton pruning would remove the model and transcript evidence spares it.
        "final_transcript_singleton_rescued": "false",
        "final_transcript_singleton_rescue_reason": criteria_reason,
        "final_transcript_singleton_rescue_criteria_met": "true" if criteria_met else "false",
        "final_transcript_singleton_rescue_criteria_reason": criteria_reason,
        "final_transcript_singleton_rescue_applied": "false",
        "final_transcript_singleton_rescue_applied_reason": ".",
        "final_transcript_singleton_rescue_min_score_scaled": min_score_scaled,
        "final_transcript_singleton_rescue_min_supporting_source_count": min_supporting_source_count,
        "final_transcript_singleton_rescue_min_scoring_eligible_count": min_scoring_eligible_count,
        "final_transcript_singleton_rescue_min_cdna_context_fraction": min_cdna_context_fraction,
        "final_transcript_singleton_rescue_require_stranded_support": "true" if require_stranded_support else "false",
        "final_transcript_singleton_rescue_block_on_cdna_context_concern": "true" if block_on_cdna_context_concern else "false",
        "final_transcript_singleton_rescue_cdna_context_eligible": "true" if context_eligible else "false",
        "final_transcript_singleton_rescue_cdna_context_reason": context_reason,
        "final_transcript_singleton_rescue_cdna_context_fraction_state": context_fraction_state,
        "final_transcript_singleton_rescue_cdna_context_fraction": context_fraction,
        "final_transcript_singleton_rescue_cdna_context_concern_class": context_class,
        "final_transcript_singleton_rescue_score_scaled": score_scaled if transcript_row else ".",
        "final_transcript_singleton_rescue_supporting_source_count": supporting_source_count if transcript_row else ".",
        "final_transcript_singleton_rescue_scoring_eligible_count": scoring_eligible_count if transcript_row else ".",
        "final_transcript_singleton_rescue_stranded_scoring_eligible_count": stranded_scoring_count if transcript_row else ".",
        "final_transcript_singleton_rescue_unstranded_scoring_eligible_count": unstranded_scoring_count if transcript_row else ".",
        "final_transcript_singleton_rescue_support_type": support_type if transcript_row else ".",
        "final_transcript_singleton_rescue_status": status,
        "final_transcript_singleton_rescue_flags": flags,
    }

SEQUENCE_INTEGRITY_STANDARD_PRUNING_FLAGS = frozenset(
    {
        "internal_stop_codon",
        "missing_start_codon",
        "missing_stop_codon",
        "all_splices_noncanonical",
    }
)


def sequence_integrity_pruning_fields(row: dict[str, object]) -> dict[str, object]:
    """Classify whitelisted sequence-integrity concerns for final pruning.

    A standard concern is any internal stop, missing start, missing stop, or
    all-splices-noncanonical flag.  The severe class is restricted to an
    internal stop, both terminal codons missing together, or an entirely
    noncanonical splice chain.
    """
    raw_flags = str(row.get("sequence_integrity_flags", "."))
    flags = {flag for flag in raw_flags.split(",") if flag and flag != "."}
    reasons = sorted(flags & SEQUENCE_INTEGRITY_STANDARD_PRUNING_FLAGS)
    severe_reasons: list[str] = []
    if "internal_stop_codon" in flags:
        severe_reasons.append("internal_stop_codon")
    if {"missing_start_codon", "missing_stop_codon"}.issubset(flags):
        severe_reasons.append("missing_start_and_stop_codon")
    if "all_splices_noncanonical" in flags:
        severe_reasons.append("all_splices_noncanonical")
    concern = bool(reasons)
    severe = bool(severe_reasons)
    if severe:
        concern_class = "severe"
    elif concern:
        concern_class = "standard"
    else:
        concern_class = "none"
    return {
        "sequence_integrity_pruning_concern": "true" if concern else "false",
        "sequence_integrity_pruning_concern_class": concern_class,
        "sequence_integrity_pruning_reasons": ",".join(reasons) if reasons else ".",
        "sequence_integrity_severe_pruning_concern": "true" if severe else "false",
        "sequence_integrity_severe_pruning_reasons": (
            ",".join(severe_reasons) if severe_reasons else "."
        ),
    }


def internal_stop_blocks_evidence_rescue(
    sequence_fields: dict[str, object],
    *,
    allow_internal_stop_evidence_rescue: bool,
) -> bool:
    """Return whether the resolved policy blocks evidence override for this model."""
    reasons = {
        reason
        for field in (
            "sequence_integrity_pruning_reasons",
            "sequence_integrity_severe_pruning_reasons",
        )
        for reason in str(sequence_fields.get(field, ".")).split(",")
        if reason and reason != "."
    }
    return (
        "internal_stop_codon" in reasons
        and not allow_internal_stop_evidence_rescue
    )


def internal_stop_evidence_audit_fields(
    sequence_fields: dict[str, object],
    *,
    allow_internal_stop_evidence_rescue: bool,
    protein_criteria_met: bool,
    transcript_criteria_met: bool,
    evidence_override_candidate: bool,
) -> dict[str, object]:
    """Describe evidence criteria and the internal-stop override policy separately."""
    has_internal_stop = "internal_stop_codon" in {
        reason
        for field in (
            "sequence_integrity_pruning_reasons",
            "sequence_integrity_severe_pruning_reasons",
        )
        for reason in str(sequence_fields.get(field, ".")).split(",")
        if reason and reason != "."
    }
    routes = []
    if protein_criteria_met:
        routes.append("protein_evidence")
    if transcript_criteria_met:
        routes.append("transcript_evidence")
    criteria_met = bool(routes)
    blocked = bool(
        evidence_override_candidate
        and has_internal_stop
        and criteria_met
        and not allow_internal_stop_evidence_rescue
    )
    return {
        "final_model_has_internal_stop_codon": "true" if has_internal_stop else "false",
        "final_allow_internal_stop_evidence_rescue": (
            "true" if allow_internal_stop_evidence_rescue else "false"
        ),
        "final_internal_stop_evidence_criteria_met": "true" if criteria_met else "false",
        "final_internal_stop_evidence_criteria_routes": ",".join(routes) if routes else ".",
        "final_internal_stop_evidence_rescue_blocked": "true" if blocked else "false",
        "final_internal_stop_evidence_rescue_block_reason": (
            "internal_stop_codon_override_disabled" if blocked else "."
        ),
        "final_internal_stop_evidence_rescue_applied": "false",
        "final_internal_stop_evidence_rescue_applied_route": ".",
    }


def merge_internal_stop_evidence_audit_fields(
    existing: dict[str, object],
    current_stage: dict[str, object],
) -> dict[str, object]:
    """Merge internal-stop evidence provenance without erasing earlier stages."""

    def is_true(row: dict[str, object], field: str) -> bool:
        return str(row.get(field, "false")).lower() == "true"

    def merged_values(field: str) -> str:
        values: list[str] = []
        for row in (existing, current_stage):
            for value in str(row.get(field, ".")).split(","):
                if value and value != "." and value not in values:
                    values.append(value)
        return ",".join(values) if values else "."

    has_internal_stop = any(
        is_true(row, "final_model_has_internal_stop_codon")
        for row in (existing, current_stage)
    )
    criteria_met = any(
        is_true(row, "final_internal_stop_evidence_criteria_met")
        for row in (existing, current_stage)
    )
    rescue_blocked = any(
        is_true(row, "final_internal_stop_evidence_rescue_blocked")
        for row in (existing, current_stage)
    )
    rescue_applied = any(
        is_true(row, "final_internal_stop_evidence_rescue_applied")
        for row in (existing, current_stage)
    )
    configured_allow = current_stage.get(
        "final_allow_internal_stop_evidence_rescue",
        existing.get("final_allow_internal_stop_evidence_rescue", "false"),
    )
    return {
        "final_model_has_internal_stop_codon": str(has_internal_stop).lower(),
        "final_allow_internal_stop_evidence_rescue": configured_allow,
        "final_internal_stop_evidence_criteria_met": str(criteria_met).lower(),
        "final_internal_stop_evidence_criteria_routes": (
            merged_values("final_internal_stop_evidence_criteria_routes")
            if criteria_met else "."
        ),
        "final_internal_stop_evidence_rescue_blocked": str(rescue_blocked).lower(),
        "final_internal_stop_evidence_rescue_block_reason": (
            merged_values("final_internal_stop_evidence_rescue_block_reason")
            if rescue_blocked else "."
        ),
        "final_internal_stop_evidence_rescue_applied": str(rescue_applied).lower(),
        "final_internal_stop_evidence_rescue_applied_route": (
            merged_values("final_internal_stop_evidence_rescue_applied_route")
            if rescue_applied else "."
        ),
    }


def count_retained_internal_stop_evidence_overrides(
    retained_rows: list[dict[str, object]],
) -> int:
    """Count retained internal-stop rows whose final provenance records an override."""
    return sum(
        str(row.get("final_model_has_internal_stop_codon", "false")).lower() == "true"
        and str(
            row.get("final_internal_stop_evidence_rescue_applied", "false")
        ).lower() == "true"
        for row in retained_rows
    )


def final_prune_singleton_selected_models(
    prune_result: dict[str, object],
    subloci: list[Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    *,
    prune_level: str,
    prune_sources: set[str],
    pruning_level_name: str,
    rescue_monosublocus_overlap: bool = False,
    stage4_result: Optional[dict[str, object]] = None,
    rescue_min_holder_member_count: int = 2,
    rescue_min_overlap_bp: int = 150,
    rescue_min_self_overlap_fraction: float = 0.80,
    rescue_min_reciprocal_overlap_fraction: float = 0.80,
    monosublocus_overlap_rescue_evidence_policy: str = "adequacy-conditioned",
    protein_adequacy_passed: bool = False,
    transcript_adequacy_passed: bool = False,
    rescue_protein_supported_singletons: bool = False,
    protein_evidence_rows_by_transcript: Optional[dict[str, dict[str, object]]] = None,
    rescue_min_protein_evidence_score_scaled: float = 0.75,
    rescue_min_protein_exact_chain_source_count: int = 1,
    rescue_transcript_supported_singletons: bool = False,
    transcript_evidence_rows_by_transcript: Optional[dict[str, dict[str, object]]] = None,
    rescue_min_transcript_evidence_score_scaled: float = 0.75,
    rescue_min_transcript_supporting_source_count: int = 1,
    rescue_min_transcript_scoring_eligible_count: int = 1,
    rescue_min_transcript_cdna_context_fraction: float = 0.30,
    transcript_rescue_require_stranded_support: bool = False,
    transcript_rescue_block_on_cdna_context_concern: bool = True,
    allow_internal_stop_evidence_rescue: bool = False,
    progress: Optional[ProgressLogger] = None,
) -> dict[str, object]:
    """Apply singleton-sublocus pruning to a current final-pruning result."""
    sub_by_model_id: dict[str, Sublocus] = {m.output_transcript_id: sub for sub in subloci for m in sub.models}
    final_models: list[CodingModel] = list(prune_result.get("retained_models", []))  # type: ignore[arg-type]
    retained_rows_in: list[dict[str, object]] = list(prune_result.get("retained_rows", []))  # type: ignore[arg-type]
    rows_by_id = {str(row.get("output_transcript_id")): row for row in retained_rows_in}
    info_by_id: dict[str, dict[str, object]] = dict(prune_result.get("retained_selection_info_by_transcript", {}))  # type: ignore[arg-type]
    protein_evidence_rows_by_transcript = protein_evidence_rows_by_transcript or {}
    transcript_evidence_rows_by_transcript = transcript_evidence_rows_by_transcript or {}

    retained_models: list[CodingModel] = []
    retained_rows: list[dict[str, object]] = []
    retained_info_by_id: dict[str, dict[str, object]] = {}
    removed_rows: list[dict[str, object]] = list(prune_result.get("removed_rows", []))  # type: ignore[arg-type]
    source_scope = ",".join(sorted(prune_sources)) if prune_sources else "ALL"
    (
        rescue_holder_member_count_by_id,
        rescue_monosublocus_by_model_id,
        rescue_edges_by_model_id,
    ) = build_monosublocus_rescue_context(stage4_result if rescue_monosublocus_overlap else None)
    rescue_count = 0
    protein_rescue_count = 0
    transcript_rescue_count = 0
    sequence_candidate_count = 0
    sequence_pruned_count = 0
    sequence_rescued_count = 0
    overlap_geometry_count = 0
    overlap_bp_geometry_candidate_count = 0
    overlap_reciprocal_rejected_count = 0
    overlap_models_with_same_source_edge_count = 0
    overlap_policy_denied_count = 0
    overlap_rescued_by_strong_local_count = 0
    overlap_rescued_without_adequate_evidence_count = 0
    internal_stop_evidence_rescue_blocked_count = 0
    internal_stop_evidence_rescue_applied_count = 0

    if progress:
        progress.log(
            f"Final selected pruning: singleton level={prune_level}; source_scope={source_scope}; "
            f"models={len(final_models):,}; rescue_monosublocus_overlap={str(rescue_monosublocus_overlap).lower()}; "
            f"rescue_protein_supported_singletons={str(rescue_protein_supported_singletons).lower()}; "
            f"rescue_transcript_supported_singletons={str(rescue_transcript_supported_singletons).lower()}; "
            f"monosublocus_overlap_rescue_evidence_policy={monosublocus_overlap_rescue_evidence_policy}; "
            f"rescue_min_reciprocal_overlap_fraction={rescue_min_reciprocal_overlap_fraction:.6f}; "
            f"protein_adequacy_passed={str(protein_adequacy_passed).lower()}; "
            f"transcript_adequacy_passed={str(transcript_adequacy_passed).lower()}"
            f"; allow_internal_stop_evidence_rescue={str(allow_internal_stop_evidence_rescue).lower()}"
        )

    for model in final_models:
        model_id = model.output_transcript_id
        sub = sub_by_model_id.get(model_id)
        row = dict(rows_by_id.get(model_id, stage_model_base_row(model, sub_by_model_id, model_structural_rows_by_transcript)))
        row.update({k: row.get(k, ".") for k in FINAL_PRUNE_EXTRA_COLUMNS})
        info = dict(info_by_id.get(model_id, row))
        source_match = (not prune_sources) or (model.source_name in prune_sources)
        is_singleton = bool(sub and len(sub.models) == 1)
        structural_row = model_structural_rows_by_transcript.get(model_id, {})
        sequence_fields = sequence_integrity_pruning_fields({**structural_row, **row})
        row.update(sequence_fields)
        intrinsic_components = str(structural_row.get("intrinsic_penalty_components", "."))
        intrinsic_penalty = value_as_float(structural_row.get("intrinsic_penalty_score", 0.0), 0.0)
        intrinsic_candidate = intrinsic_components != "." and intrinsic_penalty > 0
        sequence_concern = sequence_fields["sequence_integrity_pruning_concern"] == "true"
        sequence_candidate = False
        reason = "."
        should_remove = False
        if is_singleton and source_match and prune_level != "none":
            if prune_level == "all":
                should_remove = True
                reason = "singleton_sublocus_all"
            elif prune_level == "intrinsic" and intrinsic_candidate:
                should_remove = True
                reason = "singleton_sublocus_intrinsic_penalty"
            elif prune_level == "sequence_integrity" and sequence_concern:
                should_remove = True
                sequence_candidate = True
                reason = "singleton_sublocus_sequence_integrity"
            elif prune_level == "intrinsic_or_sequence_integrity":
                if intrinsic_candidate and sequence_concern:
                    should_remove = True
                    sequence_candidate = True
                    reason = "singleton_sublocus_intrinsic_and_sequence_integrity"
                elif intrinsic_candidate:
                    should_remove = True
                    reason = "singleton_sublocus_intrinsic_penalty"
                elif sequence_concern:
                    should_remove = True
                    sequence_candidate = True
                    reason = "singleton_sublocus_sequence_integrity"
        if sequence_candidate:
            sequence_candidate_count += 1

        rescue_fields: dict[str, object] = {
            "final_singleton_rescue_enabled": "true" if rescue_monosublocus_overlap else "false",
            "final_singleton_rescued": "false",
            "final_singleton_rescue_reason": ".",
            "final_singleton_rescue_monosublocus_id": ".",
            "final_singleton_rescue_partner_model_id": ".",
            "final_singleton_rescue_partner_source": ".",
            "final_singleton_rescue_overlap_bp": ".",
            "final_singleton_rescue_self_overlap_fraction": ".",
            "final_singleton_rescue_partner_overlap_fraction": ".",
            "final_singleton_rescue_reciprocal_overlap_fraction": ".",
            "final_singleton_rescue_min_reciprocal_overlap_fraction": rescue_min_reciprocal_overlap_fraction,
            "final_singleton_rescue_holder_member_count": ".",
            "final_singleton_overlap_rescue_evidence_policy": monosublocus_overlap_rescue_evidence_policy,
            "final_singleton_overlap_rescue_geometry_passed": "false",
            "final_singleton_overlap_rescue_geometry_reason": ".",
            "final_singleton_overlap_rescue_overlap_bp_geometry_candidate": "false",
            "final_singleton_overlap_rescue_reciprocal_rejected": "false",
            "final_singleton_overlap_rescue_same_source_edge": "false",
            "final_singleton_overlap_rescue_any_adequacy_passed": (
                "true" if protein_adequacy_passed or transcript_adequacy_passed else "false"
            ),
            "final_singleton_overlap_rescue_protein_adequacy_passed": "true" if protein_adequacy_passed else "false",
            "final_singleton_overlap_rescue_transcript_adequacy_passed": "true" if transcript_adequacy_passed else "false",
            "final_singleton_overlap_rescue_strong_local_protein": "false",
            "final_singleton_overlap_rescue_strong_local_transcript": "false",
            "final_singleton_overlap_rescue_decision": ".",
            "final_singleton_overlap_rescue_denial_reason": ".",
        }
        protein_rescue_fields = singleton_protein_rescue_fields(
            model_id,
            protein_evidence_rows_by_transcript,
            rescue_protein_supported_singletons=rescue_protein_supported_singletons,
            min_score_scaled=rescue_min_protein_evidence_score_scaled,
            min_exact_chain_source_count=rescue_min_protein_exact_chain_source_count,
        )
        transcript_rescue_fields = singleton_transcript_rescue_fields(
            model_id,
            transcript_evidence_rows_by_transcript,
            rescue_transcript_supported_singletons=rescue_transcript_supported_singletons,
            min_score_scaled=rescue_min_transcript_evidence_score_scaled,
            min_supporting_source_count=rescue_min_transcript_supporting_source_count,
            min_scoring_eligible_count=rescue_min_transcript_scoring_eligible_count,
            min_cdna_context_fraction=rescue_min_transcript_cdna_context_fraction,
            require_stranded_support=transcript_rescue_require_stranded_support,
            block_on_cdna_context_concern=transcript_rescue_block_on_cdna_context_concern,
        )
        protein_criteria_met = str(
            protein_rescue_fields.get("final_protein_singleton_rescue_criteria_met", "false")
        ) == "true"
        transcript_criteria_met = str(
            transcript_rescue_fields.get("final_transcript_singleton_rescue_criteria_met", "false")
        ) == "true"
        internal_stop_policy_blocks = internal_stop_blocks_evidence_rescue(
            sequence_fields,
            allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
        )
        internal_stop_audit = internal_stop_evidence_audit_fields(
            sequence_fields,
            allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
            protein_criteria_met=protein_criteria_met,
            transcript_criteria_met=transcript_criteria_met,
            evidence_override_candidate=should_remove,
        )
        if internal_stop_audit["final_internal_stop_evidence_rescue_blocked"] == "true":
            internal_stop_evidence_rescue_blocked_count += 1

        if should_remove and rescue_monosublocus_overlap:
            rescue_hit, geometry_audit = best_singleton_monosublocus_rescue(
                model,
                holder_member_count_by_id=rescue_holder_member_count_by_id,
                monosublocus_by_model_id=rescue_monosublocus_by_model_id,
                edges_by_model_id=rescue_edges_by_model_id,
                min_holder_member_count=rescue_min_holder_member_count,
                min_overlap_bp=rescue_min_overlap_bp,
                min_self_overlap_fraction=rescue_min_self_overlap_fraction,
                min_reciprocal_overlap_fraction=rescue_min_reciprocal_overlap_fraction,
            )
            rescue_fields.update(geometry_audit)
            if str(geometry_audit.get("final_singleton_overlap_rescue_overlap_bp_geometry_candidate")) == "true":
                overlap_bp_geometry_candidate_count += 1
            if str(geometry_audit.get("final_singleton_overlap_rescue_reciprocal_rejected")) == "true":
                overlap_reciprocal_rejected_count += 1
            if str(geometry_audit.get("final_singleton_overlap_rescue_same_source_edge")) == "true":
                overlap_models_with_same_source_edge_count += 1
            if rescue_hit is not None:
                overlap_geometry_count += 1
                raw_strong_local_protein = (
                    rescue_protein_supported_singletons
                    and protein_criteria_met
                )
                raw_strong_local_transcript = (
                    rescue_transcript_supported_singletons
                    and transcript_criteria_met
                )
                strong_local_protein = raw_strong_local_protein and not internal_stop_policy_blocks
                strong_local_transcript = raw_strong_local_transcript and not internal_stop_policy_blocks
                any_adequacy_passed = protein_adequacy_passed or transcript_adequacy_passed
                rescue_fields.update(rescue_hit)
                rescue_fields.update({
                    "final_singleton_overlap_rescue_geometry_passed": "true",
                    "final_singleton_overlap_rescue_strong_local_protein": "true" if strong_local_protein else "false",
                    "final_singleton_overlap_rescue_strong_local_transcript": "true" if strong_local_transcript else "false",
                })
                rescue_allowed = monosublocus_overlap_rescue_evidence_policy == "unconditional"
                if monosublocus_overlap_rescue_evidence_policy == "adequacy-conditioned":
                    rescue_allowed = strong_local_protein or strong_local_transcript or not any_adequacy_passed
                if rescue_allowed:
                    should_remove = False
                    reason = "singleton_sublocus_rescued_by_monosublocus_overlap"
                    rescue_fields["final_singleton_overlap_rescue_decision"] = (
                        "rescued_strong_local_evidence"
                        if strong_local_protein or strong_local_transcript
                        else (
                            "rescued_no_adequate_evidence"
                            if monosublocus_overlap_rescue_evidence_policy == "adequacy-conditioned"
                            else "rescued_unconditional"
                        )
                    )
                    rescue_count += 1
                    if strong_local_protein or strong_local_transcript:
                        overlap_rescued_by_strong_local_count += 1
                        if internal_stop_policy_blocks is False and str(
                            internal_stop_audit["final_model_has_internal_stop_codon"]
                        ) == "true" and allow_internal_stop_evidence_rescue:
                            applied_routes = []
                            if strong_local_protein:
                                applied_routes.append("monosublocus_overlap_protein_evidence")
                            if strong_local_transcript:
                                applied_routes.append("monosublocus_overlap_transcript_evidence")
                            internal_stop_audit.update({
                                "final_internal_stop_evidence_rescue_applied": "true",
                                "final_internal_stop_evidence_rescue_applied_route": ",".join(applied_routes),
                            })
                            internal_stop_evidence_rescue_applied_count += 1
                    elif monosublocus_overlap_rescue_evidence_policy == "adequacy-conditioned":
                        overlap_rescued_without_adequate_evidence_count += 1
                else:
                    rescue_fields.update({
                        "final_singleton_rescued": "false",
                        "final_singleton_rescue_reason": ".",
                        "final_singleton_overlap_rescue_decision": "denied",
                        "final_singleton_overlap_rescue_denial_reason": "adequate_evidence_without_strong_local_support",
                    })
                    overlap_policy_denied_count += 1
            else:
                rescue_fields["final_singleton_overlap_rescue_decision"] = "geometry_failed"

        if should_remove and rescue_protein_supported_singletons and protein_criteria_met and not internal_stop_policy_blocks:
            should_remove = False
            reason = "singleton_sublocus_rescued_by_protein_evidence"
            protein_rescue_reason = str(protein_rescue_fields.get("final_protein_singleton_rescue_criteria_reason", "."))
            rescue_fields.update({
                "final_singleton_rescue_enabled": rescue_fields.get("final_singleton_rescue_enabled", "false"),
                "final_singleton_rescued": "true",
                "final_singleton_rescue_reason": "protein_evidence",
            })
            protein_rescue_fields.update({
                "final_protein_singleton_rescued": "true",
                "final_protein_singleton_rescue_applied": "true",
                "final_protein_singleton_rescue_applied_reason": protein_rescue_reason,
            })
            protein_rescue_count += 1
            if str(internal_stop_audit["final_model_has_internal_stop_codon"]) == "true":
                internal_stop_audit.update({
                    "final_internal_stop_evidence_rescue_applied": "true",
                    "final_internal_stop_evidence_rescue_applied_route": "protein_evidence",
                })
                internal_stop_evidence_rescue_applied_count += 1

        if should_remove and rescue_transcript_supported_singletons and transcript_criteria_met and not internal_stop_policy_blocks:
            should_remove = False
            reason = "singleton_sublocus_rescued_by_transcript_evidence"
            transcript_rescue_reason = str(transcript_rescue_fields.get("final_transcript_singleton_rescue_criteria_reason", "."))
            rescue_fields.update({
                "final_singleton_rescue_enabled": rescue_fields.get("final_singleton_rescue_enabled", "false"),
                "final_singleton_rescued": "true",
                "final_singleton_rescue_reason": "transcript_evidence",
            })
            transcript_rescue_fields.update({
                "final_transcript_singleton_rescued": "true",
                "final_transcript_singleton_rescue_applied": "true",
                "final_transcript_singleton_rescue_applied_reason": transcript_rescue_reason,
            })
            transcript_rescue_count += 1
            if str(internal_stop_audit["final_model_has_internal_stop_codon"]) == "true":
                internal_stop_audit.update({
                    "final_internal_stop_evidence_rescue_applied": "true",
                    "final_internal_stop_evidence_rescue_applied_route": "transcript_evidence",
                })
                internal_stop_evidence_rescue_applied_count += 1

        if sequence_candidate:
            if should_remove:
                sequence_pruned_count += 1
            else:
                sequence_rescued_count += 1

        prune_fields = {
            "final_prune_level": pruning_level_name,
            "final_singleton_prune_level": prune_level,
            "final_singleton_prune_source_scope": source_scope,
            "final_singleton_pruned": "true" if should_remove else "false",
            "final_singleton_prune_reason": reason,
            "final_singleton_prune_sublocus_model_count": len(sub.models) if sub else ".",
            "final_singleton_prune_sublocus_source_count": len({m.source_name for m in sub.models}) if sub else ".",
            **sequence_fields,
            "final_sequence_integrity_singleton_candidate": "true" if sequence_candidate else "false",
            "final_sequence_integrity_singleton_rescued": "true" if sequence_candidate and not should_remove else "false",
            "final_sequence_integrity_singleton_pruned": "true" if sequence_candidate and should_remove else "false",
            "final_sequence_integrity_singleton_prune_reason": (
                str(sequence_fields.get("sequence_integrity_pruning_reasons", "."))
                if sequence_candidate else "."
            ),
            **rescue_fields,
            **protein_rescue_fields,
            **transcript_rescue_fields,
            **internal_stop_audit,
        }
        out_row = {**row, **prune_fields}
        if should_remove:
            append_prune_class(out_row, "singleton_sublocus")
            if sequence_candidate:
                append_prune_class(out_row, "sequence_integrity_singleton")
            removed_rows.append(out_row)
        else:
            retained_models.append(model)
            retained_rows.append(out_row)
            retained_info_by_id[model_id] = {**info, **prune_fields, "final_prune_classes": out_row.get("final_prune_classes", ".")}

    if progress:
        progress.log(
            f"Final selected singleton pruning complete: retained={len(retained_models):,}; "
            f"newly_removed={len(removed_rows) - len(prune_result.get('removed_rows', [])):,}; "
            f"rescued_by_monosublocus_overlap={rescue_count:,}; "
            f"rescued_by_protein_evidence={protein_rescue_count:,}; "
            f"rescued_by_transcript_evidence={transcript_rescue_count:,}; "
            f"sequence_candidates={sequence_candidate_count:,}; "
            f"sequence_pruned={sequence_pruned_count:,}; sequence_rescued={sequence_rescued_count:,}; "
            f"overlap_geometry_candidates={overlap_geometry_count:,}; "
            f"overlap_bp_geometry_candidates={overlap_bp_geometry_candidate_count:,}; "
            f"overlap_reciprocal_rejected={overlap_reciprocal_rejected_count:,}; "
            f"models_with_same_source_edge={overlap_models_with_same_source_edge_count:,}; "
            f"overlap_policy_denied={overlap_policy_denied_count:,}"
        )

    return {
        "initial_retained_models": list(prune_result.get("initial_retained_models", [])),
        "initial_retained_rows": [dict(row) for row in prune_result.get("initial_retained_rows", [])],
        "retained_models": retained_models,
        "retained_rows": retained_rows,
        "removed_rows": removed_rows,
        "retained_selection_info_by_transcript": retained_info_by_id,
        "shadow_pair_rows": list(prune_result.get("shadow_pair_rows", [])),
        "shadow_candidate_pairs_before_singleton": int(prune_result.get("shadow_candidate_pairs_before_singleton", 0)),
        "shadow_candidate_pairs_after_singleton": int(prune_result.get("shadow_candidate_pairs_after_singleton", 0)),
        "shadow_candidate_pairs_skipped_after_singleton": int(prune_result.get("shadow_candidate_pairs_skipped_after_singleton", 0)),
        "singleton_rescued_by_monosublocus_overlap_count": rescue_count,
        "singleton_rescued_by_protein_evidence_count": protein_rescue_count,
        "singleton_rescued_by_transcript_evidence_count": transcript_rescue_count,
        "sequence_integrity_singleton_candidate_count": sequence_candidate_count,
        "sequence_integrity_singleton_pruned_count": sequence_pruned_count,
        "sequence_integrity_singleton_rescued_count": sequence_rescued_count,
        "singleton_overlap_rescue_geometry_candidate_count": overlap_geometry_count,
        "singleton_overlap_rescue_overlap_bp_geometry_candidate_count": overlap_bp_geometry_candidate_count,
        "singleton_overlap_rescue_reciprocal_rejected_count": overlap_reciprocal_rejected_count,
        "singleton_overlap_rescue_models_with_same_source_edge_count": overlap_models_with_same_source_edge_count,
        "singleton_overlap_rescue_policy_denied_count": overlap_policy_denied_count,
        "singleton_overlap_rescued_by_strong_local_evidence_count": overlap_rescued_by_strong_local_count,
        "singleton_overlap_rescued_without_adequate_evidence_count": overlap_rescued_without_adequate_evidence_count,
        "internal_stop_evidence_rescue_blocked_count": internal_stop_evidence_rescue_blocked_count,
        "internal_stop_evidence_rescue_applied_count": internal_stop_evidence_rescue_applied_count,
    }


def _low_score_evidence_adequacy(
    *,
    protein_evidence_supplied: bool,
    transcript_evidence_supplied: bool,
    protein_adequacy_required: bool,
    transcript_adequacy_required: bool,
    protein_overlap_model_fraction: float,
    protein_supported_model_fraction: float,
    transcript_scoring_eligible_model_fraction: float,
    transcript_supported_model_fraction: float,
    min_protein_overlap_model_fraction: float,
    min_protein_supported_model_fraction: float,
    min_transcript_scoring_eligible_model_fraction: float,
    min_transcript_supported_model_fraction: float,
) -> tuple[bool, bool, bool]:
    """Return (overall_external_adequacy, protein_adequacy, transcript_adequacy)."""
    protein_adequacy_passed = False
    if protein_evidence_supplied:
        protein_adequacy_passed = (
            (not protein_adequacy_required)
            or (
                protein_overlap_model_fraction + 1e-15 >= min_protein_overlap_model_fraction
                and protein_supported_model_fraction + 1e-15 >= min_protein_supported_model_fraction
            )
        )
    transcript_adequacy_passed = False
    if transcript_evidence_supplied:
        transcript_adequacy_passed = (
            (not transcript_adequacy_required)
            or (
                transcript_scoring_eligible_model_fraction + 1e-15 >= min_transcript_scoring_eligible_model_fraction
                and transcript_supported_model_fraction + 1e-15 >= min_transcript_supported_model_fraction
            )
        )
    return protein_adequacy_passed or transcript_adequacy_passed, protein_adequacy_passed, transcript_adequacy_passed


def low_score_singleton_fields(
    *,
    enabled: bool,
    protein_evidence_supplied: bool,
    transcript_evidence_supplied: bool,
    protein_adequacy_required: bool,
    transcript_adequacy_required: bool,
    external_adequacy_passed: bool,
    protein_adequacy_passed: bool,
    transcript_adequacy_passed: bool,
    protein_overlap_model_fraction: float,
    protein_supported_model_fraction: float,
    transcript_scoring_eligible_model_fraction: float,
    transcript_supported_model_fraction: float,
    min_protein_overlap_model_fraction: float,
    min_protein_supported_model_fraction: float,
    min_transcript_scoring_eligible_model_fraction: float,
    min_transcript_supported_model_fraction: float,
    max_model_structural_support_score_scaled: float,
    max_protein_evidence_score_scaled: float,
    max_transcript_evidence_score_scaled: float,
    max_model_selection_score_scaled: float,
    model_structural_support_score_scaled: float,
    protein_evidence_score_scaled: float,
    transcript_evidence_score_scaled: float,
    model_selection_score_scaled: float,
    use_transcript_cdna_context: bool = True,
    cdna_context_pruning_concern: bool = False,
    cdna_context_pruning_concern_class: str = ".",
    cdna_context_pruning_supporting_flags: str = ".",
    prune_route: str = ".",
    pruned: bool = False,
    reason: str = ".",
) -> dict[str, object]:
    return {
        "final_low_score_singleton_prune_enabled": "true" if enabled else "false",
        "final_low_score_singleton_pruned": "true" if pruned else "false",
        "final_low_score_singleton_prune_reason": reason,
        "final_low_score_singleton_prune_route": prune_route,
        "final_low_score_singleton_use_transcript_cdna_context": "true" if use_transcript_cdna_context else "false",
        "final_low_score_singleton_cdna_context_pruning_concern": "true" if cdna_context_pruning_concern else "false",
        "final_low_score_singleton_cdna_context_pruning_concern_class": cdna_context_pruning_concern_class,
        "final_low_score_singleton_cdna_context_pruning_supporting_flags": cdna_context_pruning_supporting_flags,
        "final_low_score_singleton_external_evidence_adequacy_passed": "true" if external_adequacy_passed else "false",
        "final_low_score_singleton_protein_evidence_supplied": "true" if protein_evidence_supplied else "false",
        "final_low_score_singleton_protein_adequacy_required": "true" if protein_adequacy_required else "false",
        "final_low_score_singleton_protein_adequacy_passed": "true" if protein_adequacy_passed else "false",
        "final_low_score_singleton_protein_overlap_model_fraction": protein_overlap_model_fraction,
        "final_low_score_singleton_protein_supported_model_fraction": protein_supported_model_fraction,
        "final_low_score_singleton_min_protein_overlap_model_fraction": min_protein_overlap_model_fraction,
        "final_low_score_singleton_min_protein_supported_model_fraction": min_protein_supported_model_fraction,
        "final_low_score_singleton_transcript_evidence_supplied": "true" if transcript_evidence_supplied else "false",
        "final_low_score_singleton_transcript_adequacy_required": "true" if transcript_adequacy_required else "false",
        "final_low_score_singleton_transcript_adequacy_passed": "true" if transcript_adequacy_passed else "false",
        "final_low_score_singleton_transcript_scoring_eligible_model_fraction": transcript_scoring_eligible_model_fraction,
        "final_low_score_singleton_transcript_supported_model_fraction": transcript_supported_model_fraction,
        "final_low_score_singleton_min_transcript_scoring_eligible_model_fraction": min_transcript_scoring_eligible_model_fraction,
        "final_low_score_singleton_min_transcript_supported_model_fraction": min_transcript_supported_model_fraction,
        "final_low_score_singleton_max_model_structural_support_score_scaled": max_model_structural_support_score_scaled,
        "final_low_score_singleton_max_protein_evidence_score_scaled": max_protein_evidence_score_scaled,
        "final_low_score_singleton_max_transcript_evidence_score_scaled": max_transcript_evidence_score_scaled,
        "final_low_score_singleton_max_model_selection_score_scaled": max_model_selection_score_scaled,
        "final_low_score_singleton_model_structural_support_score_scaled": model_structural_support_score_scaled,
        "final_low_score_singleton_protein_evidence_score_scaled": protein_evidence_score_scaled,
        "final_low_score_singleton_transcript_evidence_score_scaled": transcript_evidence_score_scaled,
        "final_low_score_singleton_model_selection_score_scaled": model_selection_score_scaled,
    }


def final_prune_low_score_singletons(
    prune_result: dict[str, object],
    subloci: list[Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    *,
    enabled: bool,
    protein_evidence_supplied: bool,
    transcript_evidence_supplied: bool,
    protein_adequacy_required: bool,
    transcript_adequacy_required: bool,
    protein_overlap_model_fraction: float,
    protein_supported_model_fraction: float,
    transcript_scoring_eligible_model_fraction: float,
    transcript_supported_model_fraction: float,
    min_protein_overlap_model_fraction: float,
    min_protein_supported_model_fraction: float,
    min_transcript_scoring_eligible_model_fraction: float,
    min_transcript_supported_model_fraction: float,
    max_model_structural_support_score_scaled: float,
    max_protein_evidence_score_scaled: float,
    max_transcript_evidence_score_scaled: float,
    max_model_selection_score_scaled: float,
    pruning_level_name: str,
    use_transcript_cdna_context: bool = True,
    progress: Optional[ProgressLogger] = None,
) -> dict[str, object]:
    """Optionally remove selected singleton-sublocus models by either low scores or cDNA context.

    The conventional route requires all applicable score thresholds.  The
    cDNA-context route is deliberately separate: it can remove a model with a
    high transcript score when the transcript evidence itself shows that the
    CDS occupies an outlyingly small fraction of a compatible cDNA.  That route
    requires transcript adequacy, a low model score, a low protein score when
    protein evidence is supplied, and an explicit precomputed context concern.
    Missing or uncalculable cDNA context is neutral.
    """
    sub_by_model_id: dict[str, Sublocus] = {m.output_transcript_id: sub for sub in subloci for m in sub.models}
    final_models: list[CodingModel] = list(prune_result.get("retained_models", []))  # type: ignore[arg-type]
    retained_rows_in: list[dict[str, object]] = list(prune_result.get("retained_rows", []))  # type: ignore[arg-type]
    rows_by_id = {str(row.get("output_transcript_id")): row for row in retained_rows_in}
    info_by_id: dict[str, dict[str, object]] = dict(prune_result.get("retained_selection_info_by_transcript", {}))  # type: ignore[arg-type]

    retained_models: list[CodingModel] = []
    retained_rows: list[dict[str, object]] = []
    retained_info_by_id: dict[str, dict[str, object]] = {}
    removed_rows: list[dict[str, object]] = list(prune_result.get("removed_rows", []))  # type: ignore[arg-type]

    external_adequacy_passed, protein_adequacy_passed, transcript_adequacy_passed = _low_score_evidence_adequacy(
        protein_evidence_supplied=protein_evidence_supplied,
        transcript_evidence_supplied=transcript_evidence_supplied,
        protein_adequacy_required=protein_adequacy_required,
        transcript_adequacy_required=transcript_adequacy_required,
        protein_overlap_model_fraction=protein_overlap_model_fraction,
        protein_supported_model_fraction=protein_supported_model_fraction,
        transcript_scoring_eligible_model_fraction=transcript_scoring_eligible_model_fraction,
        transcript_supported_model_fraction=transcript_supported_model_fraction,
        min_protein_overlap_model_fraction=min_protein_overlap_model_fraction,
        min_protein_supported_model_fraction=min_protein_supported_model_fraction,
        min_transcript_scoring_eligible_model_fraction=min_transcript_scoring_eligible_model_fraction,
        min_transcript_supported_model_fraction=min_transcript_supported_model_fraction,
    )
    candidates = 0
    removed_count = 0
    removed_score_route = 0
    removed_context_route = 0
    skipped_protein_rescue = 0
    skipped_transcript_rescue = 0
    protein_rescue_applied_here = 0
    transcript_rescue_applied_here = 0
    skipped_adequacy = 0
    skipped_thresholds = 0
    threshold_reason_counts: Counter[str] = Counter()

    if progress:
        progress.log(
            f"Final selected pruning: low-score singleton pruning enabled={str(enabled).lower()}; "
            f"use_transcript_cdna_context={str(use_transcript_cdna_context).lower()}; "
            f"external_adequacy_passed={str(external_adequacy_passed).lower()}; "
            f"protein_supplied={str(protein_evidence_supplied).lower()}; protein_adequacy_passed={str(protein_adequacy_passed).lower()}; "
            f"protein_overlap_model_fraction={protein_overlap_model_fraction:.6f}; "
            f"protein_supported_model_fraction={protein_supported_model_fraction:.6f}; "
            f"transcript_supplied={str(transcript_evidence_supplied).lower()}; transcript_adequacy_passed={str(transcript_adequacy_passed).lower()}; "
            f"transcript_scoring_eligible_model_fraction={transcript_scoring_eligible_model_fraction:.6f}; "
            f"transcript_supported_model_fraction={transcript_supported_model_fraction:.6f}"
        )

    for model in final_models:
        model_id = model.output_transcript_id
        sub = sub_by_model_id.get(model_id)
        row = dict(rows_by_id.get(model_id, stage_model_base_row(model, sub_by_model_id, model_structural_rows_by_transcript)))
        row.update({k: row.get(k, ".") for k in FINAL_PRUNE_EXTRA_COLUMNS})
        info = dict(info_by_id.get(model_id, row))
        is_singleton = bool(sub and len(sub.models) == 1)
        model_score = value_as_float(row.get("model_structural_support_score_scaled", "."), 0.0)
        protein_score = value_as_float(row.get("protein_evidence_score_scaled", "."), 0.0)
        transcript_score = value_as_float(row.get("transcript_evidence_score_scaled", "."), 0.0)
        model_selection_score = value_as_float(row.get("model_selection_score_scaled", "."), 0.0)
        context_concern = str(row.get("transcript_cdna_context_pruning_concern", "false")) == "true"
        context_class = str(row.get("transcript_cdna_context_pruning_concern_class", "."))
        context_flags = str(row.get("transcript_cdna_context_pruning_supporting_flags", "."))

        should_remove = False
        reason = "."
        prune_route = "."
        if enabled and is_singleton:
            candidates += 1
            protein_rescue_criteria_met = (
                str(row.get("final_protein_singleton_rescue_criteria_met", "false"))
                == "true"
            )
            transcript_rescue_criteria_met = (
                str(row.get("final_transcript_singleton_rescue_criteria_met", "false"))
                == "true"
            )
            if protein_rescue_criteria_met:
                skipped_protein_rescue += 1
                reason = "protected_by_protein_evidence_rescue"
                if str(row.get("final_protein_singleton_rescue_applied", "false")) != "true":
                    protein_rescue_applied_here += 1
                applied_reason = str(
                    row.get("final_protein_singleton_rescue_criteria_reason", ".")
                )
                row.update({
                    "final_singleton_rescued": "true",
                    "final_singleton_rescue_reason": "protein_evidence",
                    "final_protein_singleton_rescued": "true",
                    "final_protein_singleton_rescue_applied": "true",
                    "final_protein_singleton_rescue_applied_reason": applied_reason,
                })
            elif transcript_rescue_criteria_met:
                skipped_transcript_rescue += 1
                reason = "protected_by_transcript_evidence_rescue"
                if str(row.get("final_transcript_singleton_rescue_applied", "false")) != "true":
                    transcript_rescue_applied_here += 1
                applied_reason = str(
                    row.get("final_transcript_singleton_rescue_criteria_reason", ".")
                )
                row.update({
                    "final_singleton_rescued": "true",
                    "final_singleton_rescue_reason": "transcript_evidence",
                    "final_transcript_singleton_rescued": "true",
                    "final_transcript_singleton_rescue_applied": "true",
                    "final_transcript_singleton_rescue_applied_reason": applied_reason,
                })
            elif not external_adequacy_passed:
                skipped_adequacy += 1
                reason = "external_evidence_adequacy_failed"
            else:
                threshold_reasons: list[str] = []
                required_threshold_count = 2
                model_low = model_score <= max_model_structural_support_score_scaled + 1e-15
                protein_low = (not protein_evidence_supplied) or protein_score <= max_protein_evidence_score_scaled + 1e-15
                transcript_low = (not transcript_evidence_supplied) or transcript_score <= max_transcript_evidence_score_scaled + 1e-15
                model_selection_low = model_selection_score <= max_model_selection_score_scaled + 1e-15
                if model_low:
                    threshold_reasons.append("model_structural_support_score_scaled")
                if protein_evidence_supplied:
                    required_threshold_count += 1
                    if protein_low:
                        threshold_reasons.append("protein_evidence_score_scaled")
                if transcript_evidence_supplied:
                    required_threshold_count += 1
                    if transcript_low:
                        threshold_reasons.append("transcript_evidence_score_scaled")
                if model_selection_low:
                    threshold_reasons.append("model_selection_score_scaled")

                if len(threshold_reasons) == required_threshold_count:
                    should_remove = True
                    prune_route = "score_thresholds"
                    reason = "+".join(threshold_reasons)
                    removed_score_route += 1
                    threshold_reason_counts[reason] += 1
                else:
                    context_route_passed = bool(
                        use_transcript_cdna_context
                        and transcript_evidence_supplied
                        and transcript_adequacy_passed
                        and context_concern
                        and model_low
                        and protein_low
                    )
                    if context_route_passed:
                        should_remove = True
                        prune_route = "transcript_cdna_context"
                        context_reason_parts = ["model_structural_support_score_scaled"]
                        if protein_evidence_supplied:
                            context_reason_parts.append("protein_evidence_score_scaled")
                        context_reason_parts.append("transcript_cdna_context_pruning_concern")
                        reason = "+".join(context_reason_parts)
                        removed_context_route += 1
                        threshold_reason_counts[reason] += 1
                    else:
                        skipped_thresholds += 1
                        reason = "score_and_context_thresholds_not_met"

        low_score_fields = low_score_singleton_fields(
            enabled=enabled,
            protein_evidence_supplied=protein_evidence_supplied,
            transcript_evidence_supplied=transcript_evidence_supplied,
            protein_adequacy_required=protein_adequacy_required,
            transcript_adequacy_required=transcript_adequacy_required,
            external_adequacy_passed=external_adequacy_passed,
            protein_adequacy_passed=protein_adequacy_passed,
            transcript_adequacy_passed=transcript_adequacy_passed,
            protein_overlap_model_fraction=protein_overlap_model_fraction,
            protein_supported_model_fraction=protein_supported_model_fraction,
            transcript_scoring_eligible_model_fraction=transcript_scoring_eligible_model_fraction,
            transcript_supported_model_fraction=transcript_supported_model_fraction,
            min_protein_overlap_model_fraction=min_protein_overlap_model_fraction,
            min_protein_supported_model_fraction=min_protein_supported_model_fraction,
            min_transcript_scoring_eligible_model_fraction=min_transcript_scoring_eligible_model_fraction,
            min_transcript_supported_model_fraction=min_transcript_supported_model_fraction,
            max_model_structural_support_score_scaled=max_model_structural_support_score_scaled,
            max_protein_evidence_score_scaled=max_protein_evidence_score_scaled,
            max_transcript_evidence_score_scaled=max_transcript_evidence_score_scaled,
            max_model_selection_score_scaled=max_model_selection_score_scaled,
            model_structural_support_score_scaled=model_score,
            protein_evidence_score_scaled=protein_score,
            transcript_evidence_score_scaled=transcript_score,
            model_selection_score_scaled=model_selection_score,
            use_transcript_cdna_context=use_transcript_cdna_context,
            cdna_context_pruning_concern=context_concern,
            cdna_context_pruning_concern_class=context_class,
            cdna_context_pruning_supporting_flags=context_flags,
            prune_route=prune_route,
            pruned=should_remove,
            reason=reason if (enabled and is_singleton) else ".",
        )
        out_row = {**row, **low_score_fields, "final_prune_level": pruning_level_name}
        if should_remove:
            append_prune_class(out_row, "low_score_singleton")
            if prune_route == "transcript_cdna_context":
                append_prune_class(out_row, "transcript_cdna_context_singleton")
            removed_rows.append(out_row)
            removed_count += 1
        else:
            retained_models.append(model)
            retained_rows.append(out_row)
            retained_info_by_id[model_id] = {
                **info,
                **row,
                **low_score_fields,
                "final_prune_classes": out_row.get("final_prune_classes", "."),
            }

    if progress:
        progress.log(
            f"Final selected low-score singleton pruning complete: candidates={candidates:,}; removed={removed_count:,}; "
            f"removed_score_route={removed_score_route:,}; removed_context_route={removed_context_route:,}; "
            f"skipped_protein_rescue={skipped_protein_rescue:,}; skipped_transcript_rescue={skipped_transcript_rescue:,}; "
            f"skipped_adequacy={skipped_adequacy:,}; skipped_thresholds={skipped_thresholds:,}"
        )

    return {
        "initial_retained_models": list(prune_result.get("initial_retained_models", [])),
        "initial_retained_rows": [dict(row) for row in prune_result.get("initial_retained_rows", [])],
        "retained_models": retained_models,
        "retained_rows": retained_rows,
        "removed_rows": removed_rows,
        "retained_selection_info_by_transcript": retained_info_by_id,
        "shadow_pair_rows": list(prune_result.get("shadow_pair_rows", [])),
        "shadow_candidate_pairs_before_singleton": int(prune_result.get("shadow_candidate_pairs_before_singleton", 0)),
        "shadow_candidate_pairs_after_singleton": int(prune_result.get("shadow_candidate_pairs_after_singleton", 0)),
        "shadow_candidate_pairs_skipped_after_singleton": int(prune_result.get("shadow_candidate_pairs_skipped_after_singleton", 0)),
        "singleton_rescued_by_monosublocus_overlap_count": int(prune_result.get("singleton_rescued_by_monosublocus_overlap_count", 0)),
        "singleton_rescued_by_protein_evidence_count": int(prune_result.get("singleton_rescued_by_protein_evidence_count", 0)) + protein_rescue_applied_here,
        "singleton_rescued_by_transcript_evidence_count": int(prune_result.get("singleton_rescued_by_transcript_evidence_count", 0)) + transcript_rescue_applied_here,
        "singleton_overlap_rescue_geometry_candidate_count": int(prune_result.get("singleton_overlap_rescue_geometry_candidate_count", 0)),
        "singleton_overlap_rescue_overlap_bp_geometry_candidate_count": int(prune_result.get("singleton_overlap_rescue_overlap_bp_geometry_candidate_count", 0)),
        "singleton_overlap_rescue_reciprocal_rejected_count": int(prune_result.get("singleton_overlap_rescue_reciprocal_rejected_count", 0)),
        "singleton_overlap_rescue_models_with_same_source_edge_count": int(prune_result.get("singleton_overlap_rescue_models_with_same_source_edge_count", 0)),
        "singleton_overlap_rescue_policy_denied_count": int(prune_result.get("singleton_overlap_rescue_policy_denied_count", 0)),
        "singleton_overlap_rescued_by_strong_local_evidence_count": int(prune_result.get("singleton_overlap_rescued_by_strong_local_evidence_count", 0)),
        "singleton_overlap_rescued_without_adequate_evidence_count": int(prune_result.get("singleton_overlap_rescued_without_adequate_evidence_count", 0)),
        "internal_stop_evidence_rescue_blocked_count": int(prune_result.get("internal_stop_evidence_rescue_blocked_count", 0)),
        "internal_stop_evidence_rescue_applied_count": int(prune_result.get("internal_stop_evidence_rescue_applied_count", 0)),
        "low_score_singleton_pruning_enabled": enabled,
        "low_score_singleton_use_transcript_cdna_context": use_transcript_cdna_context,
        "low_score_singleton_external_adequacy_passed": external_adequacy_passed,
        "low_score_singleton_protein_evidence_supplied": protein_evidence_supplied,
        "low_score_singleton_transcript_evidence_supplied": transcript_evidence_supplied,
        "low_score_singleton_protein_adequacy_required": protein_adequacy_required,
        "low_score_singleton_transcript_adequacy_required": transcript_adequacy_required,
        "low_score_singleton_protein_adequacy_passed": protein_adequacy_passed,
        "low_score_singleton_transcript_adequacy_passed": transcript_adequacy_passed,
        "low_score_singleton_candidates": candidates,
        "low_score_singleton_removed_count": removed_count,
        "low_score_singleton_removed_score_route": removed_score_route,
        "low_score_singleton_removed_context_route": removed_context_route,
        "low_score_singleton_skipped_protein_rescue": skipped_protein_rescue,
        "low_score_singleton_skipped_transcript_rescue": skipped_transcript_rescue,
        "low_score_singleton_skipped_adequacy": skipped_adequacy,
        "low_score_singleton_skipped_thresholds": skipped_thresholds,
        "low_score_singleton_threshold_reason_counts": dict(threshold_reason_counts),
        "low_score_singleton_protein_overlap_model_fraction": protein_overlap_model_fraction,
        "low_score_singleton_protein_supported_model_fraction": protein_supported_model_fraction,
        "low_score_singleton_transcript_scoring_eligible_model_fraction": transcript_scoring_eligible_model_fraction,
        "low_score_singleton_transcript_supported_model_fraction": transcript_supported_model_fraction,
        "low_score_singleton_min_protein_overlap_model_fraction": min_protein_overlap_model_fraction,
        "low_score_singleton_min_protein_supported_model_fraction": min_protein_supported_model_fraction,
        "low_score_singleton_min_transcript_scoring_eligible_model_fraction": min_transcript_scoring_eligible_model_fraction,
        "low_score_singleton_min_transcript_supported_model_fraction": min_transcript_supported_model_fraction,
    }


def low_score_nonsingleton_fields(
    *,
    enabled: bool,
    scope: str,
    protein_evidence_supplied: bool,
    transcript_evidence_supplied: bool,
    protein_adequacy_required: bool,
    transcript_adequacy_required: bool,
    external_adequacy_passed: bool,
    protein_adequacy_passed: bool,
    transcript_adequacy_passed: bool,
    protein_overlap_model_fraction: float,
    protein_supported_model_fraction: float,
    transcript_scoring_eligible_model_fraction: float,
    transcript_supported_model_fraction: float,
    min_protein_overlap_model_fraction: float,
    min_protein_supported_model_fraction: float,
    min_transcript_scoring_eligible_model_fraction: float,
    min_transcript_supported_model_fraction: float,
    max_model_structural_support_score_scaled: float,
    max_protein_evidence_score_scaled: float,
    max_transcript_evidence_score_scaled: float,
    max_model_selection_score_scaled: float,
    model_structural_support_score_scaled: float,
    protein_evidence_score_scaled: float,
    transcript_evidence_score_scaled: float,
    model_selection_score_scaled: float,
    policy: str = ".",
    run_model_source_count: object = ".",
    sublocus_source_count: object = ".",
    source_depth_fraction: object = ".",
    effective_model_confidence_threshold: object = ".",
    model_confidence_tier: str = ".",
    use_transcript_cdna_context: bool = True,
    cdna_context_pruning_concern: bool = False,
    cdna_context_pruning_concern_class: str = ".",
    cdna_context_pruning_supporting_flags: str = ".",
    prune_route: str = ".",
    pruned: bool = False,
    reason: str = ".",
) -> dict[str, object]:
    return {
        "final_low_score_nonsingleton_prune_enabled": "true" if enabled else "false",
        "final_low_score_nonsingleton_prune_scope": scope if enabled else ".",
        "final_low_score_nonsingleton_pruned": "true" if pruned else "false",
        "final_low_score_nonsingleton_prune_reason": reason,
        "final_low_score_nonsingleton_policy": policy,
        "final_low_score_nonsingleton_run_model_source_count": run_model_source_count,
        "final_low_score_nonsingleton_sublocus_source_count": sublocus_source_count,
        "final_low_score_nonsingleton_source_depth_fraction": source_depth_fraction,
        "final_low_score_nonsingleton_effective_model_confidence_threshold": effective_model_confidence_threshold,
        "final_low_score_nonsingleton_model_confidence_tier": model_confidence_tier,
        "final_low_score_nonsingleton_prune_route": prune_route,
        "final_low_score_nonsingleton_use_transcript_cdna_context": "true" if use_transcript_cdna_context else "false",
        "final_low_score_nonsingleton_cdna_context_pruning_concern": "true" if cdna_context_pruning_concern else "false",
        "final_low_score_nonsingleton_cdna_context_pruning_concern_class": cdna_context_pruning_concern_class,
        "final_low_score_nonsingleton_cdna_context_pruning_supporting_flags": cdna_context_pruning_supporting_flags,
        "final_low_score_nonsingleton_external_evidence_adequacy_passed": "true" if external_adequacy_passed else "false",
        "final_low_score_nonsingleton_protein_evidence_supplied": "true" if protein_evidence_supplied else "false",
        "final_low_score_nonsingleton_protein_adequacy_required": "true" if protein_adequacy_required else "false",
        "final_low_score_nonsingleton_protein_adequacy_passed": "true" if protein_adequacy_passed else "false",
        "final_low_score_nonsingleton_protein_overlap_model_fraction": protein_overlap_model_fraction,
        "final_low_score_nonsingleton_protein_supported_model_fraction": protein_supported_model_fraction,
        "final_low_score_nonsingleton_min_protein_overlap_model_fraction": min_protein_overlap_model_fraction,
        "final_low_score_nonsingleton_min_protein_supported_model_fraction": min_protein_supported_model_fraction,
        "final_low_score_nonsingleton_transcript_evidence_supplied": "true" if transcript_evidence_supplied else "false",
        "final_low_score_nonsingleton_transcript_adequacy_required": "true" if transcript_adequacy_required else "false",
        "final_low_score_nonsingleton_transcript_adequacy_passed": "true" if transcript_adequacy_passed else "false",
        "final_low_score_nonsingleton_transcript_scoring_eligible_model_fraction": transcript_scoring_eligible_model_fraction,
        "final_low_score_nonsingleton_transcript_supported_model_fraction": transcript_supported_model_fraction,
        "final_low_score_nonsingleton_min_transcript_scoring_eligible_model_fraction": min_transcript_scoring_eligible_model_fraction,
        "final_low_score_nonsingleton_min_transcript_supported_model_fraction": min_transcript_supported_model_fraction,
        "final_low_score_nonsingleton_max_model_structural_support_score_scaled": max_model_structural_support_score_scaled,
        "final_low_score_nonsingleton_max_protein_evidence_score_scaled": max_protein_evidence_score_scaled,
        "final_low_score_nonsingleton_max_transcript_evidence_score_scaled": max_transcript_evidence_score_scaled,
        "final_low_score_nonsingleton_max_model_selection_score_scaled": max_model_selection_score_scaled,
        "final_low_score_nonsingleton_model_structural_support_score_scaled": model_structural_support_score_scaled,
        "final_low_score_nonsingleton_protein_evidence_score_scaled": protein_evidence_score_scaled,
        "final_low_score_nonsingleton_transcript_evidence_score_scaled": transcript_evidence_score_scaled,
        "final_low_score_nonsingleton_model_selection_score_scaled": model_selection_score_scaled,
    }


def model_confidence_nonsingleton_threshold(
    pruning_level_name: str,
    run_model_source_count: int,
    sublocus_source_count: int,
) -> tuple[Optional[float], str, float]:
    """Return the canonical fallback threshold, tier, and factual source depth.

    The model-source guard and depth escalation are deliberately defined here
    rather than through user-configurable evidence-aware thresholds.
    """
    source_depth_fraction = (
        sublocus_source_count / run_model_source_count
        if run_model_source_count > 0
        else 0.0
    )
    base_threshold = MODEL_CONFIDENCE_NONSINGLETON_THRESHOLDS.get(pruning_level_name)
    if base_threshold is None:
        return None, "inactive_pruning_level", source_depth_fraction
    if run_model_source_count < 3:
        return None, "insufficient_model_sources", source_depth_fraction

    if pruning_level_name in {"high", "very-high"}:
        if (
            run_model_source_count >= 7
            and sublocus_source_count == 2
            and source_depth_fraction <= 0.30 + 1e-15
        ):
            return (
                MODEL_CONFIDENCE_VERY_LOW_DEPTH_7PLUS_THRESHOLD,
                "very_low_depth_7plus",
                source_depth_fraction,
            )
        if (
            run_model_source_count >= 5
            and sublocus_source_count in {2, 3}
            and source_depth_fraction <= 0.50 + 1e-15
        ):
            return (
                MODEL_CONFIDENCE_LOW_DEPTH_5PLUS_THRESHOLD,
                "low_depth_5plus",
                source_depth_fraction,
            )
    return base_threshold, "base", source_depth_fraction


def final_prune_low_score_nonsingletons(
    prune_result: dict[str, object],
    subloci: list[Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    *,
    enabled: bool,
    scope: str,
    protein_evidence_supplied: bool,
    transcript_evidence_supplied: bool,
    protein_adequacy_required: bool,
    transcript_adequacy_required: bool,
    protein_overlap_model_fraction: float,
    protein_supported_model_fraction: float,
    transcript_scoring_eligible_model_fraction: float,
    transcript_supported_model_fraction: float,
    min_protein_overlap_model_fraction: float,
    min_protein_supported_model_fraction: float,
    min_transcript_scoring_eligible_model_fraction: float,
    min_transcript_supported_model_fraction: float,
    max_model_structural_support_score_scaled: float,
    max_protein_evidence_score_scaled: float,
    max_transcript_evidence_score_scaled: float,
    max_model_selection_score_scaled: float,
    pruning_level_name: str,
    run_model_source_count: int = 0,
    use_transcript_cdna_context: bool = True,
    progress: Optional[ProgressLogger] = None,
) -> dict[str, object]:
    """Prune nonsingletons using adequate evidence or model confidence fallback."""
    if scope not in {"backfill", "all"}:
        raise ValueError(f"Unsupported low-score nonsingleton pruning scope: {scope!r}")

    sub_by_model_id: dict[str, Sublocus] = {m.output_transcript_id: sub for sub in subloci for m in sub.models}
    final_models: list[CodingModel] = list(prune_result.get("retained_models", []))  # type: ignore[arg-type]
    retained_rows_in: list[dict[str, object]] = list(prune_result.get("retained_rows", []))  # type: ignore[arg-type]
    rows_by_id = {str(row.get("output_transcript_id")): row for row in retained_rows_in}
    info_by_id: dict[str, dict[str, object]] = dict(prune_result.get("retained_selection_info_by_transcript", {}))  # type: ignore[arg-type]

    retained_models: list[CodingModel] = []
    retained_rows: list[dict[str, object]] = []
    retained_info_by_id: dict[str, dict[str, object]] = {}
    removed_rows: list[dict[str, object]] = list(prune_result.get("removed_rows", []))  # type: ignore[arg-type]

    external_adequacy_passed, protein_adequacy_passed, transcript_adequacy_passed = _low_score_evidence_adequacy(
        protein_evidence_supplied=protein_evidence_supplied,
        transcript_evidence_supplied=transcript_evidence_supplied,
        protein_adequacy_required=protein_adequacy_required,
        transcript_adequacy_required=transcript_adequacy_required,
        protein_overlap_model_fraction=protein_overlap_model_fraction,
        protein_supported_model_fraction=protein_supported_model_fraction,
        transcript_scoring_eligible_model_fraction=transcript_scoring_eligible_model_fraction,
        transcript_supported_model_fraction=transcript_supported_model_fraction,
        min_protein_overlap_model_fraction=min_protein_overlap_model_fraction,
        min_protein_supported_model_fraction=min_protein_supported_model_fraction,
        min_transcript_scoring_eligible_model_fraction=min_transcript_scoring_eligible_model_fraction,
        min_transcript_supported_model_fraction=min_transcript_supported_model_fraction,
    )
    candidates = 0
    removed_count = 0
    removed_score_route = 0
    removed_context_route = 0
    removed_evidence_aware_score_route = 0
    removed_evidence_aware_context_route = 0
    removed_model_confidence_score_route = 0
    removed_model_confidence_context_route = 0
    skipped_scope = 0
    skipped_adequacy = 0
    skipped_model_source_guard = 0
    skipped_thresholds = 0
    removed_from_backfill = 0
    removed_from_nonbackfill = 0
    threshold_reason_counts: Counter[str] = Counter()
    policy_candidate_counts: Counter[str] = Counter()
    model_confidence_tier_counts: Counter[str] = Counter()
    policy = "evidence_aware" if external_adequacy_passed else "model_confidence_fallback"

    if progress:
        progress.log(
            f"Final selected pruning: low-score nonsingleton pruning enabled={str(enabled).lower()}; "
            f"scope={scope}; use_transcript_cdna_context={str(use_transcript_cdna_context).lower()}; "
            f"policy={policy}; run_model_source_count={run_model_source_count:,}; "
            f"external_adequacy_passed={str(external_adequacy_passed).lower()}; "
            f"protein_supplied={str(protein_evidence_supplied).lower()}; protein_adequacy_passed={str(protein_adequacy_passed).lower()}; "
            f"protein_overlap_model_fraction={protein_overlap_model_fraction:.6f}; "
            f"protein_supported_model_fraction={protein_supported_model_fraction:.6f}; "
            f"transcript_supplied={str(transcript_evidence_supplied).lower()}; transcript_adequacy_passed={str(transcript_adequacy_passed).lower()}; "
            f"transcript_scoring_eligible_model_fraction={transcript_scoring_eligible_model_fraction:.6f}; "
            f"transcript_supported_model_fraction={transcript_supported_model_fraction:.6f}"
        )

    for model in final_models:
        model_id = model.output_transcript_id
        sub = sub_by_model_id.get(model_id)
        row = dict(rows_by_id.get(model_id, stage_model_base_row(model, sub_by_model_id, model_structural_rows_by_transcript)))
        row.update({k: row.get(k, ".") for k in FINAL_PRUNE_EXTRA_COLUMNS})
        info = dict(info_by_id.get(model_id, row))
        is_nonsingleton = bool(sub and len(sub.models) > 1)
        is_backfill = str(row.get("final_selection_stage", info.get("final_selection_stage", "."))) == "backfill" or str(row.get("backfill_selected", info.get("backfill_selected", "."))) == "true"
        in_scope = is_nonsingleton and (scope == "all" or is_backfill)
        model_score = value_as_float(row.get("model_structural_support_score_scaled", "."), 0.0)
        protein_score = value_as_float(row.get("protein_evidence_score_scaled", "."), 0.0)
        transcript_score = value_as_float(row.get("transcript_evidence_score_scaled", "."), 0.0)
        model_selection_score = value_as_float(row.get("model_selection_score_scaled", "."), 0.0)
        context_concern = str(row.get("transcript_cdna_context_pruning_concern", "false")) == "true"
        context_class = str(row.get("transcript_cdna_context_pruning_concern_class", "."))
        context_flags = str(row.get("transcript_cdna_context_pruning_supporting_flags", "."))
        sublocus_source_count = len({member.source_name for member in sub.models}) if sub else 0
        fallback_threshold, fallback_tier, source_depth_fraction = model_confidence_nonsingleton_threshold(
            pruning_level_name,
            run_model_source_count,
            sublocus_source_count,
        )

        should_remove = False
        reason = "."
        prune_route = "."
        if enabled and is_nonsingleton:
            policy_candidate_counts[policy] += 1
            if policy == "model_confidence_fallback":
                model_confidence_tier_counts[fallback_tier] += 1
            if not in_scope:
                skipped_scope += 1
                reason = "outside_scope"
            else:
                candidates += 1
                if policy == "evidence_aware":
                    threshold_reasons: list[str] = []
                    required_threshold_count = 2
                    model_low = model_score <= max_model_structural_support_score_scaled + 1e-15
                    protein_low = (not protein_evidence_supplied) or protein_score <= max_protein_evidence_score_scaled + 1e-15
                    transcript_low = (not transcript_evidence_supplied) or transcript_score <= max_transcript_evidence_score_scaled + 1e-15
                    model_selection_low = model_selection_score <= max_model_selection_score_scaled + 1e-15
                    if model_low:
                        threshold_reasons.append("model_structural_support_score_scaled")
                    if protein_evidence_supplied:
                        required_threshold_count += 1
                        if protein_low:
                            threshold_reasons.append("protein_evidence_score_scaled")
                    if transcript_evidence_supplied:
                        required_threshold_count += 1
                        if transcript_low:
                            threshold_reasons.append("transcript_evidence_score_scaled")
                    if model_selection_low:
                        threshold_reasons.append("model_selection_score_scaled")

                    if len(threshold_reasons) == required_threshold_count:
                        should_remove = True
                        prune_route = "evidence_aware_score_thresholds"
                        reason = "+".join(threshold_reasons)
                        removed_score_route += 1
                        removed_evidence_aware_score_route += 1
                        threshold_reason_counts[reason] += 1
                    else:
                        context_route_passed = bool(
                            use_transcript_cdna_context
                            and transcript_evidence_supplied
                            and transcript_adequacy_passed
                            and context_concern
                            and model_low
                            and protein_low
                        )
                        if context_route_passed:
                            should_remove = True
                            prune_route = "evidence_aware_transcript_cdna_context"
                            context_reason_parts = ["model_structural_support_score_scaled"]
                            if protein_evidence_supplied:
                                context_reason_parts.append("protein_evidence_score_scaled")
                            context_reason_parts.append("transcript_cdna_context_pruning_concern")
                            reason = "+".join(context_reason_parts)
                            removed_context_route += 1
                            removed_evidence_aware_context_route += 1
                            threshold_reason_counts[reason] += 1
                        else:
                            skipped_thresholds += 1
                            reason = "score_and_context_thresholds_not_met"
                elif fallback_threshold is None:
                    skipped_model_source_guard += 1
                    reason = fallback_tier
                else:
                    model_low = model_score <= fallback_threshold + 1e-15
                    protein_low = (
                        (not protein_evidence_supplied)
                        or protein_score <= max_protein_evidence_score_scaled + 1e-15
                    )
                    transcript_low = (
                        (not transcript_evidence_supplied)
                        or transcript_score <= max_transcript_evidence_score_scaled + 1e-15
                    )
                    if model_low and protein_low and transcript_low:
                        should_remove = True
                        prune_route = "model_confidence_score_thresholds"
                        reason_parts = ["model_structural_support_score_scaled"]
                        if protein_evidence_supplied:
                            reason_parts.append("protein_evidence_score_scaled")
                        if transcript_evidence_supplied:
                            reason_parts.append("transcript_evidence_score_scaled")
                        reason = "+".join(reason_parts)
                        removed_score_route += 1
                        removed_model_confidence_score_route += 1
                        threshold_reason_counts[reason] += 1
                    else:
                        context_route_passed = bool(
                            use_transcript_cdna_context
                            and transcript_evidence_supplied
                            and context_concern
                            and model_low
                            and protein_low
                        )
                        if context_route_passed:
                            should_remove = True
                            prune_route = "model_confidence_transcript_cdna_context"
                            context_reason_parts = ["model_structural_support_score_scaled"]
                            if protein_evidence_supplied:
                                context_reason_parts.append("protein_evidence_score_scaled")
                            context_reason_parts.append("transcript_cdna_context_pruning_concern")
                            reason = "+".join(context_reason_parts)
                            removed_context_route += 1
                            removed_model_confidence_context_route += 1
                            threshold_reason_counts[reason] += 1
                        else:
                            skipped_thresholds += 1
                            reason = "score_and_context_thresholds_not_met"

        low_score_fields = low_score_nonsingleton_fields(
            enabled=enabled,
            scope=scope,
            protein_evidence_supplied=protein_evidence_supplied,
            transcript_evidence_supplied=transcript_evidence_supplied,
            protein_adequacy_required=protein_adequacy_required,
            transcript_adequacy_required=transcript_adequacy_required,
            external_adequacy_passed=external_adequacy_passed,
            protein_adequacy_passed=protein_adequacy_passed,
            transcript_adequacy_passed=transcript_adequacy_passed,
            protein_overlap_model_fraction=protein_overlap_model_fraction,
            protein_supported_model_fraction=protein_supported_model_fraction,
            transcript_scoring_eligible_model_fraction=transcript_scoring_eligible_model_fraction,
            transcript_supported_model_fraction=transcript_supported_model_fraction,
            min_protein_overlap_model_fraction=min_protein_overlap_model_fraction,
            min_protein_supported_model_fraction=min_protein_supported_model_fraction,
            min_transcript_scoring_eligible_model_fraction=min_transcript_scoring_eligible_model_fraction,
            min_transcript_supported_model_fraction=min_transcript_supported_model_fraction,
            max_model_structural_support_score_scaled=max_model_structural_support_score_scaled,
            max_protein_evidence_score_scaled=max_protein_evidence_score_scaled,
            max_transcript_evidence_score_scaled=max_transcript_evidence_score_scaled,
            max_model_selection_score_scaled=max_model_selection_score_scaled,
            model_structural_support_score_scaled=model_score,
            protein_evidence_score_scaled=protein_score,
            transcript_evidence_score_scaled=transcript_score,
            model_selection_score_scaled=model_selection_score,
            policy=policy if is_nonsingleton else ".",
            run_model_source_count=run_model_source_count if is_nonsingleton else ".",
            sublocus_source_count=sublocus_source_count if is_nonsingleton else ".",
            source_depth_fraction=source_depth_fraction if is_nonsingleton else ".",
            effective_model_confidence_threshold=(
                fallback_threshold
                if is_nonsingleton and policy == "model_confidence_fallback" and fallback_threshold is not None
                else "."
            ),
            model_confidence_tier=(
                fallback_tier if is_nonsingleton and policy == "model_confidence_fallback" else "not_applicable"
            ),
            use_transcript_cdna_context=use_transcript_cdna_context,
            cdna_context_pruning_concern=context_concern,
            cdna_context_pruning_concern_class=context_class,
            cdna_context_pruning_supporting_flags=context_flags,
            prune_route=prune_route,
            pruned=should_remove,
            reason=reason if (enabled and is_nonsingleton) else ".",
        )
        out_row = {**row, **low_score_fields, "final_prune_level": pruning_level_name}
        if should_remove:
            append_prune_class(out_row, "low_score_nonsingleton")
            if prune_route in {
                "evidence_aware_transcript_cdna_context",
                "model_confidence_transcript_cdna_context",
            }:
                append_prune_class(out_row, "transcript_cdna_context_nonsingleton")
            removed_rows.append(out_row)
            removed_count += 1
            if is_backfill:
                removed_from_backfill += 1
            else:
                removed_from_nonbackfill += 1
        else:
            retained_models.append(model)
            retained_rows.append(out_row)
            retained_info_by_id[model_id] = {**info, **low_score_fields, "final_prune_classes": out_row.get("final_prune_classes", ".")}

    if progress:
        progress.log(
            f"Final selected low-score nonsingleton pruning complete: scope={scope}; candidates={candidates:,}; "
            f"policy={policy}; run_model_source_count={run_model_source_count:,}; "
            f"removed={removed_count:,}; removed_score_route={removed_score_route:,}; removed_context_route={removed_context_route:,}; "
            f"skipped_scope={skipped_scope:,}; "
            f"skipped_model_source_guard={skipped_model_source_guard:,}; skipped_thresholds={skipped_thresholds:,}; "
            f"removed_from_backfill={removed_from_backfill:,}; removed_from_nonbackfill={removed_from_nonbackfill:,}"
        )

    out_result = dict(prune_result)
    out_result.update({
        "initial_retained_models": list(prune_result.get("initial_retained_models", [])),
        "initial_retained_rows": [dict(row) for row in prune_result.get("initial_retained_rows", [])],
        "retained_models": retained_models,
        "retained_rows": retained_rows,
        "removed_rows": removed_rows,
        "retained_selection_info_by_transcript": retained_info_by_id,
        "shadow_pair_rows": list(prune_result.get("shadow_pair_rows", [])),
        "shadow_candidate_pairs_before_singleton": int(prune_result.get("shadow_candidate_pairs_before_singleton", 0)),
        "shadow_candidate_pairs_after_singleton": int(prune_result.get("shadow_candidate_pairs_after_singleton", 0)),
        "shadow_candidate_pairs_skipped_after_singleton": int(prune_result.get("shadow_candidate_pairs_skipped_after_singleton", 0)),
        "singleton_rescued_by_monosublocus_overlap_count": int(prune_result.get("singleton_rescued_by_monosublocus_overlap_count", 0)),
        "singleton_rescued_by_protein_evidence_count": int(prune_result.get("singleton_rescued_by_protein_evidence_count", 0)),
        "singleton_rescued_by_transcript_evidence_count": int(prune_result.get("singleton_rescued_by_transcript_evidence_count", 0)),
        "low_score_nonsingleton_pruning_enabled": enabled,
        "low_score_nonsingleton_scope": scope,
        "low_score_nonsingleton_use_transcript_cdna_context": use_transcript_cdna_context,
        "low_score_nonsingleton_policy": policy,
        "low_score_nonsingleton_run_model_source_count": run_model_source_count,
        "low_score_nonsingleton_external_adequacy_passed": external_adequacy_passed,
        "low_score_nonsingleton_protein_evidence_supplied": protein_evidence_supplied,
        "low_score_nonsingleton_transcript_evidence_supplied": transcript_evidence_supplied,
        "low_score_nonsingleton_protein_adequacy_required": protein_adequacy_required,
        "low_score_nonsingleton_transcript_adequacy_required": transcript_adequacy_required,
        "low_score_nonsingleton_protein_adequacy_passed": protein_adequacy_passed,
        "low_score_nonsingleton_transcript_adequacy_passed": transcript_adequacy_passed,
        "low_score_nonsingleton_candidates": candidates,
        "low_score_nonsingleton_removed_count": removed_count,
        "low_score_nonsingleton_removed_score_route": removed_score_route,
        "low_score_nonsingleton_removed_context_route": removed_context_route,
        "low_score_nonsingleton_removed_evidence_aware_score_route": removed_evidence_aware_score_route,
        "low_score_nonsingleton_removed_evidence_aware_context_route": removed_evidence_aware_context_route,
        "low_score_nonsingleton_removed_model_confidence_score_route": removed_model_confidence_score_route,
        "low_score_nonsingleton_removed_model_confidence_context_route": removed_model_confidence_context_route,
        "low_score_nonsingleton_skipped_scope": skipped_scope,
        "low_score_nonsingleton_skipped_adequacy": skipped_adequacy,
        "low_score_nonsingleton_skipped_model_source_guard": skipped_model_source_guard,
        "low_score_nonsingleton_skipped_thresholds": skipped_thresholds,
        "low_score_nonsingleton_removed_from_backfill": removed_from_backfill,
        "low_score_nonsingleton_removed_from_nonbackfill": removed_from_nonbackfill,
        "low_score_nonsingleton_threshold_reason_counts": dict(threshold_reason_counts),
        "low_score_nonsingleton_policy_candidate_counts": dict(policy_candidate_counts),
        "low_score_nonsingleton_model_confidence_tier_counts": dict(model_confidence_tier_counts),
        "low_score_nonsingleton_protein_overlap_model_fraction": protein_overlap_model_fraction,
        "low_score_nonsingleton_protein_supported_model_fraction": protein_supported_model_fraction,
        "low_score_nonsingleton_transcript_scoring_eligible_model_fraction": transcript_scoring_eligible_model_fraction,
        "low_score_nonsingleton_transcript_supported_model_fraction": transcript_supported_model_fraction,
        "low_score_nonsingleton_min_protein_overlap_model_fraction": min_protein_overlap_model_fraction,
        "low_score_nonsingleton_min_protein_supported_model_fraction": min_protein_supported_model_fraction,
        "low_score_nonsingleton_min_transcript_scoring_eligible_model_fraction": min_transcript_scoring_eligible_model_fraction,
        "low_score_nonsingleton_min_transcript_supported_model_fraction": min_transcript_supported_model_fraction,
    })
    return out_result


def final_prune_sequence_integrity_nonsingletons(
    prune_result: dict[str, object],
    subloci: list[Sublocus],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    *,
    enabled: bool,
    scope: str,
    pruning_level_name: str,
    protect_protein_supported: bool = False,
    protein_evidence_rows_by_transcript: Optional[dict[str, dict[str, object]]] = None,
    min_protein_evidence_score_scaled: float = 0.75,
    min_protein_exact_chain_source_count: int = 1,
    protect_transcript_supported: bool = False,
    transcript_evidence_rows_by_transcript: Optional[dict[str, dict[str, object]]] = None,
    min_transcript_evidence_score_scaled: float = 0.75,
    min_transcript_supporting_source_count: int = 1,
    min_transcript_scoring_eligible_count: int = 1,
    min_transcript_cdna_context_fraction: float = 0.30,
    transcript_require_stranded_support: bool = False,
    transcript_block_on_cdna_context_concern: bool = True,
    allow_internal_stop_evidence_rescue: bool = False,
    progress: Optional[ProgressLogger] = None,
) -> dict[str, object]:
    """Prune sequence-concerning nonsingletons in backfill or global scope.

    Preset backfill scope uses the standard concern class. All scope is a true
    superset: it retains the standard concern rule for backfilled models and
    additionally applies the severe concern class to every selected
    nonsingleton. Models meeting the configured strong protein or transcript
    singleton-rescue criteria are protected.
    """
    if scope not in {"backfill", "all"}:
        raise ValueError(f"Unsupported sequence-integrity nonsingleton pruning scope: {scope!r}")

    sub_by_model_id: dict[str, Sublocus] = {
        model.output_transcript_id: sub for sub in subloci for model in sub.models
    }
    final_models: list[CodingModel] = list(prune_result.get("retained_models", []))  # type: ignore[arg-type]
    retained_rows_in: list[dict[str, object]] = list(prune_result.get("retained_rows", []))  # type: ignore[arg-type]
    rows_by_id = {str(row.get("output_transcript_id")): row for row in retained_rows_in}
    info_by_id: dict[str, dict[str, object]] = dict(
        prune_result.get("retained_selection_info_by_transcript", {})
    )  # type: ignore[arg-type]
    protein_evidence_rows_by_transcript = protein_evidence_rows_by_transcript or {}
    transcript_evidence_rows_by_transcript = transcript_evidence_rows_by_transcript or {}

    retained_models: list[CodingModel] = []
    retained_rows: list[dict[str, object]] = []
    retained_info_by_id: dict[str, dict[str, object]] = {}
    removed_rows: list[dict[str, object]] = list(prune_result.get("removed_rows", []))  # type: ignore[arg-type]

    candidates = 0
    removed_count = 0
    skipped_scope = 0
    skipped_no_concern = 0
    protected_by_protein = 0
    protected_by_transcript = 0
    internal_stop_evidence_rescue_blocked_count = int(
        prune_result.get("internal_stop_evidence_rescue_blocked_count", 0)
    )
    internal_stop_evidence_rescue_applied_count = int(
        prune_result.get("internal_stop_evidence_rescue_applied_count", 0)
    )
    removed_from_backfill = 0
    removed_from_nonbackfill = 0
    backfill_standard_candidates = 0
    all_severe_candidates = 0
    removed_backfill_standard_route = 0
    removed_all_severe_route = 0
    reason_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()

    if progress:
        progress.log(
            "Final selected pruning: sequence-integrity nonsingleton pruning "
            f"enabled={str(enabled).lower()}; scope={scope}; models={len(final_models):,}; "
            f"protein_protection={str(protect_protein_supported).lower()}; "
            f"transcript_protection={str(protect_transcript_supported).lower()}"
            f"; allow_internal_stop_evidence_rescue={str(allow_internal_stop_evidence_rescue).lower()}"
        )

    for model in final_models:
        model_id = model.output_transcript_id
        sub = sub_by_model_id.get(model_id)
        row = dict(
            rows_by_id.get(
                model_id,
                stage_model_base_row(model, sub_by_model_id, model_structural_rows_by_transcript),
            )
        )
        row.update({key: row.get(key, ".") for key in FINAL_PRUNE_EXTRA_COLUMNS})
        info = dict(info_by_id.get(model_id, row))
        structural_row = model_structural_rows_by_transcript.get(model_id, {})
        sequence_fields = sequence_integrity_pruning_fields({**structural_row, **row})
        row.update(sequence_fields)
        is_nonsingleton = bool(sub and len(sub.models) > 1)
        is_backfill = (
            str(row.get("final_selection_stage", info.get("final_selection_stage", ".")))
            == "backfill"
            or str(row.get("backfill_selected", info.get("backfill_selected", "."))) == "true"
        )
        standard_concern = (
            sequence_fields["sequence_integrity_pruning_concern"] == "true"
        )
        severe_concern = (
            sequence_fields["sequence_integrity_severe_pruning_concern"] == "true"
        )
        prune_route = "."
        if is_nonsingleton:
            if is_backfill and standard_concern:
                prune_route = "backfill_standard"
            elif scope == "all" and severe_concern:
                prune_route = "all_severe"
        in_scope = is_nonsingleton and (scope == "all" or is_backfill)
        concern_required = prune_route != "."
        sequence_candidate = bool(enabled and concern_required)
        protein_protected = False
        transcript_protected = False
        should_remove = False
        reason = "."
        prior_internal_stop_audit = row
        internal_stop_audit = internal_stop_evidence_audit_fields(
            sequence_fields,
            allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
            protein_criteria_met=False,
            transcript_criteria_met=False,
            evidence_override_candidate=False,
        )

        if enabled and is_nonsingleton:
            if not in_scope:
                skipped_scope += 1
                reason = "outside_scope"
            elif not concern_required:
                skipped_no_concern += 1
                reason = "sequence_integrity_concern_not_met"
            else:
                candidates += 1
                if prune_route == "backfill_standard":
                    backfill_standard_candidates += 1
                elif prune_route == "all_severe":
                    all_severe_candidates += 1
                protein_fields = singleton_protein_rescue_fields(
                    model_id,
                    protein_evidence_rows_by_transcript,
                    rescue_protein_supported_singletons=protect_protein_supported,
                    min_score_scaled=min_protein_evidence_score_scaled,
                    min_exact_chain_source_count=min_protein_exact_chain_source_count,
                )
                transcript_fields = singleton_transcript_rescue_fields(
                    model_id,
                    transcript_evidence_rows_by_transcript,
                    rescue_transcript_supported_singletons=protect_transcript_supported,
                    min_score_scaled=min_transcript_evidence_score_scaled,
                    min_supporting_source_count=min_transcript_supporting_source_count,
                    min_scoring_eligible_count=min_transcript_scoring_eligible_count,
                    min_cdna_context_fraction=min_transcript_cdna_context_fraction,
                    require_stranded_support=transcript_require_stranded_support,
                    block_on_cdna_context_concern=transcript_block_on_cdna_context_concern,
                )
                protein_criteria_met = (
                    protect_protein_supported
                    and str(protein_fields.get("final_protein_singleton_rescue_criteria_met", "false"))
                    == "true"
                )
                transcript_criteria_met = (
                    protect_transcript_supported
                    and str(transcript_fields.get("final_transcript_singleton_rescue_criteria_met", "false"))
                    == "true"
                )
                internal_stop_policy_blocks = internal_stop_blocks_evidence_rescue(
                    sequence_fields,
                    allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
                )
                internal_stop_audit = internal_stop_evidence_audit_fields(
                    sequence_fields,
                    allow_internal_stop_evidence_rescue=allow_internal_stop_evidence_rescue,
                    protein_criteria_met=protein_criteria_met,
                    transcript_criteria_met=transcript_criteria_met,
                    evidence_override_candidate=True,
                )
                if internal_stop_audit["final_internal_stop_evidence_rescue_blocked"] == "true":
                    internal_stop_evidence_rescue_blocked_count += 1
                protein_protected = protein_criteria_met and not internal_stop_policy_blocks
                transcript_protected = transcript_criteria_met and not internal_stop_policy_blocks
                if protein_protected:
                    protected_by_protein += 1
                    reason = "protected_by_protein_evidence"
                    if str(internal_stop_audit["final_model_has_internal_stop_codon"]) == "true":
                        internal_stop_audit.update({
                            "final_internal_stop_evidence_rescue_applied": "true",
                            "final_internal_stop_evidence_rescue_applied_route": "protein_evidence",
                        })
                        internal_stop_evidence_rescue_applied_count += 1
                elif transcript_protected:
                    protected_by_transcript += 1
                    reason = "protected_by_transcript_evidence"
                    if str(internal_stop_audit["final_model_has_internal_stop_codon"]) == "true":
                        internal_stop_audit.update({
                            "final_internal_stop_evidence_rescue_applied": "true",
                            "final_internal_stop_evidence_rescue_applied_route": "transcript_evidence",
                        })
                        internal_stop_evidence_rescue_applied_count += 1
                else:
                    should_remove = True
                    reasons_key = (
                        "sequence_integrity_pruning_reasons"
                        if prune_route == "backfill_standard"
                        else "sequence_integrity_severe_pruning_reasons"
                    )
                    concern_reasons = str(sequence_fields.get(reasons_key, "."))
                    reason_scope = (
                        "backfill_scope"
                        if prune_route == "backfill_standard"
                        else "all_scope"
                    )
                    reason = f"sequence_integrity_{reason_scope}:{concern_reasons}"
                    reason_counts[concern_reasons] += 1
                    source_counts[model.source_name] += 1

        internal_stop_audit = merge_internal_stop_evidence_audit_fields(
            prior_internal_stop_audit,
            internal_stop_audit,
        )

        sequence_prune_fields = {
            **sequence_fields,
            "final_prune_level": pruning_level_name,
            "final_sequence_integrity_nonsingleton_prune_enabled": (
                "true" if enabled else "false"
            ),
            "final_sequence_integrity_nonsingleton_prune_scope": scope,
            "final_sequence_integrity_nonsingleton_candidate": (
                "true" if sequence_candidate else "false"
            ),
            "final_sequence_integrity_nonsingleton_pruned": (
                "true" if should_remove else "false"
            ),
            "final_sequence_integrity_nonsingleton_prune_route": (
                prune_route if enabled and is_nonsingleton else "."
            ),
            "final_sequence_integrity_nonsingleton_prune_reason": (
                reason if enabled and is_nonsingleton else "."
            ),
            "final_sequence_integrity_nonsingleton_protein_protected": (
                "true" if protein_protected else "false"
            ),
            "final_sequence_integrity_nonsingleton_transcript_protected": (
                "true" if transcript_protected else "false"
            ),
            **internal_stop_audit,
        }
        out_row = {**row, **sequence_prune_fields}
        if should_remove:
            append_prune_class(out_row, "sequence_integrity_nonsingleton")
            removed_rows.append(out_row)
            removed_count += 1
            if is_backfill:
                removed_from_backfill += 1
            else:
                removed_from_nonbackfill += 1
            if prune_route == "backfill_standard":
                removed_backfill_standard_route += 1
            elif prune_route == "all_severe":
                removed_all_severe_route += 1
        else:
            retained_models.append(model)
            retained_rows.append(out_row)
            retained_info_by_id[model_id] = {
                **info,
                **sequence_prune_fields,
                "final_prune_classes": out_row.get("final_prune_classes", "."),
            }

    if progress:
        progress.log(
            "Final selected sequence-integrity nonsingleton pruning complete: "
            f"scope={scope}; candidates={candidates:,}; removed={removed_count:,}; "
            f"protected_by_protein={protected_by_protein:,}; "
            f"protected_by_transcript={protected_by_transcript:,}; "
            f"skipped_scope={skipped_scope:,}; skipped_no_concern={skipped_no_concern:,}; "
            f"backfill_standard_candidates={backfill_standard_candidates:,}; "
            f"all_severe_candidates={all_severe_candidates:,}; "
            f"removed_backfill_standard_route={removed_backfill_standard_route:,}; "
            f"removed_all_severe_route={removed_all_severe_route:,}; "
            f"removed_from_backfill={removed_from_backfill:,}; "
            f"removed_from_nonbackfill={removed_from_nonbackfill:,}"
        )

    out_result = dict(prune_result)
    out_result.update(
        {
            "retained_models": retained_models,
            "retained_rows": retained_rows,
            "removed_rows": removed_rows,
            "retained_selection_info_by_transcript": retained_info_by_id,
            "sequence_integrity_nonsingleton_pruning_enabled": enabled,
            "sequence_integrity_nonsingleton_scope": scope,
            "sequence_integrity_nonsingleton_candidates": candidates,
            "sequence_integrity_nonsingleton_removed_count": removed_count,
            "sequence_integrity_nonsingleton_protected_by_protein": protected_by_protein,
            "sequence_integrity_nonsingleton_protected_by_transcript": protected_by_transcript,
            "internal_stop_evidence_rescue_blocked_count": internal_stop_evidence_rescue_blocked_count,
            "internal_stop_evidence_rescue_applied_count": internal_stop_evidence_rescue_applied_count,
            "sequence_integrity_nonsingleton_skipped_scope": skipped_scope,
            "sequence_integrity_nonsingleton_skipped_no_concern": skipped_no_concern,
            "sequence_integrity_nonsingleton_backfill_standard_candidates": backfill_standard_candidates,
            "sequence_integrity_nonsingleton_all_severe_candidates": all_severe_candidates,
            "sequence_integrity_nonsingleton_removed_backfill_standard_route": removed_backfill_standard_route,
            "sequence_integrity_nonsingleton_removed_all_severe_route": removed_all_severe_route,
            "sequence_integrity_nonsingleton_removed_from_backfill": removed_from_backfill,
            "sequence_integrity_nonsingleton_removed_from_nonbackfill": removed_from_nonbackfill,
            "sequence_integrity_nonsingleton_reason_counts": dict(reason_counts),
            "sequence_integrity_nonsingleton_source_counts": dict(source_counts),
        }
    )
    return out_result


def shadow_model_score_tuple(model: CodingModel, row: dict[str, object]) -> tuple[float, float, float, float, int, str]:
    return (
        value_as_float(row.get("model_selection_score_scaled")),
        value_as_float(row.get("model_selection_score")),
        value_as_float(row.get("model_completeness_score")),
        -value_as_float(row.get("intrinsic_penalty_score")),
        model.cds_length,
        model.output_transcript_id,
    )


def find_opposite_strand_shadow_candidate_pairs(
    models: list[CodingModel],
    rows_by_id: dict[str, dict[str, object]],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    *,
    min_overlap_bp: int,
    min_loser_cds_overlap_fraction: float,
    progress: Optional[ProgressLogger] = None,
    progress_label: str = "Final selected pruning shadow scan",
) -> dict[tuple[str, str], dict[str, object]]:
    """Find opposite-strand selected-model pairs meeting shadow overlap candidate thresholds.

    Candidate status is deliberately based only on CDS-exon overlap bp and
    overlap fraction of the lower-scoring model. Eligibility filters such as
    source scope, loser single-exon status, and score margin are applied later
    so the log can distinguish candidate discovery from actual pruning.
    """
    # Build a simple span-bin index by chromosome and strand to avoid all-vs-all
    # opposite-strand comparisons. CDS-exon overlap is still calculated exactly.
    bin_size = 100_000
    by_chrom_strand: dict[tuple[str, str], list[CodingModel]] = defaultdict(list)
    for model in models:
        by_chrom_strand[(model.chrom, model.strand)].append(model)

    index: dict[tuple[str, str, int], list[CodingModel]] = defaultdict(list)
    for model in models:
        first_bin = model.cds_start0 // bin_size
        last_bin = max(model.cds_end0 - 1, model.cds_start0) // bin_size
        for bin_id in range(first_bin, last_bin + 1):
            index[(model.chrom, model.strand, bin_id)].append(model)

    candidate_pairs: dict[tuple[str, str], dict[str, object]] = {}
    for chrom in sorted({m.chrom for m in models}):
        plus_models = by_chrom_strand.get((chrom, "+"), [])
        minus_models = by_chrom_strand.get((chrom, "-"), [])
        if not plus_models or not minus_models:
            continue
        minus_by_id = {m.output_transcript_id: m for m in minus_models}
        ordered_plus = sorted(plus_models, key=lambda m: (m.cds_start0, m.cds_end0, m.output_transcript_id))
        for plus_index, a in enumerate(ordered_plus, start=1):
            if progress:
                progress.every(plus_index, len(ordered_plus), f"{progress_label} {chrom}", level="basic")
            candidate_ids: set[str] = set()
            first_bin = a.cds_start0 // bin_size
            last_bin = max(a.cds_end0 - 1, a.cds_start0) // bin_size
            for bin_id in range(first_bin, last_bin + 1):
                for b in index.get((chrom, "-", bin_id), []):
                    candidate_ids.add(b.output_transcript_id)
            for b_id in sorted(candidate_ids):
                b = minus_by_id.get(b_id)
                if b is None:
                    continue
                if b.cds_end0 <= a.cds_start0 or b.cds_start0 >= a.cds_end0:
                    continue
                metrics = cds_exon_overlap_metrics(a, b)
                overlap_bp = int(metrics["cds_exon_overlap_bp"])
                if overlap_bp < min_overlap_bp:
                    continue
                row_a = rows_by_id.get(a.output_transcript_id, model_structural_rows_by_transcript.get(a.output_transcript_id, {}))
                row_b = rows_by_id.get(b.output_transcript_id, model_structural_rows_by_transcript.get(b.output_transcript_id, {}))
                score_a = shadow_model_score_tuple(a, row_a)
                score_b = shadow_model_score_tuple(b, row_b)
                if score_a >= score_b:
                    winner, loser = a, b
                    winner_row, loser_row = row_a, row_b
                    winner_frac = float(metrics["overlap_fraction_of_a"])
                    loser_frac = float(metrics["overlap_fraction_of_b"])
                else:
                    winner, loser = b, a
                    winner_row, loser_row = row_b, row_a
                    winner_frac = float(metrics["overlap_fraction_of_b"])
                    loser_frac = float(metrics["overlap_fraction_of_a"])
                if loser_frac + 1e-15 < min_loser_cds_overlap_fraction:
                    continue
                winner_scaled = value_as_float(winner_row.get("model_selection_score_scaled"))
                loser_scaled = value_as_float(loser_row.get("model_selection_score_scaled"))
                score_margin = winner_scaled - loser_scaled
                pair_key = tuple(sorted((winner.output_transcript_id, loser.output_transcript_id)))
                candidate_pairs[pair_key] = {
                    "winner_model": winner,
                    "loser_model": loser,
                    "winner_row": winner_row,
                    "loser_row": loser_row,
                    "overlap_bp": overlap_bp,
                    "winner_overlap_fraction": winner_frac,
                    "loser_overlap_fraction": loser_frac,
                    "score_margin_scaled": score_margin,
                }
    return candidate_pairs


def final_prune_opposite_strand_shadows(
    prune_result: dict[str, object],
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    *,
    prune_sources: set[str],
    require_loser_single_exon: bool,
    min_overlap_bp: int,
    min_loser_cds_overlap_fraction: float,
    min_score_margin_scaled: float,
    pruning_level_name: str,
    progress: Optional[ProgressLogger] = None,
) -> dict[str, object]:
    """Remove lower-scoring opposite-strand shadow candidates from final selected models."""
    retained_models_in: list[CodingModel] = list(prune_result.get("retained_models", []))  # type: ignore[arg-type]
    retained_rows_in: list[dict[str, object]] = list(prune_result.get("retained_rows", []))  # type: ignore[arg-type]
    rows_by_id = {str(row.get("output_transcript_id")): row for row in retained_rows_in}
    info_by_id: dict[str, dict[str, object]] = dict(prune_result.get("retained_selection_info_by_transcript", {}))  # type: ignore[arg-type]
    removed_rows: list[dict[str, object]] = list(prune_result.get("removed_rows", []))  # type: ignore[arg-type]
    source_scope = ",".join(sorted(prune_sources)) if prune_sources else "ALL"

    if progress:
        progress.log(
            f"Final selected pruning: opposite-strand shadows; source_scope={source_scope}; "
            f"require_loser_single_exon={str(require_loser_single_exon).lower()}; models={len(retained_models_in):,}"
        )

    # Compare candidate discovery before and after earlier final-pruning filters.
    # This makes it clear when broad singleton pruning has already removed one
    # member of a potential opposite-strand shadow pair before the shadow pass.
    initial_models_for_shadow: list[CodingModel] = list(prune_result.get("initial_retained_models", retained_models_in))  # type: ignore[arg-type]
    initial_rows_in: list[dict[str, object]] = list(prune_result.get("initial_retained_rows", retained_rows_in))  # type: ignore[arg-type]
    initial_rows_by_id = {str(row.get("output_transcript_id")): row for row in initial_rows_in}

    if initial_models_for_shadow:
        candidate_pairs_before_singleton = find_opposite_strand_shadow_candidate_pairs(
            initial_models_for_shadow,
            initial_rows_by_id,
            model_structural_rows_by_transcript,
            min_overlap_bp=min_overlap_bp,
            min_loser_cds_overlap_fraction=min_loser_cds_overlap_fraction,
            progress=None,
        )
    else:
        candidate_pairs_before_singleton = {}

    candidate_pairs = find_opposite_strand_shadow_candidate_pairs(
        retained_models_in,
        rows_by_id,
        model_structural_rows_by_transcript,
        min_overlap_bp=min_overlap_bp,
        min_loser_cds_overlap_fraction=min_loser_cds_overlap_fraction,
        progress=progress,
        progress_label="Final selected pruning shadow scan",
    )

    shadow_candidate_pairs_before_singleton = len(candidate_pairs_before_singleton)
    shadow_candidate_pairs_after_singleton = len(candidate_pairs)
    shadow_candidate_pairs_skipped_after_singleton = max(
        0, shadow_candidate_pairs_before_singleton - shadow_candidate_pairs_after_singleton
    )

    ordered_pairs = sorted(
        candidate_pairs.values(),
        key=lambda item: (
            -value_as_float(item["winner_row"].get("model_selection_score_scaled")),  # type: ignore[index,union-attr]
            -value_as_float(item["winner_row"].get("model_selection_score")),  # type: ignore[index,union-attr]
            -float(item["loser_overlap_fraction"]),
            -int(item["overlap_bp"]),
            item["winner_model"].output_transcript_id,  # type: ignore[union-attr]
            item["loser_model"].output_transcript_id,  # type: ignore[union-attr]
        ),
    )

    removed_model_ids: set[str] = set()
    winner_removed_ids: dict[str, list[str]] = defaultdict(list)
    shadow_pair_rows: list[dict[str, object]] = []

    for pair_index, item in enumerate(ordered_pairs, start=1):
        winner: CodingModel = item["winner_model"]  # type: ignore[assignment]
        loser: CodingModel = item["loser_model"]  # type: ignore[assignment]
        winner_row: dict[str, object] = item["winner_row"]  # type: ignore[assignment]
        loser_row: dict[str, object] = item["loser_row"]  # type: ignore[assignment]
        reason = "."
        pruned = False
        if winner.output_transcript_id in removed_model_ids:
            reason = "winner_already_removed_by_higher_scoring_shadow_pair"
        elif loser.output_transcript_id in removed_model_ids:
            reason = "loser_already_removed_by_higher_scoring_shadow_pair"
        elif prune_sources and loser.source_name not in prune_sources:
            reason = "loser_source_not_eligible"
        elif require_loser_single_exon and loser.cds_exon_count != 1:
            reason = "loser_not_single_exon"
        elif float(item["score_margin_scaled"]) + 1e-15 < min_score_margin_scaled:
            reason = "score_margin_below_threshold"
        else:
            pruned = True
            reason = "opposite_strand_shadow_lower_score"
            removed_model_ids.add(loser.output_transcript_id)
            winner_removed_ids[winner.output_transcript_id].append(loser.output_transcript_id)

        shadow_pair_rows.append({
            "shadow_pair_id": f"shadow_pair_{pair_index:06d}",
            "shadow_pruned": "true" if pruned else "false",
            "prune_reason": reason,
            "chrom": winner.chrom,
            "winner_model_id": winner.output_transcript_id,
            "winner_source": winner.source_name,
            "winner_strand": winner.strand,
            "winner_score_scaled": value_as_float(winner_row.get("model_selection_score_scaled")),
            "winner_score_raw": value_as_float(winner_row.get("model_selection_score")),
            "winner_cds_length": winner.cds_length,
            "winner_cds_exon_count": winner.cds_exon_count,
            "loser_model_id": loser.output_transcript_id,
            "loser_source": loser.source_name,
            "loser_strand": loser.strand,
            "loser_score_scaled": value_as_float(loser_row.get("model_selection_score_scaled")),
            "loser_score_raw": value_as_float(loser_row.get("model_selection_score")),
            "loser_cds_length": loser.cds_length,
            "loser_cds_exon_count": loser.cds_exon_count,
            "overlap_bp": item["overlap_bp"],
            "winner_overlap_fraction": item["winner_overlap_fraction"],
            "loser_overlap_fraction": item["loser_overlap_fraction"],
            "score_margin_scaled": item["score_margin_scaled"],
            "min_overlap_bp": min_overlap_bp,
            "min_loser_cds_overlap_fraction": min_loser_cds_overlap_fraction,
            "min_score_margin_scaled": min_score_margin_scaled,
            "require_loser_single_exon": "true" if require_loser_single_exon else "false",
            "source_scope": source_scope,
        })

    retained_models: list[CodingModel] = []
    retained_rows: list[dict[str, object]] = []
    retained_info_by_id: dict[str, dict[str, object]] = {}
    removed_rows_out: list[dict[str, object]] = list(removed_rows)

    best_pruning_pair_by_loser = {
        str(row["loser_model_id"]): row for row in shadow_pair_rows if str(row.get("shadow_pruned")) == "true"
    }

    for model in retained_models_in:
        model_id = model.output_transcript_id
        row = dict(rows_by_id.get(model_id, {}))
        info = dict(info_by_id.get(model_id, row))
        shadow_common = {
            "final_prune_level": pruning_level_name,
            "final_shadow_prune_enabled": "true",
            "final_shadow_prune_source_scope": source_scope,
            "final_shadow_prune_require_loser_single_exon": "true" if require_loser_single_exon else "false",
        }
        if model_id in removed_model_ids:
            pair_row = best_pruning_pair_by_loser.get(model_id, {})
            shadow_fields = {
                **shadow_common,
                "final_shadow_pruned": "true",
                "final_shadow_prune_reason": pair_row.get("prune_reason", "opposite_strand_shadow_lower_score"),
                "final_shadow_prune_retained_model_id": pair_row.get("winner_model_id", "."),
                "final_shadow_prune_retained_model_source": pair_row.get("winner_source", "."),
                "final_shadow_prune_overlap_bp": pair_row.get("overlap_bp", "."),
                "final_shadow_prune_retained_model_overlap_fraction": pair_row.get("winner_overlap_fraction", "."),
                "final_shadow_prune_removed_model_overlap_fraction": pair_row.get("loser_overlap_fraction", "."),
                "final_shadow_prune_score_margin_scaled": pair_row.get("score_margin_scaled", "."),
            }
            out_row = {**row, **shadow_fields}
            append_prune_class(out_row, "opposite_strand_shadow")
            removed_rows_out.append(out_row)
        else:
            removed_by_this = sorted(winner_removed_ids.get(model_id, []))
            shadow_fields = {
                **shadow_common,
                "final_shadow_pruned": "false",
                "final_shadow_prune_reason": ".",
                "final_shadow_prune_removed_model_count": len(removed_by_this),
                "final_shadow_prune_removed_model_ids": ",".join(removed_by_this) if removed_by_this else ".",
            }
            out_row = {**row, **shadow_fields}
            retained_models.append(model)
            retained_rows.append(out_row)
            retained_info_by_id[model_id] = {**info, **shadow_fields, "final_prune_classes": out_row.get("final_prune_classes", ".")}

    if progress:
        progress.log(
            f"Final selected shadow pruning complete: retained={len(retained_models):,}; "
            f"newly_removed={len(removed_model_ids):,}; candidate_pairs={len(shadow_pair_rows):,}"
        )

    # Preserve summary metadata from earlier final-pruning passes.  This is
    # important when shadow pruning runs after singleton/protein/low-score
    # pruning: row-level removals were already carried forward in removed_rows,
    # but the stdout summary also needs the earlier pass counters and adequacy
    # fields.  State fields modified by this pass are overwritten below.
    out_result = dict(prune_result)
    out_result.update({
        "initial_retained_models": list(prune_result.get("initial_retained_models", [])),
        "initial_retained_rows": [dict(row) for row in prune_result.get("initial_retained_rows", [])],
        "retained_models": retained_models,
        "retained_rows": retained_rows,
        "removed_rows": removed_rows_out,
        "retained_selection_info_by_transcript": retained_info_by_id,
        "shadow_pair_rows": shadow_pair_rows,
        "shadow_candidate_pairs_before_singleton": shadow_candidate_pairs_before_singleton,
        "shadow_candidate_pairs_after_singleton": shadow_candidate_pairs_after_singleton,
        "shadow_candidate_pairs_skipped_after_singleton": shadow_candidate_pairs_skipped_after_singleton,
    })
    return out_result


def write_final_pruned_outputs(
    prune_result: dict[str, object],
    *,
    final_pruned_selected_models_tsv: Optional[Path],
    final_pruned_selected_removed_tsv: Optional[Path],
    final_shadow_gene_pairs_tsv: Optional[Path] = None,
) -> None:
    columns = selection_stage_output_columns() + FINAL_PRUNE_EXTRA_COLUMNS
    if final_pruned_selected_models_tsv is not None:
        write_dict_tsv(prune_result["retained_rows"], final_pruned_selected_models_tsv, columns)  # type: ignore[arg-type]
    if final_pruned_selected_removed_tsv is not None:
        write_dict_tsv(prune_result["removed_rows"], final_pruned_selected_removed_tsv, columns)  # type: ignore[arg-type]
    if final_shadow_gene_pairs_tsv is not None:
        write_dict_tsv(prune_result.get("shadow_pair_rows", []), final_shadow_gene_pairs_tsv, SHADOW_PAIR_COLUMNS)  # type: ignore[arg-type]
