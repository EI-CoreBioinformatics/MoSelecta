"""Candidate-level structural-coherence penalty and region audit.

The component asks whether one candidate joins two well-supported CDS blocks
while independent annotation sources support those blocks in separate models.
Evidence is collapsed by annotation source. Human-readable anchor states are
retained for audit, while factual continuity/gene-separation evidence supplies
the independent S and C votes.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re
from typing import Iterable

from .models import CodingModel, Superlocus
from .tsvio import write_dict_tsv


Part = tuple[int, int]
AnchorCoverageCache = dict[tuple[int, Part], float]

STRUCTURAL_COHERENCE_ANCHOR_COVERAGE = 0.50
STRUCTURAL_COHERENCE_MIN_SPLIT_SOURCES = 2
STRUCTURAL_COHERENCE_MIN_ANCHOR_SOURCES = 2
STRUCTURAL_COHERENCE_BREADTH_CAP = 1.50


STRUCTURAL_COHERENCE_REGION_COLUMNS = [
    "superlocus_id", "candidate_model_id", "candidate_source", "chrom", "strand",
    "region_id", "region_start_0based", "region_end_0based", "region_length_bp",
    "left_anchor_start_0based", "left_anchor_end_0based",
    "right_anchor_start_0based", "right_anchor_end_0based",
    "first_cut_index", "last_cut_index", "number_of_candidate_boundaries_in_region",
    "connector_intron_count", "connector_internal_cds_exon_count", "boundaries_merged",
    "anchor_candidate_base_coverage_threshold",
    "bilateral_separated_source_count", "independent_spanning_source_count", "mixed_source_count",
    "bilateral_separated_sources", "independent_spanning_sources", "mixed_sources",
    "left_anchor_independent_source_count", "right_anchor_independent_source_count",
    "left_anchor_independent_sources", "right_anchor_independent_sources",
    "left_only_source_count", "right_only_source_count", "uninformative_source_count",
    "source_anchor_states",
    "bridge_hypothesis_group_id", "bridge_hypothesis_member_candidate_ids",
    "source_has_exact_chain", "source_has_exact_connector_continuity",
    "source_has_anchor_crossing", "source_has_bridge_group_continuity",
    "source_bridge_group_continuity_direct", "source_bridge_group_continuity_propagated",
    "source_has_left_support", "source_has_right_support",
    "source_left_gene_ids", "source_right_gene_ids", "source_gene_id_resolution_method",
    "source_same_gene_alternative", "source_gene_relationship", "source_state_label",
    "source_contributes_S", "source_contributes_C", "source_C_reason", "source_S_reason",
    "source_evidence_model_ids",
    "connector_intron_independent_exact_support_counts",
    "connector_donor_independent_exact_support_counts",
    "connector_acceptor_independent_exact_support_counts",
    "connector_internal_cds_exon_independent_exact_support_counts",
    "connector_internal_cds_exon_independent_anchor_support_counts",
    "candidate_only_connector_internal_cds_exon_count",
    "penalty_gate_passed", "penalty_gate_reason",
    "direction", "breadth_factor", "structural_coherence_raw_score",
    "structural_coherence_penalty_weight", "structural_coherence_penalty_score",
]


def _merged_parts(parts: Iterable[Part]) -> list[Part]:
    merged: list[list[int]] = []
    for start0, end0 in sorted(parts):
        if not merged or start0 > merged[-1][1]:
            merged.append([start0, end0])
        else:
            merged[-1][1] = max(merged[-1][1], end0)
    return [(start0, end0) for start0, end0 in merged]


def _overlap_bp(interval: Part, parts: Iterable[Part]) -> int:
    return sum(
        max(0, min(interval[1], end0) - max(interval[0], start0))
        for start0, end0 in _merged_parts(parts)
    )


def _anchor_coverage(
    model: CodingModel, anchor: Part, coverage_cache: AnchorCoverageCache | None = None,
) -> float:
    # A model/anchor pair is often evaluated for many candidate cuts and
    # regions. The cache lives only for one superlocus (see annotation loop).
    key = (id(model), anchor)
    if coverage_cache is not None and key in coverage_cache:
        return coverage_cache[key]
    length = anchor[1] - anchor[0]
    coverage = _overlap_bp(anchor, model.cds_parts) / length if length > 0 else 0.0
    if coverage_cache is not None:
        coverage_cache[key] = coverage
    return coverage


def _source_anchor_state(
    models: list[CodingModel],
    left_anchor: Part,
    right_anchor: Part,
    threshold: float,
    coverage_cache: AnchorCoverageCache | None = None,
) -> str:
    crossing = False
    left_local = False
    right_local = False
    for model in models:
        left = _anchor_coverage(model, left_anchor, coverage_cache) + 1e-15 >= threshold
        right = _anchor_coverage(model, right_anchor, coverage_cache) + 1e-15 >= threshold
        crossing = crossing or (left and right)
        left_local = left_local or (left and not right)
        right_local = right_local or (right and not left)
    if crossing and left_local and right_local:
        return "mixed"
    if crossing:
        return "crossing"
    if left_local and right_local:
        return "bilateral_separated"
    if left_local:
        return "left_only"
    if right_local:
        return "right_only"
    return "uninformative"


def _states_for_anchors(
    models_by_source: dict[str, list[CodingModel]],
    left_anchor: Part,
    right_anchor: Part,
    threshold: float,
    coverage_cache: AnchorCoverageCache | None = None,
) -> dict[str, str]:
    return {
        source: _source_anchor_state(source_models, left_anchor, right_anchor, threshold, coverage_cache)
        for source, source_models in sorted(models_by_source.items())
    }


def _source_votes(states: dict[str, str], candidate_source: str) -> dict[str, object]:
    crossing = {source for source, state in states.items() if state == "crossing"}
    separated = {source for source, state in states.items() if state == "bilateral_separated"}
    mixed = {source for source, state in states.items() if state == "mixed"}
    independent_spanning = crossing - {candidate_source}
    return {
        "bilateral_separated_source_count": len(separated),
        "independent_spanning_source_count": len(independent_spanning),
        "mixed_source_count": len(mixed),
        "bilateral_separated_sources": _format_values(separated),
        "independent_spanning_sources": _format_values(independent_spanning),
        "mixed_sources": _format_values(mixed),
        "left_only_source_count": sum(state == "left_only" for state in states.values()),
        "right_only_source_count": sum(state == "right_only" for state in states.values()),
        "uninformative_source_count": sum(state == "uninformative" for state in states.values()),
    }


def _format_values(values: Iterable[str]) -> str:
    ordered = sorted(values)
    return ",".join(ordered) if ordered else "."


def _sources_supporting_anchor(
    models_by_source: dict[str, list[CodingModel]],
    anchor: Part,
    threshold: float,
    coverage_cache: AnchorCoverageCache | None = None,
) -> set[str]:
    return {
        source
        for source, source_models in models_by_source.items()
        if any(_anchor_coverage(model, anchor, coverage_cache) + 1e-15 >= threshold for model in source_models)
    }


def _sources_with_exact_intron(models_by_source: dict[str, list[CodingModel]], intron: Part) -> set[str]:
    return {
        source
        for source, source_models in models_by_source.items()
        if any(intron in model.cds_introns for model in source_models)
    }


def _sources_with_splice_site(models_by_source: dict[str, list[CodingModel]], site0: int) -> set[str]:
    return {
        source
        for source, source_models in models_by_source.items()
        if any(site0 in model.cds_splice_sites for model in source_models)
    }


def _sources_with_exact_exon(models_by_source: dict[str, list[CodingModel]], exon: Part) -> set[str]:
    return {
        source
        for source, source_models in models_by_source.items()
        if any(exon in model.cds_parts for model in source_models)
    }


def _cut_is_region_seed(
    *,
    separated_count: int,
    independent_spanning_count: int,
    mixed_count: int,
    left_independent_count: int,
    right_independent_count: int,
    independent_exact_intron_count: int,
) -> bool:
    """Retain the accepted v2 region-discovery behavior without scoring labels."""
    if separated_count > independent_spanning_count:
        return True
    if separated_count == independent_spanning_count and separated_count > 0:
        return True
    if independent_spanning_count > 0 or mixed_count > 0:
        return False
    candidate_specific_boundary = (
        independent_exact_intron_count == 0
        and (left_independent_count == 0 or right_independent_count == 0)
    )
    return candidate_specific_boundary or min(left_independent_count, right_independent_count) > 0


def _counts_text(values: list[int]) -> str:
    return ",".join(str(value) for value in values) if values else "."


_TRANSCRIPT_SUFFIX_RE = re.compile(
    r"(?:[._-](?:t|transcript|isoform|mrna)[._-]?[0-9A-Za-z]+)$", re.IGNORECASE
)


def _gene_identity(model: CodingModel) -> tuple[str | None, str]:
    """Use source-gene metadata, never merely unequal transcript identifiers."""
    if model.input_gene_id not in {"", ".", "NA"}:
        return model.input_gene_id, "input_gene_id"
    if model.gene_record is not None:
        for key in ("original_gene_id", "gene_id", "ID"):
            value = model.gene_record.attrs.get(key)
            if value and value not in {".", "NA"}:
                return value, "gene_record"
    if model.transcript_record is not None:
        for key in ("original_gene_id", "gene_id", "Parent"):
            value = model.transcript_record.attrs.get(key)
            if value and value not in {".", "NA"} and "," not in value:
                return value, "transcript_record"
    transcript = model.input_transcript_id
    if transcript not in {"", ".", "NA"}:
        root = _TRANSCRIPT_SUFFIX_RE.sub("", transcript)
        if root and root != transcript:
            return root, "normalized_transcript_id"
    return None, "unresolved"


def _source_gene_relationship(left_models: list[CodingModel], right_models: list[CodingModel]) -> dict[str, object]:
    left = [(model, *_gene_identity(model)) for model in left_models]
    right = [(model, *_gene_identity(model)) for model in right_models]
    shared = False
    distinct = False
    missing = False
    for _left_model, left_key, _left_method in left:
        for _right_model, right_key, _right_method in right:
            if left_key is None or right_key is None:
                missing = True
            elif left_key == right_key:
                shared = True
            else:
                distinct = True
    if shared:
        relationship = "same_gene_alternative_isoforms"
    elif distinct and not missing:
        relationship = "distinct_source_genes"
    elif left_models and right_models:
        relationship = "ambiguous_missing_gene_metadata"
    else:
        relationship = "not_bilateral"
    methods = {method for _model, _key, method in left + right}
    resolution = (
        "unresolved" if "unresolved" in methods else
        "normalized_transcript_id" if "normalized_transcript_id" in methods else
        "source_record" if methods & {"gene_record", "transcript_record"} else
        "input_gene_id" if methods else "."
    )
    return {
        "relationship": relationship,
        "permits_s": relationship == "distinct_source_genes",
        "same_gene": shared,
        "left_gene_ids": {key for _model, key, _method in left if key is not None},
        "right_gene_ids": {key for _model, key, _method in right if key is not None},
        "resolution": resolution,
    }


def _contiguous_subsequence(needle: tuple[Part, ...], haystack: tuple[Part, ...]) -> bool:
    if not needle:
        return False
    width = len(needle)
    return any(haystack[index:index + width] == needle for index in range(len(haystack) - width + 1))


def _source_region_facts(
    region: dict[str, object], candidate: CodingModel, source_models: list[CodingModel], threshold: float,
    coverage_cache: AnchorCoverageCache | None = None,
) -> dict[str, object]:
    left_anchor = (int(region["left_anchor_start_0based"]), int(region["left_anchor_end_0based"]))
    right_anchor = (int(region["right_anchor_start_0based"]), int(region["right_anchor_end_0based"]))
    left_models: list[CodingModel] = []
    right_models: list[CodingModel] = []
    crossing_models: list[CodingModel] = []
    for model in source_models:
        left = _anchor_coverage(model, left_anchor, coverage_cache) + 1e-15 >= threshold
        right = _anchor_coverage(model, right_anchor, coverage_cache) + 1e-15 >= threshold
        if left and right:
            crossing_models.append(model)
        elif left:
            left_models.append(model)
        elif right:
            right_models.append(model)
    state = _source_anchor_state(source_models, left_anchor, right_anchor, threshold, coverage_cache)
    exact_chain = [model for model in source_models if model.cds_parts == candidate.cds_parts]
    first_cut = int(region["first_cut_index"])
    last_cut = int(region["last_cut_index"])
    connector_introns = candidate.cds_introns[first_cut - 1:last_cut]
    if len(connector_introns) == 1:
        exact_connector = [model for model in source_models if connector_introns[0] in model.cds_introns]
    elif len(connector_introns) > 1:
        exact_connector = [model for model in source_models if _contiguous_subsequence(connector_introns, model.cds_introns)]
    else:
        exact_connector = []
    relationship = _source_gene_relationship(left_models, right_models)
    direct_routes: list[str] = []
    if exact_chain:
        direct_routes.append("exact_complete_chain")
    if exact_connector and not exact_chain:
        direct_routes.append("exact_connector")
    if state == "crossing" and not exact_chain and not exact_connector:
        direct_routes.append("anchor_based_non_exact_continuity")
    if state == "mixed":
        direct_routes.append("mixed_source_spanning_interpretation")
    direct_continuity = bool(exact_chain or exact_connector or state in {"crossing", "mixed"})
    return {
        "state": state,
        "left_model_ids": {model.output_transcript_id for model in left_models},
        "right_model_ids": {model.output_transcript_id for model in right_models},
        "crossing_model_ids": {model.output_transcript_id for model in crossing_models},
        "exact_chain_model_ids": {model.output_transcript_id for model in exact_chain},
        "exact_connector_model_ids": {model.output_transcript_id for model in exact_connector},
        "has_exact_chain": bool(exact_chain),
        "has_exact_connector": bool(exact_connector),
        "has_anchor_crossing": bool(crossing_models),
        "has_left_support": bool(left_models),
        "has_right_support": bool(right_models),
        "relationship": relationship["relationship"],
        "same_gene": relationship["same_gene"],
        "left_gene_ids": relationship["left_gene_ids"],
        "right_gene_ids": relationship["right_gene_ids"],
        "resolution": relationship["resolution"],
        "direct_routes": direct_routes,
        "direct_continuity": direct_continuity,
        "preliminary_s": bool(relationship["permits_s"] and not direct_continuity),
        "bridge_continuity": False,
        "bridge_direct": False,
        "bridge_propagated": False,
    }


def _shorter_overlap_fraction(left: Part, right: Part) -> float:
    overlap = max(0, min(left[1], right[1]) - max(left[0], right[0]))
    denominator = min(left[1] - left[0], right[1] - right[0])
    return overlap / denominator if denominator > 0 else 0.0


def _bridge_compatible(
    left: dict[str, object], right: dict[str, object],
    candidates: dict[str, CodingModel], preliminary_s: dict[str, set[str]],
) -> bool:
    left_candidate = candidates[str(left["candidate_model_id"])]
    right_candidate = candidates[str(right["candidate_model_id"])]
    if (
        left_candidate.output_transcript_id == right_candidate.output_transcript_id
        or left_candidate.sublocus_id is None
        or right_candidate.sublocus_id is None
        or left_candidate.sublocus_id != right_candidate.sublocus_id
        or left_candidate.strand != right_candidate.strand
    ):
        return False
    left_connector = (int(left["region_start_0based"]), int(left["region_end_0based"]))
    right_connector = (int(right["region_start_0based"]), int(right["region_end_0based"]))
    nested = (
        left_connector[0] <= right_connector[0] and left_connector[1] >= right_connector[1]
        or right_connector[0] <= left_connector[0] and right_connector[1] >= left_connector[1]
    )
    shared_boundary = bool(set(left_connector) & set(right_connector))
    if _shorter_overlap_fraction(left_connector, right_connector) <= 0 or not (nested or shared_boundary):
        return False
    for side in ("left", "right"):
        left_anchor = (int(left[f"{side}_anchor_start_0based"]), int(left[f"{side}_anchor_end_0based"]))
        right_anchor = (int(right[f"{side}_anchor_start_0based"]), int(right[f"{side}_anchor_end_0based"]))
        if _shorter_overlap_fraction(left_anchor, right_anchor) < 0.50:
            return False
    left_s = preliminary_s[str(left["region_id"])]
    right_s = preliminary_s[str(right["region_id"])]
    denominator = min(len(left_s), len(right_s))
    return denominator > 0 and len(left_s & right_s) >= 2 and len(left_s & right_s) / denominator >= 2 / 3


def _assign_bridge_groups(
    regions: list[dict[str, object]], facts_by_region: dict[str, dict[str, dict[str, object]]],
    candidates: dict[str, CodingModel],
) -> None:
    preliminary_s = {
        region_id: {source for source, facts in source_facts.items() if facts["preliminary_s"]}
        for region_id, source_facts in facts_by_region.items()
    }
    by_sublocus: dict[str, list[dict[str, object]]] = defaultdict(list)
    for region in regions:
        candidate = candidates[str(region["candidate_model_id"])]
        by_sublocus[str(candidate.sublocus_id)].append(region)
        region["bridge_hypothesis_group_id"] = "."
        region["bridge_hypothesis_member_candidate_ids"] = "."
    for sublocus_id, local_regions in sorted(by_sublocus.items()):
        parent = list(range(len(local_regions)))

        def find(index: int) -> int:
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = parent[index]
            return index

        for index, left in enumerate(local_regions):
            for other in range(index + 1, len(local_regions)):
                if _bridge_compatible(left, local_regions[other], candidates, preliminary_s):
                    left_root, right_root = find(index), find(other)
                    if left_root != right_root:
                        parent[right_root] = left_root
        components: dict[int, list[dict[str, object]]] = defaultdict(list)
        for index, region in enumerate(local_regions):
            components[find(index)].append(region)
        group_number = 0
        for component in components.values():
            member_ids = {str(region["candidate_model_id"]) for region in component}
            if len(member_ids) < 2:
                continue
            group_number += 1
            group_id = f"{sublocus_id}.bridge_hypothesis{group_number:02d}"
            continuity_sources = {
                source for region in component
                for source, facts in facts_by_region[str(region["region_id"])].items()
                if facts["direct_continuity"]
            }
            for region in component:
                region_id = str(region["region_id"])
                region["bridge_hypothesis_group_id"] = group_id
                region["bridge_hypothesis_member_candidate_ids"] = _format_values(member_ids)
                for source, facts in facts_by_region[region_id].items():
                    if source not in continuity_sources:
                        continue
                    direct = bool(facts["direct_continuity"]) or any(
                        other is not region
                        and facts_by_region[str(other["region_id"])][source]["direct_continuity"]
                        and _bridge_compatible(region, other, candidates, preliminary_s)
                        for other in component
                    )
                    facts["bridge_continuity"] = True
                    facts["bridge_direct"] = direct
                    facts["bridge_propagated"] = not direct


def _serialize_source_values(facts: dict[str, dict[str, object]], key: str) -> str:
    values = []
    for source, source_facts in sorted(facts.items()):
        value = source_facts[key]
        if isinstance(value, bool):
            rendered = "1" if value else "0"
        elif isinstance(value, (set, list, tuple)):
            rendered = _format_values(str(item) for item in value)
        else:
            rendered = str(value)
        values.append(f"{source}:{rendered}")
    return ";".join(values) if values else "."


def _score_with_revised_votes(region: dict[str, object], facts: dict[str, dict[str, object]]) -> None:
    candidate_source = str(region["candidate_source"])
    separated: set[str] = set()
    spanning: set[str] = set()
    mixed = {source for source, item in facts.items() if item["state"] == "mixed"}
    for source, item in facts.items():
        continuity = bool(item["direct_continuity"] or item["bridge_continuity"])
        item["contributes_C"] = continuity and source != candidate_source
        item["contributes_S"] = bool(item["preliminary_s"] and not continuity)
        routes = list(item["direct_routes"])
        if item["bridge_continuity"]:
            routes.append("bridge_hypothesis_grouping")
        item["C_reason"] = routes[0] if item["contributes_C"] and routes else "."
        item["S_reason"] = "distinct_source_genes" if item["contributes_S"] else "."
        if item["contributes_C"]:
            spanning.add(source)
        if item["contributes_S"]:
            separated.add(source)
    region.update({
        "bilateral_separated_source_count": len(separated),
        "independent_spanning_source_count": len(spanning),
        "mixed_source_count": len(mixed),
        "bilateral_separated_sources": _format_values(separated),
        "independent_spanning_sources": _format_values(spanning),
        "mixed_sources": _format_values(mixed),
        "source_anchor_states": ";".join(f"{source}:{item['state']}" for source, item in sorted(facts.items())),
        "source_has_exact_chain": _serialize_source_values(facts, "has_exact_chain"),
        "source_has_exact_connector_continuity": _serialize_source_values(facts, "has_exact_connector"),
        "source_has_anchor_crossing": _serialize_source_values(facts, "has_anchor_crossing"),
        "source_has_bridge_group_continuity": _serialize_source_values(facts, "bridge_continuity"),
        "source_bridge_group_continuity_direct": _serialize_source_values(facts, "bridge_direct"),
        "source_bridge_group_continuity_propagated": _serialize_source_values(facts, "bridge_propagated"),
        "source_has_left_support": _serialize_source_values(facts, "has_left_support"),
        "source_has_right_support": _serialize_source_values(facts, "has_right_support"),
        "source_left_gene_ids": _serialize_source_values(facts, "left_gene_ids"),
        "source_right_gene_ids": _serialize_source_values(facts, "right_gene_ids"),
        "source_gene_id_resolution_method": _serialize_source_values(facts, "resolution"),
        "source_same_gene_alternative": _serialize_source_values(facts, "same_gene"),
        "source_gene_relationship": _serialize_source_values(facts, "relationship"),
        "source_state_label": _serialize_source_values(facts, "state"),
        "source_contributes_S": _serialize_source_values(facts, "contributes_S"),
        "source_contributes_C": _serialize_source_values(facts, "contributes_C"),
        "source_C_reason": _serialize_source_values(facts, "C_reason"),
        "source_S_reason": _serialize_source_values(facts, "S_reason"),
        "source_evidence_model_ids": ";".join(
            f"{source}:left={_format_values(item['left_model_ids'])}|right={_format_values(item['right_model_ids'])}"
            f"|crossing={_format_values(item['crossing_model_ids'])}"
            f"|chain={_format_values(item['exact_chain_model_ids'])}"
            f"|connector={_format_values(item['exact_connector_model_ids'])}"
            for source, item in sorted(facts.items())
        ) or ".",
    })
    left_count = int(region["left_anchor_independent_source_count"])
    right_count = int(region["right_anchor_independent_source_count"])
    gate_parts = {
        "S>=2": len(separated) >= STRUCTURAL_COHERENCE_MIN_SPLIT_SOURCES,
        "S>C": len(separated) > len(spanning),
        "left_anchor>=2": left_count >= STRUCTURAL_COHERENCE_MIN_ANCHOR_SOURCES,
        "right_anchor>=2": right_count >= STRUCTURAL_COHERENCE_MIN_ANCHOR_SOURCES,
    }
    gate_passed = all(gate_parts.values())
    denominator = len(separated) + len(spanning)
    direction = (len(separated) - len(spanning)) / denominator if denominator else 0.0
    breadth = min(
        STRUCTURAL_COHERENCE_BREADTH_CAP,
        1.0 + 0.10 * (len(separated) - STRUCTURAL_COHERENCE_MIN_SPLIT_SOURCES),
    )
    raw = direction * breadth if gate_passed else 0.0
    region.update({
        "penalty_gate_passed": "true" if gate_passed else "false",
        "penalty_gate_reason": "eligible" if gate_passed else ",".join(name for name, passed in gate_parts.items() if not passed),
        "direction": direction,
        "breadth_factor": breadth,
        "structural_coherence_raw_score": raw,
        "structural_coherence_penalty_score": float(region["structural_coherence_penalty_weight"]) * raw,
    })


def _candidate_regions(
    superlocus: Superlocus,
    candidate: CodingModel,
    models_by_source: dict[str, list[CodingModel]],
    *,
    anchor_threshold: float,
    penalty_weight: float,
    coverage_cache: AnchorCoverageCache | None = None,
) -> list[dict[str, object]]:
    parts = candidate.cds_parts
    if len(parts) < 2:
        return []

    seed_indices: list[int] = []
    for cut_index, (left_anchor, right_anchor) in enumerate(zip(parts, parts[1:]), start=1):
        states = _states_for_anchors(models_by_source, left_anchor, right_anchor, anchor_threshold, coverage_cache)
        votes = _source_votes(states, candidate.source_name)
        left_independent = _sources_supporting_anchor(models_by_source, left_anchor, anchor_threshold, coverage_cache) - {candidate.source_name}
        right_independent = _sources_supporting_anchor(models_by_source, right_anchor, anchor_threshold, coverage_cache) - {candidate.source_name}
        intron = (left_anchor[1], right_anchor[0])
        exact_independent = _sources_with_exact_intron(models_by_source, intron) - {candidate.source_name}
        if _cut_is_region_seed(
            separated_count=int(votes["bilateral_separated_source_count"]),
            independent_spanning_count=int(votes["independent_spanning_source_count"]),
            mixed_count=int(votes["mixed_source_count"]),
            left_independent_count=len(left_independent),
            right_independent_count=len(right_independent),
            independent_exact_intron_count=len(exact_independent),
        ):
            seed_indices.append(cut_index)

    groups: list[list[int]] = []
    for cut_index in seed_indices:
        if not groups:
            groups.append([cut_index])
            continue
        previous = groups[-1][-1]
        intervening_exon = parts[cut_index - 1]
        independent_support = (
            _sources_supporting_anchor(models_by_source, intervening_exon, anchor_threshold, coverage_cache)
            - {candidate.source_name}
        )
        if cut_index == previous + 1 and not independent_support:
            groups[-1].append(cut_index)
        else:
            groups.append([cut_index])

    rows: list[dict[str, object]] = []
    for region_index, group in enumerate(groups, start=1):
        first_cut = group[0]
        last_cut = group[-1]
        left_anchor = parts[first_cut - 1]
        right_anchor = parts[last_cut]
        states = _states_for_anchors(models_by_source, left_anchor, right_anchor, anchor_threshold, coverage_cache)
        votes = _source_votes(states, candidate.source_name)
        left_independent = _sources_supporting_anchor(models_by_source, left_anchor, anchor_threshold, coverage_cache) - {candidate.source_name}
        right_independent = _sources_supporting_anchor(models_by_source, right_anchor, anchor_threshold, coverage_cache) - {candidate.source_name}
        connector_introns = [
            (parts[index - 1][1], parts[index][0])
            for index in range(first_cut, last_cut + 1)
        ]
        connector_exons = list(parts[first_cut:last_cut])
        intron_support = [
            len(_sources_with_exact_intron(models_by_source, intron) - {candidate.source_name})
            for intron in connector_introns
        ]
        donor_support = [
            len(_sources_with_splice_site(models_by_source, intron[0]) - {candidate.source_name})
            for intron in connector_introns
        ]
        acceptor_support = [
            len(_sources_with_splice_site(models_by_source, intron[1]) - {candidate.source_name})
            for intron in connector_introns
        ]
        exon_exact_support = [
            len(_sources_with_exact_exon(models_by_source, exon) - {candidate.source_name})
            for exon in connector_exons
        ]
        exon_anchor_support = [
            len(_sources_supporting_anchor(models_by_source, exon, anchor_threshold, coverage_cache) - {candidate.source_name})
            for exon in connector_exons
        ]

        separated = int(votes["bilateral_separated_source_count"])
        spanning = int(votes["independent_spanning_source_count"])
        gate_parts = {
            "S>=2": separated >= STRUCTURAL_COHERENCE_MIN_SPLIT_SOURCES,
            "S>C": separated > spanning,
            "left_anchor>=2": len(left_independent) >= STRUCTURAL_COHERENCE_MIN_ANCHOR_SOURCES,
            "right_anchor>=2": len(right_independent) >= STRUCTURAL_COHERENCE_MIN_ANCHOR_SOURCES,
        }
        gate_passed = all(gate_parts.values())
        denominator = separated + spanning
        direction = (separated - spanning) / denominator if denominator else 0.0
        breadth_factor = min(
            STRUCTURAL_COHERENCE_BREADTH_CAP,
            1.0 + 0.10 * (separated - STRUCTURAL_COHERENCE_MIN_SPLIT_SOURCES),
        )
        raw_score = direction * breadth_factor if gate_passed else 0.0
        region_id = f"{candidate.output_transcript_id}.coherence{region_index:02d}"
        rows.append({
            "superlocus_id": superlocus.superlocus_id,
            "candidate_model_id": candidate.output_transcript_id,
            "candidate_source": candidate.source_name,
            "chrom": candidate.chrom,
            "strand": candidate.strand,
            "region_id": region_id,
            "region_start_0based": left_anchor[1],
            "region_end_0based": right_anchor[0],
            "region_length_bp": right_anchor[0] - left_anchor[1],
            "left_anchor_start_0based": left_anchor[0],
            "left_anchor_end_0based": left_anchor[1],
            "right_anchor_start_0based": right_anchor[0],
            "right_anchor_end_0based": right_anchor[1],
            "first_cut_index": first_cut,
            "last_cut_index": last_cut,
            "number_of_candidate_boundaries_in_region": len(group),
            "connector_intron_count": len(connector_introns),
            "connector_internal_cds_exon_count": len(connector_exons),
            "boundaries_merged": "true" if len(group) > 1 else "false",
            "anchor_candidate_base_coverage_threshold": anchor_threshold,
            **votes,
            "left_anchor_independent_source_count": len(left_independent),
            "right_anchor_independent_source_count": len(right_independent),
            "left_anchor_independent_sources": _format_values(left_independent),
            "right_anchor_independent_sources": _format_values(right_independent),
            "source_anchor_states": ";".join(f"{source}:{state}" for source, state in sorted(states.items())),
            "connector_intron_independent_exact_support_counts": _counts_text(intron_support),
            "connector_donor_independent_exact_support_counts": _counts_text(donor_support),
            "connector_acceptor_independent_exact_support_counts": _counts_text(acceptor_support),
            "connector_internal_cds_exon_independent_exact_support_counts": _counts_text(exon_exact_support),
            "connector_internal_cds_exon_independent_anchor_support_counts": _counts_text(exon_anchor_support),
            "candidate_only_connector_internal_cds_exon_count": sum(count == 0 for count in exon_anchor_support),
            "penalty_gate_passed": "true" if gate_passed else "false",
            "penalty_gate_reason": ",".join(
                name for name, passed in gate_parts.items() if not passed
            ) if not gate_passed else "eligible",
            "direction": direction,
            "breadth_factor": breadth_factor,
            "structural_coherence_raw_score": raw_score,
            "structural_coherence_penalty_weight": penalty_weight,
            "structural_coherence_penalty_score": penalty_weight * raw_score,
        })
    return rows


def annotate_structural_coherence(
    superloci: list[Superlocus],
    model_rows: list[dict[str, object]],
    *,
    penalty_weight: float,
    anchor_threshold: float = STRUCTURAL_COHERENCE_ANCHOR_COVERAGE,
) -> list[dict[str, object]]:
    """Attach candidate summaries and return region-level audit rows."""
    if penalty_weight < 0:
        raise ValueError("structural-coherence penalty weight must be >= 0")
    if not 0 < anchor_threshold <= 1:
        raise ValueError("structural-coherence anchor threshold must be > 0 and <= 1")

    row_by_model = {str(row.get("output_transcript_id")): row for row in model_rows}
    region_rows: list[dict[str, object]] = []
    for superlocus in sorted(superloci, key=lambda item: (item.chrom, item.start0, item.end0, item.superlocus_id)):
        coverage_cache: AnchorCoverageCache = {}
        models_by_source: dict[str, list[CodingModel]] = defaultdict(list)
        for model in superlocus.models:
            models_by_source[model.source_name].append(model)
        candidates = sorted(
            superlocus.models,
            key=lambda model: (model.cds_start0, model.cds_end0, model.source_name, model.input_transcript_id),
        )
        candidate_by_id = {model.output_transcript_id: model for model in candidates}
        candidate_regions_by_id: dict[str, list[dict[str, object]]] = {}
        local_regions: list[dict[str, object]] = []
        for candidate in candidates:
            regions = _candidate_regions(
                superlocus,
                candidate,
                models_by_source,
                anchor_threshold=anchor_threshold,
                penalty_weight=penalty_weight,
                coverage_cache=coverage_cache,
            )
            candidate_regions_by_id[candidate.output_transcript_id] = regions
            local_regions.extend(regions)
        facts_by_region: dict[str, dict[str, dict[str, object]]] = {}
        for region in local_regions:
            candidate = candidate_by_id[str(region["candidate_model_id"])]
            facts_by_region[str(region["region_id"])] = {
                source: _source_region_facts(region, candidate, source_models, anchor_threshold, coverage_cache)
                for source, source_models in sorted(models_by_source.items())
            }
        _assign_bridge_groups(local_regions, facts_by_region, candidate_by_id)
        for region in local_regions:
            _score_with_revised_votes(region, facts_by_region[str(region["region_id"])])
        region_rows.extend(local_regions)
        for candidate in candidates:
            candidate_regions = candidate_regions_by_id[candidate.output_transcript_id]
            eligible = [row for row in candidate_regions if row["penalty_gate_passed"] == "true"]
            maximum = max(
                eligible,
                key=lambda row: (float(row["structural_coherence_raw_score"]), str(row["region_id"])),
                default=None,
            )
            raw_score = float(maximum["structural_coherence_raw_score"]) if maximum is not None else 0.0
            row = row_by_model.get(candidate.output_transcript_id)
            if row is None:
                continue
            row.update({
                "structural_coherence_raw_score": raw_score,
                "structural_coherence_penalty_weight": penalty_weight,
                "structural_coherence_penalty_score": penalty_weight * raw_score,
                "structural_coherence_region_count": len(eligible),
                "structural_coherence_max_region": maximum["region_id"] if maximum is not None else ".",
            })
        # Identical CDS chains have the same continuity hypothesis. Validate
        # that candidate-source exclusion and region discovery did not create
        # a silent source-specific score difference.
        by_chain: dict[tuple[str | None, tuple[Part, ...]], list[CodingModel]] = defaultdict(list)
        for candidate in candidates:
            by_chain[(candidate.sublocus_id, tuple(candidate.cds_parts))].append(candidate)
        for group in by_chain.values():
            if len(group) < 2:
                continue
            raw_values = {
                float(row_by_model[model.output_transcript_id]["structural_coherence_raw_score"])
                for model in group if model.output_transcript_id in row_by_model
            }
            if len(raw_values) > 1 and max(raw_values) - min(raw_values) > 1e-9:
                raise RuntimeError(
                    "Identical-CDS structural-coherence invariant failed for "
                    + ",".join(sorted(model.output_transcript_id for model in group))
                )
    return region_rows


def write_structural_coherence_regions_tsv(rows: list[dict[str, object]], path: Path) -> None:
    write_dict_tsv(rows, path, STRUCTURAL_COHERENCE_REGION_COLUMNS)
