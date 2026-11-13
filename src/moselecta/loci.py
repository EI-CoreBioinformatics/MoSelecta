"""Superlocus, sublocus, and source-agreement helpers."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Optional

from .gffio import sanitize_id_part
from .models import CodingModel, RelationshipEdge, SourceAutoWeight, Sublocus, Superlocus, UnionFind


def make_superlocus_id(chrom: str, start0: int, end0: int, strand: str) -> str:
    strand_label = "plus" if strand == "+" else "minus"
    return f"superlocus_{sanitize_id_part(chrom)}_{start0 + 1}_{end0}_{strand_label}"


def build_superloci(models: list[CodingModel]) -> list[Superlocus]:
    superloci: list[Superlocus] = []
    by_key: dict[tuple[str, str], list[CodingModel]] = defaultdict(list)
    for model in models:
        by_key[(model.chrom, model.strand)].append(model)

    for (chrom, strand), key_models in sorted(by_key.items()):
        key_models = sorted(key_models, key=lambda m: (m.cds_start0, m.cds_end0, m.source_name, m.input_transcript_id))
        current_models: list[CodingModel] = []
        current_start: Optional[int] = None
        current_end: Optional[int] = None

        for model in key_models:
            if current_start is None:
                current_start = model.cds_start0
                current_end = model.cds_end0
                current_models = [model]
                continue

            assert current_end is not None
            # 0-based half-open intervals overlap by >=1 bp iff next.start < current.end.
            if model.cds_start0 < current_end:
                current_end = max(current_end, model.cds_end0)
                current_models.append(model)
            else:
                superlocus_id = make_superlocus_id(chrom, current_start, current_end, strand)
                sl = Superlocus(superlocus_id, chrom, strand, current_start, current_end, current_models)
                for member in current_models:
                    member.superlocus_id = superlocus_id
                superloci.append(sl)
                current_start = model.cds_start0
                current_end = model.cds_end0
                current_models = [model]

        if current_start is not None and current_end is not None:
            superlocus_id = make_superlocus_id(chrom, current_start, current_end, strand)
            sl = Superlocus(superlocus_id, chrom, strand, current_start, current_end, current_models)
            for member in current_models:
                member.superlocus_id = superlocus_id
            superloci.append(sl)

    return sorted(superloci, key=lambda sl: (sl.chrom, sl.start0, sl.end0, sl.strand))


def superlocus_sources(sl: Superlocus) -> list[str]:
    return sorted({m.source_name for m in sl.models})


def superlocus_models_per_source(sl: Superlocus) -> dict[str, int]:
    counts: Counter[str] = Counter(m.source_name for m in sl.models)
    return dict(sorted(counts.items()))


def format_models_per_source(sl: Superlocus) -> str:
    counts = superlocus_models_per_source(sl)
    return ",".join(f"{source}:{count}" for source, count in counts.items())


def superlocus_genes_per_source(sl: Superlocus) -> dict[str, int]:
    by_source: dict[str, set[str]] = defaultdict(set)
    for model in sl.models:
        by_source[model.source_name].add(model.input_gene_id)
    return {source: len(genes) for source, genes in sorted(by_source.items())}


def format_genes_per_source(sl: Superlocus) -> str:
    counts = superlocus_genes_per_source(sl)
    return ",".join(f"{source}:{count}" for source, count in counts.items())


def superlocus_gene_count(sl: Superlocus) -> int:
    return len({(m.source_name, m.input_gene_id) for m in sl.models})


def superlocus_span_bp(sl: Superlocus) -> int:
    return sl.end0 - sl.start0


def superlocus_max_model_cds_span_bp(sl: Superlocus) -> int:
    if not sl.models:
        return 0
    return max(m.cds_span_bp for m in sl.models)


def superlocus_max_model_cds_intron_bp(sl: Superlocus) -> int:
    if not sl.models:
        return 0
    return max(m.max_cds_intron_bp for m in sl.models)


def source_signature(sl: Superlocus) -> str:
    return "+".join(superlocus_sources(sl))


def numeric_median(values: list[int]) -> float:
    if not values:
        return 0.0
    values = sorted(values)
    n = len(values)
    mid = n // 2
    if n % 2:
        return float(values[mid])
    return (values[mid - 1] + values[mid]) / 2.0


def nearest_rank_percentile(values: list[int], pct: float) -> float:
    """Return a nearest-rank percentile as an observed value."""
    if not values:
        return 0.0
    if pct < 0 or pct > 100:
        raise ValueError(f"Percentile must be between 0 and 100, got {pct}")
    values = sorted(values)
    if pct == 0:
        return float(values[0])
    rank = math.ceil((pct / 100.0) * len(values))
    index = min(max(rank - 1, 0), len(values) - 1)
    return float(values[index])


def interval_overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def reciprocal_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> tuple[int, float, float]:
    overlap = interval_overlap_bp(a_start, a_end, b_start, b_end)
    a_len = max(1, a_end - a_start)
    b_len = max(1, b_end - b_start)
    return overlap, overlap / a_len, overlap / b_len


def model_key(model: CodingModel) -> str:
    return model.output_transcript_id


def make_sublocus_id(superlocus_id: str, index: int) -> str:
    return f"{superlocus_id}.sublocus{index:03d}"


def pairwise_single_exon_chain_flag(models: list[CodingModel], min_reciprocal_overlap: float) -> bool:
    if len(models) <= 2:
        return False
    for i, a in enumerate(models):
        for b in models[i + 1:]:
            _, ro_a, ro_b = reciprocal_overlap(a.cds_start0, a.cds_end0, b.cds_start0, b.cds_end0)
            if min(ro_a, ro_b) < min_reciprocal_overlap:
                return True
    return False


def pairwise_multiexon_splice_chain_flag(models: list[CodingModel]) -> bool:
    if len(models) <= 2:
        return False
    splice_sites = {model_key(m): m.cds_splice_sites for m in models}
    for i, a in enumerate(models):
        for b in models[i + 1:]:
            if not (splice_sites[model_key(a)] & splice_sites[model_key(b)]):
                return True
    return False


def sublocus_type_for_models(models: list[CodingModel]) -> str:
    if all(m.is_single_exon for m in models):
        return "single_exon"
    if all(m.is_multi_exon for m in models):
        return "multi_exon"
    return "mixed"


def models_share_sublocus_relationship(
    a: CodingModel,
    b: CodingModel,
    *,
    single_exon_min_reciprocal_overlap: float,
) -> bool:
    """Return the canonical first-pass MoSelecta sublocus relationship.

    Keeping this predicate public lets post-selection geometry refinements test
    against the frozen final-gene partition without maintaining a second,
    subtly different definition of gene connectivity.
    """
    if a.chrom != b.chrom or a.strand != b.strand:
        return False
    if a.is_multi_exon and b.is_multi_exon:
        return bool(a.cds_splice_sites & b.cds_splice_sites)
    if a.is_single_exon and b.is_single_exon:
        _overlap, ro_a, ro_b = reciprocal_overlap(
            a.cds_start0, a.cds_end0, b.cds_start0, b.cds_end0
        )
        return min(ro_a, ro_b) >= single_exon_min_reciprocal_overlap
    return False


def build_subloci_for_superlocus(
    sl: Superlocus,
    *,
    single_exon_min_reciprocal_overlap: float,
) -> tuple[list[Sublocus], list[RelationshipEdge]]:
    """
    First-pass sublocus rules:
      - multi-exon models connect if they share >=1 exact CDS splice-site coordinate
      - single-exon models connect if reciprocal CDS-span overlap >= threshold
      - single-exon and multi-exon models are not connected in this first pass
    """
    models = sorted(sl.models, key=lambda m: (m.cds_start0, m.cds_end0, m.source_name, m.input_transcript_id))
    uf = UnionFind([model_key(m) for m in models])
    model_by_key = {model_key(m): m for m in models}
    edges: list[RelationshipEdge] = []

    # Multi-exon: splice-site index avoids all-vs-all in most loci.
    splice_index: dict[int, list[CodingModel]] = defaultdict(list)
    for model in models:
        if model.is_multi_exon:
            for site in model.cds_splice_sites:
                splice_index[site].append(model)

    seen_multi_pairs: set[tuple[str, str]] = set()
    candidate_multi_pairs: set[tuple[str, str]] = set()
    for site_models in splice_index.values():
        if len(site_models) < 2:
            continue
        site_models = sorted(site_models, key=lambda m: model_key(m))
        for i, a in enumerate(site_models):
            for b in site_models[i + 1:]:
                ka, kb = sorted([model_key(a), model_key(b)])
                candidate_multi_pairs.add((ka, kb))

    for ka, kb in sorted(candidate_multi_pairs):
        if (ka, kb) in seen_multi_pairs:
            continue
        seen_multi_pairs.add((ka, kb))
        a = model_by_key[ka]
        b = model_by_key[kb]
        shared_splice = a.cds_splice_sites & b.cds_splice_sites
        if not models_share_sublocus_relationship(
            a,
            b,
            single_exon_min_reciprocal_overlap=single_exon_min_reciprocal_overlap,
        ):
            continue
        shared_introns = set(a.cds_introns) & set(b.cds_introns)
        uf.union(ka, kb)
        edges.append(
            RelationshipEdge(
                superlocus_id=sl.superlocus_id,
                sublocus_id=".",
                model_a=a,
                model_b=b,
                edge_type="shared_cds_splice_site",
                shared_splice_site_count=len(shared_splice),
                shared_cds_intron_count=len(shared_introns),
            )
        )

    # Single-exon: sorted interval sweep.
    singles = [m for m in models if m.is_single_exon]
    singles = sorted(singles, key=lambda m: (m.cds_start0, m.cds_end0, model_key(m)))
    for i, a in enumerate(singles):
        j = i + 1
        while j < len(singles) and singles[j].cds_start0 < a.cds_end0:
            b = singles[j]
            overlap, ro_a, ro_b = reciprocal_overlap(a.cds_start0, a.cds_end0, b.cds_start0, b.cds_end0)
            if models_share_sublocus_relationship(
                a,
                b,
                single_exon_min_reciprocal_overlap=single_exon_min_reciprocal_overlap,
            ):
                uf.union(model_key(a), model_key(b))
                edges.append(
                    RelationshipEdge(
                        superlocus_id=sl.superlocus_id,
                        sublocus_id=".",
                        model_a=a,
                        model_b=b,
                        edge_type="single_single_cds_reciprocal_overlap",
                        reciprocal_overlap_a=ro_a,
                        reciprocal_overlap_b=ro_b,
                        overlap_bp=overlap,
                    )
                )
            j += 1

    components: dict[str, list[CodingModel]] = defaultdict(list)
    for model in models:
        components[uf.find(model_key(model))].append(model)

    subloci: list[Sublocus] = []
    component_models = sorted(
        components.values(),
        key=lambda ms: (min(m.cds_start0 for m in ms), max(m.cds_end0 for m in ms), min(model_key(m) for m in ms)),
    )
    key_to_sublocus: dict[str, str] = {}
    for idx, comp in enumerate(component_models, start=1):
        sub_id = make_sublocus_id(sl.superlocus_id, idx)
        start0 = min(m.cds_start0 for m in comp)
        end0 = max(m.cds_end0 for m in comp)
        stype = sublocus_type_for_models(comp)
        flags: list[str] = []
        if stype == "single_exon" and pairwise_single_exon_chain_flag(comp, single_exon_min_reciprocal_overlap):
            flags.append("single_exon_overlap_chain")
        if stype == "multi_exon" and pairwise_multiexon_splice_chain_flag(comp):
            flags.append("splice_site_chain")
        if stype == "mixed":
            flags.append("mixed_single_and_multi_exon_sublocus")
        sub = Sublocus(
            sublocus_id=sub_id,
            superlocus_id=sl.superlocus_id,
            chrom=sl.chrom,
            strand=sl.strand,
            start0=start0,
            end0=end0,
            models=sorted(comp, key=lambda m: (m.cds_start0, m.cds_end0, m.source_name, m.input_transcript_id)),
            sublocus_type=stype,
            flags=",".join(flags) if flags else ".",
        )
        for model in sub.models:
            model.sublocus_id = sub_id
            key_to_sublocus[model_key(model)] = sub_id
        subloci.append(sub)

    for edge in edges:
        # Both models are in the same component by construction.
        edge.sublocus_id = key_to_sublocus.get(model_key(edge.model_a), ".")

    return subloci, edges


def build_subloci(
    superloci: list[Superlocus],
    *,
    single_exon_min_reciprocal_overlap: float,
) -> tuple[list[Sublocus], list[RelationshipEdge]]:
    all_subloci: list[Sublocus] = []
    all_edges: list[RelationshipEdge] = []
    for sl in superloci:
        subloci, edges = build_subloci_for_superlocus(
            sl,
            single_exon_min_reciprocal_overlap=single_exon_min_reciprocal_overlap,
        )
        all_subloci.extend(subloci)
        all_edges.extend(edges)
    return all_subloci, all_edges


def coordinate_feature_key(model: CodingModel, start0: int, end0: Optional[int] = None) -> tuple:
    """Return a strand-specific coordinate feature key.

    Coordinates are 0-based half-open for intervals. For splice sites, end0 is
    omitted and the key uses the single genomic coordinate. Exact coordinate
    matching is deliberate: these are quick diagnostic agreement estimates, not
    fuzzy overlap comparisons.
    """
    if end0 is None:
        return (model.chrom, model.strand, start0)
    return (model.chrom, model.strand, start0, end0)


def safe_fraction(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def calculate_source_auto_weights(
    models: list[CodingModel],
    *,
    source_names: list[str],
    source_prior_mode: str,
    manual_source_weights: dict[str, float],
    source_prior_scale: float,
    auto_source_weight_min: float,
    auto_source_weight_max: float,
) -> list[SourceAutoWeight]:
    """Estimate source agreement weights from exact CDS coordinate features.

    The diagnostic counts use unique coordinate features per source, and feature
    support is counted by distinct source presence, not by number of models or
    isoforms. Automatic weights are only non-neutral when there are three or
    more sources; with one or two sources the auto weight is neutral (1.0).
    """
    feature_classes = [
        "cds_exons",
        "internal_cds_exons",
        "terminal_cds_exons",
        "cds_introns",
        "cds_splice_sites",
    ]
    source_features: dict[str, dict[str, set[tuple]]] = {
        source: {feature_class: set() for feature_class in feature_classes}
        for source in source_names
    }
    coding_model_counts = Counter(model.source_name for model in models)

    for model in models:
        if model.source_name not in source_features:
            source_features[model.source_name] = {feature_class: set() for feature_class in feature_classes}
        parts = sorted(model.cds_parts)
        for start0, end0 in parts:
            source_features[model.source_name]["cds_exons"].add(
                coordinate_feature_key(model, start0, end0)
            )
        if len(parts) == 1:
            start0, end0 = parts[0]
            source_features[model.source_name]["terminal_cds_exons"].add(
                coordinate_feature_key(model, start0, end0)
            )
        elif len(parts) >= 2:
            for start0, end0 in (parts[0], parts[-1]):
                source_features[model.source_name]["terminal_cds_exons"].add(
                    coordinate_feature_key(model, start0, end0)
                )
            for start0, end0 in parts[1:-1]:
                source_features[model.source_name]["internal_cds_exons"].add(
                    coordinate_feature_key(model, start0, end0)
                )
        for intron_start0, intron_end0 in model.cds_introns:
            source_features[model.source_name]["cds_introns"].add(
                coordinate_feature_key(model, intron_start0, intron_end0)
            )
        for splice_site in model.cds_splice_sites:
            source_features[model.source_name]["cds_splice_sites"].add(
                coordinate_feature_key(model, splice_site)
            )

    feature_sources: dict[str, dict[tuple, set[str]]] = {
        feature_class: defaultdict(set) for feature_class in feature_classes
    }
    for source, class_map in source_features.items():
        for feature_class, features in class_map.items():
            for feature in features:
                feature_sources[feature_class][feature].add(source)

    consensus_features: dict[str, set[tuple]] = {
        feature_class: {
            feature
            for feature, sources in feature_sources[feature_class].items()
            if len(sources) >= 2
        }
        for feature_class in feature_classes
    }

    def source_support_counts(source: str, feature_class: str) -> tuple[int, int, int, float, float]:
        features = source_features[source][feature_class]
        total = len(features)
        supported = sum(
            1 for feature in features if len(feature_sources[feature_class][feature]) >= 2
        )
        unsupported = total - supported
        supported_fraction = safe_fraction(supported, total)
        unsupported_fraction = safe_fraction(unsupported, total)
        return total, supported, unsupported, supported_fraction, unsupported_fraction

    def source_consensus_capture(source: str, feature_class: str) -> float:
        consensus = consensus_features[feature_class]
        if not consensus:
            return 0.0
        present = source_features[source][feature_class] & consensus
        return len(present) / len(consensus)

    results: list[SourceAutoWeight] = []
    enough_sources_for_auto = len(source_names) >= 3
    for source in source_names:
        (
            cds_exons,
            supported_cds_exons,
            unsupported_cds_exons,
            supported_cds_exon_fraction,
            unsupported_cds_exon_fraction,
        ) = source_support_counts(source, "cds_exons")
        (
            internal_cds_exons,
            supported_internal_cds_exons,
            unsupported_internal_cds_exons,
            supported_internal_cds_exon_fraction,
            unsupported_internal_cds_exon_fraction,
        ) = source_support_counts(source, "internal_cds_exons")
        (
            terminal_cds_exons,
            supported_terminal_cds_exons,
            unsupported_terminal_cds_exons,
            supported_terminal_cds_exon_fraction,
            unsupported_terminal_cds_exon_fraction,
        ) = source_support_counts(source, "terminal_cds_exons")
        (
            cds_introns,
            supported_cds_introns,
            unsupported_cds_introns,
            supported_cds_intron_fraction,
            unsupported_cds_intron_fraction,
        ) = source_support_counts(source, "cds_introns")
        (
            cds_splice_sites,
            supported_cds_splice_sites,
            unsupported_cds_splice_sites,
            supported_cds_splice_site_fraction,
            unsupported_cds_splice_site_fraction,
        ) = source_support_counts(source, "cds_splice_sites")

        auto_components = [
            (0.50, supported_cds_splice_site_fraction, cds_splice_sites),
            (0.30, supported_cds_intron_fraction, cds_introns),
            (0.20, supported_internal_cds_exon_fraction, internal_cds_exons),
        ]
        component_weight_sum = sum(weight for weight, _fraction, total in auto_components if total > 0)
        if component_weight_sum > 0:
            auto_agreement_score = sum(
                weight * fraction
                for weight, fraction, total in auto_components
                if total > 0
            ) / component_weight_sum
        elif cds_exons > 0:
            auto_agreement_score = supported_cds_exon_fraction
        else:
            auto_agreement_score = 1.0

        if enough_sources_for_auto:
            auto_source_weight = auto_source_weight_min + auto_agreement_score * (
                auto_source_weight_max - auto_source_weight_min
            )
            auto_source_weight_status = "estimated_from_multisource_agreement"
            auto_source_weight_reason = (
                "auto weight estimated from exact coordinate agreement across three or more input sources"
            )
        else:
            auto_source_weight = 1.0
            auto_source_weight_status = "neutral_lt3_sources"
            auto_source_weight_reason = (
                "fewer than three input sources; automatic source weight is neutral because pairwise agreement alone is not independent support"
            )

        manual_source_weight = manual_source_weights.get(source, 1.0)
        if source_prior_mode == "none":
            effective_source_weight = 0.0
        elif source_prior_mode == "auto":
            effective_source_weight = auto_source_weight
        else:
            effective_source_weight = manual_source_weight

        results.append(
            SourceAutoWeight(
                source_name=source,
                coding_models=coding_model_counts[source],
                cds_exons=cds_exons,
                supported_cds_exons=supported_cds_exons,
                unsupported_cds_exons=unsupported_cds_exons,
                supported_cds_exon_fraction=supported_cds_exon_fraction,
                unsupported_cds_exon_fraction=unsupported_cds_exon_fraction,
                internal_cds_exons=internal_cds_exons,
                supported_internal_cds_exons=supported_internal_cds_exons,
                unsupported_internal_cds_exons=unsupported_internal_cds_exons,
                supported_internal_cds_exon_fraction=supported_internal_cds_exon_fraction,
                unsupported_internal_cds_exon_fraction=unsupported_internal_cds_exon_fraction,
                terminal_cds_exons=terminal_cds_exons,
                supported_terminal_cds_exons=supported_terminal_cds_exons,
                unsupported_terminal_cds_exons=unsupported_terminal_cds_exons,
                supported_terminal_cds_exon_fraction=supported_terminal_cds_exon_fraction,
                unsupported_terminal_cds_exon_fraction=unsupported_terminal_cds_exon_fraction,
                cds_introns=cds_introns,
                supported_cds_introns=supported_cds_introns,
                unsupported_cds_introns=unsupported_cds_introns,
                supported_cds_intron_fraction=supported_cds_intron_fraction,
                unsupported_cds_intron_fraction=unsupported_cds_intron_fraction,
                cds_splice_sites=cds_splice_sites,
                supported_cds_splice_sites=supported_cds_splice_sites,
                unsupported_cds_splice_sites=unsupported_cds_splice_sites,
                supported_cds_splice_site_fraction=supported_cds_splice_site_fraction,
                unsupported_cds_splice_site_fraction=unsupported_cds_splice_site_fraction,
                consensus_cds_exon_capture=source_consensus_capture(source, "cds_exons"),
                consensus_internal_cds_exon_capture=source_consensus_capture(source, "internal_cds_exons"),
                consensus_terminal_cds_exon_capture=source_consensus_capture(source, "terminal_cds_exons"),
                consensus_cds_intron_capture=source_consensus_capture(source, "cds_introns"),
                consensus_cds_splice_site_capture=source_consensus_capture(source, "cds_splice_sites"),
                auto_agreement_score=auto_agreement_score,
                auto_source_weight=auto_source_weight,
                auto_source_weight_status=auto_source_weight_status,
                auto_source_weight_reason=auto_source_weight_reason,
                manual_source_weight=manual_source_weight,
                source_prior_mode=source_prior_mode,
                effective_source_weight=effective_source_weight,
                source_prior_score=effective_source_weight * source_prior_scale,
            )
        )
    return results
