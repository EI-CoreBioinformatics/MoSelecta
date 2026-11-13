"""Output writers and output-specific formatting helpers for MoSelecta."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

from .cds_phase import (
    RECOVERED_MISSING_SOURCE_PHASE,
    cds_phase_output_status,
    emitted_cds_phase_by_part,
)
from .five_prime_boundary import compact_five_prime_boundary_provenance
from .constants import GFF3_SOURCE, INTRINSIC_GFF_FLAG_COUNT_KEYS
from .config import value_as_positive_int
from .gffio import format_attrs, sanitize_id_part
from .loci import (
    format_genes_per_source,
    format_models_per_source,
    model_key,
    nearest_rank_percentile,
    numeric_median,
    superlocus_gene_count,
    superlocus_max_model_cds_intron_bp,
    superlocus_max_model_cds_span_bp,
    superlocus_sources,
    superlocus_span_bp,
)
from .models import (
    CodingModel,
    ExcludedModel,
    ModelIntrinsicMetrics,
    ParsedInput,
    RelationshipEdge,
    SourceAutoWeight,
    Sublocus,
    Superlocus,
    ValidationFailure,
)
from .scoring import fmt_optional_float, value_as_float
from .validation import superlocus_flags


def feature_line(
    seqid: str,
    source: str,
    feature_type: str,
    start0: int,
    end0: int,
    score: str,
    strand: str,
    phase: str,
    attrs: dict[str, object],
) -> str:
    return "\t".join(
        [
            seqid,
            source,
            feature_type,
            str(start0 + 1),
            str(end0),
            score,
            strand,
            phase,
            format_attrs(attrs),
        ]
    )


def intrinsic_model_metric_gff_attrs(row: dict[str, object]) -> dict[str, object]:
    """GFF3 attributes for intrinsic model metric flags, emitted only when flagged."""
    attrs: dict[str, object] = {}
    flag_parts: list[str] = []
    for flag_name, count_key in INTRINSIC_GFF_FLAG_COUNT_KEYS:
        count = value_as_positive_int(row.get(count_key, 0))
        if count <= 0:
            continue
        attrs[f"intrinsic_model_metric_{flag_name}_count"] = count
        flag_parts.append(f"{flag_name}:{count}")
    if not flag_parts:
        return attrs
    attrs["intrinsic_model_metric_flags"] = ",".join(flag_parts)
    penalty = row.get("intrinsic_penalty_score")
    if penalty is not None and penalty != ".":
        attrs["intrinsic_penalty_score"] = fmt_optional_float(penalty, 6) if isinstance(penalty, float) else penalty
    components = row.get("intrinsic_penalty_components")
    if components is not None and components != ".":
        attrs["intrinsic_penalty_components"] = components
    return attrs


def model_gff_scoring_attrs(
    row: dict[str, object],
    *,
    include_model_selection_components: bool = False,
) -> dict[str, object]:
    """Compact model-level scoring attributes for GFF3 gene/mRNA features.

    Intrinsic metric attributes are emitted only when at least one intrinsic flag
    is present. Low-source conflict attributes are emitted only for true conflict
    rows, so routine false/neutral rows do not clutter GFF3 output.
    """
    attrs: dict[str, object] = {}
    for key in [
        "model_selection_score",
        "model_selection_score_scaled",
        "model_selection_score_rank",
        "model_selection_order_rank",
        "exact_structure_support_score",
        "exact_cds_exon_chain_support_source_count",
    ]:
        value = row.get(key)
        if value is None or value == ".":
            continue
        attrs[key] = fmt_optional_float(value, 6) if isinstance(value, float) else value

    if value_as_float(row.get("unsupported_short_cds_exon_penalty_score")) > 0:
        for key in [
            "unsupported_short_cds_exon_penalty_score",
            "unsupported_short_cds_exon_penalty_weight",
            "unsupported_short_cds_exon_penalty_in_selection",
        ]:
            value = row.get(key)
            if value is None or value == ".":
                continue
            attrs[key] = fmt_optional_float(value, 6) if isinstance(value, float) else value

    if value_as_float(row.get("structural_coherence_penalty_score")) > 0:
        for key in [
            "structural_coherence_penalty_score",
            "structural_coherence_raw_score",
            "structural_coherence_penalty_weight",
            "structural_coherence_region_count",
            "structural_coherence_max_region",
        ]:
            value = row.get(key)
            if value is None or value == ".":
                continue
            attrs[key] = fmt_optional_float(value, 6) if isinstance(value, float) else value

    if str(row.get("sequence_integrity_status", "genome_not_configured")) != "genome_not_configured":
        for key in [
            "sequence_integrity_status",
            "sequence_integrity_penalty_score",
            "sequence_integrity_penalty_weight",
            "sequence_integrity_penalty_in_selection",
            "sequence_integrity_flags",
        ]:
            value = row.get(key)
            if value is None or value == ".":
                continue
            attrs[key] = fmt_optional_float(value, 6) if isinstance(value, float) else value

    for key in ["model_selection_tie_group", "model_selection_tie_break_reason"]:
        value = row.get(key)
        if value is not None and value != ".":
            attrs[key] = value

    if value_as_float(row.get("structural_corroboration_weight")) > 0:
        for key in [
            "corroborating_other_source_count",
            "run_possible_other_source_count",
            "structural_corroboration_score",
            "structural_corroboration_alpha",
            "structural_corroboration_weight",
            "structural_corroboration_contribution",
        ]:
            value = row.get(key)
            if value is None or value == ".":
                continue
            attrs[key] = fmt_optional_float(value, 6) if isinstance(value, float) else value

    # Include compact protein-evidence attributes only for runs where protein
    # evidence was actually configured. This keeps no-protein GFF3 outputs
    # unchanged while making final selected models easier to inspect when
    # protein scoring is active or diagnostic.
    if str(row.get("protein_evidence_status", "protein_evidence_not_configured")) != "protein_evidence_not_configured":
        for key in [
            "protein_evidence_score",
            "protein_evidence_score_scaled",
            "protein_evidence_score_weight",
            "protein_evidence_score_in_selection",
            "protein_evidence_status",
            "protein_evidence_flags",
        ]:
            value = row.get(key)
            if value is None or value == ".":
                continue
            attrs[key] = fmt_optional_float(value, 6) if isinstance(value, float) else value

    # Transcript scoring is shown compactly when transcript evidence was configured.
    # Full association and cDNA-context detail remains in the model-support TSV.
    if str(row.get("transcript_evidence_status", "transcript_evidence_not_configured")) != "transcript_evidence_not_configured":
        for key in [
            "transcript_evidence_score",
            "transcript_evidence_score_scaled",
            "transcript_evidence_score_context_adjusted",
            "transcript_evidence_score_context_adjusted_scaled",
            "transcript_evidence_score_weight",
            "transcript_evidence_score_in_selection",
            "transcript_evidence_status",
        ]:
            value = row.get(key)
            if value is None or value == ".":
                continue
            attrs[key] = fmt_optional_float(value, 6) if isinstance(value, float) else value

        context_status = str(row.get("transcript_cdna_context_status", "."))
        context_concern = str(row.get("transcript_cdna_context_pruning_concern", ".")).lower() == "true"
        if context_status == "calculated" or context_concern:
            for key in [
                "transcript_cdna_context_status",
                "transcript_cdna_context_cds_to_min_cdna_fraction",
                "transcript_cdna_context_concern_class",
                "transcript_cdna_context_pruning_concern",
            ]:
                value = row.get(key)
                if value is None or value == ".":
                    continue
                attrs[key] = fmt_optional_float(value, 6) if isinstance(value, float) else value

    if str(row.get("low_source_conflict_flag", ".")).lower() == "true":
        attrs["low_source_conflict_flag"] = "true"
        value = row.get("low_source_conflict_agreement_factor")
        if value is not None and value != ".":
            attrs["low_source_conflict_agreement_factor"] = fmt_optional_float(value, 6) if isinstance(value, float) else value
        reason = row.get("low_source_conflict_reason")
        if reason is not None and reason != ".":
            attrs["low_source_conflict_reason"] = reason

    attrs.update(intrinsic_model_metric_gff_attrs(row))

    if include_model_selection_components:
        value = row.get("model_selection_score_components")
        if value is not None and value != ".":
            attrs["model_selection_score_components"] = value
    return attrs


def model_selection_attrs(
    model: CodingModel,
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    *,
    include_model_selection_components: bool = False,
) -> dict[str, object]:
    row = model_structural_rows_by_transcript.get(model.output_transcript_id, {})
    if not row:
        return {}
    return model_gff_scoring_attrs(row, include_model_selection_components=include_model_selection_components)


def write_combined_gff3(
    superloci: list[Superlocus],
    subloci: list[Sublocus],
    path: Path,
    validation_failures: Optional[list[ValidationFailure]] = None,
    model_structural_rows_by_transcript: Optional[dict[str, dict[str, object]]] = None,
    include_model_selection_components: bool = False,
    prune_single_model_subloci: bool = False,
    scored_gff3_prune_singleton_sources: Optional[set[str]] = None,
    keep_empty_superloci: bool = False,
) -> None:
    validation_failures = validation_failures or []
    model_structural_rows_by_transcript = model_structural_rows_by_transcript or {}
    scored_gff3_prune_singleton_sources = scored_gff3_prune_singleton_sources or set()
    subloci_by_superlocus: dict[str, list[Sublocus]] = defaultdict(list)
    for sub in subloci:
        subloci_by_superlocus[sub.superlocus_id].append(sub)

    def should_prune_sublocus_from_gff3(sub: Sublocus) -> bool:
        if not prune_single_model_subloci:
            return False
        if len(sub.models) != 1:
            return False
        source = sub.models[0].source_name
        if scored_gff3_prune_singleton_sources and source not in scored_gff3_prune_singleton_sources:
            return False
        return True

    def models_for_subloci(items: list[Sublocus]) -> list[CodingModel]:
        return [model for sub in items for model in sub.models]

    def format_models_per_source_for_models(models: list[CodingModel]) -> str:
        if not models:
            return "."
        return ",".join(
            f"{source}:{count}"
            for source, count in sorted(Counter(m.source_name for m in models).items())
        )

    def format_genes_per_source_for_models(models: list[CodingModel]) -> str:
        if not models:
            return "."
        sources = sorted({m.source_name for m in models})
        return ",".join(
            f"{source}:{len({m.input_gene_id for m in models if m.source_name == source})}"
            for source in sources
        )

    def model_scoring_attrs(model: CodingModel) -> dict[str, object]:
        row = model_structural_rows_by_transcript.get(model.output_transcript_id)
        if not row:
            return {}
        return model_gff_scoring_attrs(row, include_model_selection_components=include_model_selection_components)

    with path.open("w") as out:
        out.write("##gff-version 3\n")
        for sl in superloci:
            flags = superlocus_flags(sl, validation_failures)
            original_sl_subloci = sorted(subloci_by_superlocus.get(sl.superlocus_id, []), key=lambda s: (s.start0, s.end0, s.sublocus_id))
            sl_subloci = [sub for sub in original_sl_subloci if not should_prune_sublocus_from_gff3(sub)]
            if prune_single_model_subloci and not sl_subloci and not keep_empty_superloci:
                continue
            retained_models = models_for_subloci(sl_subloci)
            sources = sorted({m.source_name for m in retained_models})
            models_per_source = format_models_per_source_for_models(retained_models)
            genes_per_source = format_genes_per_source_for_models(retained_models)
            superlocus_attrs = {
                "ID": sl.superlocus_id,
                "model_count": len(retained_models),
                "gene_count": len({(m.source_name, m.input_gene_id) for m in retained_models}),
                "source_count": len(sources),
                "sublocus_count": len(sl_subloci),
                "sources": ",".join(sources) if sources else ".",
                "models_per_source": models_per_source,
                "genes_per_source": genes_per_source,
                "superlocus_span_bp": superlocus_span_bp(sl),
                "max_model_cds_span_bp": max((m.cds_span_bp for m in retained_models), default=0),
                "max_model_cds_intron_bp": max((m.max_cds_intron_bp for m in retained_models), default=0),
                "superlocus_flags": flags,
            }
            if prune_single_model_subloci:
                pruned_subloci = [sub for sub in original_sl_subloci if should_prune_sublocus_from_gff3(sub)]
                pruned_models = models_for_subloci(pruned_subloci)
                superlocus_attrs.update({
                    "gff3_filter": "prune_singleton_subloci_from_scored_gff3",
                    "scored_gff3_prune_singleton_sources": ",".join(sorted(scored_gff3_prune_singleton_sources)) if scored_gff3_prune_singleton_sources else "ALL",
                    "original_model_count": len(sl.models),
                    "original_sublocus_count": len(original_sl_subloci),
                    "retained_model_count": len(retained_models),
                    "retained_sublocus_count": len(sl_subloci),
                    "pruned_model_count": len(pruned_models),
                    "pruned_sublocus_count": len(pruned_subloci),
                })
            out.write(
                feature_line(
                    seqid=sl.chrom,
                    source=GFF3_SOURCE,
                    feature_type="superlocus",
                    start0=sl.start0,
                    end0=sl.end0,
                    score=".",
                    strand=sl.strand,
                    phase=".",
                    attrs=superlocus_attrs,
                )
                + "\n"
            )

            for sub in sl_subloci:
                sub_sources = sorted({m.source_name for m in sub.models})
                sub_models_per_source = ",".join(
                    f"{source}:{count}"
                    for source, count in sorted(Counter(m.source_name for m in sub.models).items())
                )
                sub_genes_per_source = ",".join(
                    f"{source}:{len({m.input_gene_id for m in sub.models if m.source_name == source})}"
                    for source in sorted({m.source_name for m in sub.models})
                )
                out.write(
                    feature_line(
                        seqid=sub.chrom,
                        source=GFF3_SOURCE,
                        feature_type="sublocus",
                        start0=sub.start0,
                        end0=sub.end0,
                        score=".",
                        strand=sub.strand,
                        phase=".",
                        attrs={
                            "ID": sub.sublocus_id,
                            "Parent": sl.superlocus_id,
                            "sublocus_type": sub.sublocus_type,
                            "model_count": len(sub.models),
                            "gene_count": len({(m.source_name, m.input_gene_id) for m in sub.models}),
                            "source_count": len(sub_sources),
                            "sources": ",".join(sub_sources),
                            "models_per_source": sub_models_per_source,
                            "genes_per_source": sub_genes_per_source,
                            "sublocus_flags": sub.flags,
                        },
                    )
                    + "\n"
                )

                for model in sorted(sub.models, key=lambda m: (m.model_start0, m.model_end0, m.source_name, m.input_transcript_id)):
                    model_source = f"{GFF3_SOURCE}:{model.source_name}"
                    phase_by_part = emitted_cds_phase_by_part(model)
                    out.write(
                        feature_line(
                            seqid=model.chrom,
                            source=model_source,
                            feature_type="gene",
                            start0=model.model_start0,
                            end0=model.model_end0,
                            score=".",
                            strand=model.strand,
                            phase=".",
                            attrs={
                                **{
                                    "ID": model.output_gene_id,
                                    "Parent": sub.sublocus_id,
                                    "original_gene_id": model.input_gene_id,
                                    "original_transcript_id": model.input_transcript_id,
                                    "source_method": model.source_name,
                                    "source_file": str(model.input_file),
                                    "candidate_type": "coding_model",
                                },
                                **model_scoring_attrs(model),
                            },
                        )
                        + "\n"
                    )
                    out.write(
                        feature_line(
                            seqid=model.chrom,
                            source=model_source,
                            feature_type="mRNA",
                            start0=model.model_start0,
                            end0=model.model_end0,
                            score=model.transcript_record.score if model.transcript_record else ".",
                            strand=model.strand,
                            phase=".",
                            attrs={
                                **{
                                    "ID": model.output_transcript_id,
                                    "Parent": model.output_gene_id,
                                    "original_gene_id": model.input_gene_id,
                                    "original_transcript_id": model.input_transcript_id,
                                    "source_method": model.source_name,
                                    "source_file": str(model.input_file),
                                    "cds_phase_output_status": cds_phase_output_status(model),
                                    **(
                                        {"cds_5prime_boundary_canonicalised": provenance}
                                        if (provenance := compact_five_prime_boundary_provenance(model)) is not None
                                        else {}
                                    ),
                                },
                                **model_scoring_attrs(model),
                            },
                        )
                        + "\n"
                    )

                    generated_counts: Counter[str] = Counter()
                    for rec in model.child_records:
                        generated_counts[rec.feature_type] += 1
                        original_id = rec.id
                        if original_id:
                            child_id = f"{model.output_transcript_id}__{sanitize_id_part(original_id)}"
                        else:
                            child_id = f"{model.output_transcript_id}__{sanitize_id_part(rec.feature_type)}{generated_counts[rec.feature_type]}"
                        attrs: dict[str, object] = {
                            "ID": child_id,
                            "Parent": model.output_transcript_id,
                            "source_method": model.source_name,
                        }
                        if original_id:
                            attrs["original_id"] = original_id
                        parent_text = rec.attrs.get("Parent") or rec.attrs.get("transcript_id")
                        if parent_text:
                            attrs["original_parent"] = parent_text
                        out.write(
                            feature_line(
                                seqid=rec.seqid,
                                source=rec.source,
                                feature_type=rec.feature_type,
                                start0=rec.start0,
                                end0=rec.end0,
                                score=rec.score,
                                strand=rec.strand,
                                phase=(
                                    phase_by_part.get((rec.start0, rec.end0), rec.phase)
                                    if rec.feature_type.lower() == "cds"
                                    else rec.phase
                                ),
                                attrs=attrs,
                            )
                            + "\n"
                        )
            out.write("###\n")


def write_model_subset_gff3(
    models: list[CodingModel],
    path: Path,
    *,
    model_structural_rows_by_transcript: dict[str, dict[str, object]],
    subloci: list[Sublocus],
    selection_info_by_transcript: Optional[dict[str, dict[str, object]]] = None,
    monosublocus_rows: Optional[list[dict[str, object]]] = None,
    include_model_selection_components: bool = False,
    include_monosublocus_features: bool = False,
    selection_label: str = "sublocus_winner",
    cds_only_model_structure: bool = False,
) -> None:
    selection_info_by_transcript = selection_info_by_transcript or {}
    monosublocus_rows = monosublocus_rows or []
    sub_by_model_id = {m.output_transcript_id: sub for sub in subloci for m in sub.models}
    models_by_mono: dict[str, list[CodingModel]] = defaultdict(list)
    for model in models:
        info = selection_info_by_transcript.get(model.output_transcript_id, {})
        mono = str(info.get("monosublocus_id", "."))
        models_by_mono[mono].append(model)
    mono_rows_by_id = {str(row.get("monosublocus_id")): row for row in monosublocus_rows}

    def write_cds_only_child_features(out, model: CodingModel, model_source: str) -> None:
        phase_by_part = emitted_cds_phase_by_part(model)
        for idx, (start0, end0) in enumerate(sorted(model.cds_parts), start=1):
            exon_attrs: dict[str, object] = {
                "ID": f"{model.output_transcript_id}__cds_exon{idx}",
                "Parent": model.output_transcript_id,
                "source_method": model.source_name,
                "generated_from": "CDS",
            }
            out.write(feature_line(
                seqid=model.chrom,
                source=model_source,
                feature_type="exon",
                start0=start0,
                end0=end0,
                score=".",
                strand=model.strand,
                phase=".",
                attrs=exon_attrs,
            ) + "\n")
            cds_attrs: dict[str, object] = {
                "ID": f"{model.output_transcript_id}__CDS{idx}",
                "Parent": model.output_transcript_id,
                "source_method": model.source_name,
                "generated_from": "CDS",
            }
            out.write(feature_line(
                seqid=model.chrom,
                source=model_source,
                feature_type="CDS",
                start0=start0,
                end0=end0,
                score=".",
                strand=model.strand,
                phase=phase_by_part.get((start0, end0), "0"),
                attrs=cds_attrs,
            ) + "\n")

    def write_model(out, model: CodingModel, parent_id: Optional[str]) -> None:
        sub = sub_by_model_id.get(model.output_transcript_id)
        info = selection_info_by_transcript.get(model.output_transcript_id, {})
        base_attrs = {
            "ID": model.output_gene_id,
            "original_gene_id": model.input_gene_id,
            "original_transcript_id": model.input_transcript_id,
            "source_method": model.source_name,
            "source_file": str(model.input_file),
            "candidate_type": "coding_model",
            selection_label: "true",
        }
        if parent_id:
            base_attrs["Parent"] = parent_id
        if cds_only_model_structure:
            base_attrs["final_gff3_cds_only_structure"] = "true"
        if sub:
            base_attrs["superlocus_id"] = sub.superlocus_id
            base_attrs["sublocus_id"] = sub.sublocus_id
        for key in [
            "monosublocus_id", "monosublocus_member_count", "monosublocus_selected_count",
            "monosublocus_suppressed_count", "monosublocus_member_order_rank", "monosublocus_selected",
            "monosublocus_selection_order", "monosublocus_suppressed_by", "monosublocus_suppression_reason",
            "monosublocus_suppression_overlap_bp", "monosublocus_suppressed_model_overlap_fraction",
            "monosublocus_selected_model_overlap_fraction", "monosublocus_suppression_max_overlap_fraction",
            "final_selected", "final_selection_stage", "backfill_candidate_order_rank", "backfill_selected",
            "backfill_selection_order", "backfill_rejection_reason", "backfill_rejected_by",
            "backfill_rejected_by_source", "backfill_conflict_shared_splice_site_count",
            "backfill_conflict_overlap_bp", "backfill_candidate_overlap_fraction",
            "backfill_selected_model_overlap_fraction", "backfill_conflict_max_overlap_fraction",
            "final_prune_level", "final_prune_classes",
            "final_singleton_prune_level", "final_singleton_prune_source_scope",
            "final_singleton_pruned", "final_singleton_prune_reason",
            "final_singleton_prune_sublocus_model_count", "final_singleton_prune_sublocus_source_count",
            "final_shadow_prune_enabled", "final_shadow_pruned", "final_shadow_prune_reason",
            "final_shadow_prune_source_scope", "final_shadow_prune_require_loser_single_exon",
            "final_shadow_prune_retained_model_id", "final_shadow_prune_retained_model_source",
            "final_shadow_prune_overlap_bp", "final_shadow_prune_retained_model_overlap_fraction",
            "final_shadow_prune_removed_model_overlap_fraction", "final_shadow_prune_score_margin_scaled",
            "final_shadow_prune_removed_model_count", "final_shadow_prune_removed_model_ids",
        ]:
            value = info.get(key)
            if value is not None and value != ".":
                base_attrs[key] = fmt_optional_float(value, 6) if isinstance(value, float) else value
        base_attrs.update(model_selection_attrs(model, model_structural_rows_by_transcript, include_model_selection_components=include_model_selection_components))
        model_source = f"{GFF3_SOURCE}:{model.source_name}"
        model_feature_start0 = model.cds_start0 if cds_only_model_structure else model.model_start0
        model_feature_end0 = model.cds_end0 if cds_only_model_structure else model.model_end0
        out.write(feature_line(
            seqid=model.chrom, source=model_source, feature_type="gene", start0=model_feature_start0, end0=model_feature_end0,
            score=".", strand=model.strand, phase=".", attrs=base_attrs,
        ) + "\n")
        mrna_attrs = dict(base_attrs)
        mrna_attrs["ID"] = model.output_transcript_id
        mrna_attrs["Parent"] = model.output_gene_id
        mrna_attrs["cds_phase_output_status"] = cds_phase_output_status(model)
        boundary_provenance = compact_five_prime_boundary_provenance(model)
        if boundary_provenance is not None:
            mrna_attrs["cds_5prime_boundary_canonicalised"] = boundary_provenance
        recovery = model.cds_phase_recovery
        if recovery is not None and recovery.status == RECOVERED_MISSING_SOURCE_PHASE:
            mrna_attrs["cds_phase_recovered"] = recovery.recovered_first_phase
        out.write(feature_line(
            seqid=model.chrom, source=model_source, feature_type="mRNA", start0=model_feature_start0, end0=model_feature_end0,
            score=model.transcript_record.score if model.transcript_record else ".", strand=model.strand, phase=".", attrs=mrna_attrs,
        ) + "\n")
        if cds_only_model_structure:
            write_cds_only_child_features(out, model, model_source)
            return
        phase_by_part = emitted_cds_phase_by_part(model)
        generated_counts: Counter[str] = Counter()
        for rec in model.child_records:
            generated_counts[rec.feature_type] += 1
            original_id = rec.id
            if original_id:
                child_id = f"{model.output_transcript_id}__{sanitize_id_part(original_id)}"
            else:
                child_id = f"{model.output_transcript_id}__{sanitize_id_part(rec.feature_type)}{generated_counts[rec.feature_type]}"
            attrs: dict[str, object] = {"ID": child_id, "Parent": model.output_transcript_id, "source_method": model.source_name}
            if original_id:
                attrs["original_id"] = original_id
            parent_text = rec.attrs.get("Parent") or rec.attrs.get("transcript_id")
            if parent_text:
                attrs["original_parent"] = parent_text
            out.write(feature_line(
                seqid=rec.seqid, source=rec.source, feature_type=rec.feature_type, start0=rec.start0, end0=rec.end0,
                score=rec.score,
                strand=rec.strand,
                phase=(
                    phase_by_part.get((rec.start0, rec.end0), rec.phase)
                    if rec.feature_type.lower() == "cds"
                    else rec.phase
                ),
                attrs=attrs,
            ) + "\n")

    with path.open("w") as out:
        out.write("##gff-version 3\n")
        if include_monosublocus_features:
            for mono_id in sorted(models_by_mono, key=lambda mid: (
                str(mono_rows_by_id.get(mid, {}).get("chrom", "")),
                int(mono_rows_by_id.get(mid, {}).get("start", 0) or 0),
                int(mono_rows_by_id.get(mid, {}).get("end", 0) or 0),
                mid,
            )):
                group_models = sorted(models_by_mono[mono_id], key=lambda m: (m.cds_start0, m.cds_end0, m.output_transcript_id))
                if mono_id != "." and mono_id in mono_rows_by_id:
                    row = mono_rows_by_id[mono_id]
                    attrs = {
                        "ID": mono_id,
                        "member_count": row.get("member_count", len(group_models)),
                        "selected_count": row.get("selected_count", "."),
                        "suppressed_count": row.get("suppressed_count", "."),
                        "edge_count": row.get("edge_count", "."),
                        "min_cds_exon_overlap_fraction": row.get("min_cds_exon_overlap_fraction", "."),
                        "not_all_pairwise_overlapping": row.get("not_all_pairwise_overlapping", "."),
                    }
                    out.write(feature_line(
                        seqid=str(row.get("chrom", group_models[0].chrom)), source=GFF3_SOURCE, feature_type="monosublocus",
                        start0=int(row.get("start", group_models[0].cds_start0 + 1)) - 1, end0=int(row.get("end", group_models[0].cds_end0)),
                        score=".", strand=str(row.get("strand", group_models[0].strand)), phase=".", attrs=attrs,
                    ) + "\n")
                    parent_id = mono_id
                else:
                    parent_id = None
                for model in group_models:
                    write_model(out, model, parent_id)
        else:
            for model in sorted(models, key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.output_transcript_id)):
                write_model(out, model, None)


def write_superloci_tsv(
    superloci: list[Superlocus],
    subloci: list[Sublocus],
    path: Path,
    validation_failures: Optional[list[ValidationFailure]] = None,
) -> None:
    validation_failures = validation_failures or []
    sublocus_counts = Counter(sub.superlocus_id for sub in subloci)
    with path.open("w") as out:
        out.write(
            "superlocus_id\tchrom\tstart\tend\tstrand\tsuperlocus_span_bp\tmodel_count\tgene_count\tsource_count\tsublocus_count\tsources\tmodels_per_source\tgenes_per_source\tinput_transcript_count\tmax_model_cds_span_bp\tmax_model_cds_intron_bp\tsuperlocus_flags\n"
        )
        for sl in superloci:
            sources = superlocus_sources(sl)
            models_per_source = format_models_per_source(sl)
            genes_per_source = format_genes_per_source(sl)
            txs = sorted({f"{m.source_name}:{m.input_transcript_id}" for m in sl.models})
            flags = superlocus_flags(sl, validation_failures)
            out.write(
                "\t".join(
                    [
                        sl.superlocus_id,
                        sl.chrom,
                        str(sl.start0 + 1),
                        str(sl.end0),
                        sl.strand,
                        str(superlocus_span_bp(sl)),
                        str(len(sl.models)),
                        str(superlocus_gene_count(sl)),
                        str(len(sources)),
                        str(sublocus_counts.get(sl.superlocus_id, 0)),
                        ",".join(sources),
                        models_per_source,
                        genes_per_source,
                        str(len(txs)),
                        str(superlocus_max_model_cds_span_bp(sl)),
                        str(superlocus_max_model_cds_intron_bp(sl)),
                        flags,
                    ]
                )
                + "\n"
            )


MODEL_TSV_COLUMNS = (
    "superlocus_id",
    "sublocus_id",
    "source_name",
    "input_file",
    "input_gene_id",
    "input_transcript_id",
    "output_gene_id",
    "output_transcript_id",
    "chrom",
    "cds_start",
    "cds_end",
    "strand",
    "cds_length",
    "cds_exon_count",
    "cds_intron_count",
    "model_start",
    "model_end",
)


def model_tsv_row(model: CodingModel) -> list[str]:
    return [
        model.superlocus_id or ".",
        model.sublocus_id or ".",
        model.source_name,
        str(model.input_file),
        model.input_gene_id,
        model.input_transcript_id,
        model.output_gene_id,
        model.output_transcript_id,
        model.chrom,
        str(model.cds_start0 + 1),
        str(model.cds_end0),
        model.strand,
        str(model.cds_length),
        str(model.cds_exon_count),
        str(max(0, model.cds_exon_count - 1)),
        str(model.model_start0 + 1),
        str(model.model_end0),
    ]


def write_models_tsv(models: list[CodingModel], path: Path) -> None:
    with path.open("w") as out:
        out.write("\t".join(MODEL_TSV_COLUMNS) + "\n")
        for model in sorted(models, key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.source_name, m.input_transcript_id)):
            out.write("\t".join(model_tsv_row(model)) + "\n")


def write_members_tsv(superloci: list[Superlocus], path: Path) -> None:
    with path.open("w") as out:
        out.write(
            "superlocus_id\tsublocus_id\tsuperlocus_models_per_source\tsuperlocus_genes_per_source\tsource_name\tinput_file\tinput_gene_id\tinput_transcript_id\toutput_gene_id\toutput_transcript_id\tchrom\tcds_start\tcds_end\tstrand\tcds_length\tcds_exon_count\tcds_intron_count\tmodel_start\tmodel_end\n"
        )
        for sl in superloci:
            models_per_source = format_models_per_source(sl)
            genes_per_source = format_genes_per_source(sl)
            for model in sorted(sl.models, key=lambda m: (m.cds_start0, m.cds_end0, m.source_name, m.input_transcript_id)):
                row = model_tsv_row(model)
                out.write("\t".join([model.superlocus_id or ".", model.sublocus_id or ".", models_per_source, genes_per_source] + row[2:]) + "\n")


def write_subloci_tsv(subloci: list[Sublocus], path: Path) -> None:
    with path.open("w") as out:
        out.write(
            "sublocus_id\tsuperlocus_id\tchrom\tstart\tend\tstrand\tsublocus_type\tmodel_count\tgene_count\tsource_count\tsources\tmodels_per_source\tgenes_per_source\tsublocus_flags\n"
        )
        for sub in sorted(subloci, key=lambda s: (s.chrom, s.start0, s.end0, s.sublocus_id)):
            sources = sorted({m.source_name for m in sub.models})
            models_per_source = ",".join(
                f"{source}:{count}" for source, count in sorted(Counter(m.source_name for m in sub.models).items())
            )
            genes_per_source = ",".join(
                f"{source}:{len({m.input_gene_id for m in sub.models if m.source_name == source})}"
                for source in sources
            )
            out.write(
                "\t".join(
                    [
                        sub.sublocus_id,
                        sub.superlocus_id,
                        sub.chrom,
                        str(sub.start0 + 1),
                        str(sub.end0),
                        sub.strand,
                        sub.sublocus_type,
                        str(len(sub.models)),
                        str(len({(m.source_name, m.input_gene_id) for m in sub.models})),
                        str(len(sources)),
                        ",".join(sources),
                        models_per_source,
                        genes_per_source,
                        sub.flags,
                    ]
                )
                + "\n"
            )


def write_sublocus_members_tsv(subloci: list[Sublocus], path: Path) -> None:
    with path.open("w") as out:
        out.write(
            "sublocus_id\tsuperlocus_id\tsublocus_type\tsublocus_flags\tsource_name\tinput_file\tinput_gene_id\tinput_transcript_id\toutput_gene_id\toutput_transcript_id\tchrom\tcds_start\tcds_end\tstrand\tcds_length\tcds_exon_count\tcds_intron_count\n"
        )
        for sub in sorted(subloci, key=lambda s: (s.chrom, s.start0, s.end0, s.sublocus_id)):
            for model in sub.models:
                out.write(
                    "\t".join(
                        [
                            sub.sublocus_id,
                            sub.superlocus_id,
                            sub.sublocus_type,
                            sub.flags,
                            model.source_name,
                            str(model.input_file),
                            model.input_gene_id,
                            model.input_transcript_id,
                            model.output_gene_id,
                            model.output_transcript_id,
                            model.chrom,
                            str(model.cds_start0 + 1),
                            str(model.cds_end0),
                            model.strand,
                            str(model.cds_length),
                            str(model.cds_exon_count),
                            str(len(model.cds_introns)),
                        ]
                    )
                    + "\n"
                )


def write_relationship_edges_tsv(edges: list[RelationshipEdge], path: Path) -> None:
    with path.open("w") as out:
        out.write(
            "superlocus_id\tsublocus_id\tedge_type\tmodel_a_source\tmodel_a_gene_id\tmodel_a_transcript_id\tmodel_b_source\tmodel_b_gene_id\tmodel_b_transcript_id\tshared_splice_site_count\tshared_cds_intron_count\toverlap_bp\treciprocal_overlap_a\treciprocal_overlap_b\n"
        )
        for edge in sorted(edges, key=lambda e: (e.superlocus_id, e.sublocus_id, e.edge_type, model_key(e.model_a), model_key(e.model_b))):
            out.write(
                "\t".join(
                    [
                        edge.superlocus_id,
                        edge.sublocus_id,
                        edge.edge_type,
                        edge.model_a.source_name,
                        edge.model_a.input_gene_id,
                        edge.model_a.input_transcript_id,
                        edge.model_b.source_name,
                        edge.model_b.input_gene_id,
                        edge.model_b.input_transcript_id,
                        str(edge.shared_splice_site_count),
                        str(edge.shared_cds_intron_count),
                        "." if edge.overlap_bp is None else str(edge.overlap_bp),
                        "." if edge.reciprocal_overlap_a is None else f"{edge.reciprocal_overlap_a:.6f}",
                        "." if edge.reciprocal_overlap_b is None else f"{edge.reciprocal_overlap_b:.6f}",
                    ]
                )
                + "\n"
            )


def write_model_intrinsic_metrics_tsv(metrics: list[ModelIntrinsicMetrics], path: Path) -> None:
    with path.open("w") as out:
        out.write(
            "superlocus_id\tsublocus_id\tsource_name\tinput_file\tinput_gene_id\tinput_transcript_id\tchrom\tcds_start\tcds_end\tstrand\tcds_length\tcds_length_mod3\tcds_exon_count\tcds_intron_count\tcds_span_bp\tmin_cds_exon_bp\tmax_cds_exon_bp\tmedian_cds_exon_bp\tmin_terminal_cds_exon_bp\tmin_internal_cds_exon_bp\tmedian_terminal_cds_exon_bp\tmedian_internal_cds_exon_bp\tshort_terminal_cds_exon_count\tshort_internal_cds_exon_count\tmin_cds_intron_bp\tmax_cds_intron_bp\tmedian_cds_intron_bp\tshort_cds_intron_count\tlong_cds_intron_count\tmodel_metric_flags\n"
        )
        for metric in metrics:
            model = metric.model
            out.write(
                "\t".join(
                    [
                        model.superlocus_id or ".",
                        model.sublocus_id or ".",
                        model.source_name,
                        str(model.input_file),
                        model.input_gene_id,
                        model.input_transcript_id,
                        model.chrom,
                        str(model.cds_start0 + 1),
                        str(model.cds_end0),
                        model.strand,
                        str(metric.cds_length),
                        str(metric.cds_length_mod3),
                        str(metric.cds_exon_count),
                        str(metric.cds_intron_count),
                        str(metric.cds_span_bp),
                        str(metric.min_cds_exon_bp),
                        str(metric.max_cds_exon_bp),
                        f"{metric.median_cds_exon_bp:.3f}",
                        str(metric.min_terminal_cds_exon_bp),
                        str(metric.min_internal_cds_exon_bp),
                        f"{metric.median_terminal_cds_exon_bp:.3f}",
                        f"{metric.median_internal_cds_exon_bp:.3f}",
                        str(metric.short_terminal_cds_exon_count),
                        str(metric.short_internal_cds_exon_count),
                        str(metric.min_cds_intron_bp),
                        str(metric.max_cds_intron_bp),
                        f"{metric.median_cds_intron_bp:.3f}",
                        str(metric.short_cds_intron_count),
                        str(metric.long_cds_intron_count),
                        metric.flags,
                    ]
                )
                + "\n"
            )


def write_source_auto_weights_tsv(source_weights: list[SourceAutoWeight], path: Path) -> None:
    with path.open("w") as out:
        out.write(
            "source_name\tcoding_models\t"
            "cds_exons\tsupported_cds_exons\tunsupported_cds_exons\tsupported_cds_exon_fraction\tunsupported_cds_exon_fraction\t"
            "internal_cds_exons\tsupported_internal_cds_exons\tunsupported_internal_cds_exons\tsupported_internal_cds_exon_fraction\tunsupported_internal_cds_exon_fraction\t"
            "terminal_cds_exons\tsupported_terminal_cds_exons\tunsupported_terminal_cds_exons\tsupported_terminal_cds_exon_fraction\tunsupported_terminal_cds_exon_fraction\t"
            "cds_introns\tsupported_cds_introns\tunsupported_cds_introns\tsupported_cds_intron_fraction\tunsupported_cds_intron_fraction\t"
            "cds_splice_sites\tsupported_cds_splice_sites\tunsupported_cds_splice_sites\tsupported_cds_splice_site_fraction\tunsupported_cds_splice_site_fraction\t"
            "consensus_cds_exon_capture\tconsensus_internal_cds_exon_capture\tconsensus_terminal_cds_exon_capture\t"
            "consensus_cds_intron_capture\tconsensus_cds_splice_site_capture\t"
            "auto_agreement_score\tauto_source_weight\tauto_source_weight_status\tauto_source_weight_reason\tmanual_source_weight\tsource_prior_mode\teffective_source_weight\tsource_prior_score\n"
        )
        for item in source_weights:
            out.write(
                "\t".join(
                    [
                        item.source_name,
                        str(item.coding_models),
                        str(item.cds_exons),
                        str(item.supported_cds_exons),
                        str(item.unsupported_cds_exons),
                        f"{item.supported_cds_exon_fraction:.6f}",
                        f"{item.unsupported_cds_exon_fraction:.6f}",
                        str(item.internal_cds_exons),
                        str(item.supported_internal_cds_exons),
                        str(item.unsupported_internal_cds_exons),
                        f"{item.supported_internal_cds_exon_fraction:.6f}",
                        f"{item.unsupported_internal_cds_exon_fraction:.6f}",
                        str(item.terminal_cds_exons),
                        str(item.supported_terminal_cds_exons),
                        str(item.unsupported_terminal_cds_exons),
                        f"{item.supported_terminal_cds_exon_fraction:.6f}",
                        f"{item.unsupported_terminal_cds_exon_fraction:.6f}",
                        str(item.cds_introns),
                        str(item.supported_cds_introns),
                        str(item.unsupported_cds_introns),
                        f"{item.supported_cds_intron_fraction:.6f}",
                        f"{item.unsupported_cds_intron_fraction:.6f}",
                        str(item.cds_splice_sites),
                        str(item.supported_cds_splice_sites),
                        str(item.unsupported_cds_splice_sites),
                        f"{item.supported_cds_splice_site_fraction:.6f}",
                        f"{item.unsupported_cds_splice_site_fraction:.6f}",
                        f"{item.consensus_cds_exon_capture:.6f}",
                        f"{item.consensus_internal_cds_exon_capture:.6f}",
                        f"{item.consensus_terminal_cds_exon_capture:.6f}",
                        f"{item.consensus_cds_intron_capture:.6f}",
                        f"{item.consensus_cds_splice_site_capture:.6f}",
                        f"{item.auto_agreement_score:.6f}",
                        f"{item.auto_source_weight:.6f}",
                        item.auto_source_weight_status,
                        item.auto_source_weight_reason,
                        f"{item.manual_source_weight:.6f}",
                        item.source_prior_mode,
                        f"{item.effective_source_weight:.6f}",
                        f"{item.source_prior_score:.6f}",
                    ]
                )
                + "\n"
            )


EXCLUDED_MODEL_TSV_COLUMNS = (
    "source_name",
    "input_file",
    "input_gene_id",
    "input_transcript_id",
    "reason",
    "cds_count",
    "cds_chroms",
    "cds_strands",
    "cds_span_bp",
    "max_cds_intron_bp",
)


def write_excluded_tsv(excluded: list[ExcludedModel], path: Path) -> None:
    with path.open("w") as out:
        out.write("\t".join(EXCLUDED_MODEL_TSV_COLUMNS) + "\n")
        for item in excluded:
            out.write(
                "\t".join(
                    [
                        item.source_name,
                        str(item.input_file),
                        item.input_gene_id,
                        item.input_transcript_id,
                        item.reason,
                        str(item.cds_count),
                        ",".join(item.cds_chroms),
                        ",".join(item.cds_strands),
                        "." if item.cds_span_bp is None else str(item.cds_span_bp),
                        "." if item.max_cds_intron_bp is None else str(item.max_cds_intron_bp),
                    ]
                )
                + "\n"
            )


def write_other_features_gff3(parsed_inputs: list[ParsedInput], path: Path) -> None:
    with path.open("w") as out:
        out.write("##gff-version 3\n")
        for parsed in parsed_inputs:
            out.write(f"# Unused/non-CDS/excluded features from {parsed.input_file} source_name={parsed.source_name}\n")
            for rec in parsed.records:
                if rec.record_index not in parsed.used_record_indices:
                    out.write(rec.raw_line + "\n")


def write_input_summary_tsv(parsed_inputs: list[ParsedInput], path: Path) -> None:
    with path.open("w") as out:
        out.write(
            "source_name\tinput_file\ttotal_features\tcoding_models\texcluded_models\tunused_features\n"
        )
        for parsed in parsed_inputs:
            unused = len(parsed.records) - len(parsed.used_record_indices)
            out.write(
                "\t".join(
                    [
                        parsed.source_name,
                        str(parsed.input_file),
                        str(len(parsed.records)),
                        str(len(parsed.coding_models)),
                        str(len(parsed.excluded_models)),
                        str(unused),
                    ]
                )
                + "\n"
            )


def write_validation_failures_tsv(failures: list[ValidationFailure], path: Path) -> None:
    with path.open("w") as out:
        out.write(
            "scope\tfailure_type\taction\tvalue\tthreshold\tsuperlocus_id\tchrom\tstart\tend\tstrand\tsource_name\tinput_file\tinput_gene_id\tinput_transcript_id\tdetails\n"
        )
        for failure in failures:
            out.write(
                "\t".join(
                    [
                        failure.scope,
                        failure.failure_type,
                        failure.action,
                        str(failure.value),
                        str(failure.threshold),
                        failure.superlocus_id,
                        failure.chrom,
                        "." if failure.start0 is None else str(failure.start0 + 1),
                        "." if failure.end0 is None else str(failure.end0),
                        failure.strand,
                        failure.source_name,
                        failure.input_file,
                        failure.input_gene_id,
                        failure.input_transcript_id,
                        failure.details,
                    ]
                )
                + "\n"
            )


def collect_cds_intron_lengths_by_source(models: list[CodingModel]) -> dict[str, list[int]]:
    lengths_by_source: dict[str, list[int]] = defaultdict(list)
    for model in models:
        lengths_by_source[model.source_name].extend(model.cds_intron_lengths)
    return lengths_by_source


def format_intron_stats(lengths: list[int]) -> str:
    if not lengths:
        return "n=0"
    return (
        f"n={len(lengths):,}; "
        f"min={min(lengths):,}; "
        f"median={numeric_median(lengths):.1f}; "
        f"p99={nearest_rank_percentile(lengths, 99):.1f}; "
        f"p99.99={nearest_rank_percentile(lengths, 99.99):.1f}; "
        f"p99.999={nearest_rank_percentile(lengths, 99.999):.1f}; "
        f"max={max(lengths):,}"
    )


def resolve_long_cds_intron_cutoff(
    models: list[CodingModel],
    *,
    user_long_cds_intron_bp: Optional[int],
    percentile_value: float,
) -> tuple[Optional[float], str]:
    if user_long_cds_intron_bp is not None:
        if user_long_cds_intron_bp < 1:
            raise ValueError("--long-cds-intron-bp must be >= 1")
        return float(user_long_cds_intron_bp), "user --long-cds-intron-bp"
    all_introns: list[int] = []
    for model in models:
        all_introns.extend(model.cds_intron_lengths)
    if not all_introns:
        return None, "not applied; no CDS introns in valid input models"
    return nearest_rank_percentile(all_introns, percentile_value), f"all valid input model intron {percentile_value:g}th percentile"
