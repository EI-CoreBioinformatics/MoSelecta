"""Intrinsic model metrics and validation helpers for MoSelecta."""

from __future__ import annotations

from typing import Optional

from .loci import (
    format_genes_per_source,
    numeric_median,
    superlocus_gene_count,
    superlocus_genes_per_source,
    superlocus_span_bp,
)
from .models import (
    CodingModel,
    ExcludedModel,
    ModelIntrinsicMetrics,
    ParsedInput,
    Superlocus,
    ValidationFailure,
)


def metric_flags(
    model: CodingModel,
    *,
    min_cds_length_bp: int,
    min_terminal_cds_exon_bp: int,
    min_internal_cds_exon_bp: int,
    min_cds_intron_bp: int,
    long_cds_intron_bp: Optional[float],
) -> str:
    flags: list[str] = []
    if model.cds_length % 3 != 0:
        flags.append("cds_length_not_multiple_of_3")
    if model.cds_length < min_cds_length_bp:
        flags.append("short_cds")
    terminal_exon_lengths = model.terminal_cds_exon_lengths
    internal_exon_lengths = model.internal_cds_exon_lengths
    intron_lengths = model.cds_intron_lengths
    if terminal_exon_lengths and min(terminal_exon_lengths) < min_terminal_cds_exon_bp:
        flags.append("short_terminal_cds_exon")
    if internal_exon_lengths and min(internal_exon_lengths) < min_internal_cds_exon_bp:
        flags.append("short_internal_cds_exon")
    if intron_lengths and min(intron_lengths) < min_cds_intron_bp:
        flags.append("short_cds_intron")
    if long_cds_intron_bp is not None and intron_lengths and max(intron_lengths) > long_cds_intron_bp:
        flags.append("long_cds_intron")
    if model.is_single_exon:
        flags.append("single_exon_model")
    return ",".join(flags) if flags else "."


def calculate_model_intrinsic_metrics(
    model: CodingModel,
    *,
    min_cds_length_bp: int,
    min_terminal_cds_exon_bp: int,
    min_internal_cds_exon_bp: int,
    min_cds_intron_bp: int,
    long_cds_intron_bp: Optional[float],
) -> ModelIntrinsicMetrics:
    exon_lengths = model.cds_exon_lengths
    terminal_exon_lengths = model.terminal_cds_exon_lengths
    internal_exon_lengths = model.internal_cds_exon_lengths
    intron_lengths = model.cds_intron_lengths
    short_terminal_cds_exon_count = sum(
        1 for length in terminal_exon_lengths if length < min_terminal_cds_exon_bp
    )
    short_internal_cds_exon_count = sum(
        1 for length in internal_exon_lengths if length < min_internal_cds_exon_bp
    )
    short_cds_intron_count = sum(
        1 for length in intron_lengths if length < min_cds_intron_bp
    )
    long_cds_intron_count = (
        sum(1 for length in intron_lengths if length > long_cds_intron_bp)
        if long_cds_intron_bp is not None
        else 0
    )

    return ModelIntrinsicMetrics(
        model=model,
        cds_length=model.cds_length,
        cds_length_mod3=model.cds_length % 3,
        cds_exon_count=model.cds_exon_count,
        cds_intron_count=len(model.cds_introns),
        cds_span_bp=model.cds_span_bp,
        min_cds_exon_bp=min(exon_lengths) if exon_lengths else 0,
        max_cds_exon_bp=max(exon_lengths) if exon_lengths else 0,
        median_cds_exon_bp=numeric_median(exon_lengths),
        min_terminal_cds_exon_bp=min(terminal_exon_lengths) if terminal_exon_lengths else 0,
        min_internal_cds_exon_bp=min(internal_exon_lengths) if internal_exon_lengths else 0,
        median_terminal_cds_exon_bp=numeric_median(terminal_exon_lengths),
        median_internal_cds_exon_bp=numeric_median(internal_exon_lengths),
        short_terminal_cds_exon_count=short_terminal_cds_exon_count,
        short_internal_cds_exon_count=short_internal_cds_exon_count,
        min_cds_intron_bp=min(intron_lengths) if intron_lengths else 0,
        max_cds_intron_bp=max(intron_lengths) if intron_lengths else 0,
        median_cds_intron_bp=numeric_median(intron_lengths),
        short_cds_intron_count=short_cds_intron_count,
        long_cds_intron_count=long_cds_intron_count,
        flags=metric_flags(
            model,
            min_cds_length_bp=min_cds_length_bp,
            min_terminal_cds_exon_bp=min_terminal_cds_exon_bp,
            min_internal_cds_exon_bp=min_internal_cds_exon_bp,
            min_cds_intron_bp=min_cds_intron_bp,
            long_cds_intron_bp=long_cds_intron_bp,
        ),
    )


def calculate_all_model_intrinsic_metrics(
    models: list[CodingModel],
    *,
    min_cds_length_bp: int,
    min_terminal_cds_exon_bp: int,
    min_internal_cds_exon_bp: int,
    min_cds_intron_bp: int,
    long_cds_intron_bp: Optional[float],
) -> list[ModelIntrinsicMetrics]:
    return [
        calculate_model_intrinsic_metrics(
            model,
            min_cds_length_bp=min_cds_length_bp,
            min_terminal_cds_exon_bp=min_terminal_cds_exon_bp,
            min_internal_cds_exon_bp=min_internal_cds_exon_bp,
            min_cds_intron_bp=min_cds_intron_bp,
            long_cds_intron_bp=long_cds_intron_bp,
        )
        for model in sorted(models, key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.source_name, m.input_transcript_id))
    ]


def model_validation_failures(
    model: CodingModel,
    *,
    max_model_cds_span_bp: int,
    max_model_cds_intron_bp: int,
    action: str,
) -> list[ValidationFailure]:
    failures: list[ValidationFailure] = []
    if model.cds_span_bp > max_model_cds_span_bp:
        failures.append(
            ValidationFailure(
                scope="model",
                failure_type="model_cds_span_bp_exceeds_threshold",
                action=action,
                value=model.cds_span_bp,
                threshold=max_model_cds_span_bp,
                chrom=model.chrom,
                start0=model.cds_start0,
                end0=model.cds_end0,
                strand=model.strand,
                source_name=model.source_name,
                input_file=str(model.input_file),
                input_gene_id=model.input_gene_id,
                input_transcript_id=model.input_transcript_id,
                details=f"cds_span_bp={model.cds_span_bp}",
            )
        )
    if model.max_cds_intron_bp > max_model_cds_intron_bp:
        failures.append(
            ValidationFailure(
                scope="model",
                failure_type="model_cds_intron_bp_exceeds_threshold",
                action=action,
                value=model.max_cds_intron_bp,
                threshold=max_model_cds_intron_bp,
                chrom=model.chrom,
                start0=model.cds_start0,
                end0=model.cds_end0,
                strand=model.strand,
                source_name=model.source_name,
                input_file=str(model.input_file),
                input_gene_id=model.input_gene_id,
                input_transcript_id=model.input_transcript_id,
                details=f"max_cds_intron_bp={model.max_cds_intron_bp}",
            )
        )
    return failures


def apply_model_validation(
    parsed_inputs: list[ParsedInput],
    *,
    action: str,
    max_model_cds_span_bp: int,
    max_model_cds_intron_bp: int,
) -> list[ValidationFailure]:
    all_failures: list[ValidationFailure] = []
    if action not in {"exclude", "stop", "warn"}:
        raise ValueError(f"Unsupported model validation action: {action}")

    for parsed in parsed_inputs:
        kept: list[CodingModel] = []
        for model in parsed.coding_models:
            failures = model_validation_failures(
                model,
                max_model_cds_span_bp=max_model_cds_span_bp,
                max_model_cds_intron_bp=max_model_cds_intron_bp,
                action=action if action != "exclude" else "excluded",
            )
            all_failures.extend(failures)
            if failures and action == "exclude":
                reasons = sorted({f.failure_type for f in failures})
                parsed.excluded_models.append(
                    ExcludedModel(
                        source_name=model.source_name,
                        input_file=model.input_file,
                        input_gene_id=model.input_gene_id,
                        input_transcript_id=model.input_transcript_id,
                        reason=";".join(reasons),
                        cds_count=len(model.cds_parts),
                        cds_chroms=[model.chrom],
                        cds_strands=[model.strand],
                        cds_span_bp=model.cds_span_bp,
                        max_cds_intron_bp=model.max_cds_intron_bp,
                    )
                )
                # These records are no longer used in candidate model output, so allow
                # them to appear in the other/excluded feature file for inspection.
                parsed.used_record_indices.difference_update(model.used_record_indices)
            else:
                kept.append(model)
        if action == "exclude":
            parsed.coding_models = sorted(
                kept,
                key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.source_name, m.input_transcript_id),
            )
    return all_failures


def superlocus_validation_failures(
    sl: Superlocus,
    *,
    action: str,
    max_superlocus_span_bp: int,
    max_superlocus_gene_count: int,
    max_genes_per_source_in_superlocus: int,
) -> list[ValidationFailure]:
    failures: list[ValidationFailure] = []
    span = superlocus_span_bp(sl)
    gene_count = superlocus_gene_count(sl)
    genes_per_source = superlocus_genes_per_source(sl)
    max_genes_for_source = max(genes_per_source.values(), default=0)

    if span > max_superlocus_span_bp:
        failures.append(
            ValidationFailure(
                scope="superlocus",
                failure_type="superlocus_span_bp_exceeds_threshold",
                action=action,
                value=span,
                threshold=max_superlocus_span_bp,
                superlocus_id=sl.superlocus_id,
                chrom=sl.chrom,
                start0=sl.start0,
                end0=sl.end0,
                strand=sl.strand,
                details=f"superlocus_span_bp={span};model_count={len(sl.models)};gene_count={gene_count}",
            )
        )
    if gene_count > max_superlocus_gene_count:
        failures.append(
            ValidationFailure(
                scope="superlocus",
                failure_type="superlocus_gene_count_exceeds_threshold",
                action=action,
                value=gene_count,
                threshold=max_superlocus_gene_count,
                superlocus_id=sl.superlocus_id,
                chrom=sl.chrom,
                start0=sl.start0,
                end0=sl.end0,
                strand=sl.strand,
                details=f"gene_count={gene_count};genes_per_source={format_genes_per_source(sl)}",
            )
        )
    if max_genes_for_source > max_genes_per_source_in_superlocus:
        offending = [
            f"{source}:{count}"
            for source, count in sorted(genes_per_source.items())
            if count > max_genes_per_source_in_superlocus
        ]
        failures.append(
            ValidationFailure(
                scope="superlocus",
                failure_type="genes_per_source_in_superlocus_exceeds_threshold",
                action=action,
                value=max_genes_for_source,
                threshold=max_genes_per_source_in_superlocus,
                superlocus_id=sl.superlocus_id,
                chrom=sl.chrom,
                start0=sl.start0,
                end0=sl.end0,
                strand=sl.strand,
                details=f"offending_sources={','.join(offending)};genes_per_source={format_genes_per_source(sl)}",
            )
        )
    return failures


def collect_superlocus_validation_failures(
    superloci: list[Superlocus],
    *,
    action: str,
    max_superlocus_span_bp: int,
    max_superlocus_gene_count: int,
    max_genes_per_source_in_superlocus: int,
) -> list[ValidationFailure]:
    failures: list[ValidationFailure] = []
    for sl in superloci:
        failures.extend(
            superlocus_validation_failures(
                sl,
                action=action,
                max_superlocus_span_bp=max_superlocus_span_bp,
                max_superlocus_gene_count=max_superlocus_gene_count,
                max_genes_per_source_in_superlocus=max_genes_per_source_in_superlocus,
            )
        )
    return failures


def superlocus_flags(sl: Superlocus, failures: list[ValidationFailure]) -> str:
    flags = sorted(
        {f.failure_type for f in failures if f.scope == "superlocus" and f.superlocus_id == sl.superlocus_id}
    )
    return ",".join(flags) if flags else "."

