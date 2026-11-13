"""Writers for the post-selection MoSelecta annotation layer."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Mapping, TextIO

from .cds_phase import (
    RECOVERED_MISSING_SOURCE_PHASE,
    cds_phase_output_status,
    emitted_cds_phase_by_part,
)
from .constants import CHILD_FEATURE_ORDER, GFF3_SOURCE
from .final_annotation import FinalAnnotation, FinalGene, FinalTranscript
from .gffio import sanitize_id_part
from .intervals import merge_intervals
from .models import CodingModel, FeatureRecord
from .terminal_stop import compact_terminal_stop_provenance
from .five_prime_boundary import compact_five_prime_boundary_provenance
from .outputs import feature_line
from .scoring import fmt_optional_float
from .genome_sequence import IndexedFasta
from .utr_validation import (
    Side,
    UTRSideValidation,
    side_has_utr,
    utr_side_introns,
    validate_utr_side,
)

if TYPE_CHECKING:
    from .alternative_splice_selection import AlternativeSpliceSelectionConfig


# Separator between a parent ID and its child-feature suffix (exon1, CDS1, ...).
# Only the public final-annotation files, standardised_primary_models.gff3,
# final_annotation.gff3 and final_annotation.with_utrs.gff3, use the dot
# (MoSelecta_000001.1.exon1).  Every other GFF3 output keeps the double
# underscore.
DEFAULT_CHILD_ID_SEPARATOR = "__"
PUBLIC_CHILD_ID_SEPARATOR = "."

SUPPORTED_FINAL_CHILD_TYPES = frozenset({
    "exon",
    "cds",
    "start_codon",
    "stop_codon",
    "utr",
    "five_prime_utr",
    "three_prime_utr",
})


@dataclass(frozen=True, slots=True)
class FinalAnnotationOutputMetadata:
    """Run-state metadata for the final annotation."""

    alternative_splice_enabled: bool
    selector_evaluated: bool
    alternatives_emitted: bool
    selection_config: AlternativeSpliceSelectionConfig | None = None

    def __post_init__(self) -> None:
        if self.selector_evaluated != (self.selection_config is not None):
            raise ValueError(
                "selector_evaluated must match selection_config availability"
            )
        if self.alternatives_emitted and not (
            self.alternative_splice_enabled and self.selector_evaluated
        ):
            raise ValueError(
                "alternatives_emitted requires enabled and evaluated selection"
            )


def _ordered_child_records(
    model: CodingModel,
    *,
    supported_only: bool = False,
    records: Iterable[FeatureRecord] | None = None,
) -> list[FeatureRecord]:
    return sorted(
        (
            record
            for record in (model.child_records if records is None else records)
            if not supported_only
            or record.feature_type.lower() in SUPPORTED_FINAL_CHILD_TYPES
        ),
        key=lambda record: (
            record.start0,
            CHILD_FEATURE_ORDER.get(record.feature_type, 99),
            record.end0,
            record.line_number,
            record.record_index,
            record.feature_type,
        ),
    )


def write_source_model_children(
    out: TextIO,
    model: CodingModel,
    parent_id: str,
    *,
    supported_only: bool = False,
    records: Iterable[FeatureRecord] | None = None,
    child_id_separator: str = DEFAULT_CHILD_ID_SEPARATOR,
) -> None:
    """Write one model's real child features under an output-specific parent ID.

    Child IDs are ``<parent_id><child_id_separator><feature type><n>``.
    """
    feature_counts: Counter[str] = Counter()
    used_child_ids: set[str] = set()
    phase_by_part = emitted_cds_phase_by_part(model)
    for record in _ordered_child_records(
        model, supported_only=supported_only, records=records
    ):
        feature_counts[record.feature_type] += 1
        feature_name = sanitize_id_part(record.feature_type)
        child_id_base = (
            f"{parent_id}{child_id_separator}{feature_name}"
            f"{feature_counts[record.feature_type]}"
        )
        child_id = child_id_base
        collision_index = 2
        while child_id in used_child_ids:
            child_id = f"{child_id_base}_{collision_index}"
            collision_index += 1
        used_child_ids.add(child_id)
        attrs: dict[str, object] = {
            "ID": child_id,
            "Parent": parent_id,
        }
        if record.id:
            attrs["original_id"] = record.id
        parent_text = record.attrs.get("Parent") or record.attrs.get("transcript_id")
        if parent_text:
            attrs["original_parent"] = parent_text
        out.write(feature_line(
            seqid=record.seqid,
            source=record.source,
            feature_type=record.feature_type,
            start0=record.start0,
            end0=record.end0,
            score=record.score,
            strand=record.strand,
            phase=(
                phase_by_part.get((record.start0, record.end0), record.phase)
                if record.feature_type.lower() == "cds"
                else record.phase
            ),
            attrs=attrs,
        ) + "\n")


@dataclass(frozen=True, slots=True)
class _StandardizedTranscriptGeometry:
    start0: int
    end0: int
    child_records: tuple[FeatureRecord, ...]


def inherited_utr_validations(
    annotation: FinalAnnotation,
    *,
    minimum_intron_bp: int,
    fasta: IndexedFasta | None = None,
) -> dict[str, dict[Side, UTRSideValidation]]:
    """Assess inherited UTR sides after the coding population is frozen."""
    results: dict[str, dict[Side, UTRSideValidation]] = {}
    for gene in annotation.genes:
        for transcript in gene.transcripts:
            model = transcript.source_model
            source_exons = tuple(sorted(
                (record.start0, record.end0)
                for record in model.child_records
                if record.feature_type.lower() == "exon"
            ))
            explicit_utr_parts = [
                (record.start0, record.end0)
                for record in model.child_records
                if record.feature_type.lower() in {
                    "utr", "five_prime_utr", "three_prime_utr"
                }
            ]
            exon_parts = source_exons or tuple(
                merge_intervals(list(model.cds_parts) + explicit_utr_parts)
            )
            geometry_valid = bool(exon_parts) and all(
                start0 >= 0 and end0 > start0
                for start0, end0 in exon_parts
            ) and all(
                left[1] <= right[0]
                for left, right in zip(exon_parts, exon_parts[1:])
            )
            by_side: dict[Side, UTRSideValidation] = {}
            for side in ("five_prime", "three_prime"):
                if not side_has_utr(exon_parts, model.cds_parts, model.strand, side):
                    continue
                by_side[side] = validate_utr_side(
                    transcript_id=transcript.transcript_id,
                    evidence_origin="inherited_source_model",
                    evidence_id=model.output_transcript_id,
                    chrom=model.chrom,
                    strand=model.strand,
                    side=side,
                    introns=utr_side_introns(
                        exon_parts, model.cds_parts, model.strand, side
                    ),
                    minimum_intron_bp=minimum_intron_bp,
                    fasta=fasta,
                    geometry_valid=geometry_valid,
                    observation_mode="exceptions",
                )
            results[transcript.transcript_id] = by_side
    return results


def _filter_rejected_utr_records(
    model: CodingModel,
    records: Iterable[FeatureRecord],
    validations: Mapping[Side, UTRSideValidation],
) -> list[FeatureRecord]:
    rejected = {
        side for side, result in validations.items() if result.status == "rejected"
    }
    if not rejected:
        return list(records)
    reject_low = (
        (model.strand == "+" and "five_prime" in rejected)
        or (model.strand == "-" and "three_prime" in rejected)
    )
    reject_high = (
        (model.strand == "+" and "three_prime" in rejected)
        or (model.strand == "-" and "five_prime" in rejected)
    )
    kept: list[FeatureRecord] = []
    for record in records:
        feature_type = record.feature_type.lower()
        if feature_type not in {"exon", "utr", "five_prime_utr", "three_prime_utr"}:
            kept.append(record)
            continue
        if (
            (feature_type == "five_prime_utr" and "five_prime" in rejected)
            or (feature_type == "three_prime_utr" and "three_prime" in rejected)
        ):
            continue
        start0 = max(record.start0, model.cds_start0) if reject_low else record.start0
        end0 = min(record.end0, model.cds_end0) if reject_high else record.end0
        if end0 <= start0:
            continue
        kept.append(
            record if (start0, end0) == (record.start0, record.end0)
            else replace(
                record,
                start=start0 + 1,
                end=end0,
                raw_line="",
            )
        )
    return kept


def _standardized_transcript_geometry(
    model: CodingModel,
    validations: Mapping[Side, UTRSideValidation] | None = None,
) -> _StandardizedTranscriptGeometry:
    """Build output bounds solely from supported emitted structure."""
    metrics = model.sequence_integrity
    records = [
        record
        for record in model.child_records
        if record.feature_type.lower() in SUPPORTED_FINAL_CHILD_TYPES
        and not (
            record.feature_type.lower() == "start_codon"
            and metrics is not None
            and metrics.start_codon_status == "missing"
        )
        and not (
            record.feature_type.lower() == "stop_codon"
            and metrics is not None
            and metrics.stop_codon_status == "missing"
        )
    ]
    if validations:
        records = _filter_rejected_utr_records(model, records, validations)
    source_exons = [
        (record.start0, record.end0)
        for record in records
        if record.feature_type.lower() == "exon"
    ]
    emitted_structure = [
        (record.start0, record.end0)
        for record in records
        if record.feature_type.lower() in {
            "cds",
            "utr",
            "five_prime_utr",
            "three_prime_utr",
            "start_codon",
            "stop_codon",
        }
    ]
    exon_parts = tuple(
        merge_intervals(source_exons + emitted_structure + list(model.cds_parts))
    )
    if not exon_parts:
        raise ValueError(
            f"Final transcript {model.output_transcript_id!r} has no recognized "
            "structural geometry"
        )
    exact_source_exons = {
        (record.start0, record.end0): record
        for record in records
        if record.feature_type.lower() == "exon"
    }
    template = next(
        (
            record
            for record in records
            if record.feature_type.lower()
            in {"exon", "cds", "utr", "five_prime_utr", "three_prime_utr"}
        ),
        None,
    )
    if template is None:
        raise ValueError(
            f"Final transcript {model.output_transcript_id!r} has no record template"
        )
    standardized_exons: list[FeatureRecord] = []
    for part in exon_parts:
        source_exon = exact_source_exons.get(part)
        standardized_exons.append(
            source_exon
            if source_exon is not None
            else replace(
                template,
                feature_type="exon",
                start=part[0] + 1,
                end=part[1],
                score=".",
                phase=".",
                attrs={
                    key: value for key, value in template.attrs.items() if key != "ID"
                },
                raw_line="",
            )
        )
    non_exons = [
        record for record in records if record.feature_type.lower() != "exon"
    ]
    return _StandardizedTranscriptGeometry(
        start0=min(start0 for start0, _end0 in exon_parts),
        end0=max(end0 for _start0, end0 in exon_parts),
        child_records=tuple((*standardized_exons, *non_exons)),
    )


def write_cds_only_transcript_children(
    out: TextIO,
    model: CodingModel,
    parent_id: str,
    *,
    child_id_separator: str = DEFAULT_CHILD_ID_SEPARATOR,
) -> None:
    """Write exon/CDS pairs reconstructed strictly from one model's CDS parts.

    Child IDs are ``<parent_id><child_id_separator>cds_exon<n>`` and
    ``<parent_id><child_id_separator>CDS<n>``.
    """
    model_source = f"{GFF3_SOURCE}:{model.source_name}"
    phase_by_part = emitted_cds_phase_by_part(model)
    for index, (start0, end0) in enumerate(sorted(model.cds_parts), start=1):
        common_attrs: dict[str, object] = {
            "Parent": parent_id,
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
            attrs={
                "ID": f"{parent_id}{child_id_separator}cds_exon{index}",
                **common_attrs,
            },
        ) + "\n")
        out.write(feature_line(
            seqid=model.chrom,
            source=model_source,
            feature_type="CDS",
            start0=start0,
            end0=end0,
            score=".",
            strand=model.strand,
            phase=phase_by_part[(start0, end0)],
            attrs={
                "ID": f"{parent_id}{child_id_separator}CDS{index}",
                **common_attrs,
            },
        ) + "\n")


def public_transcript_attributes(
    gene: FinalGene,
    transcript: FinalTranscript,
    *,
    include_alternative_metadata: bool,
    cds_only_structure: bool = False,
) -> dict[str, object]:
    """Serialize immutable public-transcript metadata for all GFF3 writers."""
    model = transcript.source_model
    transcript_attrs: dict[str, object] = {
        "ID": transcript.transcript_id,
        "Parent": gene.gene_id,
        "moselecta_transcript_role": transcript.transcript_role,
        "original_source": model.source_name,
        "original_gene_id": model.input_gene_id,
        "original_transcript_id": model.input_transcript_id,
        "moselecta_internal_gene_id": model.output_gene_id,
        "moselecta_internal_transcript_id": model.output_transcript_id,
    }
    if cds_only_structure:
        transcript_attrs["final_gff3_cds_only_structure"] = "true"
    summary = transcript.annotation_summary
    for key, value in (
        (
            "model_selection_score_scaled",
            summary.model_selection_score_scaled,
        ),
        (
            "model_structural_support_score_scaled",
            summary.model_structural_support_score_scaled,
        ),
        ("protein_evidence_score_scaled", summary.protein_evidence_score_scaled),
        ("transcript_evidence_score_scaled", summary.transcript_evidence_score_scaled),
    ):
        if value is not None:
            transcript_attrs[key] = fmt_optional_float(value, 6)
    selected_sequence_status = summary.sequence_integrity_status
    current_sequence_status = (
        model.sequence_integrity.status
        if model.sequence_integrity is not None
        else selected_sequence_status
    )
    if current_sequence_status is not None:
        transcript_attrs["sequence_integrity_status"] = current_sequence_status
    if model.cds_terminus_refined is not None:
        transcript_attrs["cds_terminus_refined"] = model.cds_terminus_refined
    recovery = model.cds_phase_recovery
    if recovery is not None and recovery.status == RECOVERED_MISSING_SOURCE_PHASE:
        transcript_attrs["cds_phase_recovered"] = recovery.recovered_first_phase
    transcript_attrs["cds_phase_output_status"] = cds_phase_output_status(model)
    boundary_provenance = compact_five_prime_boundary_provenance(model)
    if boundary_provenance is not None:
        transcript_attrs["cds_5prime_boundary_canonicalised"] = boundary_provenance
    canonicalisation_provenance = compact_terminal_stop_provenance(model)
    if canonicalisation_provenance is not None:
        transcript_attrs["terminal_stop_canonicalised"] = canonicalisation_provenance
    if (
        selected_sequence_status is not None
        and selected_sequence_status != current_sequence_status
    ):
        transcript_attrs[
            "selected_sequence_integrity_status"
        ] = selected_sequence_status
    if include_alternative_metadata and transcript.alternative_metadata is not None:
        metadata = transcript.alternative_metadata
        transcript_attrs.update({
            "moselecta_alt_candidate_order_rank": metadata.candidate_order_rank,
            "moselecta_alt_preset": metadata.preset,
            "moselecta_alt_support_policy": metadata.support_policy,
            "moselecta_alt_evidence_requirement": metadata.evidence_requirement,
            "moselecta_alt_support_routes": (
                ",".join(metadata.support_routes) if metadata.support_routes else "none"
            ),
            "moselecta_alt_mosaic_policy_required": str(
                metadata.mosaic_policy_required
            ).lower(),
            "moselecta_alt_reciprocal_cds_overlap": fmt_optional_float(
                metadata.reciprocal_cds_overlap, 6
            ),
            "moselecta_alt_frame_concordance": fmt_optional_float(
                metadata.frame_concordance, 6
            ),
            "moselecta_alt_primary_score_ratio": fmt_optional_float(
                metadata.primary_score_ratio, 6
            ),
        })
    return transcript_attrs


def _write_final_transcript(
    out: TextIO,
    gene: FinalGene,
    transcript: FinalTranscript,
    *,
    include_alternative_metadata: bool,
    cds_only_structure: bool,
    geometry: _StandardizedTranscriptGeometry | None,
    utr_validations: Mapping[Side, UTRSideValidation] | None = None,
) -> None:
    model = transcript.source_model
    if not cds_only_structure and geometry is None:
        raise ValueError("Standardized transcript geometry is required")
    transcript_start0 = (
        model.cds_start0 if cds_only_structure else geometry.start0
    )
    transcript_end0 = model.cds_end0 if cds_only_structure else geometry.end0
    model_source = f"{GFF3_SOURCE}:{model.source_name}"
    transcript_attrs = public_transcript_attributes(
        gene,
        transcript,
        include_alternative_metadata=include_alternative_metadata,
        cds_only_structure=cds_only_structure,
    )
    for side, result in sorted((utr_validations or {}).items()):
        transcript_attrs[f"moselecta_{side}_utr_validation"] = result.status
        if result.reason_codes:
            transcript_attrs[f"moselecta_{side}_utr_validation_reasons"] = ",".join(
                result.reason_codes
            )
    out.write(feature_line(
        seqid=model.chrom,
        source=model_source,
        feature_type="mRNA",
        start0=transcript_start0,
        end0=transcript_end0,
        score=model.transcript_record.score if model.transcript_record else ".",
        strand=model.strand,
        phase=".",
        attrs=transcript_attrs,
    ) + "\n")

    if cds_only_structure:
        write_cds_only_transcript_children(
            out,
            model,
            transcript.transcript_id,
            child_id_separator=PUBLIC_CHILD_ID_SEPARATOR,
        )
    else:
        write_source_model_children(
            out,
            model,
            transcript.transcript_id,
            supported_only=True,
            records=geometry.child_records,
            child_id_separator=PUBLIC_CHILD_ID_SEPARATOR,
        )


def alternative_splice_configuration_comments(
    config: AlternativeSpliceSelectionConfig,
    *,
    output_description: str,
) -> tuple[str, ...]:
    """Return safe ordinary GFF comments describing selector configuration."""
    thresholds = config.effective_thresholds
    overrides = ",".join(config.overridden_thresholds) or "none"
    return (
        f"# MoSelecta alternative-splice output; {output_description}",
        f"# moselecta_alternative_splice_preset_requested={config.requested_preset.value}",
        f"# moselecta_alternative_splice_support_policy={config.support_policy.value}",
        "# moselecta_alternative_splice_evidence_requirement="
        f"{config.evidence_requirement.value}",
        "# moselecta_alternative_splice_prediction_source_count="
        f"{config.prediction_source_count}",
        "# moselecta_alternative_splice_min_model_support_sources="
        f"{config.min_model_support_sources}",
        "# moselecta_alternative_splice_min_model_support_source_fraction="
        f"{config.min_model_support_source_fraction if config.min_model_support_source_fraction is not None else '.'}",
        "# moselecta_alternative_splice_effective_min_model_support_sources="
        f"{config.effective_min_model_support_sources}",
        "# moselecta_alternative_splice_effective_min_reciprocal_cds_overlap="
        f"{thresholds.min_reciprocal_cds_overlap:.6f}",
        "# moselecta_alternative_splice_effective_min_frame_concordance="
        f"{thresholds.min_frame_concordance:.6f}",
        "# moselecta_alternative_splice_effective_min_model_selection_score_scaled="
        f"{thresholds.min_model_selection_score_scaled:.6f}",
        "# moselecta_alternative_splice_effective_min_primary_score_ratio="
        f"{thresholds.min_primary_score_ratio:.6f}",
        "# moselecta_alternative_splice_max_alternatives_per_gene="
        f"{config.max_alternatives_per_gene}",
        f"# moselecta_alternative_splice_threshold_overrides={overrides}",
    )


def _write_final_annotation_gff3(
    annotation: FinalAnnotation,
    path: Path,
    *,
    header_comments: Iterable[str] = (),
    include_alternative_metadata: bool = False,
    cds_only_structure: bool = False,
    utr_validations: Mapping[str, Mapping[Side, UTRSideValidation]] | None = None,
) -> None:
    with path.open("w") as out:
        out.write("##gff-version 3\n")
        for comment in header_comments:
            out.write(comment.rstrip("\n") + "\n")
        for gene in annotation.genes:
            transcript_geometry = {
                transcript.transcript_id: (
                    None
                    if cds_only_structure
                    else _standardized_transcript_geometry(
                        transcript.source_model,
                        (utr_validations or {}).get(transcript.transcript_id),
                    )
                )
                for transcript in gene.transcripts
            }
            gene_attrs: dict[str, object] = {
                "ID": gene.gene_id,
                "moselecta_primary_transcript_id": gene.primary_transcript.transcript_id,
                "moselecta_transcript_count": len(gene.transcripts),
            }
            if cds_only_structure:
                gene_attrs["final_gff3_cds_only_structure"] = "true"
            out.write(feature_line(
                seqid=gene.chrom,
                source=GFF3_SOURCE,
                feature_type="gene",
                start0=(
                    min(
                        transcript.source_model.cds_start0
                        for transcript in gene.transcripts
                    )
                    if cds_only_structure
                    else min(
                        geometry.start0
                        for geometry in transcript_geometry.values()
                        if geometry is not None
                    )
                ),
                end0=(
                    max(
                        transcript.source_model.cds_end0
                        for transcript in gene.transcripts
                    )
                    if cds_only_structure
                    else max(
                        geometry.end0
                        for geometry in transcript_geometry.values()
                        if geometry is not None
                    )
                ),
                score=".",
                strand=gene.strand,
                phase=".",
                attrs=gene_attrs,
            ) + "\n")
            for transcript in gene.transcripts:
                _write_final_transcript(
                    out,
                    gene,
                    transcript,
                    include_alternative_metadata=include_alternative_metadata,
                    cds_only_structure=cds_only_structure,
                    geometry=transcript_geometry[transcript.transcript_id],
                    utr_validations=(utr_validations or {}).get(
                        transcript.transcript_id
                    ),
                )
            out.write("###\n")


def write_standardised_primary_gff3(
    annotation: FinalAnnotation,
    path: Path,
    *,
    cds_only_structure: bool = False,
) -> None:
    """Write the primary-only annotation using standardised IDs.

    The implementation is intentionally multi-transcript-ready even though
    Primary-only construction supplies exactly one transcript for every gene.
    """
    _write_final_annotation_gff3(
        annotation,
        path,
        cds_only_structure=cds_only_structure,
    )


def write_final_annotation_gff3(
    annotation: FinalAnnotation,
    path: Path,
    metadata: FinalAnnotationOutputMetadata,
    *,
    cds_only_structure: bool = False,
    utr_validations: Mapping[str, Mapping[Side, UTRSideValidation]] | None = None,
) -> None:
    """Write the final annotation."""
    if utr_validations is None and not cds_only_structure:
        utr_validations = inherited_utr_validations(
            annotation, minimum_intron_bp=20
        )
    missing_metadata = [
        transcript.transcript_id
        for gene in annotation.genes
        for transcript in gene.transcripts
        if transcript.transcript_role == "alternative"
        and transcript.alternative_metadata is None
    ]
    if missing_metadata:
        raise ValueError(
            "Final alternatives are missing selector output metadata: "
            + ", ".join(missing_metadata)
        )
    state_comments = (
        "# MoSelecta final annotation",
        "# moselecta_alternative_splice_enabled="
        f"{str(metadata.alternative_splice_enabled).lower()}",
        "# moselecta_alternative_splice_selector_evaluated="
        f"{str(metadata.selector_evaluated).lower()}",
        "# moselecta_alternative_splice_alternatives_emitted="
        f"{str(metadata.alternatives_emitted).lower()}",
    )
    configuration_comments = (
        alternative_splice_configuration_comments(
            metadata.selection_config,
            output_description="alternative-splice selector configuration",
        )
        if metadata.selection_config is not None
        else ()
    )
    _write_final_annotation_gff3(
        annotation,
        path,
        header_comments=state_comments + configuration_comments,
        include_alternative_metadata=True,
        cds_only_structure=cds_only_structure,
        utr_validations=utr_validations,
    )
