"""Final annotation representation for retained MoSelecta models.

The selection, scoring, evidence, and pruning pipeline continues to operate on
``CodingModel`` objects and their internal output identifiers.  This module is
the separate post-selection layer that assigns stable user-facing identifiers
to the frozen retained-primary population.

Public gene numbering is deterministic and depends only on retained primaries,
ordered by chromosome, CDS start, CDS end, then internal transcript ID.

Public identifiers share one prefix: gene ``MoSelecta_000001`` owns primary
transcript ``MoSelecta_000001.1`` and any accepted alternatives
``MoSelecta_000001.2``, ``MoSelecta_000001.3``, and so on.
:func:`public_gene_id` and :func:`public_transcript_id` are the single places
that build these IDs; nothing else should rebuild the strings by hand.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

from .models import CodingModel

if TYPE_CHECKING:
    from .alternative_splice_selection import AlternativeSpliceSelectionResult

PUBLIC_GENE_ID_PREFIX = "MoSelecta_"
PUBLIC_ID_NUMBER_WIDTH = 6


def public_gene_id(index: int) -> str:
    """Return the public gene ID for 1-based gene ``index``.

    >>> public_gene_id(1)
    'MoSelecta_000001'
    """
    if index < 1:
        raise ValueError(f"Public gene index must be >= 1, got {index}")
    return f"{PUBLIC_GENE_ID_PREFIX}{index:0{PUBLIC_ID_NUMBER_WIDTH}d}"


def public_transcript_id(gene_id: str, ordinal: int) -> str:
    """Return the public transcript ID for ``ordinal`` within ``gene_id``.

    The transcript ID is the gene ID followed by ``.`` and the 1-based
    transcript ordinal (``.1`` is always the primary).

    >>> public_transcript_id("MoSelecta_000001", 1)
    'MoSelecta_000001.1'
    >>> public_transcript_id("MoSelecta_000001", 2)
    'MoSelecta_000001.2'
    """
    if not gene_id.startswith(PUBLIC_GENE_ID_PREFIX):
        raise ValueError(
            f"Public gene ID {gene_id!r} does not start with {PUBLIC_GENE_ID_PREFIX!r}"
        )
    if ordinal < 1:
        raise ValueError(f"Public transcript ordinal must be >= 1, got {ordinal}")
    return f"{gene_id}.{ordinal}"


@dataclass(frozen=True, slots=True)
class AnnotationSummary:
    """Compact transcript-specific support metadata for public annotation."""

    model_selection_score_scaled: float | None = None
    model_structural_support_score_scaled: float | None = None
    protein_evidence_score_scaled: float | None = None
    transcript_evidence_score_scaled: float | None = None
    sequence_integrity_status: str | None = None


@dataclass(frozen=True, slots=True)
class AlternativeTranscriptMetadata:
    """Concise retention facts for an accepted alternative."""

    candidate_order_rank: int
    preset: str
    support_policy: str
    evidence_requirement: str
    support_routes: tuple[str, ...]
    mosaic_policy_required: bool
    reciprocal_cds_overlap: float
    frame_concordance: float
    primary_score_ratio: float


@dataclass(frozen=True, slots=True)
class FinalTranscript:
    """One public transcript backed by one actual selected source model."""

    transcript_id: str
    source_model: CodingModel
    transcript_role: Literal["primary", "alternative"] = "primary"
    annotation_summary: AnnotationSummary = field(default_factory=AnnotationSummary)
    alternative_metadata: AlternativeTranscriptMetadata | None = None


@dataclass(frozen=True, slots=True)
class FinalGene:
    """A public MoSelecta gene containing one or more transcript models."""

    gene_id: str
    transcripts: tuple[FinalTranscript, ...]

    def __post_init__(self) -> None:
        if not self.transcripts:
            raise ValueError(f"Final gene {self.gene_id!r} must contain at least one transcript")
        expected_primary_id = public_transcript_id(self.gene_id, 1)
        if self.transcripts[0].transcript_id != expected_primary_id:
            raise ValueError(
                f"Final gene {self.gene_id!r} must begin with primary transcript "
                f"{expected_primary_id!r}"
            )
        if self.transcripts[0].transcript_role != "primary":
            raise ValueError(
                f"Final gene {self.gene_id!r} first transcript must have primary role"
            )
        transcript_ids = [transcript.transcript_id for transcript in self.transcripts]
        if len(transcript_ids) != len(set(transcript_ids)):
            raise ValueError(f"Final gene {self.gene_id!r} contains duplicate transcript IDs")
        for index, transcript in enumerate(self.transcripts, start=1):
            expected_transcript_id = public_transcript_id(self.gene_id, index)
            if transcript.transcript_id != expected_transcript_id:
                raise ValueError(
                    f"Final gene {self.gene_id!r} transcript {index} must use public ID "
                    f"{expected_transcript_id!r}"
                )
            expected_role = "primary" if index == 1 else "alternative"
            if transcript.transcript_role != expected_role:
                raise ValueError(
                    f"Final gene {self.gene_id!r} transcript {index} must have "
                    f"{expected_role!r} role"
                )
        internal_ids = [
            transcript.source_model.output_transcript_id
            for transcript in self.transcripts
        ]
        if len(internal_ids) != len(set(internal_ids)):
            raise ValueError(f"Final gene {self.gene_id!r} contains duplicate source models")
        primary_model = self.transcripts[0].source_model
        for transcript in self.transcripts[1:]:
            model = transcript.source_model
            if model.chrom != primary_model.chrom or model.strand != primary_model.strand:
                raise ValueError(
                    f"Final gene {self.gene_id!r} contains transcripts on different "
                    "chromosomes or strands"
                )

    @property
    def primary_transcript(self) -> FinalTranscript:
        return self.transcripts[0]

    @property
    def chrom(self) -> str:
        return self.primary_transcript.source_model.chrom

    @property
    def strand(self) -> str:
        return self.primary_transcript.source_model.strand

    @property
    def start0(self) -> int:
        return min(transcript.source_model.model_start0 for transcript in self.transcripts)

    @property
    def end0(self) -> int:
        return max(transcript.source_model.model_end0 for transcript in self.transcripts)


@dataclass(frozen=True, slots=True)
class FinalAnnotation:
    """Ordered final annotation independent of internal ``CodingModel`` IDs."""

    genes: tuple[FinalGene, ...]

    def __post_init__(self) -> None:
        gene_ids = [gene.gene_id for gene in self.genes]
        if len(gene_ids) != len(set(gene_ids)):
            raise ValueError("Final annotation contains duplicate gene IDs")
        transcript_ids = [
            transcript.transcript_id
            for gene in self.genes
            for transcript in gene.transcripts
        ]
        if len(transcript_ids) != len(set(transcript_ids)):
            raise ValueError("Final annotation contains duplicate transcript IDs")

    @property
    def transcript_count(self) -> int:
        return sum(len(gene.transcripts) for gene in self.genes)

    @property
    def gene_count(self) -> int:
        return len(self.genes)

    @property
    def primary_transcript_count(self) -> int:
        return len(self.genes)

    @property
    def alternative_transcript_count(self) -> int:
        return self.transcript_count - self.primary_transcript_count

    def gene_by_primary_internal_transcript_id(self) -> dict[str, FinalGene]:
        """Map frozen primary internal transcript IDs to their public genes."""
        return {
            gene.primary_transcript.source_model.output_transcript_id: gene
            for gene in self.genes
        }


def final_primary_sort_key(model: CodingModel) -> tuple[str, int, int, str]:
    """Return the stable ordering key used for public primary gene numbering."""
    return (
        model.chrom,
        model.cds_start0,
        model.cds_end0,
        model.output_transcript_id,
    )


def _optional_finite_float(value: object) -> float | None:
    if value is None or value == ".":
        return None
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return None
    return numeric_value if math.isfinite(numeric_value) else None


def build_annotation_summary(
    model: CodingModel,
    model_score_row: Mapping[str, object] | None,
) -> AnnotationSummary:
    """Build compact public metadata using existing score availability semantics."""
    row = model_score_row or {}

    protein_status = str(
        row.get("protein_evidence_status", "protein_evidence_not_configured")
    )
    protein_score = (
        None
        if protein_status == "protein_evidence_not_configured"
        else _optional_finite_float(row.get("protein_evidence_score_scaled"))
    )

    transcript_status = str(
        row.get("transcript_evidence_status", "transcript_evidence_not_configured")
    )
    transcript_score = (
        None
        if transcript_status == "transcript_evidence_not_configured"
        else _optional_finite_float(row.get("transcript_evidence_score_scaled"))
    )

    sequence_status_value = row.get("sequence_integrity_status")
    if sequence_status_value is None or sequence_status_value == ".":
        sequence_status = (
            model.sequence_integrity.status
            if model.sequence_integrity is not None
            else None
        )
    else:
        sequence_status = str(sequence_status_value)
    if sequence_status == "genome_not_configured":
        sequence_status = None

    return AnnotationSummary(
        model_selection_score_scaled=_optional_finite_float(
            row.get("model_selection_score_scaled")
        ),
        model_structural_support_score_scaled=_optional_finite_float(
            row.get("model_structural_support_score_scaled")
        ),
        protein_evidence_score_scaled=protein_score,
        transcript_evidence_score_scaled=transcript_score,
        sequence_integrity_status=sequence_status,
    )


def build_primary_annotation(
    retained_primary_models: Iterable[CodingModel],
    model_score_rows_by_transcript: Mapping[str, Mapping[str, object]] | None = None,
) -> FinalAnnotation:
    """Build a primary-only public annotation from the frozen retained set.

    The function never mutates ``CodingModel`` objects.  Public IDs are derived
    exclusively from the supplied retained primaries and are consequently
    independent of any future alternative-splice candidate population.
    """
    retained_models = list(retained_primary_models)
    score_rows = (
        model_score_rows_by_transcript
        if model_score_rows_by_transcript is not None
        else {}
    )
    internal_ids = [model.output_transcript_id for model in retained_models]
    if len(internal_ids) != len(set(internal_ids)):
        raise ValueError("Frozen retained-primary set contains duplicate internal transcript IDs")
    if model_score_rows_by_transcript is not None:
        missing_score_ids = sorted(set(internal_ids) - set(score_rows))
        if missing_score_ids:
            raise ValueError(
                "Frozen retained-primary models are missing canonical score rows: "
                + ", ".join(missing_score_ids)
            )

    genes: list[FinalGene] = []
    for index, model in enumerate(sorted(retained_models, key=final_primary_sort_key), start=1):
        gene_id = public_gene_id(index)
        primary = FinalTranscript(
            transcript_id=public_transcript_id(gene_id, 1),
            source_model=model,
            transcript_role="primary",
            annotation_summary=build_annotation_summary(
                model,
                score_rows.get(model.output_transcript_id),
            ),
        )
        genes.append(FinalGene(gene_id=gene_id, transcripts=(primary,)))
    return FinalAnnotation(genes=tuple(genes))


def augment_annotation(
    primary_annotation: FinalAnnotation,
    selection_result: AlternativeSpliceSelectionResult,
    model_score_rows_by_transcript: Mapping[str, Mapping[str, object]],
) -> FinalAnnotation:
    """Add Stage 2B accepted representatives without changing public primaries.

    The returned immutable annotation is produced by the alternative-splice
    selector for serialization in the final public annotation. Alternative order and public
    ``.2+`` IDs are taken directly from the active greedy result so the writer
    cannot silently renumber them.
    """
    accepted_by_primary = selection_result.accepted_by_primary_internal_id()
    genes: list[FinalGene] = []
    for gene in primary_annotation.genes:
        primary = gene.primary_transcript
        primary_internal_id = primary.source_model.output_transcript_id
        accepted = accepted_by_primary.get(primary_internal_id, ())
        transcripts: list[FinalTranscript] = [primary]
        for item in accepted:
            model = item.candidate_model
            if item.public_gene_id != gene.gene_id:
                raise ValueError(
                    f"Accepted alternative {model.output_transcript_id!r} targets "
                    f"{item.public_gene_id!r}, expected {gene.gene_id!r}"
                )
            expected_id = public_transcript_id(gene.gene_id, len(transcripts) + 1)
            if item.accepted_public_transcript_id != expected_id:
                raise ValueError(
                    f"Accepted alternative {model.output_transcript_id!r} has public ID "
                    f"{item.accepted_public_transcript_id!r}, expected {expected_id!r}"
                )
            score_row = model_score_rows_by_transcript.get(model.output_transcript_id)
            if score_row is None:
                raise ValueError(
                    "Accepted alternative is missing its canonical score row: "
                    + model.output_transcript_id
                )
            transcripts.append(FinalTranscript(
                transcript_id=item.accepted_public_transcript_id,
                source_model=model,
                transcript_role="alternative",
                annotation_summary=build_annotation_summary(model, score_row),
                alternative_metadata=AlternativeTranscriptMetadata(
                    candidate_order_rank=item.alternative_candidate_order_rank,
                    preset=selection_result.config.requested_preset.value,
                    support_policy=selection_result.config.support_policy.value,
                    evidence_requirement=(
                        selection_result.config.evidence_requirement.value
                    ),
                    support_routes=item.factual_support_routes,
                    mosaic_policy_required=item.mosaic_policy_required,
                    reciprocal_cds_overlap=item.stage2b_reciprocal_cds_overlap,
                    frame_concordance=item.stage2b_frame_concordance,
                    primary_score_ratio=item.stage2b_primary_score_ratio,
                ),
            ))
        genes.append(FinalGene(gene_id=gene.gene_id, transcripts=tuple(transcripts)))
    return FinalAnnotation(genes=tuple(genes))
