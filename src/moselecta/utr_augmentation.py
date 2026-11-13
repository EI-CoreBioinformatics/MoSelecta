"""Post-selection UTR evidence, selection, reconstruction, and audit output.

This module is deliberately isolated from CDS scoring and selection.  It takes
an immutable :class:`~moselecta.final_annotation.FinalAnnotation`, treats every
selected CDS as authoritative, and returns a separate augmented representation.
"""

from __future__ import annotations

import math
import statistics
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Literal, Mapping, Sequence

from .cds_phase import emitted_cds_phase_by_part
from .constants import CHILD_FEATURE_ORDER, GFF3_SOURCE
from .annotation_outputs import (
    PUBLIC_CHILD_ID_SEPARATOR,
    public_transcript_attributes,
)
from .final_annotation import FinalAnnotation, FinalGene, FinalTranscript
from .intervals import merge_intervals
from .models import CodingModel
from .outputs import feature_line
from .transcript_evidence import TranscriptEvidenceRawTranscript
from .tsvio import write_dict_tsv
from .genome_sequence import IndexedFasta
from .utr_validation import UTRSideValidation, validate_utr_side


Interval = tuple[int, int]
StructuralFeature = tuple[str, str, str, Interval | int]
Side = Literal["five_prime", "three_prime"]
CandidateMode = Literal["utr", "completion"]
SIDES: tuple[Side, Side] = ("five_prime", "three_prime")
SINGLE_EXON_INDEX_BIN_SIZE = 10_000


@dataclass(frozen=True, slots=True)
class UTREvidence:
    evidence_id: str
    kind: Literal["transcript", "prediction"]
    chrom: str
    strand: str
    exon_parts: tuple[Interval, ...]
    source_member_counts: tuple[tuple[str, int], ...] = ()
    representative_ids: tuple[str, ...] = ()
    prediction_source: str = "."
    prediction_model_id: str = "."
    prediction_gene_id: str = "."
    prediction_cds_parts: tuple[Interval, ...] = ()
    prediction_selected: bool = False

    @property
    def start0(self) -> int:
        return self.exon_parts[0][0]

    @property
    def end0(self) -> int:
        return self.exon_parts[-1][1]

    @property
    def introns(self) -> tuple[Interval, ...]:
        return tuple(
            (left[1], right[0])
            for left, right in zip(self.exon_parts, self.exon_parts[1:])
        )

    @property
    def source_counts(self) -> dict[str, int]:
        return dict(self.source_member_counts)

    @property
    def sources(self) -> tuple[str, ...]:
        if self.kind == "prediction":
            return (self.prediction_source,)
        return tuple(source for source, count in self.source_member_counts if count > 0)


@dataclass(frozen=True, slots=True)
class UTRSettings:
    min_complete_exon_bp: int = 5
    min_intron_bp: int | None = 20
    max_intron_bp: float | None = None
    max_five_prime_exons: int | None = 4
    max_three_prime_exons: int | None = 2
    exonic_failsafe_bp: int | None = 20_000
    five_prime_exonic_override_bp: int | None = None
    three_prime_exonic_override_bp: int | None = None
    five_prime_genomic_override_bp: int | None = None
    three_prime_genomic_override_bp: int | None = None
    member_depth_cap: int = 3
    endpoint_window_min_bp: int = 20
    endpoint_window_max_bp: int = 250
    local_consensus_min_fraction: float = 0.5
    limit_min_observations: int = 1000

    def __post_init__(self) -> None:
        for name in (
            "min_complete_exon_bp", "member_depth_cap", "endpoint_window_min_bp",
            "endpoint_window_max_bp", "limit_min_observations",
        ):
            if int(getattr(self, name)) < 0:
                raise ValueError(f"{name} must be >= 0")
        if self.endpoint_window_max_bp < self.endpoint_window_min_bp:
            raise ValueError("endpoint_window_max_bp must be >= endpoint_window_min_bp")
        if not 0.0 <= self.local_consensus_min_fraction <= 1.0:
            raise ValueError("local_consensus_min_fraction must be between 0 and 1")


@dataclass(frozen=True, slots=True)
class UTRLearnedLimits:
    five_prime_exonic_bp: int | None
    three_prime_exonic_bp: int | None
    five_prime_genomic_bp: int | None
    three_prime_genomic_bp: int | None
    five_prime_exonic_source: str
    three_prime_exonic_source: str
    five_prime_genomic_source: str
    three_prime_genomic_source: str
    warnings: tuple[str, ...] = ()

    def value(self, side: Side, geometry: str) -> int | None:
        return getattr(self, f"{side}_{geometry}_bp")

    def source(self, side: Side, geometry: str) -> str:
        return getattr(self, f"{side}_{geometry}_source")


@dataclass(slots=True)
class UTRCandidateSide:
    gene_id: str
    transcript_id: str
    target: FinalTranscript
    evidence: UTREvidence
    compatibility_class: str
    side: Side
    endpoint0: int
    extension_interval: Interval
    added_exons: tuple[Interval, ...]
    added_introns: tuple[Interval, ...]
    complete_contributing_exon_lengths: tuple[int, ...]
    extension_exonic_bp: int
    extension_genomic_bp: int
    structure_key: str
    source_cds_class: str
    source_cds_overlap_bp: int
    validation_status: str = "passed"
    validation_reason_codes: tuple[str, ...] = ()
    rejection_reasons: tuple[str, ...] = ()
    multi_gene_bridge: bool = False

    @property
    def eligible_nonextension(self) -> bool:
        return not self.rejection_reasons and not self.multi_gene_bridge


@dataclass(slots=True)
class PooledSupport:
    transcript_members: Counter[str] = field(default_factory=Counter)
    prediction_models: dict[str, set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )
    evidence_ids: set[str] = field(default_factory=set)
    evidence_kinds: set[str] = field(default_factory=set)
    prediction_cds_classes: set[str] = field(default_factory=set)
    prediction_selection_classes: set[str] = field(default_factory=set)

    def add_candidate(self, candidate: UTRCandidateSide) -> None:
        evidence = candidate.evidence
        if evidence.kind == "transcript":
            self.transcript_members.update(evidence.source_counts)
        else:
            self.prediction_models[evidence.prediction_source].add(
                evidence.prediction_model_id
            )
            self.prediction_cds_classes.add(candidate.source_cds_class)
            self.prediction_selection_classes.add(
                "selected_model" if evidence.prediction_selected else "unselected_model"
            )
        self.evidence_ids.add(evidence.evidence_id)
        self.evidence_kinds.add(evidence.kind)

    def merge(self, other: "PooledSupport") -> None:
        self.transcript_members.update(other.transcript_members)
        for source, models in other.prediction_models.items():
            self.prediction_models[source].update(models)
        self.evidence_ids.update(other.evidence_ids)
        self.evidence_kinds.update(other.evidence_kinds)
        self.prediction_cds_classes.update(other.prediction_cds_classes)
        self.prediction_selection_classes.update(other.prediction_selection_classes)

    @property
    def source_count(self) -> int:
        return len(self.transcript_members) + len(self.prediction_models)

    @property
    def multi_member_transcript_source_count(self) -> int:
        return sum(count >= 2 for count in self.transcript_members.values())

    def effective_support(self, cap: int = 3) -> float:
        return float(len(self.prediction_models)) + sum(
            transcript_source_weight(count, cap)
            for count in self.transcript_members.values()
        )


@dataclass(slots=True)
class UTRStructure:
    key: str
    candidates: list[UTRCandidateSide] = field(default_factory=list)
    by_endpoint: dict[int, list[UTRCandidateSide]] = field(
        default_factory=lambda: defaultdict(list)
    )

    def add(self, candidate: UTRCandidateSide) -> None:
        self.candidates.append(candidate)
        self.by_endpoint[candidate.endpoint0].append(candidate)

    def pooled(self, endpoints: Iterable[int] | None = None) -> PooledSupport:
        selected = set(self.by_endpoint) if endpoints is None else set(endpoints)
        pooled = PooledSupport()
        for endpoint in sorted(selected):
            for candidate in self.by_endpoint[endpoint]:
                pooled.add_candidate(candidate)
        return pooled


@dataclass(frozen=True, slots=True)
class UTRSideSelection:
    side: Side
    candidate: UTRCandidateSide
    endpoint_method: str
    effective_support: float
    independent_source_count: int
    multi_member_transcript_source_count: int
    eligible_structure_count: int
    strong_competing_structure: bool
    strong_multimodality: bool
    local_support_fraction: float
    adaptive_window_bp: int
    selected_evidence_ids: tuple[str, ...]
    selected_sources: tuple[str, ...]
    evidence_class: str
    prediction_cds_classes: tuple[str, ...]
    prediction_selection_classes: tuple[str, ...]
    quality_status: str = "selected_recheck_passed"


@dataclass(frozen=True, slots=True)
class AugmentedTranscript:
    final_transcript: FinalTranscript
    gene_id: str
    exon_parts: tuple[Interval, ...]
    five_prime_utr_parts: tuple[Interval, ...]
    three_prime_utr_parts: tuple[Interval, ...]
    five_prime_selection: UTRSideSelection | None
    three_prime_selection: UTRSideSelection | None
    choices_share_evidence: bool

    @property
    def transcript_id(self) -> str:
        return self.final_transcript.transcript_id

    @property
    def model(self) -> CodingModel:
        return self.final_transcript.source_model

    @property
    def start0(self) -> int:
        return self.exon_parts[0][0]

    @property
    def end0(self) -> int:
        return self.exon_parts[-1][1]


@dataclass(frozen=True, slots=True)
class AugmentedGene:
    final_gene: FinalGene
    transcripts: tuple[AugmentedTranscript, ...]

    @property
    def start0(self) -> int:
        return min(transcript.start0 for transcript in self.transcripts)

    @property
    def end0(self) -> int:
        return max(transcript.end0 for transcript in self.transcripts)


@dataclass(frozen=True, slots=True)
class AugmentedAnnotation:
    genes: tuple[AugmentedGene, ...]

    @property
    def transcript_count(self) -> int:
        return sum(len(gene.transcripts) for gene in self.genes)


@dataclass(frozen=True, slots=True)
class UTRAuditSummary:
    transcript_count: int
    both_assigned: int
    five_prime_only: int
    three_prime_only: int
    neither_assigned: int
    transcript_only: int
    prediction_only: int
    both_evidence_kinds: int
    multi_gene_evidence_excluded: int
    authoritative_candidate_sides_rejected: int
    provisional_selected_sides_in_utr_conflicts: int
    conflict_fallback_sides: int
    conflict_no_alternative_sides: int
    conflict_affected_genes: int
    conflict_affected_transcripts: int
    conflict_resolution_rounds: int
    local_consensus_sides: int
    diffuse_sides: int
    strong_competing_structure_sides: int
    strongly_multimodal_sides: int
    transcript_evidence_records: int
    prediction_evidence_records: int
    catalogue_estimated_bytes: int
    runtime_seconds: float
    rejection_counts: tuple[tuple[str, int], ...]


@dataclass(frozen=True, slots=True)
class UTRAugmentationResult:
    annotation: AugmentedAnnotation
    limits: UTRLearnedLimits
    audit_rows: tuple[dict[str, object], ...]
    summary: UTRAuditSummary


def transcript_source_weight(member_count: int, cap: int = 3) -> float:
    return 1.0 + min(float(cap), math.log2(max(1, member_count)))


def _interval_overlap(a: Interval, b: Interval) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def _contains(outer: Interval, inner: Interval) -> bool:
    return outer[0] <= inner[0] and inner[1] <= outer[1]


def _introns(exons: Sequence[Interval]) -> tuple[Interval, ...]:
    return tuple((left[1], right[0]) for left, right in zip(exons, exons[1:]))


def _structural_features(
    chrom: str, strand: str, exons: Sequence[Interval]
) -> set[StructuralFeature]:
    """Exact exon, intron, and generic splice-site keys used by final genes.

    MoSelecta joins same-strand multi-exon models on *any* exact splice-site
    coordinate, independent of donor/acceptor role.  Exact exons protect
    terminal and single-exon structures that have no applicable splice site.
    """
    ordered = tuple(sorted(exons))
    keys: set[StructuralFeature] = {
        (chrom, strand, "exon", exon) for exon in ordered
    }
    for intron in _introns(ordered):
        keys.add((chrom, strand, "intron", intron))
        keys.add((chrom, strand, "splice", intron[0]))
        keys.add((chrom, strand, "splice", intron[1]))
    return keys


def _materialized_exons(
    model: CodingModel, *side_exon_groups: Sequence[Interval]
) -> tuple[Interval, ...]:
    return tuple(merge_intervals(
        list(model.cds_parts)
        + [exon for group in side_exon_groups for exon in group]
    ))


def _bins_for_interval(interval: Interval) -> range:
    start0, end0 = interval
    if end0 <= start0:
        return range(0)
    return range(
        start0 // SINGLE_EXON_INDEX_BIN_SIZE,
        ((end0 - 1) // SINGLE_EXON_INDEX_BIN_SIZE) + 1,
    )


def _prediction_exon_parts(model: CodingModel) -> tuple[Interval, ...]:
    exon_parts = sorted({
        (record.start0, record.end0)
        for record in model.child_records
        if record.feature_type.lower() == "exon" and record.end0 > record.start0
    })
    if exon_parts:
        return tuple(exon_parts)
    parts = list(model.cds_parts)
    parts.extend(
        (record.start0, record.end0)
        for record in model.child_records
        if record.feature_type.lower() in {
            "five_prime_utr", "three_prime_utr", "utr",
        }
        and record.end0 > record.start0
    )
    return tuple(merge_intervals(parts))


class FrozenFinalGeneStructureIndex:
    """Exact structural features of frozen final genes, never optional UTRs.

    The index includes authoritative CDS-derived structure and the ordinary
    final annotation's source-exon structure.  It is gene-owned, not
    transcript-owned: sharing among public alternatives of one gene is valid.
    """

    def __init__(self, annotation: FinalAnnotation) -> None:
        self.owners: dict[StructuralFeature, dict[str, str]] = defaultdict(dict)
        self.baseline_by_transcript: dict[str, set[StructuralFeature]] = {}
        for gene in annotation.genes:
            for transcript in gene.transcripts:
                model = transcript.source_model
                cds_keys = _structural_features(
                    model.chrom, model.strand, model.cds_parts
                )
                source_keys = _structural_features(
                    model.chrom, model.strand, _prediction_exon_parts(model)
                )
                self.baseline_by_transcript[transcript.transcript_id] = (
                    cds_keys | source_keys
                )
                for key in source_keys:
                    self.owners[key].setdefault(gene.gene_id, "ordinary_source")
                for key in cds_keys:
                    self.owners[key][gene.gene_id] = "authoritative"

    def candidate_reasons(self, candidate: UTRCandidateSide) -> tuple[str, ...]:
        model = candidate.target.source_model
        features = _structural_features(
            model.chrom, model.strand,
            _materialized_exons(model, candidate.added_exons),
        )
        new_features = features - self.baseline_by_transcript[candidate.transcript_id]
        reasons: set[str] = set()
        for key in new_features:
            for owner_gene, owner_kind in self.owners.get(key, {}).items():
                if owner_gene != candidate.gene_id:
                    prefix = (
                        "cross_final_gene_authoritative"
                        if owner_kind == "authoritative"
                        else "cross_final_gene_ordinary_source"
                    )
                    label = "splice" if key[2] == "splice" else key[2]
                    reasons.add(f"{prefix}_{label}_conflict")
        return tuple(sorted(reasons))


def _apply_frozen_gene_filter(
    candidates: Iterable[UTRCandidateSide],
    frozen: FrozenFinalGeneStructureIndex,
) -> None:
    for candidate in candidates:
        conflicts = frozen.candidate_reasons(candidate)
        if conflicts:
            candidate.rejection_reasons = tuple(sorted(
                set(candidate.rejection_reasons) | set(conflicts)
            ))


def build_transcript_utr_catalogue(
    raw_transcripts: Iterable[TranscriptEvidenceRawTranscript],
) -> list[UTREvidence]:
    """Collapse raw transcripts by exact endpoint-paired exon chain."""
    groups: dict[
        tuple[str, str, tuple[Interval, ...]],
        tuple[Counter[str], list[str]],
    ] = {}
    for transcript in raw_transcripts:
        if transcript.strand not in {"+", "-"}:
            continue
        key = (transcript.chrom, transcript.strand, transcript.exon_parts)
        if key not in groups:
            groups[key] = (Counter(), [])
        counts, examples = groups[key]
        counts[transcript.source_label] += 1
        if len(examples) < 3:
            examples.append(transcript.raw_id)
    catalogue: list[UTREvidence] = []
    for index, ((chrom, strand, exons), (counts, examples)) in enumerate(
        groups.items(), start=1
    ):
        catalogue.append(UTREvidence(
            evidence_id=f"UTRTX_{index:07d}",
            kind="transcript",
            chrom=chrom,
            strand=strand,
            exon_parts=exons,
            source_member_counts=tuple(sorted(counts.items())),
            representative_ids=tuple(examples),
        ))
    return catalogue


def build_prediction_utr_catalogue(
    models: Iterable[CodingModel],
    *,
    selected_internal_ids: set[str],
) -> list[UTREvidence]:
    catalogue: list[UTREvidence] = []
    for model in models:
        exons = _prediction_exon_parts(model)
        if not exons or model.strand not in {"+", "-"}:
            continue
        catalogue.append(UTREvidence(
            evidence_id=f"UTRPRED:{model.source_name}:{model.output_transcript_id}",
            kind="prediction",
            chrom=model.chrom,
            strand=model.strand,
            exon_parts=exons,
            representative_ids=(model.input_transcript_id,),
            prediction_source=model.source_name,
            prediction_model_id=model.output_transcript_id,
            prediction_gene_id=model.input_gene_id,
            prediction_cds_parts=tuple(model.cds_parts),
            prediction_selected=model.output_transcript_id in selected_internal_ids,
        ))
    return catalogue


class UTREvidenceIndex:
    """Local structural indexes for the three accepted compatibility classes."""

    def __init__(
        self,
        evidence: Iterable[UTREvidence],
        targets: FinalAnnotation | None = None,
    ) -> None:
        self.full: dict[tuple[str, str, Interval], list[tuple[UTREvidence, int]]] = defaultdict(list)
        self.terminal: dict[tuple[str, str, Side, int], list[tuple[UTREvidence, int]]] = defaultdict(list)
        self.single: dict[tuple[str, str, int], list[tuple[UTREvidence, int]]] = defaultdict(list)
        wanted_full: set[tuple[str, str, Interval]] | None = None
        wanted_terminal: set[tuple[str, str, Side, int]] | None = None
        wanted_single: set[tuple[str, str, int]] | None = None
        if targets is not None:
            wanted_full = set()
            wanted_terminal = set()
            wanted_single = set()
            for gene in targets.genes:
                for transcript in gene.transcripts:
                    model = transcript.source_model
                    if model.cds_introns:
                        wanted_full.add(
                            (model.chrom, model.strand, model.cds_introns[0])
                        )
                        for side in SIDES:
                            info = _terminal_cds_info(model, side)
                            if info is not None:
                                wanted_terminal.add(
                                    (model.chrom, model.strand, side, info[2])
                                )
                    else:
                        wanted_single.add((
                            model.chrom,
                            model.strand,
                            model.cds_start0 // SINGLE_EXON_INDEX_BIN_SIZE,
                        ))
        for item in evidence:
            introns = item.introns
            for offset, intron in enumerate(introns):
                key = (item.chrom, item.strand, intron)
                if wanted_full is None or key in wanted_full:
                    self.full[key].append((item, offset))
            if len(item.exon_parts) >= 2:
                last = len(item.exon_parts) - 1
                for exon_index, exon in enumerate(item.exon_parts):
                    for side in SIDES:
                        site: int | None = None
                        if side == "five_prime":
                            if item.strand == "+" and exon_index < last:
                                site = exon[1]
                            elif item.strand == "-" and exon_index > 0:
                                site = exon[0]
                        else:
                            if item.strand == "+" and exon_index > 0:
                                site = exon[0]
                            elif item.strand == "-" and exon_index < last:
                                site = exon[1]
                        if site is not None:
                            key = (item.chrom, item.strand, side, site)
                            if wanted_terminal is None or key in wanted_terminal:
                                self.terminal[key].append((item, exon_index))
            for exon_index, exon in enumerate(item.exon_parts):
                for bin_id in _bins_for_interval(exon):
                    key = (item.chrom, item.strand, bin_id)
                    if wanted_single is None or key in wanted_single:
                        self.single[key].append((item, exon_index))


def _terminal_cds_info(model: CodingModel, side: Side) -> tuple[int, Interval, int] | None:
    parts = model.sorted_cds_parts
    if len(parts) < 2:
        return None
    if side == "five_prime":
        index = 0 if model.strand == "+" else len(parts) - 1
    else:
        index = len(parts) - 1 if model.strand == "+" else 0
    exon = parts[index]
    internal_site = (
        exon[1]
        if (side == "five_prime" and model.strand == "+")
        or (side == "three_prime" and model.strand == "-")
        else exon[0]
    )
    return index, exon, internal_site


def _incomplete_cds_utr_block_reason(model: CodingModel, side: Side) -> str | None:
    """Return the explicit incomplete-CDS reason that invalidates a UTR anchor."""
    metrics = model.sequence_integrity
    if metrics is None:
        return None
    if side == "five_prime" and metrics.start_codon_status == "missing":
        return "missing_start_codon"
    if side == "three_prime" and metrics.stop_codon_status == "missing":
        return "missing_stop_codon"
    return None


def _extension_interval(model: CodingModel, evidence: UTREvidence, side: Side) -> Interval | None:
    if side == "five_prime":
        if model.strand == "+" and evidence.start0 < model.cds_start0:
            return evidence.start0, model.cds_start0
        if model.strand == "-" and evidence.end0 > model.cds_end0:
            return model.cds_end0, evidence.end0
    else:
        if model.strand == "+" and evidence.end0 > model.cds_end0:
            return model.cds_end0, evidence.end0
        if model.strand == "-" and evidence.start0 < model.cds_start0:
            return evidence.start0, model.cds_start0
    return None


def _endpoint(evidence: UTREvidence, side: Side) -> int:
    low = (side == "five_prime" and evidence.strand == "+") or (
        side == "three_prime" and evidence.strand == "-"
    )
    return evidence.start0 if low else evidence.end0


def _structure_key(added_exons: Sequence[Interval], endpoint0: int) -> str:
    normalized = list(added_exons)
    if endpoint0 == normalized[0][0]:
        normalized[0] = (-1, normalized[0][1])
    else:
        normalized[-1] = (normalized[-1][0], -1)
    return ",".join(f"{start0}-{end0}" for start0, end0 in normalized)


def _source_cds_provenance(
    evidence: UTREvidence, added_exons: Sequence[Interval]
) -> tuple[str, int]:
    if evidence.kind != "prediction":
        return ".", 0
    overlap = sum(
        _interval_overlap(added, cds)
        for added in added_exons
        for cds in evidence.prediction_cds_parts
    )
    total = sum(end0 - start0 for start0, end0 in added_exons)
    if overlap == 0:
        return "outside_source_cds", 0
    if overlap >= total:
        return "within_source_cds", overlap
    return "mixed_source_cds_and_outside", overlap


def _nonextension_reasons(
    *,
    complete_lengths: Sequence[int],
    added_introns: Sequence[Interval],
    added_exon_count: int,
    side: Side,
    settings: UTRSettings,
    apply_utr_exon_limit: bool = True,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if settings.min_complete_exon_bp > 0 and complete_lengths:
        if complete_lengths[0] < settings.min_complete_exon_bp:
            reasons.append("tiny_terminal_complete_evidence_exon")
        if len(complete_lengths) > 1 and min(complete_lengths[1:]) < settings.min_complete_exon_bp:
            reasons.append("tiny_internal_complete_evidence_exon")
    intron_lengths = [end0 - start0 for start0, end0 in added_introns]
    if settings.max_intron_bp is not None and settings.max_intron_bp > 0:
        if intron_lengths and max(intron_lengths) > settings.max_intron_bp:
            reasons.append("long_utr_intron_moselecta")
    if apply_utr_exon_limit:
        maximum = (
            settings.max_five_prime_exons
            if side == "five_prime" else settings.max_three_prime_exons
        )
        if maximum is not None and maximum > 0 and added_exon_count > maximum:
            reasons.append(f"excessive_{side}_utr_exon_count")
    return tuple(sorted(set(reasons)))


def _make_candidate(
    gene_id: str,
    transcript: FinalTranscript,
    evidence: UTREvidence,
    compatibility_class: str,
    side: Side,
    settings: UTRSettings,
    *,
    candidate_mode: CandidateMode = "utr",
    fasta: IndexedFasta | None = None,
    validation_recorder: Callable[[UTRSideValidation], None] | None = None,
) -> UTRCandidateSide | None:
    model = transcript.source_model
    extension = _extension_interval(model, evidence, side)
    if extension is None:
        return None
    contributing: list[tuple[int, Interval, Interval]] = []
    for index, exon in enumerate(evidence.exon_parts):
        start0 = max(exon[0], extension[0])
        end0 = min(exon[1], extension[1])
        if end0 > start0:
            contributing.append((index, exon, (start0, end0)))
    if not contributing:
        return None
    endpoint0 = _endpoint(evidence, side)
    low = endpoint0 == contributing[0][2][0]
    ordered_for_filter = contributing if low else list(reversed(contributing))
    complete_lengths = tuple(
        exon[1] - exon[0] for _index, exon, _clipped in ordered_for_filter
    )
    added_exons = tuple(item[2] for item in contributing)
    added_introns = tuple(
        intron for intron in evidence.introns
        if extension[0] <= intron[0] and intron[1] <= extension[1]
    )
    source_cds_class, source_cds_overlap = _source_cds_provenance(
        evidence, added_exons
    )
    validation = validate_utr_side(
        transcript_id=transcript.transcript_id,
        evidence_origin=evidence.kind,
        evidence_id=evidence.evidence_id,
        chrom=evidence.chrom,
        strand=evidence.strand,
        side=side,
        introns=added_introns,
        minimum_intron_bp=settings.min_intron_bp or 0,
        fasta=fasta,
        observation_mode="exceptions" if validation_recorder else "none",
    )
    if validation_recorder is not None:
        validation_recorder(validation)
    validation_rejections = {
        reason for reason in validation.reason_codes
        if reason in {
            "invalid_utr_geometry",
            "short_utr_intron",
            "resolved_unrecognized_motif",
        }
    }
    reasons = tuple(sorted(set(_nonextension_reasons(
        complete_lengths=complete_lengths,
        added_introns=added_introns,
        added_exon_count=len(added_exons),
        side=side,
        settings=settings,
        apply_utr_exon_limit=candidate_mode == "utr",
    )) | validation_rejections))
    return UTRCandidateSide(
        gene_id=gene_id,
        transcript_id=transcript.transcript_id,
        target=transcript,
        evidence=evidence,
        compatibility_class=compatibility_class,
        side=side,
        endpoint0=endpoint0,
        extension_interval=extension,
        added_exons=added_exons,
        added_introns=added_introns,
        complete_contributing_exon_lengths=complete_lengths,
        extension_exonic_bp=sum(end0 - start0 for start0, end0 in added_exons),
        extension_genomic_bp=extension[1] - extension[0],
        structure_key=_structure_key(added_exons, endpoint0),
        source_cds_class=source_cds_class,
        source_cds_overlap_bp=source_cds_overlap,
        validation_status=validation.status,
        validation_reason_codes=validation.reason_codes,
        rejection_reasons=reasons,
    )


def candidate_sides_for_transcript(
    gene_id: str,
    transcript: FinalTranscript,
    index: UTREvidenceIndex,
    settings: UTRSettings,
    *,
    candidate_mode: CandidateMode = "utr",
    requested_sides: Sequence[Side] = SIDES,
    fasta: IndexedFasta | None = None,
    validation_recorder: Callable[[UTRSideValidation], None] | None = None,
) -> list[UTRCandidateSide]:
    model = transcript.source_model
    matches: list[tuple[int, UTREvidence, str, tuple[Side, ...]]] = []
    seen_full: set[tuple[str, int]] = set()
    model_introns = tuple(model.cds_introns)
    if model_introns:
        for evidence, offset in index.full.get(
            (model.chrom, model.strand, model_introns[0]), []
        ):
            key = (evidence.evidence_id, offset)
            if key in seen_full:
                continue
            seen_full.add(key)
            evidence_introns = evidence.introns
            if evidence_introns[offset:offset + len(model_introns)] != model_introns:
                continue
            evidence_exons = evidence.exon_parts[offset:offset + len(model.cds_parts)]
            if len(evidence_exons) != len(model.cds_parts):
                continue
            if not all(
                _contains(evidence_exon, cds_exon)
                for evidence_exon, cds_exon in zip(evidence_exons, model.cds_parts)
            ):
                continue
            matches.append((0, evidence, "full_containment", tuple(requested_sides)))

    for side in requested_sides:
        info = _terminal_cds_info(model, side)
        if info is None:
            continue
        _model_index, model_exon, internal_site = info
        for evidence, evidence_exon_index in index.terminal.get(
            (model.chrom, model.strand, side, internal_site), []
        ):
            if _contains(evidence.exon_parts[evidence_exon_index], model_exon):
                matches.append((2, evidence, "terminal_compatibility", (side,)))

    if len(model.cds_parts) == 1:
        model_exon = model.cds_parts[0]
        bin_id = model_exon[0] // SINGLE_EXON_INDEX_BIN_SIZE
        seen_single: set[tuple[str, int]] = set()
        for evidence, evidence_exon_index in index.single.get(
            (model.chrom, model.strand, bin_id), []
        ):
            key = (evidence.evidence_id, evidence_exon_index)
            if key in seen_single:
                continue
            seen_single.add(key)
            evidence_exon = evidence.exon_parts[evidence_exon_index]
            if not _contains(evidence_exon, model_exon):
                continue
            extending_sides = tuple(
                side
                for side in requested_sides
                if _extension_interval(model, evidence, side) is not None
            )
            if candidate_mode == "utr":
                if evidence_exon[0] < model_exon[0] and model_exon[1] < evidence_exon[1]:
                    matches.append((1, evidence, "single_exon_containment", extending_sides))
            elif extending_sides:
                matches.append((1, evidence, "single_exon_one_sided_containment", extending_sides))

    candidates: list[UTRCandidateSide] = []
    seen_side: set[tuple[str, Side]] = set()
    for _rank, evidence, mode, sides in sorted(
        matches, key=lambda item: (item[0], item[1].evidence_id)
    ):
        for side in sides:
            if (
                candidate_mode == "utr"
                and _incomplete_cds_utr_block_reason(model, side) is not None
            ):
                continue
            key = (evidence.evidence_id, side)
            if key in seen_side:
                continue
            candidate = _make_candidate(
                gene_id,
                transcript,
                evidence,
                mode,
                side,
                settings,
                candidate_mode=candidate_mode,
                fasta=fasta,
                validation_recorder=validation_recorder,
            )
            if candidate is not None:
                candidates.append(candidate)
                seen_side.add(key)
    return candidates


def mark_multi_gene_bridges(candidates: Iterable[UTRCandidateSide]) -> int:
    genes_by_evidence: dict[str, set[str]] = defaultdict(set)
    candidate_list = list(candidates)
    for candidate in candidate_list:
        if candidate.validation_status != "rejected":
            genes_by_evidence[candidate.evidence.evidence_id].add(candidate.gene_id)
    bridged = {
        evidence_id for evidence_id, genes in genes_by_evidence.items()
        if len(genes) > 1
    }
    for candidate in candidate_list:
        candidate.multi_gene_bridge = candidate.evidence.evidence_id in bridged
    return len(bridged)


def _quantile(values: Sequence[int], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    weight = position - lower
    return float(ordered[lower] * (1.0 - weight) + ordered[upper] * weight)


def _structure_source_count(candidates: Sequence[UTRCandidateSide]) -> int:
    pooled = PooledSupport()
    for candidate in candidates:
        pooled.add_candidate(candidate)
    return pooled.source_count


def learn_extension_limits(
    candidates: Sequence[UTRCandidateSide], settings: UTRSettings
) -> UTRLearnedLimits:
    admissible = [candidate for candidate in candidates if candidate.eligible_nonextension]
    by_structure: dict[tuple[str, Side, str], list[UTRCandidateSide]] = defaultdict(list)
    for candidate in admissible:
        by_structure[
            (candidate.transcript_id, candidate.side, candidate.structure_key)
        ].append(candidate)
    strong_keys = {
        key for key, values in by_structure.items()
        if _structure_source_count(values) >= 2
    }
    values: dict[tuple[Side, str], list[int]] = defaultdict(list)
    strong_values: dict[tuple[Side, str], list[int]] = defaultdict(list)
    for candidate in admissible:
        key = (candidate.transcript_id, candidate.side, candidate.structure_key)
        for geometry, value in (
            ("exonic", candidate.extension_exonic_bp),
            ("genomic", candidate.extension_genomic_bp),
        ):
            values[(candidate.side, geometry)].append(value)
            if key in strong_keys:
                strong_values[(candidate.side, geometry)].append(value)

    return _limits_from_observations(values, strong_values, settings)


def _limits_from_observations(
    values: Mapping[tuple[Side, str], Sequence[int]],
    strong_values: Mapping[tuple[Side, str], Sequence[int]],
    settings: UTRSettings,
) -> UTRLearnedLimits:
    """Resolve learned/override limits from compact run-level observations."""

    warnings: list[str] = []
    results: dict[tuple[Side, str], tuple[int | None, str]] = {}
    for side in SIDES:
        for geometry in ("exonic", "genomic"):
            override = getattr(settings, f"{side}_{geometry}_override_bp")
            if override is not None:
                results[(side, geometry)] = (
                    override if override > 0 else None,
                    "explicit_override" if override > 0 else "explicitly_disabled",
                )
                continue
            side_values = list(values.get((side, geometry), ()))
            learning_values = side_values
            population = "side_specific"
            if len(learning_values) < settings.limit_min_observations:
                learning_values = (
                    list(values.get(("five_prime", geometry), ()))
                    + list(values.get(("three_prime", geometry), ()))
                )
                population = "pooled_sides"
            if len(learning_values) < settings.limit_min_observations:
                if geometry == "exonic":
                    results[(side, geometry)] = (
                        settings.exonic_failsafe_bp,
                        "sparse_data_exonic_failsafe",
                    )
                else:
                    results[(side, geometry)] = (None, "sparse_data_unavailable")
                    warnings.append(
                        f"Could not learn {side} genomic UTR extension limit: "
                        f"only {len(learning_values)} pooled admissible observations"
                    )
                continue
            all_p995 = _quantile(learning_values, 0.995)
            side_strong = list(strong_values.get((side, geometry), ()))
            if population == "pooled_sides":
                side_strong = (
                    list(strong_values.get(("five_prime", geometry), ()))
                    + list(strong_values.get(("three_prime", geometry), ()))
                )
            tail_values = (
                side_strong
                if len(side_strong) >= settings.limit_min_observations
                else learning_values
            )
            raw = max(all_p995, _quantile(tail_values, 0.999))
            learned = int(math.ceil(raw / 1000.0) * 1000)
            if geometry == "exonic" and settings.exonic_failsafe_bp is not None:
                learned = min(learned, settings.exonic_failsafe_bp)
            results[(side, geometry)] = (learned, f"learned_{population}")

    return UTRLearnedLimits(
        five_prime_exonic_bp=results[("five_prime", "exonic")][0],
        three_prime_exonic_bp=results[("three_prime", "exonic")][0],
        five_prime_genomic_bp=results[("five_prime", "genomic")][0],
        three_prime_genomic_bp=results[("three_prime", "genomic")][0],
        five_prime_exonic_source=results[("five_prime", "exonic")][1],
        three_prime_exonic_source=results[("three_prime", "exonic")][1],
        five_prime_genomic_source=results[("five_prime", "genomic")][1],
        three_prime_genomic_source=results[("three_prime", "genomic")][1],
        warnings=tuple(warnings),
    )


def _multi_gene_bridge_ids(
    annotation: FinalAnnotation,
    index: UTREvidenceIndex,
    settings: UTRSettings,
    fasta: IndexedFasta | None = None,
) -> set[str]:
    """Return compatible evidence IDs associated with multiple final genes.

    Only a first owning gene is retained for the common case, so this pass does
    not allocate one set per evidence record.
    """
    first_gene: dict[str, str] = {}
    bridged: set[str] = set()
    for gene in annotation.genes:
        for transcript in gene.transcripts:
            for candidate in candidate_sides_for_transcript(
                gene.gene_id, transcript, index, settings, fasta=fasta
            ):
                if candidate.validation_status == "rejected":
                    continue
                evidence_id = candidate.evidence.evidence_id
                owner = first_gene.setdefault(evidence_id, gene.gene_id)
                if owner != gene.gene_id:
                    bridged.add(evidence_id)
    return bridged


def _streaming_extension_limits(
    annotation: FinalAnnotation,
    index: UTREvidenceIndex,
    settings: UTRSettings,
    bridged_evidence_ids: set[str],
    frozen: FrozenFinalGeneStructureIndex | None = None,
    excluded_training_sides: set[tuple[str, Side]] | None = None,
    fasta: IndexedFasta | None = None,
) -> UTRLearnedLimits:
    """Learn limits without retaining the run-wide candidate-side table."""
    values: dict[tuple[Side, str], list[int]] = defaultdict(list)
    strong_values: dict[tuple[Side, str], list[int]] = defaultdict(list)
    for gene in annotation.genes:
        for transcript in gene.transcripts:
            candidates = candidate_sides_for_transcript(
                gene.gene_id, transcript, index, settings, fasta=fasta
            )
            if excluded_training_sides:
                candidates = [
                    candidate
                    for candidate in candidates
                    if (transcript.transcript_id, candidate.side)
                    not in excluded_training_sides
                ]
            if frozen is not None:
                _apply_frozen_gene_filter(candidates, frozen)
            candidates = [candidate for candidate in candidates
                          if candidate.evidence.evidence_id not in bridged_evidence_ids
                          and not candidate.rejection_reasons]
            by_structure: dict[tuple[Side, str], list[UTRCandidateSide]] = defaultdict(list)
            for candidate in candidates:
                by_structure[(candidate.side, candidate.structure_key)].append(candidate)
            strong_keys = {
                key for key, grouped in by_structure.items()
                if _structure_source_count(grouped) >= 2
            }
            for candidate in candidates:
                key = (candidate.side, candidate.structure_key)
                for geometry, value in (
                    ("exonic", candidate.extension_exonic_bp),
                    ("genomic", candidate.extension_genomic_bp),
                ):
                    values[(candidate.side, geometry)].append(value)
                    if key in strong_keys:
                        strong_values[(candidate.side, geometry)].append(value)
    return _limits_from_observations(values, strong_values, settings)


def learn_utr_limits_for_annotation(
    annotation: FinalAnnotation,
    transcript_evidence: Sequence[UTREvidence],
    prediction_evidence: Sequence[UTREvidence],
    settings: UTRSettings,
    *,
    excluded_training_sides: set[tuple[str, Side]] | None = None,
    fasta: IndexedFasta | None = None,
) -> UTRLearnedLimits:
    """Learn ordinary UTR limits without retaining run-wide candidates."""
    index = UTREvidenceIndex(list(transcript_evidence) + list(prediction_evidence), annotation)
    frozen = FrozenFinalGeneStructureIndex(annotation)
    bridged = _multi_gene_bridge_ids(annotation, index, settings, fasta)
    return _streaming_extension_limits(
        annotation,
        index,
        settings,
        bridged,
        frozen,
        excluded_training_sides,
        fasta,
    )


def _extension_reasons(
    candidate: UTRCandidateSide, limits: UTRLearnedLimits
) -> tuple[str, ...]:
    reasons: list[str] = []
    exonic = limits.value(candidate.side, "exonic")
    genomic = limits.value(candidate.side, "genomic")
    if exonic is not None and candidate.extension_exonic_bp > exonic:
        reasons.append("excessive_exonic_extension")
    if genomic is not None and candidate.extension_genomic_bp > genomic:
        reasons.append("excessive_genomic_extension")
    return tuple(reasons)


def _structure_rank(structure: UTRStructure, cap: int) -> tuple[float, int, int]:
    pooled = structure.pooled()
    return (
        -pooled.effective_support(cap),
        -pooled.source_count,
        -pooled.multi_member_transcript_source_count,
    )


def _endpoint_weights(
    structure: UTRStructure, endpoints: Iterable[int], cap: int
) -> dict[int, float]:
    selected = sorted(set(endpoints))
    weights = {endpoint: 0.0 for endpoint in selected}
    transcript_sources = {
        source
        for endpoint in selected
        for candidate in structure.by_endpoint[endpoint]
        if candidate.evidence.kind == "transcript"
        for source in candidate.evidence.source_counts
    }
    for source in sorted(transcript_sources):
        counts: Counter[int] = Counter()
        for endpoint in selected:
            for candidate in structure.by_endpoint[endpoint]:
                if candidate.evidence.kind == "transcript":
                    counts[endpoint] += candidate.evidence.source_counts.get(source, 0)
        total = sum(counts.values())
        if total:
            allowed = transcript_source_weight(total, cap)
            for endpoint, count in counts.items():
                weights[endpoint] += allowed * count / total
    prediction_sources = {
        candidate.evidence.prediction_source
        for endpoint in selected
        for candidate in structure.by_endpoint[endpoint]
        if candidate.evidence.kind == "prediction"
    }
    for source in sorted(prediction_sources):
        counts: Counter[int] = Counter()
        for endpoint in selected:
            models = {
                candidate.evidence.prediction_model_id
                for candidate in structure.by_endpoint[endpoint]
                if candidate.evidence.kind == "prediction"
                and candidate.evidence.prediction_source == source
            }
            counts[endpoint] = len(models)
        total = sum(counts.values())
        if total:
            for endpoint, count in counts.items():
                weights[endpoint] += count / total
    return weights


def weighted_observed_medoid(weights: Mapping[int, float]) -> int:
    if not weights:
        raise ValueError("Cannot calculate a medoid without observed endpoints")
    return min(
        weights,
        key=lambda coordinate: (
            sum(weight * abs(coordinate - other) for other, weight in weights.items()),
            -weights[coordinate],
            coordinate,
        ),
    )


def _adaptive_window(endpoints: Sequence[int], settings: UTRSettings) -> int:
    gaps = [right - left for left, right in zip(endpoints, endpoints[1:]) if right > left]
    value = math.floor(3 * statistics.median(gaps)) if gaps else settings.endpoint_window_min_bp
    return max(settings.endpoint_window_min_bp, min(settings.endpoint_window_max_bp, value))


def _bounded_windows(endpoints: Sequence[int], width: int) -> list[list[int]]:
    windows: list[list[int]] = []
    right = 0
    for left in range(len(endpoints)):
        right = max(right, left)
        while right + 1 < len(endpoints) and endpoints[right + 1] - endpoints[left] <= width:
            right += 1
        windows.append(list(endpoints[left:right + 1]))
    return windows


def _window_rank(
    structure: UTRStructure, endpoints: list[int], cap: int
) -> tuple[object, ...]:
    pooled = structure.pooled(endpoints)
    medoid = weighted_observed_medoid(_endpoint_weights(structure, endpoints, cap))
    return (
        -pooled.effective_support(cap),
        -pooled.source_count,
        -pooled.multi_member_transcript_source_count,
        endpoints[-1] - endpoints[0],
        medoid,
        endpoints[0],
        endpoints[-1],
    )


@dataclass(frozen=True, slots=True)
class _EndpointDecision:
    endpoint0: int
    method: str
    local_support_fraction: float
    window_width: int
    strong_multimodality: bool


def _select_endpoint(
    structure: UTRStructure, settings: UTRSettings
) -> _EndpointDecision:
    endpoints = sorted(structure.by_endpoint)
    width = _adaptive_window(endpoints, settings)
    windows = _bounded_windows(endpoints, width)
    winner = min(windows, key=lambda values: _window_rank(structure, values, settings.member_depth_cap))
    total_support = structure.pooled().effective_support(settings.member_depth_cap)
    local_support = structure.pooled(winner).effective_support(settings.member_depth_cap)
    fraction = local_support / total_support if total_support else 0.0
    local = len(endpoints) == 1 or fraction > settings.local_consensus_min_fraction
    chosen = winner if local else endpoints
    endpoint0 = weighted_observed_medoid(
        _endpoint_weights(structure, chosen, settings.member_depth_cap)
    )
    disjoint = [
        values for values in windows
        if values[-1] < winner[0] or values[0] > winner[-1]
    ]
    second = (
        min(disjoint, key=lambda values: _window_rank(structure, values, settings.member_depth_cap))
        if disjoint else None
    )
    second_pooled = structure.pooled(second) if second else None
    multimodal = bool(
        second
        and second_pooled is not None
        and second_pooled.effective_support(settings.member_depth_cap) >= 0.5 * local_support
        and second_pooled.source_count >= 2
        and min(abs(left - right) for left in winner for right in second) > width
    )
    return _EndpointDecision(
        endpoint0=endpoint0,
        method="local_consensus_medoid" if local else "diffuse_global_medoid",
        local_support_fraction=fraction,
        window_width=width,
        strong_multimodality=multimodal,
    )


def select_completion_candidate(
    candidates: Sequence[UTRCandidateSide],
    settings: UTRSettings,
) -> UTRSideSelection | None:
    """Select one codon-blind completion structure by ordinary support rank.

    Completion candidates have already passed completion-specific structural
    eligibility.  UTR extension-length and UTR-exon-count checks are
    intentionally absent here because part of the apparent extension may
    become CDS.  A codon search is performed only after this function returns,
    and failure never advances to a lower-ranked structure.
    """
    eligible = [candidate for candidate in candidates if candidate.eligible_nonextension]
    ranked, strong_competing = _ranked_structures(eligible, settings)
    if not ranked:
        return None
    structure = ranked[0]
    decision = _select_endpoint(structure, settings)
    proposals = sorted(
        structure.by_endpoint[decision.endpoint0],
        key=lambda item: item.evidence.evidence_id,
    )
    if not proposals:
        return None
    candidate = proposals[0]
    pooled = structure.pooled()
    endpoint_pooled = structure.pooled((candidate.endpoint0,))
    sources = sorted(
        set(endpoint_pooled.transcript_members)
        | set(endpoint_pooled.prediction_models)
    )
    kinds = endpoint_pooled.evidence_kinds
    evidence_class = (
        "both"
        if kinds == {"transcript", "prediction"}
        else next(iter(kinds))
        if kinds
        else "none"
    )
    return UTRSideSelection(
        side=candidate.side,
        candidate=candidate,
        endpoint_method=decision.method,
        effective_support=pooled.effective_support(settings.member_depth_cap),
        independent_source_count=pooled.source_count,
        multi_member_transcript_source_count=pooled.multi_member_transcript_source_count,
        eligible_structure_count=len(ranked),
        strong_competing_structure=strong_competing,
        strong_multimodality=decision.strong_multimodality,
        local_support_fraction=decision.local_support_fraction,
        adaptive_window_bp=decision.window_width,
        selected_evidence_ids=tuple(sorted(endpoint_pooled.evidence_ids)),
        selected_sources=tuple(sources),
        evidence_class=evidence_class,
        prediction_cds_classes=tuple(sorted(endpoint_pooled.prediction_cds_classes)),
        prediction_selection_classes=tuple(
            sorted(endpoint_pooled.prediction_selection_classes)
        ),
        quality_status="completion_structure_selected",
    )


def _candidate_recheck(
    candidate: UTRCandidateSide,
    settings: UTRSettings,
    limits: UTRLearnedLimits,
) -> bool:
    return not _nonextension_reasons(
        complete_lengths=candidate.complete_contributing_exon_lengths,
        added_introns=candidate.added_introns,
        added_exon_count=len(candidate.added_exons),
        side=candidate.side,
        settings=settings,
    ) and not _extension_reasons(candidate, limits)


def _selection_for_structure(
    structure: UTRStructure,
    structure_count: int,
    strong_competing: bool,
    settings: UTRSettings,
    limits: UTRLearnedLimits,
    *,
    fallback_structure: bool = False,
) -> UTRSideSelection | None:
    decision = _select_endpoint(structure, settings)
    endpoint_order = [decision.endpoint0] + [
        endpoint for endpoint in sorted(structure.by_endpoint)
        if endpoint != decision.endpoint0
    ]
    candidate: UTRCandidateSide | None = None
    quality_status = (
        "fallback_structure_after_recheck"
        if fallback_structure else "selected_recheck_passed"
    )
    for endpoint in endpoint_order:
        for proposal in sorted(
            structure.by_endpoint[endpoint],
            key=lambda item: item.evidence.evidence_id,
        ):
            if _candidate_recheck(proposal, settings, limits):
                candidate = proposal
                if endpoint != decision.endpoint0:
                    quality_status = (
                        "fallback_structure_and_endpoint_after_recheck"
                        if fallback_structure
                        else "fallback_endpoint_after_recheck"
                    )
                break
        if candidate is not None:
            break
    if candidate is None:
        return None
    pooled = structure.pooled()
    endpoint_pooled = structure.pooled((candidate.endpoint0,))
    sources = sorted(
        set(endpoint_pooled.transcript_members) | set(endpoint_pooled.prediction_models)
    )
    kinds = endpoint_pooled.evidence_kinds
    evidence_class = (
        "both" if kinds == {"transcript", "prediction"}
        else next(iter(kinds)) if kinds else "none"
    )
    return UTRSideSelection(
        side=candidate.side,
        candidate=candidate,
        endpoint_method=decision.method,
        effective_support=pooled.effective_support(settings.member_depth_cap),
        independent_source_count=pooled.source_count,
        multi_member_transcript_source_count=pooled.multi_member_transcript_source_count,
        eligible_structure_count=structure_count,
        strong_competing_structure=strong_competing,
        strong_multimodality=decision.strong_multimodality,
        local_support_fraction=decision.local_support_fraction,
        adaptive_window_bp=decision.window_width,
        selected_evidence_ids=tuple(sorted(endpoint_pooled.evidence_ids)),
        selected_sources=tuple(sources),
        evidence_class=evidence_class,
        prediction_cds_classes=tuple(sorted(endpoint_pooled.prediction_cds_classes)),
        prediction_selection_classes=tuple(
            sorted(endpoint_pooled.prediction_selection_classes)
        ),
        quality_status=quality_status,
    )


def _ranked_structure_ties(
    candidates: Sequence[UTRCandidateSide], settings: UTRSettings
) -> tuple[list[UTRStructure], bool]:
    structures: dict[str, UTRStructure] = {}
    for candidate in candidates:
        structures.setdefault(candidate.structure_key, UTRStructure(candidate.structure_key)).add(candidate)
    ranked = sorted(
        structures.values(),
        key=lambda structure: (_structure_rank(structure, settings.member_depth_cap), structure.key),
    )
    if not ranked:
        return [], False
    top_numerical = _structure_rank(ranked[0], settings.member_depth_cap)
    tied = [
        structure for structure in ranked
        if _structure_rank(structure, settings.member_depth_cap) == top_numerical
    ]
    top_support = ranked[0].pooled().effective_support(settings.member_depth_cap)
    strong_competing = any(
        structure.pooled().source_count >= 2
        and structure.pooled().effective_support(settings.member_depth_cap) >= 0.5 * top_support
        for structure in ranked[1:]
    )
    return tied, strong_competing


def _ranked_structures(
    candidates: Sequence[UTRCandidateSide], settings: UTRSettings
) -> tuple[list[UTRStructure], bool]:
    structures: dict[str, UTRStructure] = {}
    for candidate in candidates:
        structures.setdefault(
            candidate.structure_key, UTRStructure(candidate.structure_key)
        ).add(candidate)
    ranked = sorted(
        structures.values(),
        key=lambda structure: (
            _structure_rank(structure, settings.member_depth_cap), structure.key
        ),
    )
    if not ranked:
        return [], False
    top_support = ranked[0].pooled().effective_support(settings.member_depth_cap)
    strong_competing = any(
        structure.pooled().source_count >= 2
        and structure.pooled().effective_support(settings.member_depth_cap)
        >= 0.5 * top_support
        for structure in ranked[1:]
    )
    return ranked, strong_competing


def _select_sides_for_transcript(
    candidates: Sequence[UTRCandidateSide],
    settings: UTRSettings,
    limits: UTRLearnedLimits,
    *,
    blocked_signatures: Mapping[Side, set[tuple[Interval, ...]]] | None = None,
) -> tuple[UTRSideSelection | None, UTRSideSelection | None]:
    options: dict[Side, list[UTRSideSelection]] = {side: [] for side in SIDES}
    for side in SIDES:
        side_candidates = [
            candidate for candidate in candidates
            if candidate.side == side
            and candidate.eligible_nonextension
            and not _extension_reasons(candidate, limits)
            and (
                blocked_signatures is None
                or candidate.added_exons not in blocked_signatures.get(side, set())
            )
        ]
        all_structure_count = len({candidate.structure_key for candidate in side_candidates})
        ranked, strong_competing = _ranked_structures(side_candidates, settings)
        rank_index = 0
        while rank_index < len(ranked) and not options[side]:
            group_start = rank_index
            numerical_rank = _structure_rank(
                ranked[rank_index], settings.member_depth_cap
            )
            tied: list[UTRStructure] = []
            while (
                rank_index < len(ranked)
                and _structure_rank(ranked[rank_index], settings.member_depth_cap)
                == numerical_rank
            ):
                tied.append(ranked[rank_index])
                rank_index += 1
            fallback_structure = group_start > 0
            for structure in tied:
                selection = _selection_for_structure(
                    structure,
                    all_structure_count,
                    strong_competing,
                    settings,
                    limits,
                    fallback_structure=fallback_structure,
                )
                if selection is not None:
                    options[side].append(selection)
        options[side].sort(
            key=lambda selection: (
                selection.candidate.structure_key,
                selection.candidate.endpoint0,
            )
        )

    five = options["five_prime"]
    three = options["three_prime"]
    if len(five) > 1 and len(three) > 1:
        five_choice, three_choice = min(
            ((left, right) for left in five for right in three),
            key=lambda pair: (
                -len(set(pair[0].selected_evidence_ids) & set(pair[1].selected_evidence_ids)),
                pair[0].candidate.structure_key,
                pair[1].candidate.structure_key,
                pair[0].candidate.endpoint0,
                pair[1].candidate.endpoint0,
            ),
        )
        return five_choice, three_choice
    return (five[0] if five else None, three[0] if three else None)


def _audit_side_fields(
    side: Side,
    selection: UTRSideSelection | None,
    candidate_count: int,
    limits: UTRLearnedLimits,
    anchor_block_reason: str | None = None,
) -> dict[str, object]:
    prefix = "five_prime" if side == "five_prime" else "three_prime"
    if selection is None:
        return {
            f"{prefix}_assigned": "false",
            f"{prefix}_selected_structure_key": ".",
            f"{prefix}_selected_endpoint0": ".",
            f"{prefix}_endpoint_method": ".",
            f"{prefix}_effective_support": ".",
            f"{prefix}_independent_source_count": 0,
            f"{prefix}_multi_member_transcript_source_count": 0,
            f"{prefix}_selected_evidence_ids": ".",
            f"{prefix}_selected_sources": ".",
            f"{prefix}_evidence_class": ".",
            f"{prefix}_prediction_cds_provenance": ".",
            f"{prefix}_prediction_selection_provenance": ".",
            f"{prefix}_eligible_competing_structure_count": 0,
            f"{prefix}_strong_competing_structure": "false",
            f"{prefix}_strong_endpoint_multimodality": "false",
            f"{prefix}_local_support_fraction": ".",
            f"{prefix}_adaptive_window_bp": ".",
            f"{prefix}_exonic_limit_bp": limits.value(side, "exonic") or ".",
            f"{prefix}_genomic_limit_bp": limits.value(side, "genomic") or ".",
            f"{prefix}_selected_exonic_extension_bp": 0,
            f"{prefix}_selected_genomic_extension_bp": 0,
            f"{prefix}_anchor_block_reason": anchor_block_reason or ".",
            f"{prefix}_quality_status": (
                "blocked_incomplete_cds_anchor"
                if anchor_block_reason else
                "no_eligible_proposal" if candidate_count else
                "no_compatible_evidence"
            ),
        }
    candidate = selection.candidate
    return {
        f"{prefix}_assigned": "true",
        f"{prefix}_selected_structure_key": candidate.structure_key,
        f"{prefix}_selected_endpoint0": candidate.endpoint0,
        f"{prefix}_endpoint_method": selection.endpoint_method,
        f"{prefix}_effective_support": f"{selection.effective_support:.6f}",
        f"{prefix}_independent_source_count": selection.independent_source_count,
        f"{prefix}_multi_member_transcript_source_count": selection.multi_member_transcript_source_count,
        f"{prefix}_selected_evidence_ids": ",".join(selection.selected_evidence_ids) or ".",
        f"{prefix}_selected_sources": ",".join(selection.selected_sources) or ".",
        f"{prefix}_evidence_class": selection.evidence_class,
        f"{prefix}_prediction_cds_provenance": ",".join(selection.prediction_cds_classes) or ".",
        f"{prefix}_prediction_selection_provenance": ",".join(
            selection.prediction_selection_classes
        ) or ".",
        f"{prefix}_eligible_competing_structure_count": selection.eligible_structure_count,
        f"{prefix}_strong_competing_structure": str(selection.strong_competing_structure).lower(),
        f"{prefix}_strong_endpoint_multimodality": str(selection.strong_multimodality).lower(),
        f"{prefix}_local_support_fraction": f"{selection.local_support_fraction:.6f}",
        f"{prefix}_adaptive_window_bp": selection.adaptive_window_bp,
        f"{prefix}_exonic_limit_bp": limits.value(side, "exonic") or ".",
        f"{prefix}_genomic_limit_bp": limits.value(side, "genomic") or ".",
        f"{prefix}_selected_exonic_extension_bp": candidate.extension_exonic_bp,
        f"{prefix}_selected_genomic_extension_bp": candidate.extension_genomic_bp,
        f"{prefix}_anchor_block_reason": anchor_block_reason or ".",
        f"{prefix}_quality_status": selection.quality_status,
    }


def estimate_catalogue_bytes(evidence: Sequence[UTREvidence]) -> int:
    """Return a conservative shallow-plus-payload catalogue estimate."""
    total = sys.getsizeof(evidence)
    for item in evidence:
        total += sys.getsizeof(item) + sys.getsizeof(item.evidence_id)
        total += sys.getsizeof(item.exon_parts)
        total += sum(sys.getsizeof(interval) for interval in item.exon_parts)
        total += sys.getsizeof(item.source_member_counts)
        total += sum(
            sys.getsizeof(source) + sys.getsizeof(count)
            for source, count in item.source_member_counts
        )
        total += sys.getsizeof(item.representative_ids)
    return total


@dataclass(frozen=True, slots=True)
class _SelectedFeatureOwner:
    gene_id: str
    transcript_id: str
    contributing_sides: tuple[Side, ...]


def _selected_feature_entries(
    gene_id: str,
    transcript: FinalTranscript,
    selections: tuple[UTRSideSelection | None, UTRSideSelection | None],
    *,
    only_keys: set[StructuralFeature] | None = None,
) -> tuple[set[StructuralFeature], list[tuple[StructuralFeature, _SelectedFeatureOwner]]]:
    model = transcript.source_model
    base = _structural_features(model.chrom, model.strand, model.cds_parts)
    side_exons = {
        side: selection.candidate.added_exons
        for side, selection in zip(SIDES, selections)
        if selection is not None
    }
    final = _structural_features(
        model.chrom, model.strand,
        _materialized_exons(model, *side_exons.values()),
    )
    side_new = {
        side: _structural_features(
            model.chrom, model.strand, _materialized_exons(model, added)
        ) - base
        for side, added in side_exons.items()
    }
    new_keys = final - base
    entries: list[tuple[StructuralFeature, _SelectedFeatureOwner]] = []
    for key in final:
        if only_keys is not None and key not in only_keys:
            continue
        contributing = tuple(
            side for side in SIDES if key in side_new.get(side, set())
        )
        if key in new_keys and not contributing:
            # An exon spanning both sides of a single-exon CDS can exist only
            # after both independently selected extensions are materialized.
            contributing = tuple(side for side in SIDES if side in side_exons)
        entries.append((key, _SelectedFeatureOwner(
            gene_id, transcript.transcript_id, contributing
        )))
    return new_keys, entries


def _selected_global_conflicts(
    transcript_lookup: Mapping[str, tuple[str, FinalTranscript]],
    selected: Mapping[str, tuple[UTRSideSelection | None, UTRSideSelection | None]],
    frozen: FrozenFinalGeneStructureIndex,
) -> tuple[dict[tuple[str, Side], set[str]], set[tuple[str, Side]]]:
    """Find all new cross-gene exact links simultaneously, without ownership order."""
    new_keys: set[StructuralFeature] = set()
    for transcript_id, (gene_id, transcript) in transcript_lookup.items():
        keys, _ = _selected_feature_entries(
            gene_id, transcript, selected[transcript_id], only_keys=set()
        )
        new_keys.update(keys)
    if not new_keys:
        return {}, set()

    by_feature: dict[StructuralFeature, list[_SelectedFeatureOwner]] = defaultdict(list)
    for transcript_id, (gene_id, transcript) in transcript_lookup.items():
        _keys, entries = _selected_feature_entries(
            gene_id, transcript, selected[transcript_id], only_keys=new_keys
        )
        for key, owner in entries:
            by_feature[key].append(owner)

    conflicts: dict[tuple[str, Side], set[str]] = defaultdict(set)
    utr_against_utr: set[tuple[str, Side]] = set()
    for key, owners in by_feature.items():
        if len(owners) < 2 or len({owner.gene_id for owner in owners}) < 2:
            continue
        label = "splice" if key[2] == "splice" else key[2]
        for index, left in enumerate(owners):
            for right in owners[index + 1:]:
                if left.gene_id == right.gene_id:
                    continue
                if (
                    key in frozen.baseline_by_transcript[left.transcript_id]
                    and key in frozen.baseline_by_transcript[right.transcript_id]
                ):
                    continue  # This relationship already existed before UTRs.
                if not left.contributing_sides and not right.contributing_sides:
                    continue
                prefix = (
                    "cross_final_gene_selected_utr"
                    if left.contributing_sides and right.contributing_sides
                    else "cross_final_gene_authoritative"
                )
                reason = f"{prefix}_{label}_conflict"
                for owner in (left, right):
                    for side in owner.contributing_sides:
                        proposal = (owner.transcript_id, side)
                        conflicts[proposal].add(reason)
                        if prefix == "cross_final_gene_selected_utr":
                            utr_against_utr.add(proposal)
    return dict(conflicts), utr_against_utr


class _SelectedPopulationStructureIndex:
    """Track exact selected structural links, refreshing only changed features.

    The full-population scan above remains an independent reference for tests.
    Updates are applied simultaneously so neither gene can win a conflict by
    being processed first.
    """

    def __init__(
        self,
        transcript_lookup: Mapping[str, tuple[str, FinalTranscript]],
        selected: Mapping[str, tuple[UTRSideSelection | None, UTRSideSelection | None]],
        frozen: FrozenFinalGeneStructureIndex,
    ) -> None:
        self.transcript_lookup = transcript_lookup
        self.frozen = frozen
        self.by_feature: dict[StructuralFeature, dict[str, _SelectedFeatureOwner]] = {}
        self.selected_by_transcript = dict(selected)
        self.conflicts_by_feature: dict[
            StructuralFeature, dict[tuple[str, Side], set[str]]
        ] = {}
        self.utr_conflicts_by_feature: dict[
            StructuralFeature, set[tuple[str, Side]]
        ] = {}
        for transcript_id, (gene_id, transcript) in transcript_lookup.items():
            _new_keys, entries = _selected_feature_entries(
                gene_id, transcript, selected[transcript_id]
            )
            for key, owner in entries:
                self.by_feature.setdefault(key, {})[transcript_id] = owner
        for key, owners in self.by_feature.items():
            if len(owners) > 1:
                self._refresh_conflict(key)

    def _refresh_conflict(self, key: StructuralFeature) -> None:
        owners = list(self.by_feature.get(key, {}).values())
        if len(owners) < 2 or len({owner.gene_id for owner in owners}) < 2:
            self.conflicts_by_feature.pop(key, None)
            self.utr_conflicts_by_feature.pop(key, None)
            return
        label = "splice" if key[2] == "splice" else key[2]
        conflicts: dict[tuple[str, Side], set[str]] = defaultdict(set)
        utr_against_utr: set[tuple[str, Side]] = set()
        for index, left in enumerate(owners):
            for right in owners[index + 1:]:
                if left.gene_id == right.gene_id:
                    continue
                if (
                    key in self.frozen.baseline_by_transcript[left.transcript_id]
                    and key in self.frozen.baseline_by_transcript[right.transcript_id]
                ):
                    continue
                if not left.contributing_sides and not right.contributing_sides:
                    continue
                prefix = (
                    "cross_final_gene_selected_utr"
                    if left.contributing_sides and right.contributing_sides
                    else "cross_final_gene_authoritative"
                )
                reason = f"{prefix}_{label}_conflict"
                for owner in (left, right):
                    for side in owner.contributing_sides:
                        proposal = (owner.transcript_id, side)
                        conflicts[proposal].add(reason)
                        if prefix == "cross_final_gene_selected_utr":
                            utr_against_utr.add(proposal)
        if conflicts:
            self.conflicts_by_feature[key] = dict(conflicts)
            self.utr_conflicts_by_feature[key] = utr_against_utr
        else:
            self.conflicts_by_feature.pop(key, None)
            self.utr_conflicts_by_feature.pop(key, None)

    def current_conflicts(
        self,
    ) -> tuple[dict[tuple[str, Side], set[str]], set[tuple[str, Side]]]:
        conflicts: dict[tuple[str, Side], set[str]] = defaultdict(set)
        for feature_conflicts in self.conflicts_by_feature.values():
            for proposal, reasons in feature_conflicts.items():
                conflicts[proposal].update(reasons)
        utr_against_utr: set[tuple[str, Side]] = set()
        for proposals in self.utr_conflicts_by_feature.values():
            utr_against_utr.update(proposals)
        return dict(conflicts), utr_against_utr

    def replace(
        self,
        selected: Mapping[str, tuple[UTRSideSelection | None, UTRSideSelection | None]],
    ) -> None:
        """Commit one resolution round and recheck only affected exact keys."""
        touched: set[StructuralFeature] = set()
        for transcript_id in sorted(selected):
            gene_id, transcript = self.transcript_lookup[transcript_id]
            _old_new_keys, old_entries = _selected_feature_entries(
                gene_id, transcript, self.selected_by_transcript[transcript_id]
            )
            _new_keys, entries = _selected_feature_entries(
                gene_id, transcript, selected[transcript_id]
            )
            current = dict(entries)
            previous = dict(old_entries)
            for key in previous.keys() | current.keys():
                if previous.get(key) == current.get(key):
                    continue
                owners = self.by_feature.get(key)
                if owners is not None:
                    owners.pop(transcript_id, None)
                    if not owners:
                        del self.by_feature[key]
                if key in current:
                    self.by_feature.setdefault(key, {})[transcript_id] = current[key]
                touched.add(key)
            self.selected_by_transcript[transcript_id] = selected[transcript_id]
        for key in touched:
            self._refresh_conflict(key)


def augment_final_annotation(
    annotation: FinalAnnotation,
    transcript_evidence: Sequence[UTREvidence],
    prediction_evidence: Sequence[UTREvidence],
    settings: UTRSettings,
    *,
    learned_limits: UTRLearnedLimits | None = None,
    excluded_training_sides: set[tuple[str, Side]] | None = None,
    fasta: IndexedFasta | None = None,
    validation_recorder: Callable[[UTRSideValidation], None] | None = None,
) -> UTRAugmentationResult:
    started = time.monotonic()
    all_evidence = list(transcript_evidence) + list(prediction_evidence)
    index = UTREvidenceIndex(all_evidence, annotation)
    frozen = FrozenFinalGeneStructureIndex(annotation)
    bridged_evidence_ids = _multi_gene_bridge_ids(annotation, index, settings, fasta)
    limits = learned_limits or _streaming_extension_limits(
        annotation,
        index,
        settings,
        bridged_evidence_ids,
        frozen,
        excluded_training_sides,
        fasta,
    )
    rejection_counts: Counter[str] = Counter()

    transcript_lookup: dict[str, tuple[str, FinalTranscript]] = {}
    candidate_counts: dict[str, Counter[Side]] = {}
    hard_rejections: Counter[tuple[str, Side]] = Counter()
    hard_reasons: dict[tuple[str, Side], set[str]] = defaultdict(set)
    selected: dict[str, tuple[UTRSideSelection | None, UTRSideSelection | None]] = {}
    initial_selected: dict[str, tuple[UTRSideSelection | None, UTRSideSelection | None]] = {}

    # First select against the frozen gene population.  Candidate tables remain
    # transcript-local; only the compact selected pair is retained per public
    # transcript for the later global validation pass.
    for gene in annotation.genes:
        for transcript in gene.transcripts:
            transcript_id = transcript.transcript_id
            transcript_lookup[transcript_id] = (gene.gene_id, transcript)
            candidates = candidate_sides_for_transcript(
                gene.gene_id,
                transcript,
                index,
                settings,
                fasta=fasta,
                validation_recorder=validation_recorder,
            )
            for side in SIDES:
                anchor_block_reason = _incomplete_cds_utr_block_reason(
                    transcript.source_model, side
                )
                if anchor_block_reason:
                    rejection_counts[
                        f"incomplete_cds_anchor_{anchor_block_reason}"
                    ] += 1
            candidate_counts[transcript_id] = Counter(
                candidate.side for candidate in candidates
            )
            for candidate in candidates:
                candidate.multi_gene_bridge = (
                    candidate.evidence.evidence_id in bridged_evidence_ids
                )
            frozen_conflicts = {
                id(candidate): frozen.candidate_reasons(candidate)
                for candidate in candidates
            }
            had_hard_conflict = any(frozen_conflicts.values())
            pre_hard = (
                _select_sides_for_transcript(candidates, settings, limits)
                if had_hard_conflict else None
            )
            for candidate in candidates:
                conflicts = frozen_conflicts[id(candidate)]
                if conflicts:
                    if candidate.eligible_nonextension:
                        key = (transcript_id, candidate.side)
                        hard_rejections[key] += 1
                        hard_reasons[key].update(conflicts)
                    candidate.rejection_reasons = tuple(sorted(
                        set(candidate.rejection_reasons) | set(conflicts)
                    ))
                rejection_counts.update(candidate.rejection_reasons)
                if candidate.multi_gene_bridge:
                    rejection_counts["multi_final_gene_bridge"] += 1
                elif not candidate.rejection_reasons:
                    rejection_counts.update(_extension_reasons(candidate, limits))
            chosen = _select_sides_for_transcript(candidates, settings, limits)
            selected[transcript_id] = chosen
            initial_selected[transcript_id] = pre_hard or chosen

    # Simultaneously invalidate both sides of every selected UTR-vs-UTR link.
    # All affected transcripts reselect from their original local candidates
    # after the forbidden exact side geometry is removed.  The procedure has no
    # first-gene-wins behavior and terminates because every round bans at least
    # one previously selected exact side proposal.
    blocked: dict[str, dict[Side, set[tuple[Interval, ...]]]] = defaultdict(
        lambda: {side: set() for side in SIDES}
    )
    selected_conflict_reasons: dict[tuple[str, Side], set[str]] = defaultdict(set)
    first_round_utr_conflict_sides: set[tuple[str, Side]] = set()
    resolution_rounds = 0
    selected_index = _SelectedPopulationStructureIndex(transcript_lookup, selected, frozen)
    fallback_candidates: dict[str, list[UTRCandidateSide]] = {}
    while True:
        conflicts, utr_conflict_sides = selected_index.current_conflicts()
        if not conflicts:
            break
        if resolution_rounds == 0:
            first_round_utr_conflict_sides = utr_conflict_sides
        resolution_rounds += 1
        affected_transcripts: set[str] = set()
        for (transcript_id, side), reasons in sorted(conflicts.items()):
            choice = selected[transcript_id][SIDES.index(side)]
            if choice is None:
                raise RuntimeError("Cross-gene conflict had no selected UTR side")
            signature = choice.candidate.added_exons
            if signature in blocked[transcript_id][side]:
                raise RuntimeError("Cross-gene UTR resolution did not progress")
            blocked[transcript_id][side].add(signature)
            selected_conflict_reasons[(transcript_id, side)].update(reasons)
            rejection_counts.update(reasons)
            affected_transcripts.add(transcript_id)
        for transcript_id in sorted(affected_transcripts):
            candidates = fallback_candidates.get(transcript_id)
            if candidates is None:
                gene_id, transcript = transcript_lookup[transcript_id]
                candidates = candidate_sides_for_transcript(
                    gene_id, transcript, index, settings, fasta=fasta
                )
                for candidate in candidates:
                    candidate.multi_gene_bridge = (
                        candidate.evidence.evidence_id in bridged_evidence_ids
                    )
                _apply_frozen_gene_filter(candidates, frozen)
                fallback_candidates[transcript_id] = candidates
            selected[transcript_id] = _select_sides_for_transcript(
                candidates, settings, limits,
                blocked_signatures=blocked[transcript_id],
            )
        selected_index.replace({
            transcript_id: selected[transcript_id]
            for transcript_id in affected_transcripts
        })

    augmented_genes: list[AugmentedGene] = []
    audit_rows: list[dict[str, object]] = []
    coverage: Counter[str] = Counter()
    evidence_coverage: Counter[str] = Counter()
    endpoint_methods: Counter[str] = Counter()
    strong_competing_sides = 0
    multimodal_sides = 0
    for gene in annotation.genes:
        augmented_transcripts: list[AugmentedTranscript] = []
        for transcript in gene.transcripts:
            transcript_id = transcript.transcript_id
            five, three = selected[transcript_id]
            selections = [selection for selection in (five, three) if selection is not None]
            if any(
                selection.candidate.validation_status == "rejected"
                for selection in selections
            ):
                raise RuntimeError("Rejected UTR candidate reached final reconstruction")
            added = [
                interval
                for selection in selections
                for interval in selection.candidate.added_exons
            ]
            exon_parts = _materialized_exons(transcript.source_model, added)
            five_parts = five.candidate.added_exons if five is not None else ()
            three_parts = three.candidate.added_exons if three is not None else ()
            shared = bool(
                five is not None and three is not None
                and set(five.selected_evidence_ids) & set(three.selected_evidence_ids)
            )
            augmented = AugmentedTranscript(
                final_transcript=transcript,
                gene_id=gene.gene_id,
                exon_parts=exon_parts,
                five_prime_utr_parts=five_parts,
                three_prime_utr_parts=three_parts,
                five_prime_selection=five,
                three_prime_selection=three,
                choices_share_evidence=shared,
            )
            augmented_transcripts.append(augmented)
            state = (
                "both" if five and three else
                "five_prime_only" if five else
                "three_prime_only" if three else
                "neither"
            )
            coverage[state] += 1
            kinds = {selection.evidence_class for selection in selections}
            flattened: set[str] = set()
            for value in kinds:
                if value == "both":
                    flattened.update(("transcript", "prediction"))
                else:
                    flattened.add(value)
            evidence_state = (
                "both" if {"transcript", "prediction"} <= flattened else
                "transcript" if "transcript" in flattened else
                "prediction" if "prediction" in flattened else
                "none"
            )
            evidence_coverage[evidence_state] += 1
            for selection in selections:
                endpoint_methods[selection.endpoint_method] += 1
                strong_competing_sides += selection.strong_competing_structure
                multimodal_sides += selection.strong_multimodality
            side_counts = candidate_counts[transcript_id]
            audit_row: dict[str, object] = {
                "public_gene_id": gene.gene_id,
                "public_transcript_id": transcript.transcript_id,
                "transcript_role": transcript.transcript_role,
                "chrom": transcript.source_model.chrom,
                "strand": transcript.source_model.strand,
                **_audit_side_fields(
                    "five_prime", five, side_counts["five_prime"], limits,
                    _incomplete_cds_utr_block_reason(
                        transcript.source_model, "five_prime"
                    ),
                ),
                **_audit_side_fields(
                    "three_prime", three, side_counts["three_prime"], limits,
                    _incomplete_cds_utr_block_reason(
                        transcript.source_model, "three_prime"
                    ),
                ),
                "five_three_choices_share_evidence": str(shared).lower(),
            }
            for side in SIDES:
                key = (transcript_id, side)
                prefix = side
                current = selected[transcript_id][SIDES.index(side)]
                original = initial_selected[transcript_id][SIDES.index(side)]
                reasons = hard_reasons.get(key, set()) | selected_conflict_reasons.get(key, set())
                audit_row[f"{prefix}_cross_final_gene_conflict_reasons"] = (
                    ",".join(sorted(reasons)) or "."
                )
                audit_row[f"{prefix}_authoritative_candidate_rejections"] = hard_rejections[key]
                audit_row[f"{prefix}_selected_conflict_rejections"] = len(blocked[transcript_id][side])
                audit_row[f"{prefix}_initial_structure_key"] = (
                    original.candidate.structure_key if original else "."
                )
                audit_row[f"{prefix}_initial_endpoint0"] = (
                    original.candidate.endpoint0 if original else "."
                )
                if key in selected_conflict_reasons:
                    audit_row[f"{prefix}_quality_status"] = (
                        "cross_final_gene_conflict_fallback" if current
                        else "cross_final_gene_conflict_no_alternative"
                    )
                    rejection_counts[audit_row[f"{prefix}_quality_status"]] += 1
                elif original is not None and (
                    current is None
                    or current.candidate.added_exons != original.candidate.added_exons
                ):
                    audit_row[f"{prefix}_quality_status"] = (
                        "cross_final_gene_authoritative_filter_fallback" if current
                        else "cross_final_gene_authoritative_filter_no_alternative"
                    )
                    rejection_counts[audit_row[f"{prefix}_quality_status"]] += 1
            audit_rows.append(audit_row)
        augmented_genes.append(AugmentedGene(gene, tuple(augmented_transcripts)))

    affected_transcript_ids = {
        transcript_id for transcript_id, _side in
        set(hard_rejections) | set(selected_conflict_reasons)
    }
    conflict_fallback = sum(
        selected[transcript_id][SIDES.index(side)] is not None
        for transcript_id, side in selected_conflict_reasons
    )

    summary = UTRAuditSummary(
        transcript_count=annotation.transcript_count,
        both_assigned=coverage["both"],
        five_prime_only=coverage["five_prime_only"],
        three_prime_only=coverage["three_prime_only"],
        neither_assigned=coverage["neither"],
        transcript_only=evidence_coverage["transcript"],
        prediction_only=evidence_coverage["prediction"],
        both_evidence_kinds=evidence_coverage["both"],
        multi_gene_evidence_excluded=len(bridged_evidence_ids),
        authoritative_candidate_sides_rejected=sum(hard_rejections.values()),
        provisional_selected_sides_in_utr_conflicts=len(first_round_utr_conflict_sides),
        conflict_fallback_sides=conflict_fallback,
        conflict_no_alternative_sides=len(selected_conflict_reasons) - conflict_fallback,
        conflict_affected_genes=len({
            transcript_lookup[transcript_id][0]
            for transcript_id in affected_transcript_ids
        }),
        conflict_affected_transcripts=len(affected_transcript_ids),
        conflict_resolution_rounds=resolution_rounds,
        local_consensus_sides=endpoint_methods["local_consensus_medoid"],
        diffuse_sides=endpoint_methods["diffuse_global_medoid"],
        strong_competing_structure_sides=strong_competing_sides,
        strongly_multimodal_sides=multimodal_sides,
        transcript_evidence_records=len(transcript_evidence),
        prediction_evidence_records=len(prediction_evidence),
        catalogue_estimated_bytes=estimate_catalogue_bytes(all_evidence),
        runtime_seconds=time.monotonic() - started,
        rejection_counts=tuple(sorted(rejection_counts.items())),
    )
    return UTRAugmentationResult(
        annotation=AugmentedAnnotation(tuple(augmented_genes)),
        limits=limits,
        audit_rows=tuple(audit_rows),
        summary=summary,
    )


def write_utr_augmented_gff3(
    result: UTRAugmentationResult,
    path: Path,
) -> None:
    with path.open("w") as out:
        out.write("##gff-version 3\n")
        out.write("# MoSelecta post-selection UTR-augmented annotation\n")
        out.write("# moselecta_utr_selected_cds_authoritative=true\n")
        for side in SIDES:
            for geometry in ("exonic", "genomic"):
                value = result.limits.value(side, geometry)
                out.write(
                    f"# moselecta_utr_{side}_{geometry}_limit_bp="
                    f"{value if value is not None else '.'}\n"
                )
        for gene in result.annotation.genes:
            final_gene = gene.final_gene
            out.write(feature_line(
                seqid=final_gene.chrom,
                source=GFF3_SOURCE,
                feature_type="gene",
                start0=gene.start0,
                end0=gene.end0,
                score=".",
                strand=final_gene.strand,
                phase=".",
                attrs={
                    "ID": final_gene.gene_id,
                    "moselecta_primary_transcript_id": final_gene.primary_transcript.transcript_id,
                    "moselecta_transcript_count": len(gene.transcripts),
                    "moselecta_utr_augmented": "true",
                },
            ) + "\n")
            for transcript in gene.transcripts:
                model = transcript.model
                five = transcript.five_prime_selection
                three = transcript.three_prime_selection
                attrs = public_transcript_attributes(
                    final_gene,
                    transcript.final_transcript,
                    include_alternative_metadata=True,
                )
                attrs.update({
                    "moselecta_utr_augmented": "true",
                    "moselecta_five_prime_utr": "assigned" if five else "none",
                    "moselecta_three_prime_utr": "assigned" if three else "none",
                })
                if five:
                    attrs["moselecta_five_prime_utr_evidence"] = five.evidence_class
                if three:
                    attrs["moselecta_three_prime_utr_evidence"] = three.evidence_class
                out.write(feature_line(
                    seqid=model.chrom,
                    source=f"{GFF3_SOURCE}:{model.source_name}",
                    feature_type="mRNA",
                    start0=transcript.start0,
                    end0=transcript.end0,
                    score=model.transcript_record.score if model.transcript_record else ".",
                    strand=model.strand,
                    phase=".",
                    attrs=attrs,
                ) + "\n")
                features: list[tuple[int, int, str, str, dict[str, object]]] = []
                for index, (start0, end0) in enumerate(transcript.exon_parts, start=1):
                    features.append((start0, end0, "exon", ".", {
                        "ID": f"{transcript.transcript_id}{PUBLIC_CHILD_ID_SEPARATOR}exon{index}",
                        "Parent": transcript.transcript_id,
                        "generated_from": "authoritative_CDS_and_selected_UTR",
                    }))
                for feature_type, parts in (
                    ("five_prime_UTR", transcript.five_prime_utr_parts),
                    ("three_prime_UTR", transcript.three_prime_utr_parts),
                ):
                    for index, (start0, end0) in enumerate(parts, start=1):
                        features.append((start0, end0, feature_type, ".", {
                            "ID": (
                                f"{transcript.transcript_id}{PUBLIC_CHILD_ID_SEPARATOR}"
                                f"{feature_type}{index}"
                            ),
                            "Parent": transcript.transcript_id,
                            "moselecta_utr_selected": "true",
                        }))
                phases = emitted_cds_phase_by_part(model)
                for index, (start0, end0) in enumerate(sorted(model.cds_parts), start=1):
                    features.append((start0, end0, "CDS", phases[(start0, end0)], {
                        "ID": f"{transcript.transcript_id}{PUBLIC_CHILD_ID_SEPARATOR}CDS{index}",
                        "Parent": transcript.transcript_id,
                        "generated_from": "authoritative_CDS",
                    }))
                for start0, end0, feature_type, phase, feature_attrs in sorted(
                    features,
                    key=lambda item: (
                        item[0], CHILD_FEATURE_ORDER.get(item[2], 99), item[1], item[2]
                    ),
                ):
                    out.write(feature_line(
                        seqid=model.chrom,
                        source=f"{GFF3_SOURCE}:{model.source_name}",
                        feature_type=feature_type,
                        start0=start0,
                        end0=end0,
                        score=".",
                        strand=model.strand,
                        phase=phase,
                        attrs=feature_attrs,
                    ) + "\n")
            out.write("###\n")


def write_utr_audit_tsv(result: UTRAugmentationResult, path: Path) -> None:
    rows = list(result.audit_rows)
    if not rows:
        path.write_text("")
        return
    write_dict_tsv(rows, path, list(rows[0]))


def write_utr_summary_tsv(result: UTRAugmentationResult, path: Path) -> None:
    summary = result.summary
    rejection_text = ",".join(
        f"{reason}:{count}" for reason, count in summary.rejection_counts
    ) or "."
    row = {
        "final_public_transcript_count": summary.transcript_count,
        "both_sides_assigned": summary.both_assigned,
        "five_prime_only_assigned": summary.five_prime_only,
        "three_prime_only_assigned": summary.three_prime_only,
        "neither_side_assigned": summary.neither_assigned,
        "transcript_only_evidence": summary.transcript_only,
        "prediction_only_evidence": summary.prediction_only,
        "both_evidence_kinds": summary.both_evidence_kinds,
        "multi_gene_evidence_excluded": summary.multi_gene_evidence_excluded,
        "authoritative_candidate_sides_rejected": summary.authoritative_candidate_sides_rejected,
        "provisional_selected_sides_in_utr_conflicts": summary.provisional_selected_sides_in_utr_conflicts,
        "cross_final_gene_conflict_fallback_sides": summary.conflict_fallback_sides,
        "cross_final_gene_conflict_no_alternative_sides": summary.conflict_no_alternative_sides,
        "cross_final_gene_conflict_affected_genes": summary.conflict_affected_genes,
        "cross_final_gene_conflict_affected_transcripts": summary.conflict_affected_transcripts,
        "cross_final_gene_conflict_resolution_rounds": summary.conflict_resolution_rounds,
        "quality_rejection_counts": rejection_text,
        "five_prime_exonic_limit_bp": result.limits.five_prime_exonic_bp or ".",
        "three_prime_exonic_limit_bp": result.limits.three_prime_exonic_bp or ".",
        "five_prime_genomic_limit_bp": result.limits.five_prime_genomic_bp or ".",
        "three_prime_genomic_limit_bp": result.limits.three_prime_genomic_bp or ".",
        "limit_warnings": ";".join(result.limits.warnings) or ".",
        "local_consensus_sides": summary.local_consensus_sides,
        "diffuse_fallback_sides": summary.diffuse_sides,
        "strong_competing_structure_sides": summary.strong_competing_structure_sides,
        "strongly_multimodal_sides": summary.strongly_multimodal_sides,
        "transcript_evidence_records": summary.transcript_evidence_records,
        "prediction_evidence_records": summary.prediction_evidence_records,
        "catalogue_estimated_bytes": summary.catalogue_estimated_bytes,
        "utr_stage_runtime_seconds": f"{summary.runtime_seconds:.6f}",
    }
    write_dict_tsv([row], path, list(row))
