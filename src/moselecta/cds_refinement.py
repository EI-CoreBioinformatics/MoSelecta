"""Post-selection transcript-guided CDS terminus refinement.

The selected source model remains immutable.  Successful proposals are emitted
as replacement ``CodingModel`` values held by a new ``FinalAnnotation``; all
selection scores and public transcript metadata remain original provenance.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

from .cds_phase import calculate_cds_phase_chain
from .final_annotation import FinalAnnotation, FinalGene, FinalTranscript
from .genome_sequence import (
    CANONICAL_SPLICE_MOTIFS,
    STANDARD_STOP_CODONS,
    IndexedFasta,
    calculate_model_sequence_integrity,
    reverse_complement,
    sequence_integrity_metrics_from_row,
)
from .intervals import IntervalOverlapIndex, merge_intervals
from .loci import models_share_sublocus_relationship
from .models import CodingModel, FeatureRecord
from .tsvio import write_dict_tsv
from .utr_augmentation import (
    FrozenFinalGeneStructureIndex,
    Side,
    UTRCandidateSide,
    UTREvidence,
    UTREvidenceIndex,
    UTRLearnedLimits,
    UTRSettings,
    candidate_sides_for_transcript,
    mark_multi_gene_bridges,
    select_completion_candidate,
)


Interval = tuple[int, int]


@dataclass(frozen=True, slots=True)
class CDSTerminusRefinementSettings:
    downstream_start_search_bp: int = 30
    downstream_stop_search_bp: int = 100
    max_relative_extension: float = 1.0
    max_relative_truncation: float = 0.33
    min_cds_length_bp: int = 150
    single_exon_min_reciprocal_overlap: float = 0.2
    genetic_code: int = 1
    internal_stop_penalty: float = 1.0
    internal_stop_count_cap: int = 2
    phase_inconsistent_penalty: float = 0.75
    missing_stop_penalty: float = 0.50
    missing_start_penalty: float = 0.30
    noncanonical_splice_penalty: float = 1.0

    def __post_init__(self) -> None:
        for name in (
            "downstream_start_search_bp",
            "downstream_stop_search_bp",
            "min_cds_length_bp",
        ):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be >= 0")
        if (
            not math.isfinite(self.max_relative_extension)
            or self.max_relative_extension < 0.0
        ):
            raise ValueError("max_relative_extension must be >= 0")
        if (
            not math.isfinite(self.max_relative_truncation)
            or not 0.0 <= self.max_relative_truncation <= 1.0
        ):
            raise ValueError("max_relative_truncation must be between 0 and 1")


@dataclass(frozen=True, slots=True)
class _Proposal:
    side: Side
    direction: str
    cds_parts: tuple[Interval, ...]
    residual_utr_parts: tuple[Interval, ...]
    codon: str
    codon_positions: tuple[int, ...]
    bp_delta: int
    bases_searched: int
    splice_junctions_crossed: int
    orf_open_at_search_boundary: bool | None
    selection: object | None


@dataclass(frozen=True, slots=True)
class CDSTerminusRefinementSummary:
    execution_status: str
    audit_sides: int
    incomplete_five_prime_sides: int
    incomplete_three_prime_sides: int
    transcripts_with_two_incomplete_sides: int
    final_refined_sides: int
    final_failed_sides: int
    refined_sides_upstream_start: int
    refined_sides_downstream_start_fallback: int
    refined_sides_downstream_stop: int
    refined_transcripts: int
    unchanged_transcripts: int
    two_incomplete_sides_both_refined_transcripts: int
    two_incomplete_sides_start_only_refined_transcripts: int
    two_incomplete_sides_stop_only_refined_transcripts: int
    two_incomplete_sides_neither_refined_transcripts: int
    refined_sides_with_encountered_rejection_events: int
    final_outcome_side_counts: tuple[tuple[str, int], ...]
    failure_counts: tuple[tuple[str, int], ...]
    failure_transcript_counts: tuple[tuple[str, int], ...]
    encountered_rejection_event_counts: tuple[tuple[str, int], ...]
    encountered_rejection_side_counts: tuple[tuple[str, int], ...]
    encountered_rejection_transcript_counts: tuple[tuple[str, int], ...]


@dataclass(frozen=True, slots=True)
class CDSTerminusRefinementResult:
    original_annotation: FinalAnnotation
    annotation: FinalAnnotation
    audit_rows: tuple[dict[str, object], ...]
    summary: CDSTerminusRefinementSummary
    originally_incomplete_sides: frozenset[tuple[str, Side]]


def empty_cds_refinement_result(
    annotation: FinalAnnotation,
    *,
    execution_status: str,
) -> CDSTerminusRefinementResult:
    """Return an explicit zero-count result for a refinement stage not run."""
    return CDSTerminusRefinementResult(
        annotation,
        annotation,
        (),
        CDSTerminusRefinementSummary(
            execution_status=execution_status,
            audit_sides=0,
            incomplete_five_prime_sides=0,
            incomplete_three_prime_sides=0,
            transcripts_with_two_incomplete_sides=0,
            final_refined_sides=0,
            final_failed_sides=0,
            refined_sides_upstream_start=0,
            refined_sides_downstream_start_fallback=0,
            refined_sides_downstream_stop=0,
            refined_transcripts=0,
            unchanged_transcripts=0,
            two_incomplete_sides_both_refined_transcripts=0,
            two_incomplete_sides_start_only_refined_transcripts=0,
            two_incomplete_sides_stop_only_refined_transcripts=0,
            two_incomplete_sides_neither_refined_transcripts=0,
            refined_sides_with_encountered_rejection_events=0,
            final_outcome_side_counts=(),
            failure_counts=(),
            failure_transcript_counts=(),
            encountered_rejection_event_counts=(),
            encountered_rejection_side_counts=(),
            encountered_rejection_transcript_counts=(),
        ),
        frozenset(),
    )


@dataclass(frozen=True, slots=True)
class _SplicedPath:
    """Transcript-oriented intervals with O(part-count) coordinate state."""

    parts: tuple[Interval, ...]
    strand: str
    offsets: tuple[int, ...]

    @classmethod
    def from_parts(cls, parts: Sequence[Interval], strand: str) -> _SplicedPath:
        if strand not in {"+", "-"}:
            raise ValueError(f"Unsupported strand {strand!r}")
        merged = tuple(merge_intervals(parts))
        ordered = tuple(sorted(merged, reverse=strand == "-"))
        offsets = [0]
        for start0, end0 in ordered:
            if end0 <= start0:
                raise ValueError("Spliced path parts must have positive length")
            offsets.append(offsets[-1] + end0 - start0)
        return cls(ordered, strand, tuple(offsets))

    @property
    def length(self) -> int:
        return self.offsets[-1] if self.offsets else 0

    def offset_of_position(self, position: int) -> int | None:
        for index, (start0, end0) in enumerate(self.parts):
            if start0 <= position < end0:
                within = (
                    position - start0
                    if self.strand == "+"
                    else end0 - 1 - position
                )
                return self.offsets[index] + within
        return None

    def genomic_position(self, offset: int) -> int:
        if offset < 0 or offset >= self.length:
            raise IndexError(f"Spliced offset {offset} outside path length {self.length}")
        index = bisect_right(self.offsets, offset) - 1
        start0, end0 = self.parts[index]
        within = offset - self.offsets[index]
        return start0 + within if self.strand == "+" else end0 - 1 - within

    def slice_parts(self, start: int, end: int) -> tuple[Interval, ...]:
        if start < 0 or end < start or end > self.length:
            raise IndexError(f"Invalid spliced slice {start}:{end} for {self.length}")
        result: list[Interval] = []
        for index, (part_start, part_end) in enumerate(self.parts):
            overlap_start = max(start, self.offsets[index])
            overlap_end = min(end, self.offsets[index + 1])
            if overlap_start >= overlap_end:
                continue
            local_start = overlap_start - self.offsets[index]
            local_end = overlap_end - self.offsets[index]
            if self.strand == "+":
                result.append((part_start + local_start, part_start + local_end))
            else:
                result.append((part_end - local_end, part_end - local_start))
        return tuple(merge_intervals(result))

    def sequence_chunks(
        self,
        fasta: IndexedFasta,
        chrom: str,
        start: int,
        end: int,
    ) -> Iterable[tuple[int, str]]:
        """Yield oriented sequence chunks without materializing path coordinates."""
        if start < 0 or end < start or end > self.length:
            raise IndexError(f"Invalid spliced slice {start}:{end} for {self.length}")
        for index, (part_start, part_end) in enumerate(self.parts):
            overlap_start = max(start, self.offsets[index])
            overlap_end = min(end, self.offsets[index + 1])
            if overlap_start >= overlap_end:
                continue
            local_start = overlap_start - self.offsets[index]
            local_end = overlap_end - self.offsets[index]
            if self.strand == "+":
                genomic_start = part_start + local_start
                genomic_end = part_start + local_end
                sequence = fasta.fetch(chrom, genomic_start, genomic_end)
            else:
                genomic_start = part_end - local_end
                genomic_end = part_end - local_start
                sequence = reverse_complement(
                    fasta.fetch(chrom, genomic_start, genomic_end)
                )
            yield overlap_start, sequence

    def sequence(
        self, fasta: IndexedFasta, chrom: str, start: int, end: int
    ) -> str:
        return "".join(
            sequence
            for _offset, sequence in self.sequence_chunks(fasta, chrom, start, end)
        )

    def junction_count(self, start: int, end_inclusive: int) -> int:
        """Count part boundaries between two included spliced bases."""
        if end_inclusive <= start:
            return 0
        return sum(start < boundary <= end_inclusive for boundary in self.offsets[1:-1])


def _parts_from_positions(positions: Iterable[int]) -> tuple[Interval, ...]:
    ordered = sorted(set(positions))
    if not ordered:
        return ()
    parts: list[Interval] = []
    start = previous = ordered[0]
    for position in ordered[1:]:
        if position != previous + 1:
            parts.append((start, previous + 1))
            start = position
        previous = position
    parts.append((start, previous + 1))
    return tuple(parts)


def _first_phase(model: CodingModel) -> int:
    cds_records = {
        (record.start0, record.end0): record
        for record in model.child_records
        if record.feature_type == "CDS"
    }
    first_part = sorted(model.cds_parts, reverse=model.strand == "-")[0]
    try:
        value = int(cds_records[first_part].phase)
    except (KeyError, TypeError, ValueError):
        return 0
    return value if value in {0, 1, 2} else 0


def _path_has_prohibited_ambiguity(
    path: _SplicedPath,
    fasta: IndexedFasta,
    chrom: str,
    start: int,
    end: int,
) -> bool:
    previous_n = False
    for _offset, sequence in path.sequence_chunks(fasta, chrom, start, end):
        for base in sequence:
            if base not in "ACGTN":
                return True
            if base == "N":
                if previous_n:
                    return True
                previous_n = True
            else:
                previous_n = False
    return False


def _upstream_start_proposal(
    model: CodingModel,
    candidate: UTRCandidateSide,
    fasta: IndexedFasta,
    selection: object,
) -> tuple[_Proposal | None, str, bool, int, bool]:
    """Find the natural upstream ATG on the preferred supported path.

    The selected structure is finite, so its biological boundary replaces the
    former arbitrary 100-bp discovery window.  A preceding in-frame stop or a
    contiguous ``NN+`` run can shorten the accessible open segment.  The first
    ATG in 5-prime-to-3-prime order within that segment is the only proposal;
    nearer internal ATGs are never substituted if that proposal later fails a
    safety check.
    """
    evidence_path = _SplicedPath.from_parts(
        candidate.evidence.exon_parts, model.strand
    )
    current_path = _SplicedPath.from_parts(model.cds_parts, model.strand)
    if current_path.length == 0:
        return None, "candidate_path_does_not_contain_cds", False, 0, False
    boundary_index = evidence_path.offset_of_position(current_path.genomic_position(0))
    if boundary_index is None:
        return None, "candidate_path_does_not_contain_cds", False, 0, False
    prefix_parts = evidence_path.slice_parts(0, boundary_index)
    path = _SplicedPath.from_parts((*prefix_parts, *model.cds_parts), model.strand)
    current_start = sum(end0 - start0 for start0, end0 in prefix_parts)
    translated_start = current_start + _first_phase(model)
    first_aligned = translated_start % 3
    segment_start = first_aligned
    chosen: int | None = None
    saw_stop = False
    previous_n = False
    codon_chars: list[str] = []
    codon_start = first_aligned
    for chunk_offset, sequence in path.sequence_chunks(
        fasta, model.chrom, first_aligned, translated_start
    ):
        for chunk_index, base in enumerate(sequence):
            offset = chunk_offset + chunk_index
            if base == "N":
                if previous_n:
                    raw_boundary = offset + 1
                    segment_start = min(
                        translated_start,
                        first_aligned
                        + max(0, (raw_boundary - first_aligned + 2) // 3) * 3,
                    )
                    chosen = None
                    saw_stop = False
                previous_n = True
            else:
                previous_n = False
            if not codon_chars:
                codon_start = offset
            codon_chars.append(base)
            if len(codon_chars) != 3:
                continue
            if codon_start >= segment_start:
                codon = "".join(codon_chars)
                if codon in STANDARD_STOP_CODONS:
                    saw_stop = True
                    segment_start = codon_start + 3
                    chosen = None
                elif codon == "ATG" and chosen is None:
                    chosen = codon_start
            codon_chars.clear()
    accessible = translated_start - segment_start
    if chosen is None:
        return None, "upstream_start_not_found", not saw_stop, accessible, False
    if _path_has_prohibited_ambiguity(
        path, fasta, model.chrom, chosen, current_start
    ):
        return (
            None,
            "coding_extension_crosses_ambiguous_run",
            not saw_stop,
            accessible,
            True,
        )
    proposal = _Proposal(
        side="five_prime",
        direction="upstream_start",
        cds_parts=path.slice_parts(chosen, path.length),
        residual_utr_parts=path.slice_parts(0, chosen),
        codon="ATG",
        codon_positions=tuple(
            path.genomic_position(chosen + index) for index in range(3)
        ),
        bp_delta=current_start - chosen,
        bases_searched=accessible,
        splice_junctions_crossed=path.junction_count(chosen, current_start),
        orf_open_at_search_boundary=not saw_stop,
        selection=selection,
    )
    return proposal, ".", not saw_stop, accessible, True


def _downstream_start_proposal(
    model: CodingModel,
    fasta: IndexedFasta,
    limit: int,
    selection: object | None,
) -> tuple[_Proposal | None, str]:
    path = _SplicedPath.from_parts(model.cds_parts, model.strand)
    translated_start = _first_phase(model)
    maximum = min(limit, max(0, path.length - translated_start - 3))
    maximum -= maximum % 3
    for offset in range(3, maximum + 1, 3):
        codon_start = translated_start + offset
        if (
            path.sequence(fasta, model.chrom, codon_start, codon_start + 3)
            != "ATG"
        ):
            continue
        residual: tuple[Interval, ...] = ()
        candidate = getattr(selection, "candidate", None)
        if candidate is not None:
            evidence_path = _SplicedPath.from_parts(
                candidate.evidence.exon_parts, model.strand
            )
            evidence_codon_start = evidence_path.offset_of_position(
                path.genomic_position(codon_start)
            )
            if evidence_codon_start is not None:
                residual = evidence_path.slice_parts(0, evidence_codon_start)
        return _Proposal(
            side="five_prime",
            direction="downstream_start_fallback",
            cds_parts=path.slice_parts(codon_start, path.length),
            residual_utr_parts=residual,
            codon="ATG",
            codon_positions=tuple(
                path.genomic_position(codon_start + index) for index in range(3)
            ),
            bp_delta=-codon_start,
            bases_searched=offset,
            splice_junctions_crossed=path.junction_count(0, codon_start),
            orf_open_at_search_boundary=None,
            selection=selection,
        ), "."
    return None, "downstream_start_not_found"


def _downstream_stop_proposal(
    model: CodingModel,
    candidate: UTRCandidateSide,
    fasta: IndexedFasta,
    limit: int,
    selection: object,
) -> tuple[_Proposal | None, str]:
    evidence_path = _SplicedPath.from_parts(candidate.evidence.exon_parts, model.strand)
    current_path = _SplicedPath.from_parts(model.cds_parts, model.strand)
    if current_path.length == 0:
        return None, "candidate_path_does_not_contain_cds"
    boundary_index = evidence_path.offset_of_position(
        current_path.genomic_position(current_path.length - 1)
    )
    if boundary_index is None:
        return None, "candidate_path_does_not_contain_cds"
    suffix_parts = evidence_path.slice_parts(boundary_index + 1, evidence_path.length)
    path = _SplicedPath.from_parts((*model.cds_parts, *suffix_parts), model.strand)
    current_end = current_path.length - 1
    translated_start = _first_phase(model)
    for codon_start in range(translated_start, path.length - 2, 3):
        codon_end = codon_start + 3
        if codon_end <= current_end + 1:
            continue
        added = codon_end - (current_end + 1)
        if added > limit:
            break
        codon = path.sequence(fasta, model.chrom, codon_start, codon_end)
        if codon not in STANDARD_STOP_CODONS:
            continue
        if _path_has_prohibited_ambiguity(
            path, fasta, model.chrom, current_end + 1, codon_end
        ):
            return None, "coding_extension_crosses_ambiguous_run"
        return _Proposal(
            side="three_prime",
            direction="downstream_stop",
            cds_parts=path.slice_parts(0, codon_end),
            residual_utr_parts=path.slice_parts(codon_end, path.length),
            codon=codon,
            codon_positions=tuple(
                path.genomic_position(codon_start + index) for index in range(3)
            ),
            bp_delta=added,
            bases_searched=added,
            splice_junctions_crossed=path.junction_count(current_end, codon_end - 1),
            orf_open_at_search_boundary=None,
            selection=selection,
        ), "."
    return None, "downstream_stop_not_found"


def _synthetic_record(
    template: FeatureRecord,
    feature_type: str,
    part: Interval,
    phase: str,
) -> FeatureRecord:
    return replace(
        template,
        feature_type=feature_type,
        start=part[0] + 1,
        end=part[1],
        phase=phase,
        attrs={key: value for key, value in template.attrs.items() if key != "ID"},
        raw_line="",
    )


def _subtract_intervals(
    interval: Interval, blockers: Sequence[Interval]
) -> tuple[Interval, ...]:
    pieces = [interval]
    for block_start, block_end in sorted(blockers):
        next_pieces: list[Interval] = []
        for start0, end0 in pieces:
            if block_end <= start0 or end0 <= block_start:
                next_pieces.append((start0, end0))
                continue
            if start0 < block_start:
                next_pieces.append((start0, block_start))
            if block_end < end0:
                next_pieces.append((block_end, end0))
        pieces = next_pieces
    return tuple(piece for piece in pieces if piece[1] > piece[0])


def _consistent_source_utr_records(
    model: CodingModel,
    parts: Sequence[Interval],
) -> list[FeatureRecord]:
    """Preserve only source UTR fragments that remain outside the refined CDS."""
    cds_start0 = min(start0 for start0, _end0 in parts)
    cds_end0 = max(end0 for _start0, end0 in parts)
    records: list[FeatureRecord] = []
    for record in model.child_records:
        feature = record.feature_type.lower()
        if feature not in {"utr", "five_prime_utr", "three_prime_utr"}:
            continue
        for start0, end0 in _subtract_intervals(
            (record.start0, record.end0), parts
        ):
            if feature == "five_prime_utr":
                consistent = (
                    end0 <= cds_start0 if model.strand == "+" else start0 >= cds_end0
                )
            elif feature == "three_prime_utr":
                consistent = (
                    start0 >= cds_end0 if model.strand == "+" else end0 <= cds_start0
                )
            else:
                consistent = True
            if consistent:
                records.append(
                    replace(
                        record,
                        start=start0 + 1,
                        end=end0,
                        attrs={
                            key: value
                            for key, value in record.attrs.items()
                            if key != "ID"
                        },
                        raw_line="",
                    )
                )
    return records


def _codon_records(
    model: CodingModel,
    template: FeatureRecord,
    parts: Sequence[Interval],
    first_phase: int,
) -> list[FeatureRecord]:
    """Rebuild codon children at the current biological CDS termini."""
    existing = {record.feature_type.lower() for record in model.child_records}
    path = _SplicedPath.from_parts(parts, model.strand)
    records: list[FeatureRecord] = []
    for feature_type, codon_offsets in (
        ("start_codon", range(first_phase, min(first_phase + 3, path.length))),
        ("stop_codon", range(max(0, path.length - 3), path.length)),
    ):
        offsets = tuple(codon_offsets)
        if feature_type not in existing or len(offsets) != 3:
            continue
        codon_positions = tuple(path.genomic_position(offset) for offset in offsets)
        records.extend(
            _synthetic_record(template, feature_type, part, "0")
            for part in _parts_from_positions(codon_positions)
        )
    return records


def _refined_model(
    model: CodingModel,
    parts: Sequence[Interval],
    *,
    start_completed: bool,
) -> CodingModel:
    ordered = tuple(sorted(parts))
    template = next(
        (record for record in model.child_records if record.feature_type == "CDS"),
        model.child_records[0] if model.child_records else None,
    )
    if template is None:
        template = FeatureRecord(
            seqid=model.chrom,
            source=model.source_name,
            feature_type="CDS",
            start=ordered[0][0] + 1,
            end=ordered[0][1],
            score=".",
            strand=model.strand,
            phase="0",
            attrs={},
            raw_attrs=".",
            input_file=model.input_file,
            source_name=model.source_name,
            line_number=0,
            record_index=0,
            raw_line="",
        )
    first_phase = 0 if start_completed else _first_phase(model)
    phases = calculate_cds_phase_chain(ordered, model.strand, first_phase)
    source_exons = [
        (record.start0, record.end0)
        for record in model.child_records
        if record.feature_type.lower() == "exon"
    ]
    exon_parts = merge_intervals(source_exons + list(ordered))
    # Final children use an explicit supported set.  Arbitrary source
    # positional features (for example intron/TSS/TES records) are not copied
    # because their geometry can become stale after refinement.
    children: list[FeatureRecord] = []
    children.extend(_consistent_source_utr_records(model, ordered))
    children.extend(_codon_records(model, template, ordered, first_phase))
    children.extend(
        _synthetic_record(template, "exon", part, ".") for part in exon_parts
    )
    children.extend(
        _synthetic_record(template, "CDS", part, phases[part]) for part in ordered
    )
    return replace(
        model,
        cds_parts=list(ordered),
        child_records=children,
        sequence_integrity=None,
    )


def _combined_refined_model(
    original: CodingModel,
    proposals: Sequence[_Proposal],
) -> CodingModel:
    """Rebuild one model from independent five- and three-prime proposals."""
    by_side = {proposal.side: proposal for proposal in proposals}
    parts = (
        list(by_side["five_prime"].cds_parts)
        if "five_prime" in by_side
        else list(original.cds_parts)
    )
    if "three_prime" in by_side:
        stop_path = _SplicedPath.from_parts(
            by_side["three_prime"].cds_parts, original.strand
        )
        original_path = _SplicedPath.from_parts(original.cds_parts, original.strand)
        if (
            original_path.length == 0
            or stop_path.offset_of_position(
                original_path.genomic_position(original_path.length - 1)
            )
            is None
        ):
            raise ValueError("Stop proposal does not contain the original CDS terminus")
        added_stop_parts: list[Interval] = []
        for stop_part in by_side["three_prime"].cds_parts:
            added_stop_parts.extend(_subtract_intervals(stop_part, original.cds_parts))
        parts.extend(added_stop_parts)
    refined = _refined_model(
        original,
        merge_intervals(parts),
        start_completed="five_prime" in by_side,
    )
    if "five_prime" in by_side and "three_prime" in by_side:
        marker = "start_and_stop"
    elif "five_prime" in by_side:
        marker = (
            "upstream_start"
            if by_side["five_prime"].direction == "upstream_start"
            else "downstream_start"
        )
    else:
        marker = "stop"
    refined.cds_terminus_refined = marker
    return refined


def _integrity_row(
    model: CodingModel,
    fasta: IndexedFasta,
    settings: CDSTerminusRefinementSettings,
) -> dict[str, object]:
    contig_length = fasta.contig_length(model.chrom)
    if contig_length is None:
        raise ValueError(f"Genome FASTA is missing contig {model.chrom!r}")
    start0 = max(0, model.cds_start0 - 3)
    end0 = min(contig_length, model.cds_end0 + 3)
    return calculate_model_sequence_integrity(
        model,
        window=fasta.fetch(model.chrom, start0, end0),
        window_start0=start0,
        contig_length=contig_length,
        genetic_code=settings.genetic_code,
        internal_stop_penalty=settings.internal_stop_penalty,
        internal_stop_count_cap=settings.internal_stop_count_cap,
        phase_inconsistent_penalty=settings.phase_inconsistent_penalty,
        missing_stop_penalty=settings.missing_stop_penalty,
        missing_start_penalty=settings.missing_start_penalty,
        noncanonical_splice_penalty=settings.noncanonical_splice_penalty,
    )


def _validate_refined_model(
    original: CodingModel,
    proposals: Sequence[_Proposal],
    gene_id: str,
    frozen_relationships: _FrozenRelationshipIndex,
    fasta: IndexedFasta,
    settings: CDSTerminusRefinementSettings,
) -> tuple[CodingModel | None, str, str, str]:
    """Validate one independently assembled proposal set."""
    if not proposals:
        return None, "completion_not_found", "not_run", "."
    refined = _combined_refined_model(original, proposals)
    if frozen_relationships.has_new_same_strand_cross_gene_cds_overlap(
        original,
        refined,
        gene_id,
    ):
        return None, "new_same_strand_cross_gene_cds_overlap", "fail", "."
    if frozen_relationships.has_new_relationship(
        original,
        refined,
        gene_id,
    ):
        return None, "frozen_final_gene_structural_conflict", "fail", "."
    if not _new_introns_are_canonical(original, refined, fasta):
        return None, "new_noncanonical_cds_splice", "pass", "."
    row = _integrity_row(refined, fasta, settings)
    integrity_result = str(row["sequence_integrity_status"])
    completed_sides = {proposal.side for proposal in proposals}
    invalid = (
        int(row["sequence_internal_stop_count"]) > 0
        or str(row["sequence_phase_status"]) == "inconsistent"
        or (
            "five_prime" in completed_sides
            and str(row["sequence_start_codon_status"]) != "present"
        )
        or (
            "three_prime" in completed_sides
            and str(row["sequence_stop_codon_status"]) != "included_in_cds"
        )
    )
    if invalid:
        return None, "refined_cds_integrity_failed", "pass", integrity_result
    metrics = sequence_integrity_metrics_from_row(row)
    refined.sequence_integrity = metrics
    refined.child_records = [
        record
        for record in refined.child_records
        if not (
            record.feature_type.lower() == "start_codon"
            and metrics.start_codon_status != "present"
        )
        and not (
            record.feature_type.lower() == "stop_codon"
            and metrics.stop_codon_status != "included_in_cds"
        )
    ]
    return refined, ".", "pass", integrity_result


def _new_introns_are_canonical(
    original: CodingModel, refined: CodingModel, fasta: IndexedFasta
) -> bool:
    added = set(refined.cds_introns) - set(original.cds_introns)
    for start0, end0 in added:
        if refined.strand == "+":
            donor = fasta.fetch(refined.chrom, start0, min(start0 + 2, end0))
            acceptor = fasta.fetch(refined.chrom, max(start0, end0 - 2), end0)
        else:
            donor = reverse_complement(
                fasta.fetch(refined.chrom, max(start0, end0 - 2), end0)
            )
            acceptor = reverse_complement(
                fasta.fetch(refined.chrom, start0, min(start0 + 2, end0))
            )
        if f"{donor}-{acceptor}" not in CANONICAL_SPLICE_MOTIFS:
            return False
    return True


@dataclass(frozen=True, slots=True)
class _FinalModelEntry:
    gene_id: str
    transcript_id: str
    model: CodingModel


def _intervals_overlap(left: Interval, right: Interval) -> bool:
    return max(left[0], right[0]) < min(left[1], right[1])


def _new_cds_parts(
    original: CodingModel, refined: CodingModel
) -> tuple[Interval, ...]:
    """Return only coding intervals introduced by refinement."""
    pieces: list[Interval] = []
    for part in refined.cds_parts:
        pieces.extend(_subtract_intervals(part, original.cds_parts))
    return tuple(merge_intervals(pieces))


class _FrozenRelationshipIndex:
    """Spatially indexed original final models for canonical relationship checks."""

    def __init__(self, annotation: FinalAnnotation, threshold: float) -> None:
        self.threshold = threshold
        self.entries = tuple(
            _FinalModelEntry(gene.gene_id, transcript.transcript_id, transcript.source_model)
            for gene in annotation.genes
            for transcript in gene.transcripts
        )
        self.by_gene: dict[str, tuple[_FinalModelEntry, ...]] = {
            gene.gene_id: tuple(
                _FinalModelEntry(
                    gene.gene_id, transcript.transcript_id, transcript.source_model
                )
                for transcript in gene.transcripts
            )
            for gene in annotation.genes
        }
        self.overlaps = IntervalOverlapIndex(
            self.entries,
            key=lambda entry: (entry.model.chrom, entry.model.strand),
            start=lambda entry: entry.model.cds_start0,
            end=lambda entry: entry.model.cds_end0,
        )
        self._baseline_gene_pair_cache: dict[tuple[str, str], bool] = {}

    def _related(self, left: CodingModel, right: CodingModel) -> bool:
        return models_share_sublocus_relationship(
            left,
            right,
            single_exon_min_reciprocal_overlap=self.threshold,
        )

    def has_new_relationship(
        self,
        original: CodingModel,
        refined: CodingModel,
        own_gene_id: str,
    ) -> bool:
        for entry in self.overlaps.query(
            (refined.chrom, refined.strand), refined.cds_start0, refined.cds_end0
        ):
            if entry.gene_id == own_gene_id:
                continue
            if self._related(refined, entry.model) and not self._related(
                original, entry.model
            ):
                return True
        return False

    def has_new_same_strand_cross_gene_cds_overlap(
        self,
        original: CodingModel,
        refined: CodingModel,
        own_gene_id: str,
    ) -> bool:
        """Check newly coding bases against frozen same-strand CDS geometry."""
        for new_part in _new_cds_parts(original, refined):
            for entry in self.overlaps.query(
                (refined.chrom, refined.strand), new_part[0], new_part[1]
            ):
                if entry.gene_id == own_gene_id:
                    continue
                if any(
                    _intervals_overlap(new_part, other_part)
                    for other_part in entry.model.cds_parts
                ):
                    return True
        return False

    def baseline_gene_pair_related(self, left_gene: str, right_gene: str) -> bool:
        pair = tuple(sorted((left_gene, right_gene)))
        cached = self._baseline_gene_pair_cache.get(pair)
        if cached is not None:
            return cached
        related = any(
            self._related(left.model, right.model)
            for left in self.by_gene[left_gene]
            for right in self.by_gene[right_gene]
        )
        self._baseline_gene_pair_cache[pair] = related
        return related


def _simultaneous_conflicting_refinements(
    annotation: FinalAnnotation,
    replacements: dict[str, CodingModel],
    frozen: _FrozenRelationshipIndex,
) -> set[str]:
    """Find refined participants in new relationships using spatial neighbours."""
    entries = tuple(
        _FinalModelEntry(
            gene.gene_id,
            transcript.transcript_id,
            replacements.get(transcript.transcript_id, transcript.source_model),
        )
        for gene in annotation.genes
        for transcript in gene.transcripts
    )
    by_transcript = {entry.transcript_id: entry for entry in entries}
    overlaps = IntervalOverlapIndex(
        entries,
        key=lambda entry: (entry.model.chrom, entry.model.strand),
        start=lambda entry: entry.model.cds_start0,
        end=lambda entry: entry.model.cds_end0,
    )
    conflicting: set[str] = set()
    checked_pairs: set[tuple[str, str]] = set()
    for transcript_id in sorted(replacements):
        entry = by_transcript[transcript_id]
        for other in overlaps.query(
            (entry.model.chrom, entry.model.strand),
            entry.model.cds_start0,
            entry.model.cds_end0,
        ):
            if other.gene_id == entry.gene_id:
                continue
            transcript_pair = tuple(sorted((entry.transcript_id, other.transcript_id)))
            if transcript_pair in checked_pairs:
                continue
            checked_pairs.add(transcript_pair)
            if not frozen._related(entry.model, other.model):
                continue
            if frozen.baseline_gene_pair_related(entry.gene_id, other.gene_id):
                continue
            conflicting.add(entry.transcript_id)
            if other.transcript_id in replacements:
                conflicting.add(other.transcript_id)
    return conflicting


def _simultaneous_new_cds_overlap_refinements(
    annotation: FinalAnnotation,
    replacements: dict[str, CodingModel],
) -> set[str]:
    """Find refinements whose new CDS overlaps another final gene's CDS.

    Frozen-CDS overlaps are rejected during individual validation.  This pass
    covers the remaining refined-vs-refined case where two extensions enter
    the same previously non-coding interval.  Both contributors are rolled
    back, independent of annotation or proposal order.
    """
    original_by_transcript = {
        transcript.transcript_id: transcript.source_model
        for gene in annotation.genes
        for transcript in gene.transcripts
    }
    entries = tuple(
        _FinalModelEntry(
            gene.gene_id,
            transcript.transcript_id,
            replacements.get(transcript.transcript_id, transcript.source_model),
        )
        for gene in annotation.genes
        for transcript in gene.transcripts
    )
    entry_by_transcript = {entry.transcript_id: entry for entry in entries}
    overlaps = IntervalOverlapIndex(
        entries,
        key=lambda entry: (entry.model.chrom, entry.model.strand),
        start=lambda entry: entry.model.cds_start0,
        end=lambda entry: entry.model.cds_end0,
    )
    new_by_transcript = {
        transcript_id: _new_cds_parts(
            original_by_transcript[transcript_id], replacement
        )
        for transcript_id, replacement in replacements.items()
    }
    conflicting: set[str] = set()
    for transcript_id in sorted(replacements):
        entry = entry_by_transcript[transcript_id]
        for new_part in new_by_transcript[transcript_id]:
            for other in overlaps.query(
                (entry.model.chrom, entry.model.strand), new_part[0], new_part[1]
            ):
                if other.gene_id == entry.gene_id:
                    continue
                if not any(
                    _intervals_overlap(new_part, other_part)
                    for other_part in other.model.cds_parts
                ):
                    continue
                conflicting.add(transcript_id)
                if (
                    other.transcript_id in replacements
                    and any(
                        _intervals_overlap(new_part, other_new)
                        for other_new in new_by_transcript[other.transcript_id]
                    )
                ):
                    conflicting.add(other.transcript_id)
    return conflicting


def _audit_row(
    gene_id: str,
    transcript: FinalTranscript,
    side: Side,
    original_status: str,
    proposal: _Proposal | None,
    *,
    result: str,
    reason: str,
    original_length: int,
    refined_model: CodingModel | None,
    settings: CDSTerminusRefinementSettings,
    structural_check: str = "not_run",
    integrity_result: str = ".",
    selection_override: object | None = None,
    direction_override: str = ".",
    orf_open_override: bool | None = None,
    bases_searched_override: int = 0,
    upstream_failure_reason: str = ".",
    ambiguity_encountered: bool = False,
    encountered_rejections: Counter[str] | None = None,
) -> dict[str, object]:
    model = transcript.source_model
    selection = proposal.selection if proposal is not None else selection_override
    candidate = getattr(selection, "candidate", None)
    refined_length = refined_model.cds_length if refined_model is not None else original_length
    delta = proposal.bp_delta if proposal is not None else 0
    direction = proposal.direction if proposal else direction_override
    limit: int | str = (
        "."
        if direction == "upstream_start"
        else settings.downstream_start_search_bp
        if side == "five_prime"
        else settings.downstream_stop_search_bp
    )
    original_short = original_length < settings.min_cds_length_bp
    refined_short = refined_length < settings.min_cds_length_bp
    original_metrics = model.sequence_integrity
    refined_metrics = refined_model.sequence_integrity if refined_model is not None else None
    rejection_counts = encountered_rejections or Counter()
    def phase_text(value: CodingModel | None) -> str:
        if value is None:
            return "."
        phases = {
            (record.start0, record.end0): record.phase
            for record in value.child_records
            if record.feature_type == "CDS"
        }
        return ",".join(
            phases.get(part, ".")
            for part in sorted(value.cds_parts, reverse=value.strand == "-")
        ) or "."
    return {
        "gene_id": gene_id,
        "transcript_id": transcript.transcript_id,
        "source_name": model.source_name,
        "internal_transcript_id": model.output_transcript_id,
        "refinement_attempted": "true",
        "side": side,
        "original_status": original_status,
        "candidate_structure_key": candidate.structure_key if candidate else ".",
        "candidate_evidence_id": candidate.evidence.evidence_id if candidate else ".",
        "candidate_evidence_kind": candidate.evidence.kind if candidate else ".",
        "candidate_sources": (
            ",".join(getattr(selection, "selected_sources", ())) or "."
        ),
        "candidate_rank": 1 if candidate else ".",
        "candidate_effective_support": (
            f"{selection.effective_support:.6f}" if selection is not None else "."
        ),
        "completion_direction": proposal.direction if proposal else direction_override,
        "upstream_completion_failure_reason": upstream_failure_reason,
        "configured_search_limit_bp": limit,
        "transcript_bases_searched": (
            proposal.bases_searched if proposal else bases_searched_override
        ),
        "orf_open_at_search_boundary": (
            str(proposal.orf_open_at_search_boundary).lower()
            if proposal and proposal.orf_open_at_search_boundary is not None
            else str(orf_open_override).lower()
            if orf_open_override is not None
            else "."
        ),
        "codon_found": proposal.codon if proposal else ".",
        "codon_genomic_positions": (
            ",".join(str(position + 1) for position in proposal.codon_positions)
            if proposal else "."
        ),
        "codon_spliced_search_offset_bp": (
            proposal.bases_searched if proposal else "."
        ),
        "original_cds_parts": ",".join(f"{s + 1}-{e}" for s, e in model.cds_parts),
        "refined_cds_parts": (
            ",".join(f"{s + 1}-{e}" for s, e in refined_model.cds_parts)
            if refined_model is not None else "."
        ),
        "original_cds_length": original_length,
        "refined_cds_length": refined_length,
        "original_cds_phases": phase_text(model),
        "refined_cds_phases": phase_text(refined_model),
        "original_integrity_status": original_metrics.status if original_metrics else ".",
        "refined_integrity_status": refined_metrics.status if refined_metrics else integrity_result,
        "original_start_status": original_metrics.start_codon_status if original_metrics else ".",
        "refined_start_status": refined_metrics.start_codon_status if refined_metrics else ".",
        "original_stop_status": original_metrics.stop_codon_status if original_metrics else ".",
        "refined_stop_status": refined_metrics.stop_codon_status if refined_metrics else ".",
        "original_internal_stop_count": original_metrics.internal_stop_count if original_metrics else ".",
        "refined_internal_stop_count": refined_metrics.internal_stop_count if refined_metrics else ".",
        "coding_bp_added": max(0, delta),
        "coding_bp_removed": max(0, -delta),
        "relative_change": f"{abs(delta) / original_length:.6f}" if original_length else ".",
        "splice_junctions_crossed": proposal.splice_junctions_crossed if proposal else 0,
        "ambiguous_run_encountered": str(
            ambiguity_encountered
            or reason == "coding_extension_crosses_ambiguous_run"
        ).lower(),
        "original_short_cds": str(original_short).lower(),
        "refined_short_cds": str(refined_short).lower(),
        "crossed_short_cds_threshold": str(original_short != refined_short).lower(),
        "final_gene_structural_check": structural_check,
        "encountered_rejection_event_count": sum(rejection_counts.values()),
        "encountered_rejection_counts": json.dumps(
            dict(sorted(rejection_counts.items())),
            sort_keys=True,
            separators=(",", ":"),
        ),
        "integrity_result": integrity_result,
        "residual_utr_parts": (
            ",".join(f"{s + 1}-{e}" for s, e in proposal.residual_utr_parts)
            if proposal and proposal.residual_utr_parts else "."
        ),
        "residual_utr_exonic_bp": (
            sum(end0 - start0 for start0, end0 in proposal.residual_utr_parts)
            if proposal else 0
        ),
        "result": result,
        "failure_reason": reason,
        "final_outcome": "refined" if result == "refined" else reason,
    }


def _add_encountered_rejection_to_audit_row(
    row: dict[str, object], reason: str, count: int = 1
) -> None:
    values = Counter(json.loads(str(row["encountered_rejection_counts"])))
    values[reason] += count
    row["encountered_rejection_event_count"] = sum(values.values())
    row["encountered_rejection_counts"] = json.dumps(
        dict(sorted(values.items())), sort_keys=True, separators=(",", ":")
    )


def _summary_from_audit_rows(
    rows: Sequence[dict[str, object]],
    *,
    execution_status: str,
    expected_sides: set[tuple[str, Side]] | None = None,
) -> CDSTerminusRefinementSummary:
    """Derive unit-explicit terminal and encountered-event accounting."""
    side_keys = [(str(row["transcript_id"]), str(row["side"])) for row in rows]
    if len(side_keys) != len(set(side_keys)):
        raise AssertionError("CDS refinement audit contains duplicate transcript sides")
    if expected_sides is not None and set(side_keys) != expected_sides:
        raise AssertionError("CDS refinement audit does not cover every incomplete side")

    final_outcomes: Counter[str] = Counter()
    failures: Counter[str] = Counter()
    failure_transcripts: dict[str, set[str]] = defaultdict(set)
    event_counts: Counter[str] = Counter()
    event_sides: Counter[str] = Counter()
    event_transcripts: dict[str, set[str]] = defaultdict(set)
    rows_by_transcript: dict[str, list[dict[str, object]]] = defaultdict(list)
    refined_with_events = 0

    for row in rows:
        transcript_id = str(row["transcript_id"])
        rows_by_transcript[transcript_id].append(row)
        outcome = (
            "refined"
            if row["result"] == "refined"
            else str(row["failure_reason"])
        )
        if outcome == ".":
            raise AssertionError("Failed CDS refinement side lacks a final reason")
        row["final_outcome"] = outcome
        final_outcomes[outcome] += 1
        if outcome != "refined":
            failures[outcome] += 1
            failure_transcripts[outcome].add(transcript_id)

        encountered = Counter(
            json.loads(str(row["encountered_rejection_counts"]))
        )
        if row["result"] == "refined" and encountered:
            refined_with_events += 1
        for reason, count in encountered.items():
            if count <= 0:
                raise AssertionError("Encountered rejection counts must be positive")
            event_counts[reason] += count
            event_sides[reason] += 1
            event_transcripts[reason].add(transcript_id)

    final_refined = final_outcomes["refined"]
    final_failed = len(rows) - final_refined
    if sum(final_outcomes.values()) != len(rows):
        raise AssertionError("Final CDS refinement outcomes do not reconcile")
    if sum(failures.values()) != final_failed:
        raise AssertionError("Final CDS refinement failures do not reconcile")

    two_side_transcripts = {
        transcript_id: transcript_rows
        for transcript_id, transcript_rows in rows_by_transcript.items()
        if len(transcript_rows) == 2
    }
    both_refined = start_only = stop_only = neither = 0
    for transcript_rows in two_side_transcripts.values():
        refined_sides = {
            str(row["side"])
            for row in transcript_rows
            if row["result"] == "refined"
        }
        if refined_sides == {"five_prime", "three_prime"}:
            both_refined += 1
        elif refined_sides == {"five_prime"}:
            start_only += 1
        elif refined_sides == {"three_prime"}:
            stop_only += 1
        else:
            neither += 1

    refined_transcripts = sum(
        any(row["result"] == "refined" for row in transcript_rows)
        for transcript_rows in rows_by_transcript.values()
    )
    return CDSTerminusRefinementSummary(
        execution_status=execution_status,
        audit_sides=len(rows),
        incomplete_five_prime_sides=sum(
            row["side"] == "five_prime" for row in rows
        ),
        incomplete_three_prime_sides=sum(
            row["side"] == "three_prime" for row in rows
        ),
        transcripts_with_two_incomplete_sides=len(two_side_transcripts),
        final_refined_sides=final_refined,
        final_failed_sides=final_failed,
        refined_sides_upstream_start=sum(
            row["result"] == "refined"
            and row["completion_direction"] == "upstream_start"
            for row in rows
        ),
        refined_sides_downstream_start_fallback=sum(
            row["result"] == "refined"
            and row["completion_direction"] == "downstream_start_fallback"
            for row in rows
        ),
        refined_sides_downstream_stop=sum(
            row["result"] == "refined"
            and row["completion_direction"] == "downstream_stop"
            for row in rows
        ),
        refined_transcripts=refined_transcripts,
        unchanged_transcripts=len(rows_by_transcript) - refined_transcripts,
        two_incomplete_sides_both_refined_transcripts=both_refined,
        two_incomplete_sides_start_only_refined_transcripts=start_only,
        two_incomplete_sides_stop_only_refined_transcripts=stop_only,
        two_incomplete_sides_neither_refined_transcripts=neither,
        refined_sides_with_encountered_rejection_events=refined_with_events,
        final_outcome_side_counts=tuple(sorted(final_outcomes.items())),
        failure_counts=tuple(sorted(failures.items())),
        failure_transcript_counts=tuple(
            sorted(
                (reason, len(transcript_ids))
                for reason, transcript_ids in failure_transcripts.items()
            )
        ),
        encountered_rejection_event_counts=tuple(sorted(event_counts.items())),
        encountered_rejection_side_counts=tuple(sorted(event_sides.items())),
        encountered_rejection_transcript_counts=tuple(
            sorted(
                (reason, len(transcript_ids))
                for reason, transcript_ids in event_transcripts.items()
            )
        ),
    )


def refine_cds_termini(
    annotation: FinalAnnotation,
    transcript_evidence: Sequence[UTREvidence],
    prediction_evidence: Sequence[UTREvidence],
    fasta: IndexedFasta,
    settings: CDSTerminusRefinementSettings,
    *,
    utr_settings: UTRSettings | None = None,
    learned_utr_limits: UTRLearnedLimits | None = None,
) -> CDSTerminusRefinementResult:
    """Refine explicitly incomplete termini without changing selection provenance."""
    incomplete: list[tuple[FinalGene, FinalTranscript, tuple[Side, ...]]] = []
    originally_incomplete: set[tuple[str, Side]] = set()
    for gene in annotation.genes:
        for transcript in gene.transcripts:
            metrics = transcript.source_model.sequence_integrity
            if metrics is None:
                continue
            sides: list[Side] = []
            if metrics.start_codon_status == "missing":
                sides.append("five_prime")
            if metrics.stop_codon_status == "missing":
                sides.append("three_prime")
            if sides:
                incomplete.append((gene, transcript, tuple(sides)))
                originally_incomplete.update(
                    (transcript.transcript_id, side) for side in sides
                )

    if not incomplete:
        return empty_cds_refinement_result(
            annotation, execution_status="zero_incomplete_models"
        )

    incomplete_ids = {transcript.transcript_id for _gene, transcript, _sides in incomplete}
    evidence = list(transcript_evidence) + list(prediction_evidence)
    # Index against the complete frozen annotation so bridge detection also
    # sees evidence compatible with complete, non-target final genes.
    index = UTREvidenceIndex(evidence, annotation)
    frozen = FrozenFinalGeneStructureIndex(annotation)
    frozen_relationships = _FrozenRelationshipIndex(
        annotation, settings.single_exon_min_reciprocal_overlap
    )
    candidates_by_transcript: dict[str, list[UTRCandidateSide]] = defaultdict(list)
    for gene, transcript, sides in incomplete:
        candidates = candidate_sides_for_transcript(
            gene.gene_id,
            transcript,
            index,
            utr_settings or UTRSettings(),
            candidate_mode="completion",
            requested_sides=sides,
        )
        candidates_by_transcript[transcript.transcript_id].extend(candidates)
    all_candidates = [
        candidate
        for values in candidates_by_transcript.values()
        for candidate in values
    ]
    bridge_candidates = list(all_candidates)
    relevant_evidence_ids = {
        candidate.evidence.evidence_id for candidate in all_candidates
    }
    nearby_complete_ids: set[str] = set()
    for candidate in all_candidates:
        evidence_item = candidate.evidence
        for entry in frozen_relationships.overlaps.query(
            (evidence_item.chrom, evidence_item.strand),
            evidence_item.start0,
            evidence_item.end0,
        ):
            if entry.transcript_id not in incomplete_ids:
                nearby_complete_ids.add(entry.transcript_id)
    final_transcripts = {
        transcript.transcript_id: (gene.gene_id, transcript)
        for gene in annotation.genes
        for transcript in gene.transcripts
    }
    for transcript_id in sorted(nearby_complete_ids):
        gene_id, transcript = final_transcripts[transcript_id]
        bridge_candidates.extend(
            candidate
            for candidate in candidate_sides_for_transcript(
                gene_id,
                transcript,
                index,
                utr_settings or UTRSettings(),
                candidate_mode="utr",
                requested_sides=("five_prime", "three_prime"),
            )
            if candidate.evidence.evidence_id in relevant_evidence_ids
        )
    mark_multi_gene_bridges(bridge_candidates)
    for candidate in all_candidates:
        conflicts = frozen.candidate_reasons(candidate)
        if conflicts:
            candidate.rejection_reasons = tuple(
                sorted(set(candidate.rejection_reasons) | set(conflicts))
            )

    replacements: dict[str, CodingModel] = {}
    accepted_proposals: dict[str, tuple[_Proposal, ...]] = {}
    audit_rows: list[dict[str, object]] = []

    for gene, transcript, sides in incomplete:
        original = transcript.source_model
        original_length = original.cds_length
        working = original
        proposals: list[_Proposal] = []
        rejected_proposals: dict[Side, _Proposal] = {}
        side_failures: dict[Side, str] = {}
        side_selections: dict[Side, object | None] = {}
        side_directions: dict[Side, str] = {}
        side_orf_open: dict[Side, bool | None] = {}
        side_bases_searched: dict[Side, int] = {}
        side_upstream_failure: dict[Side, str] = {}
        side_ambiguity: dict[Side, bool] = {}
        side_events: dict[Side, Counter[str]] = {
            side: Counter() for side in sides
        }
        side_status = {
            "five_prime": original.sequence_integrity.start_codon_status,
            "three_prime": original.sequence_integrity.stop_codon_status,
        }
        if "five_prime" in sides:
            eligible = [
                candidate
                for candidate in candidates_by_transcript[transcript.transcript_id]
                if candidate.side == "five_prime" and not candidate.multi_gene_bridge
            ]
            selection = select_completion_candidate(eligible, utr_settings or UTRSettings())
            side_selections["five_prime"] = selection
            proposal = None
            reason = "no_compatible_completion_structure"
            structural_candidates_rejected = any(
                candidate.multi_gene_bridge
                or any("cross_final_gene" in item for item in candidate.rejection_reasons)
                for candidate in candidates_by_transcript[transcript.transcript_id]
                if candidate.side == "five_prime"
            )
            if selection is None and structural_candidates_rejected:
                reason = "frozen_final_gene_structural_conflict"
                side_events["five_prime"][reason] += sum(
                    candidate.multi_gene_bridge
                    or any(
                        "cross_final_gene" in item
                        for item in candidate.rejection_reasons
                    )
                    for candidate in candidates_by_transcript[
                        transcript.transcript_id
                    ]
                    if candidate.side == "five_prime"
                )
            upstream_atg_identified = False
            if selection is not None:
                side_directions["five_prime"] = "upstream_start"
                (
                    proposal,
                    reason,
                    upstream_open,
                    upstream_searched,
                    upstream_atg_identified,
                ) = _upstream_start_proposal(
                    working,
                    selection.candidate,
                    fasta,
                    selection,
                )
                side_orf_open["five_prime"] = upstream_open
                side_bases_searched["five_prime"] = upstream_searched
                if proposal is not None:
                    relative = (
                        proposal.bp_delta / original_length
                        if original_length else float("inf")
                    )
                    if relative > settings.max_relative_extension:
                        rejected_proposals["five_prime"] = proposal
                        proposal = None
                        reason = "cds_extension_relative_change_exceeded"
            if proposal is None:
                side_upstream_failure["five_prime"] = reason
                upstream_failure_reason = reason
                if selection is not None:
                    side_events["five_prime"][reason] += 1
                if reason == "coding_extension_crosses_ambiguous_run":
                    side_ambiguity["five_prime"] = True
                if not upstream_atg_identified:
                    proposal, fallback_reason = _downstream_start_proposal(
                        working,
                        fasta,
                        settings.downstream_start_search_bp,
                        selection,
                    )
                    side_directions["five_prime"] = "downstream_start_fallback"
                    if proposal is not None:
                        side_bases_searched["five_prime"] = proposal.bases_searched
                        reason = "."
                    else:
                        side_events["five_prime"][fallback_reason] += 1
                        reason = (
                            upstream_failure_reason
                            if upstream_failure_reason
                            == "frozen_final_gene_structural_conflict"
                            else fallback_reason
                        )
            if proposal is not None and proposal.direction == "downstream_start_fallback":
                relative = (
                    abs(proposal.bp_delta) / original_length
                    if original_length else float("inf")
                )
                if relative > settings.max_relative_truncation:
                    reason = "cds_truncation_relative_change_exceeded"
                    rejected_proposals["five_prime"] = proposal
                    proposal = None
                    side_events["five_prime"][reason] += 1
            if proposal is not None:
                proposals.append(proposal)
            else:
                side_failures["five_prime"] = reason

        if "three_prime" in sides:
            eligible = [
                candidate
                for candidate in candidates_by_transcript[transcript.transcript_id]
                if candidate.side == "three_prime" and not candidate.multi_gene_bridge
            ]
            selection = select_completion_candidate(eligible, utr_settings or UTRSettings())
            side_selections["three_prime"] = selection
            side_directions["three_prime"] = (
                "downstream_stop" if selection is not None else "."
            )
            proposal = None
            reason = "no_compatible_completion_structure"
            structural_candidates_rejected = any(
                candidate.multi_gene_bridge
                or any("cross_final_gene" in item for item in candidate.rejection_reasons)
                for candidate in candidates_by_transcript[transcript.transcript_id]
                if candidate.side == "three_prime"
            )
            if selection is None and structural_candidates_rejected:
                reason = "frozen_final_gene_structural_conflict"
                side_events["three_prime"][reason] += sum(
                    candidate.multi_gene_bridge
                    or any(
                        "cross_final_gene" in item
                        for item in candidate.rejection_reasons
                    )
                    for candidate in candidates_by_transcript[
                        transcript.transcript_id
                    ]
                    if candidate.side == "three_prime"
                )
            if selection is not None:
                proposal, reason = _downstream_stop_proposal(
                    original,
                    selection.candidate,
                    fasta,
                    settings.downstream_stop_search_bp,
                    selection,
                )
            if proposal is not None:
                relative = proposal.bp_delta / original_length if original_length else float("inf")
                if relative > settings.max_relative_extension:
                    reason = "cds_extension_relative_change_exceeded"
                    rejected_proposals["three_prime"] = proposal
                    proposal = None
                    side_events["three_prime"][reason] += 1
            if proposal is not None:
                proposals.append(proposal)
            else:
                side_failures["three_prime"] = reason
                if selection is not None and not side_events["three_prime"][reason]:
                    side_events["three_prime"][reason] += 1

        validation_reason = "."
        structural_check = "pending" if proposals else "not_run"
        integrity_result = "."
        if proposals:
            working, validation_reason, structural_check, integrity_result = (
                _validate_refined_model(
                    original,
                    proposals,
                    gene.gene_id,
                    frozen_relationships,
                    fasta,
                    settings,
                )
            )
            if validation_reason != ".":
                for proposal in proposals:
                    side_events[proposal.side][validation_reason] += 1
            if working is None and len(proposals) > 1:
                individually_valid: list[tuple[_Proposal, CodingModel, str, str]] = []
                individual_failures: dict[Side, tuple[str, str]] = {}
                for proposal in proposals:
                    candidate_model, reason, check, integrity = _validate_refined_model(
                        original,
                        (proposal,),
                        gene.gene_id,
                        frozen_relationships,
                        fasta,
                        settings,
                    )
                    if candidate_model is not None:
                        individually_valid.append(
                            (proposal, candidate_model, check, integrity)
                        )
                    else:
                        individual_failures[proposal.side] = (reason, check)
                        side_events[proposal.side][reason] += 1
                if len(individually_valid) == 1:
                    surviving, working, structural_check, integrity_result = (
                        individually_valid[0]
                    )
                    for proposal in proposals:
                        if proposal is surviving:
                            continue
                        reason, _check = individual_failures.get(
                            proposal.side,
                            (validation_reason, structural_check),
                        )
                        side_failures[proposal.side] = reason
                        rejected_proposals[proposal.side] = proposal
                    proposals = [surviving]
                    validation_reason = "."

        if proposals and validation_reason == ".":
            assert working is not None
            replacements[transcript.transcript_id] = working
            accepted_proposals[transcript.transcript_id] = tuple(proposals)
            for proposal in proposals:
                audit_rows.append(
                    _audit_row(
                        gene.gene_id,
                        transcript,
                        proposal.side,
                        side_status[proposal.side],
                        proposal,
                        result="refined",
                        reason=".",
                        original_length=original_length,
                        refined_model=working,
                        settings=settings,
                        structural_check=structural_check,
                        integrity_result=integrity_result,
                        selection_override=side_selections.get(proposal.side),
                        direction_override=side_directions.get(proposal.side, "."),
                        orf_open_override=side_orf_open.get(proposal.side),
                        bases_searched_override=side_bases_searched.get(proposal.side, 0),
                        upstream_failure_reason=side_upstream_failure.get(proposal.side, "."),
                        ambiguity_encountered=side_ambiguity.get(proposal.side, False),
                        encountered_rejections=side_events[proposal.side],
                    )
                )
            for side, reason in side_failures.items():
                audit_rows.append(
                    _audit_row(
                        gene.gene_id,
                        transcript,
                        side,
                        side_status[side],
                        rejected_proposals.get(side),
                        result="unchanged",
                        reason=reason,
                        original_length=original_length,
                        refined_model=working,
                        settings=settings,
                        selection_override=side_selections.get(side),
                        direction_override=side_directions.get(side, "."),
                        orf_open_override=side_orf_open.get(side),
                        bases_searched_override=side_bases_searched.get(side, 0),
                        upstream_failure_reason=side_upstream_failure.get(side, "."),
                        ambiguity_encountered=side_ambiguity.get(side, False),
                        encountered_rejections=side_events[side],
                    )
                )
        else:
            rollback_reason = validation_reason if validation_reason != "." else None
            for side in sides:
                proposal = next(
                    (item for item in proposals if item.side == side),
                    rejected_proposals.get(side),
                )
                reason = rollback_reason or side_failures.get(side, "completion_not_found")
                audit_rows.append(
                    _audit_row(
                        gene.gene_id,
                        transcript,
                        side,
                        side_status[side],
                        proposal,
                        result="rolled_back" if proposal else "unchanged",
                        reason=reason,
                        original_length=original_length,
                        refined_model=None,
                        settings=settings,
                        structural_check=structural_check,
                        integrity_result=integrity_result,
                        selection_override=side_selections.get(side),
                        direction_override=side_directions.get(side, "."),
                        orf_open_override=side_orf_open.get(side),
                        bases_searched_override=side_bases_searched.get(side, 0),
                        upstream_failure_reason=side_upstream_failure.get(side, "."),
                        ambiguity_encountered=side_ambiguity.get(side, False),
                        encountered_rejections=side_events[side],
                    )
                )

    simultaneous_reason = "simultaneous_final_population_structural_conflict"
    overlap_reason = "new_same_strand_cross_gene_cds_overlap"
    while True:
        structural_conflicting = _simultaneous_conflicting_refinements(
            annotation,
            replacements,
            frozen_relationships,
        )
        overlap_conflicting = _simultaneous_new_cds_overlap_refinements(
            annotation,
            replacements,
        )
        conflicting = structural_conflicting | overlap_conflicting
        if not conflicting:
            break
        for transcript_id in sorted(conflicting):
            proposals = accepted_proposals.pop(transcript_id, ())
            if transcript_id not in replacements or not proposals:
                continue
            replacements.pop(transcript_id)
            reason = (
                overlap_reason
                if transcript_id in overlap_conflicting
                else simultaneous_reason
            )
            for row in audit_rows:
                if (
                    row["transcript_id"] == transcript_id
                    and row["result"] == "refined"
                ):
                    row["result"] = "rolled_back"
                    row["failure_reason"] = reason
                    row["final_outcome"] = reason
                    row["final_gene_structural_check"] = "fail"
                    _add_encountered_rejection_to_audit_row(row, reason)

    refined_genes: list[FinalGene] = []
    for gene in annotation.genes:
        transcripts = tuple(
            replace(
                transcript,
                source_model=replacements.get(
                    transcript.transcript_id, transcript.source_model
                ),
            )
            for transcript in gene.transcripts
        )
        refined_genes.append(replace(gene, transcripts=transcripts))
    refined_annotation = FinalAnnotation(tuple(refined_genes))
    return CDSTerminusRefinementResult(
        original_annotation=annotation,
        annotation=refined_annotation,
        audit_rows=tuple(audit_rows),
        summary=_summary_from_audit_rows(
            audit_rows,
            execution_status="performed",
            expected_sides=originally_incomplete,
        ),
        originally_incomplete_sides=frozenset(originally_incomplete),
    )


def write_cds_refinement_audit_tsv(
    result: CDSTerminusRefinementResult, path: Path
) -> None:
    rows = list(result.audit_rows)
    if not rows:
        path.write_text("")
        return
    write_dict_tsv(rows, path, list(rows[0]))


def write_cds_refinement_summary_tsv(
    result: CDSTerminusRefinementResult, path: Path
) -> None:
    summary = result.summary
    row = {
        "execution_status": summary.execution_status,
        "audit_sides": summary.audit_sides,
        "incomplete_five_prime_sides": summary.incomplete_five_prime_sides,
        "incomplete_three_prime_sides": summary.incomplete_three_prime_sides,
        "transcripts_with_two_incomplete_sides": (
            summary.transcripts_with_two_incomplete_sides
        ),
        "final_refined_sides": summary.final_refined_sides,
        "final_failed_sides": summary.final_failed_sides,
        "refined_sides_upstream_start": summary.refined_sides_upstream_start,
        "refined_sides_downstream_start_fallback": (
            summary.refined_sides_downstream_start_fallback
        ),
        "refined_sides_downstream_stop": summary.refined_sides_downstream_stop,
        "refined_transcripts": summary.refined_transcripts,
        "unchanged_transcripts": summary.unchanged_transcripts,
        "two_incomplete_sides_both_refined_transcripts": (
            summary.two_incomplete_sides_both_refined_transcripts
        ),
        "two_incomplete_sides_start_only_refined_transcripts": (
            summary.two_incomplete_sides_start_only_refined_transcripts
        ),
        "two_incomplete_sides_stop_only_refined_transcripts": (
            summary.two_incomplete_sides_stop_only_refined_transcripts
        ),
        "two_incomplete_sides_neither_refined_transcripts": (
            summary.two_incomplete_sides_neither_refined_transcripts
        ),
        "refined_sides_with_encountered_rejection_events": (
            summary.refined_sides_with_encountered_rejection_events
        ),
        "final_outcome_side_counts": json.dumps(
            dict(summary.final_outcome_side_counts),
            sort_keys=True,
            separators=(",", ":"),
        ),
        "failure_counts": json.dumps(
            dict(summary.failure_counts), sort_keys=True, separators=(",", ":")
        ),
        "failure_transcript_counts": json.dumps(
            dict(summary.failure_transcript_counts),
            sort_keys=True,
            separators=(",", ":"),
        ),
        "encountered_rejection_event_counts": json.dumps(
            dict(summary.encountered_rejection_event_counts),
            sort_keys=True,
            separators=(",", ":"),
        ),
        "encountered_rejection_side_counts": json.dumps(
            dict(summary.encountered_rejection_side_counts),
            sort_keys=True,
            separators=(",", ":"),
        ),
        "encountered_rejection_transcript_counts": json.dumps(
            dict(summary.encountered_rejection_transcript_counts),
            sort_keys=True,
            separators=(",", ":"),
        ),
    }
    write_dict_tsv([row], path, list(row))
