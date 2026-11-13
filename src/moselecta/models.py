"""Core data containers used by the MoSelecta pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .constants import TRANSCRIPT_FEATURES

@dataclass(slots=True)
class FeatureRecord:
    seqid: str
    source: str
    feature_type: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    score: str
    strand: str
    phase: str
    attrs: dict[str, str]
    raw_attrs: str
    input_file: Path
    source_name: str
    line_number: int
    record_index: int
    raw_line: str

    @property
    def id(self) -> Optional[str]:
        return self.attrs.get("ID")

    @property
    def parents(self) -> list[str]:
        parent = self.attrs.get("Parent")
        if parent:
            return [p for p in parent.split(",") if p]
        # GTF fallback. These are not true Parent links, but allow grouping.
        transcript_id = self.attrs.get("transcript_id")
        gene_id = self.attrs.get("gene_id")
        if self.feature_type in TRANSCRIPT_FEATURES and gene_id:
            return [gene_id]
        if transcript_id:
            return [transcript_id]
        return []

    @property
    def start0(self) -> int:
        return self.start - 1

    @property
    def end0(self) -> int:
        return self.end


@dataclass(slots=True)
class SequenceIntegrityMetrics:
    """Compact genome-derived metrics retained on a candidate model.

    Identity/coordinate columns are reconstructed from the owning CodingModel
    when a dedicated TSV is written, avoiding a second full list of large row
    dictionaries for large genomes.
    """

    status: str
    start_codon: str
    start_codon_status: str
    stop_codon: str
    stop_codon_status: str
    internal_stop_count: int
    phase_status: str
    phase_mismatch_count: int
    noncanonical_splice_intron_count: int
    noncanonical_splice_intron_fraction: float
    noncanonical_splice_details: str
    has_ambiguous_bases: bool
    penalty_score: float
    penalty_components: str
    flags: str


@dataclass(slots=True)
class TerminalStopCanonicalisation:
    """Audit/provenance for early local terminal-stop normalisation."""

    status: str
    offset_bp: int = 0
    original_boundary: Optional[int] = None
    canonical_boundary: Optional[int] = None
    stop_codon: str = "."
    frame_confidence_basis: str = "."
    first_cds_phase: Optional[int] = None
    skip_reason: str = "."
    exon_adjusted: bool = False
    newly_introduced_same_strand_cds_overlap_count: int = 0
    newly_introduced_same_strand_cds_overlap_models: str = "."


@dataclass(frozen=True, slots=True)
class FivePrimeBoundaryCanonicalisation:
    """Audit/provenance for the bounded class-B 5' CDS correction."""

    status: str
    skip_reason: str = "."
    bases_trimmed: int = 0
    original_boundary: Optional[int] = None
    canonical_boundary: Optional[int] = None
    exon_context: str = "."
    five_prime_utr_extended: bool = False
    original_source_first_phase: str = "."
    canonical_first_phase: Optional[int] = None
    start_codon: str = "."
    stop_codon: str = "."
    stop_representation: str = "."
    original_start_status: str = "."
    canonical_start_status: str = "."
    original_stop_status: str = "."
    canonical_stop_status: str = "."
    original_internal_stop_count: int = 0
    canonical_internal_stop_count: int = 0


@dataclass(frozen=True, slots=True)
class CDSPhaseOutput:
    """Output-only decision for a model whose first CDS phase was inspected.

    This state is deliberately separate from :class:`SequenceIntegrityMetrics`:
    it controls only GFF3 phase serialization and is never a scoring input.
    """

    status: str
    first_phase: Optional[int] = None
    reason: str = "."
    source_first_phase: Optional[int] = None


@dataclass(frozen=True, slots=True)
class CDSPhaseRecovery:
    """Audit/provenance for pre-integrity biological phase recovery."""

    status: str
    original_source_first_phase: Optional[int] = None
    recovered_first_phase: Optional[int] = None
    evidence_basis: str = "."
    stop_representation: str = "."
    stop_codon: str = "."
    original_effective_phase: int = 0
    original_stop_status: str = "."
    recovered_stop_status: str = "."
    original_internal_stop_count: int = 0
    recovered_internal_stop_count: int = 0
    reason: str = "."


@dataclass(slots=True)
class CodingModel:
    source_name: str
    input_file: Path
    input_gene_id: str
    input_transcript_id: str
    output_gene_id: str
    output_transcript_id: str
    chrom: str
    strand: str
    cds_parts: list[tuple[int, int]]  # sorted, 0-based half-open
    gene_record: Optional[FeatureRecord]
    transcript_record: Optional[FeatureRecord]
    child_records: list[FeatureRecord]
    used_record_indices: set[int] = field(default_factory=set)
    superlocus_id: Optional[str] = None
    sublocus_id: Optional[str] = None
    sequence_integrity: Optional[SequenceIntegrityMetrics] = None
    cds_phase_output: Optional[CDSPhaseOutput] = None
    cds_phase_recovery: Optional[CDSPhaseRecovery] = None
    cds_terminus_refined: Optional[str] = None
    five_prime_boundary_canonicalisation: Optional[FivePrimeBoundaryCanonicalisation] = None
    terminal_stop_canonicalisation: Optional[TerminalStopCanonicalisation] = None
    _cds_length: int = field(init=False, repr=False)
    _cds_introns: tuple[tuple[int, int], ...] = field(init=False, repr=False)
    _cds_splice_sites: frozenset[int] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        # Keep only the geometry that is both compact and repeatedly reused.
        # This avoids patch-16's per-property __dict__ cache growth while not
        # recomputing intron/splice geometry in every evidence comparison.
        self._cds_length = sum(end0 - start0 for start0, end0 in self.cds_parts)
        self._cds_introns = tuple(
            (left[1], right[0])
            for left, right in zip(self.cds_parts, self.cds_parts[1:])
        )
        self._cds_splice_sites = frozenset(
            coordinate
            for intron in self._cds_introns
            for coordinate in intron
        )

    @property
    def sorted_cds_parts(self) -> list[tuple[int, int]]:
        # CDS parts are normalised once during parsing. Returning the original
        # tuple avoids duplicating every model's exon chain in a cached value.
        return self.cds_parts

    @property
    def cds_start0(self) -> int:
        return self.sorted_cds_parts[0][0]

    @property
    def cds_end0(self) -> int:
        return self.sorted_cds_parts[-1][1]

    @property
    def cds_length(self) -> int:
        return self._cds_length

    @property
    def cds_span_bp(self) -> int:
        return self.cds_end0 - self.cds_start0

    @property
    def max_cds_intron_bp(self) -> int:
        introns = self.cds_introns
        if not introns:
            return 0
        return max(e - s for s, e in introns)

    @property
    def cds_exon_count(self) -> int:
        return len(self.cds_parts)

    @property
    def cds_introns(self) -> tuple[tuple[int, int], ...]:
        return self._cds_introns

    @property
    def cds_intron_set(self) -> frozenset[tuple[int, int]]:
        return frozenset(self.cds_introns)

    @property
    def cds_splice_sites(self) -> frozenset[int]:
        return self._cds_splice_sites

    @property
    def cds_exon_lengths(self) -> list[int]:
        return [end - start for start, end in self.sorted_cds_parts]

    @property
    def terminal_cds_exons(self) -> tuple[tuple[int, int], ...]:
        if not self.sorted_cds_parts:
            return tuple()
        if len(self.sorted_cds_parts) == 1:
            return (self.sorted_cds_parts[0],)
        return (self.sorted_cds_parts[0], self.sorted_cds_parts[-1])

    @property
    def internal_cds_exons(self) -> tuple[tuple[int, int], ...]:
        if len(self.sorted_cds_parts) <= 2:
            return tuple()
        return self.sorted_cds_parts[1:-1]

    @property
    def terminal_cds_exon_lengths(self) -> list[int]:
        """Lengths of CDS exons at the two ends of the CDS chain.

        For a single-exon CDS, the single exon is terminal. Exons are considered
        in genomic order; for minus-strand models the biological 5'/3' order is
        reversed, but the two terminal CDS exons are the same genomic end exons.
        """
        parts = self.sorted_cds_parts
        if not parts:
            return []
        if len(parts) == 1:
            start, end = parts[0]
            return [end - start]
        left_start, left_end = parts[0]
        right_start, right_end = parts[-1]
        return [left_end - left_start, right_end - right_start]

    @property
    def internal_cds_exon_lengths(self) -> list[int]:
        """Lengths of non-terminal CDS exons. One- and two-exon CDS models have none."""
        parts = self.sorted_cds_parts
        if len(parts) <= 2:
            return []
        return [end - start for start, end in parts[1:-1]]

    @property
    def cds_intron_lengths(self) -> list[int]:
        return [end - start for start, end in self.cds_introns]

    @property
    def is_single_exon(self) -> bool:
        return self.cds_exon_count == 1

    @property
    def is_multi_exon(self) -> bool:
        return self.cds_exon_count > 1

    @property
    def model_start0(self) -> int:
        starts = [r.start0 for r in self.child_records]
        if self.transcript_record:
            starts.append(self.transcript_record.start0)
        if not starts:
            starts.append(self.cds_start0)
        return min(starts)

    @property
    def model_end0(self) -> int:
        ends = [r.end0 for r in self.child_records]
        if self.transcript_record:
            ends.append(self.transcript_record.end0)
        if not ends:
            ends.append(self.cds_end0)
        return max(ends)


@dataclass(slots=True)
class ExcludedModel:
    source_name: str
    input_file: Path
    input_gene_id: str
    input_transcript_id: str
    reason: str
    cds_count: int
    cds_chroms: list[str]
    cds_strands: list[str]
    cds_span_bp: Optional[int] = None
    max_cds_intron_bp: Optional[int] = None


@dataclass(slots=True)
class ValidationFailure:
    scope: str  # model or superlocus
    failure_type: str
    action: str
    value: int
    threshold: int
    superlocus_id: str = "."
    chrom: str = "."
    start0: Optional[int] = None
    end0: Optional[int] = None
    strand: str = "."
    source_name: str = "."
    input_file: str = "."
    input_gene_id: str = "."
    input_transcript_id: str = "."
    details: str = "."


@dataclass(slots=True)
class Superlocus:
    superlocus_id: str
    chrom: str
    strand: str
    start0: int
    end0: int
    models: list[CodingModel] = field(default_factory=list)


@dataclass(slots=True)
class Sublocus:
    sublocus_id: str
    superlocus_id: str
    chrom: str
    strand: str
    start0: int
    end0: int
    models: list[CodingModel] = field(default_factory=list)
    sublocus_type: str = "."  # multi_exon, single_exon, mixed, singleton
    flags: str = "."


@dataclass(slots=True)
class RelationshipEdge:
    superlocus_id: str
    sublocus_id: str
    model_a: CodingModel
    model_b: CodingModel
    edge_type: str
    shared_splice_site_count: int = 0
    shared_cds_intron_count: int = 0
    reciprocal_overlap_a: Optional[float] = None
    reciprocal_overlap_b: Optional[float] = None
    overlap_bp: Optional[int] = None


@dataclass(slots=True)
class ModelIntrinsicMetrics:
    model: CodingModel
    cds_length: int
    cds_length_mod3: int
    cds_exon_count: int
    cds_intron_count: int
    cds_span_bp: int
    min_cds_exon_bp: int
    max_cds_exon_bp: int
    median_cds_exon_bp: float
    min_terminal_cds_exon_bp: int
    min_internal_cds_exon_bp: int
    median_terminal_cds_exon_bp: float
    median_internal_cds_exon_bp: float
    short_terminal_cds_exon_count: int
    short_internal_cds_exon_count: int
    min_cds_intron_bp: int
    max_cds_intron_bp: int
    median_cds_intron_bp: float
    short_cds_intron_count: int
    long_cds_intron_count: int
    flags: str


@dataclass(slots=True)
class SourceAutoWeight:
    source_name: str
    coding_models: int
    cds_exons: int
    supported_cds_exons: int
    unsupported_cds_exons: int
    supported_cds_exon_fraction: float
    unsupported_cds_exon_fraction: float
    internal_cds_exons: int
    supported_internal_cds_exons: int
    unsupported_internal_cds_exons: int
    supported_internal_cds_exon_fraction: float
    unsupported_internal_cds_exon_fraction: float
    terminal_cds_exons: int
    supported_terminal_cds_exons: int
    unsupported_terminal_cds_exons: int
    supported_terminal_cds_exon_fraction: float
    unsupported_terminal_cds_exon_fraction: float
    cds_introns: int
    supported_cds_introns: int
    unsupported_cds_introns: int
    supported_cds_intron_fraction: float
    unsupported_cds_intron_fraction: float
    cds_splice_sites: int
    supported_cds_splice_sites: int
    unsupported_cds_splice_sites: int
    supported_cds_splice_site_fraction: float
    unsupported_cds_splice_site_fraction: float
    consensus_cds_exon_capture: float
    consensus_internal_cds_exon_capture: float
    consensus_terminal_cds_exon_capture: float
    consensus_cds_intron_capture: float
    consensus_cds_splice_site_capture: float
    auto_agreement_score: float
    auto_source_weight: float
    auto_source_weight_status: str
    auto_source_weight_reason: str
    manual_source_weight: float
    source_prior_mode: str
    effective_source_weight: float
    source_prior_score: float


@dataclass(slots=True)
class ParsedInput:
    input_file: Path
    source_name: str
    records: list[FeatureRecord]
    coding_models: list[CodingModel]
    excluded_models: list[ExcludedModel]
    used_record_indices: set[int]


class UnionFind:
    def __init__(self, items: list[str]):
        self.parent = {item: item for item in items}
        self.rank = {item: 0 for item in items}

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, a: str, b: str) -> None:
        root_a = self.find(a)
        root_b = self.find(b)
        if root_a == root_b:
            return
        if self.rank[root_a] < self.rank[root_b]:
            root_a, root_b = root_b, root_a
        self.parent[root_b] = root_a
        if self.rank[root_a] == self.rank[root_b]:
            self.rank[root_a] += 1
