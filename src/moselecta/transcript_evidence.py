"""Transcript evidence parsing, collapsing, scoring, and model support.

Transcript GFF/GTF inputs are parsed from exon features, collapsed into exact
evidence records, grouped into strand-aware evidence loci, and compared with
candidate CDS models. Transcript evidence scores can contribute to provisional
selection ranking when given a positive selection weight. Scoring-eligible
transcript evidence can also rescue singleton models and participate in
evidence-aware low-score pruning under the configured pruning settings.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .gffio import read_gff_records, sanitize_id_part
from .models import CodingModel, UnionFind
from .loci import reciprocal_overlap
from .intervals import (
    IntervalOverlapIndex,
    interval_list_overlap_bp,
    interval_list_overlap_bp_sorted,
    interval_total_bp,
)
from .tsvio import format_tsv_value, write_dict_tsv
from .protein_evidence import fmt_counts, fmt_flags, fmt_int_summary
from .scoring import transcript_cdna_context_concern_fields, weighted_metric_value


TRANSCRIPT_EVIDENCE_LARGE_SPAN_BP = 200_000
TRANSCRIPT_EVIDENCE_MANY_COLLAPSED_TRANSCRIPTS = 100
TRANSCRIPT_EVIDENCE_MANY_SOURCES = 5
TRANSCRIPT_EVIDENCE_HIGH_MEMBER_COUNT = 500


@dataclass(slots=True)
class TranscriptEvidenceRawTranscript:
    raw_id: str
    source_label: str
    source_index: int
    input_file: Path
    original_transcript_id: str
    chrom: str
    strand: str
    exon_parts: tuple[tuple[int, int], ...]
    line_count: int
    evidence_superlocus_id: Optional[str] = None
    collapsed_id: Optional[str] = None

    @property
    def sorted_exon_parts(self) -> tuple[tuple[int, int], ...]:
        return self.exon_parts

    @property
    def start0(self) -> int:
        return self.sorted_exon_parts[0][0]

    @property
    def end0(self) -> int:
        return self.sorted_exon_parts[-1][1]

    @property
    def span_bp(self) -> int:
        return self.end0 - self.start0

    @property
    def exonic_length(self) -> int:
        return sum(e - s for s, e in self.exon_parts)

    @property
    def exon_count(self) -> int:
        return len(self.exon_parts)

    @property
    def introns(self) -> tuple[tuple[int, int], ...]:
        parts = self.sorted_exon_parts
        return tuple((a[1], b[0]) for a, b in zip(parts, parts[1:]))

    @property
    def intron_set(self) -> frozenset[tuple[int, int]]:
        return frozenset(self.introns)

    @property
    def splice_sites(self) -> frozenset[int]:
        sites: set[int] = set()
        for s, e in self.introns:
            sites.add(s); sites.add(e)
        return frozenset(sites)

    @property
    def is_single_exon(self) -> bool:
        return len(self.exon_parts) == 1

    @property
    def is_multi_exon(self) -> bool:
        return len(self.exon_parts) > 1


@dataclass(slots=True)
class TranscriptEvidenceCollapsedTranscript:
    transcript_evidence_id: str
    chrom: str
    strand: str
    exon_parts: tuple[tuple[int, int], ...]
    evidence_superlocus_id: str
    collapse_key_type: str
    source_member_counts: dict[str, int] = field(default_factory=dict)
    member_count: int = 0
    representative_raw_id: str = "."
    evidence_sublocus_id: Optional[str] = None
    flags: str = "."
    span_start0: int = 0
    span_end0: int = 0
    left_terminal_outer_min0: int = 0
    left_terminal_outer_max0: int = 0
    right_terminal_outer_min0: int = 0
    right_terminal_outer_max0: int = 0
    _effective_exon_parts_cache: Optional[tuple[tuple[int, int], ...]] = field(
        default=None,
        init=False,
        repr=False,
        compare=False,
    )
    _exonic_length: int = field(init=False, repr=False)
    _introns: tuple[tuple[int, int], ...] = field(init=False, repr=False)
    _intron_set: frozenset[tuple[int, int]] = field(init=False, repr=False)
    _splice_sites: frozenset[int] = field(init=False, repr=False)
    _sources: frozenset[str] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self._exonic_length = sum(end0 - start0 for start0, end0 in self.exon_parts)
        self._introns = tuple(
            (left[1], right[0])
            for left, right in zip(self.exon_parts, self.exon_parts[1:])
        )
        self._intron_set = frozenset(self._introns)
        self._splice_sites = frozenset(
            coordinate
            for intron in self._introns
            for coordinate in intron
        )
        self._sources = frozenset(
            source for source, count in self.source_member_counts.items() if count > 0
        )

    @property
    def sorted_exon_parts(self) -> tuple[tuple[int, int], ...]:
        return self.exon_parts

    @property
    def start0(self) -> int:
        return self.span_start0 if self.span_end0 > self.span_start0 else self.sorted_exon_parts[0][0]

    @property
    def end0(self) -> int:
        return self.span_end0 if self.span_end0 > self.span_start0 else self.sorted_exon_parts[-1][1]

    @property
    def span_bp(self) -> int:
        return self.end0 - self.start0

    @property
    def exonic_length(self) -> int:
        return self._exonic_length

    @property
    def exon_count(self) -> int:
        return len(self.exon_parts)

    @property
    def introns(self) -> tuple[tuple[int, int], ...]:
        return self._introns

    @property
    def intron_set(self) -> frozenset[tuple[int, int]]:
        return self._intron_set

    @property
    def splice_sites(self) -> frozenset[int]:
        return self._splice_sites

    @property
    def intron_count(self) -> int:
        return len(self.introns)

    @property
    def is_single_exon(self) -> bool:
        return len(self.exon_parts) == 1

    @property
    def is_multi_exon(self) -> bool:
        return len(self.exon_parts) > 1

    @property
    def source_count(self) -> int:
        return sum(1 for count in self.source_member_counts.values() if count > 0)

    @property
    def sources(self) -> frozenset[str]:
        return self._sources

    @property
    def left_terminal_outer_range_bp(self) -> int:
        return self.left_terminal_outer_max0 - self.left_terminal_outer_min0

    @property
    def right_terminal_outer_range_bp(self) -> int:
        return self.right_terminal_outer_max0 - self.right_terminal_outer_min0

    @property
    def five_prime_terminal_outer_min0(self) -> Optional[int]:
        if self.strand == "+":
            return self.left_terminal_outer_min0
        if self.strand == "-":
            return self.right_terminal_outer_min0
        return None

    @property
    def five_prime_terminal_outer_max0(self) -> Optional[int]:
        if self.strand == "+":
            return self.left_terminal_outer_max0
        if self.strand == "-":
            return self.right_terminal_outer_max0
        return None

    @property
    def three_prime_terminal_outer_min0(self) -> Optional[int]:
        if self.strand == "+":
            return self.right_terminal_outer_min0
        if self.strand == "-":
            return self.left_terminal_outer_min0
        return None

    @property
    def three_prime_terminal_outer_max0(self) -> Optional[int]:
        if self.strand == "+":
            return self.right_terminal_outer_max0
        if self.strand == "-":
            return self.left_terminal_outer_max0
        return None

    @property
    def five_prime_terminal_outer_range_bp(self) -> object:
        if self.five_prime_terminal_outer_min0 is None or self.five_prime_terminal_outer_max0 is None:
            return "."
        return self.five_prime_terminal_outer_max0 - self.five_prime_terminal_outer_min0

    @property
    def three_prime_terminal_outer_range_bp(self) -> object:
        if self.three_prime_terminal_outer_min0 is None or self.three_prime_terminal_outer_max0 is None:
            return "."
        return self.three_prime_terminal_outer_max0 - self.three_prime_terminal_outer_min0


@dataclass(slots=True)
class TranscriptEvidenceSuperlocus:
    evidence_superlocus_id: str
    chrom: str
    strand: str
    start0: int
    end0: int
    raw_transcripts: list[TranscriptEvidenceRawTranscript] = field(default_factory=list)
    collapsed_transcripts: list[TranscriptEvidenceCollapsedTranscript] = field(default_factory=list)
    raw_transcript_count: int = 0
    raw_source_member_counts: dict[str, int] = field(default_factory=dict)


@dataclass(slots=True)
class TranscriptEvidenceSublocus:
    evidence_sublocus_id: str
    evidence_superlocus_id: str
    chrom: str
    strand: str
    start0: int
    end0: int
    collapsed_transcripts: list[TranscriptEvidenceCollapsedTranscript] = field(default_factory=list)
    sublocus_type: str = "."
    flags: str = "."


def split_transcript_label_list(values: Optional[list[str]]) -> list[str]:
    labels: list[str] = []
    if not values:
        return labels
    for item in values:
        for part in str(item).split(','):
            part = part.strip()
            if part:
                labels.append(sanitize_id_part(part))
    return labels


def default_transcript_source_labels(n: int) -> list[str]:
    return [f"transcript_{i:03d}" for i in range(1, n + 1)]


def transcript_source_labels(paths: list[Path], supplied: Optional[list[str]]) -> list[str]:
    labels = split_transcript_label_list(supplied)
    if labels:
        if len(labels) != len(paths):
            raise SystemExit("--evidence-transcript-source-name must provide one label per --evidence-transcripts file")
    else:
        labels = default_transcript_source_labels(len(paths))
    seen: Counter[str] = Counter()
    out: list[str] = []
    for label in labels:
        base = sanitize_id_part(label)
        seen[base] += 1
        out.append(base if seen[base] == 1 else f"{base}_{seen[base]}")
    return out


def _exon_parent_ids(rec) -> list[str]:
    parents = rec.parents
    if parents:
        return parents
    for key in ("transcript_id", "ID"):
        value = rec.attrs.get(key)
        if value:
            return [value]
    return []


def parse_transcript_evidence_file(path: Path, source_label: str, source_index: int, *, strict: bool = False) -> list[TranscriptEvidenceRawTranscript]:
    """Parse transcript evidence GFF/GTF exon features into raw transcripts.

    Only exon features are used. CDS, UTR, gene, and transcript container
    features are ignored for transcript evidence construction.
    """
    records = read_gff_records(path, source_label, strict=strict)
    exons_by_parent: dict[str, list] = defaultdict(list)
    for rec in records:
        if rec.feature_type.lower() != "exon":
            continue
        parents = _exon_parent_ids(rec)
        if not parents:
            if strict:
                raise ValueError(
                    f"{path}: exon on line {rec.line_number} has no Parent, transcript_id, or ID"
                )
            # Permissive parsing skips structurally unassignable exons rather
            # than inventing a synthetic one-exon transcript.
            continue
        for parent in parents:
            exons_by_parent[parent].append(rec)

    transcripts: list[TranscriptEvidenceRawTranscript] = []
    for ordinal, (original_id, exon_records) in enumerate(sorted(exons_by_parent.items()), start=1):
        chroms = sorted({r.seqid for r in exon_records})
        strands = sorted({r.strand for r in exon_records})
        valid_strands = {s for s in strands if s in {"+", "-"}}
        if len(chroms) != 1:
            if strict:
                raise ValueError(f"{path}: transcript evidence record {original_id} has mixed chromosomes")
            continue
        if len(valid_strands) == 1 and len(strands) == 1:
            transcript_strand = next(iter(valid_strands))
        elif len(valid_strands) == 0 and set(strands) <= {"."}:
            transcript_strand = "."
        else:
            if strict:
                raise ValueError(f"{path}: transcript evidence record {original_id} has mixed strands")
            continue
        exon_parts = tuple(sorted(set((r.start0, r.end0) for r in exon_records if r.end0 > r.start0)))
        if not exon_parts:
            continue
        raw_id = f"{source_label}|{sanitize_id_part(original_id)}|raw{ordinal:07d}"
        transcripts.append(
            TranscriptEvidenceRawTranscript(
                raw_id=raw_id,
                source_label=source_label,
                source_index=source_index,
                input_file=path,
                original_transcript_id=original_id,
                chrom=chroms[0],
                strand=transcript_strand,
                exon_parts=exon_parts,
                line_count=len(exon_records),
            )
        )
    return transcripts


def parse_transcript_evidence(paths: list[Path], labels: list[str], *, strict: bool = False) -> list[TranscriptEvidenceRawTranscript]:
    raw: list[TranscriptEvidenceRawTranscript] = []
    for idx, (path, label) in enumerate(zip(paths, labels), start=1):
        raw.extend(parse_transcript_evidence_file(path, label, idx, strict=strict))
    return raw


def make_transcript_superlocus_id(chrom: str, start0: int, end0: int, strand: str, index: int) -> str:
    if strand == "+":
        strand_label = "plus"
    elif strand == "-":
        strand_label = "minus"
    else:
        strand_label = "unstranded"
    return f"transcript_superlocus_{index:07d}_{sanitize_id_part(chrom)}_{start0 + 1}_{end0}_{strand_label}"


def build_transcript_superloci(raw_transcripts: list[TranscriptEvidenceRawTranscript]) -> list[TranscriptEvidenceSuperlocus]:
    """Build same-strand exon-overlap connected components from raw transcripts."""
    superloci: list[TranscriptEvidenceSuperlocus] = []
    by_key: dict[tuple[str, str], list[TranscriptEvidenceRawTranscript]] = defaultdict(list)
    for tx in raw_transcripts:
        by_key[(tx.chrom, tx.strand)].append(tx)
    next_index = 1
    for (chrom, strand), items in sorted(by_key.items()):
        ids = [tx.raw_id for tx in items]
        uf = UnionFind(ids)
        exons: list[tuple[int, int, str]] = []
        for tx in items:
            for start0, end0 in tx.exon_parts:
                exons.append((start0, end0, tx.raw_id))
        exons.sort(key=lambda x: (x[0], x[1], x[2]))
        active: list[tuple[int, str]] = []
        for start0, end0, raw_id in exons:
            active = [(ae, aid) for ae, aid in active if ae > start0]
            for _ae, aid in active:
                if aid != raw_id:
                    uf.union(raw_id, aid)
            active.append((end0, raw_id))

        comps: dict[str, list[TranscriptEvidenceRawTranscript]] = defaultdict(list)
        for tx in items:
            comps[uf.find(tx.raw_id)].append(tx)
        ordered = sorted(
            comps.values(),
            key=lambda comp: (min(t.start0 for t in comp), max(t.end0 for t in comp), min(t.raw_id for t in comp)),
        )
        for comp in ordered:
            start0 = min(t.start0 for t in comp)
            end0 = max(t.end0 for t in comp)
            sl_id = make_transcript_superlocus_id(chrom, start0, end0, strand, next_index)
            next_index += 1
            sl = TranscriptEvidenceSuperlocus(
                evidence_superlocus_id=sl_id,
                chrom=chrom,
                strand=strand,
                start0=start0,
                end0=end0,
                raw_transcripts=sorted(comp, key=lambda t: (t.start0, t.end0, t.source_label, t.raw_id)),
            )
            for member in sl.raw_transcripts:
                member.evidence_superlocus_id = sl_id
            superloci.append(sl)
    return sorted(superloci, key=lambda sl: (sl.chrom, sl.start0, sl.end0, sl.strand, sl.evidence_superlocus_id))


def transcript_collapse_key(tx: TranscriptEvidenceRawTranscript) -> tuple[str, tuple[tuple[int, int], ...]]:
    if tx.is_multi_exon:
        return ("exact_intron_chain", tx.introns)
    return ("exact_single_exon_interval", tx.exon_parts)


def _representative_member(members: list[TranscriptEvidenceRawTranscript]) -> TranscriptEvidenceRawTranscript:
    """Return the raw member whose structure represents a collapsed group.

    The representative structure is the most frequent exon chain, with ties
    resolved by shortest genomic span, coordinate order, and raw ID. Returning
    the member itself keeps ``representative_raw_id`` consistent with the
    reported representative exon chain.
    """
    counts: Counter[tuple[tuple[int, int], ...]] = Counter(m.exon_parts for m in members)
    return min(
        members,
        key=lambda m: (
            -counts[m.exon_parts],
            m.exon_parts[-1][1] - m.exon_parts[0][0],
            m.exon_parts,
            m.raw_id,
        ),
    )


def collapse_transcript_superlocus(sl: TranscriptEvidenceSuperlocus) -> list[TranscriptEvidenceCollapsedTranscript]:
    if sl.raw_transcript_count <= 0:
        sl.raw_transcript_count = len(sl.raw_transcripts)
    if not sl.raw_source_member_counts:
        sl.raw_source_member_counts = dict(sorted(Counter(raw.source_label for raw in sl.raw_transcripts).items()))
    groups: dict[tuple[str, tuple[tuple[int, int], ...]], list[TranscriptEvidenceRawTranscript]] = defaultdict(list)
    for raw in sl.raw_transcripts:
        groups[transcript_collapse_key(raw)].append(raw)
    collapsed: list[TranscriptEvidenceCollapsedTranscript] = []
    for idx, (key, members) in enumerate(
        sorted(groups.items(), key=lambda item: (item[1][0].start0, item[1][0].end0, item[0][0], item[0][1])),
        start=1,
    ):
        key_type, _key_parts = key
        representative = _representative_member(members)
        exon_parts = representative.exon_parts
        cid = f"{sl.evidence_superlocus_id}.collapsed{idx:05d}"
        counts = Counter(m.source_label for m in members)
        left_values = [m.exon_parts[0][0] for m in members]
        right_values = [m.exon_parts[-1][1] for m in members]
        obj = TranscriptEvidenceCollapsedTranscript(
            transcript_evidence_id=cid,
            chrom=sl.chrom,
            strand=sl.strand,
            exon_parts=exon_parts,
            evidence_superlocus_id=sl.evidence_superlocus_id,
            collapse_key_type=key_type,
            source_member_counts=dict(sorted(counts.items())),
            member_count=len(members),
            representative_raw_id=representative.raw_id,
            span_start0=min(m.start0 for m in members),
            span_end0=max(m.end0 for m in members),
            left_terminal_outer_min0=min(left_values),
            left_terminal_outer_max0=max(left_values),
            right_terminal_outer_min0=min(right_values),
            right_terminal_outer_max0=max(right_values),
        )
        obj.flags = collapsed_transcript_flags(obj)
        for raw in members:
            raw.collapsed_id = cid
        collapsed.append(obj)
    sl.collapsed_transcripts = collapsed
    return collapsed


def collapse_transcript_evidence(superloci: list[TranscriptEvidenceSuperlocus]) -> list[TranscriptEvidenceCollapsedTranscript]:
    out: list[TranscriptEvidenceCollapsedTranscript] = []
    for sl in superloci:
        out.extend(collapse_transcript_superlocus(sl))
    return out


def release_transcript_raw_transcripts(
    raw_transcripts: list[TranscriptEvidenceRawTranscript],
    superloci: list[TranscriptEvidenceSuperlocus],
) -> int:
    """Drop raw transcript records after collapse when no later raw output needs them."""
    released = len(raw_transcripts)
    for sl in superloci:
        if sl.raw_transcript_count <= 0:
            sl.raw_transcript_count = len(sl.raw_transcripts)
        if not sl.raw_source_member_counts:
            sl.raw_source_member_counts = dict(sorted(Counter(raw.source_label for raw in sl.raw_transcripts).items()))
        sl.raw_transcripts.clear()
    raw_transcripts.clear()
    return released


def collapsed_transcript_flags(tx: TranscriptEvidenceCollapsedTranscript) -> str:
    flags: list[str] = []
    flags.append("single_exon" if tx.is_single_exon else "multi_exon")
    flags.append(tx.collapse_key_type)
    if tx.source_count <= 1:
        flags.append("single_source_only")
    else:
        flags.append("multi_source_exact_chain")
    if isinstance(tx.five_prime_terminal_outer_range_bp, int) and tx.five_prime_terminal_outer_range_bp > 0:
        flags.append("variable_5prime_terminal_outer")
    if isinstance(tx.three_prime_terminal_outer_range_bp, int) and tx.three_prime_terminal_outer_range_bp > 0:
        flags.append("variable_3prime_terminal_outer")
    if tx.strand not in {"+", "-"}:
        flags.append("unstranded")
    return fmt_flags(flags)


def transcript_sublocus_type(items: list[TranscriptEvidenceCollapsedTranscript]) -> str:
    if all(item.is_single_exon for item in items):
        return "single_exon"
    if all(item.is_multi_exon for item in items):
        return "multi_exon"
    return "mixed"


def transcript_sublocus_flags(items: list[TranscriptEvidenceCollapsedTranscript], sublocus_type: str, start0: int, end0: int) -> str:
    source_count = len({source for tx in items for source in tx.sources})
    member_count = sum(tx.member_count for tx in items)
    collapsed_count = len(items)
    span_bp = end0 - start0
    flags: list[str] = []
    if source_count <= 1:
        flags.append("single_source_only")
    if sublocus_type == "mixed":
        flags.append("mixed_single_multi")
    if span_bp > TRANSCRIPT_EVIDENCE_LARGE_SPAN_BP:
        flags.append("large_span")
    if collapsed_count >= TRANSCRIPT_EVIDENCE_MANY_COLLAPSED_TRANSCRIPTS:
        flags.append("many_collapsed_transcripts")
    if source_count >= TRANSCRIPT_EVIDENCE_MANY_SOURCES:
        flags.append("many_sources")
    if member_count >= TRANSCRIPT_EVIDENCE_HIGH_MEMBER_COUNT:
        flags.append("high_member_count")
    return fmt_flags(flags)


def make_transcript_sublocus_id(superlocus_id: str, index: int) -> str:
    return f"{superlocus_id}.transcript_sublocus{index:03d}"


def build_transcript_subloci_for_superlocus(sl: TranscriptEvidenceSuperlocus, *, single_exon_min_reciprocal_overlap: float) -> list[TranscriptEvidenceSublocus]:
    items = sorted(sl.collapsed_transcripts, key=lambda a: (a.start0, a.end0, a.transcript_evidence_id))
    if not items:
        return []
    uf = UnionFind([a.transcript_evidence_id for a in items])

    splice_index: dict[int, list[TranscriptEvidenceCollapsedTranscript]] = defaultdict(list)
    for tx in items:
        if tx.is_multi_exon:
            for site in tx.splice_sites:
                splice_index[site].append(tx)
    for site_items in splice_index.values():
        site_items = sorted(site_items, key=lambda a: a.transcript_evidence_id)
        for i, a in enumerate(site_items):
            for b in site_items[i + 1:]:
                uf.union(a.transcript_evidence_id, b.transcript_evidence_id)

    singles = [a for a in items if a.is_single_exon]
    singles = sorted(singles, key=lambda a: (a.start0, a.end0, a.transcript_evidence_id))
    for i, a in enumerate(singles):
        j = i + 1
        while j < len(singles) and singles[j].start0 < a.end0:
            b = singles[j]
            _overlap, ro_a, ro_b = reciprocal_overlap(a.start0, a.end0, b.start0, b.end0)
            if min(ro_a, ro_b) >= single_exon_min_reciprocal_overlap:
                uf.union(a.transcript_evidence_id, b.transcript_evidence_id)
            j += 1

    comps: dict[str, list[TranscriptEvidenceCollapsedTranscript]] = defaultdict(list)
    for tx in items:
        comps[uf.find(tx.transcript_evidence_id)].append(tx)
    ordered = sorted(comps.values(), key=lambda comp: (min(a.start0 for a in comp), max(a.end0 for a in comp), min(a.transcript_evidence_id for a in comp)))
    out: list[TranscriptEvidenceSublocus] = []
    for idx, comp in enumerate(ordered, start=1):
        sid = make_transcript_sublocus_id(sl.evidence_superlocus_id, idx)
        start0 = min(a.start0 for a in comp); end0 = max(a.end0 for a in comp)
        stype = transcript_sublocus_type(comp)
        comp_sorted = sorted(comp, key=lambda a: a.transcript_evidence_id)
        sub = TranscriptEvidenceSublocus(
            evidence_sublocus_id=sid,
            evidence_superlocus_id=sl.evidence_superlocus_id,
            chrom=sl.chrom,
            strand=sl.strand,
            start0=start0,
            end0=end0,
            collapsed_transcripts=comp_sorted,
            sublocus_type=stype,
            flags=transcript_sublocus_flags(comp_sorted, stype, start0, end0),
        )
        for tx in sub.collapsed_transcripts:
            tx.evidence_sublocus_id = sid
        out.append(sub)
    return out


def build_transcript_subloci(superloci: list[TranscriptEvidenceSuperlocus], *, single_exon_min_reciprocal_overlap: float) -> list[TranscriptEvidenceSublocus]:
    out: list[TranscriptEvidenceSublocus] = []
    for sl in superloci:
        out.extend(build_transcript_subloci_for_superlocus(sl, single_exon_min_reciprocal_overlap=single_exon_min_reciprocal_overlap))
    return out


def _chain_text(parts: tuple[tuple[int, int], ...]) -> str:
    return ",".join(f"{s + 1}-{e}" for s, e in parts) if parts else "."


def _intron_chain_text(introns: tuple[tuple[int, int], ...]) -> str:
    return ",".join(f"{s + 1}-{e}" for s, e in introns) if introns else "."


def transcript_summary_counts(
    *,
    transcript_paths: list[Path],
    transcript_source_labels: list[str],
    raw_transcripts: list[TranscriptEvidenceRawTranscript],
    superloci: list[TranscriptEvidenceSuperlocus],
    collapsed_transcripts: list[TranscriptEvidenceCollapsedTranscript],
    subloci: list[TranscriptEvidenceSublocus],
    raw_to_collapsed_mapping_requested: bool,
) -> dict[str, object]:
    raw_by_source = Counter(t.source_label for t in raw_transcripts)
    collapsed_by_source: Counter[str] = Counter()
    collapsed_member_by_source: Counter[str] = Counter()
    for tx in collapsed_transcripts:
        for source, count in tx.source_member_counts.items():
            if count > 0:
                collapsed_by_source[source] += 1
                collapsed_member_by_source[source] += count
    sublocus_types = Counter(sub.sublocus_type for sub in subloci)
    collapse_types = Counter(tx.collapse_key_type for tx in collapsed_transcripts)
    flag_counts: Counter[str] = Counter()
    for tx in collapsed_transcripts:
        if tx.flags and tx.flags != ".":
            for flag in tx.flags.split(','):
                if flag:
                    flag_counts[flag] += 1
    source_raw_values = [raw_by_source.get(label, 0) for label in transcript_source_labels]
    source_collapsed_values = [collapsed_by_source.get(label, 0) for label in transcript_source_labels]
    source_member_values = [collapsed_member_by_source.get(label, 0) for label in transcript_source_labels]
    collapse_ratio = (len(collapsed_transcripts) / len(raw_transcripts)) if raw_transcripts else 0.0
    summary = {
        "file_count": len(transcript_paths),
        "source_count": len(transcript_source_labels),
        "source_labels": ",".join(transcript_source_labels) if transcript_source_labels else ".",
        "input_files": ",".join(str(p) for p in transcript_paths) if transcript_paths else ".",
        "raw_transcript_count": len(raw_transcripts),
        "collapsed_transcript_count": len(collapsed_transcripts),
        "collapse_ratio": collapse_ratio,
        "raw_to_collapsed_mapping_requested": raw_to_collapsed_mapping_requested,
        "superlocus_count": len(superloci),
        "sublocus_count": len(subloci),
        "single_exon_raw_transcript_count": sum(1 for t in raw_transcripts if t.is_single_exon),
        "multi_exon_raw_transcript_count": sum(1 for t in raw_transcripts if t.is_multi_exon),
        "unstranded_raw_transcript_count": sum(1 for t in raw_transcripts if t.strand == "."),
        "unstranded_single_exon_raw_transcript_count": sum(1 for t in raw_transcripts if t.strand == "." and t.is_single_exon),
        "unstranded_multi_exon_raw_transcript_count": sum(1 for t in raw_transcripts if t.strand == "." and t.is_multi_exon),
        "single_exon_collapsed_transcript_count": sum(1 for t in collapsed_transcripts if t.is_single_exon),
        "multi_exon_collapsed_transcript_count": sum(1 for t in collapsed_transcripts if t.is_multi_exon),
        "unstranded_collapsed_transcript_count": sum(1 for t in collapsed_transcripts if t.strand == "."),
        "unstranded_single_exon_collapsed_transcript_count": sum(1 for t in collapsed_transcripts if t.strand == "." and t.is_single_exon),
        "unstranded_multi_exon_collapsed_transcript_count": sum(1 for t in collapsed_transcripts if t.strand == "." and t.is_multi_exon),
        "single_exon_sublocus_count": sublocus_types.get("single_exon", 0),
        "multi_exon_sublocus_count": sublocus_types.get("multi_exon", 0),
        "mixed_sublocus_count": sublocus_types.get("mixed", 0),
        "collapse_type_counts": fmt_counts(dict(collapse_types)),
        "collapsed_flag_counts": fmt_counts(dict(flag_counts)),
        "per_source_raw_transcript_summary": fmt_int_summary(source_raw_values),
        "per_source_collapsed_transcript_summary": fmt_int_summary(source_collapsed_values),
        "per_source_member_count_summary": fmt_int_summary(source_member_values),
        "per_source_raw_transcript_counts": fmt_counts(dict(raw_by_source)),
        "per_source_collapsed_transcript_counts": fmt_counts(dict(collapsed_by_source)),
    }
    return summary


def transcript_summary_log_lines(summary: dict[str, object]) -> list[str]:
    return [
        (
            "Transcript evidence inputs: "
            f"files={int(summary['file_count']):,}; sources={int(summary['source_count']):,}; "
            f"raw_transcripts={int(summary['raw_transcript_count']):,}; "
            f"collapsed_transcripts={int(summary['collapsed_transcript_count']):,}; "
            f"collapse_ratio={float(summary['collapse_ratio']):.3f}; "
            f"raw_to_collapsed_mapping_requested={str(summary['raw_to_collapsed_mapping_requested']).lower()}"
        ),
        (
            "Transcript evidence source summaries: "
            f"raw_transcripts({summary['per_source_raw_transcript_summary']}); "
            f"collapsed_transcripts({summary['per_source_collapsed_transcript_summary']})"
        ),
        (
            "Transcript evidence loci: "
            f"superloci={int(summary['superlocus_count']):,}; subloci={int(summary['sublocus_count']):,}; "
            f"single_exon_subloci={int(summary['single_exon_sublocus_count']):,}; "
            f"multi_exon_subloci={int(summary['multi_exon_sublocus_count']):,}; "
            f"mixed_subloci={int(summary['mixed_sublocus_count']):,}"
        ),
        (
            "Transcript evidence collapsed records: "
            f"single_exon={int(summary['single_exon_collapsed_transcript_count']):,}; "
            f"multi_exon={int(summary['multi_exon_collapsed_transcript_count']):,}; "
            f"collapse_types={summary['collapse_type_counts']}; flags={summary['collapsed_flag_counts']}"
        ),
        (
            "Transcript evidence unstranded records: "
            f"raw={int(summary['unstranded_raw_transcript_count']):,}; "
            f"raw_single_exon={int(summary['unstranded_single_exon_raw_transcript_count']):,}; "
            f"raw_multi_exon={int(summary['unstranded_multi_exon_raw_transcript_count']):,}; "
            f"collapsed={int(summary['unstranded_collapsed_transcript_count']):,}; "
            f"collapsed_single_exon={int(summary['unstranded_single_exon_collapsed_transcript_count']):,}; "
            f"collapsed_multi_exon={int(summary['unstranded_multi_exon_collapsed_transcript_count']):,}"
        ),
    ]


def transcript_summary_tsv_rows(
    *,
    transcript_paths: list[Path],
    transcript_source_labels: list[str],
    raw_transcripts: list[TranscriptEvidenceRawTranscript],
    superloci: list[TranscriptEvidenceSuperlocus],
    collapsed_transcripts: list[TranscriptEvidenceCollapsedTranscript],
    subloci: list[TranscriptEvidenceSublocus],
    raw_to_collapsed_mapping_requested: bool,
    model_support_summary: Optional[dict[str, object]] = None,
) -> list[dict[str, object]]:
    summary = transcript_summary_counts(
        transcript_paths=transcript_paths,
        transcript_source_labels=transcript_source_labels,
        raw_transcripts=raw_transcripts,
        superloci=superloci,
        collapsed_transcripts=collapsed_transcripts,
        subloci=subloci,
        raw_to_collapsed_mapping_requested=raw_to_collapsed_mapping_requested,
    )
    rows: list[dict[str, object]] = []

    def add(scope: str, metric: str, value: object, details: str = ".") -> None:
        rows.append({"summary_scope": scope, "source_label": ".", "input_file": ".", "metric": metric, "value": value, "details": details})

    for metric in (
        "file_count", "source_count", "source_labels", "input_files", "raw_transcript_count",
        "collapsed_transcript_count", "collapse_ratio", "raw_to_collapsed_mapping_requested",
        "superlocus_count", "sublocus_count", "single_exon_raw_transcript_count", "multi_exon_raw_transcript_count",
        "unstranded_raw_transcript_count", "unstranded_single_exon_raw_transcript_count", "unstranded_multi_exon_raw_transcript_count",
        "single_exon_collapsed_transcript_count", "multi_exon_collapsed_transcript_count",
        "unstranded_collapsed_transcript_count", "unstranded_single_exon_collapsed_transcript_count", "unstranded_multi_exon_collapsed_transcript_count",
        "single_exon_sublocus_count", "multi_exon_sublocus_count", "mixed_sublocus_count",
        "collapse_type_counts", "collapsed_flag_counts", "per_source_raw_transcript_summary",
        "per_source_collapsed_transcript_summary", "per_source_member_count_summary",
        "per_source_raw_transcript_counts", "per_source_collapsed_transcript_counts",
    ):
        add("overall", metric, summary[metric])

    model_support_summary = model_support_summary or {}
    for metric in (
        "model_support_row_count",
        "models_with_broad_transcript_association",
        "models_with_scoring_eligible_transcript_evidence",
        "models_with_shared_intron_transcript_evidence",
        "models_with_shared_splice_site_transcript_evidence",
        "models_with_broad_only_transcript_association",
        "models_with_transcript_evidence_score_gt_0",
        "models_with_transcript_evidence_score_scaled_ge_0.75",
        "models_with_partial_or_strong_transcript_support",
        "models_with_partial_transcript_support",
        "models_with_strong_transcript_support",
        "models_with_unstranded_broad_transcript_association",
        "models_with_unstranded_scoring_eligible_transcript_evidence",
        "models_with_stranded_scoring_eligible_transcript_evidence",
        "models_with_only_unstranded_scoring_eligible_transcript_evidence",
        "unstranded_broad_association_total",
        "unstranded_scoring_eligible_total",
        "stranded_scoring_eligible_total",
        "model_status_counts",
        "model_flag_counts",
        "association_class_counts",
        "scoring_association_class_counts",
        "models_with_cdna_context_calculated",
        "models_with_cdna_context_low_cds_fraction",
        "models_with_cdna_context_very_low_cds_fraction",
        "models_with_cdna_context_pruning_concern",
        "models_with_cdna_context_severe_concern",
        "models_with_cdna_context_moderate_structural_concern",
        "models_with_cdna_context_not_best_pairwise",
    ):
        if metric in model_support_summary:
            add("model_support", metric, model_support_summary[metric])

    raw_by_source = Counter(t.source_label for t in raw_transcripts)
    collapsed_by_source: Counter[str] = Counter()
    collapsed_members_by_source: Counter[str] = Counter()
    for tx in collapsed_transcripts:
        for source, count in tx.source_member_counts.items():
            if count > 0:
                collapsed_by_source[source] += 1
                collapsed_members_by_source[source] += count
    path_by_label = {label: str(path) for label, path in zip(transcript_source_labels, transcript_paths)}
    for label in transcript_source_labels:
        for metric, value in (
            ("raw_transcript_count", raw_by_source.get(label, 0)),
            ("collapsed_transcript_count", collapsed_by_source.get(label, 0)),
            ("collapsed_member_count", collapsed_members_by_source.get(label, 0)),
        ):
            rows.append({"summary_scope": "source", "source_label": label, "input_file": path_by_label.get(label, "."), "metric": metric, "value": value, "details": "."})
    return rows


def write_transcript_summary_tsv(rows: list[dict[str, object]], path: Path) -> None:
    write_dict_tsv(rows, path, ["summary_scope", "source_label", "input_file", "metric", "value", "details"])


def write_transcript_raw_to_collapsed_tsv(raw_transcripts: list[TranscriptEvidenceRawTranscript], path: Path) -> None:
    columns = ["transcript_evidence_raw_id", "transcript_evidence_id", "transcript_evidence_superlocus_id", "transcript_source_label", "source_index", "input_file", "original_transcript_id", "chrom", "strand", "start", "end", "span_bp", "exonic_length", "exon_count", "intron_count", "exon_chain", "intron_chain"]
    with path.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for t in sorted(raw_transcripts, key=lambda x: (x.chrom, x.start0, x.end0, x.source_label, x.raw_id)):
            introns = t.introns
            row = (
                t.raw_id,
                t.collapsed_id or ".",
                t.evidence_superlocus_id or ".",
                t.source_label,
                t.source_index,
                str(t.input_file),
                t.original_transcript_id,
                t.chrom,
                t.strand,
                t.start0 + 1,
                t.end0,
                t.span_bp,
                t.exonic_length,
                t.exon_count,
                len(introns),
                _chain_text(t.exon_parts),
                _intron_chain_text(introns),
            )
            out.write("\t".join(format_tsv_value(value) for value in row) + "\n")


def _one_based_optional(value: Optional[int]) -> object:
    return "." if value is None else value + 1


def _end_optional(value: Optional[int]) -> object:
    return "." if value is None else value

def write_transcript_collapsed_transcripts_tsv(collapsed: list[TranscriptEvidenceCollapsedTranscript], path: Path) -> None:
    rows = [{
        "transcript_evidence_id": t.transcript_evidence_id,
        "transcript_evidence_superlocus_id": t.evidence_superlocus_id,
        "transcript_evidence_sublocus_id": t.evidence_sublocus_id or ".",
        "chrom": t.chrom,
        "strand": t.strand,
        "start": t.start0 + 1,
        "end": t.end0,
        "span_bp": t.span_bp,
        "exonic_length": t.exonic_length,
        "exon_count": t.exon_count,
        "intron_count": t.intron_count,
        "collapse_key_type": t.collapse_key_type,
        "member_count": t.member_count,
        "source_count": t.source_count,
        "source_member_counts": fmt_counts(t.source_member_counts),
        "representative_raw_id": t.representative_raw_id,
        "left_terminal_outer_min": t.left_terminal_outer_min0 + 1,
        "left_terminal_outer_max": t.left_terminal_outer_max0 + 1,
        "left_terminal_outer_range_bp": t.left_terminal_outer_range_bp,
        "right_terminal_outer_min": t.right_terminal_outer_min0,
        "right_terminal_outer_max": t.right_terminal_outer_max0,
        "right_terminal_outer_range_bp": t.right_terminal_outer_range_bp,
        "five_prime_terminal_outer_min": _one_based_optional(t.five_prime_terminal_outer_min0) if t.strand == "+" else _end_optional(t.five_prime_terminal_outer_min0),
        "five_prime_terminal_outer_max": _one_based_optional(t.five_prime_terminal_outer_max0) if t.strand == "+" else _end_optional(t.five_prime_terminal_outer_max0),
        "five_prime_terminal_outer_range_bp": t.five_prime_terminal_outer_range_bp,
        "three_prime_terminal_outer_min": _end_optional(t.three_prime_terminal_outer_min0) if t.strand == "+" else _one_based_optional(t.three_prime_terminal_outer_min0),
        "three_prime_terminal_outer_max": _end_optional(t.three_prime_terminal_outer_max0) if t.strand == "+" else _one_based_optional(t.three_prime_terminal_outer_max0),
        "three_prime_terminal_outer_range_bp": t.three_prime_terminal_outer_range_bp,
        "flags": t.flags,
        "exon_chain": _chain_text(t.exon_parts),
        "intron_chain": _intron_chain_text(t.introns),
    } for t in sorted(collapsed, key=lambda x: (x.chrom, x.start0, x.end0, x.transcript_evidence_id))]
    columns = [
        "transcript_evidence_id", "transcript_evidence_superlocus_id", "transcript_evidence_sublocus_id",
        "chrom", "strand", "start", "end", "span_bp", "exonic_length", "exon_count", "intron_count",
        "collapse_key_type", "member_count", "source_count", "source_member_counts", "representative_raw_id",
        "left_terminal_outer_min", "left_terminal_outer_max", "left_terminal_outer_range_bp",
        "right_terminal_outer_min", "right_terminal_outer_max", "right_terminal_outer_range_bp",
        "five_prime_terminal_outer_min", "five_prime_terminal_outer_max", "five_prime_terminal_outer_range_bp",
        "three_prime_terminal_outer_min", "three_prime_terminal_outer_max", "three_prime_terminal_outer_range_bp",
        "flags", "exon_chain", "intron_chain",
    ]
    write_dict_tsv(rows, path, columns)


def write_transcript_superloci_tsv(superloci: list[TranscriptEvidenceSuperlocus], path: Path) -> None:
    rows = []
    for sl in superloci:
        counts = sl.raw_source_member_counts or dict(Counter(raw.source_label for raw in sl.raw_transcripts))
        raw_count = sl.raw_transcript_count or len(sl.raw_transcripts)
        rows.append({
            "transcript_evidence_superlocus_id": sl.evidence_superlocus_id,
            "chrom": sl.chrom,
            "strand": sl.strand,
            "start": sl.start0 + 1,
            "end": sl.end0,
            "span_bp": sl.end0 - sl.start0,
            "raw_transcript_count": raw_count,
            "collapsed_transcript_count": len(sl.collapsed_transcripts),
            "source_count": len(counts),
            "source_member_counts": fmt_counts(dict(counts)),
        })
    write_dict_tsv(rows, path, ["transcript_evidence_superlocus_id", "chrom", "strand", "start", "end", "span_bp", "raw_transcript_count", "collapsed_transcript_count", "source_count", "source_member_counts"])


def write_transcript_subloci_tsv(subloci: list[TranscriptEvidenceSublocus], path: Path) -> None:
    rows = [{
        "transcript_evidence_sublocus_id": sub.evidence_sublocus_id,
        "transcript_evidence_superlocus_id": sub.evidence_superlocus_id,
        "chrom": sub.chrom,
        "strand": sub.strand,
        "start": sub.start0 + 1,
        "end": sub.end0,
        "span_bp": sub.end0 - sub.start0,
        "sublocus_type": sub.sublocus_type,
        "collapsed_transcript_count": len(sub.collapsed_transcripts),
        "source_count": len({s for t in sub.collapsed_transcripts for s in t.sources}),
        "member_count": sum(t.member_count for t in sub.collapsed_transcripts),
        "flags": sub.flags,
    } for sub in subloci]
    write_dict_tsv(rows, path, ["transcript_evidence_sublocus_id", "transcript_evidence_superlocus_id", "chrom", "strand", "start", "end", "span_bp", "sublocus_type", "collapsed_transcript_count", "source_count", "member_count", "flags"])


def write_transcript_sublocus_members_tsv(subloci: list[TranscriptEvidenceSublocus], path: Path) -> None:
    rows = []
    for sub in subloci:
        for tx in sub.collapsed_transcripts:
            rows.append({
                "transcript_evidence_sublocus_id": sub.evidence_sublocus_id,
                "transcript_evidence_superlocus_id": sub.evidence_superlocus_id,
                "transcript_evidence_id": tx.transcript_evidence_id,
                "chrom": tx.chrom,
                "strand": tx.strand,
                "start": tx.start0 + 1,
                "end": tx.end0,
                "span_bp": tx.span_bp,
                "exon_count": tx.exon_count,
                "intron_count": tx.intron_count,
                "member_count": tx.member_count,
                "source_count": tx.source_count,
                "source_member_counts": fmt_counts(tx.source_member_counts),
            })
    write_dict_tsv(rows, path, ["transcript_evidence_sublocus_id", "transcript_evidence_superlocus_id", "transcript_evidence_id", "chrom", "strand", "start", "end", "span_bp", "exon_count", "intron_count", "member_count", "source_count", "source_member_counts"])


# ---------------------------------------------------------------------------
# Diagnostic model-to-transcript evidence association
# ---------------------------------------------------------------------------

TRANSCRIPT_ASSOCIATION_MIN_EXONIC_OVERLAP_BP = 30
TRANSCRIPT_ASSOCIATION_MIN_MODEL_CDS_FRACTION = 0.05
TRANSCRIPT_ASSOCIATION_MIN_SPAN_OVERLAP_BP = 150
TRANSCRIPT_ASSOCIATION_MIN_SPAN_FRACTION = 0.25
TRANSCRIPT_SINGLE_EXON_SUPPORT_MIN_MODEL_CDS_FRACTION = 0.80
TRANSCRIPT_SINGLE_EXON_SUPPORT_MIN_BP = 150
TRANSCRIPT_CHAIN_STRINGENT_MAX_EXTRA_5PRIME_INTRONS = 4
TRANSCRIPT_CHAIN_STRINGENT_MAX_EXTRA_3PRIME_INTRONS = 2
# Single-exon CDS models can be supported by simple multi-exon transcripts
# where the additional transcript exon(s) plausibly represent UTR sequence.
TRANSCRIPT_SINGLE_EXON_MULTIEXON_MAX_EXTRA_5PRIME_EXONS = 4
TRANSCRIPT_SINGLE_EXON_MULTIEXON_MAX_EXTRA_3PRIME_EXONS = 2
TRANSCRIPT_SUPPORTING_SOURCE_SCORE_THRESHOLD = 0.50

# Diagnostic-only cDNA-context flags. These thresholds do not affect selection
# or pruning; they only make potentially fragment-like transcript/model
# contexts easier to audit.
TRANSCRIPT_CDNA_CONTEXT_LOW_CDS_TO_MIN_CDNA_FRACTION = 0.30
TRANSCRIPT_CDNA_CONTEXT_VERY_LOW_CDS_TO_MIN_CDNA_FRACTION = 0.15


def _model_introns(model: CodingModel) -> frozenset[tuple[int, int]]:
    return model.cds_intron_set


def _model_splice_sites(model: CodingModel) -> set[int]:
    return model.cds_splice_sites


def _f1(precision: float | None, sensitivity: float | None) -> float:
    if precision is None or sensitivity is None:
        return 0.0
    if precision + sensitivity <= 0:
        return 0.0
    return 2.0 * precision * sensitivity / (precision + sensitivity)


def _span_overlap_bp(a_start0: int, a_end0: int, b_start0: int, b_end0: int) -> int:
    return max(0, min(a_end0, b_end0) - max(a_start0, b_start0))


def _contains_interval(container: tuple[int, int], item: tuple[int, int]) -> bool:
    return container[0] <= item[0] and container[1] >= item[1]


def _collapsed_transcript_effective_exon_parts(tx: TranscriptEvidenceCollapsedTranscript) -> tuple[tuple[int, int], ...]:
    """Return cached exon parts expanded to collapsed terminal outer ranges."""
    cached = tx._effective_exon_parts_cache
    if cached is not None:
        return cached
    parts = list(tx.sorted_exon_parts)
    if not parts:
        cached = tuple()
    elif len(parts) == 1:
        cached = ((tx.left_terminal_outer_min0, tx.right_terminal_outer_max0),)
    else:
        parts[0] = (tx.left_terminal_outer_min0, parts[0][1])
        parts[-1] = (parts[-1][0], tx.right_terminal_outer_max0)
        cached = tuple(parts)
    tx._effective_exon_parts_cache = cached
    return cached


def _exon_fully_covered_by_transcript(cds_exon: tuple[int, int], tx: TranscriptEvidenceCollapsedTranscript) -> bool:
    return any(_contains_interval(tx_exon, cds_exon) for tx_exon in _collapsed_transcript_effective_exon_parts(tx))


def _internal_cds_exons(model: CodingModel) -> tuple[tuple[int, int], ...]:
    return model.internal_cds_exons


def _terminal_cds_exons(model: CodingModel) -> tuple[tuple[int, int], ...]:
    return model.terminal_cds_exons


def _internal_cds_exon_supported_by_transcript(cds_exon: tuple[int, int], tx: TranscriptEvidenceCollapsedTranscript) -> bool:
    if not _exon_fully_covered_by_transcript(cds_exon, tx):
        return False
    # Internal CDS exons should be bounded by transcript splice sites.
    return cds_exon[0] in tx.splice_sites and cds_exon[1] in tx.splice_sites


def _terminal_cds_exon_supported_by_transcript(model: CodingModel, cds_exon: tuple[int, int], tx: TranscriptEvidenceCollapsedTranscript) -> bool:
    if not _exon_fully_covered_by_transcript(cds_exon, tx):
        return False
    parts = model.sorted_cds_parts
    if len(parts) <= 1:
        return True
    if cds_exon == parts[0]:
        # Left/genomic terminal CDS exon: splice-facing boundary is the right end.
        return cds_exon[1] in tx.splice_sites
    if cds_exon == parts[-1]:
        # Right/genomic terminal CDS exon: splice-facing boundary is the left end.
        return cds_exon[0] in tx.splice_sites
    return False


def _single_exon_limited_multiexon_support_metrics(
    model: CodingModel,
    tx: TranscriptEvidenceCollapsedTranscript,
) -> dict[str, object]:
    """Return limited multi-exon transcript support metrics for a single-exon CDS model.

    A single-exon CDS can be compatible with a simple multi-exon transcript when
    one transcript exon fully covers the CDS and the additional transcript exons
    are limited enough to plausibly represent UTR structure. More complex
    multi-exon overlaps remain broad/diagnostic only, because they may indicate
    that the single-exon CDS model is a fragment of a larger transcript-supported
    gene.
    """
    default = {
        "single_exon_limited_multiexon_support": False,
        "single_exon_multiexon_extra_5prime_exons": 0,
        "single_exon_multiexon_extra_3prime_exons": 0,
    }
    if not (model.is_single_exon and tx.is_multi_exon):
        return default
    if not model.cds_parts:
        return default
    cds_exon = model.sorted_cds_parts[0]
    tx_exons = list(sorted(_collapsed_transcript_effective_exon_parts(tx)))
    containing_indices = [i for i, tx_exon in enumerate(tx_exons) if _contains_interval(tx_exon, cds_exon)]
    if not containing_indices:
        return default

    # If multiple effective exons could contain the CDS because of unusual
    # overlaps, use the one with the fewest total extra exons.
    best: tuple[int, int, int] | None = None
    strand = model.strand if model.strand in {"+", "-"} else tx.strand
    for idx in containing_indices:
        extra_left = idx
        extra_right = len(tx_exons) - idx - 1
        if strand == "-":
            extra_5, extra_3 = extra_right, extra_left
        else:
            # For model '+' or unstranded transcript evidence paired with a
            # stranded model, genomic-left is the 5' side. If strand cannot be
            # resolved, this conservative genomic-left/right convention is still
            # deterministic and keeps the association auditable.
            extra_5, extra_3 = extra_left, extra_right
        candidate = (extra_5 + extra_3, extra_5, extra_3)
        if best is None or candidate < best:
            best = candidate
    if best is None:
        return default
    _total, extra_5, extra_3 = best
    ok = (
        extra_5 <= TRANSCRIPT_SINGLE_EXON_MULTIEXON_MAX_EXTRA_5PRIME_EXONS
        and extra_3 <= TRANSCRIPT_SINGLE_EXON_MULTIEXON_MAX_EXTRA_3PRIME_EXONS
    )
    return {
        "single_exon_limited_multiexon_support": ok,
        "single_exon_multiexon_extra_5prime_exons": extra_5,
        "single_exon_multiexon_extra_3prime_exons": extra_3,
    }


def _supported_model_exons(
    model: CodingModel,
    transcripts: list[TranscriptEvidenceCollapsedTranscript],
    *,
    terminal: bool,
) -> set[tuple[int, int]]:
    exons = _terminal_cds_exons(model) if terminal else _internal_cds_exons(model)
    supported: set[tuple[int, int]] = set()
    for exon in exons:
        for tx in transcripts:
            ok = (
                _terminal_cds_exon_supported_by_transcript(model, exon, tx)
                if terminal
                else _internal_cds_exon_supported_by_transcript(exon, tx)
            )
            if ok:
                supported.add(exon)
                break
    return supported


def _contiguous_subchain_index(needle: tuple[tuple[int, int], ...], haystack: tuple[tuple[int, int], ...]) -> int | None:
    if not needle or len(needle) > len(haystack):
        return None
    n = len(needle)
    for i in range(0, len(haystack) - n + 1):
        if haystack[i:i + n] == needle:
            return i
    return None


def _chain_relation_metrics(model: CodingModel, tx: TranscriptEvidenceCollapsedTranscript) -> dict[str, object]:
    model_chain = tuple(model.cds_introns)
    tx_chain = tuple(tx.introns)
    if not model_chain or not tx_chain:
        return {
            "exact_chain_match": False,
            "model_chain_contained": False,
            "model_chain_stringent_contained": False,
            "extra_5prime_transcript_introns": 0,
            "extra_3prime_transcript_introns": 0,
        }
    exact = model_chain == tx_chain
    idx = _contiguous_subchain_index(model_chain, tx_chain)
    contained = idx is not None
    extra_5 = 0
    extra_3 = 0
    stringent = False
    if contained and idx is not None:
        extra_left = idx
        extra_right = len(tx_chain) - (idx + len(model_chain))
        strand = model.strand if model.strand in {"+", "-"} else tx.strand
        if strand == "-":
            extra_5, extra_3 = extra_right, extra_left
        else:
            extra_5, extra_3 = extra_left, extra_right
        stringent = (
            extra_5 <= TRANSCRIPT_CHAIN_STRINGENT_MAX_EXTRA_5PRIME_INTRONS
            and extra_3 <= TRANSCRIPT_CHAIN_STRINGENT_MAX_EXTRA_3PRIME_INTRONS
        )
    return {
        "exact_chain_match": exact,
        "model_chain_contained": contained,
        "model_chain_stringent_contained": stringent,
        "extra_5prime_transcript_introns": extra_5,
        "extra_3prime_transcript_introns": extra_3,
    }


def _collapsed_transcript_min_exon_parts(tx: TranscriptEvidenceCollapsedTranscript) -> tuple[tuple[int, int], ...]:
    """Return the conservative minimum-cDNA exon representation.

    This calculation is used only for the selected cDNA-context candidate, so
    retaining another per-transcript cache costs more memory than it saves.
    """
    parts = list(tx.sorted_exon_parts)
    if not parts:
        return tuple()
    elif len(parts) == 1:
        start = tx.left_terminal_outer_max0
        end = tx.right_terminal_outer_min0
        return ((start, end),) if end > start else tuple(parts)
    else:
        left_start = tx.left_terminal_outer_max0
        left_end = parts[0][1]
        right_start = parts[-1][0]
        right_end = tx.right_terminal_outer_min0
        if left_end > left_start:
            parts[0] = (left_start, left_end)
        if right_end > right_start:
            parts[-1] = (right_start, right_end)
        return tuple(parts)


def _exonic_bp_left_of(intervals: tuple[tuple[int, int], ...] | list[tuple[int, int]], boundary0: int) -> int:
    return sum(max(0, min(e, boundary0) - s) for s, e in intervals if s < boundary0)


def _exonic_bp_right_of(intervals: tuple[tuple[int, int], ...] | list[tuple[int, int]], boundary0: int) -> int:
    return sum(max(0, e - max(s, boundary0)) for s, e in intervals if e > boundary0)


def _strand_or_model_strand(model: CodingModel, tx: TranscriptEvidenceCollapsedTranscript) -> str:
    return model.strand if model.strand in {"+", "-"} else tx.strand


def _genomic_to_five_three(model: CodingModel, tx: TranscriptEvidenceCollapsedTranscript, left_value: int, right_value: int) -> tuple[int, int]:
    strand = _strand_or_model_strand(model, tx)
    if strand == "-":
        return right_value, left_value
    return left_value, right_value


def _cdna_context_default(model: CodingModel, *, status: str = "not_calculated_no_compatible_transcript") -> dict[str, object]:
    return {
        "transcript_cdna_context_status": status,
        "transcript_cdna_context_best_evidence_id": ".",
        "transcript_cdna_context_best_association_class": ".",
        "transcript_cdna_context_best_matches_best_pairwise": False,
        "transcript_cdna_context_best_source_count": 0,
        "transcript_cdna_context_best_member_count": 0,
        "transcript_cdna_context_model_cds_length_bp": model.cds_length,
        "transcript_cdna_context_cdna_length_min_bp": ".",
        "transcript_cdna_context_cdna_length_max_bp": ".",
        "transcript_cdna_context_cds_to_min_cdna_fraction": ".",
        "transcript_cdna_context_cds_to_max_cdna_fraction": ".",
        "transcript_cdna_context_extra_5prime_exon_count": ".",
        "transcript_cdna_context_extra_3prime_exon_count": ".",
        "transcript_cdna_context_extra_5prime_exonic_bp": ".",
        "transcript_cdna_context_extra_3prime_exonic_bp": ".",
        "transcript_cdna_context_extra_5prime_intron_count": ".",
        "transcript_cdna_context_extra_3prime_intron_count": ".",
        "transcript_cdna_context_flags": ".",
    }


def _single_exon_cdna_context_metrics(
    model: CodingModel,
    tx: TranscriptEvidenceCollapsedTranscript,
    metrics: dict[str, object],
) -> dict[str, object] | None:
    if not model.is_single_exon or not model.cds_parts:
        return None
    cds_exon = model.sorted_cds_parts[0]
    tx_exons_max = list(sorted(_collapsed_transcript_effective_exon_parts(tx)))
    containing_indices = [i for i, tx_exon in enumerate(tx_exons_max) if _contains_interval(tx_exon, cds_exon)]
    if not containing_indices:
        return None
    # Prefer the containing exon with the fewest outside exons, then the shortest
    # containing exon. This keeps the context deterministic for unusual records.
    best_idx = min(
        containing_indices,
        key=lambda i: ((i + (len(tx_exons_max) - i - 1)), tx_exons_max[i][1] - tx_exons_max[i][0], i),
    )
    extra_left_exons = best_idx
    extra_right_exons = len(tx_exons_max) - best_idx - 1
    extra_5_exons, extra_3_exons = _genomic_to_five_three(model, tx, extra_left_exons, extra_right_exons)
    extra_left_bp = _exonic_bp_left_of(tx_exons_max, model.cds_start0)
    extra_right_bp = _exonic_bp_right_of(tx_exons_max, model.cds_end0)
    extra_5_bp, extra_3_bp = _genomic_to_five_three(model, tx, extra_left_bp, extra_right_bp)
    return {
        "compatible_class": str(metrics.get("association_class", ".")),
        "extra_5prime_exon_count": extra_5_exons,
        "extra_3prime_exon_count": extra_3_exons,
        "extra_5prime_exonic_bp": extra_5_bp,
        "extra_3prime_exonic_bp": extra_3_bp,
        "extra_5prime_intron_count": 0,
        "extra_3prime_intron_count": 0,
        "matched_transcript_exon_indices": (best_idx,),
    }


def _multi_exon_cdna_context_metrics(
    model: CodingModel,
    tx: TranscriptEvidenceCollapsedTranscript,
    metrics: dict[str, object],
) -> dict[str, object] | None:
    if not (model.is_multi_exon and tx.is_multi_exon):
        return None
    model_chain = tuple(model.cds_introns)
    tx_chain = tuple(tx.introns)
    idx = _contiguous_subchain_index(model_chain, tx_chain)
    if idx is None:
        return None
    tx_exons_max = list(sorted(_collapsed_transcript_effective_exon_parts(tx)))
    matched_exon_count = len(model_chain) + 1
    model_exons = list(model.sorted_cds_parts)
    matched_tx_indices = tuple(range(idx, idx + matched_exon_count))
    if len(model_exons) != matched_exon_count or matched_tx_indices[-1] >= len(tx_exons_max):
        return None
    matched_tx_exons = [tx_exons_max[i] for i in matched_tx_indices]
    # A contained/exact intron chain alone is not sufficient for cDNA-context
    # diagnostics. A fragment can share the complete CDS intron chain but fail
    # to cover one or both terminal CDS exon extents. Qualification is therefore
    # based on the outermost collapsed transcript envelope, and every CDS exon
    # must be contained by its corresponding transcript exon.
    if not all(_contains_interval(tx_exon, cds_exon) for cds_exon, tx_exon in zip(model_exons, matched_tx_exons)):
        return None
    extra_left_exons = idx
    extra_right_exons = max(0, len(tx_exons_max) - (idx + matched_exon_count))
    extra_left_introns = idx
    extra_right_introns = max(0, len(tx_chain) - (idx + len(model_chain)))
    extra_5_exons, extra_3_exons = _genomic_to_five_three(model, tx, extra_left_exons, extra_right_exons)
    extra_5_introns, extra_3_introns = _genomic_to_five_three(model, tx, extra_left_introns, extra_right_introns)
    extra_left_bp = _exonic_bp_left_of(tx_exons_max, model.cds_start0)
    extra_right_bp = _exonic_bp_right_of(tx_exons_max, model.cds_end0)
    extra_5_bp, extra_3_bp = _genomic_to_five_three(model, tx, extra_left_bp, extra_right_bp)
    if bool(metrics.get("exact_intron_chain_match")):
        compatible_class = "exact_intron_chain_match"
    elif bool(metrics.get("model_intron_chain_stringent_contained")):
        compatible_class = "model_intron_chain_stringent_contained"
    else:
        compatible_class = "model_intron_chain_contained"
    return {
        "compatible_class": compatible_class,
        "extra_5prime_exon_count": extra_5_exons,
        "extra_3prime_exon_count": extra_3_exons,
        "extra_5prime_exonic_bp": extra_5_bp,
        "extra_3prime_exonic_bp": extra_3_bp,
        "extra_5prime_intron_count": extra_5_introns,
        "extra_3prime_intron_count": extra_3_introns,
        "matched_transcript_exon_indices": matched_tx_indices,
    }


def _condition_min_exon_parts_on_cds(
    model: CodingModel,
    min_parts: tuple[tuple[int, int], ...],
    max_parts: tuple[tuple[int, int], ...],
    matched_transcript_exon_indices: tuple[int, ...],
) -> tuple[tuple[tuple[int, int], ...], bool] | None:
    """Condition the collapsed minimum cDNA envelope on CDS containment.

    The collapsed minimum uses innermost observed terminal boundaries, while
    eligibility uses the outermost envelope. When the innermost envelope is
    shorter than the CDS but the outer envelope contains it, extend only the
    matched minimum-envelope exons as far as required to contain the CDS. Extra
    transcript exons and unavoidable non-CDS sequence are retained.
    """
    model_exons = list(model.sorted_cds_parts)
    if len(model_exons) != len(matched_transcript_exon_indices):
        return None
    if len(min_parts) != len(max_parts):
        return None
    conditioned = list(min_parts)
    changed = False
    for cds_exon, tx_index in zip(model_exons, matched_transcript_exon_indices):
        if tx_index < 0 or tx_index >= len(max_parts):
            return None
        max_exon = max_parts[tx_index]
        if not _contains_interval(max_exon, cds_exon):
            return None
        min_exon = conditioned[tx_index]
        adjusted = (min(min_exon[0], cds_exon[0]), max(min_exon[1], cds_exon[1]))
        if not _contains_interval(max_exon, adjusted):
            return None
        if adjusted != min_exon:
            changed = True
            conditioned[tx_index] = adjusted
    return tuple(conditioned), changed


def _validate_cdna_context_invariants(
    model: CodingModel,
    tx: TranscriptEvidenceCollapsedTranscript,
    *,
    cdna_min: int,
    cdna_max: int,
    cds_to_min: float,
    cds_to_max: float,
) -> None:
    """Raise on an internally inconsistent calculated cDNA context."""
    valid = (
        model.cds_length > 0
        and model.cds_length <= cdna_min <= cdna_max
        and 0.0 < cds_to_max <= cds_to_min <= 1.0 + 1e-12
    )
    if not valid:
        raise AssertionError(
            "invalid transcript cDNA context invariant for "
            f"model={model.output_transcript_id} evidence={tx.transcript_evidence_id}: "
            f"cds={model.cds_length}, min={cdna_min}, max={cdna_max}, "
            f"fractions={cds_to_min},{cds_to_max}"
        )


def _cdna_context_candidate_metrics(
    model: CodingModel,
    tx: TranscriptEvidenceCollapsedTranscript,
    metrics: dict[str, object],
) -> dict[str, object] | None:
    if model.is_single_exon:
        return _single_exon_cdna_context_metrics(model, tx, metrics)
    return _multi_exon_cdna_context_metrics(model, tx, metrics)


def _cdna_context_priority(model: CodingModel, tx: TranscriptEvidenceCollapsedTranscript, metrics: dict[str, object], context: dict[str, object]) -> tuple[float, float, int, int]:
    cls = str(context.get("compatible_class", "."))
    if model.is_single_exon:
        priority = {
            "single_exon_transcript_support": 100.0,
            "single_exon_limited_multiexon_support": 90.0,
        }.get(cls, 80.0)
    else:
        priority = {
            "exact_intron_chain_match": 100.0,
            "model_intron_chain_stringent_contained": 90.0,
            "model_intron_chain_contained": 80.0,
        }.get(cls, 0.0)
    score_value = metrics.get("pairwise_transcript_score")
    if score_value is None:
        score_value = metrics.get("association_score")
    if score_value is None:
        score_value = 0.0
    return (
        priority,
        float(score_value),
        tx.source_count,
        tx.member_count,
    )


def cdna_context_for_model(
    model: CodingModel,
    associations: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]],
    *,
    best_pairwise_id: str = ".",
) -> dict[str, object]:
    """Return cDNA-context metrics for the best compatible collapsed transcript.

    Eligibility is evaluated against the outermost collapsed transcript
    envelope. A single-exon CDS must be fully covered by one transcript exon; a
    multi-exon CDS intron chain must be a contiguous transcript subchain and
    every CDS exon must be contained by its corresponding transcript exon. The
    minimum cDNA envelope is then conditioned on CDS containment, so calculated
    CDS/cDNA fractions cannot exceed 1. Terminal extrema may originate from
    different members of the same exact-chain collapsed group.
    """
    best_item: tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object], dict[str, object], tuple[float, float, int, int]] | None = None
    for tx, metrics in associations:
        context = _cdna_context_candidate_metrics(model, tx, metrics)
        if context is None:
            continue
        priority = _cdna_context_priority(model, tx, metrics, context)
        if best_item is None or priority > best_item[3]:
            best_item = (tx, metrics, context, priority)
    if best_item is None:
        if not associations:
            return _cdna_context_default(model, status="not_calculated_no_associated_transcript")
        return _cdna_context_default(model, status="not_calculated_no_compatible_transcript")

    tx, _metrics, context, _priority = best_item
    raw_min_parts = _collapsed_transcript_min_exon_parts(tx)
    max_parts = _collapsed_transcript_effective_exon_parts(tx)
    matched_indices = tuple(int(i) for i in context.get("matched_transcript_exon_indices", ()))
    conditioned_result = _condition_min_exon_parts_on_cds(model, raw_min_parts, max_parts, matched_indices)
    if conditioned_result is None:
        return _cdna_context_default(model, status="not_calculated_no_compatible_transcript")
    min_parts, min_was_conditioned = conditioned_result
    cdna_min = interval_total_bp(list(min_parts))
    cdna_max = interval_total_bp(list(max_parts))
    cds_to_min = _safe_fraction(model.cds_length, cdna_min)
    cds_to_max = _safe_fraction(model.cds_length, cdna_max)
    _validate_cdna_context_invariants(
        model, tx, cdna_min=cdna_min, cdna_max=cdna_max,
        cds_to_min=cds_to_min, cds_to_max=cds_to_max,
    )
    flags: list[str] = []
    if min_was_conditioned:
        flags.append("min_cdna_conditioned_on_cds_containment")
    if cds_to_min < TRANSCRIPT_CDNA_CONTEXT_VERY_LOW_CDS_TO_MIN_CDNA_FRACTION:
        flags.append("very_low_cds_to_cdna_fraction")
    elif cds_to_min < TRANSCRIPT_CDNA_CONTEXT_LOW_CDS_TO_MIN_CDNA_FRACTION:
        flags.append("low_cds_to_cdna_fraction")
    if int(context.get("extra_5prime_exon_count") or 0) > TRANSCRIPT_SINGLE_EXON_MULTIEXON_MAX_EXTRA_5PRIME_EXONS:
        flags.append("many_extra_5prime_transcript_exons")
    if int(context.get("extra_3prime_exon_count") or 0) > TRANSCRIPT_SINGLE_EXON_MULTIEXON_MAX_EXTRA_3PRIME_EXONS:
        flags.append("many_extra_3prime_transcript_exons")
    if model.is_single_exon and tx.is_multi_exon:
        flags.append("single_exon_cds_in_multiexon_cdna_context")
    if model.is_multi_exon and bool(_metrics.get("model_intron_chain_contained")) and not bool(_metrics.get("exact_intron_chain_match")):
        flags.append("cds_chain_contained_in_longer_cdna")
    return {
        "transcript_cdna_context_status": "calculated",
        "transcript_cdna_context_best_evidence_id": tx.transcript_evidence_id,
        "transcript_cdna_context_best_association_class": str(context.get("compatible_class", ".")),
        "transcript_cdna_context_best_matches_best_pairwise": tx.transcript_evidence_id == best_pairwise_id,
        "transcript_cdna_context_best_source_count": tx.source_count,
        "transcript_cdna_context_best_member_count": tx.member_count,
        "transcript_cdna_context_model_cds_length_bp": model.cds_length,
        "transcript_cdna_context_cdna_length_min_bp": cdna_min,
        "transcript_cdna_context_cdna_length_max_bp": cdna_max,
        "transcript_cdna_context_cds_to_min_cdna_fraction": cds_to_min,
        "transcript_cdna_context_cds_to_max_cdna_fraction": cds_to_max,
        "transcript_cdna_context_extra_5prime_exon_count": int(context.get("extra_5prime_exon_count") or 0),
        "transcript_cdna_context_extra_3prime_exon_count": int(context.get("extra_3prime_exon_count") or 0),
        "transcript_cdna_context_extra_5prime_exonic_bp": int(context.get("extra_5prime_exonic_bp") or 0),
        "transcript_cdna_context_extra_3prime_exonic_bp": int(context.get("extra_3prime_exonic_bp") or 0),
        "transcript_cdna_context_extra_5prime_intron_count": int(context.get("extra_5prime_intron_count") or 0),
        "transcript_cdna_context_extra_3prime_intron_count": int(context.get("extra_3prime_intron_count") or 0),
        "transcript_cdna_context_flags": fmt_flags(flags),
    }

def _source_member_summary_for_transcripts(transcripts: list[TranscriptEvidenceCollapsedTranscript]) -> tuple[int, int]:
    sources = {source for tx in transcripts for source in tx.sources}
    member_count = sum(tx.member_count for tx in transcripts)
    return len(sources), member_count


def _mean_best_source_metric(
    associations: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]],
    metric_name: str,
) -> float:
    best_by_source: dict[str, float] = {}
    for tx, metrics in associations:
        value = float(metrics.get(metric_name) or 0.0)
        for source in tx.sources:
            if value > best_by_source.get(source, -1.0):
                best_by_source[source] = value
    return _mean(list(best_by_source.values())) if best_by_source else 0.0


def _compatible_transcript_strand(model: CodingModel, tx: TranscriptEvidenceCollapsedTranscript) -> bool:
    return tx.strand == model.strand or tx.strand == "." or model.strand == "."


def transcript_model_pair_metrics(model: CodingModel, tx: TranscriptEvidenceCollapsedTranscript) -> dict[str, object]:
    """Classify one model/transcript-evidence relationship.

    Broad associations are retained for audit and context diagnostics, while
    structurally tighter relationships are marked scoring-eligible and can
    contribute to transcript scoring, rescue, and evidence-aware pruning.
    """
    model_introns = _model_introns(model)
    tx_introns = tx.intron_set
    model_splice_sites = _model_splice_sites(model)
    tx_splice_sites = tx.splice_sites
    shared_introns = model_introns & tx_introns
    shared_splice_sites = model_splice_sites & tx_splice_sites

    exonic_overlap_bp = interval_list_overlap_bp_sorted(model.sorted_cds_parts, tx.sorted_exon_parts)
    model_cds_fraction = exonic_overlap_bp / model.cds_length if model.cds_length else 0.0
    tx_exonic_length = tx.exonic_length
    transcript_exonic_fraction = exonic_overlap_bp / tx_exonic_length if tx_exonic_length else 0.0
    span_overlap_bp = _span_overlap_bp(model.cds_start0, model.cds_end0, tx.start0, tx.end0)
    span_fraction = span_overlap_bp / model.cds_span_bp if model.cds_span_bp else 0.0
    transcript_span_contains_model = tx.start0 <= model.cds_start0 and tx.end0 >= model.cds_end0
    model_span_contains_transcript = model.cds_start0 <= tx.start0 and model.cds_end0 >= tx.end0

    exonic_threshold_met = (
        exonic_overlap_bp >= TRANSCRIPT_ASSOCIATION_MIN_EXONIC_OVERLAP_BP
        or model_cds_fraction >= TRANSCRIPT_ASSOCIATION_MIN_MODEL_CDS_FRACTION
    )
    span_threshold_met = (
        span_overlap_bp >= TRANSCRIPT_ASSOCIATION_MIN_SPAN_OVERLAP_BP
        or span_fraction >= TRANSCRIPT_ASSOCIATION_MIN_SPAN_FRACTION
        or transcript_span_contains_model
        or model_span_contains_transcript
    )
    single_model_min_bp = min(TRANSCRIPT_SINGLE_EXON_SUPPORT_MIN_BP, model.cds_length) if model.cds_length else TRANSCRIPT_SINGLE_EXON_SUPPORT_MIN_BP
    single_exon_support_met = (
        model.is_single_exon
        and tx.is_single_exon
        and model_cds_fraction >= TRANSCRIPT_SINGLE_EXON_SUPPORT_MIN_MODEL_CDS_FRACTION
        and exonic_overlap_bp >= single_model_min_bp
    )
    limited_multiexon_metrics = _single_exon_limited_multiexon_support_metrics(model, tx)
    single_exon_limited_multiexon_support_met = bool(
        limited_multiexon_metrics.get("single_exon_limited_multiexon_support")
    )

    association_class = "."
    scoring_eligible = False
    if model.is_multi_exon and tx.is_multi_exon and shared_introns:
        association_class = "shared_intron"
        scoring_eligible = True
    elif model.is_multi_exon and tx.is_multi_exon and shared_splice_sites:
        association_class = "shared_splice_site"
        scoring_eligible = True
    elif single_exon_support_met:
        association_class = "single_exon_transcript_support"
        scoring_eligible = True
    elif single_exon_limited_multiexon_support_met:
        association_class = "single_exon_limited_multiexon_support"
        scoring_eligible = True
    elif model.is_single_exon and tx.is_multi_exon and exonic_threshold_met:
        association_class = "single_exon_model_multiexon_transcript_overlap"
    elif model.is_multi_exon and tx.is_single_exon and exonic_threshold_met:
        association_class = "multiexon_model_single_exon_transcript_overlap"
    elif exonic_threshold_met:
        association_class = "model_cds_exonic_overlap"
    elif transcript_span_contains_model and span_threshold_met:
        association_class = "span_contains_model"
    elif model_span_contains_transcript and span_threshold_met:
        association_class = "model_contains_transcript_span"
    elif span_threshold_met:
        association_class = "span_overlap_only"

    model_intron_precision = len(shared_introns) / len(model_introns) if model_introns else None
    transcript_intron_sensitivity = len(shared_introns) / len(tx_introns) if tx_introns else None
    intron_f1 = _f1(model_intron_precision, transcript_intron_sensitivity)
    model_splice_precision = len(shared_splice_sites) / len(model_splice_sites) if model_splice_sites else None
    transcript_splice_sensitivity = len(shared_splice_sites) / len(tx_splice_sites) if tx_splice_sites else None
    splice_site_f1 = _f1(model_splice_precision, transcript_splice_sensitivity)

    chain_metrics = _chain_relation_metrics(model, tx)
    pairwise_components = _pairwise_transcript_score_components(model, tx, {
        "model_intron_precision": 0.0 if model_intron_precision is None else model_intron_precision,
        "splice_site_f1": splice_site_f1,
        "intron_f1": intron_f1,
        "model_cds_overlap_fraction": model_cds_fraction,
    })

    if bool(scoring_eligible):
        association_score = pairwise_components["pairwise_transcript_score"]
    elif model.is_single_exon and (tx.is_single_exon or single_exon_limited_multiexon_support_met):
        association_score = model_cds_fraction
    elif model.is_multi_exon and tx.is_multi_exon:
        association_score = 0.50 * intron_f1 + 0.35 * splice_site_f1 + 0.15 * model_cds_fraction
    else:
        association_score = 0.25 * model_cds_fraction

    flags: list[str] = []
    if tx.strand == ".":
        flags.append("unstranded_transcript_evidence")
    if model.is_single_exon and tx.is_multi_exon and exonic_threshold_met:
        flags.append("single_exon_model_overlapped_by_multiexon_transcript")
        if single_exon_limited_multiexon_support_met:
            flags.append("single_exon_limited_multiexon_support")
        elif model_cds_fraction >= TRANSCRIPT_SINGLE_EXON_SUPPORT_MIN_MODEL_CDS_FRACTION:
            flags.append("single_exon_model_possible_fragment")
    if model.is_multi_exon and tx.is_single_exon and exonic_threshold_met:
        flags.append("multiexon_model_overlapped_by_single_exon_transcript")
    if association_class.startswith("span") or association_class in {"model_contains_transcript_span"}:
        flags.append("span_only_diagnostic")
    if transcript_span_contains_model and tx.is_multi_exon and model.is_single_exon:
        flags.append("transcript_span_contains_single_exon_model")

    return {
        "association_class": association_class,
        "scoring_eligible": scoring_eligible,
        "exonic_overlap_bp": exonic_overlap_bp,
        "model_cds_overlap_fraction": model_cds_fraction,
        "transcript_exonic_overlap_fraction": transcript_exonic_fraction,
        "span_overlap_bp": span_overlap_bp,
        "span_overlap_fraction_of_model_cds_span": span_fraction,
        "transcript_span_contains_model": transcript_span_contains_model,
        "model_span_contains_transcript": model_span_contains_transcript,
        "shared_intron_count": len(shared_introns),
        "shared_splice_site_count": len(shared_splice_sites),
        "model_intron_precision": 0.0 if model_intron_precision is None else model_intron_precision,
        "transcript_intron_sensitivity": 0.0 if transcript_intron_sensitivity is None else transcript_intron_sensitivity,
        "intron_f1": intron_f1,
        "model_splice_site_precision": 0.0 if model_splice_precision is None else model_splice_precision,
        "transcript_splice_site_sensitivity": 0.0 if transcript_splice_sensitivity is None else transcript_splice_sensitivity,
        "splice_site_f1": splice_site_f1,
        "exact_intron_chain_match": bool(chain_metrics["exact_chain_match"]),
        "model_intron_chain_contained": bool(chain_metrics["model_chain_contained"]),
        "model_intron_chain_stringent_contained": bool(chain_metrics["model_chain_stringent_contained"]),
        "extra_5prime_transcript_introns": int(chain_metrics["extra_5prime_transcript_introns"]),
        "extra_3prime_transcript_introns": int(chain_metrics["extra_3prime_transcript_introns"]),
        "single_exon_limited_multiexon_support": single_exon_limited_multiexon_support_met,
        "single_exon_multiexon_extra_5prime_exons": int(limited_multiexon_metrics.get("single_exon_multiexon_extra_5prime_exons") or 0),
        "single_exon_multiexon_extra_3prime_exons": int(limited_multiexon_metrics.get("single_exon_multiexon_extra_3prime_exons") or 0),
        "association_score": association_score,
        "pairwise_transcript_score": pairwise_components["pairwise_transcript_score"],
        "pairwise_model_cds_intron_support": pairwise_components["model_cds_intron_support"],
        "pairwise_splice_site_f1": pairwise_components["splice_site_f1"],
        "pairwise_full_intron_chain_f1": pairwise_components["full_intron_chain_f1"],
        "pairwise_internal_cds_exon_compatibility": pairwise_components["internal_cds_exon_compatibility"],
        "pairwise_terminal_cds_exon_compatibility": pairwise_components["terminal_cds_exon_compatibility"],
        "pairwise_model_cds_exonic_coverage": pairwise_components["model_cds_exonic_coverage"],
        "flags": fmt_flags(flags),
    }


def associated_transcripts_for_model(
    model: CodingModel,
    superlocus_index: IntervalOverlapIndex[str, TranscriptEvidenceSuperlocus],
) -> list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]]:
    """Return broad transcript associations using an indexed locus lookup."""
    out: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]] = []
    for sl in superlocus_index.query(model.chrom, model.cds_start0, model.cds_end0):
        for tx in sl.collapsed_transcripts:
            if tx.end0 <= model.cds_start0 or tx.start0 >= model.cds_end0:
                continue
            if not _compatible_transcript_strand(model, tx):
                continue
            metrics = transcript_model_pair_metrics(model, tx)
            if metrics["association_class"] != ".":
                out.append((tx, metrics))
    return out


def evaluate_coherent_exact_transcript_support(
    model: CodingModel,
    associations: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]],
    *,
    min_cdna_context_fraction: float,
    require_stranded_support: bool,
    block_on_cdna_context_concern: bool,
    cdna_context_neutral_fraction: float,
    cdna_context_severe_fraction: float,
) -> dict[str, object]:
    """Evaluate Contract-v1 exact transcript support per coherent group.

    Exact-chain, scoring, strand and cDNA-context requirements are evaluated on
    each individual collapsed transcript evidence group.  Aggregate model rows
    are deliberately not used to assemble a passing conjunction from unrelated
    records.
    """
    exact_items: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]] = []
    qualifying: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object], dict[str, object]]] = []
    context_states: list[str] = []
    exact_scoring_ok = False
    exact_stranded_ok = False
    for tx, metrics in associations:
        if not bool(metrics.get("exact_intron_chain_match")):
            continue
        exact_items.append((tx, metrics))
        scoring_eligible = bool(metrics.get("scoring_eligible"))
        stranded = tx.strand != "."
        stranded_ok = stranded or not require_stranded_support
        exact_scoring_ok = exact_scoring_ok or scoring_eligible
        exact_stranded_ok = exact_stranded_ok or stranded_ok
        context = cdna_context_for_model(model, [(tx, metrics)])
        concern = transcript_cdna_context_concern_fields(
            context,
            neutral_fraction=cdna_context_neutral_fraction,
            severe_fraction=cdna_context_severe_fraction,
        )
        context_status = str(context.get("transcript_cdna_context_status", "."))
        context_evaluable = str(concern.get("transcript_cdna_context_evaluable", "false")) == "true"
        if context_evaluable:
            fraction = float(context.get("transcript_cdna_context_cds_to_min_cdna_fraction", 0.0) or 0.0)
            fraction_ok = fraction + 1e-15 >= min_cdna_context_fraction
            concern_blocks = bool(
                block_on_cdna_context_concern
                and str(concern.get("transcript_cdna_context_pruning_concern", "false")) == "true"
            )
            context_ok = fraction_ok and not concern_blocks
            if not fraction_ok:
                context_state = "evaluated_below_minimum"
            elif concern_blocks:
                context_state = "evaluated_incompatible_context_concern"
            else:
                context_state = "evaluated_acceptable"
        else:
            # This is the existing MoSelecta rescue/scoring interpretation:
            # genuinely uncalculated context is neutral and the exact coherent
            # association must still satisfy every other requirement.
            context_ok = True
            context_state = f"{context_status}_neutral"
        context_states.append(context_state)
        detail = {
            "transcript_evidence_id": tx.transcript_evidence_id,
            "sources": sorted(tx.sources),
            "member_count": tx.member_count,
            "scoring_eligible": scoring_eligible,
            "stranded": stranded,
            "stranded_ok": stranded_ok,
            "cdna_context_status": context_status,
            "cdna_context_state": context_state,
            "cdna_context_ok": context_ok,
        }
        if scoring_eligible and stranded_ok and context_ok:
            qualifying.append((tx, metrics, detail))

    qualifying_ids = {tx.transcript_evidence_id for tx, _metrics, _detail in qualifying}
    qualifying_sources = {source for tx, _metrics, _detail in qualifying for source in tx.sources}
    exact_sources = {source for tx, _metrics in exact_items for source in tx.sources}
    if qualifying:
        status = "qualifying_exact_coherent_transcript_support"
    elif exact_items:
        status = "exact_chain_present_but_coherent_requirements_failed"
    elif associations:
        status = "associated_transcript_evidence_without_exact_chain"
    else:
        status = "no_associated_transcript_evidence"
    return {
        "exact_coherent_transcript_support": bool(qualifying),
        "exact_coherent_transcript_support_status": status,
        "exact_transcript_chain_group_count": len(exact_items),
        "exact_transcript_chain_source_count": len(exact_sources),
        "exact_transcript_chain_sources": ",".join(sorted(exact_sources)) if exact_sources else ".",
        "qualifying_exact_transcript_group_count": len(qualifying_ids),
        "qualifying_exact_transcript_source_count": len(qualifying_sources),
        "qualifying_exact_transcript_sources": ",".join(sorted(qualifying_sources)) if qualifying_sources else ".",
        "qualifying_exact_transcript_group_ids": ",".join(sorted(qualifying_ids)) if qualifying_ids else ".",
        "exact_transcript_scoring_eligibility_status": (
            "satisfied" if exact_scoring_ok else
            "failed" if exact_items else
            "no_exact_group"
        ),
        "exact_transcript_stranded_requirement": "required" if require_stranded_support else "not_required",
        "exact_transcript_stranded_status": (
            "satisfied" if exact_stranded_ok else
            "failed" if exact_items and require_stranded_support else
            "not_applicable"
        ),
        "exact_transcript_cdna_context_status": (
            ",".join(sorted(set(context_states))) if context_states else "."
        ),
    }


def _fmt_id_list(values: list[str] | set[str], *, max_items: int = 25) -> str:
    ordered = sorted(values)
    if not ordered:
        return "."
    if len(ordered) <= max_items:
        return ",".join(ordered)
    return ",".join(ordered[:max_items]) + f",...(+{len(ordered) - max_items})"


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def scale_transcript_evidence_score(score: float) -> float:
    """Return the current 0-1 transcript evidence score scale.

    The first transcript score is a bounded weighted mean of 0-1 components, so
    the scaled value is currently identical to the raw score. Keeping this helper
    mirrors protein evidence and leaves room for future calibration.
    """
    if score <= 0:
        return 0.0
    return min(score, 1.0)


def _safe_fraction(numerator: int | float, denominator: int | float) -> float:
    return float(numerator) / float(denominator) if denominator else 0.0


def _log_member_support_scaled(member_count: int, *, cap: int = 10) -> float:
    """Return capped/log-scaled member support in [0, 1]."""
    import math
    if member_count <= 0:
        return 0.0
    return min(1.0, math.log1p(member_count) / math.log1p(cap))


def _single_tx_internal_exon_compatibility(model: CodingModel, tx: TranscriptEvidenceCollapsedTranscript) -> Optional[float]:
    internal = _internal_cds_exons(model)
    if not internal:
        return None
    supported = sum(1 for exon in internal if _internal_cds_exon_supported_by_transcript(exon, tx))
    return supported / len(internal)


def _single_tx_terminal_exon_compatibility(model: CodingModel, tx: TranscriptEvidenceCollapsedTranscript) -> float:
    terminal = _terminal_cds_exons(model)
    if not terminal:
        return 0.0
    supported = sum(1 for exon in terminal if _terminal_cds_exon_supported_by_transcript(model, exon, tx))
    return supported / len(terminal)


def _pairwise_transcript_score_components(
    model: CodingModel,
    tx: TranscriptEvidenceCollapsedTranscript,
    metrics: dict[str, object],
) -> dict[str, Optional[float]]:
    """Return scoring-ready pairwise transcript components for one association."""
    model_cds_intron_support = float(metrics.get("model_intron_precision") or 0.0)
    splice_site_f1 = float(metrics.get("splice_site_f1") or 0.0)
    full_intron_chain_f1 = float(metrics.get("intron_f1") or 0.0)
    internal_compat = _single_tx_internal_exon_compatibility(model, tx)
    terminal_compat = _single_tx_terminal_exon_compatibility(model, tx)
    model_cds_exonic_coverage = float(metrics.get("model_cds_overlap_fraction") or 0.0)
    if model.is_multi_exon:
        pairwise = weighted_metric_value([
            ("model_cds_intron_support", 0.30, model_cds_intron_support),
            ("splice_site_f1", 0.20, splice_site_f1),
            ("full_intron_chain_f1", 0.20, full_intron_chain_f1),
            ("internal_cds_exon_compatibility", 0.15, internal_compat),
            ("terminal_cds_exon_compatibility", 0.10, terminal_compat),
            ("model_cds_exonic_coverage", 0.05, model_cds_exonic_coverage),
        ])
        if pairwise is None:
            pairwise = 0.0
    else:
        pairwise = 0.85 * model_cds_exonic_coverage + 0.15 * terminal_compat
    return {
        "pairwise_transcript_score": max(0.0, min(pairwise, 1.0)),
        "model_cds_intron_support": model_cds_intron_support,
        "splice_site_f1": splice_site_f1,
        "full_intron_chain_f1": full_intron_chain_f1,
        "internal_cds_exon_compatibility": internal_compat,
        "terminal_cds_exon_compatibility": terminal_compat,
        "model_cds_exonic_coverage": model_cds_exonic_coverage,
    }


def _best_per_source_pairwise_components(
    scoring_associations: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]],
) -> dict[str, object]:
    """Choose the best scoring-eligible transcript per source and average components."""
    best_by_source: dict[str, tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]] = {}
    for tx, metrics in scoring_associations:
        score = float(metrics.get("pairwise_transcript_score") or 0.0)
        for source in tx.sources:
            prev = best_by_source.get(source)
            if prev is None or score > float(prev[1].get("pairwise_transcript_score") or 0.0):
                best_by_source[source] = (tx, metrics)
    best_items = list(best_by_source.values())
    out: dict[str, object] = {
        "source_count": len(best_by_source),
        "mean_best_source_pairwise_score": _mean([float(m.get("pairwise_transcript_score") or 0.0) for _tx, m in best_items]),
    }
    component_metric_keys = {
        "model_cds_intron_support": "pairwise_model_cds_intron_support",
        "splice_site_f1": "pairwise_splice_site_f1",
        "full_intron_chain_f1": "pairwise_full_intron_chain_f1",
        "internal_cds_exon_compatibility": "pairwise_internal_cds_exon_compatibility",
        "terminal_cds_exon_compatibility": "pairwise_terminal_cds_exon_compatibility",
        "model_cds_exonic_coverage": "pairwise_model_cds_exonic_coverage",
    }
    for output_key, metric_key in component_metric_keys.items():
        if output_key == "internal_cds_exon_compatibility":
            # Internal CDS exon compatibility is not applicable when a model
            # has no internal CDS exons. Exclude unavailable source values from
            # the mean and preserve None when every selected source is N/A.
            values = [
                float(m[metric_key])
                for _tx, m in best_items
                if m.get(metric_key) is not None
            ]
            out[f"mean_best_source_{output_key}"] = _mean(values) if values else None
        else:
            out[f"mean_best_source_{output_key}"] = _mean(
                [float(m.get(metric_key) or 0.0) for _tx, m in best_items]
            )
    return out


def _weighted_mean_scoring_eligible_pairwise_score(
    scoring_associations: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]],
    total_sources: int,
) -> float:
    """Return weighted mean pairwise score across all scoring-eligible records.

    Weighting favours independent source support while allowing capped member
    count to contribute as a secondary confidence signal.
    """
    weighted_sum = 0.0
    weight_sum = 0.0
    for tx, metrics in scoring_associations:
        source_fraction = _safe_fraction(tx.source_count, total_sources)
        member_scaled = _log_member_support_scaled(tx.member_count)
        weight = 0.80 * source_fraction + 0.20 * member_scaled
        if weight <= 0:
            continue
        weighted_sum += float(metrics.get("pairwise_transcript_score") or 0.0) * weight
        weight_sum += weight
    return weighted_sum / weight_sum if weight_sum > 0 else 0.0


def _transcript_chain_source_score(
    *,
    exact_source_count: int,
    stringent_source_count: int,
    contained_source_count: int,
    total_sources: int,
) -> float:
    exact_fraction = _safe_fraction(exact_source_count, total_sources)
    stringent_fraction = _safe_fraction(stringent_source_count, total_sources)
    contained_fraction = _safe_fraction(contained_source_count, total_sources)
    score = (
        exact_fraction
        + 0.75 * max(0.0, stringent_fraction - exact_fraction)
        + 0.35 * max(0.0, contained_fraction - stringent_fraction)
    )
    return max(0.0, min(score, 1.0))


def _transcript_score_status(score: float, *, has_associations: bool, has_scoring: bool) -> str:
    if not has_associations:
        return "no_associated_transcript_evidence"
    if not has_scoring:
        return "transcript_broad_association_only"
    if score >= 0.75:
        return "strong_transcript_support"
    if score >= 0.50:
        return "partial_transcript_support"
    if score > 0:
        return "weak_transcript_support"
    return "transcript_scoring_eligible_no_score"


def _fmt_score_components(parts: list[tuple[str, float, float]]) -> str:
    if not parts:
        return "."
    return ",".join(f"+{name}:{value:.6f}*{weight:.3f}" for name, value, weight in parts)


def transcript_model_support_flags(
    *,
    total_sources: int,
    associations: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]],
    scoring_associations: list[tuple[TranscriptEvidenceCollapsedTranscript, dict[str, object]]],
) -> str:
    flags: list[str] = []
    if total_sources <= 0:
        flags.append("no_transcript_evidence_configured")
        return fmt_flags(flags)
    if not associations:
        flags.append("no_associated_transcript_evidence")
        return fmt_flags(flags)
    if scoring_associations:
        flags.append("has_scoring_eligible_transcript_evidence")
    else:
        flags.append("broad_association_only")
    classes = Counter(str(m.get("association_class", ".")) for _tx, m in associations)
    for cls in sorted(classes):
        if cls and cls != ".":
            flags.append(f"transcript_assoc_{cls}")
    if any(str(m.get("flags", ".")).find("possible_fragment") >= 0 for _tx, m in associations):
        flags.append("possible_fragment_signal")
    if any(tx.strand == "." for tx, _m in associations):
        flags.append("has_unstranded_transcript_association")
    return fmt_flags(flags)


def calculate_transcript_model_support(
    models: list[CodingModel],
    transcript_superloci: list[TranscriptEvidenceSuperlocus],
    transcript_source_labels: list[str],
) -> list[dict[str, object]]:
    """Calculate model-to-transcript support summaries and evidence scores.

    Broad associations are retained for diagnostics, while scoring-eligible
    associations contribute to ``transcript_evidence_score``. Downstream code
    can use that score in model selection ranking, singleton rescue, and
    evidence-aware low-score pruning according to the active configuration.
    """
    superlocus_index = IntervalOverlapIndex(
        transcript_superloci,
        key=lambda sl: sl.chrom,
        start=lambda sl: sl.start0,
        end=lambda sl: sl.end0,
    )
    total_sources = len(transcript_source_labels)
    rows: list[dict[str, object]] = []
    for model in sorted(models, key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.output_transcript_id)):
        associations = associated_transcripts_for_model(model, superlocus_index)
        scoring_associations = [(tx, metrics) for tx, metrics in associations if bool(metrics.get("scoring_eligible"))]
        broad_sources = sorted({s for tx, _m in associations for s in tx.sources})
        scoring_sources = sorted({s for tx, _m in scoring_associations for s in tx.sources})
        class_counts = Counter(str(m.get("association_class", ".")) for _tx, m in associations)
        scoring_class_counts = Counter(str(m.get("association_class", ".")) for _tx, m in scoring_associations)
        evidence_super_ids = {tx.evidence_superlocus_id for tx, _m in associations}
        evidence_sub_ids = {tx.evidence_sublocus_id for tx, _m in associations if tx.evidence_sublocus_id}
        scoring_super_ids = {tx.evidence_superlocus_id for tx, _m in scoring_associations}
        scoring_sub_ids = {tx.evidence_sublocus_id for tx, _m in scoring_associations if tx.evidence_sublocus_id}
        best_tx: TranscriptEvidenceCollapsedTranscript | None = None
        best_metrics: dict[str, object] | None = None
        for tx, metrics in associations:
            if best_metrics is None or float(metrics.get("association_score") or 0.0) > float(best_metrics.get("association_score") or 0.0):
                best_tx = tx
                best_metrics = metrics
        best_scoring_tx: TranscriptEvidenceCollapsedTranscript | None = None
        best_scoring_metrics: dict[str, object] | None = None
        for tx, metrics in scoring_associations:
            if best_scoring_metrics is None or float(metrics.get("association_score") or 0.0) > float(best_scoring_metrics.get("association_score") or 0.0):
                best_scoring_tx = tx
                best_scoring_metrics = metrics

        scoring_member_count = sum(tx.member_count for tx, _m in scoring_associations)
        broad_member_count = sum(tx.member_count for tx, _m in associations)
        union_exonic_bp = interval_list_overlap_bp(list(model.cds_parts), [part for tx, _m in associations for part in tx.exon_parts]) if associations else 0
        union_scoring_exonic_bp = interval_list_overlap_bp(list(model.cds_parts), [part for tx, _m in scoring_associations for part in tx.exon_parts]) if scoring_associations else 0
        broad_model_cds_fraction = union_exonic_bp / model.cds_length if model.cds_length else 0.0
        scoring_model_cds_fraction = union_scoring_exonic_bp / model.cds_length if model.cds_length else 0.0
        model_introns = _model_introns(model)
        model_splice_sites = _model_splice_sites(model)
        tx_introns_union: set[tuple[int, int]] = set()
        tx_splice_union: set[int] = set()
        scoring_tx_introns_union: set[tuple[int, int]] = set()
        scoring_tx_splice_union: set[int] = set()
        for tx, _metrics in associations:
            tx_introns_union.update(tx.introns)
            tx_splice_union.update(tx.splice_sites)
        for tx, _metrics in scoring_associations:
            scoring_tx_introns_union.update(tx.introns)
            scoring_tx_splice_union.update(tx.splice_sites)
        supported_introns = model_introns & tx_introns_union
        supported_splice_sites = model_splice_sites & tx_splice_union
        scoring_supported_introns = model_introns & scoring_tx_introns_union
        scoring_supported_splice_sites = model_splice_sites & scoring_tx_splice_union

        association_txs = [tx for tx, _metrics in associations]
        scoring_association_txs = [tx for tx, _metrics in scoring_associations]
        internal_cds_exons = _internal_cds_exons(model)
        terminal_cds_exons = _terminal_cds_exons(model)
        supported_internal_exons = _supported_model_exons(model, association_txs, terminal=False)
        supported_terminal_exons = _supported_model_exons(model, association_txs, terminal=True)
        scoring_supported_internal_exons = _supported_model_exons(model, scoring_association_txs, terminal=False)
        scoring_supported_terminal_exons = _supported_model_exons(model, scoring_association_txs, terminal=True)

        exact_chain_txs = [tx for tx, metrics in associations if bool(metrics.get("exact_intron_chain_match"))]
        contained_chain_txs = [tx for tx, metrics in associations if bool(metrics.get("model_intron_chain_contained"))]
        stringent_contained_chain_txs = [tx for tx, metrics in associations if bool(metrics.get("model_intron_chain_stringent_contained"))]
        exact_chain_source_count, exact_chain_member_count = _source_member_summary_for_transcripts(exact_chain_txs)
        contained_chain_source_count, contained_chain_member_count = _source_member_summary_for_transcripts(contained_chain_txs)
        stringent_chain_source_count, stringent_chain_member_count = _source_member_summary_for_transcripts(stringent_contained_chain_txs)
        best_intron_chain_f1 = max((float(metrics.get("intron_f1") or 0.0) for _tx, metrics in associations), default=0.0)
        mean_best_source_intron_chain_f1 = _mean_best_source_metric(associations, "intron_f1")

        limited_multiexon_txs = [tx for tx, metrics in associations if bool(metrics.get("single_exon_limited_multiexon_support"))]
        limited_multiexon_source_count, limited_multiexon_member_count = _source_member_summary_for_transcripts(limited_multiexon_txs)
        unstranded_association_txs = [tx for tx, _metrics in associations if tx.strand == "."]
        unstranded_scoring_txs = [tx for tx, _metrics in scoring_associations if tx.strand == "."]
        stranded_scoring_txs = [tx for tx, _metrics in scoring_associations if tx.strand != "."]
        unstranded_broad_source_count, unstranded_broad_member_count = _source_member_summary_for_transcripts(unstranded_association_txs)
        unstranded_scoring_source_count, unstranded_scoring_member_count = _source_member_summary_for_transcripts(unstranded_scoring_txs)
        stranded_scoring_source_count, stranded_scoring_member_count = _source_member_summary_for_transcripts(stranded_scoring_txs)

        best_pairwise_tx: TranscriptEvidenceCollapsedTranscript | None = None
        best_pairwise_metrics: dict[str, object] | None = None
        for tx, metrics in scoring_associations:
            if best_pairwise_metrics is None or float(metrics.get("pairwise_transcript_score") or 0.0) > float(best_pairwise_metrics.get("pairwise_transcript_score") or 0.0):
                best_pairwise_tx = tx
                best_pairwise_metrics = metrics
        cdna_context = cdna_context_for_model(
            model,
            associations,
            best_pairwise_id=best_pairwise_tx.transcript_evidence_id if best_pairwise_tx else ".",
        )
        best_source_components = _best_per_source_pairwise_components(scoring_associations)
        mean_best_source_pairwise_score = float(best_source_components.get("mean_best_source_pairwise_score") or 0.0)
        raw_mean_internal_compatibility = best_source_components.get("mean_best_source_internal_cds_exon_compatibility")
        mean_best_source_internal_compatibility = (
            None
            if not internal_cds_exons
            else float(raw_mean_internal_compatibility) if raw_mean_internal_compatibility is not None else 0.0
        )
        weighted_mean_pairwise_score = _weighted_mean_scoring_eligible_pairwise_score(scoring_associations, total_sources)
        supporting_sources = sorted({
            source
            for tx, metrics in scoring_associations
            if float(metrics.get("pairwise_transcript_score") or 0.0) >= TRANSCRIPT_SUPPORTING_SOURCE_SCORE_THRESHOLD
            for source in tx.sources
        })
        supporting_source_count = len(supporting_sources)
        supporting_source_fraction = _safe_fraction(supporting_source_count, total_sources)
        member_support_scaled = _log_member_support_scaled(scoring_member_count)
        chain_source_score = _transcript_chain_source_score(
            exact_source_count=exact_chain_source_count,
            stringent_source_count=stringent_chain_source_count,
            contained_source_count=contained_chain_source_count,
            total_sources=total_sources,
        )
        support_confidence = 0.80 * supporting_source_fraction + 0.20 * member_support_scaled
        if not scoring_associations:
            transcript_score = 0.0
            score_type = "none"
            score_components = "."
        elif model.is_multi_exon:
            transcript_score = (
                0.40 * mean_best_source_pairwise_score
                + 0.25 * weighted_mean_pairwise_score
                + 0.25 * chain_source_score
                + 0.10 * support_confidence
            )
            score_type = "multi_exon"
            score_components = _fmt_score_components([
                ("transcript_mean_best_source_pairwise_score", mean_best_source_pairwise_score, 0.40),
                ("transcript_weighted_mean_scoring_eligible_pairwise_score", weighted_mean_pairwise_score, 0.25),
                ("transcript_chain_source_score", chain_source_score, 0.25),
                ("transcript_support_confidence", support_confidence, 0.10),
            ])
        else:
            transcript_score = (
                0.45 * mean_best_source_pairwise_score
                + 0.35 * weighted_mean_pairwise_score
                + 0.15 * supporting_source_fraction
                + 0.05 * member_support_scaled
            )
            score_type = "single_exon"
            score_components = _fmt_score_components([
                ("transcript_mean_best_source_pairwise_score", mean_best_source_pairwise_score, 0.45),
                ("transcript_weighted_mean_scoring_eligible_pairwise_score", weighted_mean_pairwise_score, 0.35),
                ("transcript_supporting_source_fraction", supporting_source_fraction, 0.15),
                ("transcript_member_support_scaled", member_support_scaled, 0.05),
            ])
        transcript_score = scale_transcript_evidence_score(transcript_score)

        flags = transcript_model_support_flags(total_sources=total_sources, associations=associations, scoring_associations=scoring_associations)
        status = (
            "transcript_evidence_not_configured" if total_sources <= 0 else
            _transcript_score_status(transcript_score, has_associations=bool(associations), has_scoring=bool(scoring_associations))
        )
        score_flag = {
            "strong_transcript_support": "strong_transcript_support",
            "partial_transcript_support": "partial_transcript_support",
            "weak_transcript_support": "weak_transcript_support",
        }.get(status)
        if score_flag:
            flags = fmt_flags(([] if flags == "." else flags.split(',')) + [score_flag])
        row = {
            "superlocus_id": model.superlocus_id or ".",
            "sublocus_id": model.sublocus_id or ".",
            "source_name": model.source_name,
            "input_gene_id": model.input_gene_id,
            "input_transcript_id": model.input_transcript_id,
            "output_gene_id": model.output_gene_id,
            "output_transcript_id": model.output_transcript_id,
            "chrom": model.chrom,
            "strand": model.strand,
            "cds_start": model.cds_start0 + 1,
            "cds_end": model.cds_end0,
            "cds_length": model.cds_length,
            "cds_exon_count": model.cds_exon_count,
            "transcript_evidence_source_total": total_sources,
            "transcript_evidence_score": transcript_score,
            "transcript_evidence_score_scaled": scale_transcript_evidence_score(transcript_score),
            "transcript_evidence_score_type": score_type,
            "transcript_evidence_score_components": score_components,
            "transcript_evidence_superlocus_count": len(evidence_super_ids),
            "transcript_evidence_sublocus_count": len(evidence_sub_ids),
            "transcript_scoring_evidence_superlocus_count": len(scoring_super_ids),
            "transcript_scoring_evidence_sublocus_count": len(scoring_sub_ids),
            "transcript_broad_association_count": len(associations),
            "transcript_broad_association_source_count": len(broad_sources),
            "transcript_broad_association_sources": ",".join(broad_sources) if broad_sources else ".",
            "transcript_broad_association_member_count": broad_member_count,
            "transcript_scoring_eligible_count": len(scoring_associations),
            "transcript_scoring_eligible_source_count": len(scoring_sources),
            "transcript_scoring_eligible_sources": ",".join(scoring_sources) if scoring_sources else ".",
            "transcript_scoring_eligible_member_count": scoring_member_count,
            "transcript_unstranded_broad_association_count": len(unstranded_association_txs),
            "transcript_unstranded_broad_association_source_count": unstranded_broad_source_count,
            "transcript_unstranded_broad_association_member_count": unstranded_broad_member_count,
            "transcript_unstranded_scoring_eligible_count": len(unstranded_scoring_txs),
            "transcript_unstranded_scoring_eligible_source_count": unstranded_scoring_source_count,
            "transcript_unstranded_scoring_eligible_member_count": unstranded_scoring_member_count,
            "transcript_stranded_scoring_eligible_count": len(stranded_scoring_txs),
            "transcript_stranded_scoring_eligible_source_count": stranded_scoring_source_count,
            "transcript_stranded_scoring_eligible_member_count": stranded_scoring_member_count,
            "transcript_best_pairwise_score": float(best_pairwise_metrics.get("pairwise_transcript_score") or 0.0) if best_pairwise_metrics else 0.0,
            "transcript_best_pairwise_id": best_pairwise_tx.transcript_evidence_id if best_pairwise_tx else ".",
            "transcript_best_pairwise_source_count": best_pairwise_tx.source_count if best_pairwise_tx else 0,
            "transcript_best_pairwise_member_count": best_pairwise_tx.member_count if best_pairwise_tx else 0,
            **cdna_context,
            "transcript_mean_best_source_pairwise_score": mean_best_source_pairwise_score,
            "transcript_weighted_mean_scoring_eligible_pairwise_score": weighted_mean_pairwise_score,
            "transcript_supporting_source_count": supporting_source_count,
            "transcript_supporting_source_fraction": supporting_source_fraction,
            "transcript_supporting_sources": ",".join(supporting_sources) if supporting_sources else ".",
            "transcript_member_support_scaled": member_support_scaled,
            "transcript_chain_source_score": chain_source_score,
            "transcript_support_confidence": support_confidence,
            "transcript_mean_best_source_model_cds_intron_support": float(best_source_components.get("mean_best_source_model_cds_intron_support") or 0.0),
            "transcript_mean_best_source_splice_site_f1": float(best_source_components.get("mean_best_source_splice_site_f1") or 0.0),
            "transcript_mean_best_source_full_intron_chain_f1": float(best_source_components.get("mean_best_source_full_intron_chain_f1") or 0.0),
            "transcript_mean_best_source_internal_cds_exon_compatibility": mean_best_source_internal_compatibility,
            "transcript_mean_best_source_terminal_cds_exon_compatibility": float(best_source_components.get("mean_best_source_terminal_cds_exon_compatibility") or 0.0),
            "transcript_mean_best_source_model_cds_exonic_coverage": float(best_source_components.get("mean_best_source_model_cds_exonic_coverage") or 0.0),
            "transcript_association_class_counts": fmt_counts(dict(class_counts)),
            "transcript_scoring_association_class_counts": fmt_counts(dict(scoring_class_counts)),
            "transcript_shared_intron_association_count": class_counts.get("shared_intron", 0),
            "transcript_shared_splice_site_association_count": class_counts.get("shared_splice_site", 0),
            "transcript_single_exon_support_association_count": class_counts.get("single_exon_transcript_support", 0),
            "transcript_single_exon_limited_multiexon_support_count": class_counts.get("single_exon_limited_multiexon_support", 0),
            "transcript_single_exon_limited_multiexon_support_source_count": limited_multiexon_source_count,
            "transcript_single_exon_limited_multiexon_support_member_count": limited_multiexon_member_count,
            "transcript_single_exon_model_multiexon_overlap_count": class_counts.get("single_exon_model_multiexon_transcript_overlap", 0),
            "transcript_multiexon_model_single_exon_overlap_count": class_counts.get("multiexon_model_single_exon_transcript_overlap", 0),
            "transcript_model_cds_exonic_overlap_count": class_counts.get("model_cds_exonic_overlap", 0),
            "transcript_span_only_association_count": class_counts.get("span_contains_model", 0) + class_counts.get("model_contains_transcript_span", 0) + class_counts.get("span_overlap_only", 0),
            "transcript_model_cds_intron_count": len(model_introns),
            "transcript_supported_cds_intron_count": len(supported_introns),
            "transcript_supported_cds_intron_fraction": len(supported_introns) / len(model_introns) if model_introns else 0.0,
            "transcript_scoring_supported_cds_intron_count": len(scoring_supported_introns),
            "transcript_scoring_supported_cds_intron_fraction": len(scoring_supported_introns) / len(model_introns) if model_introns else 0.0,
            "transcript_model_cds_splice_site_count": len(model_splice_sites),
            "transcript_supported_cds_splice_site_count": len(supported_splice_sites),
            "transcript_supported_cds_splice_site_fraction": len(supported_splice_sites) / len(model_splice_sites) if model_splice_sites else 0.0,
            "transcript_scoring_supported_cds_splice_site_count": len(scoring_supported_splice_sites),
            "transcript_scoring_supported_cds_splice_site_fraction": len(scoring_supported_splice_sites) / len(model_splice_sites) if model_splice_sites else 0.0,
            "transcript_model_internal_cds_exon_count": len(internal_cds_exons),
            "transcript_supported_internal_cds_exon_count": len(supported_internal_exons),
            "transcript_supported_internal_cds_exon_fraction": len(supported_internal_exons) / len(internal_cds_exons) if internal_cds_exons else None,
            "transcript_scoring_supported_internal_cds_exon_count": len(scoring_supported_internal_exons),
            "transcript_scoring_supported_internal_cds_exon_fraction": len(scoring_supported_internal_exons) / len(internal_cds_exons) if internal_cds_exons else None,
            "transcript_model_terminal_cds_exon_count": len(terminal_cds_exons),
            "transcript_supported_terminal_cds_exon_count": len(supported_terminal_exons),
            "transcript_supported_terminal_cds_exon_fraction": len(supported_terminal_exons) / len(terminal_cds_exons) if terminal_cds_exons else 0.0,
            "transcript_scoring_supported_terminal_cds_exon_count": len(scoring_supported_terminal_exons),
            "transcript_scoring_supported_terminal_cds_exon_fraction": len(scoring_supported_terminal_exons) / len(terminal_cds_exons) if terminal_cds_exons else 0.0,
            "transcript_exact_intron_chain_match_count": len(exact_chain_txs),
            "transcript_exact_intron_chain_match_source_count": exact_chain_source_count,
            "transcript_exact_intron_chain_match_member_count": exact_chain_member_count,
            "transcript_model_intron_chain_contained_count": len(contained_chain_txs),
            "transcript_model_intron_chain_contained_source_count": contained_chain_source_count,
            "transcript_model_intron_chain_contained_member_count": contained_chain_member_count,
            "transcript_model_intron_chain_stringent_contained_count": len(stringent_contained_chain_txs),
            "transcript_model_intron_chain_stringent_contained_source_count": stringent_chain_source_count,
            "transcript_model_intron_chain_stringent_contained_member_count": stringent_chain_member_count,
            "transcript_best_intron_chain_f1": best_intron_chain_f1,
            "transcript_mean_best_source_intron_chain_f1": mean_best_source_intron_chain_f1,
            "transcript_supported_cds_base_bp": union_exonic_bp,
            "transcript_supported_cds_base_fraction": broad_model_cds_fraction,
            "transcript_scoring_supported_cds_base_bp": union_scoring_exonic_bp,
            "transcript_scoring_supported_cds_base_fraction": scoring_model_cds_fraction,
            "transcript_best_association_id": best_tx.transcript_evidence_id if best_tx else ".",
            "transcript_best_association_class": str(best_metrics.get("association_class")) if best_metrics else ".",
            "transcript_best_association_score": float(best_metrics.get("association_score") or 0.0) if best_metrics else 0.0,
            "transcript_best_association_source_count": best_tx.source_count if best_tx else 0,
            "transcript_best_association_member_count": best_tx.member_count if best_tx else 0,
            "transcript_best_association_model_cds_overlap_fraction": float(best_metrics.get("model_cds_overlap_fraction") or 0.0) if best_metrics else 0.0,
            "transcript_best_association_shared_intron_count": int(best_metrics.get("shared_intron_count") or 0) if best_metrics else 0,
            "transcript_best_association_shared_splice_site_count": int(best_metrics.get("shared_splice_site_count") or 0) if best_metrics else 0,
            "transcript_best_scoring_association_id": best_scoring_tx.transcript_evidence_id if best_scoring_tx else ".",
            "transcript_best_scoring_association_class": str(best_scoring_metrics.get("association_class")) if best_scoring_metrics else ".",
            "transcript_best_scoring_association_score": float(best_scoring_metrics.get("association_score") or 0.0) if best_scoring_metrics else 0.0,
            "transcript_best_scoring_association_source_count": best_scoring_tx.source_count if best_scoring_tx else 0,
            "transcript_best_scoring_association_member_count": best_scoring_tx.member_count if best_scoring_tx else 0,
            "transcript_associated_evidence_ids": _fmt_id_list({tx.transcript_evidence_id for tx, _m in associations}),
            "transcript_scoring_evidence_ids": _fmt_id_list({tx.transcript_evidence_id for tx, _m in scoring_associations}),
            "transcript_evidence_status": status,
            "transcript_evidence_flags": flags,
        }
        rows.append(row)
    return rows


TRANSCRIPT_MODEL_SUPPORT_COLUMNS = [
    "superlocus_id", "sublocus_id", "source_name", "input_gene_id", "input_transcript_id", "output_gene_id", "output_transcript_id",
    "chrom", "strand", "cds_start", "cds_end", "cds_length", "cds_exon_count",
    "transcript_evidence_source_total", "transcript_evidence_score", "transcript_evidence_score_scaled", "transcript_evidence_score_type", "transcript_evidence_score_components",
    "transcript_evidence_superlocus_count", "transcript_evidence_sublocus_count",
    "transcript_scoring_evidence_superlocus_count", "transcript_scoring_evidence_sublocus_count",
    "transcript_broad_association_count", "transcript_broad_association_source_count", "transcript_broad_association_sources", "transcript_broad_association_member_count",
    "transcript_scoring_eligible_count", "transcript_scoring_eligible_source_count", "transcript_scoring_eligible_sources", "transcript_scoring_eligible_member_count",
    "transcript_unstranded_broad_association_count", "transcript_unstranded_broad_association_source_count", "transcript_unstranded_broad_association_member_count",
    "transcript_unstranded_scoring_eligible_count", "transcript_unstranded_scoring_eligible_source_count", "transcript_unstranded_scoring_eligible_member_count",
    "transcript_stranded_scoring_eligible_count", "transcript_stranded_scoring_eligible_source_count", "transcript_stranded_scoring_eligible_member_count",
    "transcript_best_pairwise_score", "transcript_best_pairwise_id", "transcript_best_pairwise_source_count", "transcript_best_pairwise_member_count",
    "transcript_cdna_context_status", "transcript_cdna_context_best_evidence_id", "transcript_cdna_context_best_association_class", "transcript_cdna_context_best_matches_best_pairwise",
    "transcript_cdna_context_best_source_count", "transcript_cdna_context_best_member_count", "transcript_cdna_context_model_cds_length_bp",
    "transcript_cdna_context_cdna_length_min_bp", "transcript_cdna_context_cdna_length_max_bp",
    "transcript_cdna_context_cds_to_min_cdna_fraction", "transcript_cdna_context_cds_to_max_cdna_fraction",
    "transcript_cdna_context_extra_5prime_exon_count", "transcript_cdna_context_extra_3prime_exon_count",
    "transcript_cdna_context_extra_5prime_exonic_bp", "transcript_cdna_context_extra_3prime_exonic_bp",
    "transcript_cdna_context_extra_5prime_intron_count", "transcript_cdna_context_extra_3prime_intron_count", "transcript_cdna_context_flags",
    "transcript_cdna_context_evaluable", "transcript_cdna_context_concern_score", "transcript_cdna_context_concern_class",
    "transcript_cdna_context_pruning_concern", "transcript_cdna_context_pruning_concern_class", "transcript_cdna_context_pruning_supporting_flags",
    "transcript_cdna_context_neutral_fraction", "transcript_cdna_context_severe_fraction",
    "transcript_mean_best_source_pairwise_score", "transcript_weighted_mean_scoring_eligible_pairwise_score",
    "transcript_supporting_source_count", "transcript_supporting_source_fraction", "transcript_supporting_sources", "transcript_member_support_scaled",
    "transcript_chain_source_score", "transcript_support_confidence",
    "transcript_mean_best_source_model_cds_intron_support", "transcript_mean_best_source_splice_site_f1", "transcript_mean_best_source_full_intron_chain_f1",
    "transcript_mean_best_source_internal_cds_exon_compatibility", "transcript_mean_best_source_terminal_cds_exon_compatibility", "transcript_mean_best_source_model_cds_exonic_coverage",
    "transcript_association_class_counts", "transcript_scoring_association_class_counts",
    "transcript_shared_intron_association_count", "transcript_shared_splice_site_association_count", "transcript_single_exon_support_association_count",
    "transcript_single_exon_limited_multiexon_support_count", "transcript_single_exon_limited_multiexon_support_source_count", "transcript_single_exon_limited_multiexon_support_member_count",
    "transcript_single_exon_model_multiexon_overlap_count", "transcript_multiexon_model_single_exon_overlap_count", "transcript_model_cds_exonic_overlap_count", "transcript_span_only_association_count",
    "transcript_model_cds_intron_count", "transcript_supported_cds_intron_count", "transcript_supported_cds_intron_fraction", "transcript_scoring_supported_cds_intron_count", "transcript_scoring_supported_cds_intron_fraction",
    "transcript_model_cds_splice_site_count", "transcript_supported_cds_splice_site_count", "transcript_supported_cds_splice_site_fraction", "transcript_scoring_supported_cds_splice_site_count", "transcript_scoring_supported_cds_splice_site_fraction",
    "transcript_model_internal_cds_exon_count", "transcript_supported_internal_cds_exon_count", "transcript_supported_internal_cds_exon_fraction", "transcript_scoring_supported_internal_cds_exon_count", "transcript_scoring_supported_internal_cds_exon_fraction",
    "transcript_model_terminal_cds_exon_count", "transcript_supported_terminal_cds_exon_count", "transcript_supported_terminal_cds_exon_fraction", "transcript_scoring_supported_terminal_cds_exon_count", "transcript_scoring_supported_terminal_cds_exon_fraction",
    "transcript_exact_intron_chain_match_count", "transcript_exact_intron_chain_match_source_count", "transcript_exact_intron_chain_match_member_count",
    "transcript_model_intron_chain_contained_count", "transcript_model_intron_chain_contained_source_count", "transcript_model_intron_chain_contained_member_count",
    "transcript_model_intron_chain_stringent_contained_count", "transcript_model_intron_chain_stringent_contained_source_count", "transcript_model_intron_chain_stringent_contained_member_count",
    "transcript_best_intron_chain_f1", "transcript_mean_best_source_intron_chain_f1",
    "transcript_supported_cds_base_bp", "transcript_supported_cds_base_fraction", "transcript_scoring_supported_cds_base_bp", "transcript_scoring_supported_cds_base_fraction",
    "transcript_best_association_id", "transcript_best_association_class", "transcript_best_association_score", "transcript_best_association_source_count", "transcript_best_association_member_count", "transcript_best_association_model_cds_overlap_fraction", "transcript_best_association_shared_intron_count", "transcript_best_association_shared_splice_site_count",
    "transcript_best_scoring_association_id", "transcript_best_scoring_association_class", "transcript_best_scoring_association_score", "transcript_best_scoring_association_source_count", "transcript_best_scoring_association_member_count",
    "transcript_associated_evidence_ids", "transcript_scoring_evidence_ids", "transcript_evidence_status", "transcript_evidence_flags",
]


def transcript_model_support_summary_counts(rows: list[dict[str, object]]) -> dict[str, object]:
    status_counts = Counter(str(row.get("transcript_evidence_status", ".")) for row in rows)
    flag_counts: Counter[str] = Counter()
    class_counts: Counter[str] = Counter()
    scoring_class_counts: Counter[str] = Counter()
    for row in rows:
        flags = str(row.get("transcript_evidence_flags", "."))
        if flags and flags != ".":
            for flag in flags.split(','):
                if flag:
                    flag_counts[flag] += 1
        for item in str(row.get("transcript_association_class_counts", ".")).split(','):
            if ':' in item:
                key, value = item.rsplit(':', 1)
                try:
                    class_counts[key] += int(value)
                except ValueError:
                    pass
        for item in str(row.get("transcript_scoring_association_class_counts", ".")).split(','):
            if ':' in item:
                key, value = item.rsplit(':', 1)
                try:
                    scoring_class_counts[key] += int(value)
                except ValueError:
                    pass
    return {
        "model_support_row_count": len(rows),
        "models_with_broad_transcript_association": sum(1 for row in rows if int(row.get("transcript_broad_association_count", 0) or 0) > 0),
        "models_with_scoring_eligible_transcript_evidence": sum(1 for row in rows if int(row.get("transcript_scoring_eligible_count", 0) or 0) > 0),
        "models_with_shared_intron_transcript_evidence": sum(1 for row in rows if int(row.get("transcript_shared_intron_association_count", 0) or 0) > 0),
        "models_with_shared_splice_site_transcript_evidence": sum(1 for row in rows if int(row.get("transcript_shared_splice_site_association_count", 0) or 0) > 0),
        "models_with_broad_only_transcript_association": sum(1 for row in rows if str(row.get("transcript_evidence_status", ".")) == "transcript_broad_association_only"),
        "models_with_transcript_evidence_score_gt_0": sum(1 for row in rows if float(row.get("transcript_evidence_score") or 0.0) > 0),
        "models_with_partial_or_strong_transcript_support": sum(
            1 for row in rows
            if str(row.get("transcript_evidence_status", ".")) in {"partial_transcript_support", "strong_transcript_support"}
        ),
        "models_with_strong_transcript_support": sum(1 for row in rows if str(row.get("transcript_evidence_status", ".")) == "strong_transcript_support"),
        "models_with_partial_transcript_support": sum(1 for row in rows if str(row.get("transcript_evidence_status", ".")) == "partial_transcript_support"),
        "models_with_transcript_evidence_score_scaled_ge_0.75": sum(1 for row in rows if float(row.get("transcript_evidence_score_scaled") or 0.0) >= 0.75),
        "models_with_unstranded_broad_transcript_association": sum(1 for row in rows if int(row.get("transcript_unstranded_broad_association_count", 0) or 0) > 0),
        "models_with_unstranded_scoring_eligible_transcript_evidence": sum(1 for row in rows if int(row.get("transcript_unstranded_scoring_eligible_count", 0) or 0) > 0),
        "models_with_stranded_scoring_eligible_transcript_evidence": sum(1 for row in rows if int(row.get("transcript_stranded_scoring_eligible_count", 0) or 0) > 0),
        "models_with_only_unstranded_scoring_eligible_transcript_evidence": sum(
            1 for row in rows
            if int(row.get("transcript_unstranded_scoring_eligible_count", 0) or 0) > 0
            and int(row.get("transcript_stranded_scoring_eligible_count", 0) or 0) == 0
        ),
        "models_with_cdna_context_calculated": sum(1 for row in rows if str(row.get("transcript_cdna_context_status", ".")) == "calculated"),
        "models_with_cdna_context_low_cds_fraction": sum(
            1 for row in rows
            if "low_cds_to_cdna_fraction" in str(row.get("transcript_cdna_context_flags", "."))
        ),
        "models_with_cdna_context_very_low_cds_fraction": sum(
            1 for row in rows
            if "very_low_cds_to_cdna_fraction" in str(row.get("transcript_cdna_context_flags", "."))
        ),
        "models_with_cdna_context_not_best_pairwise": sum(
            1 for row in rows
            if str(row.get("transcript_cdna_context_status", ".")) == "calculated"
            and str(row.get("transcript_cdna_context_best_matches_best_pairwise", "False")) not in {"True", "true", "1"}
        ),
        "models_with_cdna_context_pruning_concern": sum(
            1 for row in rows if str(row.get("transcript_cdna_context_pruning_concern", "false")) == "true"
        ),
        "models_with_cdna_context_severe_concern": sum(
            1 for row in rows if str(row.get("transcript_cdna_context_pruning_concern_class", ".")) == "severe_low_fraction"
        ),
        "models_with_cdna_context_moderate_structural_concern": sum(
            1 for row in rows if str(row.get("transcript_cdna_context_pruning_concern_class", ".")) == "moderate_low_fraction_with_structural_context"
        ),
        "unstranded_broad_association_total": sum(int(row.get("transcript_unstranded_broad_association_count", 0) or 0) for row in rows),
        "unstranded_scoring_eligible_total": sum(int(row.get("transcript_unstranded_scoring_eligible_count", 0) or 0) for row in rows),
        "stranded_scoring_eligible_total": sum(int(row.get("transcript_stranded_scoring_eligible_count", 0) or 0) for row in rows),
        "model_status_counts": fmt_counts(dict(status_counts)),
        "model_flag_counts": fmt_counts(dict(flag_counts)),
        "association_class_counts": fmt_counts(dict(class_counts)),
        "scoring_association_class_counts": fmt_counts(dict(scoring_class_counts)),
    }


def transcript_model_support_log_lines(summary: dict[str, object]) -> list[str]:
    return [
        (
            "Transcript evidence model associations: "
            f"models={int(summary['model_support_row_count']):,}; "
            f"broad_association={int(summary['models_with_broad_transcript_association']):,}; "
            f"scoring_eligible={int(summary['models_with_scoring_eligible_transcript_evidence']):,}; "
            f"shared_intron={int(summary['models_with_shared_intron_transcript_evidence']):,}; "
            f"shared_splice_site={int(summary['models_with_shared_splice_site_transcript_evidence']):,}; "
            f"broad_only={int(summary['models_with_broad_only_transcript_association']):,}; "
            f"partial_or_strong={int(summary.get('models_with_partial_or_strong_transcript_support', 0)):,}; "
            f"score_gt_0={int(summary.get('models_with_transcript_evidence_score_gt_0', 0)):,}; "
            f"score_ge_0.75={int(summary.get('models_with_transcript_evidence_score_scaled_ge_0.75', 0)):,}"
        ),
        (
            "Transcript evidence association classes: "
            f"all={summary['association_class_counts']}; scoring={summary['scoring_association_class_counts']}"
        ),
        (
            "Transcript evidence unstranded model associations: "
            f"models_broad={int(summary.get('models_with_unstranded_broad_transcript_association', 0)):,}; "
            f"models_scoring={int(summary.get('models_with_unstranded_scoring_eligible_transcript_evidence', 0)):,}; "
            f"models_only_unstranded_scoring={int(summary.get('models_with_only_unstranded_scoring_eligible_transcript_evidence', 0)):,}; "
            f"unstranded_broad_records={int(summary.get('unstranded_broad_association_total', 0)):,}; "
            f"unstranded_scoring_records={int(summary.get('unstranded_scoring_eligible_total', 0)):,}; "
            f"stranded_scoring_records={int(summary.get('stranded_scoring_eligible_total', 0)):,}"
        ),
        (
            "Transcript cDNA-context diagnostics: "
            f"calculated={int(summary.get('models_with_cdna_context_calculated', 0)):,}; "
            f"low_cds_fraction={int(summary.get('models_with_cdna_context_low_cds_fraction', 0)):,}; "
            f"very_low_cds_fraction={int(summary.get('models_with_cdna_context_very_low_cds_fraction', 0)):,}; "
            f"context_best_not_pairwise_best={int(summary.get('models_with_cdna_context_not_best_pairwise', 0)):,}"
        ),
    ]


def write_transcript_model_support_tsv(rows: list[dict[str, object]], path: Path) -> None:
    write_dict_tsv(rows, path, TRANSCRIPT_MODEL_SUPPORT_COLUMNS)
