"""Protein evidence parsing, exact-chain collapsing, loci, and model-support metrics.

Protein evidence metrics are written to audit TSVs and may optionally contribute
to final model selection scoring when a positive effective protein-evidence
selection weight is active.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .gffio import read_gff_records, sanitize_id_part
from .intervals import (
    IntervalOverlapIndex,
    interval_list_overlap_bp,
    interval_list_overlap_bp_sorted,
    interval_prf_sorted,
)
from .loci import reciprocal_overlap
from .models import CodingModel, UnionFind
from .scoring import weighted_metric_value
from .tsvio import format_tsv_value, write_dict_tsv


# Diagnostic flag thresholds for protein evidence audit tables. These are
# intentionally conservative and currently affect reporting only, not selection.
PROTEIN_EVIDENCE_SHORT_CDS_BP = 90
PROTEIN_EVIDENCE_LARGE_SPAN_BP = 100_000
PROTEIN_EVIDENCE_MANY_COLLAPSED_ALIGNMENTS = 50
PROTEIN_EVIDENCE_MANY_SOURCES = 5
PROTEIN_EVIDENCE_HIGH_MEMBER_COUNT = 100
PROTEIN_EVIDENCE_STRONG_SCORE = 0.75


def fmt_flags(flags: list[str] | set[str]) -> str:
    return ",".join(sorted(set(flags))) if flags else "."


def scale_protein_evidence_score(score: float) -> float:
    """Return a review-friendly 0-1 scaled protein evidence score.

    The current protein_evidence_score is already a weighted mean of 0-1
    components. This helper clamps it for robust downstream reporting and keeps
    a stable column available if future protein scores use a different raw scale.
    """
    if score <= 0:
        return 0.0
    return min(float(score), 1.0)


def weighted_metric_value_with_components(components: list[tuple[str, float, Optional[float]]]) -> tuple[float, str]:
    """Return weighted mean and transparent component text for protein evidence scoring.

    The protein evidence score is currently a diagnostic weighted average whose
    configured components sum to 1.0 for both single- and multi-exon models.
    Component text intentionally mirrors model_selection_score_components
    so benchmark/audit outputs show the exact values and weights used.
    """
    available = [(name, weight, value) for name, weight, value in components if value is not None]
    if not available:
        return 0.0, "."
    total_weight = sum(weight for _name, weight, _value in available)
    if total_weight <= 0:
        return 0.0, "."
    score = sum(weight * float(value) for _name, weight, value in available) / total_weight
    component_text = ",".join(f"+{name}:{float(value):.6f}*{weight:.3f}" for name, weight, value in available)
    if abs(total_weight - 1.0) > 1e-12:
        component_text += f",/available_weight:{total_weight:.3f}"
    return score, component_text


@dataclass(slots=True)
class ProteinEvidenceRawAlignment:
    raw_id: str
    source_label: str
    source_index: int
    input_file: Path
    original_alignment_id: str
    chrom: str
    strand: str
    cds_parts: tuple[tuple[int, int], ...]
    line_count: int
    evidence_superlocus_id: Optional[str] = None
    collapsed_id: Optional[str] = None

    @property
    def sorted_cds_parts(self) -> tuple[tuple[int, int], ...]:
        return self.cds_parts

    @property
    def cds_start0(self) -> int:
        return self.sorted_cds_parts[0][0]

    @property
    def cds_end0(self) -> int:
        return self.sorted_cds_parts[-1][1]

    @property
    def cds_length(self) -> int:
        return sum(e - s for s, e in self.sorted_cds_parts)

    @property
    def cds_exon_count(self) -> int:
        return len(self.cds_parts)

    @property
    def cds_introns(self) -> tuple[tuple[int, int], ...]:
        return tuple((a[1], b[0]) for a, b in zip(self.sorted_cds_parts, self.sorted_cds_parts[1:]))

    @property
    def cds_splice_sites(self) -> frozenset[int]:
        sites: set[int] = set()
        for s, e in self.cds_introns:
            sites.add(s); sites.add(e)
        return frozenset(sites)

    @property
    def is_single_exon(self) -> bool:
        return len(self.cds_parts) == 1

    @property
    def is_multi_exon(self) -> bool:
        return len(self.cds_parts) > 1


@dataclass(slots=True)
class ProteinEvidenceCollapsedAlignment:
    collapsed_id: str
    chrom: str
    strand: str
    cds_parts: tuple[tuple[int, int], ...]
    evidence_superlocus_id: str
    source_member_counts: dict[str, int] = field(default_factory=dict)
    member_count: int = 0
    representative_raw_id: str = "."
    evidence_sublocus_id: Optional[str] = None
    flags: str = "."
    _cds_length: int = field(init=False, repr=False)
    _cds_introns: tuple[tuple[int, int], ...] = field(init=False, repr=False)
    _cds_splice_sites: frozenset[int] = field(init=False, repr=False)
    _sources: frozenset[str] = field(init=False, repr=False)

    def __post_init__(self) -> None:
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
        self._sources = frozenset(
            source for source, count in self.source_member_counts.items() if count > 0
        )

    @property
    def sorted_cds_parts(self) -> tuple[tuple[int, int], ...]:
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
    def is_single_exon(self) -> bool:
        return len(self.cds_parts) == 1

    @property
    def is_multi_exon(self) -> bool:
        return len(self.cds_parts) > 1

    @property
    def source_count(self) -> int:
        return sum(1 for count in self.source_member_counts.values() if count > 0)

    @property
    def sources(self) -> frozenset[str]:
        return self._sources


@dataclass(slots=True)
class ProteinEvidenceSuperlocus:
    evidence_superlocus_id: str
    chrom: str
    strand: str
    start0: int
    end0: int
    raw_alignments: list[ProteinEvidenceRawAlignment] = field(default_factory=list)
    collapsed_alignments: list[ProteinEvidenceCollapsedAlignment] = field(default_factory=list)
    raw_alignment_count: int = 0
    raw_source_member_counts: dict[str, int] = field(default_factory=dict)


@dataclass(slots=True)
class ProteinEvidenceSublocus:
    evidence_sublocus_id: str
    evidence_superlocus_id: str
    chrom: str
    strand: str
    start0: int
    end0: int
    collapsed_alignments: list[ProteinEvidenceCollapsedAlignment] = field(default_factory=list)
    sublocus_type: str = "."
    flags: str = "."


def split_path_list(values: Optional[list[object]]) -> list[Path]:
    paths: list[Path] = []
    if not values:
        return paths
    for item in values:
        for part in str(item).split(','):
            part = part.strip()
            if part:
                paths.append(Path(part))
    return paths


def split_label_list(values: Optional[list[str]]) -> list[str]:
    labels: list[str] = []
    if not values:
        return labels
    for item in values:
        for part in str(item).split(','):
            part = part.strip()
            if part:
                labels.append(sanitize_id_part(part))
    return labels


def default_protein_source_labels(n: int) -> list[str]:
    return [f"protein_{i:03d}" for i in range(1, n + 1)]


def protein_source_labels(paths: list[Path], supplied: Optional[list[str]]) -> list[str]:
    labels = split_label_list(supplied)
    if labels:
        if len(labels) != len(paths):
            raise SystemExit("--evidence-protein-source-name must provide one label per --evidence-proteins file")
    else:
        labels = default_protein_source_labels(len(paths))
    seen: Counter[str] = Counter()
    out: list[str] = []
    for label in labels:
        base = sanitize_id_part(label)
        seen[base] += 1
        out.append(base if seen[base] == 1 else f"{base}_{seen[base]}")
    return out


def parse_protein_evidence_file(path: Path, source_label: str, source_index: int, *, strict: bool = False) -> list[ProteinEvidenceRawAlignment]:
    """Parse protein evidence GFF/GTF CDS features into raw alignments.

    v0.2.1-dev expects CDS features. Features are grouped by Parent,
    transcript_id, protein_id, Target, or ID, in that order.
    """
    records = read_gff_records(path, source_label, strict=strict)
    cds_by_parent: dict[str, list] = defaultdict(list)
    for rec in records:
        if rec.feature_type != "CDS":
            continue
        parents = rec.parents
        if not parents:
            for key in ("transcript_id", "protein_id", "Target", "ID"):
                value = rec.attrs.get(key)
                if value:
                    parents = [value.split()[0]]
                    break
        if not parents:
            parents = [f"unparented_CDS_line_{rec.line_number}"]
        for parent in parents:
            cds_by_parent[parent].append(rec)

    alignments: list[ProteinEvidenceRawAlignment] = []
    for ordinal, (original_id, cds_records) in enumerate(sorted(cds_by_parent.items()), start=1):
        chroms = sorted({r.seqid for r in cds_records})
        strands = sorted({r.strand for r in cds_records})
        valid_strands = {s for s in strands if s in {"+", "-"}}
        if len(chroms) != 1 or len(valid_strands) != 1 or len(strands) != 1:
            # Diagnostic v0.2.1-dev skips malformed evidence alignments rather than
            # failing the whole gene-model run.
            if strict:
                raise ValueError(f"{path}: protein evidence alignment {original_id} has mixed chromosomes/strands")
            continue
        cds_parts = tuple(sorted((r.start0, r.end0) for r in cds_records if r.end0 > r.start0))
        if not cds_parts:
            continue
        raw_id = f"{source_label}|{sanitize_id_part(original_id)}|raw{ordinal:07d}"
        alignments.append(
            ProteinEvidenceRawAlignment(
                raw_id=raw_id,
                source_label=source_label,
                source_index=source_index,
                input_file=path,
                original_alignment_id=original_id,
                chrom=chroms[0],
                strand=next(iter(valid_strands)),
                cds_parts=cds_parts,
                line_count=len(cds_records),
            )
        )
    return alignments


def parse_protein_evidence(paths: list[Path], labels: list[str], *, strict: bool = False) -> list[ProteinEvidenceRawAlignment]:
    raw: list[ProteinEvidenceRawAlignment] = []
    for idx, (path, label) in enumerate(zip(paths, labels), start=1):
        raw.extend(parse_protein_evidence_file(path, label, idx, strict=strict))
    return raw


def make_protein_superlocus_id(chrom: str, start0: int, end0: int, strand: str, index: int) -> str:
    strand_label = "plus" if strand == "+" else "minus"
    return f"protein_superlocus_{index:07d}_{sanitize_id_part(chrom)}_{start0 + 1}_{end0}_{strand_label}"


def build_protein_superloci(raw_alignments: list[ProteinEvidenceRawAlignment]) -> list[ProteinEvidenceSuperlocus]:
    superloci: list[ProteinEvidenceSuperlocus] = []
    by_key: dict[tuple[str, str], list[ProteinEvidenceRawAlignment]] = defaultdict(list)
    for aln in raw_alignments:
        by_key[(aln.chrom, aln.strand)].append(aln)
    next_index = 1
    for (chrom, strand), items in sorted(by_key.items()):
        items = sorted(items, key=lambda a: (a.cds_start0, a.cds_end0, a.source_label, a.raw_id))
        cur: list[ProteinEvidenceRawAlignment] = []
        cur_start: Optional[int] = None
        cur_end: Optional[int] = None
        for aln in items:
            if cur_start is None:
                cur_start = aln.cds_start0; cur_end = aln.cds_end0; cur = [aln]
                continue
            assert cur_end is not None
            if aln.cds_start0 < cur_end:
                cur_end = max(cur_end, aln.cds_end0); cur.append(aln)
            else:
                sl_id = make_protein_superlocus_id(chrom, cur_start, cur_end, strand, next_index); next_index += 1
                sl = ProteinEvidenceSuperlocus(sl_id, chrom, strand, cur_start, cur_end, cur)
                for member in cur:
                    member.evidence_superlocus_id = sl_id
                superloci.append(sl)
                cur_start = aln.cds_start0; cur_end = aln.cds_end0; cur = [aln]
        if cur_start is not None and cur_end is not None:
            sl_id = make_protein_superlocus_id(chrom, cur_start, cur_end, strand, next_index); next_index += 1
            sl = ProteinEvidenceSuperlocus(sl_id, chrom, strand, cur_start, cur_end, cur)
            for member in cur:
                member.evidence_superlocus_id = sl_id
            superloci.append(sl)
    return sorted(superloci, key=lambda sl: (sl.chrom, sl.start0, sl.end0, sl.strand, sl.evidence_superlocus_id))


def chain_key(aln: ProteinEvidenceRawAlignment | ProteinEvidenceCollapsedAlignment | CodingModel) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(aln.cds_parts))


def collapse_protein_superlocus(sl: ProteinEvidenceSuperlocus) -> list[ProteinEvidenceCollapsedAlignment]:
    if sl.raw_alignment_count <= 0:
        sl.raw_alignment_count = len(sl.raw_alignments)
    if not sl.raw_source_member_counts:
        sl.raw_source_member_counts = dict(sorted(Counter(raw.source_label for raw in sl.raw_alignments).items()))
    groups: dict[tuple[tuple[int, int], ...], list[ProteinEvidenceRawAlignment]] = defaultdict(list)
    for raw in sl.raw_alignments:
        groups[chain_key(raw)].append(raw)
    collapsed: list[ProteinEvidenceCollapsedAlignment] = []
    for idx, (parts, members) in enumerate(
        sorted(groups.items(), key=lambda item: (item[0][0][0], item[0][-1][1], len(item[0]), item[0])), start=1
    ):
        cid = f"{sl.evidence_superlocus_id}.collapsed{idx:05d}"
        counts = Counter(m.source_label for m in members)
        obj = ProteinEvidenceCollapsedAlignment(
            collapsed_id=cid,
            chrom=sl.chrom,
            strand=sl.strand,
            cds_parts=parts,
            evidence_superlocus_id=sl.evidence_superlocus_id,
            source_member_counts=dict(sorted(counts.items())),
            member_count=len(members),
            representative_raw_id=sorted(m.raw_id for m in members)[0],
        )
        obj.flags = collapsed_alignment_flags(obj)
        for raw in members:
            raw.collapsed_id = cid
        collapsed.append(obj)
    sl.collapsed_alignments = collapsed
    return collapsed


def collapse_protein_alignments(superloci: list[ProteinEvidenceSuperlocus]) -> list[ProteinEvidenceCollapsedAlignment]:
    out: list[ProteinEvidenceCollapsedAlignment] = []
    for sl in superloci:
        out.extend(collapse_protein_superlocus(sl))
    return out


def release_protein_raw_alignments(
    raw_alignments: list[ProteinEvidenceRawAlignment],
    superloci: list[ProteinEvidenceSuperlocus],
) -> int:
    """Drop raw protein evidence after collapse when no later raw output needs it.

    Counts required by summary/superlocus outputs are retained on each
    superlocus, while the much larger raw objects and duplicate member lists are
    released before transcript evidence and model scoring are processed.
    """
    released = len(raw_alignments)
    for sl in superloci:
        if sl.raw_alignment_count <= 0:
            sl.raw_alignment_count = len(sl.raw_alignments)
        if not sl.raw_source_member_counts:
            sl.raw_source_member_counts = dict(sorted(Counter(raw.source_label for raw in sl.raw_alignments).items()))
        sl.raw_alignments.clear()
    raw_alignments.clear()
    return released


def collapsed_alignment_flags(aln: ProteinEvidenceCollapsedAlignment) -> str:
    flags: list[str] = []
    flags.append("single_exon" if aln.is_single_exon else "multi_exon")
    if aln.source_count <= 1:
        flags.append("single_source_only")
    else:
        flags.append("multi_source_exact_chain")
    if aln.cds_length < PROTEIN_EVIDENCE_SHORT_CDS_BP:
        flags.append("short_cds")
    if aln.cds_length % 3 != 0:
        flags.append("cds_length_not_multiple_of_3")
    return fmt_flags(flags)


def protein_sublocus_flags(items: list[ProteinEvidenceCollapsedAlignment], sublocus_type: str, start0: int, end0: int) -> str:
    source_count = len({source for aln in items for source in aln.sources})
    member_count = sum(aln.member_count for aln in items)
    collapsed_count = len(items)
    span_bp = end0 - start0
    flags: list[str] = []
    if source_count <= 1:
        flags.append("single_source_only")
    if sublocus_type == "mixed":
        flags.append("mixed_single_multi")
    if span_bp > PROTEIN_EVIDENCE_LARGE_SPAN_BP:
        flags.append("large_span")
    if collapsed_count >= PROTEIN_EVIDENCE_MANY_COLLAPSED_ALIGNMENTS:
        flags.append("many_collapsed_alignments")
    if source_count >= PROTEIN_EVIDENCE_MANY_SOURCES:
        flags.append("many_sources")
    if member_count >= PROTEIN_EVIDENCE_HIGH_MEMBER_COUNT:
        flags.append("high_member_count")
    return fmt_flags(flags)


def protein_model_support_flags(
    *,
    protein_source_labels: list[str],
    alns: list[ProteinEvidenceCollapsedAlignment],
    exact: list[ProteinEvidenceCollapsedAlignment],
    supporting_sources: list[str],
    evidence_score: float,
    support_score_threshold: float,
    feature_counts: dict[str, tuple[int, int, Optional[float]]],
    supported_cds_fraction: float,
) -> str:
    flags: list[str] = []
    if not protein_source_labels:
        flags.append("no_protein_evidence_configured")
        return fmt_flags(flags)
    if not alns:
        flags.append("no_overlapping_protein_evidence")
        return fmt_flags(flags)
    if exact:
        flags.append("exact_chain_supported")
    if evidence_score >= PROTEIN_EVIDENCE_STRONG_SCORE:
        flags.append("strong_protein_support")
    elif evidence_score >= support_score_threshold:
        flags.append("partial_protein_support")
    else:
        flags.append("weak_protein_support")
    if len(supporting_sources) == 1:
        flags.append("single_source_support_only")
    elif len(supporting_sources) > 1:
        flags.append("multi_source_support")
    if feature_counts["cds_introns"][0] > 0:
        flags.append("protein_supports_introns")
    if feature_counts["cds_splice_sites"][0] > 0:
        flags.append("protein_supports_splice_sites")
    if feature_counts["internal_cds_exons"][0] > 0:
        flags.append("protein_supports_internal_exons")
    if feature_counts["terminal_cds_exons"][0] > 0:
        flags.append("protein_supports_terminal_exons")
    if supported_cds_fraction > 0 and not any(
        feature_counts[name][0] > 0 for name in ("cds_introns", "cds_splice_sites", "internal_cds_exons", "terminal_cds_exons")
    ):
        flags.append("protein_supports_cds_bases_only")
    return fmt_flags(flags)


def make_protein_sublocus_id(superlocus_id: str, index: int) -> str:
    return f"{superlocus_id}.protein_sublocus{index:03d}"


def protein_sublocus_type(items: list[ProteinEvidenceCollapsedAlignment]) -> str:
    if all(item.is_single_exon for item in items):
        return "single_exon"
    if all(item.is_multi_exon for item in items):
        return "multi_exon"
    return "mixed"


def build_protein_subloci_for_superlocus(sl: ProteinEvidenceSuperlocus, *, single_exon_min_reciprocal_overlap: float) -> list[ProteinEvidenceSublocus]:
    items = sorted(sl.collapsed_alignments, key=lambda a: (a.cds_start0, a.cds_end0, a.collapsed_id))
    if not items:
        return []
    uf = UnionFind([a.collapsed_id for a in items])

    splice_index: dict[int, list[ProteinEvidenceCollapsedAlignment]] = defaultdict(list)
    for aln in items:
        if aln.is_multi_exon:
            for site in aln.cds_splice_sites:
                splice_index[site].append(aln)
    for site_items in splice_index.values():
        site_items = sorted(site_items, key=lambda a: a.collapsed_id)
        for i, a in enumerate(site_items):
            for b in site_items[i + 1:]:
                uf.union(a.collapsed_id, b.collapsed_id)

    singles = [a for a in items if a.is_single_exon]
    singles = sorted(singles, key=lambda a: (a.cds_start0, a.cds_end0, a.collapsed_id))
    for i, a in enumerate(singles):
        j = i + 1
        while j < len(singles) and singles[j].cds_start0 < a.cds_end0:
            b = singles[j]
            _overlap, ro_a, ro_b = reciprocal_overlap(a.cds_start0, a.cds_end0, b.cds_start0, b.cds_end0)
            if min(ro_a, ro_b) >= single_exon_min_reciprocal_overlap:
                uf.union(a.collapsed_id, b.collapsed_id)
            j += 1

    comps: dict[str, list[ProteinEvidenceCollapsedAlignment]] = defaultdict(list)
    for aln in items:
        comps[uf.find(aln.collapsed_id)].append(aln)
    ordered = sorted(comps.values(), key=lambda comp: (min(a.cds_start0 for a in comp), max(a.cds_end0 for a in comp), min(a.collapsed_id for a in comp)))
    out: list[ProteinEvidenceSublocus] = []
    for idx, comp in enumerate(ordered, start=1):
        sid = make_protein_sublocus_id(sl.evidence_superlocus_id, idx)
        start0 = min(a.cds_start0 for a in comp); end0 = max(a.cds_end0 for a in comp)
        stype = protein_sublocus_type(comp)
        comp_sorted = sorted(comp, key=lambda a: a.collapsed_id)
        sub = ProteinEvidenceSublocus(
            sid,
            sl.evidence_superlocus_id,
            sl.chrom,
            sl.strand,
            start0,
            end0,
            comp_sorted,
            stype,
            protein_sublocus_flags(comp_sorted, stype, start0, end0),
        )
        for aln in sub.collapsed_alignments:
            aln.evidence_sublocus_id = sid
        out.append(sub)
    return out


def build_protein_subloci(superloci: list[ProteinEvidenceSuperlocus], *, single_exon_min_reciprocal_overlap: float) -> list[ProteinEvidenceSublocus]:
    out: list[ProteinEvidenceSublocus] = []
    for sl in superloci:
        out.extend(build_protein_subloci_for_superlocus(sl, single_exon_min_reciprocal_overlap=single_exon_min_reciprocal_overlap))
    return out


def terminal_exons_with_side(
    obj: ProteinEvidenceCollapsedAlignment | CodingModel,
) -> tuple[tuple[str, str, int, int, str], ...]:
    """Return compact biological terminal-exon records without object caches."""
    parts = obj.sorted_cds_parts
    if not parts:
        return tuple()
    if len(parts) == 1:
        start0, end0 = parts[0]
        return ((obj.chrom, obj.strand, start0, end0, "single"),)
    if obj.strand == "+":
        items = ((parts[0], "5_prime"), (parts[-1], "3_prime"))
    else:
        items = ((parts[-1], "5_prime"), (parts[0], "3_prime"))
    return tuple((obj.chrom, obj.strand, s, e, side) for (s, e), side in items)


def terminal_match(
    a: tuple[str, str, int, int, str],
    b: tuple[str, str, int, int, str],
    wobble_bp: int,
) -> bool:
    a_chrom, a_strand, a_start, a_end, a_side = a
    b_chrom, b_strand, b_start, b_end, b_side = b
    if a_chrom != b_chrom or a_strand != b_strand or a_side != b_side:
        return False
    if a_side != "3_prime" or wobble_bp <= 0:
        return a_start == b_start and a_end == b_end
    if a_strand == "+":
        return a_start == b_start and abs(a_end - b_end) <= wobble_bp
    if a_strand == "-":
        return a_end == b_end and abs(a_start - b_start) <= wobble_bp
    return a_start == b_start and a_end == b_end


def terminal_prf(query: ProteinEvidenceCollapsedAlignment | CodingModel, target: ProteinEvidenceCollapsedAlignment | CodingModel, wobble_bp: int) -> dict[str, Optional[float]]:
    q = terminal_exons_with_side(query); t = terminal_exons_with_side(target)
    q_supported = sum(1 for qf in q if any(terminal_match(qf, tf, wobble_bp) for tf in t))
    t_supported = sum(1 for tf in t if any(terminal_match(qf, tf, wobble_bp) for qf in q))
    return prf_counts(q_supported, len(q), t_supported, len(t))


def prf_counts(query_supported: int, query_total: int, target_supported: int, target_total: int) -> dict[str, Optional[float]]:
    precision = query_supported / query_total if query_total else None
    sensitivity = target_supported / target_total if target_total else None
    if precision is None or sensitivity is None:
        f1 = None
    elif precision + sensitivity <= 0:
        f1 = 0.0
    else:
        f1 = 2 * precision * sensitivity / (precision + sensitivity)
    return {"precision": precision, "sensitivity": sensitivity, "f1": f1}


def exact_prf(query_features: set | frozenset, target_features: set | frozenset) -> dict[str, Optional[float]]:
    shared = len(query_features & target_features)
    return prf_counts(shared, len(query_features), shared, len(target_features))


def model_feature_sets(model: CodingModel) -> dict[str, frozenset]:
    parts = model.sorted_cds_parts
    internal = parts[1:-1] if len(parts) > 2 else tuple()
    terminal = parts if len(parts) == 1 else ((parts[0], parts[-1]) if parts else tuple())
    return {
        "cds_introns": frozenset(model.cds_introns),
        "cds_splice_sites": model.cds_splice_sites,
        "internal_cds_exons": frozenset(internal),
        "terminal_cds_exons": frozenset(terminal),
        "cds_exons": frozenset(parts),
    }


def evidence_feature_sets(aln: ProteinEvidenceCollapsedAlignment) -> dict[str, frozenset]:
    parts = aln.sorted_cds_parts
    internal = parts[1:-1] if len(parts) > 2 else tuple()
    terminal = parts if len(parts) == 1 else ((parts[0], parts[-1]) if parts else tuple())
    return {
        "cds_introns": aln.cds_intron_set,
        "cds_splice_sites": aln.cds_splice_sites,
        "internal_cds_exons": frozenset(internal),
        "terminal_cds_exons": frozenset(terminal),
        "cds_exons": frozenset(parts),
    }


def chains_match_with_3prime_wobble(model: CodingModel, aln: ProteinEvidenceCollapsedAlignment, wobble_bp: int) -> bool:
    if model.chrom != aln.chrom or model.strand != aln.strand:
        return False
    m_parts = model.sorted_cds_parts; a_parts = aln.sorted_cds_parts
    if len(m_parts) != len(a_parts) or not m_parts:
        return False
    three_prime_index = len(m_parts) - 1 if model.strand == "+" else 0
    for idx, ((ms, me), (as_, ae)) in enumerate(zip(m_parts, a_parts)):
        if idx == three_prime_index and wobble_bp > 0:
            if model.strand == "+":
                if ms != as_ or abs(me - ae) > wobble_bp:
                    return False
            else:
                if me != ae or abs(ms - as_) > wobble_bp:
                    return False
        elif ms != as_ or me != ae:
            return False
    return True


def pairwise_model_protein_metrics(
    model: CodingModel,
    aln: ProteinEvidenceCollapsedAlignment,
    *,
    wobble_bp: int,
    model_features: Optional[dict[str, frozenset]] = None,
    evidence_features: Optional[dict[str, frozenset]] = None,
) -> dict[str, Optional[float]]:
    base = interval_prf_sorted(model.sorted_cds_parts, aln.sorted_cds_parts)
    mf = model_features if model_features is not None else model_feature_sets(model)
    ef = evidence_features if evidence_features is not None else evidence_feature_sets(aln)
    intron = exact_prf(mf["cds_introns"], ef["cds_introns"])
    splice = exact_prf(mf["cds_splice_sites"], ef["cds_splice_sites"])
    if mf["internal_cds_exons"]:
        internal = exact_prf(mf["internal_cds_exons"], ef["internal_cds_exons"])
        # Internal-exon support is model-applicable here. An alignment with no
        # internal exons provides genuine zero support rather than making the
        # model's internal-exon F1 unavailable.
        if internal["f1"] is None:
            internal["f1"] = 0.0
    else:
        internal = {"precision": None, "sensitivity": None, "f1": None}
    terminal = terminal_prf(model, aln, wobble_bp)
    # Mixed single-/multi-exon comparisons use CDS-base F1 for the aggregate
    # protein score. Terminal-exon metrics remain populated for diagnostics,
    # but must not introduce an asymmetric penalty depending on which side is
    # the single-exon object.
    if model.is_single_exon or aln.is_single_exon:
        score = base["f1"] if base["f1"] is not None else 0.0
    else:
        score = weighted_metric_value([
            ("cds_intron_f1", 0.30, intron["f1"]),
            ("cds_splice_site_f1", 0.25, splice["f1"]),
            ("internal_cds_exon_f1", 0.20, internal["f1"]),
            ("cds_base_f1", 0.20, base["f1"]),
            ("terminal_cds_exon_f1", 0.05, terminal["f1"]),
        ]) or 0.0
    return {
        "score": score,
        "cds_base_precision": base["precision"], "cds_base_sensitivity": base["sensitivity"], "cds_base_f1": base["f1"],
        "cds_intron_precision": intron["precision"], "cds_intron_sensitivity": intron["sensitivity"], "cds_intron_f1": intron["f1"],
        "cds_splice_site_precision": splice["precision"], "cds_splice_site_sensitivity": splice["sensitivity"], "cds_splice_site_f1": splice["f1"],
        "internal_cds_exon_precision": internal["precision"], "internal_cds_exon_sensitivity": internal["sensitivity"], "internal_cds_exon_f1": internal["f1"],
        "terminal_cds_exon_precision": terminal["precision"], "terminal_cds_exon_sensitivity": terminal["sensitivity"], "terminal_cds_exon_f1": terminal["f1"],
    }


def fmt_counts(counts: dict[str, int]) -> str:
    return ",".join(f"{k}:{v}" for k, v in sorted(counts.items())) if counts else "."


def median_int(values: list[int]) -> float:
    if not values:
        return 0.0
    ordered = sorted(int(v) for v in values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[mid])
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def fmt_int_summary(values: list[int]) -> str:
    if not values:
        return "min=0,median=0,max=0"
    return f"min={min(values):,},median={median_int(values):,.1f},max={max(values):,}"


def protein_summary_counts(
    *,
    protein_paths: list[Path],
    protein_source_labels: list[str],
    raw_alignments: list[ProteinEvidenceRawAlignment],
    superloci: list[ProteinEvidenceSuperlocus],
    collapsed_alignments: list[ProteinEvidenceCollapsedAlignment],
    subloci: list[ProteinEvidenceSublocus],
    model_support_rows: list[dict[str, object]],
    raw_to_collapsed_mapping_requested: bool,
) -> dict[str, object]:
    """Return compact protein-evidence summary counts for logs and TSV output."""
    raw_by_source = Counter(a.source_label for a in raw_alignments)
    collapsed_by_source: Counter[str] = Counter()
    collapsed_member_by_source: Counter[str] = Counter()
    for aln in collapsed_alignments:
        for source, count in aln.source_member_counts.items():
            if count > 0:
                collapsed_by_source[source] += 1
                collapsed_member_by_source[source] += count
    sublocus_types = Counter(sub.sublocus_type for sub in subloci)
    model_status = Counter(str(row.get("protein_evidence_status", ".")) for row in model_support_rows)
    flag_counts: Counter[str] = Counter()
    for row in model_support_rows:
        flags = str(row.get("protein_evidence_flags", "."))
        if flags and flags != ".":
            for flag in flags.split(','):
                if flag:
                    flag_counts[flag] += 1
    source_raw_values = [raw_by_source.get(label, 0) for label in protein_source_labels]
    source_collapsed_values = [collapsed_by_source.get(label, 0) for label in protein_source_labels]
    source_member_values = [collapsed_member_by_source.get(label, 0) for label in protein_source_labels]
    collapse_ratio = (len(collapsed_alignments) / len(raw_alignments)) if raw_alignments else 0.0
    rows_total = len(model_support_rows)
    models_with_overlap = sum(1 for row in model_support_rows if int(row.get("protein_overlapping_collapsed_alignment_count", 0) or 0) > 0)
    models_with_exact = sum(1 for row in model_support_rows if str(row.get("protein_exact_cds_chain_supported", "false")) == "true")
    models_with_supporting_sources = sum(1 for row in model_support_rows if int(row.get("protein_supporting_source_count", 0) or 0) > 0)
    models_with_multi_source_support = flag_counts.get("multi_source_support", 0)
    models_with_strong_support = flag_counts.get("strong_protein_support", 0)
    models_with_partial_support = flag_counts.get("partial_protein_support", 0)
    models_with_weak_support = flag_counts.get("weak_protein_support", 0)
    summary = {
        "file_count": len(protein_paths),
        "source_count": len(protein_source_labels),
        "source_labels": ",".join(protein_source_labels) if protein_source_labels else ".",
        "input_files": ",".join(str(p) for p in protein_paths) if protein_paths else ".",
        "raw_alignment_count": len(raw_alignments),
        "collapsed_alignment_count": len(collapsed_alignments),
        "collapse_ratio": collapse_ratio,
        "raw_to_collapsed_mapping_requested": raw_to_collapsed_mapping_requested,
        "superlocus_count": len(superloci),
        "sublocus_count": len(subloci),
        "single_exon_sublocus_count": sublocus_types.get("single_exon", 0),
        "multi_exon_sublocus_count": sublocus_types.get("multi_exon", 0),
        "mixed_sublocus_count": sublocus_types.get("mixed", 0),
        "model_support_row_count": rows_total,
        "models_with_overlapping_protein_evidence": models_with_overlap,
        "models_with_exact_chain_support": models_with_exact,
        "models_with_supporting_protein_source": models_with_supporting_sources,
        "models_with_multi_source_support": models_with_multi_source_support,
        "models_with_strong_protein_support": models_with_strong_support,
        "models_with_partial_protein_support": models_with_partial_support,
        "models_with_weak_protein_support": models_with_weak_support,
        "models_with_no_overlapping_protein_evidence": model_status.get("no_overlapping_protein_evidence", 0),
        "model_status_counts": fmt_counts(dict(model_status)),
        "model_flag_counts": fmt_counts(dict(flag_counts)),
        "per_source_raw_alignment_summary": fmt_int_summary(source_raw_values),
        "per_source_collapsed_alignment_summary": fmt_int_summary(source_collapsed_values),
        "per_source_member_count_summary": fmt_int_summary(source_member_values),
        "per_source_raw_alignment_counts": fmt_counts(dict(raw_by_source)),
        "per_source_collapsed_alignment_counts": fmt_counts(dict(collapsed_by_source)),
    }
    return summary


def protein_summary_log_lines(summary: dict[str, object]) -> list[str]:
    return [
        (
            "Protein evidence inputs: "
            f"files={int(summary['file_count']):,}; sources={int(summary['source_count']):,}; "
            f"raw_alignments={int(summary['raw_alignment_count']):,}; "
            f"collapsed_alignments={int(summary['collapsed_alignment_count']):,}; "
            f"collapse_ratio={float(summary['collapse_ratio']):.3f}; "
            f"raw_to_collapsed_mapping_requested={str(summary['raw_to_collapsed_mapping_requested']).lower()}"
        ),
        (
            "Protein evidence source summaries: "
            f"raw_alignments({summary['per_source_raw_alignment_summary']}); "
            f"collapsed_alignments({summary['per_source_collapsed_alignment_summary']})"
        ),
        (
            "Protein evidence loci: "
            f"superloci={int(summary['superlocus_count']):,}; subloci={int(summary['sublocus_count']):,}; "
            f"single_exon_subloci={int(summary['single_exon_sublocus_count']):,}; "
            f"multi_exon_subloci={int(summary['multi_exon_sublocus_count']):,}; "
            f"mixed_subloci={int(summary['mixed_sublocus_count']):,}"
        ),
        (
            "Protein evidence model associations: "
            f"models={int(summary['model_support_row_count']):,}; "
            f"with_overlap={int(summary['models_with_overlapping_protein_evidence']):,}; "
            f"exact_chain={int(summary['models_with_exact_chain_support']):,}; "
            f"supporting_source={int(summary['models_with_supporting_protein_source']):,}; "
            f"multi_source_support={int(summary['models_with_multi_source_support']):,}; "
            f"strong={int(summary['models_with_strong_protein_support']):,}; "
            f"partial={int(summary['models_with_partial_protein_support']):,}; "
            f"weak={int(summary['models_with_weak_protein_support']):,}; "
            f"no_overlap={int(summary['models_with_no_overlapping_protein_evidence']):,}"
        ),
    ]


def protein_summary_tsv_rows(
    *,
    protein_paths: list[Path],
    protein_source_labels: list[str],
    raw_alignments: list[ProteinEvidenceRawAlignment],
    superloci: list[ProteinEvidenceSuperlocus],
    collapsed_alignments: list[ProteinEvidenceCollapsedAlignment],
    subloci: list[ProteinEvidenceSublocus],
    model_support_rows: list[dict[str, object]],
    raw_to_collapsed_mapping_requested: bool,
) -> list[dict[str, object]]:
    summary = protein_summary_counts(
        protein_paths=protein_paths,
        protein_source_labels=protein_source_labels,
        raw_alignments=raw_alignments,
        superloci=superloci,
        collapsed_alignments=collapsed_alignments,
        subloci=subloci,
        model_support_rows=model_support_rows,
        raw_to_collapsed_mapping_requested=raw_to_collapsed_mapping_requested,
    )
    rows: list[dict[str, object]] = []

    def add(scope: str, metric: str, value: object, details: str = ".") -> None:
        rows.append({"summary_scope": scope, "source_label": ".", "input_file": ".", "metric": metric, "value": value, "details": details})

    for metric in (
        "file_count", "source_count", "source_labels", "input_files", "raw_alignment_count",
        "collapsed_alignment_count", "collapse_ratio", "raw_to_collapsed_mapping_requested",
        "superlocus_count", "sublocus_count", "single_exon_sublocus_count", "multi_exon_sublocus_count",
        "mixed_sublocus_count", "model_support_row_count", "models_with_overlapping_protein_evidence",
        "models_with_exact_chain_support", "models_with_supporting_protein_source", "models_with_multi_source_support",
        "models_with_strong_protein_support", "models_with_partial_protein_support", "models_with_weak_protein_support",
        "models_with_no_overlapping_protein_evidence", "model_status_counts", "model_flag_counts",
        "per_source_raw_alignment_summary", "per_source_collapsed_alignment_summary", "per_source_member_count_summary",
        "per_source_raw_alignment_counts", "per_source_collapsed_alignment_counts",
    ):
        add("overall", metric, summary[metric])

    raw_by_source = Counter(a.source_label for a in raw_alignments)
    collapsed_by_source: Counter[str] = Counter()
    collapsed_members_by_source: Counter[str] = Counter()
    exact_chain_multi_source_by_source: Counter[str] = Counter()
    for aln in collapsed_alignments:
        for source, count in aln.source_member_counts.items():
            if count > 0:
                collapsed_by_source[source] += 1
                collapsed_members_by_source[source] += count
                if aln.source_count > 1:
                    exact_chain_multi_source_by_source[source] += 1
    support_by_source: Counter[str] = Counter()
    overlap_by_source: Counter[str] = Counter()
    for row in model_support_rows:
        for source in str(row.get("protein_overlapping_sources", ".")).split(','):
            if source and source != ".":
                overlap_by_source[source] += 1
        for source in str(row.get("protein_supporting_sources", ".")).split(','):
            if source and source != ".":
                support_by_source[source] += 1
        # Exact-chain source names are not written explicitly yet, so exact-by-source
        # is approximated by collapsed IDs only in the overall exact count.
    path_by_label = {label: str(path) for label, path in zip(protein_source_labels, protein_paths)}
    for label in protein_source_labels:
        detail = f"source_index={protein_source_labels.index(label) + 1}"
        rows.extend([
            {"summary_scope": "source", "source_label": label, "input_file": path_by_label.get(label, "."), "metric": "raw_alignment_count", "value": raw_by_source.get(label, 0), "details": detail},
            {"summary_scope": "source", "source_label": label, "input_file": path_by_label.get(label, "."), "metric": "collapsed_alignment_count", "value": collapsed_by_source.get(label, 0), "details": detail},
            {"summary_scope": "source", "source_label": label, "input_file": path_by_label.get(label, "."), "metric": "member_count_in_collapsed_alignments", "value": collapsed_members_by_source.get(label, 0), "details": detail},
            {"summary_scope": "source", "source_label": label, "input_file": path_by_label.get(label, "."), "metric": "multi_source_exact_chain_alignment_count", "value": exact_chain_multi_source_by_source.get(label, 0), "details": detail},
            {"summary_scope": "source", "source_label": label, "input_file": path_by_label.get(label, "."), "metric": "models_with_overlap", "value": overlap_by_source.get(label, 0), "details": detail},
            {"summary_scope": "source", "source_label": label, "input_file": path_by_label.get(label, "."), "metric": "models_with_supporting_source", "value": support_by_source.get(label, 0), "details": detail},
        ])
    return rows


def write_protein_summary_tsv(rows: list[dict[str, object]], path: Path) -> None:
    write_dict_tsv(rows, path, ["summary_scope", "source_label", "input_file", "metric", "value", "details"])


def associated_evidence_for_model(
    model: CodingModel,
    superlocus_index: IntervalOverlapIndex[tuple[str, str], ProteinEvidenceSuperlocus],
) -> list[ProteinEvidenceCollapsedAlignment]:
    """Return overlapping collapsed protein alignments using an indexed locus lookup."""
    out: list[ProteinEvidenceCollapsedAlignment] = []
    for sl in superlocus_index.query(
        (model.chrom, model.strand),
        model.cds_start0,
        model.cds_end0,
    ):
        for aln in sl.collapsed_alignments:
            if interval_list_overlap_bp_sorted(model.sorted_cds_parts, aln.sorted_cds_parts) > 0:
                out.append(aln)
    return out


def calculate_protein_model_support(
    models: list[CodingModel],
    protein_superloci: list[ProteinEvidenceSuperlocus],
    protein_source_labels: list[str],
    *,
    terminal_cds_3prime_wobble_bp: int,
    support_score_threshold: float = 0.50,
) -> list[dict[str, object]]:
    superlocus_index = IntervalOverlapIndex(
        protein_superloci,
        key=lambda sl: (sl.chrom, sl.strand),
        start=lambda sl: sl.start0,
        end=lambda sl: sl.end0,
    )
    total_sources = len(protein_source_labels)
    rows: list[dict[str, object]] = []
    for model in sorted(models, key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.output_transcript_id)):
        alns = associated_evidence_for_model(model, superlocus_index)
        evidence_super_ids = {a.evidence_superlocus_id for a in alns}
        evidence_sub_ids = {a.evidence_sublocus_id for a in alns if a.evidence_sublocus_id}
        overlapping_sources = sorted({s for a in alns for s in a.sources})
        exact = [a for a in alns if chains_match_with_3prime_wobble(model, a, terminal_cds_3prime_wobble_bp)]
        exact_sources = sorted({s for a in exact for s in a.sources})
        exact_member_count = sum(a.member_count for a in exact)

        # Feature support: a model feature is supported if any overlapping protein source/alignment has that exact feature.
        mf = model_feature_sets(model)
        # Keep alignment geometry only for this model. Patch 16 attached large
        # feature-set dictionaries to every model/alignment for the whole run;
        # this bounded cache preserves reuse within a model while allowing the
        # memory to be reclaimed on the next iteration.
        aln_feature_sets = {a.collapsed_id: evidence_feature_sets(a) for a in alns}
        union_ef: dict[str, set[tuple]] = {k: set() for k in mf}
        for aln in alns:
            ef = aln_feature_sets[aln.collapsed_id]
            for k in union_ef:
                union_ef[k].update(ef[k])
        feature_counts: dict[str, tuple[int, int, Optional[float]]] = {}
        for k in mf:
            total = len(mf[k])
            supported = len(mf[k] & union_ef[k])
            if total:
                fraction = supported / total
            elif k == "internal_cds_exons":
                fraction = None
            else:
                fraction = 0.0
            feature_counts[k] = (supported, total, fraction)
        supported_cds_bp = 0
        if alns:
            supported_cds_bp = interval_list_overlap_bp(list(model.cds_parts), [part for aln in alns for part in aln.cds_parts])
        supported_cds_fraction = supported_cds_bp / model.cds_length if model.cds_length else 0.0

        # Source-collapsed best pairwise metrics.
        best_by_source: dict[str, tuple[ProteinEvidenceCollapsedAlignment, dict[str, Optional[float]]]] = {}
        best_alignment: Optional[ProteinEvidenceCollapsedAlignment] = None
        best_alignment_metrics: Optional[dict[str, Optional[float]]] = None
        for aln in alns:
            metrics = pairwise_model_protein_metrics(
                model,
                aln,
                wobble_bp=terminal_cds_3prime_wobble_bp,
                model_features=mf,
                evidence_features=aln_feature_sets[aln.collapsed_id],
            )
            if best_alignment_metrics is None or float(metrics.get("score") or 0.0) > float(best_alignment_metrics.get("score") or 0.0):
                best_alignment = aln; best_alignment_metrics = metrics
            for source in aln.sources:
                old = best_by_source.get(source)
                if old is None or float(metrics.get("score") or 0.0) > float(old[1].get("score") or 0.0):
                    best_by_source[source] = (aln, metrics)

        def mean_metric(name: str, *, model_applicable: bool = True) -> Optional[float]:
            if not model_applicable:
                return None
            vals = [float(m[name]) for _aln, m in best_by_source.values() if m.get(name) is not None]
            return sum(vals) / len(vals) if vals else 0.0

        has_internal_cds_exons = feature_counts["internal_cds_exons"][1] > 0

        best_score = float(best_alignment_metrics.get("score") or 0.0) if best_alignment_metrics else 0.0
        mean_best_source_score = mean_metric("score")
        supporting_sources = sorted(s for s, (_aln, m) in best_by_source.items() if float(m.get("score") or 0.0) >= support_score_threshold)
        supporting_source_fraction = len(supporting_sources) / total_sources if total_sources else 0.0
        exact_source_fraction = len(exact_sources) / total_sources if total_sources else 0.0

        if model.is_single_exon:
            evidence_score, score_components = weighted_metric_value_with_components([
                ("protein_cds_base_f1", 0.70, mean_metric("cds_base_f1")),
                ("protein_exact_cds_chain_source_fraction", 0.15, exact_source_fraction),
                ("protein_supporting_source_fraction", 0.15, supporting_source_fraction),
            ])
            score_type = "single_exon"
        else:
            evidence_score, score_components = weighted_metric_value_with_components([
                ("protein_cds_intron_f1", 0.25, mean_metric("cds_intron_f1")),
                ("protein_cds_splice_site_f1", 0.20, mean_metric("cds_splice_site_f1")),
                ("protein_internal_cds_exon_f1", 0.20, mean_metric("internal_cds_exon_f1", model_applicable=has_internal_cds_exons)),
                ("protein_cds_base_f1", 0.15, mean_metric("cds_base_f1")),
                ("protein_exact_cds_chain_source_fraction", 0.10, exact_source_fraction),
                ("protein_supporting_source_fraction", 0.10, supporting_source_fraction),
            ])
            score_type = "multi_exon"
        if not protein_source_labels:
            status = "protein_evidence_not_configured"
        elif not alns:
            status = "no_overlapping_protein_evidence"
        elif supporting_sources:
            status = "protein_evidence_supporting"
        else:
            status = "protein_evidence_overlap_low_similarity"

        evidence_flags = protein_model_support_flags(
            protein_source_labels=protein_source_labels,
            alns=alns,
            exact=exact,
            supporting_sources=supporting_sources,
            evidence_score=evidence_score,
            support_score_threshold=support_score_threshold,
            feature_counts=feature_counts,
            supported_cds_fraction=supported_cds_fraction,
        )

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
            "protein_evidence_source_total": total_sources,
            "protein_evidence_superlocus_count": len(evidence_super_ids),
            "protein_evidence_sublocus_count": len(evidence_sub_ids),
            "protein_overlapping_source_count": len(overlapping_sources),
            "protein_overlapping_sources": ",".join(overlapping_sources) if overlapping_sources else ".",
            "protein_overlapping_collapsed_alignment_count": len(alns),
            "protein_exact_cds_chain_supported": "true" if exact else "false",
            "protein_exact_cds_chain_source_count": len(exact_sources),
            "protein_exact_cds_chain_source_fraction": exact_source_fraction,
            "protein_exact_cds_chain_member_count": exact_member_count,
            "protein_exact_cds_chain_collapsed_ids": ",".join(a.collapsed_id for a in exact) if exact else ".",
            "protein_model_cds_intron_count": feature_counts["cds_introns"][1],
            "protein_supported_cds_intron_count": feature_counts["cds_introns"][0],
            "protein_supported_cds_intron_fraction": feature_counts["cds_introns"][2],
            "protein_model_cds_splice_site_count": feature_counts["cds_splice_sites"][1],
            "protein_supported_cds_splice_site_count": feature_counts["cds_splice_sites"][0],
            "protein_supported_cds_splice_site_fraction": feature_counts["cds_splice_sites"][2],
            "protein_model_internal_cds_exon_count": feature_counts["internal_cds_exons"][1],
            "protein_supported_internal_cds_exon_count": feature_counts["internal_cds_exons"][0],
            "protein_supported_internal_cds_exon_fraction": feature_counts["internal_cds_exons"][2],
            "protein_model_terminal_cds_exon_count": feature_counts["terminal_cds_exons"][1],
            "protein_supported_terminal_cds_exon_count": feature_counts["terminal_cds_exons"][0],
            "protein_supported_terminal_cds_exon_fraction": feature_counts["terminal_cds_exons"][2],
            "protein_supported_cds_base_bp": supported_cds_bp,
            "protein_supported_cds_base_fraction": supported_cds_fraction,
            "protein_best_alignment_score": best_score,
            "protein_best_alignment_id": best_alignment.collapsed_id if best_alignment else ".",
            "protein_mean_best_source_score": mean_best_source_score,
            "protein_supporting_source_count": len(supporting_sources),
            "protein_supporting_source_fraction": supporting_source_fraction,
            "protein_supporting_sources": ",".join(supporting_sources) if supporting_sources else ".",
            "protein_cds_base_precision": mean_metric("cds_base_precision"),
            "protein_cds_base_sensitivity": mean_metric("cds_base_sensitivity"),
            "protein_cds_base_f1": mean_metric("cds_base_f1"),
            "protein_cds_intron_precision": mean_metric("cds_intron_precision"),
            "protein_cds_intron_sensitivity": mean_metric("cds_intron_sensitivity"),
            "protein_cds_intron_f1": mean_metric("cds_intron_f1"),
            "protein_cds_splice_site_precision": mean_metric("cds_splice_site_precision"),
            "protein_cds_splice_site_sensitivity": mean_metric("cds_splice_site_sensitivity"),
            "protein_cds_splice_site_f1": mean_metric("cds_splice_site_f1"),
            "protein_internal_cds_exon_precision": mean_metric("internal_cds_exon_precision", model_applicable=has_internal_cds_exons),
            "protein_internal_cds_exon_sensitivity": mean_metric("internal_cds_exon_sensitivity", model_applicable=has_internal_cds_exons),
            "protein_internal_cds_exon_f1": mean_metric("internal_cds_exon_f1", model_applicable=has_internal_cds_exons),
            "protein_terminal_cds_exon_precision": mean_metric("terminal_cds_exon_precision"),
            "protein_terminal_cds_exon_sensitivity": mean_metric("terminal_cds_exon_sensitivity"),
            "protein_terminal_cds_exon_f1": mean_metric("terminal_cds_exon_f1"),
            "protein_evidence_score": evidence_score,
            "protein_evidence_score_scaled": scale_protein_evidence_score(evidence_score),
            "protein_evidence_score_type": score_type,
            "protein_evidence_score_components": score_components,
            "protein_evidence_status": status,
            "protein_evidence_flags": evidence_flags,
        }
        rows.append(row)
    return rows


PROTEIN_MODEL_SUPPORT_COLUMNS = [
    "superlocus_id", "sublocus_id", "source_name", "input_gene_id", "input_transcript_id", "output_gene_id", "output_transcript_id",
    "chrom", "strand", "cds_start", "cds_end", "cds_length", "cds_exon_count",
    "protein_evidence_source_total", "protein_evidence_superlocus_count", "protein_evidence_sublocus_count",
    "protein_overlapping_source_count", "protein_overlapping_sources", "protein_overlapping_collapsed_alignment_count",
    "protein_exact_cds_chain_supported", "protein_exact_cds_chain_source_count", "protein_exact_cds_chain_source_fraction", "protein_exact_cds_chain_member_count", "protein_exact_cds_chain_collapsed_ids",
    "protein_model_cds_intron_count", "protein_supported_cds_intron_count", "protein_supported_cds_intron_fraction",
    "protein_model_cds_splice_site_count", "protein_supported_cds_splice_site_count", "protein_supported_cds_splice_site_fraction",
    "protein_model_internal_cds_exon_count", "protein_supported_internal_cds_exon_count", "protein_supported_internal_cds_exon_fraction",
    "protein_model_terminal_cds_exon_count", "protein_supported_terminal_cds_exon_count", "protein_supported_terminal_cds_exon_fraction",
    "protein_supported_cds_base_bp", "protein_supported_cds_base_fraction",
    "protein_best_alignment_score", "protein_best_alignment_id", "protein_mean_best_source_score", "protein_supporting_source_count", "protein_supporting_source_fraction", "protein_supporting_sources",
    "protein_cds_base_precision", "protein_cds_base_sensitivity", "protein_cds_base_f1",
    "protein_cds_intron_precision", "protein_cds_intron_sensitivity", "protein_cds_intron_f1",
    "protein_cds_splice_site_precision", "protein_cds_splice_site_sensitivity", "protein_cds_splice_site_f1",
    "protein_internal_cds_exon_precision", "protein_internal_cds_exon_sensitivity", "protein_internal_cds_exon_f1",
    "protein_terminal_cds_exon_precision", "protein_terminal_cds_exon_sensitivity", "protein_terminal_cds_exon_f1",
    "protein_evidence_score", "protein_evidence_score_scaled", "protein_evidence_score_type", "protein_evidence_score_components", "protein_evidence_status", "protein_evidence_flags",
]


def write_protein_raw_to_collapsed_tsv(raw_alignments: list[ProteinEvidenceRawAlignment], path: Path) -> None:
    columns = ["protein_evidence_raw_id", "protein_evidence_collapsed_id", "protein_evidence_superlocus_id", "protein_source_label", "source_index", "input_file", "original_alignment_id", "chrom", "strand", "cds_start", "cds_end", "cds_length", "cds_exon_count", "cds_chain"]
    with path.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for a in sorted(raw_alignments, key=lambda x: (x.chrom, x.cds_start0, x.cds_end0, x.source_label, x.raw_id)):
            row = (
                a.raw_id,
                a.collapsed_id or ".",
                a.evidence_superlocus_id or ".",
                a.source_label,
                a.source_index,
                str(a.input_file),
                a.original_alignment_id,
                a.chrom,
                a.strand,
                a.cds_start0 + 1,
                a.cds_end0,
                a.cds_length,
                a.cds_exon_count,
                ",".join(f"{s+1}-{e}" for s, e in a.cds_parts),
            )
            out.write("\t".join(format_tsv_value(value) for value in row) + "\n")


def write_protein_collapsed_alignments_tsv(collapsed: list[ProteinEvidenceCollapsedAlignment], path: Path) -> None:
    rows = [{
        "protein_evidence_collapsed_id": a.collapsed_id,
        "protein_evidence_superlocus_id": a.evidence_superlocus_id,
        "protein_evidence_sublocus_id": a.evidence_sublocus_id or ".",
        "chrom": a.chrom,
        "strand": a.strand,
        "cds_start": a.cds_start0 + 1,
        "cds_end": a.cds_end0,
        "cds_length": a.cds_length,
        "cds_exon_count": a.cds_exon_count,
        "member_count": a.member_count,
        "source_count": a.source_count,
        "source_member_counts": fmt_counts(a.source_member_counts),
        "representative_raw_id": a.representative_raw_id,
        "flags": a.flags,
        "cds_chain": ",".join(f"{s+1}-{e}" for s, e in a.cds_parts),
    } for a in sorted(collapsed, key=lambda x: (x.chrom, x.cds_start0, x.cds_end0, x.collapsed_id))]
    write_dict_tsv(rows, path, ["protein_evidence_collapsed_id", "protein_evidence_superlocus_id", "protein_evidence_sublocus_id", "chrom", "strand", "cds_start", "cds_end", "cds_length", "cds_exon_count", "member_count", "source_count", "source_member_counts", "representative_raw_id", "flags", "cds_chain"])


def write_protein_superloci_tsv(superloci: list[ProteinEvidenceSuperlocus], path: Path) -> None:
    rows = []
    for sl in superloci:
        counts = sl.raw_source_member_counts or dict(Counter(raw.source_label for raw in sl.raw_alignments))
        raw_count = sl.raw_alignment_count or len(sl.raw_alignments)
        rows.append({
            "protein_evidence_superlocus_id": sl.evidence_superlocus_id,
            "chrom": sl.chrom,
            "strand": sl.strand,
            "start": sl.start0 + 1,
            "end": sl.end0,
            "span_bp": sl.end0 - sl.start0,
            "raw_alignment_count": raw_count,
            "collapsed_alignment_count": len(sl.collapsed_alignments),
            "source_count": len(counts),
            "source_member_counts": fmt_counts(dict(counts)),
        })
    write_dict_tsv(rows, path, ["protein_evidence_superlocus_id", "chrom", "strand", "start", "end", "span_bp", "raw_alignment_count", "collapsed_alignment_count", "source_count", "source_member_counts"])


def write_protein_subloci_tsv(subloci: list[ProteinEvidenceSublocus], path: Path) -> None:
    rows = [{
        "protein_evidence_sublocus_id": sub.evidence_sublocus_id,
        "protein_evidence_superlocus_id": sub.evidence_superlocus_id,
        "chrom": sub.chrom,
        "strand": sub.strand,
        "start": sub.start0 + 1,
        "end": sub.end0,
        "span_bp": sub.end0 - sub.start0,
        "sublocus_type": sub.sublocus_type,
        "collapsed_alignment_count": len(sub.collapsed_alignments),
        "source_count": len({s for a in sub.collapsed_alignments for s in a.sources}),
        "member_count": sum(a.member_count for a in sub.collapsed_alignments),
        "flags": sub.flags,
    } for sub in subloci]
    write_dict_tsv(rows, path, ["protein_evidence_sublocus_id", "protein_evidence_superlocus_id", "chrom", "strand", "start", "end", "span_bp", "sublocus_type", "collapsed_alignment_count", "source_count", "member_count", "flags"])


def write_protein_sublocus_members_tsv(subloci: list[ProteinEvidenceSublocus], path: Path) -> None:
    rows = []
    for sub in subloci:
        for aln in sub.collapsed_alignments:
            rows.append({
                "protein_evidence_sublocus_id": sub.evidence_sublocus_id,
                "protein_evidence_superlocus_id": sub.evidence_superlocus_id,
                "protein_evidence_collapsed_id": aln.collapsed_id,
                "chrom": aln.chrom,
                "strand": aln.strand,
                "cds_start": aln.cds_start0 + 1,
                "cds_end": aln.cds_end0,
                "cds_length": aln.cds_length,
                "cds_exon_count": aln.cds_exon_count,
                "member_count": aln.member_count,
                "source_count": aln.source_count,
                "source_member_counts": fmt_counts(aln.source_member_counts),
            })
    write_dict_tsv(rows, path, ["protein_evidence_sublocus_id", "protein_evidence_superlocus_id", "protein_evidence_collapsed_id", "chrom", "strand", "cds_start", "cds_end", "cds_length", "cds_exon_count", "member_count", "source_count", "source_member_counts"])


def write_protein_model_support_tsv(rows: list[dict[str, object]], path: Path) -> None:
    write_dict_tsv(rows, path, PROTEIN_MODEL_SUPPORT_COLUMNS)
