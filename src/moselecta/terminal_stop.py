"""Early, local terminal-stop boundary canonicalisation.

This module deliberately performs no stop search.  It considers only the
current translated CDS edge and the immediately adjacent genomic bases needed
for one unique +/-1--3 bp representation correction.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import replace
from pathlib import Path
from typing import Iterable, Optional, Protocol

from .genome_sequence import (
    STANDARD_START_CODONS,
    STANDARD_STOP_CODONS,
    model_cds_phases,
    reverse_complement,
)
from .intervals import IntervalOverlapIndex
from .models import (
    CodingModel,
    FeatureRecord,
    ParsedInput,
    TerminalStopCanonicalisation,
)
from .tsvio import format_tsv_value, write_dict_tsv


class FastaReader(Protocol):
    def contig_length(self, name: str) -> Optional[int]: ...
    def fetch(self, name: str, start0: int, end0: int) -> str: ...


AUDIT_COLUMNS = [
    "source_name", "input_file", "input_gene_id", "input_transcript_id",
    "output_transcript_id", "chrom", "strand", "status", "offset",
    "original_biological_3prime_cds_boundary", "canonical_biological_3prime_cds_boundary",
    "stop_codon", "frame_confidence_basis", "first_cds_phase", "skip_reason",
    "terminal_exon_adjusted", "new_same_strand_cds_overlap_count",
    "new_same_strand_cds_overlap_models",
]


def _ordered_parts(model: CodingModel) -> list[tuple[int, int]]:
    parts = list(model.sorted_cds_parts)
    return list(reversed(parts)) if model.strand == "-" else parts


def _oriented_fetch(
    fasta: FastaReader,
    model: CodingModel,
    start0: int,
    end0: int,
) -> str:
    sequence = fasta.fetch(model.chrom, start0, end0)
    return reverse_complement(sequence) if model.strand == "-" else sequence


def _spliced_slice(
    fasta: FastaReader,
    model: CodingModel,
    start: int,
    end: int,
) -> str:
    """Fetch a short transcript-oriented CDS slice by spliced offset."""
    if start < 0 or end < start or end > model.cds_length:
        raise ValueError(f"Invalid spliced CDS slice: {start}-{end}/{model.cds_length}")
    chunks: list[str] = []
    cursor = 0
    for part_start0, part_end0 in _ordered_parts(model):
        part_length = part_end0 - part_start0
        local_start = max(start, cursor)
        local_end = min(end, cursor + part_length)
        if local_start < local_end:
            offset_start = local_start - cursor
            offset_end = local_end - cursor
            if model.strand == "+":
                genomic_start0 = part_start0 + offset_start
                genomic_end0 = part_start0 + offset_end
            else:
                genomic_start0 = part_end0 - offset_end
                genomic_end0 = part_end0 - offset_start
            chunks.append(
                _oriented_fetch(fasta, model, genomic_start0, genomic_end0)
            )
        cursor += part_length
        if cursor >= end:
            break
    sequence = "".join(chunks)
    if len(sequence) != end - start:
        raise ValueError(
            f"Could not recover complete spliced CDS slice for {model.output_transcript_id}: "
            f"requested={end - start}; recovered={len(sequence)}"
        )
    return sequence


def _biological_boundary(model: CodingModel) -> int:
    """Return the GFF-style coordinate of the biological 3' CDS boundary."""
    return model.cds_end0 if model.strand == "+" else model.cds_start0 + 1


def _adjacent_interval(model: CodingModel, count: int) -> tuple[int, int]:
    if model.strand == "+":
        return model.cds_end0, model.cds_end0 + count
    return model.cds_start0 - count, model.cds_start0


def _first_phase(model: CodingModel) -> Optional[int]:
    phases = model_cds_phases(model)
    return phases.get(_ordered_parts(model)[0])


def _frame_basis(
    model: CodingModel,
    fasta: FastaReader,
) -> tuple[Optional[int], str, Optional[str]]:
    phase = _first_phase(model)
    if phase is not None:
        return phase, f"numeric_phase_{phase}", None
    if model.cds_length < 3:
        return None, ".", "uncertain_frame"
    try:
        start = _spliced_slice(fasta, model, 0, 3)
    except (KeyError, ValueError):
        return None, ".", "uncertain_frame"
    if any(base not in "ACGT" for base in start):
        return None, ".", "uncertain_frame_ambiguous_start"
    if start not in STANDARD_START_CODONS:
        return None, ".", "uncertain_frame_missing_start"
    return 0, "start_confirmed_zero_frame", None


def _terminal_exon(model: CodingModel) -> tuple[Optional[FeatureRecord], Optional[str]]:
    terminal_part = _ordered_parts(model)[-1]
    containing = [
        record for record in model.child_records
        if record.feature_type.lower() == "exon"
        and record.seqid == model.chrom
        and record.strand == model.strand
        and record.start0 <= terminal_part[0]
        and record.end0 >= terminal_part[1]
    ]
    geometries = {(record.start0, record.end0) for record in containing}
    if not containing:
        return None, "missing_terminal_exon"
    if len(geometries) != 1:
        return None, "ambiguous_terminal_exon"
    return containing[0], None


def _extension_mode(
    model: CodingModel,
    count: int,
) -> tuple[Optional[str], Optional[FeatureRecord], Optional[str]]:
    exon, exon_error = _terminal_exon(model)
    if exon is None:
        return None, None, exon_error
    old_boundary0 = model.cds_end0 if model.strand == "+" else model.cds_start0
    exon_outer0 = exon.end0 if model.strand == "+" else exon.start0
    required_outer0 = old_boundary0 + count if model.strand == "+" else old_boundary0 - count

    if model.strand == "+":
        already_inside = required_outer0 <= exon.end0
    else:
        already_inside = required_outer0 >= exon.start0
    if already_inside:
        return "cds_only", exon, None

    if exon_outer0 != old_boundary0:
        return None, exon, "splice_crossing_or_structurally_nonlocal"

    # This must be the biological terminal exon.  A downstream exon/UTR exon is
    # never followed by this local representation normalisation.
    exons = [
        record for record in model.child_records
        if record.feature_type.lower() == "exon"
        and record.seqid == model.chrom and record.strand == model.strand
    ]
    if model.strand == "+":
        has_downstream_exon = any(record.start0 >= exon.end0 for record in exons if record is not exon)
    else:
        has_downstream_exon = any(record.end0 <= exon.start0 for record in exons if record is not exon)
    if has_downstream_exon:
        return None, exon, "downstream_exon_or_structurally_nonlocal"
    return "cds_and_exon", exon, None


def _changed_record(
    record: FeatureRecord,
    *,
    start0: int,
    end0: int,
) -> FeatureRecord:
    return replace(record, start=start0 + 1, end=end0)


def _reconcile_three_prime_utr(
    record: FeatureRecord,
    *,
    strand: str,
    old_terminal: tuple[int, int],
    new_terminal: tuple[int, int],
    offset: int,
) -> Optional[FeatureRecord]:
    """Reconcile only a source biological 3' UTR affected by the CDS move."""
    if record.feature_type.lower() != "three_prime_utr":
        return record

    if strand == "+":
        old_boundary0 = old_terminal[1]
        new_boundary0 = new_terminal[1]
        if offset > 0:
            # Clip only when newly coding bases actually intersect this UTR.
            if record.start0 >= new_boundary0 or record.end0 <= old_boundary0:
                return record
            new_start0 = max(record.start0, new_boundary0)
        else:
            # A trim contributes exon bases to an explicitly adjacent source UTR.
            if record.start0 != old_boundary0:
                return record
            new_start0 = new_boundary0
        if new_start0 >= record.end0:
            return None
        return _changed_record(record, start0=new_start0, end0=record.end0)

    old_boundary0 = old_terminal[0]
    new_boundary0 = new_terminal[0]
    if offset > 0:
        if record.end0 <= new_boundary0 or record.start0 >= old_boundary0:
            return record
        new_end0 = min(record.end0, new_boundary0)
    else:
        if record.end0 != old_boundary0:
            return record
        new_end0 = new_boundary0
    if new_end0 <= record.start0:
        return None
    return _changed_record(record, start0=record.start0, end0=new_end0)


def _apply_adjustment(
    model: CodingModel,
    offset: int,
    *,
    extension_mode: Optional[str],
    terminal_exon: Optional[FeatureRecord],
) -> CodingModel:
    old_terminal = _ordered_parts(model)[-1]
    if model.strand == "+":
        new_terminal = (old_terminal[0], old_terminal[1] + offset)
    else:
        # Positive offsets extend toward smaller genomic coordinates; negative
        # offsets trim from that same biological 3' edge.
        new_terminal = (old_terminal[0] - offset, old_terminal[1])
    if new_terminal[1] <= new_terminal[0]:
        raise ValueError("terminal CDS adjustment would remove the terminal CDS part")

    replacement_by_geometry = {old_terminal: new_terminal}
    new_parts = sorted(replacement_by_geometry.get(part, part) for part in model.cds_parts)
    new_children: list[FeatureRecord] = []
    for record in model.child_records:
        geometry = (record.start0, record.end0)
        if record.feature_type == "CDS" and geometry == old_terminal:
            new_children.append(
                _changed_record(record, start0=new_terminal[0], end0=new_terminal[1])
            )
        elif (
            extension_mode == "cds_and_exon"
            and terminal_exon is not None
            and record.feature_type.lower() == "exon"
            and geometry == (terminal_exon.start0, terminal_exon.end0)
        ):
            if model.strand == "+":
                new_children.append(
                    _changed_record(record, start0=record.start0, end0=record.end0 + offset)
                )
            else:
                new_children.append(
                    _changed_record(record, start0=record.start0 - offset, end0=record.end0)
                )
        else:
            reconciled = _reconcile_three_prime_utr(
                record,
                strand=model.strand,
                old_terminal=old_terminal,
                new_terminal=new_terminal,
                offset=offset,
            )
            if reconciled is not None:
                new_children.append(reconciled)
    return replace(model, cds_parts=new_parts, child_records=new_children)


def canonicalise_model_terminal_stop(
    model: CodingModel,
    fasta: FastaReader,
) -> CodingModel:
    original_boundary = _biological_boundary(model)
    phase, frame_basis, frame_error = _frame_basis(model, fasta)
    if frame_error is not None or phase is None:
        model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
            status="skipped", original_boundary=original_boundary,
            canonical_boundary=original_boundary, frame_confidence_basis=frame_basis,
            first_cds_phase=phase, skip_reason=frame_error or "uncertain_frame",
        )
        return model

    translated_length = model.cds_length - phase
    if translated_length <= 0:
        model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
            status="skipped", original_boundary=original_boundary,
            canonical_boundary=original_boundary, frame_confidence_basis=frame_basis,
            first_cds_phase=phase, skip_reason="uncertain_frame",
        )
        return model
    remainder = translated_length % 3
    suffix_needed = min(translated_length, 3 + remainder)
    try:
        suffix = _spliced_slice(
            fasta, model, model.cds_length - suffix_needed, model.cds_length
        )
    except (KeyError, ValueError):
        model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
            status="skipped", original_boundary=original_boundary,
            canonical_boundary=original_boundary, frame_confidence_basis=frame_basis,
            first_cds_phase=phase, skip_reason="sequence_unavailable",
        )
        return model

    translated_suffix = suffix
    if phase and model.cds_length == suffix_needed:
        translated_suffix = suffix[phase:]
    if remainder == 0 and len(translated_suffix) >= 3:
        included = translated_suffix[-3:]
        if included in STANDARD_STOP_CODONS:
            model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
                status="already_canonical", original_boundary=original_boundary,
                canonical_boundary=original_boundary, stop_codon=included,
                frame_confidence_basis=frame_basis, first_cds_phase=phase,
            )
            return model

    extension_count = 3 if remainder == 0 else 3 - remainder
    adjacent_start0, adjacent_end0 = _adjacent_interval(model, extension_count)
    contig_length = fasta.contig_length(model.chrom)
    if (
        contig_length is None or adjacent_start0 < 0
        or adjacent_end0 > contig_length
    ):
        model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
            status="skipped", original_boundary=original_boundary,
            canonical_boundary=original_boundary, frame_confidence_basis=frame_basis,
            first_cds_phase=phase, skip_reason="contig_boundary",
        )
        return model
    try:
        adjacent = _oriented_fetch(
            fasta, model, adjacent_start0, adjacent_end0
        )
    except (KeyError, ValueError):
        model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
            status="skipped", original_boundary=original_boundary,
            canonical_boundary=original_boundary, frame_confidence_basis=frame_basis,
            first_cds_phase=phase, skip_reason="sequence_unavailable",
        )
        return model

    partial = translated_suffix[-remainder:] if remainder else ""
    completion_codon = partial + adjacent
    if any(base not in "ACGT" for base in completion_codon):
        model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
            status="skipped", original_boundary=original_boundary,
            canonical_boundary=original_boundary, frame_confidence_basis=frame_basis,
            first_cds_phase=phase, skip_reason="ambiguous_genomic_bases",
        )
        return model
    completion_valid = completion_codon in STANDARD_STOP_CODONS

    trim_valid = False
    trim_codon = "."
    if remainder and translated_length - remainder >= 3:
        trim_codon = translated_suffix[-(remainder + 3):-remainder]
        trim_valid = trim_codon in STANDARD_STOP_CODONS
    if completion_valid and trim_valid:
        model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
            status="skipped", original_boundary=original_boundary,
            canonical_boundary=original_boundary, frame_confidence_basis=frame_basis,
            first_cds_phase=phase,
            skip_reason="ambiguous_partial_completion_vs_overrun_trim",
        )
        return model

    if completion_valid:
        offset = extension_count
        extension_mode, terminal_exon, structural_error = _extension_mode(model, offset)
        if structural_error is not None:
            model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
                status="skipped", original_boundary=original_boundary,
                canonical_boundary=original_boundary, stop_codon=completion_codon,
                frame_confidence_basis=frame_basis, first_cds_phase=phase,
                skip_reason=structural_error,
            )
            return model
        new_model = _apply_adjustment(
            model, offset, extension_mode=extension_mode, terminal_exon=terminal_exon
        )
        new_model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
            status="canonicalised", offset_bp=offset,
            original_boundary=original_boundary,
            canonical_boundary=_biological_boundary(new_model),
            stop_codon=completion_codon, frame_confidence_basis=frame_basis,
            first_cds_phase=phase, exon_adjusted=extension_mode == "cds_and_exon",
        )
        return new_model

    if trim_valid:
        offset = -remainder
        terminal_part = _ordered_parts(model)[-1]
        if terminal_part[1] - terminal_part[0] <= remainder:
            model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
                status="skipped", original_boundary=original_boundary,
                canonical_boundary=original_boundary, stop_codon=trim_codon,
                frame_confidence_basis=frame_basis, first_cds_phase=phase,
                skip_reason="would_remove_terminal_cds_part",
            )
            return model
        new_model = _apply_adjustment(
            model, offset, extension_mode=None, terminal_exon=None
        )
        new_model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
            status="canonicalised", offset_bp=offset,
            original_boundary=original_boundary,
            canonical_boundary=_biological_boundary(new_model),
            stop_codon=trim_codon, frame_confidence_basis=frame_basis,
            first_cds_phase=phase,
        )
        return new_model

    model.terminal_stop_canonicalisation = TerminalStopCanonicalisation(
        status="skipped", original_boundary=original_boundary,
        canonical_boundary=original_boundary, frame_confidence_basis=frame_basis,
        first_cds_phase=phase, skip_reason="no_qualifying_local_stop",
    )
    return model


def _interval_overlap(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return min(left[1], right[1]) > max(left[0], right[0])


def _original_parts(model: CodingModel) -> list[tuple[int, int]]:
    audit = model.terminal_stop_canonicalisation
    if audit is None or audit.status != "canonicalised":
        return list(model.cds_parts)
    offset = audit.offset_bp
    terminal = _ordered_parts(model)[-1]
    if model.strand == "+":
        old_terminal = (terminal[0], terminal[1] - offset)
    else:
        old_terminal = (terminal[0] + offset, terminal[1])
    return sorted(old_terminal if part == terminal else part for part in model.cds_parts)


def _parts_overlap(left: Iterable[tuple[int, int]], right: Iterable[tuple[int, int]]) -> bool:
    return any(_interval_overlap(a, b) for a in left for b in right)


def annotate_new_same_strand_overlaps(models: list[CodingModel]) -> None:
    """Audit newly created cross-gene candidate overlaps without vetoing them."""
    overlap_index = IntervalOverlapIndex(
        models,
        key=lambda item: (item.chrom, item.strand),
        start=lambda item: item.cds_start0,
        end=lambda item: item.cds_end0,
    )
    for model in models:
        audit = model.terminal_stop_canonicalisation
        if audit is None or audit.status != "canonicalised":
            continue
        overlaps: list[str] = []
        for other in overlap_index.query(
            (model.chrom, model.strand), model.cds_start0, model.cds_end0
        ):
            if other is model:
                continue
            if (other.source_name, other.input_gene_id) == (model.source_name, model.input_gene_id):
                continue
            if not _parts_overlap(model.cds_parts, other.cds_parts):
                continue
            if _parts_overlap(_original_parts(model), _original_parts(other)):
                continue
            overlaps.append(other.output_transcript_id)
        unique = sorted(set(overlaps))
        audit.newly_introduced_same_strand_cds_overlap_count = len(unique)
        audit.newly_introduced_same_strand_cds_overlap_models = (
            ",".join(unique) if unique else "."
        )


def canonicalise_parsed_inputs(
    parsed_inputs: list[ParsedInput],
    fasta: FastaReader,
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for parsed in parsed_inputs:
        canonical_models: list[CodingModel] = []
        for model in parsed.coding_models:
            canonical = canonicalise_model_terminal_stop(model, fasta)
            canonical_models.append(canonical)
            audit = canonical.terminal_stop_canonicalisation
            if audit is not None:
                if audit.status == "canonicalised":
                    counts[f"canonicalised_{audit.offset_bp:+d}"] += 1
                elif audit.status == "skipped":
                    counts[f"skipped:{audit.skip_reason}"] += 1
                else:
                    counts[audit.status] += 1
        parsed.coding_models = sorted(
            canonical_models,
            key=lambda item: (item.chrom, item.cds_start0, item.cds_end0, item.input_transcript_id),
        )
    models = [model for parsed in parsed_inputs for model in parsed.coding_models]
    annotate_new_same_strand_overlaps(models)
    return counts


def terminal_stop_audit_row(model: CodingModel) -> dict[str, object]:
    audit = model.terminal_stop_canonicalisation
    return {
        "source_name": model.source_name,
        "input_file": str(model.input_file),
        "input_gene_id": model.input_gene_id,
        "input_transcript_id": model.input_transcript_id,
        "output_transcript_id": model.output_transcript_id,
        "chrom": model.chrom,
        "strand": model.strand,
        "status": audit.status if audit else "not_evaluated",
        "offset": f"{audit.offset_bp:+d}" if audit and audit.offset_bp else ".",
        "original_biological_3prime_cds_boundary": audit.original_boundary if audit else ".",
        "canonical_biological_3prime_cds_boundary": audit.canonical_boundary if audit else ".",
        "stop_codon": audit.stop_codon if audit else ".",
        "frame_confidence_basis": audit.frame_confidence_basis if audit else ".",
        "first_cds_phase": audit.first_cds_phase if audit and audit.first_cds_phase is not None else ".",
        "skip_reason": audit.skip_reason if audit else "genome_fasta_not_supplied",
        "terminal_exon_adjusted": audit.exon_adjusted if audit else False,
        "new_same_strand_cds_overlap_count": audit.newly_introduced_same_strand_cds_overlap_count if audit else 0,
        "new_same_strand_cds_overlap_models": audit.newly_introduced_same_strand_cds_overlap_models if audit else ".",
    }


def write_terminal_stop_audit_tsv(models: list[CodingModel], path: Path) -> None:
    # Stream because this audit has one row per prediction model and can be
    # large in multi-source production runs.
    with path.open("w") as out:
        out.write("\t".join(AUDIT_COLUMNS) + "\n")
        for model in models:
            row = terminal_stop_audit_row(model)
            out.write(
                "\t".join(
                    format_tsv_value(row.get(column, "."))
                    for column in AUDIT_COLUMNS
                )
                + "\n"
            )


def write_terminal_stop_actions_tsv(models: list[CodingModel], path: Path) -> None:
    """Write only models whose terminal-stop representation was modified."""
    with path.open("w") as out:
        out.write("\t".join(AUDIT_COLUMNS) + "\n")
        for model in models:
            row = terminal_stop_audit_row(model)
            if row["status"] != "canonicalised":
                continue
            out.write(
                "\t".join(
                    format_tsv_value(row.get(column, "."))
                    for column in AUDIT_COLUMNS
                )
                + "\n"
            )


def write_terminal_stop_summary_tsv(
    models: list[CodingModel], path: Path, *, execution_status: str
) -> None:
    statuses: Counter[str] = Counter()
    offsets: Counter[str] = Counter()
    skip_reasons: Counter[str] = Counter()
    exon_adjusted = 0
    new_overlap = 0
    for model in models:
        audit = model.terminal_stop_canonicalisation
        if audit is None:
            statuses["not_evaluated"] += 1
            continue
        statuses[audit.status] += 1
        if audit.status == "canonicalised":
            offsets[f"{audit.offset_bp:+d}"] += 1
            exon_adjusted += int(audit.exon_adjusted)
            new_overlap += int(audit.newly_introduced_same_strand_cds_overlap_count > 0)
        elif audit.status == "skipped":
            skip_reasons[audit.skip_reason] += 1
    row = {
        "execution_status": execution_status,
        "total_models": len(models),
        "canonicalised_models": statuses["canonicalised"],
        "already_canonical_models": statuses["already_canonical"],
        "skipped_models": statuses["skipped"],
        "canonicalised_offset_counts": json.dumps(dict(sorted(offsets.items())), separators=(",", ":")),
        "skip_reason_counts": json.dumps(dict(sorted(skip_reasons.items())), separators=(",", ":")),
        "canonicalised_with_terminal_exon_adjustment": exon_adjusted,
        "canonicalised_with_new_same_strand_overlap": new_overlap,
    }
    write_dict_tsv([row], path, list(row))


def compact_terminal_stop_provenance(model: CodingModel) -> Optional[str]:
    audit = model.terminal_stop_canonicalisation
    if audit is None or audit.status != "canonicalised":
        return None
    return (
        f"offset:{audit.offset_bp:+d}|original_boundary:{audit.original_boundary}|"
        f"canonical_boundary:{audit.canonical_boundary}|stop:{audit.stop_codon}|"
        f"frame:{audit.frame_confidence_basis}"
    )
