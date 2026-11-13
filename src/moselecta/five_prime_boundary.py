"""Bounded class-B canonicalisation of 1--2 bp 5' CDS boundary errors."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import replace
from pathlib import Path
from typing import Optional, Protocol

from .cds_phase import (
    CDSFrameEvidence,
    STANDARD_STOP_CODONS,
    analyse_cds_frame,
    biological_cds_parts,
    calculate_cds_phase_chain,
    parse_cds_phase,
    source_phase_by_part,
)
from .genome_sequence import reverse_complement
from .models import (
    CodingModel,
    FeatureRecord,
    FivePrimeBoundaryCanonicalisation,
    ParsedInput,
)
from .tsvio import format_tsv_value, write_dict_tsv


class FastaReader(Protocol):
    def contig_length(self, name: str) -> Optional[int]: ...
    def fetch(self, name: str, start0: int, end0: int) -> str: ...


AUDIT_COLUMNS = [
    "source_name", "input_file", "input_gene_id", "input_transcript_id",
    "output_transcript_id", "chrom", "strand", "status", "skip_reason",
    "bases_trimmed", "original_biological_5prime_cds_boundary",
    "canonical_biological_5prime_cds_boundary", "exon_context",
    "five_prime_utr_extended", "original_source_first_phase",
    "canonical_first_phase", "start_codon", "stop_codon",
    "stop_representation", "original_start_status", "canonical_start_status",
    "original_stop_status", "canonical_stop_status",
    "original_internal_stop_count", "canonical_internal_stop_count",
]


def _oriented_fetch(
    fasta: FastaReader, model: CodingModel, start0: int, end0: int
) -> str:
    sequence = fasta.fetch(model.chrom, start0, end0).upper()
    return reverse_complement(sequence) if model.strand == "-" else sequence


def _model_sequence(model: CodingModel, fasta: FastaReader) -> tuple[str, str]:
    sequence = "".join(
        _oriented_fetch(fasta, model, start0, end0)
        for start0, end0 in biological_cds_parts(model)
    )
    contig_length = fasta.contig_length(model.chrom)
    if contig_length is None:
        raise KeyError(model.chrom)
    if model.strand == "+":
        adjacent_start0, adjacent_end0 = model.cds_end0, model.cds_end0 + 3
    else:
        adjacent_start0, adjacent_end0 = model.cds_start0 - 3, model.cds_start0
    adjacent = ""
    if adjacent_start0 >= 0 and adjacent_end0 <= contig_length:
        adjacent = _oriented_fetch(fasta, model, adjacent_start0, adjacent_end0)
    return sequence, adjacent


def _boundary(model: CodingModel) -> int:
    first = biological_cds_parts(model)[0]
    return first[0] + 1 if model.strand == "+" else first[1]


def _stop_status(frame: CDSFrameEvidence | None) -> str:
    if frame is None:
        return "."
    if frame.stop_representation == "included":
        return "included_in_cds"
    if frame.stop_representation == "adjacent":
        return "adjacent_to_cds"
    if frame.stop_representation.startswith("ambiguous_"):
        return "ambiguous"
    return "missing"


def _audit(
    status: str,
    *,
    reason: str = ".",
    source_phase: str = ".",
    original_boundary: int | None = None,
    canonical_boundary: int | None = None,
    exon_context: str = ".",
    trimmed: int = 0,
    utr_extended: bool = False,
    original: CDSFrameEvidence | None = None,
    canonical: CDSFrameEvidence | None = None,
) -> FivePrimeBoundaryCanonicalisation:
    return FivePrimeBoundaryCanonicalisation(
        status=status,
        skip_reason=reason,
        bases_trimmed=trimmed,
        original_boundary=original_boundary,
        canonical_boundary=canonical_boundary if canonical_boundary is not None else original_boundary,
        exon_context=exon_context,
        five_prime_utr_extended=utr_extended,
        original_source_first_phase=source_phase,
        canonical_first_phase=0 if status == "canonicalised" else None,
        start_codon="ATG" if canonical is not None and canonical.definite_start else ".",
        stop_codon=(canonical or original).terminal_codon if (canonical or original) else ".",
        stop_representation=(canonical or original).stop_representation if (canonical or original) else ".",
        original_start_status="present" if original and original.definite_start else "missing",
        canonical_start_status="present" if canonical and canonical.definite_start else ".",
        original_stop_status=_stop_status(original),
        canonical_stop_status=_stop_status(canonical),
        original_internal_stop_count=original.definite_internal_stop_count if original else 0,
        canonical_internal_stop_count=canonical.definite_internal_stop_count if canonical else 0,
    )


def _source_exon_context(
    model: CodingModel, first_part: tuple[int, int]
) -> tuple[FeatureRecord | None, str, str | None]:
    containing = [
        record for record in model.child_records
        if record.feature_type.lower() == "exon"
        and record.seqid == model.chrom
        and record.strand == model.strand
        and record.start0 <= first_part[0]
        and record.end0 >= first_part[1]
    ]
    geometries = {(record.start0, record.end0) for record in containing}
    if not containing:
        return None, "missing_source_exon", "missing_source_exon"
    if len(geometries) != 1:
        return None, "ambiguous_source_exon", "ambiguous_source_exon"
    exon = containing[0]
    strictly_inside = (
        exon.start0 < first_part[0]
        if model.strand == "+"
        else first_part[1] < exon.end0
    )
    if not strictly_inside:
        return exon, "cds_at_source_exon_boundary", "cds_starts_at_exon_boundary"
    return exon, "strictly_inside_source_exon", None


def _changed_record(
    record: FeatureRecord, *, start0: int, end0: int, phase: str | None = None
) -> FeatureRecord:
    return replace(
        record, start=start0 + 1, end=end0,
        phase=record.phase if phase is None else phase,
    )


def _apply_trim(model: CodingModel, offset: int) -> tuple[CodingModel, bool]:
    first = biological_cds_parts(model)[0]
    new_first = (
        (first[0] + offset, first[1])
        if model.strand == "+"
        else (first[0], first[1] - offset)
    )
    replacements = {first: new_first}
    new_parts = sorted(replacements.get(part, part) for part in model.cds_parts)
    phases = calculate_cds_phase_chain(new_parts, model.strand, 0)
    utr_extended = False
    children: list[FeatureRecord] = []
    for record in model.child_records:
        geometry = (record.start0, record.end0)
        if record.feature_type.lower() == "cds" and geometry in replacements:
            replacement = replacements[geometry]
            children.append(
                _changed_record(record, start0=replacement[0], end0=replacement[1], phase=phases[replacement])
            )
        elif record.feature_type.lower() == "cds" and geometry in phases:
            children.append(
                _changed_record(record, start0=record.start0, end0=record.end0, phase=phases[geometry])
            )
        elif record.feature_type.lower() == "five_prime_utr":
            if model.strand == "+" and record.end0 == first[0]:
                children.append(_changed_record(record, start0=record.start0, end0=new_first[0]))
                utr_extended = True
            elif model.strand == "-" and record.start0 == first[1]:
                children.append(_changed_record(record, start0=new_first[1], end0=record.end0))
                utr_extended = True
            else:
                children.append(record)
        else:
            children.append(record)
    return replace(model, cds_parts=new_parts, child_records=children), utr_extended


def _validate_canonical_model(
    model: CodingModel,
    fasta: FastaReader,
    *,
    expected_stop: str,
    expected_representation: str,
) -> tuple[CDSFrameEvidence | None, str | None]:
    first = biological_cds_parts(model)[0]
    if first[1] <= first[0]:
        return None, "postcheck_empty_first_cds_part"
    try:
        sequence, adjacent = _model_sequence(model, fasta)
    except (KeyError, ValueError, IndexError):
        return None, "postcheck_sequence_unavailable"
    frame = analyse_cds_frame(sequence, adjacent, 0, cds_is_acgt=True)
    phases = calculate_cds_phase_chain(model.cds_parts, model.strand, 0)
    observed = source_phase_by_part(model)
    if not frame.definite_start or sequence[:3] != "ATG":
        return frame, "postcheck_start_not_atg"
    if len(sequence) % 3 != 0 and frame.stop_representation == "included":
        return frame, "postcheck_incoherent_cds_length"
    if any(observed.get(part) != phase for part, phase in phases.items()):
        return frame, "postcheck_phase_chain_inconsistent"
    if frame.definite_internal_stop_count:
        return frame, "postcheck_internal_stop"
    if (
        not frame.definite_stop
        or frame.terminal_codon != expected_stop
        or frame.stop_representation != expected_representation
    ):
        return frame, "postcheck_terminal_stop_changed"
    return frame, None


def canonicalise_model_five_prime_boundary(
    model: CodingModel, fasta: FastaReader
) -> CodingModel:
    if not model.cds_parts:
        model.five_prime_boundary_canonicalisation = _audit("skipped", reason="no_cds_parts")
        return model
    first = biological_cds_parts(model)[0]
    original_boundary = _boundary(model)
    source_phase = source_phase_by_part(model).get(first, ".")
    if source_phase != ".":
        numeric_source_phase = parse_cds_phase(source_phase) is not None
        model.five_prime_boundary_canonicalisation = _audit(
            "not_eligible_numeric_source_phase" if numeric_source_phase else "skipped",
            reason="numeric_source_phase" if numeric_source_phase else "source_first_phase_not_dot",
            source_phase=source_phase, original_boundary=original_boundary,
        )
        return model
    _exon, exon_context, exon_error = _source_exon_context(model, first)
    if exon_error is not None:
        model.five_prime_boundary_canonicalisation = _audit(
            "skipped", reason=exon_error, source_phase=source_phase,
            original_boundary=original_boundary, exon_context=exon_context,
        )
        return model
    try:
        sequence, adjacent = _model_sequence(model, fasta)
    except (KeyError, ValueError, IndexError):
        model.five_prime_boundary_canonicalisation = _audit(
            "skipped", reason="sequence_context_unavailable", source_phase=source_phase,
            original_boundary=original_boundary, exon_context=exon_context,
        )
        return model
    if not sequence or any(base not in "ACGT" for base in sequence):
        model.five_prime_boundary_canonicalisation = _audit(
            "skipped", reason="ambiguous_sequence" if sequence else "empty_cds_sequence",
            source_phase=source_phase, original_boundary=original_boundary,
            exon_context=exon_context,
        )
        return model
    frames = {
        phase: analyse_cds_frame(sequence, adjacent, phase, cds_is_acgt=True)
        for phase in (0, 1, 2)
    }
    original = frames[0]
    if original.definite_start or original.definite_stop:
        model.five_prime_boundary_canonicalisation = _audit(
            "skipped", reason="phase_zero_positively_supported", source_phase=source_phase,
            original_boundary=original_boundary, exon_context=exon_context, original=original,
        )
        return model
    stop_supported = [
        phase for phase in (1, 2)
        if frames[phase].definite_start
        and frames[phase].terminal_codon in STANDARD_STOP_CODONS
        and frames[phase].stop_representation in {"included", "adjacent"}
    ]
    clean = [phase for phase in stop_supported if not frames[phase].definite_internal_stop_count]
    if len(clean) != 1:
        reason = "offset_frame_internal_stop" if stop_supported and not clean else "no_unique_start_stop_offset"
        model.five_prime_boundary_canonicalisation = _audit(
            "skipped", reason=reason, source_phase=source_phase,
            original_boundary=original_boundary, exon_context=exon_context, original=original,
        )
        return model
    offset = clean[0]
    candidate = frames[offset]
    if first[1] - first[0] < offset + 3:
        model.five_prime_boundary_canonicalisation = _audit(
            "skipped", reason="start_codon_not_in_first_cds_part", source_phase=source_phase,
            original_boundary=original_boundary, exon_context=exon_context, original=original,
        )
        return model
    proposed, utr_extended = _apply_trim(model, offset)
    canonical, post_error = _validate_canonical_model(
        proposed, fasta, expected_stop=candidate.terminal_codon,
        expected_representation=candidate.stop_representation,
    )
    if post_error is not None:
        model.five_prime_boundary_canonicalisation = _audit(
            "skipped", reason=post_error, source_phase=source_phase,
            original_boundary=original_boundary, exon_context=exon_context, original=original,
        )
        return model
    proposed.five_prime_boundary_canonicalisation = _audit(
        "canonicalised", source_phase=source_phase, original_boundary=original_boundary,
        canonical_boundary=_boundary(proposed), exon_context=exon_context,
        trimmed=offset, utr_extended=utr_extended, original=original, canonical=canonical,
    )
    return proposed


def canonicalise_parsed_input_five_prime_boundaries(
    parsed_inputs: list[ParsedInput], fasta: FastaReader
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for parsed in parsed_inputs:
        corrected: list[CodingModel] = []
        for model in parsed.coding_models:
            result = canonicalise_model_five_prime_boundary(model, fasta)
            corrected.append(result)
            audit = result.five_prime_boundary_canonicalisation
            if audit is None:
                counts["not_evaluated"] += 1
            elif audit.status == "canonicalised":
                counts[f"canonicalised_trim_{audit.bases_trimmed}"] += 1
            elif audit.status == "skipped":
                counts[f"skipped:{audit.skip_reason}"] += 1
            else:
                counts[audit.status] += 1
        parsed.coding_models = sorted(
            corrected,
            key=lambda item: (item.chrom, item.cds_start0, item.cds_end0, item.input_transcript_id),
        )
    return counts


def five_prime_boundary_audit_row(model: CodingModel) -> dict[str, object]:
    audit = model.five_prime_boundary_canonicalisation
    return {
        "source_name": model.source_name,
        "input_file": str(model.input_file),
        "input_gene_id": model.input_gene_id,
        "input_transcript_id": model.input_transcript_id,
        "output_transcript_id": model.output_transcript_id,
        "chrom": model.chrom,
        "strand": model.strand,
        "status": audit.status if audit else "not_evaluated",
        "skip_reason": audit.skip_reason if audit else "genome_fasta_not_supplied",
        "bases_trimmed": audit.bases_trimmed if audit else 0,
        "original_biological_5prime_cds_boundary": audit.original_boundary if audit else ".",
        "canonical_biological_5prime_cds_boundary": audit.canonical_boundary if audit else ".",
        "exon_context": audit.exon_context if audit else ".",
        "five_prime_utr_extended": audit.five_prime_utr_extended if audit else False,
        "original_source_first_phase": audit.original_source_first_phase if audit else ".",
        "canonical_first_phase": audit.canonical_first_phase if audit and audit.canonical_first_phase is not None else ".",
        "start_codon": audit.start_codon if audit else ".",
        "stop_codon": audit.stop_codon if audit else ".",
        "stop_representation": audit.stop_representation if audit else ".",
        "original_start_status": audit.original_start_status if audit else ".",
        "canonical_start_status": audit.canonical_start_status if audit else ".",
        "original_stop_status": audit.original_stop_status if audit else ".",
        "canonical_stop_status": audit.canonical_stop_status if audit else ".",
        "original_internal_stop_count": audit.original_internal_stop_count if audit else 0,
        "canonical_internal_stop_count": audit.canonical_internal_stop_count if audit else 0,
    }


def write_five_prime_boundary_audit_tsv(models: list[CodingModel], path: Path) -> None:
    with path.open("w") as out:
        out.write("\t".join(AUDIT_COLUMNS) + "\n")
        for model in sorted(
            models,
            key=lambda item: (
                item.source_name, str(item.input_file), item.input_transcript_id,
                item.output_transcript_id,
            ),
        ):
            row = five_prime_boundary_audit_row(model)
            out.write("\t".join(format_tsv_value(row.get(column, ".")) for column in AUDIT_COLUMNS) + "\n")


def write_five_prime_boundary_actions_tsv(
    models: list[CodingModel], path: Path
) -> None:
    """Write only models whose biological 5-prime CDS boundary was modified."""
    with path.open("w") as out:
        out.write("\t".join(AUDIT_COLUMNS) + "\n")
        for model in sorted(
            models,
            key=lambda item: (
                item.source_name, str(item.input_file), item.input_transcript_id,
                item.output_transcript_id,
            ),
        ):
            row = five_prime_boundary_audit_row(model)
            if row["status"] != "canonicalised":
                continue
            out.write("\t".join(format_tsv_value(row.get(column, ".")) for column in AUDIT_COLUMNS) + "\n")


def write_five_prime_boundary_summary_tsv(
    models: list[CodingModel], path: Path, *, execution_status: str
) -> None:
    statuses: Counter[str] = Counter()
    trims: Counter[str] = Counter()
    skip_reasons: Counter[str] = Counter()
    utr_extensions = 0
    for model in models:
        audit = model.five_prime_boundary_canonicalisation
        if audit is None:
            statuses["not_evaluated"] += 1
            continue
        statuses[audit.status] += 1
        if audit.status == "canonicalised":
            trims[str(audit.bases_trimmed)] += 1
            utr_extensions += int(audit.five_prime_utr_extended)
        elif audit.status == "skipped":
            skip_reasons[audit.skip_reason] += 1
    row = {
        "execution_status": execution_status,
        "total_models": len(models),
        "canonicalised_models": statuses["canonicalised"],
        "skipped_models": statuses["skipped"],
        "not_evaluated_models": statuses["not_evaluated"],
        "canonicalised_bases_trimmed_counts": json.dumps(dict(sorted(trims.items())), separators=(",", ":")),
        "skip_reason_counts": json.dumps(dict(sorted(skip_reasons.items())), separators=(",", ":")),
        "canonicalised_with_utr_extension": utr_extensions,
    }
    write_dict_tsv([row], path, list(row))


def compact_five_prime_boundary_provenance(model: CodingModel) -> Optional[str]:
    audit = model.five_prime_boundary_canonicalisation
    if audit is None or audit.status != "canonicalised":
        return None
    return (
        f"trim:{audit.bases_trimmed}|original_boundary:{audit.original_boundary}|"
        f"canonical_boundary:{audit.canonical_boundary}|start:ATG|source_phase:."
    )
