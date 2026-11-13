"""Conservative pre-integrity recovery of missing biological CDS phase."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import replace
from pathlib import Path
from typing import Optional, Protocol

from .cds_phase import (
    CANONICALISED_FIVE_PRIME_BOUNDARY,
    CDSFrameEvidence,
    RECOVERED_MISSING_SOURCE_PHASE,
    analyse_cds_frame,
    biological_cds_parts,
    calculate_cds_phase_chain,
    source_first_cds_phase,
)
from .genome_sequence import reverse_complement
from .models import CDSPhaseRecovery, CodingModel, ParsedInput
from .tsvio import format_tsv_value, write_dict_tsv


class FastaReader(Protocol):
    def contig_length(self, name: str) -> Optional[int]: ...
    def fetch(self, name: str, start0: int, end0: int) -> str: ...


AUDIT_COLUMNS = [
    "source_name", "input_file", "input_gene_id", "input_transcript_id",
    "output_transcript_id", "chrom", "strand", "status", "reason",
    "original_source_first_phase", "original_effective_phase",
    "recovered_first_phase", "evidence_basis", "stop_representation",
    "stop_codon", "original_stop_status", "recovered_stop_status",
    "original_internal_stop_count", "recovered_internal_stop_count",
]


def _oriented_fetch(
    fasta: FastaReader,
    model: CodingModel,
    start0: int,
    end0: int,
) -> str:
    sequence = fasta.fetch(model.chrom, start0, end0)
    return reverse_complement(sequence) if model.strand == "-" else sequence


def model_sequence_context(
    model: CodingModel,
    fasta: FastaReader,
) -> tuple[str, str]:
    """Fetch the spliced CDS and biological 3-prime adjacent codon once."""
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
        adjacent = _oriented_fetch(
            fasta, model, adjacent_start0, adjacent_end0
        )
    return sequence, adjacent


def _recovery_not_applicable(model: CodingModel) -> bool:
    boundary = model.five_prime_boundary_canonicalisation
    return bool(
        (boundary is not None and boundary.status == "canonicalised")
        or source_first_cds_phase(model) is not None
    )


def _stop_status(frame: CDSFrameEvidence) -> str:
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
    reason: str,
    source_phase: int | None = None,
    recovered_phase: int | None = None,
    original: CDSFrameEvidence | None = None,
    recovered: CDSFrameEvidence | None = None,
) -> CDSPhaseRecovery:
    return CDSPhaseRecovery(
        status=status,
        original_source_first_phase=source_phase,
        recovered_first_phase=recovered_phase,
        evidence_basis=(
            "biological_3prime_terminal_stop"
            if recovered is not None and recovered.definite_stop
            else "."
        ),
        stop_representation=(recovered.stop_representation if recovered else "."),
        stop_codon=(recovered.terminal_codon if recovered else "."),
        original_stop_status=_stop_status(original) if original else ".",
        recovered_stop_status=_stop_status(recovered) if recovered else ".",
        original_internal_stop_count=(
            original.definite_internal_stop_count if original else 0
        ),
        recovered_internal_stop_count=(
            recovered.definite_internal_stop_count if recovered else 0
        ),
        reason=reason,
    )


def choose_recovered_phase(
    frames: dict[int, CDSFrameEvidence],
) -> tuple[int | None, str]:
    """Apply the sequence-only recovery policy to precomputed frame evidence."""
    phase_zero = frames[0]
    if phase_zero.definite_start or phase_zero.definite_stop:
        return None, "phase_zero_positively_supported"
    candidates = [phase for phase in (1, 2) if frames[phase].definite_stop]
    if len(candidates) != 1:
        return None, "no_unique_nonzero_stop_supported_frame"
    phase = candidates[0]
    candidate = frames[phase]
    if candidate.definite_internal_stop_count:
        return None, "recovered_frame_internal_stop"
    if candidate.definite_start:
        return None, "five_prime_leading_base_offset"
    return phase, "recovered"


def recover_model_cds_phase_from_context(
    model: CodingModel,
    sequence: str | None,
    adjacent: str = "",
    *,
    context_available: bool = True,
) -> CodingModel:
    """Recover phase using a sequence context already fetched by the caller."""
    boundary = model.five_prime_boundary_canonicalisation
    if boundary is not None and boundary.status == "canonicalised":
        # Phase zero was established by a geometry correction, not recovered
        # from the source's missing phase.
        return model
    source_phase = source_first_cds_phase(model)
    if source_phase is not None:
        # Numeric-source models need no per-model recovery object: the audit
        # writer can derive this immutable non-eligibility directly.  This
        # keeps the common numeric-source path memory-neutral.
        return model
    if not model.cds_parts:
        model.cds_phase_recovery = _audit(
            "skipped", reason="no_cds_parts"
        )
        return model
    if not context_available or sequence is None:
        model.cds_phase_recovery = _audit(
            "skipped", reason="sequence_context_unavailable"
        )
        return model
    if not sequence:
        model.cds_phase_recovery = _audit(
            "skipped", reason="empty_cds_sequence"
        )
        return model
    # Early biological recovery deliberately remains stricter than the later
    # output-only phase-zero standardisation.
    if any(base not in "ACGT" for base in sequence):
        model.cds_phase_recovery = _audit(
            "skipped", reason="ambiguous_terminal_or_frame_evidence"
        )
        return model
    frames = {
        phase: analyse_cds_frame(
            sequence, adjacent, phase, cds_is_acgt=True
        )
        for phase in (0, 1, 2)
    }
    if any(
        frame.possible_start != frame.definite_start
        or frame.possible_stop != frame.definite_stop
        for frame in frames.values()
    ):
        model.cds_phase_recovery = _audit(
            "skipped",
            reason="ambiguous_terminal_or_frame_evidence",
            original=frames[0],
        )
        return model
    phase, reason = choose_recovered_phase(frames)
    original = frames[0]
    if phase is None:
        model.cds_phase_recovery = _audit(
            "skipped", reason=reason, original=original
        )
        return model
    candidate = frames[phase]
    # A terminal stop used for biological recovery must be a literal,
    # unambiguous DNA codon even though the shared analyzer understands IUPAC.
    if (
        candidate.terminal_codon not in {"TAA", "TAG", "TGA"}
        or candidate.stop_representation not in {"included", "adjacent"}
    ):
        model.cds_phase_recovery = _audit(
            "skipped",
            reason="ambiguous_terminal_or_frame_evidence",
            original=original,
            recovered=candidate,
        )
        return model
    phases = calculate_cds_phase_chain(model.cds_parts, model.strand, phase)
    children = [
        replace(record, phase=phases[(record.start0, record.end0)])
        if (
            record.feature_type.lower() == "cds"
            and (record.start0, record.end0) in phases
        )
        else record
        for record in model.child_records
    ]
    recovery = _audit(
        RECOVERED_MISSING_SOURCE_PHASE,
        reason=f"recovered_phase_{phase}",
        recovered_phase=phase,
        original=original,
        recovered=candidate,
    )
    return replace(
        model,
        child_records=children,
        cds_phase_recovery=recovery,
    )


def recover_model_cds_phase(
    model: CodingModel,
    fasta: FastaReader,
) -> CodingModel:
    """Return a model with a uniquely stop-supported missing phase recovered."""
    if _recovery_not_applicable(model):
        return model
    if not model.cds_parts:
        return recover_model_cds_phase_from_context(model, None)
    try:
        sequence, adjacent = model_sequence_context(model, fasta)
    except (KeyError, ValueError, IndexError):
        return recover_model_cds_phase_from_context(
            model, None, context_available=False
        )
    return recover_model_cds_phase_from_context(model, sequence, adjacent)


def phase_recovery_count_key(model: CodingModel) -> str:
    """Return the stable summary bucket for one phase-recovery result."""
    audit = model.cds_phase_recovery
    key = (
        audit.status
        if audit is not None
        else CANONICALISED_FIVE_PRIME_BOUNDARY
        if (
            model.five_prime_boundary_canonicalisation is not None
            and model.five_prime_boundary_canonicalisation.status == "canonicalised"
        )
        else "not_eligible_numeric_source_phase"
        if source_first_cds_phase(model) is not None
        else "not_evaluated"
    )
    if audit is not None and audit.status == RECOVERED_MISSING_SOURCE_PHASE:
        return f"recovered_phase_{audit.recovered_first_phase}"
    if audit is not None and audit.status == "skipped":
        return f"skipped:{audit.reason}"
    return key


def recover_parsed_input_phases(
    parsed_inputs: list[ParsedInput],
    fasta: FastaReader,
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for parsed in parsed_inputs:
        recovered_models: list[CodingModel] = []
        for model in parsed.coding_models:
            recovered = recover_model_cds_phase(model, fasta)
            recovered_models.append(recovered)
            counts[phase_recovery_count_key(recovered)] += 1
        parsed.coding_models = recovered_models
    return counts


def phase_recovery_audit_row(model: CodingModel) -> dict[str, object]:
    audit = model.cds_phase_recovery
    source_phase = source_first_cds_phase(model)
    boundary = model.five_prime_boundary_canonicalisation
    boundary_canonicalised = boundary is not None and boundary.status == "canonicalised"
    derived_numeric = audit is None and source_phase is not None and not boundary_canonicalised
    return {
        "source_name": model.source_name,
        "input_file": str(model.input_file),
        "input_gene_id": model.input_gene_id,
        "input_transcript_id": model.input_transcript_id,
        "output_transcript_id": model.output_transcript_id,
        "chrom": model.chrom,
        "strand": model.strand,
        "status": (
            audit.status if audit else
            CANONICALISED_FIVE_PRIME_BOUNDARY if boundary_canonicalised else
            "not_eligible_numeric_source_phase" if derived_numeric else
            "not_evaluated"
        ),
        "reason": (
            audit.reason if audit else
            "five_prime_boundary_canonicalised" if boundary_canonicalised else
            "numeric_source_phase" if derived_numeric else
            "genome_fasta_not_supplied"
        ),
        "original_source_first_phase": (
            audit.original_source_first_phase
            if audit and audit.original_source_first_phase is not None else
            "." if boundary_canonicalised else
            source_phase if derived_numeric else "."
        ),
        "original_effective_phase": (
            audit.original_effective_phase if audit else
            0 if boundary_canonicalised else
            source_phase if derived_numeric else "."
        ),
        "recovered_first_phase": (
            audit.recovered_first_phase
            if audit and audit.recovered_first_phase is not None else "."
        ),
        "evidence_basis": audit.evidence_basis if audit else ".",
        "stop_representation": audit.stop_representation if audit else ".",
        "stop_codon": audit.stop_codon if audit else ".",
        "original_stop_status": audit.original_stop_status if audit else ".",
        "recovered_stop_status": audit.recovered_stop_status if audit else ".",
        "original_internal_stop_count": audit.original_internal_stop_count if audit else 0,
        "recovered_internal_stop_count": audit.recovered_internal_stop_count if audit else 0,
    }


def write_phase_recovery_audit_tsv(
    models: list[CodingModel],
    path: Path,
) -> None:
    with path.open("w") as out:
        out.write("\t".join(AUDIT_COLUMNS) + "\n")
        for model in models:
            row = phase_recovery_audit_row(model)
            out.write(
                "\t".join(
                    format_tsv_value(row.get(column, "."))
                    for column in AUDIT_COLUMNS
                )
                + "\n"
            )


def write_phase_recovery_actions_tsv(
    models: list[CodingModel], path: Path
) -> None:
    """Write only models whose missing biological first CDS phase was recovered."""
    with path.open("w") as out:
        out.write("\t".join(AUDIT_COLUMNS) + "\n")
        for model in models:
            row = phase_recovery_audit_row(model)
            if row["status"] != RECOVERED_MISSING_SOURCE_PHASE:
                continue
            out.write(
                "\t".join(
                    format_tsv_value(row.get(column, "."))
                    for column in AUDIT_COLUMNS
                )
                + "\n"
            )


def write_phase_recovery_summary_tsv(
    models: list[CodingModel], path: Path, *, execution_status: str
) -> None:
    statuses: Counter[str] = Counter()
    recovered_phases: Counter[str] = Counter()
    skip_reasons: Counter[str] = Counter()
    for model in models:
        row = phase_recovery_audit_row(model)
        status = str(row["status"])
        statuses[status] += 1
        if status == RECOVERED_MISSING_SOURCE_PHASE:
            recovered_phases[str(row["recovered_first_phase"])] += 1
        elif status == "skipped":
            skip_reasons[str(row["reason"])] += 1
    row = {
        "execution_status": execution_status,
        "total_models": len(models),
        "recovered_models": statuses[RECOVERED_MISSING_SOURCE_PHASE],
        "numeric_source_phase_models": statuses["not_eligible_numeric_source_phase"],
        "five_prime_boundary_canonicalised_models": statuses[CANONICALISED_FIVE_PRIME_BOUNDARY],
        "skipped_models": statuses["skipped"],
        "not_evaluated_models": statuses["not_evaluated"],
        "recovered_phase_counts": json.dumps(dict(sorted(recovered_phases.items())), separators=(",", ":")),
        "skip_reason_counts": json.dumps(dict(sorted(skip_reasons.items())), separators=(",", ":")),
    }
    write_dict_tsv([row], path, list(row))
