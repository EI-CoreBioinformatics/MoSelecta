"""Early exclusion of models with frame-definite internal stop codons."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from .cds_phase import (
    CANONICALISED_FIVE_PRIME_BOUNDARY,
    PRESERVED_NUMERIC_SOURCE,
    RECOVERED_MISSING_SOURCE_PHASE,
    UNRESOLVED_SOURCE_MISSING,
    analyse_cds_frame,
    source_first_cds_phase,
)
from .models import CodingModel, ParsedInput
from .phase_recovery import (
    FastaReader,
    model_sequence_context,
    phase_recovery_count_key,
    recover_model_cds_phase_from_context,
)
from .tsvio import write_dict_tsv


RESOLVED_FRAME_INTERNAL_STOP = "resolved_frame_internal_stop"
UNRESOLVED_ALL_FRAMES_INTERNAL_STOP = "unresolved_all_frames_internal_stop"
INVALID_CDS_GEOMETRY = "invalid_cds_geometry"

LEDGER_COLUMNS = [
    "model_id", "source", "input_file", "original_gene_id",
    "original_transcript_id", "seqid", "strand", "frame_resolution_status",
    "phase_provenance", "active_first_cds_phase", "exclusion_class",
    "frame0_internal_stop_count", "frame1_internal_stop_count",
    "frame2_internal_stop_count", "has_ambiguous_bases",
    "sequence_evaluability",
]


@dataclass(frozen=True, slots=True)
class DefiniteInternalStopExclusion:
    model_id: str
    source: str
    input_file: str
    original_gene_id: str
    original_transcript_id: str
    seqid: str
    strand: str
    frame_resolution_status: str
    phase_provenance: str
    active_first_cds_phase: int | None
    exclusion_class: str
    frame0_internal_stop_count: int
    frame1_internal_stop_count: int
    frame2_internal_stop_count: int
    has_ambiguous_bases: bool
    sequence_evaluability: str = "evaluable"

    def row(self) -> dict[str, object]:
        return {
            field: getattr(self, field)
            for field in LEDGER_COLUMNS
        }


@dataclass(slots=True)
class DefiniteInternalStopResult:
    execution_status: str
    total_models: int
    evaluable_models: int = 0
    not_evaluable_models: int = 0
    resolved_frame_models: int = 0
    unresolved_frame_models: int = 0
    excluded_resolved_frame_internal_stop: int = 0
    excluded_unresolved_all_frames_internal_stop: int = 0
    excluded_invalid_cds_geometry: int = 0
    retained_unresolved_with_at_least_one_clean_frame: int = 0
    exclusions: list[DefiniteInternalStopExclusion] = field(default_factory=list)
    excluded_source_counts: Counter[str] = field(default_factory=Counter)

    @property
    def excluded_models(self) -> int:
        return len(self.exclusions)


def empty_definite_internal_stop_result(
    total_models: int,
    *,
    execution_status: str = "prerequisite_unavailable:no_genome_fasta",
) -> DefiniteInternalStopResult:
    return DefiniteInternalStopResult(
        execution_status=execution_status,
        total_models=total_models,
        not_evaluable_models=total_models,
    )


def _phase_provenance(model: CodingModel, active_phase: int | None) -> str:
    boundary = model.five_prime_boundary_canonicalisation
    if boundary is not None and boundary.status == "canonicalised":
        return CANONICALISED_FIVE_PRIME_BOUNDARY
    recovery = model.cds_phase_recovery
    if recovery is not None and recovery.status == RECOVERED_MISSING_SOURCE_PHASE:
        return RECOVERED_MISSING_SOURCE_PHASE
    if active_phase is not None:
        return PRESERVED_NUMERIC_SOURCE
    return UNRESOLVED_SOURCE_MISSING


def _literal_internal_stop_count(
    sequence: str,
    adjacent: str,
    phase: int,
) -> int:
    # ``cds_is_acgt=True`` intentionally requests literal matching even when
    # the sequence contains ambiguity codes: TAR/TAN/etc. are not H3 stops.
    return analyse_cds_frame(
        sequence, adjacent, phase, cds_is_acgt=True
    ).definite_internal_stop_count


def _has_invalid_cds_geometry(model: CodingModel) -> bool:
    """Return whether CDS parts cannot form a non-overlapping splice chain."""
    previous_end0: int | None = None
    for start0, end0 in sorted(model.cds_parts):
        if start0 < 0 or end0 <= start0:
            return True
        if previous_end0 is not None and start0 < previous_end0:
            return True
        previous_end0 = end0
    return False


def recover_phases_and_exclude_definite_internal_stops(
    parsed_inputs: list[ParsedInput],
    fasta: FastaReader,
) -> tuple[Counter[str], list[CodingModel], DefiniteInternalStopResult]:
    """Recover phase and apply H3 using one fetched context per model."""
    total_models = sum(len(parsed.coding_models) for parsed in parsed_inputs)
    result = DefiniteInternalStopResult(
        execution_status="performed", total_models=total_models
    )
    phase_counts: Counter[str] = Counter()
    audit_models: list[CodingModel] = []

    for parsed in parsed_inputs:
        retained: list[CodingModel] = []
        for model in parsed.coding_models:
            sequence: str | None = None
            adjacent = ""
            context_available = True
            if not model.cds_parts:
                context_available = False
            else:
                try:
                    sequence, adjacent = model_sequence_context(model, fasta)
                except (KeyError, ValueError, IndexError):
                    context_available = False

            recovered = recover_model_cds_phase_from_context(
                model,
                sequence,
                adjacent,
                context_available=context_available,
            )
            audit_models.append(recovered)
            phase_counts[phase_recovery_count_key(recovered)] += 1

            # Overlapping CDS features do not define a coherent spliced coding
            # sequence.  Detect them before applying either H3 frame predicate
            # so literal concatenation cannot masquerade as an ordinary
            # all-frame internal-stop result.
            if _has_invalid_cds_geometry(recovered):
                active_phase = source_first_cds_phase(recovered)
                result.not_evaluable_models += 1
                result.excluded_invalid_cds_geometry += 1
                result.exclusions.append(DefiniteInternalStopExclusion(
                    model_id=recovered.output_transcript_id,
                    source=recovered.source_name,
                    input_file=str(recovered.input_file),
                    original_gene_id=recovered.input_gene_id,
                    original_transcript_id=recovered.input_transcript_id,
                    seqid=recovered.chrom,
                    strand=recovered.strand,
                    frame_resolution_status=(
                        "resolved" if active_phase is not None else "unresolved"
                    ),
                    phase_provenance=_phase_provenance(recovered, active_phase),
                    active_first_cds_phase=active_phase,
                    exclusion_class=INVALID_CDS_GEOMETRY,
                    frame0_internal_stop_count=0,
                    frame1_internal_stop_count=0,
                    frame2_internal_stop_count=0,
                    has_ambiguous_bases=(
                        any(base not in "ACGT" for base in sequence)
                        if sequence is not None
                        else False
                    ),
                    sequence_evaluability=INVALID_CDS_GEOMETRY,
                ))
                result.excluded_source_counts[recovered.source_name] += 1
                continue

            if not context_available or not sequence:
                result.not_evaluable_models += 1
                retained.append(recovered)
                continue

            result.evaluable_models += 1
            active_phase = source_first_cds_phase(recovered)
            counts: dict[int, int] = {}
            exclusion_class: str | None = None
            if active_phase is not None:
                result.resolved_frame_models += 1
                counts[active_phase] = _literal_internal_stop_count(
                    sequence, adjacent, active_phase
                )
                if counts[active_phase] > 0:
                    exclusion_class = RESOLVED_FRAME_INTERNAL_STOP
                    result.excluded_resolved_frame_internal_stop += 1
            else:
                result.unresolved_frame_models += 1
                counts = {
                    phase: _literal_internal_stop_count(
                        sequence, adjacent, phase
                    )
                    for phase in (0, 1, 2)
                }
                if all(counts[phase] > 0 for phase in (0, 1, 2)):
                    exclusion_class = UNRESOLVED_ALL_FRAMES_INTERNAL_STOP
                    result.excluded_unresolved_all_frames_internal_stop += 1
                else:
                    result.retained_unresolved_with_at_least_one_clean_frame += 1

            if exclusion_class is None:
                retained.append(recovered)
                continue

            # Populate all three ledger counts only for the small excluded set;
            # resolved retained models require only their active-frame scan.
            for phase in (0, 1, 2):
                counts.setdefault(
                    phase,
                    _literal_internal_stop_count(sequence, adjacent, phase),
                )
            result.exclusions.append(DefiniteInternalStopExclusion(
                model_id=recovered.output_transcript_id,
                source=recovered.source_name,
                input_file=str(recovered.input_file),
                original_gene_id=recovered.input_gene_id,
                original_transcript_id=recovered.input_transcript_id,
                seqid=recovered.chrom,
                strand=recovered.strand,
                frame_resolution_status=(
                    "resolved" if active_phase is not None else "unresolved"
                ),
                phase_provenance=_phase_provenance(recovered, active_phase),
                active_first_cds_phase=active_phase,
                exclusion_class=exclusion_class,
                frame0_internal_stop_count=counts[0],
                frame1_internal_stop_count=counts[1],
                frame2_internal_stop_count=counts[2],
                has_ambiguous_bases=any(base not in "ACGT" for base in sequence),
            ))
            result.excluded_source_counts[recovered.source_name] += 1
        parsed.coding_models = retained

    return phase_counts, audit_models, result


def write_definite_internal_stop_exclusions_tsv(
    result: DefiniteInternalStopResult,
    path: Path,
) -> None:
    rows = [item.row() for item in result.exclusions]
    write_dict_tsv(rows, path, LEDGER_COLUMNS)


def write_definite_internal_stop_summary_tsv(
    result: DefiniteInternalStopResult,
    path: Path,
) -> None:
    source_counts = dict(sorted(result.excluded_source_counts.items()))
    row = {
        "execution_status": result.execution_status,
        "total_models": result.total_models,
        "evaluable_models": result.evaluable_models,
        "not_evaluable_models": result.not_evaluable_models,
        "resolved_frame_models": result.resolved_frame_models,
        "unresolved_frame_models": result.unresolved_frame_models,
        "excluded_models": result.excluded_models,
        "excluded_resolved_frame_internal_stop": result.excluded_resolved_frame_internal_stop,
        "excluded_unresolved_all_frames_internal_stop": result.excluded_unresolved_all_frames_internal_stop,
        "excluded_invalid_cds_geometry": result.excluded_invalid_cds_geometry,
        "retained_unresolved_with_at_least_one_clean_frame": result.retained_unresolved_with_at_least_one_clean_frame,
        "excluded_source_counts": json.dumps(source_counts, separators=(",", ":")),
    }
    write_dict_tsv([row], path, list(row))
