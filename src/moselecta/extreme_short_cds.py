"""Early hard exclusion of candidates below the coding-length floor."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from .cds_phase import analyse_cds_frame, source_first_cds_phase
from .models import CodingModel, ParsedInput
from .phase_recovery import FastaReader, model_sequence_context
from .tsvio import write_dict_tsv


EXCLUSION_CLASS = "extreme_short_cds"
INCLUDED_STOP_SUBTRACTION_BP = 3

LEDGER_COLUMNS = [
    "model_id", "source", "input_file", "original_gene_id",
    "original_transcript_id", "seqid", "strand", "total_cds_length_bp",
    "terminal_stop_status", "included_terminal_stop_subtracted",
    "terminal_stop_subtraction_bp", "coding_length_bp",
    "minimum_coding_length_bp", "cds_exon_count", "exclusion_class",
]


@dataclass(frozen=True, slots=True)
class ExtremeShortCDSExclusion:
    model_id: str
    source: str
    input_file: str
    original_gene_id: str
    original_transcript_id: str
    seqid: str
    strand: str
    total_cds_length_bp: int
    terminal_stop_status: str
    included_terminal_stop_subtracted: bool
    terminal_stop_subtraction_bp: int
    coding_length_bp: int
    minimum_coding_length_bp: int
    cds_exon_count: int
    exclusion_class: str = EXCLUSION_CLASS

    def row(self) -> dict[str, object]:
        return {column: getattr(self, column) for column in LEDGER_COLUMNS}


@dataclass(slots=True)
class ExtremeShortCDSResult:
    execution_status: str
    minimum_coding_length_bp: int
    candidates_evaluated: int
    retained_models: int = 0
    exclusions: list[ExtremeShortCDSExclusion] = field(default_factory=list)
    excluded_source_counts: Counter[str] = field(default_factory=Counter)
    excluded_cds_exon_count_counts: Counter[int] = field(default_factory=Counter)

    @property
    def excluded_models(self) -> int:
        return len(self.exclusions)


def _terminal_stop_status(
    model: CodingModel,
    fasta: FastaReader | None,
) -> str:
    """Return the existing frame analyser's terminal-stop representation."""
    if fasta is None:
        return "not_evaluable:no_genome_fasta"
    try:
        sequence, adjacent = model_sequence_context(model, fasta)
    except (KeyError, ValueError, IndexError):
        return "not_evaluable:sequence_context_unavailable"
    phase = source_first_cds_phase(model)
    frame = analyse_cds_frame(
        sequence,
        adjacent,
        phase if phase is not None else 0,
        cds_is_acgt=True,
    )
    return {
        "included": "included_in_cds",
        "adjacent": "adjacent_to_cds",
        "ambiguous_included": "ambiguous_included",
        "ambiguous_adjacent": "ambiguous_adjacent",
        ".": "missing",
    }[frame.stop_representation]


def exclude_extreme_short_cds(
    parsed_inputs: list[ParsedInput],
    *,
    minimum_coding_length_bp: int,
    fasta: FastaReader | None = None,
) -> ExtremeShortCDSResult:
    """Remove candidates below the post-correction coding-length floor.

    Only candidates whose raw CDS length is within three bases of the floor
    need sequence inspection. Longer candidates cannot cross the threshold
    after the sole permitted adjustment (subtracting one included stop codon).
    """
    if minimum_coding_length_bp <= 0:
        raise ValueError("minimum_coding_length_bp must be > 0")
    total = sum(len(parsed.coding_models) for parsed in parsed_inputs)
    result = ExtremeShortCDSResult(
        execution_status="performed",
        minimum_coding_length_bp=minimum_coding_length_bp,
        candidates_evaluated=total,
    )

    for parsed in parsed_inputs:
        retained: list[CodingModel] = []
        for model in parsed.coding_models:
            total_length = model.cds_length
            if total_length >= minimum_coding_length_bp + INCLUDED_STOP_SUBTRACTION_BP:
                result.retained_models += 1
                retained.append(model)
                continue

            stop_status = _terminal_stop_status(model, fasta)
            subtract_stop = stop_status == "included_in_cds"
            subtraction = INCLUDED_STOP_SUBTRACTION_BP if subtract_stop else 0
            coding_length = total_length - subtraction
            if coding_length >= minimum_coding_length_bp:
                result.retained_models += 1
                retained.append(model)
                continue

            exclusion = ExtremeShortCDSExclusion(
                model_id=model.output_transcript_id,
                source=model.source_name,
                input_file=str(model.input_file),
                original_gene_id=model.input_gene_id,
                original_transcript_id=model.input_transcript_id,
                seqid=model.chrom,
                strand=model.strand,
                total_cds_length_bp=total_length,
                terminal_stop_status=stop_status,
                included_terminal_stop_subtracted=subtract_stop,
                terminal_stop_subtraction_bp=subtraction,
                coding_length_bp=coding_length,
                minimum_coding_length_bp=minimum_coding_length_bp,
                cds_exon_count=model.cds_exon_count,
            )
            result.exclusions.append(exclusion)
            result.excluded_source_counts[model.source_name] += 1
            result.excluded_cds_exon_count_counts[model.cds_exon_count] += 1
        parsed.coding_models = retained

    if result.retained_models + result.excluded_models != total:
        raise AssertionError("extreme-short CDS accounting did not reconcile")
    return result


def write_extreme_short_cds_exclusions_tsv(
    result: ExtremeShortCDSResult,
    path: Path,
) -> None:
    write_dict_tsv(
        [exclusion.row() for exclusion in result.exclusions],
        path,
        LEDGER_COLUMNS,
    )


def write_extreme_short_cds_summary_tsv(
    result: ExtremeShortCDSResult,
    path: Path,
) -> None:
    row = {
        "execution_status": result.execution_status,
        "minimum_coding_length_bp": result.minimum_coding_length_bp,
        "candidates_evaluated": result.candidates_evaluated,
        "excluded_models": result.excluded_models,
        "retained_models": result.retained_models,
        "included_terminal_stop_subtraction_bp": INCLUDED_STOP_SUBTRACTION_BP,
        "excluded_source_counts": json.dumps(
            dict(sorted(result.excluded_source_counts.items())),
            separators=(",", ":"),
        ),
        "excluded_cds_exon_count_counts": json.dumps(
            {
                str(key): value
                for key, value in sorted(
                    result.excluded_cds_exon_count_counts.items()
                )
            },
            separators=(",", ":"),
        ),
    }
    write_dict_tsv([row], path, list(row))
