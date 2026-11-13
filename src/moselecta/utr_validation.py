"""Shared, post-selection validation of biological UTR splice junctions.

The validator is deliberately independent of CDS scoring and model selection.
Callers decide whether a rejected side is suppressed from an inherited model or
excluded from UTR evidence; the scientific assessment is identical.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Sequence, TextIO

from .genome_sequence import (
    CANONICAL_SPLICE_MOTIFS,
    IndexedFasta,
    reverse_complement,
)
from .tsvio import write_dict_tsv


Interval = tuple[int, int]
Side = Literal["five_prime", "three_prime"]
ValidationStatus = Literal["passed", "rejected", "not_fully_evaluated"]
MotifStatus = Literal["passed", "rejected", "not_evaluated", "not_applicable"]
ObservationMode = Literal["all", "exceptions", "none"]

DEFINITE_REJECTION_REASONS = frozenset({
    "invalid_utr_geometry", "short_utr_intron", "resolved_unrecognized_motif",
})
ROUTINE_UNEVALUATED_REASONS = frozenset({"no_genome_fasta"})


@dataclass(frozen=True, slots=True)
class UTRJunctionObservation:
    intron_start0: int
    intron_end0: int
    intron_length: int
    motif: str
    motif_status: MotifStatus
    length_status: Literal["passed", "rejected"]
    reason_codes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class UTRSideValidation:
    transcript_id: str
    evidence_origin: str
    evidence_id: str
    chrom: str
    strand: str
    side: Side
    minimum_intron_bp: int
    status: ValidationStatus
    reason_codes: tuple[str, ...]
    reason_counts: tuple[tuple[str, int], ...]
    junction_count: int
    motif_passed_count: int
    motif_not_evaluated_count: int
    length_violation_count: int
    observations: tuple[UTRJunctionObservation, ...]

    @property
    def retained(self) -> bool:
        return self.status != "rejected"


def splice_motif(
    fasta: IndexedFasta,
    chrom: str,
    intron: Interval,
    strand: str,
) -> tuple[str, MotifStatus, str | None]:
    """Return one strand-oriented donor-acceptor motif and its assessability."""
    start0, end0 = intron
    contig_length = fasta.contig_length(chrom)
    if contig_length is None:
        return ".", "not_evaluated", "genome_contig_missing"
    if start0 < 0 or end0 > contig_length:
        return ".", "not_evaluated", "splice_context_out_of_range"
    if end0 - start0 < 2:
        return ".", "not_evaluated", "splice_context_too_short"
    try:
        if strand == "+":
            donor = fasta.fetch(chrom, start0, start0 + 2)
            acceptor = fasta.fetch(chrom, end0 - 2, end0)
        elif strand == "-":
            donor = reverse_complement(fasta.fetch(chrom, end0 - 2, end0))
            acceptor = reverse_complement(fasta.fetch(chrom, start0, start0 + 2))
        else:
            return ".", "not_evaluated", "invalid_strand"
    except (KeyError, ValueError):
        return ".", "not_evaluated", "splice_context_unavailable"
    motif = f"{donor}-{acceptor}"
    if any(base not in "ACGT" for base in donor + acceptor):
        return motif, "not_evaluated", "ambiguous_splice_sequence"
    if motif in CANONICAL_SPLICE_MOTIFS:
        return motif, "passed", None
    return motif, "rejected", "resolved_unrecognized_motif"


def validate_utr_side(
    *,
    transcript_id: str,
    evidence_origin: str,
    evidence_id: str,
    chrom: str,
    strand: str,
    side: Side,
    introns: Sequence[Interval],
    minimum_intron_bp: int,
    fasta: IndexedFasta | None = None,
    geometry_valid: bool = True,
    observation_mode: ObservationMode = "all",
) -> UTRSideValidation:
    """Validate every UTR-associated intron on one complete biological side."""
    side_reasons: set[str] = set()
    reason_counts: Counter[str] = Counter()
    observations: list[UTRJunctionObservation] = []
    if strand not in {"+", "-"} or not geometry_valid:
        side_reasons.add("invalid_utr_geometry")
        reason_counts["invalid_utr_geometry"] += 1

    ordered = tuple(sorted(introns))
    if any(start0 < 0 or end0 <= start0 for start0, end0 in ordered):
        side_reasons.add("invalid_utr_geometry")
        reason_counts["invalid_utr_geometry"] = 1
    if any(left[1] > right[0] for left, right in zip(ordered, ordered[1:])):
        side_reasons.add("invalid_utr_geometry")
        reason_counts["invalid_utr_geometry"] = 1

    motif_passed_count = 0
    motif_not_evaluated_count = 0
    length_violation_count = 0
    for start0, end0 in ordered:
        length = end0 - start0
        reasons: set[str] = set()
        length_status: Literal["passed", "rejected"] = "passed"
        if length < minimum_intron_bp:
            length_status = "rejected"
            reasons.add("short_utr_intron")
            length_violation_count += 1

        if fasta is None:
            motif = "."
            motif_status: MotifStatus = "not_evaluated"
            motif_reason = "no_genome_fasta"
        else:
            motif, motif_status, motif_reason = splice_motif(
                fasta, chrom, (start0, end0), strand
            )
        if motif_reason is not None:
            reasons.add(motif_reason)
        if motif_status == "passed":
            motif_passed_count += 1
        elif motif_status == "not_evaluated":
            motif_not_evaluated_count += 1
        side_reasons.update(reasons)
        reason_counts.update(reasons)
        retain_observation = (
            observation_mode == "all"
            or observation_mode == "exceptions" and bool(
                reasons - ROUTINE_UNEVALUATED_REASONS
            )
        )
        if retain_observation:
            observations.append(
                UTRJunctionObservation(
                    intron_start0=start0,
                    intron_end0=end0,
                    intron_length=length,
                    motif=motif,
                    motif_status=motif_status,
                    length_status=length_status,
                    reason_codes=tuple(sorted(reasons)),
                )
            )

    definite_failure = bool(DEFINITE_REJECTION_REASONS & side_reasons)
    unevaluated = motif_not_evaluated_count > 0
    status: ValidationStatus = (
        "rejected" if definite_failure
        else "not_fully_evaluated" if unevaluated
        else "passed"
    )
    return UTRSideValidation(
        transcript_id=transcript_id,
        evidence_origin=evidence_origin,
        evidence_id=evidence_id,
        chrom=chrom,
        strand=strand,
        side=side,
        minimum_intron_bp=minimum_intron_bp,
        status=status,
        reason_codes=tuple(sorted(side_reasons)),
        reason_counts=tuple(sorted(reason_counts.items())),
        junction_count=len(ordered),
        motif_passed_count=motif_passed_count,
        motif_not_evaluated_count=motif_not_evaluated_count,
        length_violation_count=length_violation_count,
        observations=tuple(observations),
    )


def utr_side_introns(
    exon_parts: Sequence[Interval],
    cds_parts: Sequence[Interval],
    strand: str,
    side: Side,
) -> tuple[Interval, ...]:
    """Return UTR-only and CDS/UTR-transition introns for one biological side."""
    exons = tuple(sorted(exon_parts))
    cds = tuple(sorted(cds_parts))
    if not exons or not cds:
        return ()
    low_side = (side == "five_prime" and strand == "+") or (
        side == "three_prime" and strand == "-"
    )
    boundary = cds[0][0] if low_side else cds[-1][1]
    introns: list[Interval] = []
    for left, right in zip(exons, exons[1:]):
        intron = (left[1], right[0])
        if low_side and left[1] <= boundary:
            introns.append(intron)
        elif not low_side and right[0] >= boundary:
            introns.append(intron)
    return tuple(introns)


def side_has_utr(
    exon_parts: Sequence[Interval],
    cds_parts: Sequence[Interval],
    strand: str,
    side: Side,
) -> bool:
    cds_start0 = min(part[0] for part in cds_parts)
    cds_end0 = max(part[1] for part in cds_parts)
    low_side = (side == "five_prime" and strand == "+") or (
        side == "three_prime" and strand == "-"
    )
    if low_side:
        return any(start0 < cds_start0 for start0, _end0 in exon_parts)
    return any(end0 > cds_end0 for _start0, end0 in exon_parts)


UTR_VALIDATION_AUDIT_COLUMNS = (
    "transcript_id", "evidence_origin", "evidence_id", "chrom", "strand",
    "utr_side", "intron_start", "intron_end", "intron_length",
    "minimum_intron_bp", "splice_motif", "motif_status", "length_status",
    "side_status", "exception_scope", "side_reason_codes",
    "intron_reason_codes", "reason_codes", "final_action",
)

UTR_VALIDATION_SUMMARY_COLUMNS = (
    "evidence_origin", "sides_assessed", "sides_accepted", "sides_rejected",
    "sides_not_fully_evaluated", "sides_without_splice_junction",
    "junctions_assessed", "motif_checks_passed", "motif_checks_not_evaluated",
    "length_violations", "individual_rejection_violations",
    "five_prime_sides_assessed", "five_prime_sides_rejected",
    "three_prime_sides_assessed", "three_prime_sides_rejected",
    "rejection_reason_counts", "unevaluable_reason_counts",
)


def validation_audit_rows(
    validation: UTRSideValidation,
    *,
    rejected_action: str,
) -> list[dict[str, object]]:
    """Return only consequential exception rows for one side decision."""
    if validation.status == "passed":
        return []
    if (
        validation.status == "not_fully_evaluated"
        and set(validation.reason_codes) <= ROUTINE_UNEVALUATED_REASONS
    ):
        return []
    action = (
        rejected_action if validation.status == "rejected"
        else "retain_unevaluated" if validation.status == "not_fully_evaluated"
        else "retain_side"
    )
    observations: list[UTRJunctionObservation | None] = []
    for observation in validation.observations:
        reasons = set(observation.reason_codes)
        if validation.status == "rejected":
            if reasons & DEFINITE_REJECTION_REASONS:
                observations.append(observation)
        elif reasons - ROUTINE_UNEVALUATED_REASONS:
            observations.append(observation)
    geometry_only = (
        "invalid_utr_geometry" in validation.reason_codes
        and not any(
            observation is not None
            and set(observation.reason_codes) & DEFINITE_REJECTION_REASONS
            for observation in observations
        )
    )
    if geometry_only or not observations:
        observations.append(None)
    rows: list[dict[str, object]] = []
    for observation in observations:
        side_reasons = ",".join(validation.reason_codes) or "."
        intron_reasons = (
            ",".join(observation.reason_codes) if observation else "."
        ) or "."
        rows.append({
            "transcript_id": validation.transcript_id,
            "evidence_origin": validation.evidence_origin,
            "evidence_id": validation.evidence_id,
            "chrom": validation.chrom,
            "strand": validation.strand,
            "utr_side": validation.side,
            "intron_start": observation.intron_start0 + 1 if observation else ".",
            "intron_end": observation.intron_end0 if observation else ".",
            "intron_length": observation.intron_length if observation else 0,
            "minimum_intron_bp": validation.minimum_intron_bp,
            "splice_motif": observation.motif if observation else ".",
            "motif_status": observation.motif_status if observation else "not_applicable",
            "length_status": observation.length_status if observation else "not_applicable",
            "side_status": validation.status,
            "exception_scope": "junction" if observation else "side",
            "side_reason_codes": side_reasons,
            "intron_reason_codes": intron_reasons,
            "reason_codes": side_reasons,
            "final_action": action,
        })
    return rows


def write_utr_validation_audit(
    rows: Sequence[dict[str, object]], path: Path
) -> None:
    write_dict_tsv(list(rows), path, list(UTR_VALIDATION_AUDIT_COLUMNS))


class UTRValidationRecorder:
    """Bounded-memory counters plus an optional streaming exception audit."""

    def __init__(self, detail_path: Path | None = None) -> None:
        self._counts: Counter[tuple[str, str]] = Counter()
        self._reason_counts: Counter[tuple[str, str]] = Counter()
        self._handle: TextIO | None = None
        self._writer: csv.DictWriter | None = None
        self.detail_rows_written = 0
        if detail_path is not None:
            detail_path.parent.mkdir(parents=True, exist_ok=True)
            self._handle = detail_path.open("w", newline="")
            self._writer = csv.DictWriter(
                self._handle,
                fieldnames=list(UTR_VALIDATION_AUDIT_COLUMNS),
                delimiter="\t",
                lineterminator="\n",
            )
            self._writer.writeheader()

    @property
    def detail_enabled(self) -> bool:
        return self._writer is not None

    def record(
        self,
        validation: UTRSideValidation,
        *,
        rejected_action: str,
    ) -> None:
        origins = ("all", validation.evidence_origin)
        reason_counts = dict(validation.reason_counts)
        for origin in origins:
            self._counts[(origin, "sides_assessed")] += 1
            self._counts[(origin, f"sides_{validation.status}")] += 1
            self._counts[(origin, f"{validation.side}_sides_assessed")] += 1
            if validation.status == "rejected":
                self._counts[(origin, f"{validation.side}_sides_rejected")] += 1
            if validation.junction_count == 0:
                self._counts[(origin, "sides_without_splice_junction")] += 1
            self._counts[(origin, "junctions_assessed")] += validation.junction_count
            self._counts[(origin, "motif_checks_passed")] += (
                validation.motif_passed_count
            )
            self._counts[(origin, "motif_checks_not_evaluated")] += (
                validation.motif_not_evaluated_count
            )
            self._counts[(origin, "length_violations")] += (
                validation.length_violation_count
            )
            for reason, count in reason_counts.items():
                self._reason_counts[(origin, reason)] += count
                if reason in DEFINITE_REJECTION_REASONS:
                    self._counts[(origin, "individual_rejection_violations")] += count
        if self._writer is not None:
            for row in validation_audit_rows(
                validation, rejected_action=rejected_action
            ):
                self._writer.writerow(row)
                self.detail_rows_written += 1

    def summary_rows(self) -> list[dict[str, object]]:
        origins = sorted(
            {origin for origin, _metric in self._counts if origin != "all"}
        )
        rows: list[dict[str, object]] = []
        for origin in ("all", *origins):
            rejection_reasons = {
                reason: count
                for (item_origin, reason), count in self._reason_counts.items()
                if item_origin == origin and reason in DEFINITE_REJECTION_REASONS
            }
            unevaluable_reasons = {
                reason: count
                for (item_origin, reason), count in self._reason_counts.items()
                if item_origin == origin and reason not in DEFINITE_REJECTION_REASONS
            }
            rows.append({
                "evidence_origin": origin,
                "sides_assessed": self._counts[(origin, "sides_assessed")],
                "sides_accepted": self._counts[(origin, "sides_passed")],
                "sides_rejected": self._counts[(origin, "sides_rejected")],
                "sides_not_fully_evaluated": self._counts[
                    (origin, "sides_not_fully_evaluated")
                ],
                "sides_without_splice_junction": self._counts[
                    (origin, "sides_without_splice_junction")
                ],
                "junctions_assessed": self._counts[(origin, "junctions_assessed")],
                "motif_checks_passed": self._counts[(origin, "motif_checks_passed")],
                "motif_checks_not_evaluated": self._counts[
                    (origin, "motif_checks_not_evaluated")
                ],
                "length_violations": self._counts[(origin, "length_violations")],
                "individual_rejection_violations": self._counts[
                    (origin, "individual_rejection_violations")
                ],
                "five_prime_sides_assessed": self._counts[
                    (origin, "five_prime_sides_assessed")
                ],
                "five_prime_sides_rejected": self._counts[
                    (origin, "five_prime_sides_rejected")
                ],
                "three_prime_sides_assessed": self._counts[
                    (origin, "three_prime_sides_assessed")
                ],
                "three_prime_sides_rejected": self._counts[
                    (origin, "three_prime_sides_rejected")
                ],
                "rejection_reason_counts": _compact_counts(rejection_reasons),
                "unevaluable_reason_counts": _compact_counts(unevaluable_reasons),
            })
        return rows

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None
            self._writer = None


def _compact_counts(counts: dict[str, int]) -> str:
    return ",".join(f"{key}:{counts[key]}" for key in sorted(counts)) or "."


def write_utr_validation_summary(
    recorder: UTRValidationRecorder,
    path: Path,
) -> None:
    write_dict_tsv(
        recorder.summary_rows(), path, list(UTR_VALIDATION_SUMMARY_COLUMNS)
    )
