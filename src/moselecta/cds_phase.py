"""Shared CDS-frame analysis and phase propagation policy.

The small frame-evidence helpers are shared by early biological phase recovery
and the later output-only phase-zero decision.  Only the former may alter the
effective biological phase used by sequence integrity.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache

from .models import CodingModel


Interval = tuple[int, int]

PRESERVED_NUMERIC_SOURCE = "preserved_numeric_source"
STANDARDISED_SOURCE_MISSING_TO_ZERO = "standardised_source_missing_to_zero"
UNRESOLVED_SOURCE_MISSING = "unresolved_source_missing"
REFINED_NUMERIC_PHASE = "refined_numeric_phase"
RECOVERED_MISSING_SOURCE_PHASE = "recovered_missing_source_phase"
CANONICALISED_FIVE_PRIME_BOUNDARY = "canonicalised_five_prime_boundary"

STANDARD_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
STANDARD_START_CODONS = frozenset({"ATG"})
PLAIN_DNA_BASES = frozenset("ACGT")

IUPAC_BASES: dict[str, frozenset[str]] = {
    "A": frozenset("A"), "C": frozenset("C"),
    "G": frozenset("G"), "T": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"),
    "S": frozenset("CG"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"),
    "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


@dataclass(frozen=True, slots=True)
class CDSFrameEvidence:
    phase: int
    definite_start: bool
    possible_start: bool
    definite_stop: bool
    possible_stop: bool
    stop_representation: str
    terminal_codon: str
    definite_internal_stop_count: int
    possible_internal_stop_count: int
    possible_included_stop: bool = False
    literal_adjacent_stop: bool = False

    @property
    def definitely_supported(self) -> bool:
        return self.definite_start or self.definite_stop

    @property
    def possibly_supported(self) -> bool:
        return self.possible_start or self.possible_stop

    @property
    def supported_for_output_competition(self) -> bool:
        """Return support relevant to output-only competing-frame vetoes.

        Included stop evidence retains the existing ambiguity-aware behaviour.
        Adjacent stop evidence is positive here only for a literal canonical
        stop; an ambiguous neighbouring codon must not veto otherwise definite
        phase-zero output merely because one expansion could be a stop.
        """
        return (
            self.possible_start
            or self.possible_included_stop
            or self.literal_adjacent_stop
        )


@lru_cache(maxsize=4096)
def _cached_ambiguous_codon_membership(
    codon: str,
    targets: frozenset[str],
) -> tuple[bool, bool]:
    """Return ``(possible, definite)`` for one non-plain codon."""
    possibilities = tuple(IUPAC_BASES.get(base, frozenset()) for base in codon)
    if any(not values for values in possibilities):
        return False, False
    possible = False
    definite = True
    for first in possibilities[0]:
        for second in possibilities[1]:
            for third in possibilities[2]:
                matches = first + second + third in targets
                possible = possible or matches
                definite = definite and matches
    return possible, definite


def codon_membership(
    codon: str,
    targets: frozenset[str],
) -> tuple[bool, bool]:
    """Return ``(possible, definite)`` membership under IUPAC semantics.

    Ordinary DNA codons use a direct set lookup. Only codons containing an
    ambiguity code enter the cached expansion path.
    """
    if len(codon) != 3:
        return False, False
    normalized = codon.upper()
    if all(base in PLAIN_DNA_BASES for base in normalized):
        matches = normalized in targets
        return matches, matches
    return _cached_ambiguous_codon_membership(normalized, targets)


def codon_possibly_in(codon: str, targets: frozenset[str]) -> bool:
    """Return whether at least one IUPAC expansion is in ``targets``."""
    return codon_membership(codon, targets)[0]


def codon_definitely_in(codon: str, targets: frozenset[str]) -> bool:
    """Return whether every possible IUPAC expansion is in ``targets``."""
    return codon_membership(codon, targets)[1]


def analyse_cds_frame(
    cds_sequence: str,
    adjacent_codon: str,
    phase: int,
    *,
    cds_is_acgt: bool = False,
) -> CDSFrameEvidence:
    """Summarize definite and ambiguity-sensitive evidence for one frame."""
    translated_length = max(0, len(cds_sequence) - phase)
    complete_length = (translated_length // 3) * 3
    terminal_aligned = translated_length % 3 == 0
    included = (
        cds_sequence[phase + complete_length - 3 : phase + complete_length]
        if terminal_aligned and complete_length else ""
    )

    def sequence_membership(
        codon: str,
        targets: frozenset[str],
    ) -> tuple[bool, bool]:
        if cds_is_acgt and len(codon) == 3:
            matches = codon in targets
            return matches, matches
        return codon_membership(codon, targets)

    possible_included, definite_included = sequence_membership(
        included, STANDARD_STOP_CODONS
    )
    possible_adjacent_match, definite_adjacent_match = codon_membership(
        adjacent_codon, STANDARD_STOP_CODONS
    )
    definite_adjacent = bool(
        terminal_aligned
        and not definite_included
        and definite_adjacent_match
    )
    possible_adjacent = bool(
        terminal_aligned
        and not possible_included
        and possible_adjacent_match
    )
    if definite_included:
        representation = "included"
        terminal_codon = included
    elif definite_adjacent:
        representation = "adjacent"
        terminal_codon = adjacent_codon
    elif possible_included:
        representation = "ambiguous_included"
        terminal_codon = included
    elif possible_adjacent:
        representation = "ambiguous_adjacent"
        terminal_codon = adjacent_codon
    else:
        representation = "."
        terminal_codon = "."
    internal_length = complete_length - 3 if possible_included else complete_length
    definite_internal_stop_count = 0
    possible_internal_stop_count = 0
    for offset in range(0, internal_length, 3):
        codon = cds_sequence[phase + offset : phase + offset + 3]
        possible, definite = sequence_membership(codon, STANDARD_STOP_CODONS)
        possible_internal_stop_count += int(possible)
        definite_internal_stop_count += int(definite)
    first_codon = cds_sequence[phase : phase + 3]
    possible_start, definite_start = sequence_membership(
        first_codon, STANDARD_START_CODONS
    )
    return CDSFrameEvidence(
        phase=phase,
        definite_start=definite_start,
        possible_start=possible_start,
        definite_stop=definite_included or definite_adjacent,
        possible_stop=possible_included or possible_adjacent,
        stop_representation=representation,
        terminal_codon=terminal_codon,
        definite_internal_stop_count=definite_internal_stop_count,
        possible_internal_stop_count=possible_internal_stop_count,
        possible_included_stop=possible_included,
        literal_adjacent_stop=bool(
            terminal_aligned and adjacent_codon in STANDARD_STOP_CODONS
        ),
    )


def parse_cds_phase(value: object) -> int | None:
    """Return a valid GFF3 CDS phase, leaving missing/invalid phase unresolved."""
    try:
        phase = int(str(value))
    except (TypeError, ValueError):
        return None
    return phase if phase in {0, 1, 2} else None


def source_phase_by_part(model: CodingModel) -> dict[Interval, str]:
    """Return the preferred source phase text for each normalized CDS part."""
    phases: dict[Interval, str] = {}
    for record in model.child_records:
        if record.feature_type.lower() != "cds":
            continue
        part = (record.start0, record.end0)
        if part not in model.cds_parts:
            continue
        current = phases.get(part)
        if current is None or (
            parse_cds_phase(current) is None
            and parse_cds_phase(record.phase) is not None
        ):
            phases[part] = record.phase
    return phases


def biological_cds_parts(model: CodingModel) -> list[Interval]:
    return sorted(model.cds_parts, reverse=model.strand == "-")


def source_first_cds_phase(model: CodingModel) -> int | None:
    ordered = biological_cds_parts(model)
    if not ordered:
        return None
    return parse_cds_phase(source_phase_by_part(model).get(ordered[0], "."))


def calculate_cds_phase_chain(
    parts: Sequence[Interval],
    strand: str,
    first_phase: int,
) -> dict[Interval, str]:
    """Propagate one established first phase in transcript order."""
    if first_phase not in {0, 1, 2}:
        raise ValueError("first CDS phase must be 0, 1, or 2")
    phases: dict[Interval, str] = {}
    cumulative = -first_phase
    for index, part in enumerate(sorted(parts, reverse=strand == "-")):
        phase = first_phase if index == 0 else (3 - (cumulative % 3)) % 3
        phases[part] = str(phase)
        cumulative += part[1] - part[0]
    return phases


def emitted_cds_phase_by_part(model: CodingModel) -> dict[Interval, str]:
    """Return phases for GFF output without inventing an unresolved frame.

    A numeric first phase in the model is authoritative.  A missing first
    phase becomes numeric only when the separately attached output-only
    decision explicitly establishes phase zero.  Otherwise source phase text
    is retained part by part.
    """
    source_phases = source_phase_by_part(model)
    first_phase = source_first_cds_phase(model)
    if first_phase is None:
        decision = model.cds_phase_output
        if (
            decision is not None
            and decision.status == STANDARDISED_SOURCE_MISSING_TO_ZERO
            and decision.first_phase == 0
        ):
            first_phase = 0
    if first_phase is not None:
        return calculate_cds_phase_chain(model.cds_parts, model.strand, first_phase)
    return {part: source_phases.get(part, ".") for part in model.cds_parts}


def cds_phase_output_status(model: CodingModel) -> str:
    """Describe how the first emitted CDS phase was obtained."""
    first_phase = source_first_cds_phase(model)
    if model.cds_terminus_refined is not None and first_phase is not None:
        return REFINED_NUMERIC_PHASE
    boundary = model.five_prime_boundary_canonicalisation
    if boundary is not None and boundary.status == "canonicalised" and first_phase == 0:
        return CANONICALISED_FIVE_PRIME_BOUNDARY
    recovery = model.cds_phase_recovery
    if (
        recovery is not None
        and recovery.status == RECOVERED_MISSING_SOURCE_PHASE
        and first_phase is not None
    ):
        return RECOVERED_MISSING_SOURCE_PHASE
    if first_phase is not None:
        decision = model.cds_phase_output
        if decision is not None and decision.source_first_phase is None:
            return REFINED_NUMERIC_PHASE
        return PRESERVED_NUMERIC_SOURCE
    decision = model.cds_phase_output
    if decision is not None:
        return decision.status
    return UNRESOLVED_SOURCE_MISSING


def cds_phase_output_first_phase(model: CodingModel) -> int | None:
    """Return the first phase that a writer will emit, if established."""
    ordered = biological_cds_parts(model)
    if not ordered:
        return None
    return parse_cds_phase(emitted_cds_phase_by_part(model).get(ordered[0], "."))


def cds_phase_output_reason(model: CodingModel) -> str:
    boundary = model.five_prime_boundary_canonicalisation
    if boundary is not None and boundary.status == "canonicalised":
        return "five_prime_boundary_canonicalised"
    decision = model.cds_phase_output
    if decision is not None:
        return decision.reason
    recovery = model.cds_phase_recovery
    if recovery is not None and recovery.status == RECOVERED_MISSING_SOURCE_PHASE:
        return recovery.reason
    return "numeric_source_phase" if source_first_cds_phase(model) is not None else "not_evaluated"
