"""Historical synthetic characterization of CDS terminus completion.

QUARANTINE NOTE: all completion helpers in this module reproduce the old,
pre-implementation bounded-search experiment.  They deliberately materialize
per-base transcript coordinates and accept a historical upstream distance
limit.  Neither representation nor limit describes the current production
algorithm.  Current behavior is covered by ``test_cds_terminus_refinement.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

from moselecta.final_annotation import FinalTranscript
from moselecta.models import CodingModel
from moselecta.utr_augmentation import (
    UTRCandidateSide,
    UTREvidence,
    UTRLearnedLimits,
    UTRSettings,
    _extension_reasons,
    _ranked_structures,
)


Interval = tuple[int, int]
STOPS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class PrototypeCompletion:
    success: bool
    direction: str
    reason: str
    codon: str = "."
    codon_positions: tuple[int, ...] = ()
    cds_parts: tuple[Interval, ...] = ()
    utr_parts: tuple[Interval, ...] = ()
    cds_bases_delta: int = 0
    transcript_bases_searched: int = 0
    splice_junctions_crossed: int = 0


def _historical_materialized_transcript_positions(
    exons: tuple[Interval, ...], strand: str = "+"
) -> list[int]:
    """Reproduce the obsolete prototype's per-base coordinate materialization."""
    ordered = sorted(exons, reverse=strand == "-")
    positions: list[int] = []
    for start0, end0 in ordered:
        positions.extend(
            range(end0 - 1, start0 - 1, -1)
            if strand == "-" else range(start0, end0)
        )
    return positions


def _genome_for_transcript(
    exons: tuple[Interval, ...], transcript_sequence: str, strand: str = "+"
) -> str:
    positions = _historical_materialized_transcript_positions(exons, strand)
    assert len(positions) == len(transcript_sequence)
    genome = ["C"] * (max(positions) + 1)
    for position, base in zip(positions, transcript_sequence):
        genome[position] = base if strand == "+" else base.translate(COMPLEMENT)
    return "".join(genome)


def _oriented_sequence(genome: str, positions: list[int], strand: str) -> str:
    sequence = "".join(genome[position] for position in positions)
    return sequence if strand == "+" else sequence.translate(COMPLEMENT)


def _parts_from_positions(positions: list[int]) -> tuple[Interval, ...]:
    if not positions:
        return ()
    ordered = sorted(set(positions))
    parts: list[Interval] = []
    start = previous = ordered[0]
    for position in ordered[1:]:
        if position != previous + 1:
            parts.append((start, previous + 1))
            start = position
        previous = position
    parts.append((start, previous + 1))
    return tuple(parts)


def _junction_count(positions: list[int]) -> int:
    return sum(abs(right - left) != 1 for left, right in zip(positions, positions[1:]))


def _coding_bounds(
    exons: tuple[Interval, ...], cds_parts: tuple[Interval, ...], strand: str
) -> tuple[list[int], int, int]:
    transcript_positions = _historical_materialized_transcript_positions(
        exons, strand
    )
    cds_positions = _historical_materialized_transcript_positions(cds_parts, strand)
    start_index = transcript_positions.index(cds_positions[0])
    end_index = transcript_positions.index(cds_positions[-1])
    assert transcript_positions[start_index : end_index + 1] == cds_positions
    return transcript_positions, start_index, end_index


def _contains_internal_stop(sequence: str) -> bool:
    complete = len(sequence) // 3 * 3
    codons = [sequence[index : index + 3] for index in range(0, complete, 3)]
    terminal_is_stop = len(sequence) % 3 == 0 and bool(codons) and codons[-1] in STOPS
    internal = codons[:-1] if terminal_is_stop else codons
    return any(codon in STOPS for codon in internal)


def _historical_bounded_upstream_start(
    genome: str,
    exons: tuple[Interval, ...],
    cds_parts: tuple[Interval, ...],
    *,
    historical_limit: int,
    strand: str = "+",
) -> PrototypeCompletion:
    positions, current_start, current_end = _coding_bounds(exons, cds_parts, strand)
    sequence = _oriented_sequence(genome, positions, strand)
    maximum = min(historical_limit, current_start)
    searched = 0
    for distance in range(3, maximum + 1, 3):
        searched = distance
        codon_start = current_start - distance
        codon = sequence[codon_start : codon_start + 3]
        if codon in STOPS:
            return PrototypeCompletion(
                False, "upstream_start", "upstream_in_frame_stop_barrier",
                transcript_bases_searched=searched,
            )
        if codon == "ATG":
            new_positions = positions[codon_start : current_end + 1]
            added_positions = positions[codon_start:current_start]
            return PrototypeCompletion(
                True,
                "upstream_start",
                "completed",
                codon,
                tuple(positions[codon_start : codon_start + 3]),
                _parts_from_positions(new_positions),
                _parts_from_positions(positions[:codon_start]),
                len(added_positions),
                searched,
                _junction_count(added_positions + positions[current_start:current_start + 1]),
            )
    return PrototypeCompletion(
        False,
        "upstream_start",
        "start_not_found_within_limit",
        transcript_bases_searched=searched,
    )


def _downstream_start(
    genome: str,
    exons: tuple[Interval, ...],
    cds_parts: tuple[Interval, ...],
    *,
    limit: int,
    strand: str = "+",
) -> PrototypeCompletion:
    positions, current_start, current_end = _coding_bounds(exons, cds_parts, strand)
    sequence = _oriented_sequence(genome, positions, strand)
    maximum = min(limit, current_end - current_start - 2)
    searched = 0
    for offset in range(3, maximum + 1, 3):
        searched = offset
        codon_start = current_start + offset
        codon = sequence[codon_start : codon_start + 3]
        if codon != "ATG":
            continue
        corrected_sequence = sequence[codon_start : current_end + 1]
        if _contains_internal_stop(corrected_sequence):
            return PrototypeCompletion(
                False,
                "downstream_start_fallback",
                "proposed_cds_has_internal_stop",
                codon,
                tuple(positions[codon_start : codon_start + 3]),
                transcript_bases_searched=searched,
            )
        removed_positions = positions[current_start:codon_start]
        return PrototypeCompletion(
            True,
            "downstream_start_fallback",
            "completed",
            codon,
            tuple(positions[codon_start : codon_start + 3]),
            _parts_from_positions(positions[codon_start : current_end + 1]),
            _parts_from_positions(positions[:codon_start]),
            -len(removed_positions),
            searched,
            _junction_count(removed_positions + positions[codon_start:codon_start + 1]),
        )
    return PrototypeCompletion(
        False,
        "downstream_start_fallback",
        "start_not_found_within_limit",
        transcript_bases_searched=searched,
    )


def _historical_bounded_complete_start(
    genome: str,
    exons: tuple[Interval, ...],
    cds_parts: tuple[Interval, ...],
    *,
    historical_upstream_limit: int,
    downstream_limit: int,
    strand: str = "+",
) -> PrototypeCompletion:
    upstream = _historical_bounded_upstream_start(
        genome,
        exons,
        cds_parts,
        historical_limit=historical_upstream_limit,
        strand=strand,
    )
    if upstream.success:
        return upstream
    return _downstream_start(
        genome, exons, cds_parts, limit=downstream_limit, strand=strand
    )


def _downstream_stop(
    genome: str,
    exons: tuple[Interval, ...],
    cds_parts: tuple[Interval, ...],
    *,
    limit: int,
    strand: str = "+",
) -> PrototypeCompletion:
    positions, current_start, current_end = _coding_bounds(exons, cds_parts, strand)
    sequence = _oriented_sequence(genome, positions, strand)
    coding_and_extension = sequence[current_start:]
    old_length = current_end - current_start + 1
    first_codon_start = old_length if old_length % 3 == 0 else old_length - old_length % 3
    searched = 0
    for codon_start in range(first_codon_start, len(coding_and_extension) - 2, 3):
        codon_end = codon_start + 3
        added = codon_end - old_length
        if added <= 0:
            continue
        if added > limit:
            break
        searched = added
        codon = coding_and_extension[codon_start:codon_end]
        if codon not in STOPS:
            continue
        corrected_sequence = coding_and_extension[:codon_end]
        if _contains_internal_stop(corrected_sequence):
            return PrototypeCompletion(
                False,
                "downstream_stop",
                "proposed_cds_has_internal_stop",
                codon,
                tuple(positions[current_start + codon_start : current_start + codon_end]),
                transcript_bases_searched=searched,
            )
        new_end = current_start + codon_end
        added_positions = positions[current_end + 1:new_end]
        return PrototypeCompletion(
            True,
            "downstream_stop",
            "completed",
            codon,
            tuple(positions[current_start + codon_start : new_end]),
            _parts_from_positions(positions[current_start:new_end]),
            _parts_from_positions(positions[new_end:]),
            len(added_positions),
            searched,
            _junction_count(positions[current_end:new_end]),
        )
    return PrototypeCompletion(
        False,
        "downstream_stop",
        "stop_not_found_within_limit",
        transcript_bases_searched=searched,
    )


def _sequence(length: int, *overrides: tuple[int, str]) -> str:
    bases = ["C"] * length
    for start, value in overrides:
        bases[start : start + len(value)] = value
    return "".join(bases)


def test_upstream_start_within_50_bp_and_residual_sequence_becomes_utr() -> None:
    exons = ((0, 90),)
    sequence = _sequence(90, (24, "ATG"), (87, "TAA"))
    result = _historical_bounded_upstream_start(
        _genome_for_transcript(exons, sequence),
        exons,
        ((30, 90),),
        historical_limit=50,
    )
    assert result.success and result.codon_positions == (24, 25, 26)
    assert result.cds_parts == ((24, 90),)
    assert result.utr_parts == ((0, 24),)
    assert result.cds_bases_delta == 6


@pytest.mark.parametrize(
    ("atg_offset", "limit", "reason"),
    [
        (30, 20, "start_not_found_within_limit"),
        (None, 50, "start_not_found_within_limit"),
    ],
)
def test_upstream_start_beyond_limit_or_absent(
    atg_offset: int | None, limit: int, reason: str
) -> None:
    overrides = [(87, "TAA")]
    if atg_offset is not None:
        overrides.append((30 - atg_offset, "ATG"))
    exons = ((0, 90),)
    result = _historical_bounded_upstream_start(
        _genome_for_transcript(exons, _sequence(90, *overrides)),
        exons,
        ((30, 90),),
        historical_limit=limit,
    )
    assert not result.success and result.reason == reason


def test_nearest_upstream_atg_is_preferred_when_several_are_in_frame() -> None:
    exons = ((0, 90),)
    sequence = _sequence(90, (18, "ATG"), (24, "ATG"), (87, "TAA"))
    result = _historical_bounded_upstream_start(
        _genome_for_transcript(exons, sequence),
        exons,
        ((30, 90),),
        historical_limit=50,
    )
    assert result.codon_positions == (24, 25, 26)


def test_upstream_start_codon_can_cross_a_splice_junction() -> None:
    exons = ((0, 5), (20, 61))
    sequence = _sequence(46, (4, "ATG"), (43, "TAA"))
    result = _historical_bounded_upstream_start(
        _genome_for_transcript(exons, sequence),
        exons,
        ((22, 61),),
        historical_limit=20,
    )
    assert result.success
    assert result.codon_positions == (4, 20, 21)
    assert result.cds_parts == ((4, 5), (20, 61))
    assert result.splice_junctions_crossed == 1


def test_upstream_in_frame_stop_is_a_barrier_to_more_distant_atg() -> None:
    exons = ((0, 90),)
    sequence = _sequence(90, (24, "ATG"), (27, "TGA"), (87, "TAA"))
    result = _historical_bounded_upstream_start(
        _genome_for_transcript(exons, sequence),
        exons,
        ((30, 90),),
        historical_limit=50,
    )
    assert not result.success
    assert result.reason == "upstream_in_frame_stop_barrier"


@pytest.mark.parametrize("fallback_limit", [15, 30, 50])
def test_downstream_start_fallback_within_exercised_windows(fallback_limit: int) -> None:
    exons = ((0, 90),)
    sequence = _sequence(90, (33, "ATG"), (87, "TAA"))
    result = _historical_bounded_complete_start(
        _genome_for_transcript(exons, sequence),
        exons,
        ((30, 90),),
        historical_upstream_limit=30,
        downstream_limit=fallback_limit,
    )
    assert result.success and result.direction == "downstream_start_fallback"
    assert result.cds_parts == ((33, 90),)
    assert result.utr_parts == ((0, 33),)
    assert result.cds_bases_delta == -3


def test_downstream_start_beyond_fallback_window_is_not_used() -> None:
    exons = ((0, 120),)
    sequence = _sequence(120, (63, "ATG"), (117, "TAA"))
    result = _historical_bounded_complete_start(
        _genome_for_transcript(exons, sequence),
        exons,
        ((30, 120),),
        historical_upstream_limit=20,
        downstream_limit=30,
    )
    assert not result.success


def test_upstream_start_wins_when_upstream_and_downstream_are_valid() -> None:
    exons = ((0, 90),)
    sequence = _sequence(90, (24, "ATG"), (33, "ATG"), (87, "TAA"))
    result = _historical_bounded_complete_start(
        _genome_for_transcript(exons, sequence),
        exons,
        ((30, 90),),
        historical_upstream_limit=50,
        downstream_limit=30,
    )
    assert result.success and result.direction == "upstream_start"
    assert result.codon_positions == (24, 25, 26)


def test_downstream_start_can_cross_existing_splice_structure() -> None:
    exons = ((0, 5), (20, 60))
    sequence = _sequence(45, (3, "ATG"), (42, "TAA"))
    result = _downstream_start(
        _genome_for_transcript(exons, sequence),
        exons,
        exons,
        limit=15,
    )
    assert result.success
    assert result.codon_positions == (3, 4, 20)
    assert result.cds_parts == ((3, 5), (20, 60))


def test_downstream_start_proposal_with_internal_stop_fails_locally() -> None:
    exons = ((0, 90),)
    sequence = _sequence(90, (33, "ATG"), (39, "TGA"), (87, "TAA"))
    result = _downstream_start(
        _genome_for_transcript(exons, sequence), exons, ((30, 90),), limit=30
    )
    assert not result.success
    assert result.reason == "proposed_cds_has_internal_stop"


@pytest.mark.parametrize("limit", [50, 100, 200])
def test_downstream_stop_within_exercised_windows(limit: int) -> None:
    exons = ((0, 120),)
    sequence = _sequence(120, (30, "TAA"))
    result = _downstream_stop(
        _genome_for_transcript(exons, sequence), exons, ((0, 30),), limit=limit
    )
    assert result.success and result.codon == "TAA"
    assert result.cds_parts == ((0, 33),)


@pytest.mark.parametrize(
    ("stop_offset", "limit"),
    [(60, 50), (None, 100)],
)
def test_downstream_stop_beyond_limit_or_absent(
    stop_offset: int | None, limit: int
) -> None:
    exons = ((0, 150),)
    overrides = [] if stop_offset is None else [(30 + stop_offset, "TAG")]
    result = _downstream_stop(
        _genome_for_transcript(exons, _sequence(150, *overrides)),
        exons,
        ((0, 30),),
        limit=limit,
    )
    assert not result.success


def test_first_in_frame_downstream_stop_is_selected_and_residual_is_utr() -> None:
    exons = ((0, 90),)
    sequence = _sequence(90, (33, "TAG"), (42, "TGA"))
    result = _downstream_stop(
        _genome_for_transcript(exons, sequence), exons, ((0, 30),), limit=50
    )
    assert result.success and result.codon_positions == (33, 34, 35)
    assert result.cds_parts == ((0, 36),)
    assert result.utr_parts == ((36, 90),)


def test_downstream_stop_can_cross_a_splice_junction() -> None:
    exons = ((0, 31), (50, 80))
    sequence = _sequence(61, (30, "TGA"))
    result = _downstream_stop(
        _genome_for_transcript(exons, sequence), exons, ((0, 30),), limit=20
    )
    assert result.success
    assert result.codon_positions == (30, 50, 51)
    assert result.cds_parts == ((0, 31), (50, 52))
    assert result.splice_junctions_crossed == 1


@pytest.mark.parametrize(
    ("start_at", "stop_at", "expected"),
    [
        (24, 90, (True, True)),
        (24, None, (True, False)),
        (None, 90, (False, True)),
        (None, None, (False, False)),
    ],
)
def test_missing_both_sides_are_independently_completable(
    start_at: int | None,
    stop_at: int | None,
    expected: tuple[bool, bool],
) -> None:
    exons = ((0, 120),)
    overrides: list[tuple[int, str]] = []
    if start_at is not None:
        overrides.append((start_at, "ATG"))
    if stop_at is not None:
        overrides.append((stop_at, "TAA"))
    genome = _genome_for_transcript(exons, _sequence(120, *overrides))
    start = _historical_bounded_complete_start(
        genome,
        exons,
        ((30, 90),),
        historical_upstream_limit=50,
        downstream_limit=15,
    )
    start_parts = start.cds_parts if start.success else ((30, 90),)
    stop = _downstream_stop(genome, exons, start_parts, limit=50)
    assert (start.success, stop.success) == expected


def _ranking_candidate(
    model: CodingModel,
    *,
    evidence_id: str,
    structure_key: str,
    exon_parts: tuple[Interval, ...],
    source_count: int,
    extension_bp: int,
) -> UTRCandidateSide:
    evidence = UTREvidence(
        evidence_id=evidence_id,
        kind="transcript",
        chrom="chr1",
        strand="+",
        exon_parts=exon_parts,
        source_member_counts=tuple((f"source{index}", 1) for index in range(source_count)),
    )
    return UTRCandidateSide(
        gene_id="MoSelecta_000001",
        transcript_id="MoSelecta_000001.1",
        target=FinalTranscript("MoSelecta_000001.1", model),
        evidence=evidence,
        compatibility_class="characterization",
        side="five_prime",
        endpoint0=exon_parts[0][0],
        extension_interval=(exon_parts[0][0], 100),
        added_exons=exon_parts,
        added_introns=(),
        complete_contributing_exon_lengths=tuple(end - start for start, end in exon_parts),
        extension_exonic_bp=extension_bp,
        extension_genomic_bp=extension_bp,
        structure_key=structure_key,
        source_cds_class=".",
        source_cds_overlap_bp=0,
    )


def _ranking_model() -> CodingModel:
    return CodingModel(
        source_name="prediction",
        input_file=Path("prediction.gff3"),
        input_gene_id="gene",
        input_transcript_id="transcript",
        output_gene_id="out_gene",
        output_transcript_id="out_transcript",
        chrom="chr1",
        strand="+",
        cds_parts=[(100, 160)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def test_current_support_ranking_can_choose_structure_before_codon_search() -> None:
    model = _ranking_model()
    preferred = _ranking_candidate(
        model,
        evidence_id="preferred",
        structure_key="preferred_structure",
        exon_parts=((0, 160),),
        source_count=3,
        extension_bp=100,
    )
    weaker = _ranking_candidate(
        model,
        evidence_id="weaker",
        structure_key="weaker_structure",
        exon_parts=((40, 160),),
        source_count=1,
        extension_bp=60,
    )
    ranked, _strong_competing = _ranked_structures(
        [weaker, preferred], UTRSettings()
    )
    assert ranked[0].key == "preferred_structure"


def test_preferred_structure_can_fail_when_lower_ranked_structure_would_complete() -> None:
    model = _ranking_model()
    preferred = _ranking_candidate(
        model,
        evidence_id="preferred",
        structure_key="preferred_structure",
        exon_parts=((0, 160),),
        source_count=3,
        extension_bp=100,
    )
    lower = _ranking_candidate(
        model,
        evidence_id="lower",
        structure_key="lower_structure",
        exon_parts=((40, 160),),
        source_count=1,
        extension_bp=60,
    )
    ranked, _ = _ranked_structures([lower, preferred], UTRSettings())
    assert ranked[0].key == "preferred_structure"

    preferred_sequence = _sequence(160, (157, "TAA"))
    lower_sequence = _sequence(120, (54, "ATG"), (117, "TAA"))
    preferred_result = _historical_bounded_upstream_start(
        _genome_for_transcript(((0, 160),), preferred_sequence),
        ((0, 160),),
        ((100, 160),),
        historical_limit=100,
    )
    lower_result = _historical_bounded_upstream_start(
        _genome_for_transcript(((40, 160),), lower_sequence),
        ((40, 160),),
        ((100, 160),),
        historical_limit=100,
    )
    assert not preferred_result.success and lower_result.success


def test_different_splice_paths_can_produce_different_completed_cds_geometry() -> None:
    unspliced_exons = ((0, 90),)
    spliced_exons = ((0, 5), (20, 61))
    unspliced = _historical_bounded_upstream_start(
        _genome_for_transcript(unspliced_exons, _sequence(90, (24, "ATG"))),
        unspliced_exons,
        ((30, 90),),
        historical_limit=50,
    )
    spliced = _historical_bounded_upstream_start(
        _genome_for_transcript(spliced_exons, _sequence(46, (4, "ATG"))),
        spliced_exons,
        ((22, 61),),
        historical_limit=20,
    )
    assert unspliced.cds_parts == ((24, 90),)
    assert spliced.cds_parts == ((4, 5), (20, 61))


def test_apparent_utr_length_can_be_mostly_reassigned_to_cds() -> None:
    model = _ranking_model()
    candidate = _ranking_candidate(
        model,
        evidence_id="long-apparent-utr",
        structure_key="long_structure",
        exon_parts=((0, 160),),
        source_count=3,
        extension_bp=100,
    )
    limits = UTRLearnedLimits(
        five_prime_exonic_bp=50,
        three_prime_exonic_bp=50,
        five_prime_genomic_bp=None,
        three_prime_genomic_bp=None,
        five_prime_exonic_source="learned",
        three_prime_exonic_source="learned",
        five_prime_genomic_source="unavailable",
        three_prime_genomic_source="unavailable",
    )
    assert _extension_reasons(candidate, limits) == ("excessive_exonic_extension",)

    sequence = _sequence(160, (40, "ATG"), (157, "TAA"))
    completion = _historical_bounded_upstream_start(
        _genome_for_transcript(((0, 160),), sequence),
        ((0, 160),),
        ((100, 160),),
        historical_limit=100,
    )
    assert completion.success
    assert completion.cds_bases_delta == 60
    assert sum(end - start for start, end in completion.utr_parts) == 40


def test_minus_strand_spliced_coordinate_mapping_is_symmetric() -> None:
    exons = ((20, 61), (80, 85))
    sequence = _sequence(46, (4, "ATG"), (43, "TAA"))
    genome = _genome_for_transcript(exons, sequence, strand="-")
    result = _historical_bounded_upstream_start(
        genome,
        exons,
        ((20, 59),),
        historical_limit=20,
        strand="-",
    )
    assert result.success
    assert result.codon_positions == (80, 60, 59)
    assert result.cds_parts == ((20, 61), (80, 81))
