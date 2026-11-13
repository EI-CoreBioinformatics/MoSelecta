from __future__ import annotations

from itertools import product

import pytest

from moselecta import cds_phase
from moselecta.cds_phase import (
    IUPAC_BASES,
    STANDARD_START_CODONS,
    STANDARD_STOP_CODONS,
    analyse_cds_frame,
    codon_definitely_in,
    codon_membership,
    codon_possibly_in,
)


def _reference_membership(
    codon: str,
    targets: frozenset[str],
) -> tuple[bool, bool]:
    if len(codon) != 3:
        return False, False
    choices = [IUPAC_BASES.get(base.upper(), frozenset()) for base in codon]
    if any(not values for values in choices):
        return False, False
    expanded = {
        "".join(bases)
        for bases in product(*choices)
    }
    return bool(expanded & targets), bool(expanded) and expanded <= targets


@pytest.mark.parametrize(
    ("codon", "targets", "expected"),
    [
        ("ATG", STANDARD_START_CODONS, (True, True)),
        ("TAA", STANDARD_STOP_CODONS, (True, True)),
        ("TAG", STANDARD_STOP_CODONS, (True, True)),
        ("TGA", STANDARD_STOP_CODONS, (True, True)),
        ("GCC", STANDARD_STOP_CODONS, (False, False)),
        ("TAR", STANDARD_STOP_CODONS, (True, True)),
        ("TAN", STANDARD_STOP_CODONS, (True, False)),
        ("NCC", STANDARD_STOP_CODONS, (False, False)),
    ],
)
def test_codon_membership_distinguishes_possible_and_definite(
    codon: str,
    targets: frozenset[str],
    expected: tuple[bool, bool],
) -> None:
    assert codon_membership(codon, targets) == expected
    assert codon_possibly_in(codon, targets) is expected[0]
    assert codon_definitely_in(codon, targets) is expected[1]


def test_plain_acgt_fast_path_matches_general_iupac_semantics() -> None:
    for bases in product("ACGT", repeat=3):
        codon = "".join(bases)
        for targets in (STANDARD_START_CODONS, STANDARD_STOP_CODONS):
            assert codon_membership(codon, targets) == _reference_membership(
                codon, targets
            )


def test_ambiguity_cache_is_deterministic_and_plain_codons_bypass_it() -> None:
    cds_phase._cached_ambiguous_codon_membership.cache_clear()
    assert codon_membership("TAN", STANDARD_STOP_CODONS) == (True, False)
    after_first = cds_phase._cached_ambiguous_codon_membership.cache_info()
    assert codon_membership("TAN", STANDARD_STOP_CODONS) == (True, False)
    after_second = cds_phase._cached_ambiguous_codon_membership.cache_info()
    assert after_first.misses == 1
    assert after_second.misses == 1
    assert after_second.hits == after_first.hits + 1

    assert codon_membership("TAA", STANDARD_STOP_CODONS) == (True, True)
    after_plain = cds_phase._cached_ambiguous_codon_membership.cache_info()
    assert after_plain == after_second


def test_acgt_frame_fast_path_is_semantically_identical() -> None:
    sequence = "ATG" + "GCC" * 20 + "TAA"
    for phase in (0, 1, 2):
        assert analyse_cds_frame(
            sequence, "TGA", phase, cds_is_acgt=True
        ) == analyse_cds_frame(sequence, "TGA", phase)


def test_output_competition_distinguishes_included_from_adjacent_ambiguity() -> None:
    ambiguous_adjacent = analyse_cds_frame("ATGCCCCCCC", "TAN", 1)
    assert not ambiguous_adjacent.possible_included_stop
    assert not ambiguous_adjacent.literal_adjacent_stop
    assert not ambiguous_adjacent.supported_for_output_competition

    literal_adjacent = analyse_cds_frame("ATGCCCCCCC", "TAA", 1)
    assert literal_adjacent.literal_adjacent_stop
    assert literal_adjacent.supported_for_output_competition

    ambiguous_included = analyse_cds_frame("ATGCCCCTAN", "TAR", 1)
    assert ambiguous_included.possible_included_stop
    assert not ambiguous_included.literal_adjacent_stop
    assert ambiguous_included.supported_for_output_competition


def test_acgt_frame_analysis_never_expands_ambiguity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def unexpected_expansion(*_args: object, **_kwargs: object) -> tuple[bool, bool]:
        raise AssertionError("ordinary ACGT codon entered ambiguity expansion")

    monkeypatch.setattr(
        cds_phase, "_cached_ambiguous_codon_membership", unexpected_expansion
    )
    evidence = analyse_cds_frame(
        "ATG" + "GCC" * 20 + "TAA",
        "CCC",
        0,
        cds_is_acgt=True,
    )
    assert evidence.definite_start
    assert evidence.definite_stop
    assert evidence.definite_internal_stop_count == 0
