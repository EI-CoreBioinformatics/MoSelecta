from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.transcript_evidence import (
    build_transcript_subloci,
    build_transcript_superloci,
    collapse_transcript_evidence,
    parse_transcript_evidence,
    transcript_source_labels,
)


def write_gff(path: Path) -> None:
    path.write_text(
        "##gff-version 3\n"
        # tx1 and tx2 share exact intron chain but have different terminal outer boundaries.
        "chr1\tRNA\tmRNA\t100\t400\t.\t+\t.\tID=tx1\n"
        "chr1\tRNA\texon\t100\t150\t.\t+\t.\tParent=tx1\n"
        "chr1\tRNA\texon\t250\t400\t.\t+\t.\tParent=tx1\n"
        "chr1\tRNA\tmRNA\t90\t410\t.\t+\t.\tID=tx2\n"
        "chr1\tRNA\texon\t90\t150\t.\t+\t.\tParent=tx2\n"
        "chr1\tRNA\texon\t250\t410\t.\t+\t.\tParent=tx2\n"
        # tx3 is a distinct single-exon transcript.
        "chr1\tRNA\tmRNA\t1000\t1100\t.\t+\t.\tID=tx3\n"
        "chr1\tRNA\texon\t1000\t1100\t.\t+\t.\tParent=tx3\n"
    )


def test_transcript_evidence_exact_intron_chain_collapse(tmp_path: Path) -> None:
    gff = tmp_path / "tx.gff3"
    write_gff(gff)
    labels = transcript_source_labels([gff], None)
    assert labels == ["transcript_001"]

    raw = parse_transcript_evidence([gff], labels, strict=True)
    assert len(raw) == 3

    superloci = build_transcript_superloci(raw)
    assert len(superloci) == 2

    collapsed = collapse_transcript_evidence(superloci)
    assert len(collapsed) == 2
    multi = [tx for tx in collapsed if tx.is_multi_exon][0]
    assert multi.member_count == 2
    assert multi.collapse_key_type == "exact_intron_chain"
    assert multi.left_terminal_outer_min0 == 89
    assert multi.left_terminal_outer_max0 == 99
    assert multi.right_terminal_outer_min0 == 400
    assert multi.right_terminal_outer_max0 == 410

    subloci = build_transcript_subloci(superloci, single_exon_min_reciprocal_overlap=0.8)
    assert len(subloci) == 2


def test_transcript_evidence_accepts_unstranded_gtf_exons_under_strict(tmp_path: Path) -> None:
    gtf = tmp_path / "unstranded.gtf"
    gtf.write_text(
        'chr1\tStringTie\ttranscript\t100\t300\t1000\t.\t.\tgene_id "g1"; transcript_id "tx1";\n'
        'chr1\tStringTie\texon\t100\t150\t1000\t.\t.\tgene_id "g1"; transcript_id "tx1"; exon_number "1";\n'
        'chr1\tStringTie\texon\t250\t300\t1000\t.\t.\tgene_id "g1"; transcript_id "tx1"; exon_number "2";\n'
    )
    labels = transcript_source_labels([gtf], None)
    raw = parse_transcript_evidence([gtf], labels, strict=True)
    assert len(raw) == 1
    assert raw[0].strand == "."

    superloci = build_transcript_superloci(raw)
    assert len(superloci) == 1
    assert superloci[0].strand == "."
    assert "unstranded" in superloci[0].evidence_superlocus_id

    collapsed = collapse_transcript_evidence(superloci)
    assert len(collapsed) == 1
    assert collapsed[0].strand == "."
    assert collapsed[0].five_prime_terminal_outer_min0 is None
    assert collapsed[0].three_prime_terminal_outer_min0 is None


def test_transcript_model_support_separates_scoring_from_broad_association(tmp_path: Path) -> None:
    from moselecta.models import CodingModel
    from moselecta.transcript_evidence import calculate_transcript_model_support

    gff = tmp_path / "support.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tRNA\tmRNA\t100\t300\t.\t+\t.\tID=tx_spliced\n"
        "chr1\tRNA\texon\t100\t150\t.\t+\t.\tParent=tx_spliced\n"
        "chr1\tRNA\texon\t250\t300\t.\t+\t.\tParent=tx_spliced\n"
        "chr1\tRNA\tmRNA\t1000\t1300\t.\t+\t.\tID=tx_single\n"
        "chr1\tRNA\texon\t1000\t1300\t.\t+\t.\tParent=tx_single\n"
        "chr1\tRNA\tmRNA\t2000\t2750\t.\t+\t.\tID=tx_complex\n"
        "chr1\tRNA\texon\t2000\t2050\t.\t+\t.\tParent=tx_complex\n"
        "chr1\tRNA\texon\t2100\t2150\t.\t+\t.\tParent=tx_complex\n"
        "chr1\tRNA\texon\t2200\t2250\t.\t+\t.\tParent=tx_complex\n"
        "chr1\tRNA\texon\t2300\t2350\t.\t+\t.\tParent=tx_complex\n"
        "chr1\tRNA\texon\t2400\t2450\t.\t+\t.\tParent=tx_complex\n"
        "chr1\tRNA\texon\t2500\t2550\t.\t+\t.\tParent=tx_complex\n"
        "chr1\tRNA\texon\t2600\t2650\t.\t+\t.\tParent=tx_complex\n"
        "chr1\tRNA\texon\t2700\t2750\t.\t+\t.\tParent=tx_complex\n"
    )
    labels = transcript_source_labels([gff], None)
    raw = parse_transcript_evidence([gff], labels, strict=True)
    superloci = build_transcript_superloci(raw)
    collapse_transcript_evidence(superloci)
    build_transcript_subloci(superloci, single_exon_min_reciprocal_overlap=0.8)

    model_spliced = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g1", input_transcript_id="m1",
        output_gene_id="g1", output_transcript_id="m1", chrom="chr1", strand="+",
        cds_parts=[(99, 150), (249, 300)], gene_record=None, transcript_record=None, child_records=[],
    )
    model_single_inside_simple_multiexon = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g2", input_transcript_id="m2",
        output_gene_id="g2", output_transcript_id="m2", chrom="chr1", strand="+",
        cds_parts=[(110, 140)], gene_record=None, transcript_record=None, child_records=[],
    )
    model_single_inside_complex_multiexon = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g3", input_transcript_id="m3",
        output_gene_id="g3", output_transcript_id="m3", chrom="chr1", strand="+",
        cds_parts=[(2010, 2040)], gene_record=None, transcript_record=None, child_records=[],
    )
    rows = calculate_transcript_model_support([model_spliced, model_single_inside_simple_multiexon, model_single_inside_complex_multiexon], superloci, labels)
    by_id = {row["output_transcript_id"]: row for row in rows}

    assert by_id["m1"]["transcript_scoring_eligible_count"] == 1
    assert by_id["m1"]["transcript_best_scoring_association_class"] == "shared_intron"
    assert by_id["m2"]["transcript_broad_association_count"] == 1
    assert by_id["m2"]["transcript_scoring_eligible_count"] == 1
    assert by_id["m2"]["transcript_best_scoring_association_class"] == "single_exon_limited_multiexon_support"
    assert by_id["m2"]["transcript_single_exon_limited_multiexon_support_count"] == 1
    assert by_id["m3"]["transcript_broad_association_count"] == 1
    assert by_id["m3"]["transcript_scoring_eligible_count"] == 0
    assert by_id["m3"]["transcript_best_association_class"] == "single_exon_model_multiexon_transcript_overlap"


def test_transcript_model_support_chain_and_exon_compatibility(tmp_path: Path) -> None:
    from moselecta.models import CodingModel
    from moselecta.transcript_evidence import calculate_transcript_model_support

    gff = tmp_path / "chain.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # Exact CDS-compatible transcript with UTR-like terminal extension.
        "chr1\tRNA\tmRNA\t90\t500\t.\t+\t.\tID=tx_exact\n"
        "chr1\tRNA\texon\t90\t150\t.\t+\t.\tParent=tx_exact\n"
        "chr1\tRNA\texon\t250\t300\t.\t+\t.\tParent=tx_exact\n"
        "chr1\tRNA\texon\t400\t500\t.\t+\t.\tParent=tx_exact\n"
        # Longer transcript containing the model intron chain with one extra intron at each end.
        "chr1\tRNA\tmRNA\t50\t700\t.\t+\t.\tID=tx_long\n"
        "chr1\tRNA\texon\t50\t80\t.\t+\t.\tParent=tx_long\n"
        "chr1\tRNA\texon\t90\t150\t.\t+\t.\tParent=tx_long\n"
        "chr1\tRNA\texon\t250\t300\t.\t+\t.\tParent=tx_long\n"
        "chr1\tRNA\texon\t400\t500\t.\t+\t.\tParent=tx_long\n"
        "chr1\tRNA\texon\t650\t700\t.\t+\t.\tParent=tx_long\n"
    )
    labels = transcript_source_labels([gff], None)
    raw = parse_transcript_evidence([gff], labels, strict=True)
    superloci = build_transcript_superloci(raw)
    collapse_transcript_evidence(superloci)
    build_transcript_subloci(superloci, single_exon_min_reciprocal_overlap=0.8)

    model = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g1", input_transcript_id="m1",
        output_gene_id="g1", output_transcript_id="m1", chrom="chr1", strand="+",
        cds_parts=[(99, 150), (249, 300), (399, 450)], gene_record=None, transcript_record=None, child_records=[],
    )
    row = calculate_transcript_model_support([model], superloci, labels)[0]

    assert row["transcript_exact_intron_chain_match_count"] == 1
    assert row["transcript_model_intron_chain_contained_count"] == 2
    assert row["transcript_model_intron_chain_stringent_contained_count"] == 2
    assert row["transcript_best_intron_chain_f1"] == 1.0
    assert row["transcript_model_internal_cds_exon_count"] == 1
    assert row["transcript_supported_internal_cds_exon_count"] == 1
    assert row["transcript_model_terminal_cds_exon_count"] == 2
    assert row["transcript_supported_terminal_cds_exon_count"] == 2
    assert row["transcript_scoring_supported_terminal_cds_exon_count"] == 2


def test_transcript_cdna_context_diagnostics_for_best_compatible_transcript(tmp_path: Path) -> None:
    from moselecta.models import CodingModel
    from moselecta.transcript_evidence import calculate_transcript_model_support

    gff = tmp_path / "cdna_context.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # A multi-exon cDNA that can legitimately cover a single-exon CDS but
        # also provides cDNA-context diagnostics for fragment-risk review.
        "chr1\tRNA\tmRNA\t100\t300\t.\t+\t.\tID=tx_spliced\n"
        "chr1\tRNA\texon\t100\t150\t.\t+\t.\tParent=tx_spliced\n"
        "chr1\tRNA\texon\t250\t300\t.\t+\t.\tParent=tx_spliced\n"
        # A shared-splice-only transcript for a separate multi-exon model; this
        # should not get cDNA-context values because the chain is not contained.
        "chr1\tRNA\tmRNA\t1000\t1300\t.\t+\t.\tID=tx_shared_splice\n"
        "chr1\tRNA\texon\t1000\t1050\t.\t+\t.\tParent=tx_shared_splice\n"
        "chr1\tRNA\texon\t1150\t1200\t.\t+\t.\tParent=tx_shared_splice\n"
        "chr1\tRNA\texon\t1250\t1300\t.\t+\t.\tParent=tx_shared_splice\n"
    )
    labels = transcript_source_labels([gff], None)
    raw = parse_transcript_evidence([gff], labels, strict=True)
    superloci = build_transcript_superloci(raw)
    collapse_transcript_evidence(superloci)
    build_transcript_subloci(superloci, single_exon_min_reciprocal_overlap=0.8)

    single_cds = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g1", input_transcript_id="m1",
        output_gene_id="g1", output_transcript_id="single", chrom="chr1", strand="+",
        cds_parts=[(109, 139)], gene_record=None, transcript_record=None, child_records=[],
    )
    shared_splice_only = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g2", input_transcript_id="m2",
        output_gene_id="g2", output_transcript_id="shared", chrom="chr1", strand="+",
        cds_parts=[(999, 1050), (1159, 1180)], gene_record=None, transcript_record=None, child_records=[],
    )
    rows = calculate_transcript_model_support([single_cds, shared_splice_only], superloci, labels)
    by_id = {row["output_transcript_id"]: row for row in rows}

    assert by_id["single"]["transcript_cdna_context_status"] == "calculated"
    assert by_id["single"]["transcript_cdna_context_best_association_class"] == "single_exon_limited_multiexon_support"
    assert by_id["single"]["transcript_cdna_context_extra_3prime_exon_count"] == 1
    assert by_id["single"]["transcript_cdna_context_cds_to_min_cdna_fraction"] < 0.31
    assert "low_cds_to_cdna_fraction" in by_id["single"]["transcript_cdna_context_flags"]

    assert by_id["shared"]["transcript_best_scoring_association_class"] == "shared_splice_site"
    assert by_id["shared"]["transcript_cdna_context_status"] == "not_calculated_no_compatible_transcript"
    assert by_id["shared"]["transcript_cdna_context_best_evidence_id"] == "."


def test_transcript_cdna_context_requires_full_cds_exonic_containment(tmp_path: Path) -> None:
    from moselecta.models import CodingModel
    from moselecta.transcript_evidence import calculate_transcript_model_support

    gff = tmp_path / "cdna_context_fragment.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # The transcript has the same intron chain as the model but is a cDNA
        # fragment: its terminal exon extents do not cover the full CDS.  This
        # should not be used as a cDNA-context denominator.
        "chr1\tRNA\tmRNA\t120\t380\t.\t+\t.\tID=tx_fragment\n"
        "chr1\tRNA\texon\t120\t180\t.\t+\t.\tParent=tx_fragment\n"
        "chr1\tRNA\texon\t300\t380\t.\t+\t.\tParent=tx_fragment\n"
    )
    labels = transcript_source_labels([gff], None)
    raw = parse_transcript_evidence([gff], labels, strict=True)
    superloci = build_transcript_superloci(raw)
    collapse_transcript_evidence(superloci)
    build_transcript_subloci(superloci, single_exon_min_reciprocal_overlap=0.8)

    model = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g1", input_transcript_id="m1",
        output_gene_id="g1", output_transcript_id="fragment_case", chrom="chr1", strand="+",
        cds_parts=[(99, 180), (299, 450)], gene_record=None, transcript_record=None, child_records=[],
    )
    row = calculate_transcript_model_support([model], superloci, labels)[0]

    assert row["transcript_exact_intron_chain_match_count"] == 1
    assert row["transcript_cdna_context_status"] == "not_calculated_no_compatible_transcript"
    assert row["transcript_cdna_context_best_evidence_id"] == "."
    assert row["transcript_cdna_context_cds_to_min_cdna_fraction"] == "."



def _assert_calculated_cdna_context_invariants(row: dict[str, object]) -> None:
    if row["transcript_cdna_context_status"] != "calculated":
        return
    cds_length = int(row["transcript_cdna_context_model_cds_length_bp"])
    cdna_min = int(row["transcript_cdna_context_cdna_length_min_bp"])
    cdna_max = int(row["transcript_cdna_context_cdna_length_max_bp"])
    cds_to_min = float(row["transcript_cdna_context_cds_to_min_cdna_fraction"])
    cds_to_max = float(row["transcript_cdna_context_cds_to_max_cdna_fraction"])
    assert 0 < cds_length <= cdna_min <= cdna_max
    assert 0.0 < cds_to_max <= cds_to_min <= 1.0


def _calculate_one_model_row(gff: Path, model) -> tuple[dict[str, object], list, list]:
    from moselecta.transcript_evidence import calculate_transcript_model_support

    labels = transcript_source_labels([gff], None)
    raw = parse_transcript_evidence([gff], labels, strict=True)
    superloci = build_transcript_superloci(raw)
    collapse_transcript_evidence(superloci)
    build_transcript_subloci(superloci, single_exon_min_reciprocal_overlap=0.8)
    row = calculate_transcript_model_support([model], superloci, labels)[0]
    return row, raw, superloci


def test_transcript_cdna_context_falls_back_to_another_compatible_collapsed_transcript(tmp_path: Path) -> None:
    from moselecta.models import CodingModel

    gff = tmp_path / "cdna_context_fallback.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # Exact-chain fragment: expected to remain pairwise-best but must not be
        # selected for cDNA context because it misses both CDS termini.
        "chr1\tRNA\tmRNA\t120\t380\t.\t+\t.\tID=tx_fragment\n"
        "chr1\tRNA\texon\t120\t180\t.\t+\t.\tParent=tx_fragment\n"
        "chr1\tRNA\texon\t300\t380\t.\t+\t.\tParent=tx_fragment\n"
        # Longer-chain transcript containing and fully covering the model CDS.
        "chr1\tRNA\tmRNA\t50\t550\t.\t+\t.\tID=tx_compatible\n"
        "chr1\tRNA\texon\t50\t80\t.\t+\t.\tParent=tx_compatible\n"
        "chr1\tRNA\texon\t100\t180\t.\t+\t.\tParent=tx_compatible\n"
        "chr1\tRNA\texon\t300\t450\t.\t+\t.\tParent=tx_compatible\n"
        "chr1\tRNA\texon\t500\t550\t.\t+\t.\tParent=tx_compatible\n"
    )
    model = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g1", input_transcript_id="m1",
        output_gene_id="g1", output_transcript_id="fallback", chrom="chr1", strand="+",
        cds_parts=[(99, 180), (299, 450)], gene_record=None, transcript_record=None, child_records=[],
    )
    row, raw, _superloci = _calculate_one_model_row(gff, model)
    collapsed_by_original = {tx.original_transcript_id: tx.collapsed_id for tx in raw}

    assert row["transcript_best_pairwise_id"] == collapsed_by_original["tx_fragment"]
    assert row["transcript_cdna_context_status"] == "calculated"
    assert row["transcript_cdna_context_best_evidence_id"] == collapsed_by_original["tx_compatible"]
    assert row["transcript_cdna_context_best_matches_best_pairwise"] is False
    _assert_calculated_cdna_context_invariants(row)


def test_transcript_cdna_context_conditions_minimum_envelope_on_cds_containment(tmp_path: Path) -> None:
    from moselecta.models import CodingModel

    gff = tmp_path / "cdna_context_conditioned_min.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tRNA\tmRNA\t100\t450\t.\t+\t.\tID=tx_full\n"
        "chr1\tRNA\texon\t100\t200\t.\t+\t.\tParent=tx_full\n"
        "chr1\tRNA\texon\t300\t450\t.\t+\t.\tParent=tx_full\n"
        "chr1\tRNA\tmRNA\t150\t350\t.\t+\t.\tID=tx_short\n"
        "chr1\tRNA\texon\t150\t200\t.\t+\t.\tParent=tx_short\n"
        "chr1\tRNA\texon\t300\t350\t.\t+\t.\tParent=tx_short\n"
    )
    model = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g1", input_transcript_id="m1",
        output_gene_id="g1", output_transcript_id="conditioned", chrom="chr1", strand="+",
        cds_parts=[(99, 200), (299, 400)], gene_record=None, transcript_record=None, child_records=[],
    )
    row, _raw, _superloci = _calculate_one_model_row(gff, model)

    assert row["transcript_cdna_context_status"] == "calculated"
    assert row["transcript_cdna_context_cdna_length_min_bp"] == model.cds_length == 202
    assert row["transcript_cdna_context_cdna_length_max_bp"] == 252
    assert row["transcript_cdna_context_cds_to_min_cdna_fraction"] == 1.0
    assert "min_cdna_conditioned_on_cds_containment" in row["transcript_cdna_context_flags"]
    _assert_calculated_cdna_context_invariants(row)


def test_transcript_cdna_context_accepts_complementary_terminal_members(tmp_path: Path) -> None:
    from moselecta.models import CodingModel

    gff = tmp_path / "cdna_context_complementary.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # Neither member covers both CDS termini, but their exact-chain collapsed
        # outer envelope does. This is intentionally accepted by MoSelecta.
        "chr1\tRNA\tmRNA\t100\t350\t.\t+\t.\tID=tx_left\n"
        "chr1\tRNA\texon\t100\t200\t.\t+\t.\tParent=tx_left\n"
        "chr1\tRNA\texon\t300\t350\t.\t+\t.\tParent=tx_left\n"
        "chr1\tRNA\tmRNA\t150\t450\t.\t+\t.\tID=tx_right\n"
        "chr1\tRNA\texon\t150\t200\t.\t+\t.\tParent=tx_right\n"
        "chr1\tRNA\texon\t300\t450\t.\t+\t.\tParent=tx_right\n"
    )
    model = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g1", input_transcript_id="m1",
        output_gene_id="g1", output_transcript_id="complementary", chrom="chr1", strand="+",
        cds_parts=[(99, 200), (299, 450)], gene_record=None, transcript_record=None, child_records=[],
    )
    row, _raw, superloci = _calculate_one_model_row(gff, model)
    collapsed = [tx for sl in superloci for tx in sl.collapsed_transcripts]

    assert len(collapsed) == 1
    assert collapsed[0].member_count == 2
    assert row["transcript_cdna_context_status"] == "calculated"
    assert row["transcript_cdna_context_cdna_length_min_bp"] == model.cds_length == 252
    assert row["transcript_cdna_context_cdna_length_max_bp"] == 252
    assert row["transcript_cdna_context_cds_to_min_cdna_fraction"] == 1.0
    assert "min_cdna_conditioned_on_cds_containment" in row["transcript_cdna_context_flags"]
    _assert_calculated_cdna_context_invariants(row)


def test_transcript_cdna_context_conditioned_minimum_retains_obligatory_extra_exons(tmp_path: Path) -> None:
    from moselecta.models import CodingModel

    gff = tmp_path / "cdna_context_obligatory_extra.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tRNA\tmRNA\t100\t650\t.\t+\t.\tID=tx_full\n"
        "chr1\tRNA\texon\t100\t300\t.\t+\t.\tParent=tx_full\n"
        "chr1\tRNA\texon\t400\t500\t.\t+\t.\tParent=tx_full\n"
        "chr1\tRNA\texon\t600\t650\t.\t+\t.\tParent=tx_full\n"
        "chr1\tRNA\tmRNA\t290\t650\t.\t+\t.\tID=tx_short\n"
        "chr1\tRNA\texon\t290\t300\t.\t+\t.\tParent=tx_short\n"
        "chr1\tRNA\texon\t400\t500\t.\t+\t.\tParent=tx_short\n"
        "chr1\tRNA\texon\t600\t650\t.\t+\t.\tParent=tx_short\n"
    )
    model = CodingModel(
        source_name="A", input_file=gff, input_gene_id="g1", input_transcript_id="m1",
        output_gene_id="g1", output_transcript_id="obligatory", chrom="chr1", strand="+",
        cds_parts=[(99, 300), (399, 449)], gene_record=None, transcript_record=None, child_records=[],
    )
    row, _raw, _superloci = _calculate_one_model_row(gff, model)

    assert model.cds_length == 251
    # Scalar max(observed_min, CDS length) would be 251. Coordinate conditioning
    # retains the full matched second exon and the obligatory downstream exon.
    assert row["transcript_cdna_context_cdna_length_min_bp"] == 353
    assert row["transcript_cdna_context_cdna_length_max_bp"] == 353
    assert row["transcript_cdna_context_extra_3prime_exon_count"] == 1
    assert row["transcript_cdna_context_extra_3prime_exonic_bp"] == 102
    assert "min_cdna_conditioned_on_cds_containment" in row["transcript_cdna_context_flags"]
    _assert_calculated_cdna_context_invariants(row)


def test_transcript_cdna_context_minus_strand_has_same_lengths_and_swapped_terminal_orientation(tmp_path: Path) -> None:
    from moselecta.models import CodingModel
    from moselecta.transcript_evidence import calculate_transcript_model_support

    gff = tmp_path / "cdna_context_strands.gff3"
    records = ["##gff-version 3\n"]
    for chrom, strand, prefix in (("chrP", "+", "p"), ("chrM", "-", "m")):
        records.extend([
            f"{chrom}\tRNA\tmRNA\t100\t650\t.\t{strand}\t.\tID={prefix}_full\n",
            f"{chrom}\tRNA\texon\t100\t300\t.\t{strand}\t.\tParent={prefix}_full\n",
            f"{chrom}\tRNA\texon\t400\t500\t.\t{strand}\t.\tParent={prefix}_full\n",
            f"{chrom}\tRNA\texon\t600\t650\t.\t{strand}\t.\tParent={prefix}_full\n",
            f"{chrom}\tRNA\tmRNA\t290\t650\t.\t{strand}\t.\tID={prefix}_short\n",
            f"{chrom}\tRNA\texon\t290\t300\t.\t{strand}\t.\tParent={prefix}_short\n",
            f"{chrom}\tRNA\texon\t400\t500\t.\t{strand}\t.\tParent={prefix}_short\n",
            f"{chrom}\tRNA\texon\t600\t650\t.\t{strand}\t.\tParent={prefix}_short\n",
        ])
    gff.write_text("".join(records))
    models = [
        CodingModel(
            source_name="A", input_file=gff, input_gene_id="gp", input_transcript_id="mp",
            output_gene_id="gp", output_transcript_id="plus", chrom="chrP", strand="+",
            cds_parts=[(99, 300), (399, 449)], gene_record=None, transcript_record=None, child_records=[],
        ),
        CodingModel(
            source_name="A", input_file=gff, input_gene_id="gm", input_transcript_id="mm",
            output_gene_id="gm", output_transcript_id="minus", chrom="chrM", strand="-",
            cds_parts=[(99, 300), (399, 449)], gene_record=None, transcript_record=None, child_records=[],
        ),
    ]
    labels = transcript_source_labels([gff], None)
    raw = parse_transcript_evidence([gff], labels, strict=True)
    superloci = build_transcript_superloci(raw)
    collapse_transcript_evidence(superloci)
    build_transcript_subloci(superloci, single_exon_min_reciprocal_overlap=0.8)
    rows = calculate_transcript_model_support(models, superloci, labels)
    by_id = {row["output_transcript_id"]: row for row in rows}

    assert by_id["plus"]["transcript_cdna_context_cdna_length_min_bp"] == 353
    assert by_id["minus"]["transcript_cdna_context_cdna_length_min_bp"] == 353
    assert by_id["plus"]["transcript_cdna_context_extra_3prime_exon_count"] == 1
    assert by_id["plus"]["transcript_cdna_context_extra_5prime_exon_count"] == 0
    assert by_id["minus"]["transcript_cdna_context_extra_5prime_exon_count"] == 1
    assert by_id["minus"]["transcript_cdna_context_extra_3prime_exon_count"] == 0
    assert by_id["plus"]["transcript_cdna_context_extra_3prime_exonic_bp"] == 102
    assert by_id["minus"]["transcript_cdna_context_extra_5prime_exonic_bp"] == 102
    for row in rows:
        _assert_calculated_cdna_context_invariants(row)


def test_every_calculated_cdna_context_row_satisfies_fraction_invariants(tmp_path: Path) -> None:
    from moselecta.models import CodingModel
    from moselecta.transcript_evidence import calculate_transcript_model_support

    gff = tmp_path / "cdna_context_invariants.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tRNA\tmRNA\t100\t450\t.\t+\t.\tID=tx1\n"
        "chr1\tRNA\texon\t100\t200\t.\t+\t.\tParent=tx1\n"
        "chr1\tRNA\texon\t300\t450\t.\t+\t.\tParent=tx1\n"
        "chr1\tRNA\tmRNA\t150\t350\t.\t+\t.\tID=tx2\n"
        "chr1\tRNA\texon\t150\t200\t.\t+\t.\tParent=tx2\n"
        "chr1\tRNA\texon\t300\t350\t.\t+\t.\tParent=tx2\n"
        "chr1\tRNA\tmRNA\t1000\t1300\t.\t+\t.\tID=tx_single_context\n"
        "chr1\tRNA\texon\t1000\t1100\t.\t+\t.\tParent=tx_single_context\n"
        "chr1\tRNA\texon\t1200\t1300\t.\t+\t.\tParent=tx_single_context\n"
    )
    models = [
        CodingModel(
            source_name="A", input_file=gff, input_gene_id="g1", input_transcript_id="m1",
            output_gene_id="g1", output_transcript_id="multi", chrom="chr1", strand="+",
            cds_parts=[(99, 200), (299, 400)], gene_record=None, transcript_record=None, child_records=[],
        ),
        CodingModel(
            source_name="A", input_file=gff, input_gene_id="g2", input_transcript_id="m2",
            output_gene_id="g2", output_transcript_id="single", chrom="chr1", strand="+",
            cds_parts=[(1009, 1050)], gene_record=None, transcript_record=None, child_records=[],
        ),
    ]
    labels = transcript_source_labels([gff], None)
    raw = parse_transcript_evidence([gff], labels, strict=True)
    superloci = build_transcript_superloci(raw)
    collapse_transcript_evidence(superloci)
    build_transcript_subloci(superloci, single_exon_min_reciprocal_overlap=0.8)
    rows = calculate_transcript_model_support(models, superloci, labels)

    assert all(row["transcript_cdna_context_status"] == "calculated" for row in rows)
    for row in rows:
        _assert_calculated_cdna_context_invariants(row)


def test_strict_transcript_parsing_rejects_unparented_exons_and_permissive_mode_skips_them(tmp_path: Path) -> None:
    gff = tmp_path / "unparented.gff3"
    gff.write_text("##gff-version 3\nchr1\tRNA\texon\t100\t200\t.\t+\t.\t.\n")
    labels = transcript_source_labels([gff], None)

    with pytest.raises(ValueError, match="has no Parent, transcript_id, or ID"):
        parse_transcript_evidence([gff], labels, strict=True)
    assert parse_transcript_evidence([gff], labels, strict=False) == []


def test_collapsed_representative_raw_id_matches_representative_exon_chain(tmp_path: Path) -> None:
    gff = tmp_path / "representative.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # Lexicographically first ID is the uncommon longer structure.
        "chr1\tRNA\tmRNA\t90\t410\t.\t+\t.\tID=txA\n"
        "chr1\tRNA\texon\t90\t150\t.\t+\t.\tParent=txA\n"
        "chr1\tRNA\texon\t250\t410\t.\t+\t.\tParent=txA\n"
        "chr1\tRNA\tmRNA\t100\t400\t.\t+\t.\tID=txB\n"
        "chr1\tRNA\texon\t100\t150\t.\t+\t.\tParent=txB\n"
        "chr1\tRNA\texon\t250\t400\t.\t+\t.\tParent=txB\n"
        "chr1\tRNA\tmRNA\t100\t400\t.\t+\t.\tID=txC\n"
        "chr1\tRNA\texon\t100\t150\t.\t+\t.\tParent=txC\n"
        "chr1\tRNA\texon\t250\t400\t.\t+\t.\tParent=txC\n"
    )
    labels = transcript_source_labels([gff], None)
    raw = parse_transcript_evidence([gff], labels, strict=True)
    superloci = build_transcript_superloci(raw)
    collapsed = collapse_transcript_evidence(superloci)
    assert len(collapsed) == 1
    representative = collapsed[0]
    raw_by_id = {tx.raw_id: tx for tx in raw}

    assert representative.exon_parts == ((99, 150), (249, 400))
    assert raw_by_id[representative.representative_raw_id].exon_parts == representative.exon_parts
    assert raw_by_id[representative.representative_raw_id].original_transcript_id == "txB"
