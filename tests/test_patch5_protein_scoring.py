from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.models import CodingModel
from moselecta.protein_evidence import (
    ProteinEvidenceCollapsedAlignment,
    pairwise_model_protein_metrics,
)


def make_model(parts: tuple[tuple[int, int], ...]) -> CodingModel:
    return CodingModel(
        source_name="model",
        input_file=Path("model.gff3"),
        input_gene_id="g1",
        input_transcript_id="m1",
        output_gene_id="g1",
        output_transcript_id="m1",
        chrom="chr1",
        strand="+",
        cds_parts=list(parts),
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def make_alignment(parts: tuple[tuple[int, int], ...]) -> ProteinEvidenceCollapsedAlignment:
    return ProteinEvidenceCollapsedAlignment(
        collapsed_id="p1",
        chrom="chr1",
        strand="+",
        cds_parts=parts,
        evidence_superlocus_id="psl1",
        source_member_counts={"protein": 1},
        member_count=1,
    )


@pytest.mark.parametrize(
    "model_parts,alignment_parts",
    [
        (((100, 150),), ((100, 150),)),
        (((100, 150),), ((100, 150), (250, 320))),
        (((100, 150), (250, 320)), ((100, 320),)),
    ],
)
def test_single_exon_on_either_side_uses_cds_base_f1_for_aggregate_score(
    model_parts: tuple[tuple[int, int], ...],
    alignment_parts: tuple[tuple[int, int], ...],
) -> None:
    metrics = pairwise_model_protein_metrics(
        make_model(model_parts),
        make_alignment(alignment_parts),
        wobble_bp=3,
    )
    assert metrics["score"] == pytest.approx(metrics["cds_base_f1"])


def test_mixed_single_multi_score_is_symmetric_and_terminal_diagnostics_remain_populated() -> None:
    multi_model_single_alignment = pairwise_model_protein_metrics(
        make_model(((100, 150), (250, 320))),
        make_alignment(((100, 320),)),
        wobble_bp=3,
    )
    single_model_multi_alignment = pairwise_model_protein_metrics(
        make_model(((100, 320),)),
        make_alignment(((100, 150), (250, 320))),
        wobble_bp=3,
    )

    assert multi_model_single_alignment["score"] == pytest.approx(
        single_model_multi_alignment["score"]
    )
    assert multi_model_single_alignment["score"] == pytest.approx(
        multi_model_single_alignment["cds_base_f1"]
    )
    assert multi_model_single_alignment["terminal_cds_exon_f1"] == 0.0
    assert single_model_multi_alignment["terminal_cds_exon_f1"] == 0.0


def test_multi_exon_pair_retains_weighted_aggregate_scoring() -> None:
    metrics = pairwise_model_protein_metrics(
        make_model(((100, 150), (250, 320))),
        make_alignment(((100, 150), (260, 320))),
        wobble_bp=3,
    )
    assert metrics["score"] != pytest.approx(metrics["cds_base_f1"])
