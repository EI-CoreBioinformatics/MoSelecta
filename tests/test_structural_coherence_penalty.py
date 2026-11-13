from pathlib import Path

import pytest

from moselecta.constants import MODEL_SELECTION_DEFAULT_WEIGHTS
from moselecta.config import build_arg_parser, output_keys_for_profile
from moselecta.models import CodingModel, Superlocus
from moselecta.scoring import apply_model_selection_scores
from moselecta.structural_coherence import annotate_structural_coherence


A = (100, 200)
B = (500, 600)


def model(model_id: str, source: str, parts: tuple[tuple[int, int], ...]) -> CodingModel:
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"{model_id}.gene",
        input_transcript_id=model_id,
        output_gene_id=f"{model_id}.gene.out",
        output_transcript_id=model_id,
        chrom="chr1",
        strand="+",
        cds_parts=list(parts),
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def run(models: list[CodingModel], candidate_id: str = "candidate", weight: float = 0.30):
    superlocus = Superlocus("super1", "chr1", "+", 100, 600, models)
    rows = [{"output_transcript_id": item.output_transcript_id} for item in models]
    regions = annotate_structural_coherence([superlocus], rows, penalty_weight=weight)
    return {str(row["output_transcript_id"]): row for row in rows}[candidate_id], [
        row for row in regions if row["candidate_model_id"] == candidate_id
    ]


def ratio_models(separated: int, spanning: int) -> list[CodingModel]:
    models = [model("candidate", "Candidate", (A, B))]
    for index in range(separated):
        source = f"Split{index + 1}"
        models.extend([
            model(f"{source}.left", source, (A,)),
            model(f"{source}.right", source, (B,)),
        ])
    for index in range(spanning):
        source = f"Span{index + 1}"
        models.append(model(f"{source}.cross", source, (A, B)))
    return models


@pytest.mark.parametrize(
    ("separated", "spanning", "expected_raw"),
    [
        (2, 0, 1.0),
        (2, 1, 1.0 / 3.0),
        (3, 0, 1.1),
        (3, 1, 0.55),
        (4, 0, 1.2),
        (4, 2, 0.4),
    ],
)
def test_requested_source_ratio_cases(separated: int, spanning: int, expected_raw: float) -> None:
    summary, regions = run(ratio_models(separated, spanning))
    eligible = [row for row in regions if row["penalty_gate_passed"] == "true"]

    assert len(eligible) == 1
    assert eligible[0]["bilateral_separated_source_count"] == separated
    assert eligible[0]["independent_spanning_source_count"] == spanning
    assert summary["structural_coherence_raw_score"] == pytest.approx(expected_raw)
    assert summary["structural_coherence_penalty_score"] == pytest.approx(0.30 * expected_raw)


def test_duplicate_isoforms_from_one_source_count_once() -> None:
    models = ratio_models(2, 0)
    models.extend([
        model("Split1.left.duplicate", "Split1", (A,)),
        model("Split1.right.duplicate", "Split1", (B,)),
    ])
    summary, regions = run(models)

    assert summary["structural_coherence_raw_score"] == pytest.approx(1.0)
    assert regions[0]["bilateral_separated_source_count"] == 2


def test_mixed_source_is_audited_and_contributes_one_continuity_vote() -> None:
    models = ratio_models(2, 0)
    models.extend([
        model("mixed.cross", "Mixed", (A, B)),
        model("mixed.left", "Mixed", (A,)),
        model("mixed.right", "Mixed", (B,)),
    ])
    summary, regions = run(models)
    region = regions[0]

    assert summary["structural_coherence_raw_score"] == pytest.approx(1.0 / 3.0)
    assert region["bilateral_separated_source_count"] == 2
    assert region["independent_spanning_source_count"] == 1
    assert region["mixed_source_count"] == 1
    assert region["mixed_sources"] == "Mixed"
    assert "Mixed:1" in region["source_contributes_C"]
    assert "Mixed:0" in region["source_contributes_S"]


def test_one_sided_support_fails_anchor_gate() -> None:
    models = [
        model("candidate", "Candidate", (A, B)),
        model("left1", "Left1", (A,)),
        model("left2", "Left2", (A,)),
        model("right1", "Right1", (B,)),
    ]
    summary, regions = run(models)

    assert summary["structural_coherence_raw_score"] == 0.0
    assert summary["structural_coherence_region_count"] == 0
    assert regions[0]["penalty_gate_passed"] == "false"
    assert "right_anchor>=2" in str(regions[0]["penalty_gate_reason"])


def test_equal_split_and_span_counts_fail_comparative_gate() -> None:
    summary, regions = run(ratio_models(2, 2))

    assert summary["structural_coherence_raw_score"] == 0.0
    assert summary["structural_coherence_penalty_score"] == 0.0
    assert regions[0]["penalty_gate_passed"] == "false"
    assert "S>C" in str(regions[0]["penalty_gate_reason"])


def test_different_splice_long_model_counts_as_independent_span() -> None:
    models = ratio_models(2, 0)
    models.append(model("alternate.cross", "Alternate", ((90, 210), (490, 610))))
    summary, regions = run(models)
    region = regions[0]

    assert summary["structural_coherence_raw_score"] == pytest.approx(1.0 / 3.0)
    assert region["independent_spanning_source_count"] == 1
    assert region["connector_intron_independent_exact_support_counts"] == "0"


def test_multiple_bridge_regions_use_maximum_not_sum() -> None:
    c = (900, 1000)
    d = (1300, 1400)
    models = [model("candidate", "Candidate", (A, B, c, d))]
    for source in ("AB1", "AB2"):
        models.extend([model(f"{source}.left", source, (A,)), model(f"{source}.right", source, (B,))])
    for source in ("CD1", "CD2", "CD3"):
        models.extend([model(f"{source}.left", source, (c,)), model(f"{source}.right", source, (d,))])
    models.append(model("middle.cross", "Middle", (B, c)))

    summary, regions = run(models)
    eligible = [row for row in regions if row["penalty_gate_passed"] == "true"]

    assert len(eligible) == 2
    assert summary["structural_coherence_region_count"] == 2
    assert summary["structural_coherence_raw_score"] == pytest.approx(1.1)
    assert summary["structural_coherence_penalty_score"] == pytest.approx(0.33)


def test_breadth_factor_caps_at_one_point_five() -> None:
    summary, regions = run(ratio_models(8, 0))

    assert regions[0]["breadth_factor"] == pytest.approx(1.5)
    assert summary["structural_coherence_raw_score"] == pytest.approx(1.5)
    assert summary["structural_coherence_penalty_score"] == pytest.approx(0.45)


def test_candidate_specific_internal_block_merges_into_one_outer_anchor_region() -> None:
    internal = (300, 350)
    models = [model("candidate", "Candidate", (A, internal, B))]
    for source in ("Split1", "Split2"):
        models.extend([
            model(f"{source}.left", source, (A,)),
            model(f"{source}.right", source, (B,)),
        ])

    summary, regions = run(models)

    assert len(regions) == 1
    assert regions[0]["boundaries_merged"] == "true"
    assert regions[0]["number_of_candidate_boundaries_in_region"] == 2
    assert regions[0]["candidate_only_connector_internal_cds_exon_count"] == 1
    assert summary["structural_coherence_raw_score"] == pytest.approx(1.0)


def test_configuration_default_and_output_profile() -> None:
    parser = build_arg_parser(full_help=True)
    args = parser.parse_args(["--inputs", "input.gff3"])

    assert args.selection_penalty_structural_coherence == pytest.approx(0.30)
    assert args.structural_coherence_regions_tsv is None
    assert "structural_coherence_regions_tsv" not in output_keys_for_profile("standard")
    assert "structural_coherence_regions_tsv" in output_keys_for_profile("full")
    assert "structural_coherence_regions_tsv" not in output_keys_for_profile("minimal")


def test_model_score_subtracts_configured_weight_times_raw_score() -> None:
    weights = dict(MODEL_SELECTION_DEFAULT_WEIGHTS)
    base = {
        "output_transcript_id": "candidate",
        "source_name": "Candidate",
        "sublocus_id": "sub1",
        "structural_coherence_region_count": 1,
        "structural_coherence_max_region": "candidate.coherence01",
    }
    unpenalized = dict(base, structural_coherence_raw_score=0.0)
    penalized = dict(base, structural_coherence_raw_score=1.0)
    rows = [unpenalized, penalized]

    apply_model_selection_scores(
        rows,
        source_prior_scores={},
        selection_weights=weights,
        tie_epsilon=1e-12,
    )

    assert unpenalized["model_structural_support_score"] - penalized["model_structural_support_score"] == pytest.approx(0.30)
    assert penalized["structural_coherence_penalty_score"] == pytest.approx(0.30)
    assert "-structural_coherence_raw_score:1.000000*0.300" in str(
        penalized["model_structural_support_score_components"]
    )
