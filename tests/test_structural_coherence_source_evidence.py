"""Production regressions for revised independent structural source evidence."""

import csv
from pathlib import Path

import pytest

from moselecta.models import CodingModel, Superlocus
from moselecta.structural_coherence import (
    _anchor_coverage,
    _bridge_compatible,
    _gene_identity,
    _source_region_facts,
    annotate_structural_coherence,
    write_structural_coherence_regions_tsv,
)
from moselecta import structural_coherence


A = (100, 200)
B = (500, 600)


def model(
    model_id: str, source: str, parts: tuple[tuple[int, int], ...],
    *, gene_id: str | None = None, transcript_id: str | None = None,
) -> CodingModel:
    item = CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"{model_id}.gene" if gene_id is None else gene_id,
        input_transcript_id=model_id if transcript_id is None else transcript_id,
        output_gene_id=f"{model_id}.gene.out",
        output_transcript_id=model_id,
        chrom="chr1", strand="+", cds_parts=list(parts),
        gene_record=None, transcript_record=None, child_records=[],
    )
    item.superlocus_id = "super1"
    item.sublocus_id = "super1.sublocus001"
    return item


def run(models: list[CodingModel]) -> tuple[dict[str, dict[str, object]], dict[str, list[dict[str, object]]]]:
    superlocus = Superlocus("super1", "chr1", "+", 0, 1000, models)
    rows = [{"output_transcript_id": item.output_transcript_id, "source_prior_score": item.source_name} for item in models]
    regions = annotate_structural_coherence([superlocus], rows, penalty_weight=0.30)
    by_model = {str(row["output_transcript_id"]): row for row in rows}
    by_candidate: dict[str, list[dict[str, object]]] = {item.output_transcript_id: [] for item in models}
    for row in regions:
        by_candidate[str(row["candidate_model_id"])].append(row)
    return by_model, by_candidate


def split_models(*sources: str, left: tuple[int, int] = A, right: tuple[int, int] = B) -> list[CodingModel]:
    return [item for source in sources for item in (
        model(f"{source}.left", source, (left,), gene_id=f"{source}.left.gene"),
        model(f"{source}.right", source, (right,), gene_id=f"{source}.right.gene"),
    )]


def test_anchor_coverage_cache_reuses_model_anchor_result_without_changing_value(monkeypatch: pytest.MonkeyPatch) -> None:
    first = model("first", "First", (A, B))
    second = model("second", "Second", ((150, 200), B))
    expected = {
        (id(item), anchor): _anchor_coverage(item, anchor)
        for item in (first, second) for anchor in (A, B)
    }
    calls = 0
    original = structural_coherence._overlap_bp

    def counted(*args: object) -> int:
        nonlocal calls
        calls += 1
        return original(*args)

    monkeypatch.setattr(structural_coherence, "_overlap_bp", counted)
    cache = {}
    for _ in range(3):
        for item in (first, second):
            for anchor in (A, B):
                assert _anchor_coverage(item, anchor, cache) == expected[(id(item), anchor)]
    assert calls == 4
    assert len(cache) == 4


def test_bombus_case1_exact_annevo_connector_and_same_gene_braker_alternatives() -> None:
    # Synthetic equivalent of NC_063271.1:14233066-14239691.
    models = [model("Tiberius.g1724.t1", "Tiberius", (A, B))]
    models.extend(split_models("Split"))
    models.extend([
        model("Braker.g1861.t1", "Braker", (A,), gene_id="Braker.g1861"),
        model("Braker.g1861.t2", "Braker", (B,), gene_id="Braker.g1861"),
        model("Annevo.connector", "Annevo", ((190, 200), (500, 510))),
    ])
    scores, regions = run(models)
    row = regions["Tiberius.g1724.t1"][0]
    assert row["bilateral_separated_source_count"] == 1
    assert row["independent_spanning_source_count"] == 1
    assert row["source_has_exact_connector_continuity"].find("Annevo:1") >= 0
    assert row["source_has_anchor_crossing"].find("Annevo:0") >= 0
    assert "Braker:1" in row["source_same_gene_alternative"]
    assert "Braker:0" in row["source_contributes_S"]
    assert scores["Tiberius.g1724.t1"]["structural_coherence_penalty_score"] == 0.0


def test_bombus_case2_identical_annevo_galba_chain_and_mixed_continuity() -> None:
    # Synthetic equivalent of NC_063271.1:6623174-6637958.
    models = [
        model("Annevo.g129.t1", "Annevo", (A, B)),
        model("Galba.g2021.t1", "Galba", (A, B), gene_id="Galba.g2021"),
        model("Galba.g2021.t2", "Galba", (A,), gene_id="Galba.g2021"),
        model("Galba.g2021.t3", "Galba", (B,), gene_id="Galba.g2021"),
    ] + split_models("Split1", "Split2")
    scores, regions = run(models)
    annevo = scores["Annevo.g129.t1"]
    galba = scores["Galba.g2021.t1"]
    assert annevo["structural_coherence_raw_score"] == pytest.approx(1 / 3)
    assert annevo["structural_coherence_raw_score"] == galba["structural_coherence_raw_score"]
    assert annevo["structural_coherence_penalty_score"] == galba["structural_coherence_penalty_score"]
    assert annevo["source_prior_score"] != galba["source_prior_score"]
    row = regions["Annevo.g129.t1"][0]
    assert "Galba:mixed" in row["source_state_label"]
    assert "Galba:1" in row["source_contributes_C"]
    assert "Galba:0" in row["source_contributes_S"]


def test_bombus_case3_nested_annevo_tiberius_bridge_is_symmetric() -> None:
    # Synthetic equivalent of NC_063275.1:3193795-3196356.
    wide_left = (100, 400)
    models = [
        model("Annevo.g432.t1", "Annevo", (wide_left, B)),
        model("Tiberius.g3682.t1", "Tiberius", (A, B)),
    ] + split_models("Split1", "Split2", left=wide_left)
    scores, regions = run(models)
    wide = regions["Annevo.g432.t1"][0]
    narrow = regions["Tiberius.g3682.t1"][0]
    assert wide["bridge_hypothesis_group_id"] == narrow["bridge_hypothesis_group_id"] != "."
    assert wide["bilateral_separated_source_count"] == narrow["bilateral_separated_source_count"] == 2
    assert wide["independent_spanning_source_count"] == narrow["independent_spanning_source_count"] == 1
    assert "Tiberius:1" in wide["source_has_bridge_group_continuity"]
    assert "Annevo:1" in narrow["source_has_bridge_group_continuity"]
    assert scores["Annevo.g432.t1"]["structural_coherence_penalty_score"] == pytest.approx(0.10)
    assert scores["Annevo.g432.t1"]["structural_coherence_penalty_score"] == scores["Tiberius.g3682.t1"]["structural_coherence_penalty_score"]


def test_exact_simple_connector_below_half_anchor_coverage_is_one_C_not_veto() -> None:
    models = [model("candidate", "Candidate", (A, B))] + split_models("Split1", "Split2")
    models.append(model("connector", "Connector", ((190, 200), (500, 510))))
    scores, regions = run(models)
    row = regions["candidate"][0]
    assert row["bilateral_separated_source_count"] == 2
    assert row["independent_spanning_source_count"] == 1
    assert row["source_C_reason"].find("Connector:exact_connector") >= 0
    assert scores["candidate"]["structural_coherence_penalty_score"] == pytest.approx(0.10)


def test_genuine_distinct_gene_separation_and_same_gene_suppression() -> None:
    models = [model("candidate", "Candidate", (A, B))]
    models.extend(split_models("Distinct"))
    models.extend([
        model("Same.t1", "Same", (A,), gene_id="Same.gene"),
        model("Same.t2", "Same", (B,), gene_id="Same.gene"),
    ])
    scores, regions = run(models)
    row = regions["candidate"][0]
    assert "Distinct:1" in row["source_contributes_S"]
    assert "Same:0" in row["source_contributes_S"]
    assert row["bilateral_separated_source_count"] == 1
    assert scores["candidate"]["structural_coherence_raw_score"] == 0.0


def test_unresolved_identity_suppresses_S_and_recognizable_transcript_roots_allow_it() -> None:
    models = [model("candidate", "Candidate", (A, B))] + split_models("Distinct")
    models.extend([
        model("Unknown.left", "Unknown", (A,), gene_id=".", transcript_id="alpha"),
        model("Unknown.right", "Unknown", (B,), gene_id=".", transcript_id="beta"),
    ])
    scores, regions = run(models)
    assert scores["candidate"]["structural_coherence_penalty_score"] == 0.0
    assert "Unknown:unresolved" in regions["candidate"][0]["source_gene_id_resolution_method"]
    assert "Unknown:0" in regions["candidate"][0]["source_contributes_S"]
    assert _gene_identity(models[-1]) == (None, "unresolved")
    fallback_models = [model("candidate", "Candidate", (A, B))] + split_models("Distinct")
    fallback_models.extend([
        model("Fallback.left", "Fallback", (A,), gene_id=".", transcript_id="leftGene.t1"),
        model("Fallback.right", "Fallback", (B,), gene_id=".", transcript_id="rightGene.t1"),
    ])
    scores, regions = run(fallback_models)
    assert scores["candidate"]["structural_coherence_penalty_score"] == pytest.approx(0.30)
    assert "Fallback:normalized_transcript_id" in regions["candidate"][0]["source_gene_id_resolution_method"]
    assert "Fallback:1" in regions["candidate"][0]["source_contributes_S"]


def test_same_source_identical_isoforms_count_once_and_receive_equal_penalty() -> None:
    models = [
        model("candidate.t1", "Candidate", (A, B), gene_id="Candidate.gene"),
        model("candidate.t2", "Candidate", (A, B), gene_id="Candidate.gene"),
    ] + split_models("Split1", "Split2")
    scores, regions = run(models)
    assert scores["candidate.t1"]["structural_coherence_raw_score"] == scores["candidate.t2"]["structural_coherence_raw_score"]
    assert scores["candidate.t1"]["structural_coherence_penalty_score"] == pytest.approx(0.30)
    assert regions["candidate.t1"][0]["bilateral_separated_source_count"] == 2
    assert regions["candidate.t1"][0]["independent_spanning_source_count"] == 0


def test_multifeature_connector_requires_all_introns_in_one_model() -> None:
    candidate = model("candidate", "Candidate", (A, (300, 350), B))
    region = {
        "left_anchor_start_0based": 100, "left_anchor_end_0based": 200,
        "right_anchor_start_0based": 500, "right_anchor_end_0based": 600,
        "first_cut_index": 1, "last_cut_index": 2,
    }
    partial = model("partial", "Partial", ((190, 200), (300, 350)))
    full = model("full", "Full", ((190, 200), (300, 350), (500, 510)))
    partial_facts = _source_region_facts(region, candidate, [partial], 0.50)
    full_facts = _source_region_facts(region, candidate, [full], 0.50)
    assert not partial_facts["has_exact_connector"]
    assert not partial_facts["direct_continuity"]
    assert full_facts["has_exact_connector"]
    assert not full_facts["has_exact_chain"]
    assert full_facts["direct_routes"] == ["exact_connector"]


def test_component_only_bridge_vote_is_distinguished_from_direct_compatibility() -> None:
    # A/B and B/C are compatible, but A/C have overlapping non-nested
    # connectors without a shared boundary. C crosses neither A nor B's
    # long right anchor, but supports its own hypothesis directly.
    models = [
        model("A", "Annevo", (A, (600, 900))),
        model("B", "Braker", ((100, 300), (600, 900))),
        model("C", "Galba", ((100, 300), (650, 750))),
    ] + split_models("Split1", "Split2", left=(100, 300), right=(600, 900))
    _scores, regions = run(models)
    a, b, c = regions["A"][0], regions["B"][0], regions["C"][0]
    assert a["bridge_hypothesis_group_id"] == b["bridge_hypothesis_group_id"] == c["bridge_hypothesis_group_id"]
    candidates = {item.output_transcript_id: item for item in models}
    prelim = {str(row["region_id"]): {"Split1", "Split2"} for row in (a, b, c)}
    assert _bridge_compatible(a, b, candidates, prelim)
    assert _bridge_compatible(b, c, candidates, prelim)
    assert not _bridge_compatible(a, c, candidates, prelim)
    assert "Galba:1" in a["source_contributes_C"]
    assert "Galba:0" in a["source_bridge_group_continuity_direct"]
    assert "Galba:1" in a["source_bridge_group_continuity_propagated"]


def test_region_tsv_serializes_source_votes_and_supporting_models(tmp_path: Path) -> None:
    models = [model("candidate", "Candidate", (A, B))] + split_models("Split1", "Split2")
    models.append(model("connector", "Connector", ((190, 200), (500, 510))))
    _scores, regions = run(models)
    path = tmp_path / "regions.tsv"
    write_structural_coherence_regions_tsv(regions["candidate"], path)
    with path.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert "Connector:1" in row["source_has_exact_connector_continuity"]
    assert "Connector:1" in row["source_contributes_C"]
    assert "Connector:exact_connector" in row["source_C_reason"]
    assert "Connector:left=.|right=.|crossing=.|chain=.|connector=connector" in row["source_evidence_model_ids"]
