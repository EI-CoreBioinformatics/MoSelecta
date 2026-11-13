from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import moselecta.utr_augmentation as utr_augmentation_module

from moselecta.annotation_outputs import (
    FinalAnnotationOutputMetadata,
    inherited_utr_validations,
    write_final_annotation_gff3,
)
from moselecta.genome_sequence import IndexedFasta
from moselecta.final_annotation import (
    AlternativeTranscriptMetadata,
    AnnotationSummary,
    FinalAnnotation,
    FinalGene,
    FinalTranscript,
    public_gene_id,
    public_transcript_id,
)
from moselecta.gffio import parse_attributes
from moselecta.models import CodingModel, FeatureRecord
from moselecta.utr_augmentation import (
    PooledSupport,
    FrozenFinalGeneStructureIndex,
    UTRCandidateSide,
    UTREvidence,
    UTREvidenceIndex,
    UTRSideSelection,
    UTRSettings,
    UTRStructure,
    _SelectedPopulationStructureIndex,
    _bounded_windows,
    _limits_from_observations,
    _select_endpoint,
    _select_sides_for_transcript,
    _selected_global_conflicts,
    _structural_features,
    _structure_source_count,
    augment_final_annotation,
    build_prediction_utr_catalogue,
    candidate_sides_for_transcript,
    learn_extension_limits,
    transcript_source_weight,
    weighted_observed_medoid,
    write_utr_augmented_gff3,
)
from moselecta.utr_validation import UTRValidationRecorder


def _record(
    feature_type: str,
    start0: int,
    end0: int,
    *,
    strand: str = "+",
    phase: str = ".",
    parent: str = "tx",
    source: str = "prediction",
) -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1",
        source=source,
        feature_type=feature_type,
        start=start0 + 1,
        end=end0,
        score=".",
        strand=strand,
        phase=phase,
        attrs={"Parent": parent},
        raw_attrs=".",
        input_file=Path("input.gff3"),
        source_name=source,
        line_number=1,
        record_index=1,
        raw_line="",
    )


def _model(
    name: str,
    *,
    parts: list[tuple[int, int]] | None = None,
    strand: str = "+",
    source: str = "prediction",
    exons: list[tuple[int, int]] | None = None,
    phases: list[str] | None = None,
) -> CodingModel:
    cds = parts or [(100, 150), (200, 250)]
    children: list[FeatureRecord] = []
    for start0, end0 in exons or []:
        children.append(_record("exon", start0, end0, strand=strand, source=source))
    for index, (start0, end0) in enumerate(cds):
        children.append(_record(
            "CDS", start0, end0, strand=strand,
            phase=(phases or ["0"] * len(cds))[index], source=source,
        ))
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"g_{name}",
        input_transcript_id=f"t_{name}",
        output_gene_id=f"og_{name}",
        output_transcript_id=f"ot_{name}",
        chrom="chr1",
        strand=strand,
        cds_parts=cds,
        gene_record=None,
        transcript_record=None,
        child_records=children,
    )


def _annotation(*models: CodingModel) -> FinalAnnotation:
    genes = []
    for gene_index, model in enumerate(models, start=1):
        gene_id = public_gene_id(gene_index)
        genes.append(FinalGene(gene_id, (FinalTranscript(public_transcript_id(gene_id, 1), model),)))
    return FinalAnnotation(tuple(genes))


def _multi_isoform_annotation(primary: CodingModel, alternative: CodingModel) -> FinalAnnotation:
    gene_id = "MoSelecta_000001"
    return FinalAnnotation((FinalGene(gene_id, (
        FinalTranscript(public_transcript_id(gene_id, 1), primary),
        FinalTranscript(public_transcript_id(gene_id, 2), alternative, transcript_role="alternative"),
    )),))


def _tx(
    evidence_id: str,
    exons: tuple[tuple[int, int], ...],
    *,
    strand: str = "+",
    counts: tuple[tuple[str, int], ...] = (("rna", 1),),
) -> UTREvidence:
    return UTREvidence(
        evidence_id=evidence_id,
        kind="transcript",
        chrom="chr1",
        strand=strand,
        exon_parts=exons,
        source_member_counts=counts,
    )


def _pred(
    evidence_id: str,
    exons: tuple[tuple[int, int], ...],
    *,
    source: str = "pred",
    model_id: str | None = None,
    strand: str = "+",
    cds: tuple[tuple[int, int], ...] = (),
    selected: bool = False,
) -> UTREvidence:
    return UTREvidence(
        evidence_id=evidence_id,
        kind="prediction",
        chrom="chr1",
        strand=strand,
        exon_parts=exons,
        prediction_source=source,
        prediction_model_id=model_id or evidence_id,
        prediction_cds_parts=cds,
        prediction_selected=selected,
    )


def _settings(**overrides: object) -> UTRSettings:
    values: dict[str, object] = {
        "five_prime_exonic_override_bp": 100_000,
        "three_prime_exonic_override_bp": 100_000,
        "five_prime_genomic_override_bp": 100_000,
        "three_prime_genomic_override_bp": 100_000,
    }
    values.update(overrides)
    return UTRSettings(**values)


def _result(
    model: CodingModel,
    evidence: list[UTREvidence],
    *,
    settings: UTRSettings | None = None,
):
    transcripts = [item for item in evidence if item.kind == "transcript"]
    predictions = [item for item in evidence if item.kind == "prediction"]
    return augment_final_annotation(
        _annotation(model), transcripts, predictions, settings or _settings()
    )


def _write_sequence_fasta(tmp_path: Path, edits: dict[tuple[int, int], str]) -> Path:
    sequence = list("A" * 400)
    for (start0, end0), value in edits.items():
        sequence[start0:end0] = value
    path = tmp_path / "genome.fa"
    path.write_text(">chr1\n" + "".join(sequence) + "\n")
    return path


def test_inherited_rejection_is_side_specific_and_preserves_cds(tmp_path: Path) -> None:
    model = _model(
        "inherited",
        exons=[(50, 70), (80, 150), (200, 270)],
    )
    annotation = _annotation(model)
    validations = inherited_utr_validations(annotation, minimum_intron_bp=20)
    by_side = validations["MoSelecta_000001.1"]
    assert by_side["five_prime"].status == "rejected"
    assert by_side["three_prime"].status == "passed"

    output = tmp_path / "final.gff3"
    write_final_annotation_gff3(
        annotation,
        output,
        FinalAnnotationOutputMetadata(False, False, False),
        utr_validations=validations,
    )
    lines = [line.split("\t") for line in output.read_text().splitlines() if not line.startswith("#")]
    exons = {(int(row[3]), int(row[4])) for row in lines if row[2] == "exon"}
    cds = {(int(row[3]), int(row[4])) for row in lines if row[2] == "CDS"}
    assert exons == {(101, 150), (201, 270)}
    assert cds == {(101, 150), (201, 250)}
    assert "moselecta_five_prime_utr_validation=rejected" in output.read_text()


def test_inherited_explicit_utr_records_without_exons_are_validated(
    tmp_path: Path,
) -> None:
    model = _model("explicit_utr")
    model.child_records.extend([
        _record("five_prime_UTR", 50, 70),
        _record("five_prime_UTR", 80, 100),
    ])
    annotation = _annotation(model)
    validations = inherited_utr_validations(annotation, minimum_intron_bp=20)
    five_prime = validations["MoSelecta_000001.1"]["five_prime"]
    assert five_prime.status == "rejected"
    assert "short_utr_intron" in five_prime.reason_codes

    output = tmp_path / "explicit_utr.gff3"
    write_final_annotation_gff3(
        annotation,
        output,
        FinalAnnotationOutputMetadata(False, False, False),
        utr_validations=validations,
    )
    rows = [
        row.split("\t")
        for row in output.read_text().splitlines()
        if not row.startswith("#")
    ]
    assert not any(
        row[2].lower() in {"utr", "five_prime_utr"} and int(row[4]) <= 100
        for row in rows
    )
    assert {(int(row[3]), int(row[4])) for row in rows if row[2] == "CDS"} == {
        (101, 150),
        (201, 250),
    }


def test_transcript_candidate_motif_rejected_before_selection(tmp_path: Path) -> None:
    model = _model("target")
    evidence = _tx("rna_bad", ((20, 50), (100, 150), (200, 250)))
    validations = []
    fasta_path = _write_sequence_fasta(
        tmp_path, {(50, 52): "GG", (98, 100): "TT"}
    )
    with IndexedFasta(fasta_path) as fasta:
        result = augment_final_annotation(
            _annotation(model), [evidence], [], _settings(), fasta=fasta,
            validation_recorder=validations.append,
        )
    transcript = result.annotation.genes[0].transcripts[0]
    assert transcript.five_prime_selection is None
    rejected = [
        validation for validation in validations
        if validation.evidence_id == "rna_bad" and validation.side == "five_prime"
    ]
    assert len(rejected) == 1
    assert rejected[0].evidence_origin == "transcript"
    assert rejected[0].status == "rejected"
    assert rejected[0].reason_codes == ("resolved_unrecognized_motif",)
    assert dict(result.summary.rejection_counts)["resolved_unrecognized_motif"] >= 1


def test_validation_detail_profile_does_not_change_utr_decisions(
    tmp_path: Path,
) -> None:
    model = _model("target")
    evidence = [
        _tx("bad", ((20, 50), (100, 150), (200, 250))),
        _tx("good", ((80, 150), (200, 270))),
    ]
    fasta_path = _write_sequence_fasta(
        tmp_path,
        {(50, 52): "GG", (98, 100): "TT"},
    )
    standard = UTRValidationRecorder()
    full = UTRValidationRecorder(tmp_path / "exceptions.tsv")
    with IndexedFasta(fasta_path) as fasta:
        standard_result = augment_final_annotation(
            _annotation(model), evidence, [], _settings(), fasta=fasta,
            validation_recorder=lambda value: standard.record(
                value, rejected_action="exclude_candidate_side"
            ),
        )
    with IndexedFasta(fasta_path) as fasta:
        full_result = augment_final_annotation(
            _annotation(model), evidence, [], _settings(), fasta=fasta,
            validation_recorder=lambda value: full.record(
                value, rejected_action="exclude_candidate_side"
            ),
        )
    standard.close()
    full.close()

    standard_tx = standard_result.annotation.genes[0].transcripts[0]
    full_tx = full_result.annotation.genes[0].transcripts[0]
    assert standard_tx.exon_parts == full_tx.exon_parts
    assert standard_tx.five_prime_utr_parts == full_tx.five_prime_utr_parts
    assert standard_tx.three_prime_utr_parts == full_tx.three_prime_utr_parts
    assert standard.summary_rows() == full.summary_rows()
    assert full.detail_rows_written > 0


@pytest.mark.parametrize(
    ("strand", "exons", "five", "three"),
    [
        ("+", ((80, 150), (200, 250)), True, False),
        ("+", ((100, 150), (200, 270)), False, True),
        ("+", ((80, 150), (200, 270)), True, True),
        ("-", ((100, 150), (200, 270)), True, False),
        ("-", ((80, 150), (200, 250)), False, True),
        ("-", ((80, 150), (200, 270)), True, True),
    ],
)
def test_biological_sides_are_strand_aware(
    strand: str,
    exons: tuple[tuple[int, int], ...],
    five: bool,
    three: bool,
) -> None:
    model = _model("target", strand=strand)
    transcript = _annotation(model).genes[0].transcripts[0]
    candidates = candidate_sides_for_transcript(
        "MoSelecta_000001", transcript,
        UTREvidenceIndex([_tx("e", exons, strand=strand)]), _settings(),
    )
    sides = {candidate.side for candidate in candidates}
    assert ("five_prime" in sides) is five
    assert ("three_prime" in sides) is three


def test_compatibility_classes_and_plain_overlap_rejection() -> None:
    multi = _model("multi")
    transcript = _annotation(multi).genes[0].transcripts[0]
    evidence = [
        _tx("full", ((80, 150), (200, 270))),
        _tx("terminal", ((80, 150), (200, 220))),
        _tx("overlap_only", ((90, 140), (210, 260))),
    ]
    candidates = candidate_sides_for_transcript(
        "MoSelecta_000001", transcript, UTREvidenceIndex(evidence), _settings()
    )
    modes = {(item.evidence.evidence_id, item.side): item.compatibility_class for item in candidates}
    assert modes[("full", "five_prime")] == "full_containment"
    assert modes[("terminal", "five_prime")] == "terminal_compatibility"
    assert not any(key[0] == "overlap_only" for key in modes)

    single = _model("single", parts=[(100, 200)])
    single_tx = _annotation(single).genes[0].transcripts[0]
    single_candidates = candidate_sides_for_transcript(
        "MoSelecta_000001", single_tx,
        UTREvidenceIndex([_tx("strict", ((80, 220),))]), _settings(),
    )
    assert {item.compatibility_class for item in single_candidates} == {
        "single_exon_containment"
    }


def test_complete_exon_filter_uses_source_exon_not_clipped_utr() -> None:
    model = _model("single", parts=[(100, 396)])
    transcript = _annotation(model).genes[0].transcripts[0]
    passing = candidate_sides_for_transcript(
        "MoSelecta_000001", transcript,
        UTREvidenceIndex([_tx("large", ((96, 400),))]), _settings(),
    )
    assert passing
    assert all(candidate.extension_exonic_bp == 4 for candidate in passing)
    assert all(not candidate.rejection_reasons for candidate in passing)

    multi = _model("multi")
    multi_tx = _annotation(multi).genes[0].transcripts[0]
    failing = candidate_sides_for_transcript(
        "MoSelecta_000001", multi_tx,
        UTREvidenceIndex([_tx("tiny", ((80, 84), (100, 150), (200, 250)))]),
        _settings(),
    )
    five = next(candidate for candidate in failing if candidate.side == "five_prime")
    assert "tiny_terminal_complete_evidence_exon" in five.rejection_reasons

    internal = candidate_sides_for_transcript(
        "MoSelecta_000001", multi_tx,
        UTREvidenceIndex([_tx(
            "tiny_internal", ((60, 70), (80, 84), (100, 150), (200, 250))
        )]),
        _settings(min_intron_bp=None),
    )
    internal_five = next(
        candidate for candidate in internal if candidate.side == "five_prime"
    )
    assert "tiny_internal_complete_evidence_exon" in internal_five.rejection_reasons


@pytest.mark.parametrize(
    ("exons", "setting", "reason"),
    [
        (((80, 90), (100, 150), (200, 250)), {"min_intron_bp": 20}, "short_utr_intron"),
        (((20, 40), (100, 150), (200, 250)), {"max_intron_bp": 50}, "long_utr_intron_moselecta"),
    ],
)
def test_new_utr_intron_safeguards(
    exons: tuple[tuple[int, int], ...], setting: dict[str, int], reason: str
) -> None:
    model = _model("target")
    transcript = _annotation(model).genes[0].transcripts[0]
    candidates = candidate_sides_for_transcript(
        "MoSelecta_000001", transcript,
        UTREvidenceIndex([_tx("e", exons)]), _settings(**setting),
    )
    five = next(candidate for candidate in candidates if candidate.side == "five_prime")
    assert reason in five.rejection_reasons


def test_utr_exon_complexity_limits_are_side_specific() -> None:
    model = _model("target")
    transcript = _annotation(model).genes[0].transcripts[0]
    accepted_five = _tx("f4", ((0, 10), (20, 30), (40, 50), (60, 150), (200, 250)))
    rejected_five = _tx("f5", ((0, 10), (20, 30), (40, 50), (60, 70), (80, 150), (200, 250)))
    accepted_three = _tx("t2", ((100, 150), (200, 260), (270, 280)))
    rejected_three = _tx("t3", ((100, 150), (200, 260), (270, 280), (290, 300)))
    candidates = candidate_sides_for_transcript(
        "MoSelecta_000001", transcript,
        UTREvidenceIndex([accepted_five, rejected_five, accepted_three, rejected_three]),
        _settings(min_intron_bp=None),
    )
    by_id_side = {(c.evidence.evidence_id, c.side): c for c in candidates}
    assert not by_id_side[("f4", "five_prime")].rejection_reasons
    assert "excessive_five_prime_utr_exon_count" in by_id_side[("f5", "five_prime")].rejection_reasons
    assert not by_id_side[("t2", "three_prime")].rejection_reasons
    assert "excessive_three_prime_utr_exon_count" in by_id_side[("t3", "three_prime")].rejection_reasons


def test_same_topology_endpoints_group_and_splice_topology_separates() -> None:
    model = _model("target")
    transcript = _annotation(model).genes[0].transcripts[0]
    evidence = [
        _tx("a", ((80, 150), (200, 250))),
        _tx("b", ((70, 150), (200, 250))),
        _tx("c", ((60, 70), (80, 150), (200, 250))),
    ]
    candidates = candidate_sides_for_transcript(
        "MoSelecta_000001", transcript, UTREvidenceIndex(evidence),
        _settings(min_intron_bp=None),
    )
    keys = {
        candidate.evidence.evidence_id: candidate.structure_key
        for candidate in candidates if candidate.side == "five_prime"
    }
    assert keys["a"] == keys["b"]
    assert keys["c"] != keys["a"]


def _dummy_candidate(
    evidence: UTREvidence,
    *,
    endpoint: int = 80,
    side: str = "five_prime",
    structure: str = "-1-100",
    extension: int = 20,
) -> UTRCandidateSide:
    model = _model("dummy")
    transcript = _annotation(model).genes[0].transcripts[0]
    return UTRCandidateSide(
        gene_id="MoSelecta_000001",
        transcript_id=transcript.transcript_id,
        target=transcript,
        evidence=evidence,
        compatibility_class="full_containment",
        side=side,  # type: ignore[arg-type]
        endpoint0=endpoint,
        extension_interval=(endpoint, 100) if endpoint < 100 else (250, endpoint),
        added_exons=((endpoint, 100),) if endpoint < 100 else ((250, endpoint),),
        added_introns=(),
        complete_contributing_exon_lengths=(max(5, extension),),
        extension_exonic_bp=extension,
        extension_genomic_bp=extension,
        structure_key=structure,
        source_cds_class=".",
        source_cds_overlap_bp=0,
    )


def test_cap_three_support_order_and_prediction_source_cap() -> None:
    def support(counts: list[int]) -> float:
        pooled = PooledSupport()
        for index, count in enumerate(counts):
            pooled.add_candidate(_dummy_candidate(_tx(
                f"e{index}", ((80, 150), (200, 250)), counts=((f"s{index}", count),)
            )))
        return pooled.effective_support(3)

    assert support([42, 1, 1, 1]) > support([2, 1, 1, 1, 1]) > support([45])
    assert transcript_source_weight(1, 3) == 1
    assert transcript_source_weight(8, 3) == 4
    assert transcript_source_weight(64, 3) == 4

    pooled = PooledSupport()
    for index in range(10):
        pooled.add_candidate(_dummy_candidate(_pred(
            f"p{index}", ((80, 150), (200, 250)), source="one_source",
            model_id=f"m{index}",
        )))
    assert pooled.effective_support(3) == 1
    assert pooled.source_count == 1


def test_source_pooling_occurs_across_endpoints() -> None:
    structure = UTRStructure("k")
    structure.add(_dummy_candidate(_tx("a", ((80, 150),), counts=(("rna", 2),)), endpoint=80))
    structure.add(_dummy_candidate(_tx("b", ((70, 150),), counts=(("rna", 2),)), endpoint=70))
    assert structure.pooled().effective_support(3) == 3
    assert structure.pooled((80,)).effective_support(3) == 2


def test_bounded_endpoint_windows_do_not_gap_chain() -> None:
    assert _bounded_windows([0, 20, 40, 60], 40) == [
        [0, 20, 40], [20, 40, 60], [40, 60], [60]
    ]


def test_endpoint_majority_diffuse_fallback_and_observed_medoid() -> None:
    compact = UTRStructure("k")
    for endpoint, source in ((80, "a"), (85, "b"), (300, "c")):
        compact.add(_dummy_candidate(_tx(
            f"{source}", ((endpoint, 150),), counts=((source, 1),)
        ), endpoint=endpoint))
    decision = _select_endpoint(compact, _settings())
    assert decision.method == "local_consensus_medoid"
    assert decision.endpoint0 in {80, 85}

    half = UTRStructure("half")
    for endpoint, source in ((80, "a"), (85, "b"), (300, "c"), (305, "d")):
        half.add(_dummy_candidate(_tx(
            f"{source}", ((endpoint, 150),), counts=((source, 1),)
        ), endpoint=endpoint))
    diffuse = _select_endpoint(half, _settings())
    assert diffuse.method == "diffuse_global_medoid"
    assert diffuse.endpoint0 in {80, 85, 300, 305}
    assert diffuse.strong_multimodality
    assert weighted_observed_medoid({0: 1, 10: 1}) == 0


def test_prediction_endpoint_normalization_prevents_vote_multiplication() -> None:
    structure = UTRStructure("k")
    for index in range(10):
        structure.add(_dummy_candidate(_pred(
            f"many{index}", ((80, 150),), source="pred", model_id=f"m{index}"
        ), endpoint=80))
    structure.add(_dummy_candidate(_pred(
        "other", ((400, 450),), source="other", model_id="other"
    ), endpoint=400))
    decision = _select_endpoint(structure, _settings())
    assert decision.method == "diffuse_global_medoid"
    assert decision.endpoint0 == 80  # equal source votes, lower-coordinate medoid tie


def test_learned_limits_side_specific_pool_sparse_and_overrides() -> None:
    candidates = [
        _dummy_candidate(_tx(f"f{i}", ((80, 150),)), extension=i + 1)
        for i in range(1000)
    ]
    candidates += [
        _dummy_candidate(
            _tx(f"t{i}", ((200, 270),)), endpoint=270,
            side="three_prime", structure="250--1", extension=(i + 1) * 2,
        )
        for i in range(1000)
    ]
    limits = learn_extension_limits(candidates, UTRSettings())
    assert limits.five_prime_exonic_bp == 1000
    assert limits.three_prime_exonic_bp == 2000

    sparse = learn_extension_limits(candidates[:10], UTRSettings())
    assert sparse.five_prime_exonic_bp == 20_000
    assert sparse.five_prime_genomic_bp is None
    assert sparse.warnings

    overridden = learn_extension_limits([], UTRSettings(
        five_prime_exonic_override_bp=123,
        five_prime_genomic_override_bp=0,
    ))
    assert overridden.five_prime_exonic_bp == 123
    assert overridden.five_prime_genomic_bp is None

    pooled_candidates = candidates[:600] + candidates[1000:1600]
    pooled_limits = learn_extension_limits(pooled_candidates, UTRSettings())
    assert pooled_limits.five_prime_exonic_source == "learned_pooled_sides"
    assert pooled_limits.three_prime_exonic_source == "learned_pooled_sides"


def test_multi_source_p999_term_protects_a_supported_tail() -> None:
    all_values = [100] * 1995 + [10_000] * 5
    multi_source_values = [100] * 995 + [10_000] * 5
    values = {
        (side, geometry): list(all_values)
        for side in ("five_prime", "three_prime")
        for geometry in ("exonic", "genomic")
    }
    strong = {
        (side, geometry): list(multi_source_values)
        for side in ("five_prime", "three_prime")
        for geometry in ("exonic", "genomic")
    }
    limits = _limits_from_observations(values, strong, UTRSettings())
    assert limits.five_prime_exonic_bp == 10_000
    assert limits.five_prime_genomic_bp == 10_000


def test_multi_final_gene_bridge_rejected_but_same_gene_isoforms_can_share() -> None:
    first = _model("first")
    second = _model("second")
    evidence = _tx("bridge", ((80, 150), (200, 270)))
    two_genes = augment_final_annotation(
        _annotation(first, second), [evidence], [], _settings()
    )
    assert two_genes.summary.multi_gene_evidence_excluded == 1
    assert two_genes.summary.neither_assigned == 2

    same_gene = augment_final_annotation(
        _multi_isoform_annotation(first, second), [evidence], [], _settings()
    )
    assert same_gene.summary.multi_gene_evidence_excluded == 0
    assert same_gene.summary.both_assigned == 2


@pytest.mark.parametrize("rejected_side", ["five_prime", "three_prime"])
def test_invalid_transcript_side_cannot_create_multi_gene_bridge(
    monkeypatch: pytest.MonkeyPatch,
    rejected_side: str,
) -> None:
    """Bridge ownership is established only after side validation."""
    annotation = _annotation(_model("first"), _model("second"))
    evidence = _tx("shared", ((80, 150), (200, 270)))

    def synthetic_candidates(gene_id, transcript, index, settings, **kwargs):
        del transcript, index, settings, kwargs
        if gene_id == "MoSelecta_000001":
            return [SimpleNamespace(
                evidence=evidence,
                side="three_prime" if rejected_side == "five_prime" else "five_prime",
                validation_status="passed",
                rejection_reasons=(),
            )]
        return [SimpleNamespace(
            evidence=evidence,
            side=rejected_side,
            validation_status="rejected",
            rejection_reasons=("resolved_unrecognized_motif",),
        )]

    monkeypatch.setattr(
        utr_augmentation_module,
        "candidate_sides_for_transcript",
        synthetic_candidates,
    )
    bridged = utr_augmentation_module._multi_gene_bridge_ids(
        annotation, UTREvidenceIndex([]), _settings()
    )
    assert bridged == set()


def test_prediction_catalogue_uses_selected_and_unselected_normalized_exons() -> None:
    selected = _model("selected", source="A", exons=[(80, 150), (200, 270)])
    unselected = _model("unselected", source="B", exons=[(70, 150), (200, 280)])
    catalogue = build_prediction_utr_catalogue(
        [selected, unselected], selected_internal_ids={selected.output_transcript_id}
    )
    assert [item.prediction_selected for item in catalogue] == [True, False]
    assert catalogue[0].exon_parts == ((80, 150), (200, 270))
    result = augment_final_annotation(_annotation(selected), [], catalogue, _settings())
    row = result.audit_rows[0]
    assert "unselected_model" in {
        row["five_prime_prediction_selection_provenance"],
        row["three_prime_prediction_selection_provenance"],
    }


def test_prediction_cds_and_non_cds_extension_provenance_is_audited() -> None:
    target = _model("target")
    evidence = [
        _pred("cds", ((80, 150), (200, 250)), cds=((80, 150),)),
        _pred("noncds", ((100, 150), (200, 270)), source="other", cds=((100, 150), (200, 250))),
    ]
    result = _result(target, evidence)
    row = result.audit_rows[0]
    assert row["five_prime_prediction_cds_provenance"] == "within_source_cds"
    assert row["three_prime_prediction_cds_provenance"] == "outside_source_cds"


def test_structure_ranking_support_breadth_depth_and_lexical_ties() -> None:
    limits = learn_extension_limits([], _settings())
    candidates = [
        _dummy_candidate(_tx("a1", ((80, 150),), counts=(("A", 1),)), structure="a"),
        _dummy_candidate(_tx("a2", ((80, 150),), counts=(("B", 1),)), structure="a"),
        _dummy_candidate(_tx("b1", ((70, 150),), counts=(("A", 2),)), endpoint=70, structure="b"),
    ]
    five, _ = _select_sides_for_transcript(candidates, _settings(), limits)
    assert five is not None and five.candidate.structure_key == "a"

    lexical = [
        _dummy_candidate(_tx("x", ((80, 150),), counts=(("A", 1),)), structure="z"),
        _dummy_candidate(_tx("y", ((70, 150),), counts=(("B", 1),)), endpoint=70, structure="a"),
    ]
    five, _ = _select_sides_for_transcript(lexical, _settings(), limits)
    assert five is not None and five.candidate.structure_key == "a"


def test_exact_tie_cross_side_coherence_cannot_override_stronger_structure() -> None:
    limits = learn_extension_limits([], _settings())
    shared = _tx("shared", ((80, 150), (200, 270)), counts=(("A", 1),))
    candidates = [
        _dummy_candidate(shared, structure="five_shared"),
        _dummy_candidate(_tx("five_strong_a", ((70, 150),), counts=(("B", 1),)), endpoint=70, structure="five_strong"),
        _dummy_candidate(_tx("five_strong_b", ((70, 150),), counts=(("C", 1),)), endpoint=70, structure="five_strong"),
        _dummy_candidate(shared, endpoint=270, side="three_prime", structure="three_shared"),
    ]
    five, three = _select_sides_for_transcript(candidates, _settings(), limits)
    assert five is not None and five.candidate.structure_key == "five_strong"
    assert three is not None and three.candidate.structure_key == "three_shared"


def test_exact_numerical_cross_side_tie_prefers_shared_evidence() -> None:
    limits = learn_extension_limits([], _settings())
    shared = _tx("shared", ((80, 150), (200, 270)), counts=(("A", 1),))
    candidates = [
        _dummy_candidate(shared, structure="z_five_shared"),
        _dummy_candidate(_tx("f_lexical", ((70, 150),), counts=(("B", 1),)), endpoint=70, structure="a_five"),
        _dummy_candidate(shared, endpoint=270, side="three_prime", structure="z_three_shared"),
        _dummy_candidate(_tx("t_lexical", ((200, 280),), counts=(("C", 1),)), endpoint=280, side="three_prime", structure="a_three"),
    ]
    five, three = _select_sides_for_transcript(candidates, _settings(), limits)
    assert five is not None and five.candidate.evidence.evidence_id == "shared"
    assert three is not None and three.candidate.evidence.evidence_id == "shared"


def test_final_recheck_falls_through_to_next_ranked_structure(monkeypatch: pytest.MonkeyPatch) -> None:
    import moselecta.utr_augmentation as utr

    limits = learn_extension_limits([], _settings())
    strongest = _dummy_candidate(
        _tx("bad", ((80, 150),), counts=(("A", 4),)), structure="bad"
    )
    fallback = _dummy_candidate(
        _tx("good", ((70, 150),), counts=(("B", 1),)), endpoint=70, structure="good"
    )
    monkeypatch.setattr(
        utr,
        "_candidate_recheck",
        lambda candidate, settings, learned: candidate.evidence.evidence_id != "bad",
    )
    five, _ = _select_sides_for_transcript(
        [strongest, fallback], _settings(), limits
    )
    assert five is not None
    assert five.candidate.evidence.evidence_id == "good"
    assert five.quality_status == "fallback_structure_after_recheck"


def test_reconstruction_preserves_ids_cds_phases_and_gene_span(tmp_path: Path) -> None:
    model = _model("target", phases=["1", "2"])
    annotation = _annotation(model)
    result = augment_final_annotation(
        annotation, [_tx("both", ((80, 150), (200, 270)), counts=(("rna", 2),))],
        [], _settings(),
    )
    assert result.annotation.genes[0].start0 == 80
    assert result.annotation.genes[0].end0 == 270
    augmented = result.annotation.genes[0].transcripts[0]
    assert augmented.transcript_id == annotation.genes[0].transcripts[0].transcript_id
    assert augmented.exon_parts == ((80, 150), (200, 270))

    path = tmp_path / "utr.gff3"
    write_utr_augmented_gff3(result, path)
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows.append((fields, parse_attributes(fields[8])))
    assert [fields[2] for fields, _ in rows].count("mRNA") == 1
    assert [fields[2] for fields, _ in rows].count("five_prime_UTR") == 1
    assert [fields[2] for fields, _ in rows].count("three_prime_UTR") == 1
    cds = [(int(fields[3]) - 1, int(fields[4]), fields[7]) for fields, _ in rows if fields[2] == "CDS"]
    assert cds == [(100, 150, "1"), (200, 250, "2")]


def test_no_usable_evidence_emits_one_cds_derived_transcript() -> None:
    model = _model("target")
    result = _result(model, [])
    augmented = result.annotation.genes[0].transcripts[0]
    assert result.summary.neither_assigned == 1
    assert augmented.exon_parts == tuple(model.cds_parts)
    assert augmented.five_prime_utr_parts == ()
    assert augmented.three_prime_utr_parts == ()


def test_cross_final_gene_exact_cds_exon_is_hard_filtered_with_fallback() -> None:
    first = _model("first", parts=[(100, 150), (180, 200)])
    second = _model("second", parts=[(220, 260)])
    bad = _tx("bad", ((100, 150), (180, 200), (220, 260)), counts=(("strong", 8),))
    good = _tx("good", ((100, 150), (180, 200), (225, 270)), counts=(("weak", 1),))
    result = augment_final_annotation(_annotation(first, second), [bad, good], [], _settings())
    rows = {row["public_gene_id"]: row for row in result.audit_rows}
    assert result.summary.multi_gene_evidence_excluded == 0
    assert result.summary.authoritative_candidate_sides_rejected >= 1
    assert rows["MoSelecta_000001"]["three_prime_selected_evidence_ids"] == "good"
    assert rows["MoSelecta_000001"]["three_prime_quality_status"] == (
        "cross_final_gene_authoritative_filter_fallback"
    )
    assert "cross_final_gene_authoritative_exon_conflict" in rows[
        "MoSelecta_000001"
    ]["three_prime_cross_final_gene_conflict_reasons"]
    assert result.annotation.genes[0].transcripts[0].three_prime_utr_parts == ((225, 270),)
    assert result.annotation.genes[1].transcripts[0].exon_parts == ((220, 260),)
    assert tuple(first.cds_parts) == ((100, 150), (180, 200))
    assert tuple(second.cds_parts) == ((220, 260),)


def test_cross_final_gene_exact_cds_exon_without_fallback_leaves_side_unassigned() -> None:
    first = _model("first", parts=[(100, 150), (180, 200)])
    second = _model("second", parts=[(220, 260)])
    bad = _tx("bad", ((100, 150), (180, 200), (220, 260)))
    result = augment_final_annotation(_annotation(first, second), [bad], [], _settings())
    assert result.annotation.genes[0].transcripts[0].three_prime_selection is None
    assert result.audit_rows[0]["three_prime_quality_status"] == (
        "cross_final_gene_authoritative_filter_no_alternative"
    )


def test_cross_final_gene_ordinary_non_cds_source_exon_is_protected() -> None:
    first = _model("first", parts=[(100, 150), (180, 200)])
    second = _model(
        "second", parts=[(230, 260), (300, 340)],
        exons=[(220, 260), (300, 340)],
    )
    evidence = _tx("source_exon", ((100, 150), (180, 200), (220, 260)))
    result = augment_final_annotation(_annotation(first, second), [evidence], [], _settings())
    assert result.summary.multi_gene_evidence_excluded == 0
    assert result.annotation.genes[0].transcripts[0].three_prime_selection is None
    assert "cross_final_gene_ordinary_source_exon_conflict" in (
        result.audit_rows[0]["three_prime_cross_final_gene_conflict_reasons"]
    )


@pytest.mark.parametrize(
    ("other_parts", "reason"),
    [
        ([(205, 240), (280, 300)], "cross_final_gene_authoritative_intron_conflict"),
        ([(205, 240), (290, 310)], "cross_final_gene_authoritative_splice_conflict"),
    ],
)
def test_cross_final_gene_authoritative_intron_and_splice_are_hard_filtered(
    other_parts: list[tuple[int, int]], reason: str,
) -> None:
    first = _model("first", parts=[(100, 150), (180, 200)])
    second = _model("second", parts=other_parts)
    evidence = _tx("bad", ((100, 150), (180, 200), (220, 240), (280, 300)))
    result = augment_final_annotation(_annotation(first, second), [evidence], [], _settings())
    assert result.summary.multi_gene_evidence_excluded == 0
    assert result.annotation.genes[0].transcripts[0].three_prime_selection is None
    assert reason in result.audit_rows[0]["three_prime_cross_final_gene_conflict_reasons"]


def test_independent_selected_utrs_with_exact_exon_conflict_both_fall_back() -> None:
    first = _model("first", parts=[(100, 150), (180, 200)])
    second = _model("second", parts=[(300, 340), (370, 400)])
    evidence = [
        _tx("bad_a", ((100, 150), (180, 200), (220, 260)), counts=(("A", 8),)),
        _tx("bad_b", ((220, 260), (300, 340), (370, 400)), counts=(("B", 8),)),
        _tx("good_a", ((100, 150), (180, 200), (220, 245)), counts=(("C", 1),)),
        _tx("good_b", ((260, 280), (300, 340), (370, 400)), counts=(("D", 1),)),
    ]
    result = augment_final_annotation(_annotation(first, second), evidence, [], _settings())
    assert result.summary.multi_gene_evidence_excluded == 0
    assert result.summary.provisional_selected_sides_in_utr_conflicts == 2
    assert result.summary.conflict_fallback_sides == 2
    assert result.summary.conflict_no_alternative_sides == 0
    assert result.summary.conflict_resolution_rounds == 1
    assert result.audit_rows[0]["three_prime_selected_evidence_ids"] == "good_a"
    assert result.audit_rows[1]["five_prime_selected_evidence_ids"] == "good_b"
    assert all(
        row[side + "_quality_status"] == "cross_final_gene_conflict_fallback"
        for row, side in zip(result.audit_rows, ("three_prime", "five_prime"))
    )


@pytest.mark.parametrize(
    ("second_utr", "reason"),
    [
        (((225, 240), (280, 305)), "cross_final_gene_selected_utr_intron_conflict"),
        (((225, 240), (290, 305)), "cross_final_gene_selected_utr_splice_conflict"),
    ],
)
def test_independent_selected_utrs_with_intron_or_splice_conflict_lose_both_sides(
    second_utr: tuple[tuple[int, int], ...], reason: str,
) -> None:
    first = _model("first", parts=[(100, 150), (180, 200)])
    second = _model("second", parts=[(330, 350), (380, 410)])
    evidence = [
        _tx("a", ((100, 150), (180, 200), (220, 240), (280, 300))),
        _tx("b", second_utr + ((330, 350), (380, 410))),
    ]
    result = augment_final_annotation(_annotation(first, second), evidence, [], _settings())
    assert result.summary.multi_gene_evidence_excluded == 0
    assert result.summary.provisional_selected_sides_in_utr_conflicts == 2
    assert result.summary.conflict_no_alternative_sides == 2
    assert result.annotation.genes[0].transcripts[0].three_prime_selection is None
    assert result.annotation.genes[1].transcripts[0].five_prime_selection is None
    assert reason in result.audit_rows[0]["three_prime_cross_final_gene_conflict_reasons"]
    assert reason in result.audit_rows[1]["five_prime_cross_final_gene_conflict_reasons"]


def test_opposite_strand_coincident_utr_exons_are_allowed() -> None:
    first = _model("first", parts=[(100, 150), (180, 200)], strand="+")
    second = _model("second", parts=[(300, 340), (370, 400)], strand="-")
    evidence = [
        _tx("plus", ((100, 150), (180, 200), (220, 260)), strand="+"),
        _tx("minus", ((220, 260), (300, 340), (370, 400)), strand="-"),
    ]
    result = augment_final_annotation(_annotation(first, second), evidence, [], _settings())
    assert result.summary.conflict_resolution_rounds == 0
    assert result.annotation.genes[0].transcripts[0].three_prime_utr_parts == ((220, 260),)
    assert result.annotation.genes[1].transcripts[0].three_prime_utr_parts == ((220, 260),)


def test_same_strand_utr_cds_overlap_without_exact_link_is_allowed() -> None:
    first = _model("first", parts=[(100, 150), (180, 200)])
    second = _model("second", parts=[(230, 270)])
    evidence = _tx("overlap", ((100, 150), (180, 200), (220, 250)))
    result = augment_final_annotation(_annotation(first, second), [evidence], [], _settings())
    assert result.summary.authoritative_candidate_sides_rejected == 0
    assert result.summary.conflict_resolution_rounds == 0
    assert result.annotation.genes[0].transcripts[0].three_prime_utr_parts == ((220, 250),)


def test_iterative_cross_gene_fallback_is_order_independent() -> None:
    models = [
        _model("a", parts=[(100, 150), (180, 200)]),
        _model("b", parts=[(300, 340), (370, 400)]),
        _model("c", parts=[(500, 540), (570, 600)]),
    ]
    evidence = [
        _tx("a_bad", ((100, 150), (180, 200), (220, 260)), counts=(("A", 8),)),
        _tx("a_fallback_one", ((100, 150), (180, 200), (260, 280)), counts=(("B", 4),)),
        _tx("a_fallback_two", ((100, 150), (180, 200), (220, 245)), counts=(("C", 1),)),
        _tx("b_bad", ((220, 260), (300, 340), (370, 400)), counts=(("D", 8),)),
        _tx("c_bad", ((260, 280), (500, 540), (570, 600)), counts=(("E", 8),)),
    ]
    genes = tuple(
        FinalGene(
            public_gene_id(index),
            (FinalTranscript(public_transcript_id(public_gene_id(index), 1), model),),
        )
        for index, model in enumerate(models, start=1)
    )

    def run(ordered: tuple[FinalGene, ...]):
        return augment_final_annotation(FinalAnnotation(ordered), evidence, [], _settings())

    forward, reverse = run(genes), run(tuple(reversed(genes)))
    assert forward.summary.conflict_resolution_rounds == 2
    assert reverse.summary.conflict_resolution_rounds == 2
    assert forward.summary.conflict_fallback_sides == 1
    assert forward.summary.conflict_no_alternative_sides == 2
    decisions = lambda result: {
        row["public_gene_id"]: (
            row["five_prime_selected_evidence_ids"],
            row["three_prime_selected_evidence_ids"],
            row["five_prime_quality_status"],
            row["three_prime_quality_status"],
        )
        for row in result.audit_rows
    }
    assert decisions(forward) == decisions(reverse)
    assert decisions(forward)["MoSelecta_000001"][1] == "a_fallback_two"
    assert decisions(forward)["MoSelecta_000002"][0] == "."
    assert decisions(forward)["MoSelecta_000003"][0] == "."


def test_incremental_selected_index_matches_full_scan_after_simultaneous_updates() -> None:
    models = (
        _model("a", parts=[(100, 150), (180, 200)]),
        _model("b", parts=[(300, 340), (370, 400)]),
        _model("c", parts=[(500, 540), (570, 600)]),
        _model("d", parts=[(700, 740), (770, 800)]),
    )
    annotation = _annotation(*models)
    evidence = [
        _tx("a_bad", ((100, 150), (180, 200), (220, 260))),
        _tx("a_middle", ((100, 150), (180, 200), (260, 280))),
        _tx("a_later", ((100, 150), (180, 200), (280, 295))),
        _tx("b_bad", ((220, 260), (300, 340), (370, 400))),
        _tx("c_bad", ((260, 280), (500, 540), (570, 600))),
        _tx("d_bad", ((280, 295), (700, 740), (770, 800))),
    ]
    lookup = {
        transcript.transcript_id: (gene.gene_id, transcript)
        for gene in annotation.genes for transcript in gene.transcripts
    }
    evidence_index = UTREvidenceIndex(evidence, annotation)

    def selection(transcript_id: str, evidence_id: str) -> UTRSideSelection:
        gene_id, transcript = lookup[transcript_id]
        candidate = next(candidate for candidate in candidate_sides_for_transcript(
            gene_id, transcript, evidence_index, _settings()
        ) if candidate.evidence.evidence_id == evidence_id)
        return UTRSideSelection(
            side=candidate.side, candidate=candidate, endpoint_method="test",
            effective_support=1.0, independent_source_count=1,
            multi_member_transcript_source_count=0, eligible_structure_count=1,
            strong_competing_structure=False, strong_multimodality=False,
            local_support_fraction=1.0, adaptive_window_bp=20,
            selected_evidence_ids=(evidence_id,), selected_sources=("rna",),
            evidence_class="transcript", prediction_cds_classes=(),
            prediction_selection_classes=(),
        )

    ids = list(lookup)
    selected = {
        ids[0]: (None, selection(ids[0], "a_bad")),
        ids[1]: (selection(ids[1], "b_bad"), None),
        ids[2]: (selection(ids[2], "c_bad"), None),
        ids[3]: (selection(ids[3], "d_bad"), None),
    }
    frozen = FrozenFinalGeneStructureIndex(annotation)
    incremental = _SelectedPopulationStructureIndex(lookup, selected, frozen)

    def assert_matches_reference() -> None:
        assert incremental.current_conflicts() == _selected_global_conflicts(
            lookup, selected, frozen
        )

    assert_matches_reference()
    for changes in (
        {ids[0]: (None, selection(ids[0], "a_middle")), ids[1]: (None, None)},
        {ids[0]: (None, selection(ids[0], "a_later")), ids[2]: (None, None)},
        {ids[0]: (None, None), ids[3]: (None, None)},
    ):
        selected.update(changes)
        incremental.replace(changes)
        assert_matches_reference()


def test_three_round_fallback_chain_is_order_independent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import moselecta.utr_augmentation as utr_module

    models = (
        _model("a", parts=[(100, 150), (180, 200)]),
        _model("b", parts=[(300, 340), (370, 400)]),
        _model("c", parts=[(500, 540), (570, 600)]),
        _model("d", parts=[(700, 740), (770, 800)]),
    )
    evidence = [
        _tx("a_bad", ((100, 150), (180, 200), (220, 260)), counts=(("A", 8),)),
        _tx("a_second", ((100, 150), (180, 200), (260, 280)), counts=(("B", 4),)),
        _tx("a_third", ((100, 150), (180, 200), (280, 295)), counts=(("C", 2),)),
        _tx("a_safe", ((100, 150), (180, 200), (220, 245)), counts=(("D", 1),)),
        _tx("b_bad", ((220, 260), (300, 340), (370, 400)), counts=(("E", 8),)),
        _tx("c_bad", ((260, 280), (500, 540), (570, 600)), counts=(("F", 8),)),
        _tx("d_bad", ((280, 295), (700, 740), (770, 800)), counts=(("G", 8),)),
    ]
    genes = _annotation(*models).genes
    generation_counts: dict[str, int] = {}
    original_generate = utr_module.candidate_sides_for_transcript

    def counted_generate(*args, **kwargs):
        transcript_id = args[1].transcript_id
        generation_counts[transcript_id] = generation_counts.get(transcript_id, 0) + 1
        return original_generate(*args, **kwargs)

    monkeypatch.setattr(utr_module, "candidate_sides_for_transcript", counted_generate)
    forward = augment_final_annotation(FinalAnnotation(genes), evidence, [], _settings())
    assert generation_counts["MoSelecta_000001.1"] == 4
    generation_counts.clear()
    reverse = augment_final_annotation(FinalAnnotation(tuple(reversed(genes))), evidence, [], _settings())
    assert generation_counts["MoSelecta_000001.1"] == 4
    assert forward.summary.conflict_resolution_rounds == 3
    assert reverse.summary.conflict_resolution_rounds == 3
    assert forward.summary.conflict_fallback_sides == 1
    assert forward.summary.conflict_no_alternative_sides == 3
    decisions = lambda result: {
        row["public_gene_id"]: (
            row["five_prime_selected_evidence_ids"],
            row["three_prime_selected_evidence_ids"],
            row["five_prime_quality_status"],
            row["three_prime_quality_status"],
        )
        for row in result.audit_rows
    }
    assert decisions(forward) == decisions(reverse)
    assert decisions(forward)["MoSelecta_000001"][1] == "a_safe"


def test_structure_source_count_counts_distinct_sources_once() -> None:
    model = _model("breadth")
    annotation = _annotation(model)
    evidence = [
        _tx("rna_one", ((100, 150), (200, 250), (270, 300)), counts=(("rna", 8),)),
        _tx("rna_two", ((100, 150), (200, 250), (275, 300)), counts=(("rna", 1),)),
        _pred("prediction", ((100, 150), (200, 250), (270, 300)), source="pred"),
    ]
    gene = annotation.genes[0]
    candidates = candidate_sides_for_transcript(
        gene.gene_id, gene.transcripts[0], UTREvidenceIndex(evidence, annotation), _settings()
    )
    assert _structure_source_count(candidates) == 2


def test_exact_structural_keys_separate_type_chromosome_and_strand() -> None:
    exons = ((100, 150), (200, 250))
    keys = _structural_features("chr1", "+", exons)
    assert keys.isdisjoint(_structural_features("chr1", "-", exons))
    assert keys.isdisjoint(_structural_features("chr2", "+", exons))
    assert ("chr1", "+", "exon", (100, 150)) in keys
    assert ("chr1", "+", "intron", (150, 200)) in keys
    assert ("chr1", "+", "splice", 150) in keys
    assert ("chr1", "+", "exon", (150, 200)) not in keys


def test_cli_flag_adds_only_separate_products_and_ignores_cds_only_writer_mode(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    input_a = tmp_path / "a.gff3"
    input_b = tmp_path / "b.gff3"
    template = (
        "##gff-version 3\n"
        "chr1\t{source}\tgene\t81\t270\t.\t+\t.\tID=g{source}\n"
        "chr1\t{source}\tmRNA\t81\t270\t.\t+\t.\tID=t{source};Parent=g{source}\n"
        "chr1\t{source}\texon\t81\t150\t.\t+\t.\tParent=t{source}\n"
        "chr1\t{source}\texon\t201\t270\t.\t+\t.\tParent=t{source}\n"
        "chr1\t{source}\tCDS\t101\t150\t.\t+\t1\tParent=t{source}\n"
        "chr1\t{source}\tCDS\t201\t250\t.\t+\t2\tParent=t{source}\n"
    )
    input_a.write_text(template.format(source="A"))
    input_b.write_text(template.format(source="B"))
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")

    def run(prefix: Path, *, add_utrs: bool) -> dict[str, Path]:
        command = [
            sys.executable, "-m", "moselecta",
            "--inputs", str(input_a), str(input_b),
            "--source-name", "A", "B",
            "--out-prefix", str(prefix),
            "--pruning-level", "none",
            "--output-profile", "minimal",
            "--progress-log-level", "none",
            "--final-gff3-cds-only-structure",
        ]
        if add_utrs:
            command.append("--add-utrs")
        subprocess.run(command, cwd=root, env=env, check=True, capture_output=True, text=True)
        return {
            path.name.removeprefix(prefix.name + "."): path
            for path in tmp_path.glob(prefix.name + ".*")
        }

    disabled = run(tmp_path / "disabled", add_utrs=False)
    enabled = run(tmp_path / "enabled", add_utrs=True)
    new_suffixes = {
        "final_annotation.with_utrs.gff3",
        "utr_augmentation.summary.tsv",
    }
    assert set(enabled) - set(disabled) == new_suffixes
    for suffix in set(disabled) - {
        "run_configuration.json", "utr_junction_validation.summary.tsv",
    }:
        assert disabled[suffix].read_bytes() == enabled[suffix].read_bytes(), suffix
    assert (
        enabled["utr_junction_validation.summary.tsv"].read_text()
        != disabled["utr_junction_validation.summary.tsv"].read_text()
    )

    ordinary = enabled["final_annotation.gff3"].read_text()
    augmented = enabled["final_annotation.with_utrs.gff3"].read_text()
    assert "five_prime_UTR" not in ordinary
    assert "three_prime_UTR" not in ordinary
    assert "five_prime_UTR" in augmented
    assert "three_prime_UTR" in augmented

    def identities_and_cds(path: Path):
        identities: list[tuple[str, str]] = []
        cds: list[tuple[str, str, str, str]] = []
        for line in path.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            attrs = parse_attributes(fields[8])
            if fields[2] in {"gene", "mRNA"}:
                identities.append((fields[2], attrs["ID"]))
            if fields[2] == "CDS":
                cds.append((attrs["Parent"], fields[3], fields[4], fields[7]))
        return identities, cds

    assert identities_and_cds(enabled["final_annotation.gff3"]) == identities_and_cds(
        enabled["final_annotation.with_utrs.gff3"]
    )


def test_augmented_mrnas_reuse_authoritative_cds_metadata_per_transcript(
    tmp_path: Path,
) -> None:
    primary = _model("primary", phases=["1", "2"])
    alternative = _model(
        "alternative",
        parts=[(100, 150), (220, 250)],
        phases=["1", "0"],
    )
    gene_id = "MoSelecta_000001"
    alt_metadata = AlternativeTranscriptMetadata(
        candidate_order_rank=1,
        preset="balanced",
        support_policy="test",
        evidence_requirement="none",
        support_routes=("structural",),
        mosaic_policy_required=False,
        reciprocal_cds_overlap=0.8,
        frame_concordance=1.0,
        primary_score_ratio=0.9,
    )
    annotation = FinalAnnotation((FinalGene(gene_id, (
        FinalTranscript(
            public_transcript_id(gene_id, 1),
            primary,
            annotation_summary=AnnotationSummary(
                model_selection_score_scaled=0.71,
                model_structural_support_score_scaled=0.66,
                protein_evidence_score_scaled=0.95,
                transcript_evidence_score_scaled=0.49,
                sequence_integrity_status="evaluated_no_concern",
            ),
        ),
        FinalTranscript(
            public_transcript_id(gene_id, 2),
            alternative,
            transcript_role="alternative",
            annotation_summary=AnnotationSummary(
                model_selection_score_scaled=0.42,
                model_structural_support_score_scaled=0.33,
            ),
            alternative_metadata=alt_metadata,
        ),
    )),))
    evidence = [_tx("extension", ((80, 150), (200, 250), (270, 290)))]
    result = augment_final_annotation(annotation, evidence, [], _settings())
    ordinary_path = tmp_path / "ordinary.gff3"
    augmented_path = tmp_path / "augmented.gff3"
    write_final_annotation_gff3(
        annotation,
        ordinary_path,
        FinalAnnotationOutputMetadata(False, False, False),
    )
    write_utr_augmented_gff3(result, augmented_path)

    def rows(path: Path, feature_type: str):
        return [
            (fields, attrs)
            for fields, attrs in (
                (line.split("\t"), parse_attributes(line.split("\t")[8]))
                for line in path.read_text().splitlines()
                if line and not line.startswith("#")
            )
            if fields[2] == feature_type
        ]

    ordinary_mrnas = {attrs["ID"]: (fields, attrs) for fields, attrs in rows(ordinary_path, "mRNA")}
    augmented_mrnas = {attrs["ID"]: (fields, attrs) for fields, attrs in rows(augmented_path, "mRNA")}
    compact_keys = {
        "model_selection_score_scaled",
        "model_structural_support_score_scaled",
        "protein_evidence_score_scaled",
        "transcript_evidence_score_scaled",
        "sequence_integrity_status",
    }
    for transcript_id in ordinary_mrnas:
        ordinary_attrs = ordinary_mrnas[transcript_id][1]
        augmented_attrs = augmented_mrnas[transcript_id][1]
        assert {
            key: ordinary_attrs[key] for key in compact_keys if key in ordinary_attrs
        } == {
            key: augmented_attrs[key] for key in compact_keys if key in augmented_attrs
        }
    assert augmented_mrnas[public_transcript_id(gene_id, 1)][1]["protein_evidence_score_scaled"] == "0.950000"
    assert "protein_evidence_score_scaled" not in augmented_mrnas[public_transcript_id(gene_id, 2)][1]
    assert "transcript_evidence_score_scaled" not in augmented_mrnas[public_transcript_id(gene_id, 2)][1]
    assert "sequence_integrity_status" not in augmented_mrnas[public_transcript_id(gene_id, 2)][1]
    assert "moselecta_utr_augmented" not in ordinary_mrnas[public_transcript_id(gene_id, 1)][1]
    assert augmented_mrnas[public_transcript_id(gene_id, 1)][1]["moselecta_utr_augmented"] == "true"
    assert augmented_mrnas[public_transcript_id(gene_id, 2)][1]["moselecta_alt_preset"] == "balanced"

    ordinary_cds = [
        (attrs["Parent"], fields[3], fields[4], fields[7])
        for fields, attrs in rows(ordinary_path, "CDS")
    ]
    augmented_cds = [
        (attrs["Parent"], fields[3], fields[4], fields[7])
        for fields, attrs in rows(augmented_path, "CDS")
    ]
    assert augmented_cds == ordinary_cds
    for feature_type in ("exon", "five_prime_UTR", "three_prime_UTR", "CDS"):
        for parent in augmented_mrnas:
            ids = [
                attrs["ID"] for _fields, attrs in rows(augmented_path, feature_type)
                if attrs["Parent"] == parent
            ]
            assert ids == [f"{parent}.{feature_type}{i}" for i in range(1, len(ids) + 1)]
    ordinary_primary_fields = ordinary_mrnas[public_transcript_id(gene_id, 1)][0]
    augmented_primary_fields = augmented_mrnas[public_transcript_id(gene_id, 1)][0]
    assert ordinary_primary_fields[3:5] == ["101", "250"]
    assert augmented_primary_fields[3:5] == ["81", "290"]
