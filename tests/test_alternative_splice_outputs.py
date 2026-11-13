from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.alternative_splice_outputs import (
    build_alternative_splice_audit_candidates,
    write_alternative_splice_candidate_audit_gff3,
)
from moselecta.alternative_splice_selection import (
    AlternativeSpliceSelectionConfig,
    select_alternative_splice_candidates,
)
from moselecta.annotation_outputs import (
    FinalAnnotationOutputMetadata,
    write_final_annotation_gff3,
    write_standardised_primary_gff3,
)
from moselecta.final_annotation import (
    AlternativeTranscriptMetadata,
    FinalAnnotation,
    FinalGene,
    FinalTranscript,
    augment_annotation,
    build_primary_annotation,
    public_transcript_id,
)
from moselecta.gffio import parse_attributes, parse_input
from moselecta.models import CodingModel, FeatureRecord


def _write_final(
    annotation: FinalAnnotation,
    config: AlternativeSpliceSelectionConfig,
    path: Path,
    *,
    cds_only_structure: bool = False,
) -> None:
    write_final_annotation_gff3(
        annotation,
        path,
        FinalAnnotationOutputMetadata(True, True, True, config),
        cds_only_structure=cds_only_structure,
    )


def _record(
    feature_type: str,
    start0: int,
    end0: int,
    *,
    source: str,
    strand: str,
    phase: str = ".",
    parent: str = "input_tx",
    record_index: int = 0,
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
        attrs={"Parent": parent, "ID": f"{parent}.{feature_type}{record_index}"},
        raw_attrs=".",
        input_file=Path(f"{source}.gff3"),
        source_name=source,
        line_number=record_index + 1,
        record_index=record_index,
        raw_line="",
    )


def _model(
    name: str,
    parts: list[tuple[int, int]],
    *,
    source: str,
    strand: str = "+",
    model_start0: int | None = None,
    model_end0: int | None = None,
) -> CodingModel:
    transcript_id = f"input_{name}"
    transcript_start = parts[0][0] if model_start0 is None else model_start0
    transcript_end = parts[-1][1] if model_end0 is None else model_end0
    transcript = _record(
        "mRNA",
        transcript_start,
        transcript_end,
        source=source,
        strand=strand,
        parent=f"gene_{name}",
    )
    children: list[FeatureRecord] = []
    if transcript_start < parts[0][0]:
        children.append(_record(
            "five_prime_UTR",
            transcript_start,
            parts[0][0],
            source=source,
            strand=strand,
            parent=transcript_id,
            record_index=1,
        ))
    for index, (start0, end0) in enumerate(parts, start=1):
        children.extend((
            _record(
                "exon", start0, end0, source=source, strand=strand,
                parent=transcript_id, record_index=index * 2,
            ),
            _record(
                "CDS", start0, end0, source=source, strand=strand,
                phase=str((index - 1) % 3), parent=transcript_id,
                record_index=index * 2 + 1,
            ),
        ))
    if transcript_end > parts[-1][1]:
        children.append(_record(
            "three_prime_UTR",
            parts[-1][1],
            transcript_end,
            source=source,
            strand=strand,
            parent=transcript_id,
            record_index=100,
        ))
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"input_gene_{name}",
        input_transcript_id=transcript_id,
        output_gene_id=f"internal_gene_{name}",
        output_transcript_id=f"internal_transcript_{name}",
        chrom="chr1",
        strand=strand,
        cds_parts=parts,
        gene_record=None,
        transcript_record=transcript,
        child_records=children,
        superlocus_id="super1",
        sublocus_id="sub1",
    )


def _score_row(score: float = 0.9) -> dict[str, object]:
    return {
        "model_selection_score": score,
        "model_selection_score_scaled": score,
        "model_structural_support_score_scaled": 0.8,
        "protein_evidence_status": "protein_evidence_not_configured",
        "transcript_evidence_status": "transcript_evidence_not_configured",
        "sequence_integrity_status": "genome_not_configured",
    }


def _candidate_row(
    candidate: CodingModel,
    primary: CodingModel | None,
    *,
    rank: int | str = 1,
    ownership: str = "unique_owner",
    score: float = 0.9,
    exact_support: bool = True,
    full_mosaic: bool = False,
    exact_equivalent: bool = False,
    **overrides: object,
) -> dict[str, object]:
    primary_id = primary.output_transcript_id if primary is not None else "."
    row: dict[str, object] = {
        "alternative_candidate_assessment_id": (
            f"assessment:{candidate.output_transcript_id}:{primary_id}"
        ),
        "candidate_output_transcript_id": candidate.output_transcript_id,
        "primary_output_transcript_id": primary_id,
        "ownership_status": ownership,
        "structural_candidate": primary is not None,
        "alternative_candidate_order_rank": rank,
        "exact_cds_equivalent_to_primary": exact_equivalent,
        "same_intron_chain_terminal_variant": False,
        "frame_concordance_status": "calculated",
        "sequence_integrity_default_exclusion": False,
        "sequence_integrity_status": "unavailable",
        "intrinsic_default_exclusion": False,
        "exact_multisource_model_support": exact_support,
        "protein_exact_cds_chain_support": False,
        "exact_coherent_transcript_support": False,
        "full_structural_mosaic": full_mosaic,
        "mosaic_cds_introns_status": "complete",
        "mosaic_cds_splice_sites_status": "complete",
        "mosaic_internal_cds_exons_status": "complete",
        "mosaic_terminal_cds_exons_status": "complete",
        "mosaic_cds_base_consensus_status": "complete",
        "same_frame_overlap_fraction": 0.99,
        "candidate_model_selection_score": score,
        "candidate_model_selection_score_scaled": score,
        "primary_model_selection_score": 1.0,
    }
    row.update(overrides)
    return row


def _selection(
    primary: CodingModel,
    candidates: list[CodingModel],
    rows: list[dict[str, object]],
    config: AlternativeSpliceSelectionConfig,
):
    score_rows = {primary.output_transcript_id: _score_row(1.0)}
    score_rows.update({
        candidate.output_transcript_id: _score_row(
            float(row.get("candidate_model_selection_score_scaled", 0.9))
        )
        for candidate, row in zip(candidates, rows)
    })
    annotation = build_primary_annotation([primary], score_rows)
    result = select_alternative_splice_candidates(
        candidate_rows=rows,
        all_scored_models=[primary, *candidates],
        primary_annotation=annotation,
        config=config,
        terminal_cds_3prime_wobble_bp=3,
    )
    return result, annotation, score_rows


def _feature_rows(path: Path) -> list[tuple[list[str], dict[str, str]]]:
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows.append((fields, parse_attributes(fields[8])))
    return rows


@pytest.mark.parametrize("strand", ("+", "-"))
@pytest.mark.parametrize(
    ("alt_start", "alt_end"),
    ((70, 500), (100, 530), (70, 530)),
)
def test_final_annotation_preserves_primary_and_writes_union_span_and_alt_metadata(
    tmp_path: Path,
    strand: str,
    alt_start: int,
    alt_end: int,
) -> None:
    primary = _model(
        "primary", [(100, 180), (260, 340), (420, 500)], source="Annevo",
        strand=strand, model_start0=100, model_end0=500,
    )
    alternative = _model(
        "alternative", [(100, 180), (260, 330), (410, 500)], source="Braker",
        strand=strand, model_start0=alt_start, model_end0=alt_end,
    )
    row = _candidate_row(alternative, primary)
    config = AlternativeSpliceSelectionConfig.resolve()
    result, annotation, score_rows = _selection(
        primary, [alternative], [row], config
    )
    augmented = augment_annotation(annotation, result, score_rows)
    path = tmp_path / f"final-{strand}-{alt_start}-{alt_end}.gff3"
    _write_final(augmented, config, path)
    rows = _feature_rows(path)
    gene_fields, gene_attrs = next(item for item in rows if item[0][2] == "gene")
    transcripts = [item for item in rows if item[0][2] == "mRNA"]

    assert gene_fields[3:5] == [str(min(100, alt_start) + 1), str(max(500, alt_end))]
    assert gene_attrs["moselecta_transcript_count"] == "2"
    assert [attrs["ID"] for _fields, attrs in transcripts] == [
        "MoSelecta_000001.1", "MoSelecta_000001.2"
    ]
    primary_fields, primary_attrs = transcripts[0]
    alt_fields, alt_attrs = transcripts[1]
    assert primary_fields[1] == "MoSelecta:Annevo"
    assert primary_fields[3:5] == ["101", "500"]
    assert primary_attrs["moselecta_transcript_role"] == "primary"
    assert not any(key.startswith("moselecta_alt_") for key in primary_attrs)
    assert alt_fields[1] == "MoSelecta:Braker"
    assert alt_fields[3:5] == [str(alt_start + 1), str(alt_end)]
    assert alt_attrs["moselecta_transcript_role"] == "alternative"
    assert alt_attrs["original_source"] == "Braker"
    assert alt_attrs["moselecta_alt_candidate_order_rank"] == "1"
    assert alt_attrs["moselecta_alt_preset"] == "intermediate"
    assert alt_attrs["moselecta_alt_support_policy"] == "exact"
    assert alt_attrs["moselecta_alt_evidence_requirement"] == "none"
    assert alt_attrs["moselecta_alt_support_routes"] == "exact_multisource_model_support"
    assert alt_attrs["moselecta_alt_frame_concordance"] == "0.990000"
    assert alt_attrs["model_selection_score_scaled"] == "0.900000"
    alt_children = [
        (fields, attrs) for fields, attrs in rows
        if attrs.get("Parent") == "MoSelecta_000001.2"
    ]
    child_types = {fields[2] for fields, _attrs in alt_children}
    assert child_types >= {"exon", "CDS"}
    if alt_start < 100:
        assert "five_prime_UTR" in child_types
    if alt_end > 500:
        assert "three_prime_UTR" in child_types
    assert [fields[7] for fields, _attrs in alt_children if fields[2] == "CDS"] == (
        ["0", "1", "0"] if strand == "+" else ["1", "2", "2"]
    )
    text = path.read_text()
    assert "# moselecta_alternative_splice_enabled=true\n" in text
    assert "# moselecta_alternative_splice_selector_evaluated=true\n" in text
    assert "# moselecta_alternative_splice_alternatives_emitted=true\n" in text
    assert "# moselecta_alternative_splice_evidence_requirement=none\n" in text


def test_primary_only_final_metadata_and_feature_content(tmp_path: Path) -> None:
    primary = _model("primary", [(0, 90), (180, 270)], source="P")
    annotation = build_primary_annotation(
        [primary], {primary.output_transcript_id: _score_row()}
    )
    standardised = tmp_path / "standardised.gff3"
    final = tmp_path / "final.gff3"
    write_standardised_primary_gff3(annotation, standardised)
    write_final_annotation_gff3(
        annotation,
        final,
        FinalAnnotationOutputMetadata(False, False, False),
    )
    final_text = final.read_text()
    assert "# moselecta_alternative_splice_enabled=false\n" in final_text
    assert "# moselecta_alternative_splice_selector_evaluated=false\n" in final_text
    assert "# moselecta_alternative_splice_alternatives_emitted=false\n" in final_text
    assert "moselecta_alternative_splice_preset_requested" not in final_text
    standardised_features = [
        line for line in standardised.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    final_features = [
        line for line in final_text.splitlines()
        if line and not line.startswith("#")
    ]
    assert final_features == standardised_features


def test_public_primary_outputs_honour_cds_only_structure(tmp_path: Path) -> None:
    primary = _model(
        "primary",
        [(100, 180), (260, 340), (420, 500)],
        source="P",
        model_start0=70,
        model_end0=540,
    )
    annotation = build_primary_annotation(
        [primary], {primary.output_transcript_id: _score_row()}
    )
    standardised = tmp_path / "standardised.cds-only.gff3"
    final = tmp_path / "final.cds-only.gff3"
    write_standardised_primary_gff3(
        annotation, standardised, cds_only_structure=True
    )
    write_final_annotation_gff3(
        annotation,
        final,
        FinalAnnotationOutputMetadata(False, False, False),
        cds_only_structure=True,
    )

    for path in (standardised, final):
        rows = _feature_rows(path)
        gene_fields, gene_attrs = next(item for item in rows if item[0][2] == "gene")
        mrna_fields, mrna_attrs = next(item for item in rows if item[0][2] == "mRNA")
        children = [item for item in rows if item[1].get("Parent") == mrna_attrs["ID"]]

        assert gene_fields[3:5] == ["101", "500"]
        assert mrna_fields[3:5] == ["101", "500"]
        assert gene_attrs["ID"] == "MoSelecta_000001"
        assert mrna_attrs["ID"] == "MoSelecta_000001.1"
        assert gene_attrs["final_gff3_cds_only_structure"] == "true"
        assert mrna_attrs["final_gff3_cds_only_structure"] == "true"
        assert [fields[2] for fields, _attrs in children] == [
            "exon", "CDS", "exon", "CDS", "exon", "CDS"
        ]
        cds = [fields for fields, _attrs in children if fields[2] == "CDS"]
        assert [(fields[3], fields[4]) for fields in cds] == [
            ("101", "180"), ("261", "340"), ("421", "500")
        ]
        assert [fields[7] for fields in cds] == ["0", "1", "2"]


@pytest.mark.parametrize("strand", ("+", "-"))
def test_public_final_with_alternatives_honours_cds_only_structure(
    tmp_path: Path,
    strand: str,
) -> None:
    primary = _model(
        "primary", [(100, 180), (260, 340), (420, 500)], source="P",
        strand=strand, model_start0=60, model_end0=540,
    )
    alternative = _model(
        "alternative", [(80, 180), (260, 330), (410, 520)], source="A",
        strand=strand, model_start0=40, model_end0=570,
    )
    config = AlternativeSpliceSelectionConfig.resolve()
    augmented = FinalAnnotation((FinalGene(
        gene_id="MoSelecta_000001",
        transcripts=(
            FinalTranscript(
                transcript_id="MoSelecta_000001.1",
                source_model=primary,
            ),
            FinalTranscript(
                transcript_id="MoSelecta_000001.2",
                source_model=alternative,
                transcript_role="alternative",
                alternative_metadata=AlternativeTranscriptMetadata(
                    candidate_order_rank=1,
                    preset="intermediate",
                    support_policy="exact",
                    evidence_requirement="none",
                    support_routes=("exact_multisource_model_support",),
                    mosaic_policy_required=False,
                    reciprocal_cds_overlap=0.9,
                    frame_concordance=0.9,
                    primary_score_ratio=0.9,
                ),
            ),
        ),
    ),)
    )
    default_path = tmp_path / f"default-{strand}.gff3"
    cds_only_path = tmp_path / f"cds-only-{strand}.gff3"
    _write_final(augmented, config, default_path)
    _write_final(
        augmented, config, cds_only_path, cds_only_structure=True
    )

    default_rows = _feature_rows(default_path)
    rows = _feature_rows(cds_only_path)
    default_transcript_ids = [
        attrs["ID"] for fields, attrs in default_rows if fields[2] == "mRNA"
    ]
    transcript_rows = [item for item in rows if item[0][2] == "mRNA"]
    assert [attrs["ID"] for _fields, attrs in transcript_rows] == default_transcript_ids
    assert default_transcript_ids == [
        "MoSelecta_000001.1", "MoSelecta_000001.2"
    ]

    gene_fields, gene_attrs = next(item for item in rows if item[0][2] == "gene")
    assert gene_fields[3:5] == ["81", "520"]
    assert gene_attrs["moselecta_transcript_count"] == "2"
    assert [fields[3:5] for fields, _attrs in transcript_rows] == [
        ["101", "500"], ["81", "520"]
    ]
    for fields, attrs in transcript_rows:
        children = [
            item for item in rows if item[1].get("Parent") == attrs["ID"]
        ]
        assert {child_fields[2] for child_fields, _child_attrs in children} == {
            "exon", "CDS"
        }
        assert all(
            child_attrs.get("generated_from") == "CDS"
            for _child_fields, child_attrs in children
        )
        phases = [
            child_fields[7]
            for child_fields, _child_attrs in children
            if child_fields[2] == "CDS"
        ]
        if attrs["ID"].endswith(".1"):
            assert phases == (["0", "1", "2"] if strand == "+" else ["1", "0", "2"])
        else:
            assert phases == (["0", "2", "1"] if strand == "+" else ["2", "0", "2"])


def test_final_gff_records_model_source_consensus_configuration(
    tmp_path: Path,
) -> None:
    primary = _model("primary", [(0, 90), (180, 270)], source="P")
    annotation = build_primary_annotation(
        [primary], {primary.output_transcript_id: _score_row()}
    )
    config = AlternativeSpliceSelectionConfig.resolve(
        min_model_support_sources=2,
        min_model_support_source_fraction=0.50,
    ).for_prediction_source_count(5)
    path = tmp_path / "source-consensus.gff3"
    write_final_annotation_gff3(
        annotation,
        path,
        FinalAnnotationOutputMetadata(True, True, False, config),
    )
    text = path.read_text()
    assert "# moselecta_alternative_splice_prediction_source_count=5\n" in text
    assert "# moselecta_alternative_splice_min_model_support_sources=2\n" in text
    assert (
        "# moselecta_alternative_splice_min_model_support_source_fraction=0.5\n"
        in text
    )
    assert (
        "# moselecta_alternative_splice_effective_min_model_support_sources=3\n"
        in text
    )
def test_audit_population_hierarchy_and_namespace_are_truthful(tmp_path: Path) -> None:
    primary = _model("primary", [(0, 90), (180, 270), (360, 450)], source="P")
    accepted = _model("accepted", [(0, 90), (180, 240), (330, 450)], source="A")
    exact = _model("exact", list(primary.cds_parts), source="B")
    threshold = _model("threshold", [(0, 90), (180, 246), (336, 450)], source="C")
    zero = _model("zero", [(700, 790), (880, 970)], source="D")
    ambiguous = _model("ambiguous", [(0, 90), (180, 252), (342, 450)], source="E")
    rows = [
        _candidate_row(accepted, primary, rank=1),
        _candidate_row(exact, primary, rank=2, exact_equivalent=True),
        _candidate_row(threshold, primary, rank=3, score=0.4),
        _candidate_row(zero, None, rank=".", ownership="zero_owner"),
        _candidate_row(ambiguous, primary, rank=4, ownership="ambiguous_multiple_owners"),
    ]
    config = AlternativeSpliceSelectionConfig.resolve()
    result, annotation, score_rows = _selection(
        primary, [accepted, exact, threshold, zero, ambiguous], rows, config
    )
    path = tmp_path / "audit.gff3"
    audit = write_alternative_splice_candidate_audit_gff3(
        candidate_rows=rows,
        selection_result=result,
        all_scored_models=[primary, accepted, exact, threshold, zero, ambiguous],
        primary_annotation=annotation,
        path=path,
    )
    assert [item.model for item in audit] == [accepted, exact, threshold]
    feature_rows = _feature_rows(path)
    transcripts = [item for item in feature_rows if item[0][2] == "mRNA"]
    assert [attrs["ID"] for _fields, attrs in transcripts] == [
        "MoSelecta_alt_candidate_000000001",
        "MoSelecta_alt_candidate_000000002",
        "MoSelecta_alt_candidate_000000003",
    ]
    assert all("Parent" not in attrs for _fields, attrs in transcripts)
    assert all(fields[1].startswith("MoSelecta:Audit:") for fields, _attrs in transcripts)
    assert transcripts[0][1]["moselecta_alt_public_transcript_id"] == (
        "MoSelecta_000001.2"
    )
    assert transcripts[1][1]["moselecta_alt_decision"] == "rejected_exact_equivalent"
    assert transcripts[2][1]["moselecta_alt_decision"] == (
        "rejected_absolute_model_selection_score"
    )
    ids = {attrs["ID"] for _fields, attrs in feature_rows}
    assert len(ids) == len(feature_rows)
    assert all(
        attrs.get("Parent", "MoSelecta_alt_candidate_").startswith(
            "MoSelecta_alt_candidate_"
        )
        for _fields, attrs in feature_rows
    )
    audit_children = [attrs for _fields, attrs in feature_rows if "Parent" in attrs]
    assert audit_children
    assert all(
        attrs["ID"].startswith(attrs["Parent"] + "__") for attrs in audit_children
    )
    parsed_audit = parse_input(path, "Audit", strict=True)
    assert len(parsed_audit.coding_models) == 3
    final_path = tmp_path / "final.gff3"
    _write_final(
        augment_annotation(annotation, result, score_rows), config, final_path
    )
    final_ids = {
        attrs["ID"] for _fields, attrs in _feature_rows(final_path)
    }
    assert ids.isdisjoint(final_ids)


def test_audit_ids_are_stable_across_presets_and_support_policies(tmp_path: Path) -> None:
    primary = _model("primary", [(0, 90), (180, 270), (360, 450)], source="P")
    exact = _model("exact", [(0, 90), (180, 240), (330, 450)], source="A")
    mosaic = _model("mosaic", [(0, 90), (180, 246), (336, 450)], source="B")
    rows = [
        _candidate_row(exact, primary, rank=1),
        _candidate_row(
            mosaic, primary, rank=2, exact_support=False, full_mosaic=True
        ),
    ]
    identity_maps = []
    policy_outputs: dict[str, tuple[object, FinalAnnotation, Path]] = {}
    for preset in ("permissive", "intermediate", "stringent"):
        for policy in ("exact", "exact-plus-mosaic"):
            config = AlternativeSpliceSelectionConfig.resolve(
                preset=preset, support_policy=policy
            )
            result, annotation, score_rows = _selection(
                primary, [exact, mosaic], rows, config
            )
            audit = build_alternative_splice_audit_candidates(
                candidate_rows=rows,
                selection_result=result,
                all_scored_models=[primary, exact, mosaic],
                primary_annotation=annotation,
            )
            identity_maps.append({
                item.model.output_transcript_id: item.audit_id for item in audit
            })
            if preset == "intermediate":
                final = tmp_path / f"{policy}.gff3"
                augmented = augment_annotation(annotation, result, score_rows)
                _write_final(augmented, config, final)
                audit_path = tmp_path / f"{policy}.audit.gff3"
                write_alternative_splice_candidate_audit_gff3(
                    candidate_rows=rows,
                    selection_result=result,
                    all_scored_models=[primary, exact, mosaic],
                    primary_annotation=annotation,
                    path=audit_path,
                )
                policy_outputs[policy] = (result, augmented, audit_path)

    assert all(identity_map == identity_maps[0] for identity_map in identity_maps)
    reordered_config = AlternativeSpliceSelectionConfig.resolve()
    reordered_result, reordered_annotation, _ = _selection(
        primary, [mosaic, exact], list(reversed(rows)), reordered_config
    )
    reordered_audit = build_alternative_splice_audit_candidates(
        candidate_rows=reversed(rows),
        selection_result=reordered_result,
        all_scored_models=[mosaic, primary, exact],
        primary_annotation=reordered_annotation,
    )
    assert {
        item.model.output_transcript_id: item.audit_id for item in reordered_audit
    } == identity_maps[0]
    exact_result, exact_annotation, exact_audit = policy_outputs["exact"]
    mosaic_result, mosaic_annotation, mosaic_audit = policy_outputs["exact-plus-mosaic"]
    assert exact_annotation.genes[0].gene_id == mosaic_annotation.genes[0].gene_id
    assert (
        exact_annotation.genes[0].primary_transcript.transcript_id
        == mosaic_annotation.genes[0].primary_transcript.transcript_id
    )
    assert len(exact_annotation.genes[0].transcripts) == 2
    assert len(mosaic_annotation.genes[0].transcripts) == 3
    exact_mosaic_attrs = next(
        attrs for fields, attrs in _feature_rows(exact_audit)
        if fields[2] == "mRNA"
        and attrs["moselecta_internal_transcript_id"] == mosaic.output_transcript_id
    )
    relaxed_mosaic_attrs = next(
        attrs for fields, attrs in _feature_rows(mosaic_audit)
        if fields[2] == "mRNA"
        and attrs["moselecta_internal_transcript_id"] == mosaic.output_transcript_id
    )
    assert exact_mosaic_attrs["moselecta_alt_support_routes"] == "full_structural_mosaic"
    assert exact_mosaic_attrs["moselecta_alt_support_policy_eligible"] == "false"
    assert exact_mosaic_attrs["moselecta_alt_decision"] == "rejected_no_qualifying_support"
    assert relaxed_mosaic_attrs["moselecta_alt_support_routes"] == "full_structural_mosaic"
    assert relaxed_mosaic_attrs["moselecta_alt_support_policy_eligible"] == "true"
    assert relaxed_mosaic_attrs["moselecta_alt_mosaic_policy_required"] == "true"
    assert relaxed_mosaic_attrs["moselecta_alt_decision"] == "accepted"
    assert len(exact_result.accepted_alternatives) == 1
    assert len(mosaic_result.accepted_alternatives) == 2


@pytest.mark.parametrize(
    "requirement",
    ("none", "protein-or-transcript", "protein", "transcript", "protein-and-transcript"),
)
def test_audit_ids_are_stable_across_evidence_requirements(
    tmp_path: Path,
    requirement: str,
) -> None:
    primary = _model("primary", [(0, 90), (180, 270), (360, 450)], source="P")
    candidate = _model("candidate", [(0, 90), (180, 240), (330, 450)], source="A")
    row = _candidate_row(
        candidate,
        primary,
        protein_exact_cds_chain_support=True,
        exact_coherent_transcript_support=True,
    )
    config = AlternativeSpliceSelectionConfig.resolve(
        evidence_requirement=requirement
    )
    result, annotation, _ = _selection(primary, [candidate], [row], config)
    audit = build_alternative_splice_audit_candidates(
        candidate_rows=[row],
        selection_result=result,
        all_scored_models=[primary, candidate],
        primary_annotation=annotation,
    )
    assert audit[0].audit_id == "MoSelecta_alt_candidate_000000001"


def test_audit_reports_external_evidence_requirement_failure(tmp_path: Path) -> None:
    primary = _model("primary", [(0, 90), (180, 270), (360, 450)], source="P")
    candidate = _model("candidate", [(0, 90), (180, 240), (330, 450)], source="A")
    row = _candidate_row(candidate, primary)
    config = AlternativeSpliceSelectionConfig.resolve(evidence_requirement="protein")
    result, annotation, _ = _selection(primary, [candidate], [row], config)
    path = tmp_path / "external.gff3"
    write_alternative_splice_candidate_audit_gff3(
        candidate_rows=[row],
        selection_result=result,
        all_scored_models=[primary, candidate],
        primary_annotation=annotation,
        path=path,
    )
    attrs = next(
        attrs for fields, attrs in _feature_rows(path) if fields[2] == "mRNA"
    )
    assert attrs["moselecta_alt_evidence_requirement"] == "protein"
    assert attrs["moselecta_alt_support_policy_eligible"] == "true"
    assert attrs["moselecta_alt_external_evidence_requirement_met"] == "false"
    assert attrs["moselecta_alt_external_evidence_requirement_routes"] == "none"
    assert attrs["moselecta_alt_prethreshold_eligible"] == "false"
    assert attrs["moselecta_alt_decision"] == "rejected_external_evidence_requirement"


def test_audit_exposes_base_union_without_relaxed_chain_containment(tmp_path: Path) -> None:
    primary = _model("primary", [(0, 100), (200, 300), (400, 500)], source="P")
    candidate = _model("candidate", [(0, 100), (200, 250), (450, 500)], source="A")
    row = _candidate_row(candidate, primary)
    config = AlternativeSpliceSelectionConfig.resolve()
    result, annotation, _score_rows = _selection(primary, [candidate], [row], config)
    path = tmp_path / "containment.gff3"
    write_alternative_splice_candidate_audit_gff3(
        candidate_rows=[row], selection_result=result,
        all_scored_models=[primary, candidate], primary_annotation=annotation, path=path,
    )
    attrs = next(attrs for fields, attrs in _feature_rows(path) if fields[2] == "mRNA")
    assert attrs["moselecta_alt_cds_base_union_containment"] == "candidate_in_primary"
    assert attrs["moselecta_alt_relaxed_chain_containment"] == "none"


def test_audit_exposes_greedy_conflict_and_cap_rejection(tmp_path: Path) -> None:
    primary = _model("primary", [(0, 90), (180, 270), (360, 450)], source="P")
    first = _model("first", [(0, 90), (180, 240), (330, 450)], source="A")
    equivalent = _model("equivalent", list(first.cds_parts), source="B")
    beyond_cap = _model("beyond_cap", [(0, 90), (180, 246), (336, 450)], source="C")
    rows = [
        _candidate_row(first, primary, rank=1),
        _candidate_row(equivalent, primary, rank=2),
        _candidate_row(beyond_cap, primary, rank=3),
    ]
    config = AlternativeSpliceSelectionConfig.resolve(max_alternatives_per_gene=1)
    result, annotation, _score_rows = _selection(
        primary, [first, equivalent, beyond_cap], rows, config
    )
    path = tmp_path / "greedy-cap.gff3"
    write_alternative_splice_candidate_audit_gff3(
        candidate_rows=rows, selection_result=result,
        all_scored_models=[primary, first, equivalent, beyond_cap],
        primary_annotation=annotation, path=path,
    )
    by_internal = {
        attrs["moselecta_internal_transcript_id"]: attrs
        for fields, attrs in _feature_rows(path) if fields[2] == "mRNA"
    }
    conflict = by_internal[equivalent.output_transcript_id]
    assert conflict["moselecta_alt_decision"] == "rejected_greedy_exact_equivalent"
    assert conflict["moselecta_alt_conflict_transcript"] == "MoSelecta_000001.2"
    assert conflict["moselecta_alt_conflict_type"] == "exact_equivalent"
    assert by_internal[beyond_cap.output_transcript_id]["moselecta_alt_decision"] == (
        "rejected_max_alternatives_per_gene"
    )


def test_audit_exposes_existing_sequence_and_intrinsic_concerns(tmp_path: Path) -> None:
    primary = _model("primary", [(0, 90), (180, 270), (360, 450)], source="P")
    sequence = _model("sequence", [(0, 90), (180, 240), (330, 450)], source="A")
    intrinsic = _model("intrinsic", [(0, 90), (180, 246), (336, 450)], source="B")
    rows = [
        _candidate_row(
            sequence,
            primary,
            rank=1,
            sequence_integrity_default_exclusion=True,
            sequence_integrity_status="evaluated_concern",
            sequence_internal_stop_concern=True,
        ),
        _candidate_row(
            intrinsic,
            primary,
            rank=2,
            intrinsic_default_exclusion=True,
            intrinsic_short_internal_cds_exon=1,
        ),
    ]
    config = AlternativeSpliceSelectionConfig.resolve()
    result, annotation, _score_rows = _selection(
        primary, [sequence, intrinsic], rows, config
    )
    path = tmp_path / "concerns.gff3"
    write_alternative_splice_candidate_audit_gff3(
        candidate_rows=rows, selection_result=result,
        all_scored_models=[primary, sequence, intrinsic],
        primary_annotation=annotation, path=path,
    )
    by_internal = {
        attrs["moselecta_internal_transcript_id"]: attrs
        for fields, attrs in _feature_rows(path) if fields[2] == "mRNA"
    }
    sequence_attrs = by_internal[sequence.output_transcript_id]
    assert sequence_attrs["moselecta_alt_decision"] == "rejected_sequence_integrity"
    assert sequence_attrs["moselecta_alt_sequence_status"] == "evaluated_concern"
    assert sequence_attrs["moselecta_alt_sequence_concerns"] == (
        "sequence_internal_stop_concern"
    )
    intrinsic_attrs = by_internal[intrinsic.output_transcript_id]
    assert intrinsic_attrs["moselecta_alt_decision"] == "rejected_intrinsic_concern"
    assert intrinsic_attrs["moselecta_alt_intrinsic_concerns"] == (
        "intrinsic_short_internal_cds_exon"
    )


def test_standardised_writer_ignores_optional_alternative_metadata(tmp_path: Path) -> None:
    primary = _model("primary", [(0, 90), (180, 270), (360, 450)], source="P")
    alternative = _model("alternative", [(0, 90), (180, 240), (330, 450)], source="A")
    primary_annotation = build_primary_annotation([primary], {primary.output_transcript_id: _score_row()})
    primary_tx = primary_annotation.genes[0].primary_transcript
    without = FinalAnnotation(genes=(FinalGene(
        primary_annotation.genes[0].gene_id,
        (primary_tx, FinalTranscript("MoSelecta_000001.2", alternative, "alternative")),
    ),))
    with_metadata = FinalAnnotation(genes=(FinalGene(
        primary_annotation.genes[0].gene_id,
        (
            primary_tx,
            FinalTranscript(
                "MoSelecta_000001.2",
                alternative,
                "alternative",
                alternative_metadata=AlternativeTranscriptMetadata(
                    1, "intermediate", "exact", "none",
                    ("exact_multisource_model_support",),
                    False, 0.9, 0.99, 0.9,
                ),
            ),
        ),
    ),))
    without_path = tmp_path / "without.gff3"
    with_path = tmp_path / "with.gff3"
    write_standardised_primary_gff3(without, without_path)
    write_standardised_primary_gff3(with_metadata, with_path)
    assert without_path.read_bytes() == with_path.read_bytes()
    assert "moselecta_alt_" not in with_path.read_text()


def test_final_annotation_preserves_numeric_order_for_many_alternatives(tmp_path: Path) -> None:
    primary = _model("primary", [(0, 100), (200, 300), (400, 500)], source="P")
    candidates = [
        _model(
            f"candidate_{index:02d}",
            [(0, 100), (200, 270 + index), (370 + index, 500)],
            source=f"S{index:02d}",
        )
        for index in range(1, 21)
    ]
    rows = [
        _candidate_row(candidate, primary, rank=index)
        for index, candidate in enumerate(candidates, start=1)
    ]
    config = AlternativeSpliceSelectionConfig.resolve(max_alternatives_per_gene=50)
    result, annotation, score_rows = _selection(primary, candidates, rows, config)
    augmented = augment_annotation(annotation, result, score_rows)
    path = tmp_path / "many.gff3"
    _write_final(augmented, config, path)
    transcript_ids = [
        attrs["ID"] for fields, attrs in _feature_rows(path) if fields[2] == "mRNA"
    ]
    assert transcript_ids == [public_transcript_id("MoSelecta_000001", index) for index in range(1, 22)]
    gene_attrs = next(attrs for fields, attrs in _feature_rows(path) if fields[2] == "gene")
    assert gene_attrs["moselecta_transcript_count"] == "21"
