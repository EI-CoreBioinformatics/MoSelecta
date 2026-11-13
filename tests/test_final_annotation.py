from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

from moselecta.annotation_outputs import write_standardised_primary_gff3
from moselecta.final_annotation import (
    AnnotationSummary,
    FinalAnnotation,
    FinalGene,
    FinalTranscript,
    build_primary_annotation,
    public_gene_id,
    public_transcript_id,
)
from moselecta.gffio import parse_attributes
from moselecta.models import CodingModel, FeatureRecord, SequenceIntegrityMetrics


def _record(
    feature_type: str,
    start0: int,
    end0: int,
    *,
    source: str,
    strand: str = "+",
    phase: str = ".",
    record_id: str | None = None,
    parent: str | None = None,
    record_index: int = 0,
) -> FeatureRecord:
    attrs: dict[str, str] = {}
    if record_id is not None:
        attrs["ID"] = record_id
    if parent is not None:
        attrs["Parent"] = parent
    return FeatureRecord(
        seqid="chr1",
        source=source,
        feature_type=feature_type,
        start=start0 + 1,
        end=end0,
        score=".",
        strand=strand,
        phase=phase,
        attrs=attrs,
        raw_attrs=".",
        input_file=Path(f"{source}.gff3"),
        source_name=source,
        line_number=record_index + 1,
        record_index=record_index,
        raw_line="",
    )


def _model(
    name: str,
    *,
    source: str,
    chrom: str = "chr1",
    strand: str = "+",
    parts: list[tuple[int, int]] | None = None,
    utr_extension: int = 0,
) -> CodingModel:
    cds_parts = parts or [(100, 160), (220, 280)]
    transcript_id = f"input_{name}"
    transcript_record = _record(
        "mRNA",
        cds_parts[0][0] - utr_extension,
        cds_parts[-1][1] + utr_extension,
        source=source,
        strand=strand,
        record_id=transcript_id,
        parent=f"input_gene_{name}",
    )
    transcript_record.seqid = chrom
    children: list[FeatureRecord] = []
    if utr_extension:
        utr = _record(
            "five_prime_UTR",
            cds_parts[0][0] - utr_extension,
            cds_parts[0][0],
            source=source,
            strand=strand,
            parent=transcript_id,
            record_index=1,
        )
        utr.seqid = chrom
        children.append(utr)
    for index, (start0, end0) in enumerate(cds_parts, start=1):
        exon = _record(
            "exon", start0, end0, source=source, strand=strand,
            record_id=f"{transcript_id}.exon{index}", parent=transcript_id,
            record_index=index * 2,
        )
        cds = _record(
            "CDS", start0, end0, source=source, strand=strand, phase="0",
            record_id=f"{transcript_id}.cds{index}", parent=transcript_id,
            record_index=index * 2 + 1,
        )
        exon.seqid = chrom
        cds.seqid = chrom
        children.extend((exon, cds))
    model = CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"input_gene_{name}",
        input_transcript_id=transcript_id,
        output_gene_id=f"internal_gene_{name}",
        output_transcript_id=f"internal_transcript_{name}",
        chrom=chrom,
        strand=strand,
        cds_parts=cds_parts,
        gene_record=None,
        transcript_record=transcript_record,
        child_records=children,
    )
    return model


def _feature_rows(path: Path) -> list[tuple[list[str], dict[str, str]]]:
    rows: list[tuple[list[str], dict[str, str]]] = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows.append((fields, parse_attributes(fields[8])))
    return rows


def _assert_standardized_parent_child_invariants(
    rows: list[tuple[list[str], dict[str, str]]],
) -> None:
    by_id = {attrs["ID"]: fields for fields, attrs in rows if "ID" in attrs}
    children: dict[str, list[list[str]]] = {}
    for fields, attrs in rows:
        if "Parent" in attrs:
            children.setdefault(attrs["Parent"], []).append(fields)
    for fields, attrs in rows:
        if fields[2] != "mRNA":
            continue
        gene = by_id[attrs["Parent"]]
        assert int(gene[3]) <= int(fields[3]) <= int(fields[4]) <= int(gene[4])
        transcript_children = children[attrs["ID"]]
        assert all(
            int(fields[3]) <= int(child[3]) <= int(child[4]) <= int(fields[4])
            for child in transcript_children
        )
        exons = [child for child in transcript_children if child[2] == "exon"]
        assert exons
        assert (int(fields[3]), int(fields[4])) == (
            min(int(exon[3]) for exon in exons),
            max(int(exon[4]) for exon in exons),
        )
        for child in transcript_children:
            if child[2].lower() not in {
                "cds",
                "utr",
                "five_prime_utr",
                "three_prime_utr",
                "start_codon",
                "stop_codon",
            }:
                continue
            assert any(
                int(exon[3]) <= int(child[3]) <= int(child[4]) <= int(exon[4])
                for exon in exons
            )


def _score_row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "model_selection_score_scaled": 0.328745,
        "model_structural_support_score_scaled": 0.287654,
        "protein_evidence_score_scaled": 0.536445,
        "protein_evidence_status": "protein_evidence_supporting",
        "transcript_evidence_score_scaled": 0.414817,
        "transcript_evidence_status": "partial_transcript_support",
        "sequence_integrity_status": "evaluated_no_concern",
    }
    row.update(overrides)
    return row


def _primary_annotation_with_scores(model: CodingModel) -> FinalAnnotation:
    return build_primary_annotation(
        [model],
        {model.output_transcript_id: _score_row()},
    )


def test_public_id_helpers_share_one_prefix() -> None:
    assert public_gene_id(1) == "MoSelecta_000001"
    assert public_gene_id(123456) == "MoSelecta_123456"
    assert public_transcript_id("MoSelecta_000001", 1) == "MoSelecta_000001.1"
    assert public_transcript_id("MoSelecta_000042", 3) == "MoSelecta_000042.3"
    with pytest.raises(ValueError):
        public_gene_id(0)
    with pytest.raises(ValueError):
        public_transcript_id("MoSelecta_000001", 0)
    with pytest.raises(ValueError):
        public_transcript_id("Other_000001", 1)


def test_final_gene_rejects_transcript_ids_from_another_gene() -> None:
    model = _model("tx", source="src")
    with pytest.raises(ValueError, match="must begin with primary transcript"):
        FinalGene("MoSelecta_000001", (FinalTranscript("MoSelecta_000002.1", model),))


def test_primary_annotation_ids_are_deterministic_and_do_not_mutate_models() -> None:
    later = _model("later", source="B", parts=[(500, 560), (620, 680)])
    first_tie_b = _model("tie_b", source="B", parts=[(100, 160), (220, 280)])
    first_tie_a = _model("tie_a", source="A", parts=[(100, 160), (220, 280)])
    before = {
        id(model): (model.output_gene_id, model.output_transcript_id)
        for model in (later, first_tie_b, first_tie_a)
    }

    first = build_primary_annotation([later, first_tie_b, first_tie_a])
    repeated = build_primary_annotation([first_tie_a, later, first_tie_b])

    assert [gene.gene_id for gene in first.genes] == [
        "MoSelecta_000001",
        "MoSelecta_000002",
        "MoSelecta_000003",
    ]
    assert [
        gene.primary_transcript.source_model.output_transcript_id
        for gene in first.genes
    ] == [
        "internal_transcript_tie_a",
        "internal_transcript_tie_b",
        "internal_transcript_later",
    ]
    assert [
        (gene.gene_id, gene.primary_transcript.transcript_id)
        for gene in first.genes
    ] == [
        (gene.gene_id, gene.primary_transcript.transcript_id)
        for gene in repeated.genes
    ]
    assert before == {
        id(model): (model.output_gene_id, model.output_transcript_id)
        for model in (later, first_tie_b, first_tie_a)
    }


def test_future_transcript_addition_cannot_change_existing_primary_public_ids(
    tmp_path: Path,
) -> None:
    primary_model = _model("primary", source="Annevo")
    alternative_model = _model(
        "alternative",
        source="Tiberius",
        parts=[(100, 160), (240, 300)],
    )
    primary_only = build_primary_annotation([primary_model])
    original_gene = primary_only.genes[0]
    original_primary = original_gene.primary_transcript

    expanded_gene = FinalGene(
        gene_id=original_gene.gene_id,
        transcripts=(
            original_primary,
            FinalTranscript(
                transcript_id=public_transcript_id(original_gene.gene_id, 2),
                source_model=alternative_model,
                transcript_role="alternative",
            ),
        ),
    )
    expanded = FinalAnnotation(genes=(expanded_gene,))

    assert expanded.genes[0].gene_id == original_gene.gene_id
    assert expanded.genes[0].primary_transcript is original_primary
    assert expanded.genes[0].primary_transcript.transcript_id == "MoSelecta_000001.1"
    assert expanded.genes[0].primary_transcript.source_model is primary_model

    path = tmp_path / "expanded.gff3"
    write_standardised_primary_gff3(expanded, path)
    rows = _feature_rows(path)
    gene_attrs = next(attrs for fields, attrs in rows if fields[2] == "gene")
    transcript_attrs = [attrs for fields, attrs in rows if fields[2] == "mRNA"]
    assert gene_attrs["moselecta_transcript_count"] == "2"
    assert [attrs["moselecta_transcript_role"] for attrs in transcript_attrs] == [
        "primary",
        "alternative",
    ]


def test_standardised_gff_uses_only_selected_representative_provenance(tmp_path: Path) -> None:
    selected = _model("selected", source="Annevo", utr_extension=25)
    unretained_equivalent = _model("equivalent", source="Tiberius", utr_extension=50)
    assert selected.cds_parts == unretained_equivalent.cds_parts
    annotation = _primary_annotation_with_scores(selected)
    path = tmp_path / "standardised.gff3"

    write_standardised_primary_gff3(annotation, path)
    rows = _feature_rows(path)
    gene_fields, gene_attrs = next(row for row in rows if row[0][2] == "gene")
    transcript_fields, transcript_attrs = next(row for row in rows if row[0][2] == "mRNA")

    assert gene_fields[3:5] == ["76", "280"]
    assert gene_attrs == {
        "ID": "MoSelecta_000001",
        "moselecta_primary_transcript_id": "MoSelecta_000001.1",
        "moselecta_transcript_count": "1",
    }
    assert transcript_fields[3:5] == ["76", "280"]
    assert transcript_attrs == {
        "ID": "MoSelecta_000001.1",
        "Parent": "MoSelecta_000001",
        "moselecta_transcript_role": "primary",
        "original_source": "Annevo",
        "original_gene_id": "input_gene_selected",
        "original_transcript_id": "input_selected",
        "moselecta_internal_gene_id": "internal_gene_selected",
        "moselecta_internal_transcript_id": "internal_transcript_selected",
        "model_selection_score_scaled": "0.328745",
        "model_structural_support_score_scaled": "0.287654",
        "protein_evidence_score_scaled": "0.536445",
        "transcript_evidence_score_scaled": "0.414817",
        "sequence_integrity_status": "evaluated_no_concern",
        "cds_phase_output_status": "preserved_numeric_source",
    }
    assert "Tiberius" not in path.read_text()
    assert "internal_transcript_equivalent" not in path.read_text()


def test_standardised_gff_ids_parents_and_source_features_are_valid(tmp_path: Path) -> None:
    models = [
        _model("one", source="Annevo", utr_extension=20),
        _model("two", source="Helixer", parts=[(500, 560), (620, 680)]),
    ]
    annotation = build_primary_annotation(
        models,
        {model.output_transcript_id: _score_row() for model in models},
    )
    path = tmp_path / "standardised.gff3"
    write_standardised_primary_gff3(annotation, path)
    rows = _feature_rows(path)

    ids = [attrs["ID"] for _fields, attrs in rows]
    assert len(ids) == len(set(ids))
    child_counts: dict[tuple[str, str], int] = {}
    for fields, attrs in rows:
        if fields[2] in {"gene", "mRNA"}:
            continue
        key = (attrs["Parent"], fields[2])
        child_counts[key] = child_counts.get(key, 0) + 1
        assert attrs["ID"] == f"{attrs['Parent']}.{fields[2]}{child_counts[key]}"
        assert "__" not in attrs["ID"]
    id_set = set(ids)
    for _fields, attrs in rows:
        if "Parent" in attrs:
            assert attrs["Parent"] in id_set
    assert sum(fields[2] == "gene" for fields, _attrs in rows) == 2
    assert sum(fields[2] == "mRNA" for fields, _attrs in rows) == 2
    cds_rows = [fields for fields, _attrs in rows if fields[2] == "CDS"]
    assert [fields[7] for fields in cds_rows] == ["0", "0", "0", "0"]
    assert {fields[1] for fields in cds_rows} == {"Annevo", "Helixer"}
    assert {
        fields[1] for fields, _attrs in rows if fields[2] == "gene"
    } == {"MoSelecta"}
    assert {
        fields[1] for fields, _attrs in rows if fields[2] == "mRNA"
    } == {"MoSelecta:Annevo", "MoSelecta:Helixer"}


def test_final_gff_preserves_only_explicitly_supported_child_types(
    tmp_path: Path,
) -> None:
    model = _model("children", source="Prediction", utr_extension=20)
    transcript_id = model.input_transcript_id
    model.child_records.extend([
        _record(
            "start_codon", 100, 103, source="Prediction", phase="0",
            parent=transcript_id, record_index=20,
        ),
        _record(
            "stop_codon", 277, 280, source="Prediction", phase="0",
            parent=transcript_id, record_index=21,
        ),
        _record(
            "intron", 160, 220, source="Prediction",
            parent=transcript_id, record_index=22,
        ),
        _record(
            "transcription_start_site", 79, 80, source="Prediction",
            parent=transcript_id, record_index=23,
        ),
        _record(
            "transcription_end_site", 280, 281, source="Prediction",
            parent=transcript_id, record_index=24,
        ),
        _record(
            "source_translation_hint", 140, 143, source="Prediction",
            parent=transcript_id, record_index=25,
        ),
    ])
    path = tmp_path / "supported_children.gff3"
    write_standardised_primary_gff3(_primary_annotation_with_scores(model), path)
    rows = _feature_rows(path)
    child_types = {
        fields[2]
        for fields, _attrs in rows
        if fields[2] not in {"gene", "mRNA"}
    }
    assert child_types == {
        "five_prime_UTR",
        "exon",
        "CDS",
        "start_codon",
        "stop_codon",
    }
    assert all(int(fields[3]) <= int(fields[4]) for fields, _attrs in rows)


def test_standardised_gff_omits_only_unavailable_support_values(tmp_path: Path) -> None:
    unavailable = _model("unavailable", source="Annevo")
    configured_zero = _model(
        "configured_zero",
        source="Tiberius",
        parts=[(500, 560), (620, 680)],
    )
    annotation = build_primary_annotation(
        [unavailable, configured_zero],
        {
            unavailable.output_transcript_id: _score_row(
                protein_evidence_score_scaled=0.0,
                protein_evidence_status="protein_evidence_not_configured",
                transcript_evidence_score_scaled=0.0,
                transcript_evidence_status="transcript_evidence_not_configured",
                sequence_integrity_status="genome_not_configured",
            ),
            configured_zero.output_transcript_id: _score_row(
                protein_evidence_score_scaled=0.0,
                protein_evidence_status="no_overlapping_protein_evidence",
                transcript_evidence_score_scaled=0.0,
                transcript_evidence_status="weak_transcript_support",
                sequence_integrity_status="not_evaluable_missing_contig",
            ),
        },
    )
    path = tmp_path / "availability.gff3"
    write_standardised_primary_gff3(annotation, path)
    transcripts = [attrs for fields, attrs in _feature_rows(path) if fields[2] == "mRNA"]

    unavailable_attrs, configured_attrs = transcripts
    assert "protein_evidence_score_scaled" not in unavailable_attrs
    assert "transcript_evidence_score_scaled" not in unavailable_attrs
    assert "sequence_integrity_status" not in unavailable_attrs
    assert unavailable_attrs["model_selection_score_scaled"] == "0.328745"
    assert unavailable_attrs["model_structural_support_score_scaled"] == "0.287654"

    assert configured_attrs["protein_evidence_score_scaled"] == "0.000000"
    assert configured_attrs["transcript_evidence_score_scaled"] == "0.000000"
    assert configured_attrs["sequence_integrity_status"] == "not_evaluable_missing_contig"


def test_standardised_gff_escapes_representative_provenance_attributes(tmp_path: Path) -> None:
    model = _model("escaped", source="Annevo")
    model.input_gene_id = "gene; odd=1"
    model.input_transcript_id = "transcript value,%"
    model.output_gene_id = "internal gene;1"
    model.output_transcript_id = "internal transcript;1"
    annotation = _primary_annotation_with_scores(model)
    path = tmp_path / "escaped.gff3"

    write_standardised_primary_gff3(annotation, path)
    raw_text = path.read_text()
    transcript_attrs = next(
        attrs for fields, attrs in _feature_rows(path) if fields[2] == "mRNA"
    )
    assert transcript_attrs["original_gene_id"] == "gene; odd=1"
    assert transcript_attrs["original_transcript_id"] == "transcript value,%"
    assert transcript_attrs["moselecta_internal_gene_id"] == "internal gene;1"
    assert transcript_attrs["moselecta_internal_transcript_id"] == "internal transcript;1"
    assert "gene%3B%20odd%3D1" in raw_text
    assert "transcript%20value%2C%25" in raw_text


def test_final_gene_span_is_union_of_transcript_source_models() -> None:
    primary = _model("primary", source="Annevo", utr_extension=10)
    future_transcript = _model(
        "future", source="Tiberius", parts=[(80, 160), (220, 320)], utr_extension=20,
    )
    gene = FinalGene(
        gene_id="MoSelecta_000001",
        transcripts=(
            FinalTranscript("MoSelecta_000001.1", primary),
            FinalTranscript(
                "MoSelecta_000001.2",
                future_transcript,
                transcript_role="alternative",
            ),
        ),
    )
    assert gene.start0 == 60
    assert gene.end0 == 340
    annotation = FinalAnnotation((gene,))
    assert annotation.gene_count == 1
    assert annotation.transcript_count == 2
    assert annotation.primary_transcript_count == 1
    assert annotation.alternative_transcript_count == 1


def _integrity_metrics(
    *, start_status: str = "present", stop_status: str = "included_in_cds"
) -> SequenceIntegrityMetrics:
    missing = start_status == "missing" or stop_status == "missing"
    return SequenceIntegrityMetrics(
        status="evaluated_with_concern" if missing else "evaluated_no_concern",
        start_codon="CCC" if start_status == "missing" else "ATG",
        start_codon_status=start_status,
        stop_codon="CCC" if stop_status == "missing" else "TAA",
        stop_codon_status=stop_status,
        internal_stop_count=0,
        phase_status="consistent",
        phase_mismatch_count=0,
        noncanonical_splice_intron_count=0,
        noncanonical_splice_intron_fraction=0.0,
        noncanonical_splice_details=".",
        has_ambiguous_bases=False,
        penalty_score=0.0,
        penalty_components=".",
        flags=".",
    )


def test_standardized_edge_intron_bounds_match_emitted_exons(
    tmp_path: Path,
) -> None:
    cases = (
        ("gene_003011", "+", (0, 808), (568, 808), (0, 568), (569, 808)),
        ("gene_009062", "-", (0, 346), (0, 81), (81, 346), (1, 81)),
        ("gene_009077", "+", (0, 536), (136, 536), (0, 155), (137, 536)),
        ("gene_009110", "-", (881, 1257), (881, 1228), (1228, 1257), (882, 1228)),
    )
    models: list[CodingModel] = []
    for name, strand, source_span, exon_span, intron_span, _expected in cases:
        model = _model(name, source="Prediction", strand=strand, parts=[exon_span])
        model.transcript_record = _record(
            "mRNA",
            *source_span,
            source="Prediction",
            strand=strand,
            record_id=f"{name}.t1",
            parent=name,
        )
        model.child_records.append(
            _record(
                "intron",
                *intron_span,
                source="Prediction",
                strand=strand,
                parent=f"{name}.t1",
                record_index=99,
            )
        )
        if name == "gene_009077":
            model.cds_terminus_refined = "upstream_start"
        models.append(model)

    before = sum(
        (model.model_start0, model.model_end0)
        != (
            min(record.start0 for record in model.child_records if record.feature_type == "exon"),
            max(record.end0 for record in model.child_records if record.feature_type == "exon"),
        )
        for model in models
    )
    assert before == 4

    annotation = build_primary_annotation(models)
    path = tmp_path / "edge_introns.gff3"
    write_standardised_primary_gff3(annotation, path)
    rows = _feature_rows(path)
    transcript_rows = {
        attrs["original_gene_id"]: (fields, attrs)
        for fields, attrs in rows
        if fields[2] == "mRNA"
    }
    for name, _strand, _source, _exon, _intron, expected in cases:
        fields, attrs = transcript_rows[f"input_gene_{name}"]
        assert tuple(map(int, fields[3:5])) == expected
        if name == "gene_009077":
            assert attrs["cds_terminus_refined"] == "upstream_start"

    children_by_parent: dict[str, list[list[str]]] = {}
    for fields, attrs in rows:
        if "Parent" in attrs and fields[2] != "mRNA":
            children_by_parent.setdefault(attrs["Parent"], []).append(fields)
    after = 0
    for fields, attrs in rows:
        if fields[2] != "mRNA":
            continue
        exons = [
            child
            for child in children_by_parent[attrs["ID"]]
            if child[2] == "exon"
        ]
        exon_bounds = (min(int(row[3]) for row in exons), max(int(row[4]) for row in exons))
        after += tuple(map(int, fields[3:5])) != exon_bounds
    assert after == 0


def test_standardized_output_parent_child_and_multitranscript_bounds(
    tmp_path: Path,
) -> None:
    left = _model(
        "left",
        source="Prediction",
        parts=[(100, 150), (200, 250)],
        utr_extension=20,
    )
    right = _model(
        "right",
        source="Prediction",
        parts=[(300, 350), (400, 460)],
    )
    left.child_records.append(
        _record(
            "intron", 150, 200, source="Prediction", parent=left.input_transcript_id
        )
    )
    gene_id = "MoSelecta_000001"
    annotation = FinalAnnotation((FinalGene(gene_id, (
        FinalTranscript(public_transcript_id(gene_id, 1), left),
        FinalTranscript(public_transcript_id(gene_id, 2), right, transcript_role="alternative"),
    )),))
    path = tmp_path / "parents.gff3"
    write_standardised_primary_gff3(annotation, path)
    rows = _feature_rows(path)
    gene = next(fields for fields, _attrs in rows if fields[2] == "gene")
    assert gene[3:5] == ["81", "460"]

    by_id = {attrs["ID"]: fields for fields, attrs in rows if "ID" in attrs}
    children: dict[str, list[list[str]]] = {}
    for fields, attrs in rows:
        if "Parent" in attrs:
            children.setdefault(attrs["Parent"], []).append(fields)
    for fields, attrs in rows:
        if fields[2] != "mRNA":
            continue
        gene_fields = by_id[attrs["Parent"]]
        assert int(gene_fields[3]) <= int(fields[3]) <= int(fields[4]) <= int(gene_fields[4])
        transcript_children = children[attrs["ID"]]
        assert all(
            int(fields[3]) <= int(child[3]) <= int(child[4]) <= int(fields[4])
            for child in transcript_children
        )
        exons = [child for child in transcript_children if child[2] == "exon"]
        assert (int(fields[3]), int(fields[4])) == (
            min(int(exon[3]) for exon in exons),
            max(int(exon[4]) for exon in exons),
        )
    _assert_standardized_parent_child_invariants(rows)
    left_exons = [
        child for child in children[public_transcript_id(gene_id, 1)] if child[2] == "exon"
    ]
    assert [(row[3], row[4]) for row in left_exons] == [("81", "150"), ("201", "250")]


@pytest.mark.parametrize(
    ("start_status", "stop_status", "expected"),
    (
        ("missing", "included_in_cds", {"stop_codon"}),
        ("present", "missing", {"start_codon"}),
        ("present", "included_in_cds", {"start_codon", "stop_codon"}),
    ),
)
@pytest.mark.parametrize("strand", ("+", "-"))
def test_standardized_output_suppresses_only_explicitly_missing_codon_children(
    tmp_path: Path,
    start_status: str,
    stop_status: str,
    expected: set[str],
    strand: str,
) -> None:
    model = _model(
        "codons", source="Prediction", parts=[(100, 160)], strand=strand
    )
    model.child_records.extend((
        _record(
            "start_codon", 100, 103, source="Prediction", phase="0", strand=strand
        ),
        _record(
            "stop_codon", 157, 160, source="Prediction", phase="0", strand=strand
        ),
    ))
    model.sequence_integrity = _integrity_metrics(
        start_status=start_status, stop_status=stop_status
    )
    path = tmp_path / f"codons_{strand}_{start_status}_{stop_status}.gff3"
    write_standardised_primary_gff3(build_primary_annotation([model]), path)
    emitted = {
        fields[2]
        for fields, _attrs in _feature_rows(path)
        if fields[2] in {"start_codon", "stop_codon"}
    }
    assert emitted == expected


def test_standardised_gff_is_byte_deterministic(tmp_path: Path) -> None:
    models = [
        _model("later", source="Tiberius", parts=[(500, 560), (620, 680)]),
        _model("first", source="Annevo"),
    ]
    score_rows = {
        model.output_transcript_id: _score_row() for model in models
    }
    first = build_primary_annotation(models, score_rows)
    repeated = build_primary_annotation(reversed(models), score_rows)
    first_path = tmp_path / "first.gff3"
    repeated_path = tmp_path / "repeated.gff3"
    write_standardised_primary_gff3(first, first_path)
    write_standardised_primary_gff3(repeated, repeated_path)
    assert first_path.read_bytes() == repeated_path.read_bytes()


def test_annotation_summary_is_typed_and_immutable() -> None:
    summary = AnnotationSummary(model_selection_score_scaled=0.5)
    with pytest.raises((AttributeError, TypeError)):
        summary.model_selection_score_scaled = 0.6  # type: ignore[misc]


def test_duplicate_frozen_internal_transcript_ids_are_rejected() -> None:
    first = _model("duplicate", source="A")
    second = _model("duplicate", source="B", parts=[(500, 560), (620, 680)])
    with pytest.raises(ValueError, match="duplicate internal transcript IDs"):
        build_primary_annotation([first, second])


def test_primary_annotation_rejects_missing_requested_score_rows() -> None:
    model = _model("missing_score", source="Annevo")
    with pytest.raises(ValueError, match="missing canonical score rows"):
        build_primary_annotation([model], {})


@pytest.mark.parametrize("pruning_level", ("none", "auto"))
def test_cli_standardised_output_does_not_change_existing_scientific_files(
    tmp_path: Path,
    pruning_level: str,
) -> None:
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    inputs = [
        root / "tests" / "data" / "alternative_splice_smoke_a.gff3",
        root / "tests" / "data" / "alternative_splice_smoke_b.gff3",
    ]

    def run(prefix: Path, *, enabled: bool) -> None:
        command = [
            sys.executable, "-m", "moselecta",
            "--inputs", *(str(path) for path in inputs),
            "--source-name", "A", "B",
            "--out-prefix", str(prefix),
            "--pruning-level", pruning_level,
            "--output-profile", "minimal",
            "--progress-log-level", "none",
            "--write-output", "final_selected_models_tsv,final_pruned_selected_models_tsv",
        ]
        if enabled:
            command.extend([
                "--write-output", "standardised_primary_models_gff3",
            ])
        else:
            command.extend([
                "--skip-output", "standardised_primary_models_gff3",
            ])
        subprocess.run(
            command,
            cwd=root,
            env=env,
            check=True,
            text=True,
            capture_output=True,
        )

    disabled_prefix = tmp_path / "disabled"
    enabled_prefix = tmp_path / "enabled"
    run(disabled_prefix, enabled=False)
    run(enabled_prefix, enabled=True)
    disabled_files = {
        path.name.removeprefix("disabled."): path
        for path in tmp_path.glob("disabled.*")
    }
    enabled_files = {
        path.name.removeprefix("enabled."): path
        for path in tmp_path.glob("enabled.*")
    }
    assert set(enabled_files) - set(disabled_files) == {
        "standardised_primary_models.gff3"
    }
    for suffix, path in disabled_files.items():
        if suffix == "run_configuration.json":
            continue
        assert path.read_bytes() == enabled_files[suffix].read_bytes(), suffix

    public_rows = _feature_rows(enabled_files["standardised_primary_models.gff3"])
    public_transcripts = [
        attrs["moselecta_internal_transcript_id"]
        for fields, attrs in public_rows
        if fields[2] == "mRNA"
    ]
    legacy_suffix = (
        "final_selected_models.tsv"
        if pruning_level == "none"
        else "final_selected_models.pruned.tsv"
    )
    with enabled_files[legacy_suffix].open() as handle:
        legacy_transcripts = [
            row["output_transcript_id"]
            for row in csv.DictReader(handle, delimiter="\t")
        ]
    assert public_transcripts == legacy_transcripts


def test_cli_cds_only_flag_changes_only_final_gff3_serialization(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    inputs = [
        root / "tests" / "data" / "alternative_splice_smoke_a.gff3",
        root / "tests" / "data" / "alternative_splice_smoke_b.gff3",
    ]

    def run(prefix: Path, *, cds_only: bool) -> dict[str, Path]:
        command = [
            sys.executable, "-m", "moselecta",
            "--inputs", *(str(path) for path in inputs),
            "--source-name", "A", "B",
            "--out-prefix", str(prefix),
            "--pruning-level", "auto",
            "--output-profile", "standard",
            "--progress-log-level", "none",
        ]
        if cds_only:
            command.append("--final-gff3-cds-only-structure")
        subprocess.run(
            command,
            cwd=root,
            env=env,
            check=True,
            text=True,
            capture_output=True,
        )
        return {
            path.name.removeprefix(prefix.name + "."): path
            for path in tmp_path.glob(prefix.name + ".*")
        }

    default_files = run(tmp_path / "default", cds_only=False)
    cds_only_files = run(tmp_path / "cds_only", cds_only=True)
    assert set(default_files) == set(cds_only_files)
    affected_gff3 = {
        "final_selected_models.gff3",
        "final_selected_models.pruned.gff3",
        "standardised_primary_models.gff3",
        "final_annotation.gff3",
    }
    for suffix in set(default_files) - affected_gff3 - {"run_configuration.json"}:
        assert default_files[suffix].read_bytes() == cds_only_files[suffix].read_bytes(), suffix

    for suffix in ("standardised_primary_models.gff3", "final_annotation.gff3"):
        default_rows = _feature_rows(default_files[suffix])
        cds_only_rows = _feature_rows(cds_only_files[suffix])
        _assert_standardized_parent_child_invariants(default_rows)
        default_transcripts = {
            attrs["ID"]: attrs["moselecta_internal_transcript_id"]
            for fields, attrs in default_rows if fields[2] == "mRNA"
        }
        cds_only_transcripts = {
            attrs["ID"]: attrs["moselecta_internal_transcript_id"]
            for fields, attrs in cds_only_rows if fields[2] == "mRNA"
        }
        assert cds_only_transcripts == default_transcripts
        for internal_id in default_transcripts.values():
            assert internal_id.endswith("__mRNA")
        for fields, attrs in default_rows:
            if fields[2] == "mRNA":
                assert attrs["moselecta_internal_gene_id"].endswith("__gene")
        for rows in (default_rows, cds_only_rows):
            for fields, attrs in rows:
                if fields[2] not in {"gene", "mRNA"}:
                    assert attrs["ID"].startswith(attrs["Parent"] + ".")
                    assert "__" not in attrs["ID"]
        assert {
            attrs["ID"] for fields, attrs in cds_only_rows if fields[2] == "gene"
        } == {
            attrs["ID"] for fields, attrs in default_rows if fields[2] == "gene"
        }

        def cds_geometry(rows: list[tuple[list[str], dict[str, str]]]):
            return sorted(
                (attrs["Parent"], fields[3], fields[4], fields[6], fields[7])
                for fields, attrs in rows if fields[2] == "CDS"
            )

        assert cds_geometry(cds_only_rows) == cds_geometry(default_rows)
        assert {
            fields[2] for fields, _attrs in cds_only_rows
        } == {"gene", "mRNA", "exon", "CDS"}
        for parent in cds_only_transcripts:
            exon_ids = [
                attrs["ID"] for fields, attrs in cds_only_rows
                if fields[2] == "exon" and attrs["Parent"] == parent
            ]
            cds_ids = [
                attrs["ID"] for fields, attrs in cds_only_rows
                if fields[2] == "CDS" and attrs["Parent"] == parent
            ]
            assert exon_ids == [f"{parent}.cds_exon{i}" for i in range(1, len(exon_ids) + 1)]
            assert cds_ids == [f"{parent}.CDS{i}" for i in range(1, len(cds_ids) + 1)]


@pytest.mark.parametrize("cds_only", (False, True))
def test_cli_only_main_public_gff3_use_dot_separated_child_ids(
    tmp_path: Path,
    cds_only: bool,
) -> None:
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    prefix = tmp_path / "run"
    command = [
        sys.executable, "-m", "moselecta",
        "--inputs",
        str(root / "tests" / "data" / "alternative_splice_smoke_a.gff3"),
        str(root / "tests" / "data" / "alternative_splice_smoke_b.gff3"),
        "--source-name", "A", "B",
        "--out-prefix", str(prefix),
        "--pruning-level", "auto",
        "--output-profile", "debug",
        "--progress-log-level", "none",
    ]
    if cds_only:
        command.append("--final-gff3-cds-only-structure")
    subprocess.run(command, cwd=root, env=env, check=True, text=True, capture_output=True)

    main_public = {"standardised_primary_models.gff3", "final_annotation.gff3"}
    dot_separated = main_public | {"final_annotation.with_utrs.gff3"}
    gff3_files = {
        path.name.removeprefix("run."): path for path in tmp_path.glob("run.*.gff3")
    }
    assert main_public <= set(gff3_files)
    assert {"final_selected_models.gff3", "selected_models.gff3"} <= set(gff3_files)

    def child_rows(path: Path) -> list[dict[str, str]]:
        return [
            attrs for fields, attrs in _feature_rows(path)
            if fields[2] in {"exon", "CDS"} and "Parent" in attrs and "ID" in attrs
        ]

    for suffix, path in gff3_files.items():
        children = child_rows(path)
        if suffix in dot_separated:
            assert children, suffix
            for attrs in children:
                assert attrs["ID"].startswith(attrs["Parent"] + "."), suffix
                assert "__" not in attrs["ID"], suffix
        else:
            for attrs in children:
                assert attrs["ID"].startswith(attrs["Parent"] + "__"), suffix
                assert attrs["Parent"].endswith("__mRNA"), suffix
    assert any(child_rows(gff3_files["final_selected_models.gff3"]))


@pytest.mark.parametrize("stop_after", ("scored-candidates", "monosublocus", "backfill"))
def test_cli_standardised_primary_output_is_final_stage_only(
    tmp_path: Path,
    stop_after: str,
) -> None:
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    output_path = tmp_path / f"{stop_after}.public.gff3"
    completed = subprocess.run(
        [
            sys.executable, "-m", "moselecta",
            "--inputs",
            str(root / "tests" / "data" / "alternative_splice_smoke_a.gff3"),
            str(root / "tests" / "data" / "alternative_splice_smoke_b.gff3"),
            "--source-name", "A", "B",
            "--out-prefix", str(tmp_path / stop_after),
            "--stop-after", stop_after,
            "--output-profile", "minimal",
            "--standardised-primary-models-gff3", str(output_path),
            "--progress-log-level", "none",
        ],
        cwd=root,
        env=env,
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0
    assert not output_path.exists()
