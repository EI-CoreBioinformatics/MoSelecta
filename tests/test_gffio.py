from __future__ import annotations

from pathlib import Path
from typing import get_type_hints

import pytest

from moselecta.gffio import first_parent, make_unique_source_names, parse_attributes, parse_input, read_gff_records


def test_gffio_optional_annotations_resolve() -> None:
    assert get_type_hints(first_parent)["return"] == str | None


def test_parse_attributes_supports_gff3_and_gtf() -> None:
    assert parse_attributes("ID=tx1;Parent=gene1;Note=a%20b") == {"ID": "tx1", "Parent": "gene1", "Note": "a b"}
    assert parse_attributes('gene_id "gene1"; transcript_id "tx1";') == {"gene_id": "gene1", "transcript_id": "tx1"}


def test_make_unique_source_names_sanitizes_and_deduplicates(tmp_path: Path) -> None:
    inputs = [tmp_path / "a.gff3", tmp_path / "b.gff3", tmp_path / "c.gff3"]
    assert make_unique_source_names(inputs, ["my source", "my/source", "my_source"]) == ["my_source", "my_source_2", "my_source_3"]


def test_read_gff_records_strict_rejects_malformed_line(tmp_path: Path) -> None:
    gff = tmp_path / "bad.gff3"
    gff.write_text("chr1\tsrc\tgene\t1\t10\t.\t+\t.\n")
    with pytest.raises(ValueError):
        read_gff_records(gff, "src", strict=True)


def test_parse_input_extracts_coding_model_from_gtf_style_records(tmp_path: Path) -> None:
    gtf = tmp_path / "model.gtf"
    gtf.write_text(
        'chr1\tsrc\ttranscript\t100\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tsrc\tCDS\t100\t150\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tsrc\tCDS\t200\t300\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
    )
    parsed = parse_input(gtf, "GTF", strict=True)
    assert len(parsed.coding_models) == 1
    model = parsed.coding_models[0]
    assert model.input_gene_id == "g1"
    assert model.input_transcript_id == "t1"
    assert model.cds_parts == [(99, 150), (199, 300)]
    assert not parsed.excluded_models
