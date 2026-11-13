from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

from moselecta.models import CodingModel, ExcludedModel
from moselecta.outputs import write_excluded_tsv, write_models_tsv
from moselecta.protein_evidence import terminal_prf
from moselecta.transcript_evidence import (
    TranscriptEvidenceCollapsedTranscript,
    _cdna_context_priority,
    write_transcript_collapsed_transcripts_tsv,
)


def package_root() -> Path:
    return Path(__file__).resolve().parents[1]


def run_cli(args: list[str]) -> subprocess.CompletedProcess[str]:
    root = package_root()
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    return subprocess.run(
        [sys.executable, "-m", "moselecta", *args],
        cwd=root,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def simple_model(name: str = "tx1", parts: tuple[tuple[int, int], ...] = ((9, 30),)) -> CodingModel:
    return CodingModel(
        source_name="SourceA",
        input_file=Path("input.gff3"),
        input_gene_id=f"g_{name}",
        input_transcript_id=name,
        output_gene_id=f"out_g_{name}",
        output_transcript_id=f"out_{name}",
        chrom="chr1",
        strand="+",
        cds_parts=parts,
        gene_record=None,
        transcript_record=None,
        child_records=[],
        superlocus_id="sl1",
        sublocus_id="sub1",
    )


def test_model_and_excluded_tsv_headers_match_data_rows(tmp_path: Path) -> None:
    models_path = tmp_path / "models.tsv"
    write_models_tsv([simple_model()], models_path)
    model_rows = list(csv.reader(models_path.open(), delimiter="\t"))
    assert len(model_rows[0]) == 17
    assert model_rows[0][1] == "sublocus_id"
    assert len(model_rows[1]) == len(model_rows[0])
    assert model_rows[1][1] == "sub1"

    excluded_path = tmp_path / "excluded.tsv"
    write_excluded_tsv(
        [
            ExcludedModel(
                source_name="SourceA",
                input_file=Path("input.gff3"),
                input_gene_id="g_bad",
                input_transcript_id="tx_bad",
                reason="test_exclusion",
                cds_count=2,
                cds_chroms=["chr1"],
                cds_strands=["+"],
                cds_span_bp=101,
                max_cds_intron_bp=50,
            )
        ],
        excluded_path,
    )
    excluded_rows = list(csv.reader(excluded_path.open(), delimiter="\t"))
    assert excluded_rows[0][-2:] == ["cds_span_bp", "max_cds_intron_bp"]
    assert len(excluded_rows[1]) == len(excluded_rows[0]) == 10
    assert excluded_rows[1][-2:] == ["101", "50"]


def test_minus_strand_five_prime_collapsed_coordinates_are_end_coordinates(tmp_path: Path) -> None:
    plus = TranscriptEvidenceCollapsedTranscript(
        transcript_evidence_id="plus",
        chrom="chr1",
        strand="+",
        exon_parts=((100, 120), (200, 220)),
        evidence_superlocus_id="tesl1",
        collapse_key_type="intron_chain",
        source_member_counts={"txsrc": 1},
        member_count=1,
        representative_raw_id="raw_plus",
        span_start0=100,
        span_end0=220,
        left_terminal_outer_min0=100,
        left_terminal_outer_max0=105,
        right_terminal_outer_min0=215,
        right_terminal_outer_max0=220,
    )
    minus = TranscriptEvidenceCollapsedTranscript(
        transcript_evidence_id="minus",
        chrom="chr1",
        strand="-",
        exon_parts=((100, 120), (200, 220)),
        evidence_superlocus_id="tesl2",
        collapse_key_type="intron_chain",
        source_member_counts={"txsrc": 1},
        member_count=1,
        representative_raw_id="raw_minus",
        span_start0=100,
        span_end0=220,
        left_terminal_outer_min0=100,
        left_terminal_outer_max0=105,
        right_terminal_outer_min0=215,
        right_terminal_outer_max0=220,
    )
    path = tmp_path / "collapsed.tsv"
    write_transcript_collapsed_transcripts_tsv([plus, minus], path)
    rows = {row["transcript_evidence_id"]: row for row in csv.DictReader(path.open(), delimiter="\t")}
    assert rows["plus"]["five_prime_terminal_outer_min"] == "101"
    assert rows["plus"]["five_prime_terminal_outer_max"] == "106"
    assert rows["minus"]["five_prime_terminal_outer_min"] == "215"
    assert rows["minus"]["five_prime_terminal_outer_max"] == "220"


def test_cdna_context_priority_preserves_genuine_zero_pairwise_score() -> None:
    model = simple_model()
    tx = TranscriptEvidenceCollapsedTranscript(
        transcript_evidence_id="txe1",
        chrom="chr1",
        strand="+",
        exon_parts=((9, 30),),
        evidence_superlocus_id="tesl1",
        collapse_key_type="single_exon_interval",
        source_member_counts={"txsrc": 1},
        member_count=1,
        representative_raw_id="raw1",
        span_start0=9,
        span_end0=30,
        left_terminal_outer_min0=9,
        left_terminal_outer_max0=9,
        right_terminal_outer_min0=30,
        right_terminal_outer_max0=30,
    )
    priority = _cdna_context_priority(
        model,
        tx,
        {"pairwise_transcript_score": 0.0, "association_score": 0.75},
        {"compatible_class": "single_exon_transcript_support"},
    )
    assert priority[1] == 0.0


def test_single_exon_terminal_prf_is_consistent_with_terminal_feature_count() -> None:
    query = simple_model("q", ((9, 30),))
    exact = simple_model("t", ((9, 30),))
    shifted = simple_model("u", ((10, 31),))
    exact_metrics = terminal_prf(query, exact, wobble_bp=3)
    shifted_metrics = terminal_prf(query, shifted, wobble_bp=3)
    assert exact_metrics == {"precision": 1.0, "sensitivity": 1.0, "f1": 1.0}
    assert shifted_metrics == {"precision": 0.0, "sensitivity": 0.0, "f1": 0.0}


@pytest.mark.parametrize(
    "option,value,expected",
    [
        ("--single-exon-min-reciprocal-overlap", "1.1", "must be between 0 and 1"),
        ("--long-cds-intron-percentile", "101", "must be between 0 and 100"),
    ],
)
def test_cli_rejects_invalid_clustering_and_percentile_ranges(
    tmp_path: Path, option: str, value: str, expected: str
) -> None:
    result = run_cli(
        [
            "--inputs",
            str(tmp_path / "not_read.gff3"),
            "--source-name",
            "A",
            "--out-prefix",
            str(tmp_path / "out"),
            option,
            value,
            "--progress-log-level",
            "none",
        ]
    )
    assert result.returncode != 0
    assert expected in result.stderr
    assert "Traceback" not in result.stderr


def test_superlocus_stop_exits_before_sublocus_and_final_outputs(tmp_path: Path) -> None:
    gff = tmp_path / "model.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=g1\n"
        "chr1\ttest\tmRNA\t100\t300\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\tCDS\t100\t150\t.\t+\t0\tParent=t1\n"
        "chr1\ttest\tCDS\t250\t300\t.\t+\t0\tParent=t1\n"
    )
    prefix = tmp_path / "out"
    result = run_cli(
        [
            "--inputs",
            str(gff),
            "--source-name",
            "A",
            "--out-prefix",
            str(prefix),
            "--superlocus-validation-action",
            "stop",
            "--max-superlocus-span-bp",
            "50",
            "--progress-log-level",
            "none",
        ]
    )
    assert result.returncode == 2
    assert "Superlocus validation failed after clustering" in result.stdout
    assert Path(str(prefix) + ".validation_failures.tsv").exists()
    assert Path(str(prefix) + ".superloci.tsv").exists()
    assert Path(str(prefix) + ".superlocus_members.tsv").exists()
    assert not Path(str(prefix) + ".subloci.tsv").exists()
    assert not Path(str(prefix) + ".final_selected_models.tsv").exists()
