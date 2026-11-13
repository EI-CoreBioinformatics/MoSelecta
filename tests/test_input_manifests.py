from __future__ import annotations

from pathlib import Path
import os
import subprocess
import sys

import pytest

from moselecta.input_manifests import resolve_input_family, split_input_tokens


def _file(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("test\n")
    return path


def _resolve(
    values,
    labels=None,
    *,
    family="model",
    input_option="--inputs",
    label_option="--source-name",
    comma=False,
    prefix=None,
):
    return resolve_input_family(
        values,
        labels,
        family_name=family,
        input_option=input_option,
        label_option=label_option,
        comma_separated_inputs=comma,
        comma_separated_labels=comma,
        fallback_prefix=prefix,
    )


def test_model_manifest_expands_relative_paths_labels_comments_and_order(
    tmp_path: Path,
) -> None:
    direct_a = _file(tmp_path / "direct_a.gff3")
    direct_z = _file(tmp_path / "direct_z.gff3")
    b = _file(tmp_path / "annotations" / "b file.gff3")
    c = _file(tmp_path / "annotations" / "c.gff3")
    manifest = tmp_path / "manifests" / "models.tsv"
    manifest.parent.mkdir()
    manifest.write_text(
        "  # comment\n\n"
        "../annotations/b file.gff3\tBRAKER\n"
        "../annotations/c.gff3\tHelixer\n"
    )

    resolved = _resolve([direct_a, f"@{manifest}", direct_z])

    assert resolved.paths == [direct_a, b.resolve(), c.resolve(), direct_z]
    assert resolved.labels == ("direct_a", "BRAKER", "Helixer", "direct_z")
    assert resolved.entries[1].origin == f"manifest:{manifest.resolve()}:3"
    assert resolved.expansion_messages == (
        f"Expanded @{manifest.resolve()}: 2 model inputs",
    )


def test_explicit_model_labels_apply_to_complete_expansion_and_override_manifest(
    tmp_path: Path,
) -> None:
    paths = [_file(tmp_path / f"{name}.gff3") for name in "AbcdZ"]
    manifest = tmp_path / "models.tsv"
    manifest.write_text(f"{paths[1]}\tignored1\n{paths[2]}\tignored2\n{paths[3]}\n")

    resolved = _resolve(
        [paths[0], f"@{manifest}", paths[4]],
        ["A", "B", "C", "D", "Z"],
    )
    assert resolved.labels == ("A", "B", "C", "D", "Z")

    with pytest.raises(SystemExit, match=r"--source-name supplied 2 labels for 5 expanded model inputs"):
        _resolve([paths[0], f"@{manifest}", paths[4]], ["A", "Z"])


@pytest.mark.parametrize(
    ("family", "input_option", "label_option", "prefix"),
    [
        ("transcript evidence", "--evidence-transcripts", "--evidence-transcript-source-name", "transcript"),
        ("protein evidence", "--evidence-proteins", "--evidence-protein-source-name", "protein"),
    ],
)
def test_evidence_manifest_label_precedence_and_fallback(
    tmp_path: Path,
    family: str,
    input_option: str,
    label_option: str,
    prefix: str,
) -> None:
    paths = [_file(tmp_path / f"{prefix}_{index}.gff3") for index in range(1, 4)]
    manifest = tmp_path / f"{prefix}.tsv"
    manifest.write_text(f"{paths[0]}\tlabel_one\n{paths[1]}\n{paths[2]}\tlabel_three\n")
    kwargs = dict(
        family=family,
        input_option=input_option,
        label_option=label_option,
        comma=True,
        prefix=prefix,
    )

    resolved = _resolve([f"@{manifest}"], **kwargs)
    assert resolved.labels == ("label_one", f"{prefix}_002", "label_three")
    overridden = _resolve(
        [f"@{manifest}"], ["CLI_1,CLI_2", "CLI_3"], **kwargs
    )
    assert overridden.labels == ("CLI_1", "CLI_2", "CLI_3")
    with pytest.raises(SystemExit, match=r"supplied 1 labels for 3 expanded"):
        _resolve([f"@{manifest}"], ["partial"], **kwargs)


def test_multiple_manifests_and_existing_comma_syntax_preserve_order(tmp_path: Path) -> None:
    paths = [_file(tmp_path / f"e{index}.gtf") for index in range(1, 6)]
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    first.write_text(f"{paths[1]}\n{paths[2]}\n")
    second.write_text(f"{paths[3]}\n")
    resolved = _resolve(
        [f"{paths[0]},@{first}", f"@{second},{paths[4]}"],
        family="transcript evidence",
        input_option="--evidence-transcripts",
        label_option="--evidence-transcript-source-name",
        comma=True,
        prefix="transcript",
    )
    assert resolved.paths == paths
    assert resolved.labels is None
    assert split_input_tokens([f"{paths[0]},{paths[1]}", str(paths[2])], comma_separated=True) == [
        str(paths[0]), str(paths[1]), str(paths[2])
    ]


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("\tlabel\n", r":1: empty input path"),
        ("file\tlabel\textra\n", r":1: expected one or two tab-separated columns"),
        ("@nested.tsv\n", r":1: nested input manifests are not supported"),
        ("missing.gff3\n", r":1: referenced input does not exist"),
    ],
)
def test_manifest_row_errors_report_file_and_line(
    tmp_path: Path, contents: str, message: str
) -> None:
    manifest = tmp_path / "bad.tsv"
    manifest.write_text(contents)
    with pytest.raises(SystemExit, match=message):
        _resolve([f"@{manifest}"])


def test_missing_manifest_and_duplicate_resolved_paths_are_rejected(tmp_path: Path) -> None:
    missing = tmp_path / "missing.tsv"
    with pytest.raises(SystemExit, match=r"Input manifest does not exist"):
        _resolve([f"@{missing}"])

    item = _file(tmp_path / "inputs" / "same.gff3")
    manifest = tmp_path / "duplicates.tsv"
    manifest.write_text("inputs/../inputs/same.gff3\n")
    with pytest.raises(SystemExit, match=r"Duplicate model input.*first: direct.*manifest:"):
        _resolve([item, f"@{manifest}"])


def test_direct_only_inputs_retain_existing_paths_and_no_synthetic_labels(tmp_path: Path) -> None:
    paths = [_file(tmp_path / "a.gff3"), _file(tmp_path / "b.gff3")]
    resolved = _resolve(paths)
    assert resolved.paths == paths
    assert resolved.labels is None


def test_cli_manifest_and_direct_specifications_produce_equivalent_outputs(
    tmp_path: Path,
) -> None:
    model_template = (
        "##gff-version 3\n"
        "chr1\t{source}\tgene\t81\t290\t.\t+\t.\tID=g{source}\n"
        "chr1\t{source}\tmRNA\t81\t290\t.\t+\t.\tID=t{source};Parent=g{source}\n"
        "chr1\t{source}\texon\t81\t150\t.\t+\t.\tParent=t{source}\n"
        "chr1\t{source}\texon\t201\t290\t.\t+\t.\tParent=t{source}\n"
        "chr1\t{source}\tCDS\t101\t150\t.\t+\t1\tParent=t{source}\n"
        "chr1\t{source}\tCDS\t201\t250\t.\t+\t2\tParent=t{source}\n"
    )
    model_a = tmp_path / "a.gff3"
    model_b = tmp_path / "b.gff3"
    model_a.write_text(model_template.format(source="A"))
    model_b.write_text(model_template.format(source="B"))
    protein = tmp_path / "protein.gff3"
    protein.write_text(
        "##gff-version 3\n"
        "chr1\tP\tCDS\t101\t150\t.\t+\t1\tParent=p1\n"
        "chr1\tP\tCDS\t201\t250\t.\t+\t2\tParent=p1\n"
    )
    transcript = tmp_path / "transcript.gtf"
    transcript.write_text(
        'chr1\tT\texon\t81\t150\t.\t+\t.\tgene_id "tg"; transcript_id "tt";\n'
        'chr1\tT\texon\t201\t290\t.\t+\t.\tgene_id "tg"; transcript_id "tt";\n'
    )
    manifests = tmp_path / "manifests"
    manifests.mkdir()
    (manifests / "models.tsv").write_text("../a.gff3\tA\n../b.gff3\tB\n")
    (manifests / "proteins.tsv").write_text("../protein.gff3\tProteinSource\n")
    (manifests / "transcripts.tsv").write_text("../transcript.gtf\tTranscriptSource\n")
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")

    common = [
        "--pruning-level", "none",
        "--output-profile", "standard",
        "--write-output", "protein_evidence_model_support_tsv,transcript_evidence_model_support_tsv",
        "--progress-log-level", "basic",
        "--add-utrs",
    ]

    def run(prefix: Path, inputs: list[str]) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, "-m", "moselecta", *inputs, "--out-prefix", str(prefix), *common],
            cwd=tmp_path / "elsewhere",
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )

    (tmp_path / "elsewhere").mkdir()
    direct = run(tmp_path / "direct", [
        "--inputs", str(model_a), str(model_b),
        "--source-name", "A", "B",
        "--evidence-proteins", str(protein),
        "--evidence-protein-source-name", "ProteinSource",
        "--evidence-transcripts", str(transcript),
        "--evidence-transcript-source-name", "TranscriptSource",
    ])
    manifest = run(tmp_path / "manifest", [
        "--inputs", f"@{manifests / 'models.tsv'}",
        "--evidence-proteins", f"@{manifests / 'proteins.tsv'}",
        "--evidence-transcripts", f"@{manifests / 'transcripts.tsv'}",
    ])
    for suffix in (
        "final_annotation.gff3",
        "final_annotation.with_utrs.gff3",
        "utr_augmentation.audit.tsv",
        "protein_evidence_model_support.tsv",
        "transcript_evidence_model_support.tsv",
    ):
        assert (tmp_path / f"direct.{suffix}").read_bytes() == (
            tmp_path / f"manifest.{suffix}"
        ).read_bytes(), suffix
    assert "Expanded @" not in direct.stderr
    assert "2 model inputs" in manifest.stderr
    assert "1 protein evidence inputs" in manifest.stderr
    assert "1 transcript evidence inputs" in manifest.stderr

    mixed_manifest = manifests / "one_model.tsv"
    mixed_manifest.write_text("../b.gff3\tignored\n")
    mixed = run(tmp_path / "mixed", [
        "--inputs", str(model_a), f"@{mixed_manifest}",
        "--source-name", "A", "B",
        "--evidence-proteins", str(protein),
        "--evidence-protein-source-name", "ProteinSource",
        "--evidence-transcripts", str(transcript),
        "--evidence-transcript-source-name", "TranscriptSource",
    ])
    assert mixed.returncode == 0
    assert (tmp_path / "mixed.final_annotation.with_utrs.gff3").read_bytes() == (
        tmp_path / "direct.final_annotation.with_utrs.gff3"
    ).read_bytes()
