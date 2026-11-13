from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
from pathlib import Path

import moselecta


def package_root() -> Path:
    return Path(__file__).resolve().parents[1]


def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    root = package_root()
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    return subprocess.run(
        args,
        cwd=root,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )


def test_module_version_smoke() -> None:
    result = run_command([sys.executable, "-m", "moselecta", "--version"])
    assert result.stdout.strip() == f"MoSelecta {moselecta.__version__}"


def test_package_version_smoke() -> None:
    result = run_command(
        [
            sys.executable,
            "-c",
            "import moselecta; print(moselecta.__version__)",
        ]
    )
    assert result.stdout.strip() == moselecta.__version__


def test_missing_input_reports_concise_cli_error_without_traceback(tmp_path: Path) -> None:
    missing = tmp_path / "missing.gff3"
    root = package_root()
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "moselecta",
            "--inputs",
            str(missing),
            "--out-prefix",
            str(tmp_path / "out"),
        ],
        cwd=root,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert result.returncode == 2
    assert "MoSelecta input/output error:" in result.stderr
    # Match the file name, not the full path: Windows paths are printed with
    # doubled backslashes in OSError messages, so str(missing) would not match.
    assert missing.name in result.stderr
    assert "Traceback" not in result.stderr


def test_project_and_runtime_versions_match() -> None:
    pyproject = (package_root() / "pyproject.toml").read_text()
    match = re.search(r'(?m)^version = "([^"]+)"$', pyproject)
    assert match is not None
    assert match.group(1) == moselecta.__version__ == "0.3.0"


def test_minimal_two_source_run(tmp_path: Path) -> None:
    gff_a = tmp_path / "a.gff3"
    gff_b = tmp_path / "b.gff3"
    gff_a.write_text(
        "##gff-version 3\n"
        "chr1\tsrcA\tgene\t100\t300\t.\t+\t.\tID=geneA\n"
        "chr1\tsrcA\tmRNA\t100\t300\t.\t+\t.\tID=txA;Parent=geneA\n"
        "chr1\tsrcA\texon\t100\t150\t.\t+\t.\tParent=txA\n"
        "chr1\tsrcA\texon\t200\t300\t.\t+\t.\tParent=txA\n"
        "chr1\tsrcA\tCDS\t100\t150\t.\t+\t0\tID=cdsA1;Parent=txA\n"
        "chr1\tsrcA\tCDS\t200\t300\t.\t+\t0\tID=cdsA2;Parent=txA\n"
    )
    gff_b.write_text(
        "##gff-version 3\n"
        "chr1\tsrcB\tgene\t100\t300\t.\t+\t.\tID=geneB\n"
        "chr1\tsrcB\tmRNA\t100\t300\t.\t+\t.\tID=txB;Parent=geneB\n"
        "chr1\tsrcB\texon\t100\t150\t.\t+\t.\tParent=txB\n"
        "chr1\tsrcB\texon\t200\t300\t.\t+\t.\tParent=txB\n"
        "chr1\tsrcB\tCDS\t100\t150\t.\t+\t0\tID=cdsB1;Parent=txB\n"
        "chr1\tsrcB\tCDS\t200\t300\t.\t+\t0\tID=cdsB2;Parent=txB\n"
    )
    out_prefix = tmp_path / "out"
    result = run_command([
        sys.executable,
        "-m",
        "moselecta",
        "--inputs",
        str(gff_a),
        str(gff_b),
        "--source-name",
        "A",
        "B",
        "--out-prefix",
        str(out_prefix),
        "--progress-log-level",
        "none",
        "--output-profile",
        "minimal",
    ])
    assert result.returncode == 0
    assert "Protein evidence selection weight is >0" not in result.stdout
    assert "Transcript evidence selection weight is >0" not in result.stdout
    assert "stage\tduration\tcumulative\tcurrent_rss_mib\tpeak_rss_mib" in result.stdout
    assert not (tmp_path / "out.final_selected_models.tsv").exists()
    assert not (tmp_path / "out.final_selected_models.pruned.tsv").exists()
    assert (tmp_path / "out.pruning_scenario_summary.tsv").exists()
    assert (tmp_path / "out.auto_pruning_diagnostics.tsv").exists()
    assert not (tmp_path / "out.standardised_primary_models.gff3").exists()
    assert (tmp_path / "out.final_annotation.gff3").exists()
    assert (tmp_path / "out.run_configuration.json").exists()
    assert "Automatic pruning assessment" in result.stdout
    assert "effective=medium" in result.stdout

    report_prefix = tmp_path / "report_only"
    report_only = run_command([
        sys.executable,
        "-m",
        "moselecta",
        "--inputs",
        str(gff_a),
        str(gff_b),
        "--source-name",
        "A",
        "B",
        "--out-prefix",
        str(report_prefix),
        "--pruning-level",
        "none",
        "--auto-pruning-report-only",
        "--progress-log-level",
        "none",
        "--output-profile",
        "standard",
        "--write-output",
        "final_selected_models_tsv",
    ])
    assert "Pruning level: requested=none; effective=none" in report_only.stdout
    assert (tmp_path / "report_only.pruning_scenario_summary.tsv").exists()
    assert (tmp_path / "report_only.auto_pruning_diagnostics.tsv").exists()
    assert (
        len((tmp_path / "report_only.final_selected_models.tsv").read_text().splitlines())
        == len((tmp_path / "report_only.final_selected_models.pruned.tsv").read_text().splitlines())
    )


def test_retained_internal_stop_is_counted_and_always_warned(tmp_path: Path) -> None:
    gff = tmp_path / "internal_stop.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=gene1\n"
        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\tsrc\texon\t1\t12\t.\t+\t.\tParent=tx1\n"
        "chr1\tsrc\tCDS\t1\t12\t.\t+\t.\tID=cds1;Parent=tx1\n"
    )
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nATGTAAAAATAA\n")
    result = run_command([
        sys.executable,
        "-m",
        "moselecta",
        "--inputs",
        str(gff),
        "--source-name",
        "A",
        "--genome-fasta",
        str(fasta),
        "--out-prefix",
        str(tmp_path / "out"),
        "--pruning-level",
        "none",
        "--output-profile",
        "minimal",
        "--minimum-coding-length-bp",
        "1",
        "--progress-log-level",
        "none",
    ])
    assert "selected_pre_pruning=1/1; final_retained=1/1" in result.stdout
    assert "retained_via_evidence_override=0; override_enabled=false" in result.stdout
    assert "WARNING: Final annotation contains retained internal-stop CDS models (1)." in result.stderr
    assert "genetic code" in result.stderr
    assert "--allow-internal-stop-evidence-rescue" in result.stderr


def test_terminal_stop_is_canonicalised_before_pipeline_scoring(tmp_path: Path) -> None:
    gff = tmp_path / "adjacent_stop.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=gene1\n"
        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\tsrc\texon\t1\t12\t.\t+\t.\tParent=tx1\n"
        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "chr1\tsrc\tstop_codon\t10\t12\t.\t+\t0\tParent=tx1\n"
    )
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nATGCCCGGGTAA\n")
    out_prefix = tmp_path / "canonical"
    completed = run_command([
        sys.executable, "-m", "moselecta", "--inputs", str(gff),
        "--source-name", "A", "--genome-fasta", str(fasta),
        "--out-prefix", str(out_prefix), "--pruning-level", "none",
        "--minimum-coding-length-bp", "1",
        "--output-profile", "standard", "--progress-log-level", "none",
        "--write-output", "sequence_integrity_model_metrics_tsv",
    ])
    assert completed.returncode == 0

    with (tmp_path / "canonical.terminal_stop_canonicalisation.actions.tsv").open(
        newline=""
    ) as handle:
        audit = next(csv.DictReader(handle, delimiter="\t"))
    assert audit["status"] == "canonicalised"
    assert audit["offset"] == "+3"
    assert audit["stop_codon"] == "TAA"

    with (tmp_path / "canonical.sequence_integrity_model_metrics.tsv").open(
        newline=""
    ) as handle:
        integrity = next(csv.DictReader(handle, delimiter="\t"))
    assert integrity["sequence_stop_codon_status"] == "included_in_cds"

    final_gff = (tmp_path / "canonical.final_annotation.gff3").read_text()
    assert "\tCDS\t1\t12\t" in final_gff
    assert "\tstop_codon\t10\t12\t" in final_gff
    assert "terminal_stop_canonicalised=offset:%2B3" in final_gff

    with (tmp_path / "canonical.cds_terminus_refinement.summary.tsv").open(
        newline=""
    ) as handle:
        refinement = next(csv.DictReader(handle, delimiter="\t"))
    assert refinement["execution_status"] == "zero_incomplete_models"


def test_auto_scenarios_and_effective_pruning_share_internal_stop_override(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "internal_stop.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=gene1\n"
        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\tsrc\texon\t1\t12\t.\t+\t.\tParent=tx1\n"
        "chr1\tsrc\tCDS\t1\t12\t.\t+\t.\tID=cds1;Parent=tx1\n"
    )
    protein = tmp_path / "protein.gff3"
    protein.write_text(
        "##gff-version 3\n"
        "chr1\tprotein\tCDS\t1\t12\t.\t+\t0\tID=p1.cds;Parent=p1\n"
    )
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nATGTAAAAATAA\n")

    results: dict[bool, tuple[subprocess.CompletedProcess[str], dict[str, str], dict[str, str]]] = {}
    for allow_override in (False, True):
        out_prefix = tmp_path / ("allow" if allow_override else "block")
        command = [
            sys.executable,
            "-m",
            "moselecta",
            "--inputs",
            str(gff),
            "--source-name",
            "A",
            "--genome-fasta",
            str(fasta),
            "--evidence-proteins",
            str(protein),
            "--out-prefix",
            str(out_prefix),
            "--pruning-level",
            "auto",
            "--output-profile",
            "standard",
            "--minimum-coding-length-bp",
            "1",
            "--progress-log-level",
            "none",
        ]
        if allow_override:
            command.append("--allow-internal-stop-evidence-rescue")
        completed = run_command(command)
        with (tmp_path / f"{out_prefix.name}.pruning_scenario_summary.tsv").open(
            newline=""
        ) as handle:
            medium = next(
                row for row in csv.DictReader(handle, delimiter="\t")
                if row["level"] == "medium"
            )
        retained_path = tmp_path / f"{out_prefix.name}.final_selected_models.pruned.tsv"
        with retained_path.open(newline="") as handle:
            retained_rows = list(csv.DictReader(handle, delimiter="\t"))
        results[allow_override] = (
            completed,
            medium,
            retained_rows[0] if retained_rows else {},
        )

    blocked_run, blocked_medium, blocked_row = results[False]
    assert blocked_medium["retained_models"] == "0"
    assert blocked_row == {}
    assert "retained_via_evidence_override=0; override_enabled=false" in blocked_run.stdout

    allowed_run, allowed_medium, allowed_row = results[True]
    assert allowed_medium["retained_models"] == "1"
    assert allowed_row["final_internal_stop_evidence_criteria_met"] == "true"
    assert allowed_row["final_internal_stop_evidence_rescue_applied"] == "true"
    assert allowed_row["final_internal_stop_evidence_rescue_applied_route"] == "protein_evidence"
    assert "retained_via_evidence_override=1; override_enabled=true" in allowed_run.stdout


def test_compare_module_smoke() -> None:
    result = run_command([sys.executable, "-m", "moselecta.compare", "-h"])
    assert result.returncode == 0
    assert "Compare MoSelecta output files" in result.stdout


def test_expert_help_reports_current_nonsingleton_scopes() -> None:
    result = run_command([sys.executable, "-m", "moselecta", "--expert-help"])
    help_text = " ".join(result.stdout.split())
    assert "The revised medium/high/very-high presets use all" in help_text
    assert "Medium uses backfill; high, very-high, and explicit" in help_text
    assert "automatic preset pruning uses backfill" not in help_text
    assert "medium/high presets use backfill" not in help_text.lower()


def test_help_levels_are_help_only_and_cumulative() -> None:
    basic = run_command([sys.executable, "-m", "moselecta", "--help"])
    full = run_command([sys.executable, "-m", "moselecta", "--full-help"])
    expert = run_command([sys.executable, "-m", "moselecta", "--expert-help"])
    assert "--inputs" in basic.stdout
    assert "--sequence-penalty-missing-start" not in basic.stdout
    assert "--sequence-penalty-missing-start" in full.stdout
    assert "--strict-r1-backfill-replacement" not in full.stdout
    assert "--strict-r1-backfill-replacement" in expert.stdout
    assert len(basic.stdout) < len(full.stdout) < len(expert.stdout)


def test_auto_high_stdout_uses_effective_high_configuration(tmp_path: Path) -> None:
    inputs: list[str] = []
    for source in ("A", "B", "C"):
        gff = tmp_path / f"{source}.gff3"
        lines = ["##gff-version 3"]
        for index in range(100):
            start = 1_000 + index * 1_000
            end = start + 149
            gene_id = f"{source}_gene_{index:03d}"
            transcript_id = f"{source}_tx_{index:03d}"
            lines.extend([
                f"chr1\t{source}\tgene\t{start}\t{end}\t.\t+\t.\tID={gene_id}",
                f"chr1\t{source}\tmRNA\t{start}\t{end}\t.\t+\t.\tID={transcript_id};Parent={gene_id}",
                f"chr1\t{source}\texon\t{start}\t{end}\t.\t+\t.\tParent={transcript_id}",
                f"chr1\t{source}\tCDS\t{start}\t{end}\t.\t+\t0\tID={transcript_id}.cds;Parent={transcript_id}",
            ])
        gff.write_text("\n".join(lines) + "\n")
        inputs.append(str(gff))

    out_prefix = tmp_path / "auto_high"
    result = run_command([
        sys.executable,
        "-m",
        "moselecta",
        "--inputs",
        *inputs,
        "--source-name",
        "A",
        "B",
        "C",
        "--out-prefix",
        str(out_prefix),
        "--pruning-level",
        "auto",
        "--output-profile",
        "minimal",
        "--progress-log-level",
        "none",
    ])
    stdout = result.stdout
    assert "Pruning level: requested=auto; effective=high" in stdout
    assert "max_model_score_scaled=0.350" in stdout
    assert "max_model_selection_score_scaled=0.350" in stdout
    assert "Low-score nonsingleton pruning: enabled=true; scope=all; policy=model_confidence_fallback" in stdout
    assert "run_model_source_count=3" in stdout
    assert "nonsingleton_policy=model_confidence_fallback" in stdout
    assert "Sequence-integrity nonsingleton pruning: enabled=false; scope=all" in stdout
