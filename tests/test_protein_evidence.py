from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def package_root() -> Path:
    return Path(__file__).resolve().parents[1]


def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    root = package_root()
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    return subprocess.run(args, cwd=root, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)


def test_protein_evidence_outputs_and_three_prime_wobble(tmp_path: Path) -> None:
    a = tmp_path / "a.gff3"
    b = tmp_path / "b.gff3"
    p = tmp_path / "protein.gff3"
    a.write_text(
        "##gff-version 3\n"
        "chr1\tsrcA\tgene\t100\t400\t.\t+\t.\tID=geneA\n"
        "chr1\tsrcA\tmRNA\t100\t400\t.\t+\t.\tID=txA;Parent=geneA\n"
        "chr1\tsrcA\tCDS\t100\t150\t.\t+\t0\tID=cdsA1;Parent=txA\n"
        "chr1\tsrcA\tCDS\t250\t400\t.\t+\t0\tID=cdsA2;Parent=txA\n"
    )
    b.write_text(
        "##gff-version 3\n"
        "chr1\tsrcB\tgene\t100\t397\t.\t+\t.\tID=geneB\n"
        "chr1\tsrcB\tmRNA\t100\t397\t.\t+\t.\tID=txB;Parent=geneB\n"
        "chr1\tsrcB\tCDS\t100\t150\t.\t+\t0\tID=cdsB1;Parent=txB\n"
        "chr1\tsrcB\tCDS\t250\t397\t.\t+\t0\tID=cdsB2;Parent=txB\n"
    )
    p.write_text(
        "##gff-version 3\n"
        "chr1\tprot\tCDS\t100\t150\t.\t+\t0\tID=p1.cds1;Parent=p1\n"
        "chr1\tprot\tCDS\t250\t397\t.\t+\t0\tID=p1.cds2;Parent=p1\n"
    )
    out = tmp_path / "out"
    result = run_command([
        sys.executable,
        "-m",
        "moselecta",
        "--inputs", str(a), str(b),
        "--source-name", "A", "B",
        "--evidence-proteins", str(p),
        "--out-prefix", str(out),
        "--progress-log-level", "basic",
        "--output-profile", "standard",
        "--write-output", "protein_evidence_model_support_tsv,model_selection_scores_tsv",
    ])
    assert result.returncode == 0
    support = tmp_path / "out.protein_evidence_model_support.tsv"
    assert support.exists()
    text = support.read_text()
    assert "protein_exact_cds_chain_supported" in text
    assert "true" in text
    assert "+protein_cds_intron_f1:" in text
    assert "+protein_exact_cds_chain_source_fraction:" in text
    assert "diagnostic_v0.2.0_dev_source_collapsed_mean_best_f1_plus_exact_and_source_fraction" not in text
    assert "protein_evidence_flags" in text
    assert "exact_chain_supported" in text
    assert "strong_protein_support" in text
    summary = tmp_path / "out.protein_evidence_summary.tsv"
    assert summary.exists()
    summary_text = summary.read_text()
    assert "raw_alignment_count" in summary_text
    assert "collapsed_alignment_count" in summary_text
    assert "models_with_exact_chain_support" in summary_text
    assert "per_source_raw_alignment_summary" in summary_text
    assert "Protein evidence inputs:" in result.stderr
    assert "Protein evidence loci:" in result.stderr
    assert "Protein evidence model associations:" in result.stderr
    selection_scores = tmp_path / "out.model_selection_scores.tsv"
    assert selection_scores.exists()
    selection_text = selection_scores.read_text()
    assert "model_structural_support_score" in selection_text.splitlines()[0]
    assert "protein_evidence_score_weight" in selection_text.splitlines()[0]
    assert "protein_evidence_score_scaled" in selection_text.splitlines()[0]
    assert "protein_evidence_score_in_selection" in selection_text.splitlines()[0]
    assert "	3.000000	true	" in selection_text
    score_summary = tmp_path / "out.model_selection_score_summary.tsv"
    assert score_summary.exists()
    score_summary_text = score_summary.read_text()
    assert "effective_protein_evidence_weight" in score_summary_text
    assert "rank_changes" in score_summary_text
    final_pruned_header = (tmp_path / "out.final_selected_models.pruned.tsv").read_text().splitlines()[0].split("\t")
    assert "cds_location" in final_pruned_header
    assert "model_structural_support_score" in final_pruned_header
    assert "model_structural_support_score_scaled" in final_pruned_header
    assert "protein_evidence_score" in final_pruned_header
    assert "protein_evidence_score_scaled" in final_pruned_header
    final_pruned_text = (tmp_path / "out.final_selected_models.pruned.tsv").read_text()
    assert "chr1:100..397" in final_pruned_text or "chr1:100..400" in final_pruned_text
    # In v0.2.4-dev patch 2, standard keeps the model-level protein audit
    # but no longer writes protein evidence construction internals by default.
    collapsed = tmp_path / "out.protein_evidence_collapsed_alignments.tsv"
    assert not collapsed.exists()
    assert not (tmp_path / "out.protein_evidence_subloci.tsv").exists()
    assert not (tmp_path / "out.protein_evidence_raw_to_collapsed.tsv").exists()

    out_explicit = tmp_path / "out_explicit"
    result_explicit = run_command([
        sys.executable,
        "-m",
        "moselecta",
        "--inputs", str(a), str(b),
        "--source-name", "A", "B",
        "--evidence-proteins", str(p),
        "--out-prefix", str(out_explicit),
        "--progress-log-level", "none",
        "--output-profile", "minimal",
        "--write-output", "protein_evidence_collapsed_alignments_tsv,protein_evidence_subloci_tsv",
    ])
    assert result_explicit.returncode == 0
    collapsed_explicit = tmp_path / "out_explicit.protein_evidence_collapsed_alignments.tsv"
    assert collapsed_explicit.exists()
    collapsed_text = collapsed_explicit.read_text()
    assert "flags" in collapsed_text.splitlines()[0]
    assert "single_source_only" in collapsed_text
    assert "multi_exon" in collapsed_text
    subloci_text = (tmp_path / "out_explicit.protein_evidence_subloci.tsv").read_text()
    assert "flags" in subloci_text.splitlines()[0]
    assert "single_source_only" in subloci_text
    assert not (tmp_path / "out_explicit.protein_evidence_raw_to_collapsed.tsv").exists()



def test_protein_evidence_weight_enters_model_selection_scores(tmp_path: Path) -> None:
    a = tmp_path / "a.gff3"
    b = tmp_path / "b.gff3"
    p = tmp_path / "protein.gff3"
    a.write_text(
        "##gff-version 3\n"
        "chr1\tsrcA\tgene\t100\t400\t.\t+\t.\tID=geneA\n"
        "chr1\tsrcA\tmRNA\t100\t400\t.\t+\t.\tID=txA;Parent=geneA\n"
        "chr1\tsrcA\tCDS\t100\t150\t.\t+\t0\tID=cdsA1;Parent=txA\n"
        "chr1\tsrcA\tCDS\t250\t400\t.\t+\t0\tID=cdsA2;Parent=txA\n"
    )
    b.write_text(
        "##gff-version 3\n"
        "chr1\tsrcB\tgene\t100\t397\t.\t+\t.\tID=geneB\n"
        "chr1\tsrcB\tmRNA\t100\t397\t.\t+\t.\tID=txB;Parent=geneB\n"
        "chr1\tsrcB\tCDS\t100\t150\t.\t+\t0\tID=cdsB1;Parent=txB\n"
        "chr1\tsrcB\tCDS\t250\t397\t.\t+\t0\tID=cdsB2;Parent=txB\n"
    )
    p.write_text(
        "##gff-version 3\n"
        "chr1\tprot\tCDS\t100\t150\t.\t+\t0\tID=p1.cds1;Parent=p1\n"
        "chr1\tprot\tCDS\t250\t397\t.\t+\t0\tID=p1.cds2;Parent=p1\n"
    )
    out = tmp_path / "weighted"
    result = run_command([
        sys.executable,
        "-m",
        "moselecta",
        "--inputs", str(a), str(b),
        "--source-name", "A", "B",
        "--evidence-proteins", str(p),
        "--selection-weight-protein-evidence", "0.5",
        "--out-prefix", str(out),
        "--progress-log-level", "none",
        "--output-profile", "standard",
        "--write-output", "model_selection_scores_tsv",
    ])
    assert result.returncode == 0
    scores = (tmp_path / "weighted.model_selection_scores.tsv").read_text().splitlines()
    header = scores[0].split("\t")
    assert "model_structural_support_score" in header
    assert "model_selection_score" in header
    assert "protein_evidence_score" in header
    assert "rank_changed_by_evidence" in header
    assert any("+protein_evidence_score:" in line for line in scores[1:])
    assert any("\ttrue\t" in line or line.endswith("\ttrue") for line in scores[1:])


def test_protein_supported_singleton_rescue_from_high_pruning(tmp_path: Path) -> None:
    a = tmp_path / "single.gff3"
    p = tmp_path / "protein.gff3"
    a.write_text(
        "##gff-version 3\n"
        "chr1\tsrcA\tgene\t100\t400\t.\t+\t.\tID=geneA\n"
        "chr1\tsrcA\tmRNA\t100\t400\t.\t+\t.\tID=txA;Parent=geneA\n"
        "chr1\tsrcA\tCDS\t100\t150\t.\t+\t0\tID=cdsA1;Parent=txA\n"
        "chr1\tsrcA\tCDS\t250\t400\t.\t+\t0\tID=cdsA2;Parent=txA\n"
    )
    p.write_text(
        "##gff-version 3\n"
        "chr1\tprot\tCDS\t100\t150\t.\t+\t0\tID=p1.cds1;Parent=p1\n"
        "chr1\tprot\tCDS\t250\t400\t.\t+\t0\tID=p1.cds2;Parent=p1\n"
    )
    out = tmp_path / "rescue"
    result = run_command([
        sys.executable,
        "-m",
        "moselecta",
        "--inputs", str(a),
        "--source-name", "A",
        "--evidence-proteins", str(p),
        "--pruning-level", "high",
        "--out-prefix", str(out),
        "--progress-log-level", "none",
        "--output-profile", "standard",
    ])
    assert result.returncode == 0
    pruned = (tmp_path / "rescue.final_selected_models.pruned.tsv").read_text().splitlines()
    assert len(pruned) == 2
    header = pruned[0].split("\t")
    row = dict(zip(header, pruned[1].split("\t")))
    assert row["final_singleton_pruned"] == "false"
    assert row["final_singleton_prune_reason"] == "singleton_sublocus_rescued_by_protein_evidence"
    assert row["final_protein_singleton_rescued"] == "true"
    assert "protein_exact_cds_chain_support" in row["final_protein_singleton_rescue_reason"]
    removed = (tmp_path / "rescue.final_selected_models.pruned.removed.tsv").read_text().splitlines()
    assert len(removed) == 1


def test_low_score_singleton_pruning_is_optional_and_audited(tmp_path: Path) -> None:
    a = tmp_path / "single_low.gff3"
    a.write_text(
        "##gff-version 3\n"
        "chr1\tsrcA\tgene\t100\t400\t.\t+\t.\tID=geneA\n"
        "chr1\tsrcA\tmRNA\t100\t400\t.\t+\t.\tID=txA;Parent=geneA\n"
        "chr1\tsrcA\tCDS\t100\t150\t.\t+\t0\tID=cdsA1;Parent=txA\n"
        "chr1\tsrcA\tCDS\t250\t400\t.\t+\t0\tID=cdsA2;Parent=txA\n"
    )
    p = tmp_path / "protein_far.gff3"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tprot\tgene\t1000\t1200\t.\t+\t.\tID=pgene1\n"
        "chr1\tprot\tmRNA\t1000\t1200\t.\t+\t.\tID=ptx1;Parent=pgene1\n"
        "chr1\tprot\tCDS\t1000\t1050\t.\t+\t0\tID=pcds1;Parent=ptx1\n"
        "chr1\tprot\tCDS\t1100\t1200\t.\t+\t0\tID=pcds2;Parent=ptx1\n"
    )
    out = tmp_path / "low"
    result = run_command([
        sys.executable,
        "-m",
        "moselecta",
        "--inputs", str(a),
        "--source-name", "A",
        "--evidence-proteins", str(p),
        "--pruning-level", "none",
        "--final-prune-low-score-singletons",
        "--no-low-score-singleton-require-protein-evidence-adequacy",
        "--out-prefix", str(out),
        "--progress-log-level", "none",
        "--output-profile", "standard",
    ])
    assert result.returncode == 0
    assert "Evidence-aware low-score singleton pruning: enabled=true" in result.stdout
    retained = (tmp_path / "low.final_selected_models.pruned.tsv").read_text().splitlines()
    removed = (tmp_path / "low.final_selected_models.pruned.removed.tsv").read_text().splitlines()
    assert len(retained) == 1
    assert len(removed) == 2
    header = removed[0].split("\t")
    row = dict(zip(header, removed[1].split("\t")))
    assert row["final_low_score_singleton_pruned"] == "true"
    assert row["final_low_score_singleton_prune_enabled"] == "true"
    assert row["final_low_score_singleton_protein_adequacy_required"] == "false"
    assert row["final_prune_classes"] == "low_score_singleton"
    assert "model_structural_support_score_scaled" in row["final_low_score_singleton_prune_reason"]


def test_low_score_nonsingleton_uses_scores_not_singleton_rescue_flag() -> None:
    from moselecta.models import CodingModel, Sublocus
    from moselecta.selection import final_prune_low_score_nonsingletons

    m1 = CodingModel(
        source_name="A",
        input_file=Path("a.gff3"),
        input_gene_id="g1",
        input_transcript_id="tx1",
        output_gene_id="MoSelecta_gene1",
        output_transcript_id="tx1",
        chrom="chr1",
        strand="+",
        cds_parts=[(99, 150), (249, 400)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )
    m2 = CodingModel(
        source_name="B",
        input_file=Path("b.gff3"),
        input_gene_id="g2",
        input_transcript_id="tx2",
        output_gene_id="MoSelecta_gene2",
        output_transcript_id="tx2",
        chrom="chr1",
        strand="+",
        cds_parts=[(99, 150), (260, 400)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )
    sub = Sublocus(
        sublocus_id="sub001",
        superlocus_id="sl001",
        chrom="chr1",
        strand="+",
        start0=99,
        end0=400,
        models=[m1, m2],
    )
    low_row = {
        "model_structural_support_score_scaled": "0.000000",
        "protein_evidence_score_scaled": "0.000000",
        "model_selection_score_scaled": "0.000000",
        "final_protein_singleton_rescue_applied": "true",
    }
    prune_result = {
        "initial_retained_models": [m1, m2],
        "initial_retained_rows": [],
        "retained_models": [m1, m2],
        "retained_rows": [
            {"output_transcript_id": "tx1", "source_name": "A", **low_row},
            {
                "output_transcript_id": "tx2",
                "source_name": "B",
                **low_row,
                "protein_evidence_score_scaled": "0.900000",
                "final_protein_singleton_rescue_applied": "false",
            },
        ],
        "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }
    out = final_prune_low_score_nonsingletons(
        prune_result,
        [sub],
        {},
        enabled=True,
        scope="all",
        protein_evidence_supplied=True,
        transcript_evidence_supplied=False,
        protein_adequacy_required=False,
        transcript_adequacy_required=True,
        protein_overlap_model_fraction=0.0,
        protein_supported_model_fraction=0.0,
        transcript_scoring_eligible_model_fraction=0.0,
        transcript_supported_model_fraction=0.0,
        min_protein_overlap_model_fraction=0.75,
        min_protein_supported_model_fraction=0.50,
        min_transcript_scoring_eligible_model_fraction=0.75,
        min_transcript_supported_model_fraction=0.50,
        max_model_structural_support_score_scaled=0.10,
        max_protein_evidence_score_scaled=0.10,
        max_transcript_evidence_score_scaled=0.10,
        max_model_selection_score_scaled=0.10,
        pruning_level_name="none",
    )
    assert out["retained_models"] == [m2]
    assert len(out["removed_rows"]) == 1
    assert out["low_score_nonsingleton_removed_count"] == 1
    assert out["low_score_nonsingleton_candidates"] == 2
    row = out["removed_rows"][0]
    assert row["final_low_score_nonsingleton_pruned"] == "true"
    assert row["final_low_score_nonsingleton_prune_enabled"] == "true"
    assert row["final_low_score_nonsingleton_prune_scope"] == "all"
    assert row["final_prune_classes"] == "low_score_nonsingleton"
    assert row["final_low_score_nonsingleton_policy"] == "evidence_aware"
    assert row["final_low_score_nonsingleton_prune_route"] == "evidence_aware_score_thresholds"
    assert "model_structural_support_score_scaled" in row["final_low_score_nonsingleton_prune_reason"]
    retained = out["retained_rows"][0]
    assert retained["output_transcript_id"] == "tx2"
    assert retained["final_low_score_nonsingleton_prune_reason"] == "score_and_context_thresholds_not_met"
    assert "low_score_nonsingleton_skipped_protein_rescue" not in out



def test_low_score_nonsingleton_backfill_scope_skips_nonbackfill_models() -> None:
    from moselecta.models import CodingModel, Sublocus
    from moselecta.selection import final_prune_low_score_nonsingletons

    m1 = CodingModel(
        source_name="A",
        input_file=Path("a.gff3"),
        input_gene_id="g1",
        input_transcript_id="tx1",
        output_gene_id="MoSelecta_gene1",
        output_transcript_id="tx1",
        chrom="chr1",
        strand="+",
        cds_parts=[(99, 150), (249, 400)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )
    m2 = CodingModel(
        source_name="B",
        input_file=Path("b.gff3"),
        input_gene_id="g2",
        input_transcript_id="tx2",
        output_gene_id="MoSelecta_gene2",
        output_transcript_id="tx2",
        chrom="chr1",
        strand="+",
        cds_parts=[(99, 150), (260, 400)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )
    sub = Sublocus(
        sublocus_id="sub001",
        superlocus_id="sl001",
        chrom="chr1",
        strand="+",
        start0=99,
        end0=400,
        models=[m1, m2],
    )
    base_row = {
        "model_structural_support_score_scaled": "0.000000",
        "protein_evidence_score_scaled": "0.000000",
        "model_selection_score_scaled": "0.000000",
        "final_protein_singleton_rescue_applied": "false",
    }
    prune_result = {
        "initial_retained_models": [m1, m2],
        "initial_retained_rows": [],
        "retained_models": [m1, m2],
        "retained_rows": [
            {"output_transcript_id": "tx1", "source_name": "A", "final_selection_stage": "monosublocus", "backfill_selected": ".", **base_row},
            {"output_transcript_id": "tx2", "source_name": "B", "final_selection_stage": "backfill", "backfill_selected": "true", **base_row},
        ],
        "removed_rows": [],
        "retained_selection_info_by_transcript": {},
    }
    out = final_prune_low_score_nonsingletons(
        prune_result,
        [sub],
        {},
        enabled=True,
        scope="backfill",
        protein_evidence_supplied=True,
        transcript_evidence_supplied=False,
        protein_adequacy_required=False,
        transcript_adequacy_required=True,
        protein_overlap_model_fraction=0.0,
        protein_supported_model_fraction=0.0,
        transcript_scoring_eligible_model_fraction=0.0,
        transcript_supported_model_fraction=0.0,
        min_protein_overlap_model_fraction=0.75,
        min_protein_supported_model_fraction=0.50,
        min_transcript_scoring_eligible_model_fraction=0.75,
        min_transcript_supported_model_fraction=0.50,
        max_model_structural_support_score_scaled=0.20,
        max_protein_evidence_score_scaled=0.20,
        max_transcript_evidence_score_scaled=0.20,
        max_model_selection_score_scaled=0.20,
        pruning_level_name="mid",
    )
    assert len(out["retained_models"]) == 1
    assert len(out["removed_rows"]) == 1
    assert out["low_score_nonsingleton_candidates"] == 1
    assert out["low_score_nonsingleton_skipped_scope"] == 1
    assert out["low_score_nonsingleton_removed_from_backfill"] == 1
    assert out["low_score_nonsingleton_removed_from_nonbackfill"] == 0
    retained_row = out["retained_rows"][0]
    removed_row = out["removed_rows"][0]
    assert retained_row["output_transcript_id"] == "tx1"
    assert retained_row["final_low_score_nonsingleton_prune_reason"] == "outside_scope"
    assert removed_row["output_transcript_id"] == "tx2"
    assert removed_row["final_low_score_nonsingleton_prune_scope"] == "backfill"
