from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
INPUTS = [
    ROOT / "tests" / "data" / "alternative_splice_smoke_a.gff3",
    ROOT / "tests" / "data" / "alternative_splice_smoke_b.gff3",
]


def _run(prefix: Path, *extra: str, profile: str = "minimal") -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "moselecta",
            "--inputs",
            *(str(path) for path in INPUTS),
            "--source-name",
            "A",
            "B",
            "--out-prefix",
            str(prefix),
            "--pruning-level",
            "auto",
            "--output-profile",
            profile,
            "--progress-log-level",
            "none",
            *extra,
        ],
        cwd=ROOT,
        env=env,
        check=True,
        text=True,
        capture_output=True,
    )


def _suffixes(prefix: Path) -> set[str]:
    return {
        path.name.removeprefix(prefix.name + ".")
        for path in prefix.parent.glob(prefix.name + ".*")
    }


def test_disabling_strict_r1_effectively_disables_default_r1me(
    tmp_path: Path,
) -> None:
    prefix = tmp_path / "no_r1"
    _run(prefix, "--no-strict-r1-backfill-replacement")
    configuration = json.loads(
        (tmp_path / "no_r1.run_configuration.json").read_text()
    )
    selection = configuration["selection"]
    assert selection["strict_r1_backfill_replacement"] is False
    assert selection["strict_r1_backfill_replacement_configured"] is False
    assert selection["strict_r1_backfill_replacement_explicit"] is True
    assert selection["strict_r1_allow_multiexon_blockers"] is False
    assert selection["strict_r1_allow_multiexon_blockers_configured"] is True
    assert selection["strict_r1_allow_multiexon_blockers_explicit"] is False


def _ids(path: Path) -> set[str]:
    with path.open(newline="") as handle:
        return {
            row["output_transcript_id"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


@pytest.mark.parametrize("profile", ("minimal", "standard", "full", "debug"))
def test_profile_file_semantics_and_run_configuration(
    tmp_path: Path, profile: str
) -> None:
    prefix = tmp_path / profile
    _run(prefix, profile=profile)
    suffixes = _suffixes(prefix)
    for suffix in (
        "terminal_stop_canonicalisation.summary.tsv",
        "cds_5prime_boundary_canonicalisation.summary.tsv",
        "cds_phase_recovery.summary.tsv",
        "definite_internal_stop_exclusion.summary.tsv",
        "extreme_short_cds_exclusion.summary.tsv",
        "cds_terminus_refinement.summary.tsv",
        "alternative_splice_run_summary.tsv",
        "final_annotation.gff3",
        "run_configuration.json",
    ):
        assert suffix in suffixes

    action_suffixes = {
        "terminal_stop_canonicalisation.actions.tsv",
        "cds_5prime_boundary_canonicalisation.actions.tsv",
        "cds_phase_recovery.actions.tsv",
        "definite_internal_stop_exclusions.tsv",
        "extreme_short_cds_exclusions.tsv",
    }
    audit_suffixes = {
        "terminal_stop_canonicalisation.audit.tsv",
        "cds_5prime_boundary_canonicalisation.audit.tsv",
        "cds_phase_recovery.audit.tsv",
    }
    assert action_suffixes.issubset(suffixes) is (profile != "minimal")
    assert audit_suffixes.issubset(suffixes) is (profile in {"full", "debug"})
    assert ("cds_terminus_refinement.audit.tsv" in suffixes) is (
        profile != "minimal"
    )
    assert ("standardised_primary_models.gff3" in suffixes) is (
        profile != "minimal"
    )
    assert ("model_selection_scores.tsv" in suffixes) is (
        profile in {"full", "debug"}
    )
    assert ("structural_corroboration.tsv" in suffixes) is (
        profile in {"full", "debug"}
    )
    assert ("selected_models.tsv" in suffixes) is (profile != "minimal")
    assert "selected_models.gff3" not in suffixes or profile in {"full", "debug"}

    configuration = json.loads((tmp_path / f"{profile}.run_configuration.json").read_text())
    assert configuration["moselecta"]["version"] == "0.3.0"
    assert configuration["output"]["profile"] == profile
    assert configuration["pruning"]["requested_level"] == "auto"
    assert configuration["definite_internal_stop_exclusion"] == {
        "requested_enabled": True,
        "effective_enabled": False,
        "evidence_rescue": False,
        "execution_status": "prerequisite_unavailable:no_genome_fasta",
        "resolved_frame_policy": "exclude_if_active_frame_has_literal_internal_stop",
        "unresolved_frame_policy": "exclude_if_all_three_frames_have_literal_internal_stop",
    }
    assert configuration["extreme_short_cds_exclusion"] == {
        "default_enabled": True,
        "effective_enabled": True,
        "evidence_rescue": False,
        "included_terminal_stop_detection": "unavailable_without_genome_fasta",
        "included_terminal_stop_subtraction_bp": 3,
        "minimum_coding_length_bp": 30,
        "no_genome_policy": "subtract_nothing_unless_included_stop_established",
    }
    assert configuration["pruning"]["effective_level"] in {
        "medium", "high", "very-high"
    }
    assert configuration["output"]["enabled_output_keys"] == sorted(
        configuration["output"]["enabled_output_keys"]
    )


def test_profile_choice_does_not_change_final_biology_and_pruning_partition(
    tmp_path: Path,
) -> None:
    prefixes = {name: tmp_path / name for name in ("minimal", "standard", "full")}
    for profile, prefix in prefixes.items():
        _run(prefix, profile=profile)
    final_bytes = {
        (tmp_path / f"{profile}.final_annotation.gff3").read_bytes()
        for profile in prefixes
    }
    assert len(final_bytes) == 1

    retained = _ids(tmp_path / "standard.final_selected_models.pruned.tsv")
    removed = _ids(tmp_path / "standard.final_selected_models.pruned.removed.tsv")
    pruning_input = _ids(tmp_path / "full.final_selected_models.tsv")
    assert retained.isdisjoint(removed)
    assert retained | removed == pruning_input


@pytest.mark.parametrize("profile", ("minimal", "standard"))
def test_stop_after_scored_candidates_forces_principal_products(
    tmp_path: Path, profile: str
) -> None:
    prefix = tmp_path / profile
    _run(prefix, "--stop-after", "scored-candidates", profile=profile)
    suffixes = _suffixes(prefix)
    assert {
        "coding_models_with_superloci.gff3",
        "model_selection_scores.tsv",
        "model_structural_support.tsv",
        "model_selection_score_summary.tsv",
    }.issubset(suffixes)
    assert "final_annotation.gff3" not in suffixes
    assert "coding_models_with_superloci.no_singleton_subloci.gff3" not in suffixes


def test_explicit_output_precedence_and_skip_all_are_absolute(tmp_path: Path) -> None:
    added = tmp_path / "added"
    _run(
        added,
        "--write-output",
        "model_selection_scores_tsv",
        profile="minimal",
    )
    assert (tmp_path / "added.model_selection_scores.tsv").exists()

    corroboration = tmp_path / "corroboration"
    _run(
        corroboration,
        "--write-output",
        "structural_corroboration_tsv",
        profile="standard",
    )
    report = tmp_path / "corroboration.structural_corroboration.tsv"
    assert report.exists()
    with report.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert {
        "candidate_model_id",
        "candidate_source",
        "run_model_source_count",
        "sublocus_model_count",
        "corroborating_other_source_count",
        "corroborating_other_sources",
        "exact_cds_exon_chain_support_source_count",
        "structural_corroboration_linear_score",
        "structural_corroboration_log_score",
        "structural_corroboration_saturating_alpha_8_score",
    }.issubset(rows[0])
    assert {row["run_model_source_count"] for row in rows} == {"2"}

    skipped = tmp_path / "skipped"
    _run(
        skipped,
        "--write-output",
        "model_selection_scores_tsv",
        "--skip-output",
        "model_selection_scores_tsv",
        profile="full",
    )
    assert not (tmp_path / "skipped.model_selection_scores.tsv").exists()

    feature_skip = tmp_path / "feature_skip"
    _run(
        feature_skip,
        "--add-utrs",
        "--skip-output",
        "utr_augmentation_audit_tsv",
        profile="standard",
    )
    assert (tmp_path / "feature_skip.final_annotation.with_utrs.gff3").exists()
    assert not (tmp_path / "feature_skip.utr_augmentation.audit.tsv").exists()

    stage_skip = tmp_path / "stage_skip"
    _run(
        stage_skip,
        "--stop-after",
        "scored-candidates",
        "--skip-output",
        "combined_gff3",
        profile="minimal",
    )
    assert not (tmp_path / "stage_skip.coding_models_with_superloci.gff3").exists()
    assert (tmp_path / "stage_skip.model_selection_scores.tsv").exists()

    nothing = tmp_path / "nothing"
    _run(nothing, "--add-utrs", "--skip-output", "all", profile="full")
    assert _suffixes(nothing) == set()


def test_minimal_utr_and_explicit_refinement_audit_controls(tmp_path: Path) -> None:
    prefix = tmp_path / "minimal_utr"
    _run(prefix, "--add-utrs", profile="minimal")
    suffixes = _suffixes(prefix)
    assert "final_annotation.with_utrs.gff3" in suffixes
    assert "utr_augmentation.summary.tsv" in suffixes
    assert "utr_augmentation.audit.tsv" not in suffixes
    assert "cds_terminus_refinement.summary.tsv" in suffixes
    assert "cds_terminus_refinement.audit.tsv" not in suffixes

    explicit = tmp_path / "explicit"
    _run(
        explicit,
        "--add-utrs",
        "--write-output",
        "utr_augmentation_audit_tsv,cds_terminus_refinement_audit_tsv",
        profile="minimal",
    )
    assert (tmp_path / "explicit.utr_augmentation.audit.tsv").exists()
    assert (tmp_path / "explicit.cds_terminus_refinement.audit.tsv").exists()


def test_run_configuration_is_byte_deterministic_for_identical_invocation(
    tmp_path: Path,
) -> None:
    prefix = tmp_path / "repeat"
    _run(prefix, profile="standard")
    first = (tmp_path / "repeat.run_configuration.json").read_bytes()
    _run(prefix, profile="standard")
    second = (tmp_path / "repeat.run_configuration.json").read_bytes()
    assert first == second

