from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from moselecta.alternative_splice import AlternativeSpliceReportPolicy
from moselecta.console_reporting import (
    format_counter,
    log_section,
    summarise_score_rows_by_source,
)
from moselecta.models import CodingModel, FeatureRecord
from moselecta.pipeline_stages import (
    AlternativeSpliceConfigurationError,
    AlternativeSpliceSelectionSettings,
    alternative_splice_evidence_requirement_available,
    build_standardised_primary_annotation,
    evaluate_alternative_splice_stage,
    resolve_alternative_splice_execution_plan,
)


def _record(start0: int, end0: int, source: str) -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1",
        source=source,
        feature_type="CDS",
        start=start0 + 1,
        end=end0,
        score=".",
        strand="+",
        phase="0",
        attrs={"Parent": "tx"},
        raw_attrs="Parent=tx",
        input_file=Path(f"{source}.gff3"),
        source_name=source,
        line_number=1,
        record_index=1,
        raw_line="",
    )


def _model(
    name: str,
    parts: list[tuple[int, int]],
    source: str,
) -> CodingModel:
    model = CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"input_gene_{name}",
        input_transcript_id=f"input_transcript_{name}",
        output_gene_id=f"internal_gene_{name}",
        output_transcript_id=name,
        chrom="chr1",
        strand="+",
        cds_parts=parts,
        gene_record=None,
        transcript_record=None,
        child_records=[_record(start0, end0, source) for start0, end0 in parts],
    )
    model.superlocus_id = "super1"
    model.sublocus_id = "sub1"
    return model


def _score_row(model: CodingModel, score: float) -> dict[str, object]:
    return {
        "output_transcript_id": model.output_transcript_id,
        "model_selection_score": score,
        "model_selection_score_scaled": score,
        "model_selection_order_rank": 1,
        "model_completeness_score": 1.0,
        "intrinsic_penalty_score": 0.0,
        "consensus_similarity_score": 1.0,
        "structural_support_score": 1.0,
        "source_prior_score": 1.0,
        "cds_length": model.cds_length,
        "model_structural_support_score": 1.0,
        "model_structural_support_score_scaled": 1.0,
        "pairwise_similarity_score": 1.0,
        "protein_evidence_status": "protein_evidence_not_configured",
        "protein_evidence_score_scaled": 0.0,
        "transcript_evidence_status": "transcript_evidence_not_configured",
        "transcript_evidence_score_scaled": 0.0,
        "sequence_integrity_status": "genome_not_configured",
        "cds_introns_total": len(model.cds_introns),
        "cds_introns_supported_by_other_source": 0,
        "cds_intron_support_source_min": 1 if model.cds_introns else ".",
        "cds_splice_sites_total": len(model.cds_splice_sites),
        "cds_splice_sites_supported_by_other_source": 0,
        "cds_splice_site_support_source_min": (
            1 if model.cds_splice_sites else "."
        ),
        "internal_cds_exons_total": len(model.internal_cds_exons),
        "internal_cds_exons_supported_by_other_source": 0,
        "internal_cds_exon_support_source_min": (
            1 if model.internal_cds_exons else "."
        ),
        "terminal_cds_exons_total": len(model.terminal_cds_exons),
        "terminal_cds_exons_supported_by_other_source": 0,
        "terminal_cds_exon_support_source_min": (
            1 if model.terminal_cds_exons else "."
        ),
        "cds_base_consensus_precision": 0.0,
        "intrinsic_penalty_cds_length_not_multiple_of_3_count": 0,
        "intrinsic_penalty_short_cds_count": 0,
        "intrinsic_penalty_short_cds_intron_count": 0,
        "intrinsic_penalty_short_internal_cds_exon_count": 0,
        "intrinsic_penalty_short_terminal_cds_exon_count": 0,
        "intrinsic_penalty_long_cds_intron_count": 0,
    }


def _report_policy() -> AlternativeSpliceReportPolicy:
    return AlternativeSpliceReportPolicy(
        terminal_cds_3prime_wobble_bp=3,
        protein_terminal_cds_3prime_wobble_bp=3,
        transcript_min_cdna_context_fraction=0.30,
        transcript_require_stranded_support=False,
        transcript_block_on_cdna_context_concern=True,
        transcript_cdna_context_neutral_fraction=0.30,
        transcript_cdna_context_severe_fraction=0.15,
    )


def _selection_settings() -> AlternativeSpliceSelectionSettings:
    return AlternativeSpliceSelectionSettings(
        preset="intermediate",
        support_policy="exact",
        evidence_requirement="none",
        max_alternatives_per_gene=50,
        min_model_support_sources=2,
        min_model_support_source_fraction=None,
        min_reciprocal_cds_overlap=None,
        min_frame_concordance=None,
        min_model_selection_score_scaled=None,
        min_primary_score_ratio=None,
    )


@pytest.mark.parametrize(
    ("no_alternatives", "report", "expected"),
    (
        (False, False, (True, True, False)),
        (True, False, (False, False, False)),
        (True, True, (False, True, True)),
    ),
)
def test_execution_plan_preserves_final_selector_states(
    no_alternatives: bool,
    report: bool,
    expected: tuple[bool, bool, bool],
) -> None:
    plan = resolve_alternative_splice_execution_plan(
        stop_after="final",
        no_alternative_splice=no_alternatives,
        alternative_splice_report=report,
        explicitly_requested_outputs=(),
        enabled_outputs=(),
        raw_args=(),
    )
    assert (plan.enabled, plan.analysis_requested, plan.tsv_requested) == expected


def test_execution_plan_preserves_stop_and_inactive_tuning_semantics() -> None:
    stopped = resolve_alternative_splice_execution_plan(
        stop_after="backfill",
        no_alternative_splice=False,
        alternative_splice_report=False,
        explicitly_requested_outputs=(),
        enabled_outputs=(),
        raw_args=(),
    )
    assert not stopped.analysis_requested

    primary_only = resolve_alternative_splice_execution_plan(
        stop_after="final",
        no_alternative_splice=True,
        alternative_splice_report=False,
        explicitly_requested_outputs=(),
        enabled_outputs=(),
        raw_args=("--alternative-splice-preset=stringent",),
    )
    assert primary_only.supplied_inactive_tuning_options == (
        "--alternative-splice-preset",
    )


@pytest.mark.parametrize(
    ("requirement", "protein", "transcript", "available"),
    (
        ("none", False, False, True),
        ("protein", True, False, True),
        ("protein", False, True, False),
        ("transcript", False, True, True),
        ("protein-or-transcript", True, False, True),
        ("protein-or-transcript", False, False, False),
        ("protein-and-transcript", True, True, True),
        ("protein-and-transcript", True, False, False),
    ),
)
def test_evidence_requirement_availability_is_unchanged(
    requirement: str,
    protein: bool,
    transcript: bool,
    available: bool,
) -> None:
    assert alternative_splice_evidence_requirement_available(
        requirement,
        protein_configured=protein,
        transcript_configured=transcript,
    ) is available


def test_frozen_primary_stage_prefers_pruned_models_and_stabilises_public_ids() -> None:
    first = _model("first", [(100, 160), (220, 280)], "A")
    second = _model("second", [(400, 460), (520, 580)], "B")
    rows = {
        first.output_transcript_id: _score_row(first, 1.0),
        second.output_transcript_id: _score_row(second, 0.9),
    }
    frozen = build_standardised_primary_annotation(
        final_selected_models=[first, second],
        final_pruned_models=[second],
        model_score_rows_by_transcript=rows,
    )
    assert frozen.models == (second,)
    assert frozen.annotation.genes[0].gene_id == "MoSelecta_000001"
    assert frozen.annotation.genes[0].primary_transcript.transcript_id.endswith(".1")


def test_alternative_stage_keeps_selection_and_emission_as_separate_states() -> None:
    primary = _model(
        "primary",
        [(0, 90), (180, 270), (360, 450)],
        "A",
    )
    alternatives = [
        _model(
            f"alternative_{source}",
            [(0, 90), (180, 240), (330, 450)],
            source,
        )
        for source in ("B", "C", "D")
    ]
    models = [primary, *alternatives]
    rows = {
        model.output_transcript_id: _score_row(
            model, 1.0 if model is primary else 0.8
        )
        for model in models
    }
    frozen = build_standardised_primary_annotation(
        final_selected_models=[primary],
        final_pruned_models=None,
        model_score_rows_by_transcript=rows,
    )
    common = {
        "frozen_primary": frozen,
        "all_scored_models": models,
        "model_score_rows_by_transcript": rows,
        "protein_rows_by_transcript": {},
        "protein_collapsed_alignments": [],
        "transcript_rows_by_transcript": {},
        "transcript_superloci": [],
        "report_policy": _report_policy(),
        "run_metadata": {},
        "selection_settings": _selection_settings(),
        "terminal_cds_3prime_wobble_bp": 3,
    }
    emitted = evaluate_alternative_splice_stage(
        **common,
        alternatives_emitted=True,
    )
    diagnostic = evaluate_alternative_splice_stage(
        **common,
        alternatives_emitted=False,
    )

    accepted_emitted = tuple(
        item.candidate_model.output_transcript_id
        for item in emitted.selection.accepted_alternatives
    )
    accepted_diagnostic = tuple(
        item.candidate_model.output_transcript_id
        for item in diagnostic.selection.accepted_alternatives
    )
    assert accepted_emitted == accepted_diagnostic == ("alternative_B",)
    assert emitted.final_annotation.transcript_count == 2
    assert diagnostic.final_annotation is frozen.annotation
    assert diagnostic.report["run_rows"][0][
        "alternative_splice_alternatives_emitted"
    ] is False


def test_selector_configuration_is_resolved_after_candidate_facts(
    monkeypatch,
) -> None:
    primary = _model("primary", [(0, 90), (180, 270)], "A")
    rows = {primary.output_transcript_id: _score_row(primary, 1.0)}
    frozen = build_standardised_primary_annotation(
        final_selected_models=[primary],
        final_pruned_models=None,
        model_score_rows_by_transcript=rows,
    )
    calls: list[str] = []

    def fake_analyse(**kwargs):
        calls.append("candidate_facts")
        return {"candidate_rows": [], "gene_rows": [], "run_rows": []}

    monkeypatch.setattr(
        "moselecta.pipeline_stages.analyze_alternative_splice_report",
        fake_analyse,
    )
    invalid_settings = _selection_settings()
    invalid_settings = AlternativeSpliceSelectionSettings(
        preset=invalid_settings.preset,
        support_policy=invalid_settings.support_policy,
        evidence_requirement=invalid_settings.evidence_requirement,
        max_alternatives_per_gene=-1,
        min_model_support_sources=invalid_settings.min_model_support_sources,
        min_model_support_source_fraction=None,
        min_reciprocal_cds_overlap=None,
        min_frame_concordance=None,
        min_model_selection_score_scaled=None,
        min_primary_score_ratio=None,
    )
    with pytest.raises(AlternativeSpliceConfigurationError):
        evaluate_alternative_splice_stage(
            frozen_primary=frozen,
            all_scored_models=[primary],
            model_score_rows_by_transcript=rows,
            protein_rows_by_transcript={},
            protein_collapsed_alignments=[],
            transcript_rows_by_transcript={},
            transcript_superloci=[],
            report_policy=_report_policy(),
            run_metadata={},
            selection_settings=invalid_settings,
            terminal_cds_3prime_wobble_bp=3,
            alternatives_emitted=True,
        )
    assert calls == ["candidate_facts"]


def test_console_reporting_helpers_preserve_text_and_summary(capsys) -> None:
    log_section("Example")
    assert capsys.readouterr().out == "\n=======\nExample\n=======\n"
    assert format_counter(Counter({"B": 1, "A": 2})) == (
        "A=2, B=1"
    )
    assert summarise_score_rows_by_source([
        {
            "source_name": "A",
            "model_selection_score_scaled": 0.5,
            "model_selection_score": -0.25,
        },
        {
            "source_name": "A",
            "model_selection_score_scaled": 1.0,
            "model_selection_score": 0.75,
        },
    ]) == {
        "A": {"count": 2.0, "scaled_sum": 1.5, "raw_sum": 0.5}
    }
