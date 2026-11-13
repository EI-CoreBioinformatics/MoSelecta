from __future__ import annotations

from pathlib import Path

from moselecta.models import CodingModel, Sublocus
from moselecta.selection import (
    final_prune_sequence_integrity_nonsingletons,
    final_prune_singleton_selected_models,
    sequence_integrity_pruning_fields,
)


def _model(model_id: str, source: str, start0: int = 99) -> CodingModel:
    return CodingModel(
        source_name=source,
        input_file=Path(f"{source}.gff3"),
        input_gene_id=f"g_{model_id}",
        input_transcript_id=model_id,
        output_gene_id=f"gene_{model_id}",
        output_transcript_id=model_id,
        chrom="chr1",
        strand="+",
        cds_parts=[(start0, start0 + 51), (start0 + 150, start0 + 301)],
        gene_record=None,
        transcript_record=None,
        child_records=[],
    )


def _prune_result(models: list[CodingModel], rows: list[dict[str, object]]) -> dict[str, object]:
    return {
        "initial_retained_models": list(models),
        "initial_retained_rows": list(rows),
        "retained_models": list(models),
        "retained_rows": list(rows),
        "removed_rows": [],
        "retained_selection_info_by_transcript": {
            str(row["output_transcript_id"]): dict(row) for row in rows
        },
    }


def test_sequence_integrity_concern_classes_are_whitelisted() -> None:
    fields = sequence_integrity_pruning_fields(
        {
            "sequence_integrity_flags": (
                "missing_start_codon,noncanonical_splice,ambiguous_sequence"
            )
        }
    )
    assert fields["sequence_integrity_pruning_concern"] == "true"
    assert fields["sequence_integrity_pruning_concern_class"] == "standard"
    assert fields["sequence_integrity_pruning_reasons"] == "missing_start_codon"
    assert fields["sequence_integrity_severe_pruning_concern"] == "false"

    severe = sequence_integrity_pruning_fields(
        {"sequence_integrity_flags": "missing_start_codon,missing_stop_codon"}
    )
    assert severe["sequence_integrity_pruning_concern_class"] == "severe"
    assert severe["sequence_integrity_severe_pruning_concern"] == "true"
    assert severe["sequence_integrity_severe_pruning_reasons"] == "missing_start_and_stop_codon"


def test_sequence_integrity_singleton_is_pruned_and_can_be_transcript_rescued() -> None:
    model = _model("tx1", "Galba")
    sub = Sublocus("sub1", "sl1", "chr1", "+", model.cds_start0, model.cds_end0, [model])
    row = {"output_transcript_id": "tx1", "source_name": "Galba", "final_prune_classes": "."}
    structural = {
        "tx1": {
            "intrinsic_penalty_components": ".",
            "intrinsic_penalty_score": 0.0,
            "sequence_integrity_flags": "missing_start_codon",
        }
    }
    removed = final_prune_singleton_selected_models(
        _prune_result([model], [row]),
        [sub],
        structural,
        prune_level="sequence_integrity",
        prune_sources=set(),
        pruning_level_name="mid",
    )
    assert removed["retained_models"] == []
    assert len(removed["removed_rows"]) == 1
    removed_row = removed["removed_rows"][0]
    assert removed_row["final_sequence_integrity_singleton_candidate"] == "true"
    assert removed_row["final_sequence_integrity_singleton_pruned"] == "true"
    assert "sequence_integrity_singleton" in str(removed_row["final_prune_classes"])

    transcript_rows = {
        "tx1": {
            "transcript_evidence_score_scaled": 0.90,
            "transcript_supporting_source_count": 2,
            "transcript_scoring_eligible_count": 2,
            "transcript_stranded_scoring_eligible_count": 2,
            "transcript_unstranded_scoring_eligible_count": 0,
            "transcript_cdna_context_evaluable": "true",
            "transcript_cdna_context_cds_to_min_cdna_fraction": 0.50,
        }
    }
    rescued = final_prune_singleton_selected_models(
        _prune_result([model], [row]),
        [sub],
        structural,
        prune_level="intrinsic_or_sequence_integrity",
        prune_sources=set(),
        pruning_level_name="mid",
        rescue_transcript_supported_singletons=True,
        transcript_evidence_rows_by_transcript=transcript_rows,
    )
    assert rescued["retained_models"] == [model]
    rescued_row = rescued["retained_rows"][0]
    assert rescued_row["final_sequence_integrity_singleton_candidate"] == "true"
    assert rescued_row["final_sequence_integrity_singleton_rescued"] == "true"
    assert rescued_row["final_transcript_singleton_rescue_applied"] == "true"
    assert rescued_row["final_transcript_singleton_rescue_min_cdna_context_fraction"] == 0.30


def test_sequence_integrity_nonsingleton_backfill_scope_uses_standard_concerns() -> None:
    primary = _model("primary", "Annevo")
    backfill = _model("backfill", "Galba", 500)
    sub = Sublocus("sub1", "sl1", "chr1", "+", 99, 801, [primary, backfill])
    rows = [
        {
            "output_transcript_id": "primary",
            "source_name": "Annevo",
            "final_selection_stage": "monosublocus",
            "sequence_integrity_flags": ".",
        },
        {
            "output_transcript_id": "backfill",
            "source_name": "Galba",
            "final_selection_stage": "backfill",
            "backfill_selected": "true",
            "sequence_integrity_flags": "missing_start_codon",
        },
    ]
    structural = {
        "primary": {"sequence_integrity_flags": "."},
        "backfill": {"sequence_integrity_flags": "missing_start_codon"},
    }
    out = final_prune_sequence_integrity_nonsingletons(
        _prune_result([primary, backfill], rows),
        [sub],
        structural,
        enabled=True,
        scope="backfill",
        pruning_level_name="mid",
    )
    assert out["retained_models"] == [primary]
    assert len(out["removed_rows"]) == 1
    removed = out["removed_rows"][0]
    assert removed["output_transcript_id"] == "backfill"
    assert removed["final_sequence_integrity_nonsingleton_pruned"] == "true"
    assert removed["final_prune_classes"] == "sequence_integrity_nonsingleton"


def test_sequence_integrity_nonsingleton_transcript_protection_uses_three_state_context() -> None:
    primary = _model("primary_context", "Annevo")
    backfill = _model("backfill_context", "Galba", 500)
    sub = Sublocus("sub_context", "sl_context", "chr1", "+", 99, 801, [primary, backfill])
    rows = [
        {
            "output_transcript_id": primary.output_transcript_id,
            "source_name": primary.source_name,
            "final_selection_stage": "monosublocus",
            "sequence_integrity_flags": ".",
        },
        {
            "output_transcript_id": backfill.output_transcript_id,
            "source_name": backfill.source_name,
            "final_selection_stage": "backfill",
            "backfill_selected": "true",
            "sequence_integrity_flags": "missing_start_codon",
        },
    ]
    structural = {
        primary.output_transcript_id: {"sequence_integrity_flags": "."},
        backfill.output_transcript_id: {"sequence_integrity_flags": "missing_start_codon"},
    }
    common_transcript = {
        "transcript_evidence_score_scaled": 0.90,
        "transcript_supporting_source_count": 2,
        "transcript_scoring_eligible_count": 2,
        "transcript_stranded_scoring_eligible_count": 2,
        "transcript_unstranded_scoring_eligible_count": 0,
        "transcript_cdna_context_pruning_concern": "false",
    }

    neutral = final_prune_sequence_integrity_nonsingletons(
        _prune_result([primary, backfill], rows),
        [sub],
        structural,
        enabled=True,
        scope="backfill",
        pruning_level_name="medium",
        protect_transcript_supported=True,
        transcript_evidence_rows_by_transcript={
            backfill.output_transcript_id: {
                **common_transcript,
                "transcript_cdna_context_evaluable": "false",
                "transcript_cdna_context_cds_to_min_cdna_fraction": ".",
            }
        },
    )
    assert neutral["retained_models"] == [primary, backfill]
    assert neutral["sequence_integrity_nonsingleton_protected_by_transcript"] == 1

    below = final_prune_sequence_integrity_nonsingletons(
        _prune_result([primary, backfill], rows),
        [sub],
        structural,
        enabled=True,
        scope="backfill",
        pruning_level_name="medium",
        protect_transcript_supported=True,
        transcript_evidence_rows_by_transcript={
            backfill.output_transcript_id: {
                **common_transcript,
                "transcript_cdna_context_evaluable": "true",
                "transcript_cdna_context_cds_to_min_cdna_fraction": 0.299,
            }
        },
    )
    assert below["retained_models"] == [primary]
    assert below["removed_rows"][0]["output_transcript_id"] == backfill.output_transcript_id


def test_sequence_integrity_nonsingleton_all_scope_keeps_backfill_standard_rule_and_adds_severe_global_rule() -> None:
    standard = _model("standard", "Galba")
    severe = _model("severe", "Helixer", 500)
    protected = _model("protected", "Tiberius", 900)
    backfill_standard = _model("backfill_standard", "Galba", 1300)
    sub = Sublocus(
        "sub1",
        "sl1",
        "chr1",
        "+",
        99,
        1601,
        [standard, severe, protected, backfill_standard],
    )
    rows = [
        {
            "output_transcript_id": "standard",
            "source_name": "Galba",
            "final_selection_stage": "monosublocus",
            "sequence_integrity_flags": "missing_start_codon",
        },
        {
            "output_transcript_id": "severe",
            "source_name": "Helixer",
            "final_selection_stage": "monosublocus",
            "sequence_integrity_flags": "internal_stop_codon",
        },
        {
            "output_transcript_id": "protected",
            "source_name": "Tiberius",
            "final_selection_stage": "monosublocus",
            "sequence_integrity_flags": "all_splices_noncanonical",
        },
        {
            "output_transcript_id": "backfill_standard",
            "source_name": "Galba",
            "final_selection_stage": "backfill",
            "backfill_selected": "true",
            "sequence_integrity_flags": "missing_start_codon",
        },
    ]
    structural = {
        str(row["output_transcript_id"]): {
            "sequence_integrity_flags": row["sequence_integrity_flags"]
        }
        for row in rows
    }
    protein_rows = {
        "protected": {
            "protein_evidence_score_scaled": 0.90,
            "protein_exact_cds_chain_source_count": 0,
        }
    }
    out = final_prune_sequence_integrity_nonsingletons(
        _prune_result([standard, severe, protected, backfill_standard], rows),
        [sub],
        structural,
        enabled=True,
        scope="all",
        pruning_level_name="mid",
        protect_protein_supported=True,
        protein_evidence_rows_by_transcript=protein_rows,
    )
    retained_ids = {model.output_transcript_id for model in out["retained_models"]}
    assert retained_ids == {"standard", "protected"}
    assert len(out["removed_rows"]) == 2
    removed_by_id = {
        row["output_transcript_id"]: row for row in out["removed_rows"]
    }
    assert removed_by_id["severe"]["final_sequence_integrity_nonsingleton_prune_route"] == "all_severe"
    assert removed_by_id["backfill_standard"]["final_sequence_integrity_nonsingleton_prune_route"] == "backfill_standard"
    assert out["sequence_integrity_nonsingleton_removed_backfill_standard_route"] == 1
    assert out["sequence_integrity_nonsingleton_removed_all_severe_route"] == 1
    protected_row = next(
        row for row in out["retained_rows"] if row["output_transcript_id"] == "protected"
    )
    assert protected_row["final_sequence_integrity_nonsingleton_protein_protected"] == "true"
