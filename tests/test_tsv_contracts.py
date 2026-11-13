from __future__ import annotations

from pathlib import Path

import pytest

from moselecta.constants import FINAL_PRUNE_EXTRA_COLUMNS
from moselecta.selection import (
    selection_stage_output_columns,
    write_backfill_outputs,
    write_final_pruned_outputs,
)
from moselecta.tsvio import format_tsv_value, write_dict_tsv


@pytest.mark.parametrize(("value", "expected"), ((True, "true"), (False, "false")))
def test_tsv_boole_use_lowercase_convention(value: bool, expected: str) -> None:
    assert format_tsv_value(value) == expected


def test_write_dict_tsv_serializes_boolean_cells_lowercase(tmp_path: Path) -> None:
    path = tmp_path / "booleans.tsv"
    write_dict_tsv([{"yes": True, "no": False}], path, ["yes", "no"])
    assert path.read_text().splitlines() == ["yes\tno", "true\tfalse"]


def test_selection_outputs_share_one_header_contract(tmp_path: Path) -> None:
    expected = selection_stage_output_columns()
    assert len(expected) == len(set(expected))

    paths = {
        "backfill_candidates_tsv": tmp_path / "candidates.tsv",
        "backfilled_models_tsv": tmp_path / "backfilled.tsv",
        "backfill_rejected_models_tsv": tmp_path / "rejected.tsv",
        "final_selected_models_tsv": tmp_path / "selected.tsv",
    }
    write_backfill_outputs(
        {
            "backfill_candidate_rows": [],
            "backfilled_rows": [],
            "backfill_rejected_rows": [],
            "final_rows": [],
        },
        **paths,
    )
    for path in paths.values():
        assert path.read_text().rstrip("\n").split("\t") == expected


def test_final_prune_headers_extend_selection_contract_exactly(tmp_path: Path) -> None:
    retained = tmp_path / "retained.tsv"
    removed = tmp_path / "removed.tsv"
    expected = selection_stage_output_columns() + FINAL_PRUNE_EXTRA_COLUMNS
    assert len(expected) == len(set(expected))

    write_final_pruned_outputs(
        {"retained_rows": [], "removed_rows": []},
        final_pruned_selected_models_tsv=retained,
        final_pruned_selected_removed_tsv=removed,
    )
    assert retained.read_text().rstrip("\n").split("\t") == expected
    assert removed.read_text().rstrip("\n").split("\t") == expected


def test_group_c_audit_columns_use_descriptive_edge_terminology() -> None:
    assert "final_singleton_overlap_rescue_overlap_bp_geometry_candidate" in FINAL_PRUNE_EXTRA_COLUMNS
    assert "final_singleton_overlap_rescue_same_source_edge" in FINAL_PRUNE_EXTRA_COLUMNS
    assert "final_singleton_overlap_rescue_preexisting_geometry_candidate" not in FINAL_PRUNE_EXTRA_COLUMNS
    assert "final_singleton_overlap_rescue_same_source_rejected" not in FINAL_PRUNE_EXTRA_COLUMNS


def test_internal_stop_evidence_override_audit_columns_are_in_final_prune_contract() -> None:
    expected = {
        "final_model_has_internal_stop_codon",
        "final_allow_internal_stop_evidence_rescue",
        "final_internal_stop_evidence_criteria_met",
        "final_internal_stop_evidence_criteria_routes",
        "final_internal_stop_evidence_rescue_blocked",
        "final_internal_stop_evidence_rescue_block_reason",
        "final_internal_stop_evidence_rescue_applied",
        "final_internal_stop_evidence_rescue_applied_route",
    }
    assert expected <= set(FINAL_PRUNE_EXTRA_COLUMNS)
