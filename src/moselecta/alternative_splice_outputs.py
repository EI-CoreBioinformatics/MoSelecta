"""Diagnostic GFF output for existing alternative-splice decisions.

This module never evaluates biological eligibility. It joins the immutable
selector result to alternative-splice fact rows and actual candidate models for
browser-oriented serialization.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

from .alternative_splice_selection import (
    AlternativeCandidateDecision,
    AlternativeSpliceSelectionResult,
)
from .annotation_outputs import (
    alternative_splice_configuration_comments,
    write_source_model_children,
)
from .constants import GFF3_SOURCE
from .final_annotation import FinalAnnotation
from .models import CodingModel
from .outputs import feature_line
from .scoring import fmt_optional_float


AUDIT_CANDIDATE_ID_PREFIX = "MoSelecta_alt_candidate_"


@dataclass(frozen=True, slots=True)
class AlternativeSpliceAuditCandidate:
    """One uniquely owned structural candidate prepared for deterministic output."""

    audit_id: str
    public_gene_id: str
    public_primary_transcript_id: str
    preset: str
    support_policy: str
    evidence_requirement: str
    model: CodingModel
    fact_row: Mapping[str, object]
    decision: AlternativeCandidateDecision


def _truth(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "passed"}


def _optional_float(value: object) -> float | None:
    if value in {None, ".", ""}:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _boolean_text(value: bool) -> str:
    return str(value).lower()


def _containment_value(
    candidate_in_primary: bool,
    primary_in_candidate: bool,
) -> str:
    if candidate_in_primary and primary_in_candidate:
        return "bidirectional"
    if candidate_in_primary:
        return "candidate_in_primary"
    if primary_in_candidate:
        return "primary_in_candidate"
    return "none"


def _concern_names(
    row: Mapping[str, object],
    names: Iterable[str],
) -> tuple[str, ...]:
    return tuple(name for name in names if _truth(row.get(name)))


def build_alternative_splice_audit_candidates(
    *,
    candidate_rows: Iterable[Mapping[str, object]],
    selection_result: AlternativeSpliceSelectionResult,
    all_scored_models: Iterable[CodingModel],
    primary_annotation: FinalAnnotation,
) -> tuple[AlternativeSpliceAuditCandidate, ...]:
    """Build stable audit identities for all uniquely owned structural candidates.

    IDs depend only on frozen public gene order, the candidate rank,
    and internal transcript ID.  Preset thresholds, support policy, decisions,
    and active alternative ordinals never participate in identity assignment.
    """
    model_by_internal_id = {
        model.output_transcript_id: model for model in all_scored_models
    }
    decision_by_assessment_id = selection_result.decisions_by_assessment_id()
    public_by_primary_internal_id = {
        internal_id: gene
        for internal_id, gene in (
            primary_annotation.gene_by_primary_internal_transcript_id().items()
        )
    }
    public_gene_order = {
        gene.gene_id: index for index, gene in enumerate(primary_annotation.genes)
    }

    prepared: list[
        tuple[int, int, str, str, str, Mapping[str, object], AlternativeCandidateDecision]
    ] = []
    seen_candidate_ids: set[str] = set()
    for original_row in candidate_rows:
        row = dict(original_row)
        if not _truth(row.get("structural_candidate")):
            continue
        if str(row.get("ownership_status", ".")) != "unique_owner":
            continue
        candidate_id = str(row.get("candidate_output_transcript_id", "."))
        primary_id = str(row.get("primary_output_transcript_id", "."))
        assessment_id = str(row.get("alternative_candidate_assessment_id", "."))
        if "." in {candidate_id, primary_id, assessment_id}:
            raise ValueError("Uniquely owned audit candidate has an incomplete identity")
        if candidate_id in seen_candidate_ids:
            raise ValueError(
                "Uniquely owned audit population contains duplicate candidate: "
                + candidate_id
            )
        seen_candidate_ids.add(candidate_id)
        model = model_by_internal_id.get(candidate_id)
        if model is None:
            raise ValueError("Audit candidate model is unavailable: " + candidate_id)
        gene = public_by_primary_internal_id.get(primary_id)
        if gene is None:
            raise ValueError("Audit candidate primary is not in FinalAnnotation: " + primary_id)
        decision = decision_by_assessment_id.get(assessment_id)
        if decision is None:
            raise ValueError("Audit candidate selector decision is unavailable: " + assessment_id)
        if (
            decision.candidate_internal_transcript_id != candidate_id
            or decision.primary_internal_transcript_id != primary_id
        ):
            raise ValueError("Audit candidate fact/decision identities disagree: " + assessment_id)
        rank_value = row.get("alternative_candidate_order_rank")
        try:
            rank = int(rank_value)
        except (TypeError, ValueError) as exc:
            raise ValueError("Audit candidate has no deterministic candidate rank: " + candidate_id) from exc
        prepared.append((
            public_gene_order[gene.gene_id],
            rank,
            candidate_id,
            gene.gene_id,
            gene.primary_transcript.transcript_id,
            row,
            decision,
        ))

    prepared.sort(key=lambda item: item[:3])
    audit_candidates = tuple(
        AlternativeSpliceAuditCandidate(
            audit_id=f"{AUDIT_CANDIDATE_ID_PREFIX}{index:09d}",
            public_gene_id=public_gene_id,
            public_primary_transcript_id=public_primary_id,
            preset=selection_result.config.requested_preset.value,
            support_policy=selection_result.config.support_policy.value,
            evidence_requirement=selection_result.config.evidence_requirement.value,
            model=model_by_internal_id[candidate_id],
            fact_row=row,
            decision=decision,
        )
        for index, (
            _gene_order,
            _rank,
            candidate_id,
            public_gene_id,
            public_primary_id,
            row,
            decision,
        ) in enumerate(prepared, start=1)
    )

    accepted_ids = {
        accepted.candidate_model.output_transcript_id
        for accepted in selection_result.accepted_alternatives
    }
    missing_accepted = sorted(accepted_ids - seen_candidate_ids)
    if missing_accepted:
        raise ValueError(
            "Accepted alternatives are missing from the structural audit population: "
            + ", ".join(missing_accepted)
        )
    return audit_candidates


def _audit_attrs(candidate: AlternativeSpliceAuditCandidate) -> dict[str, object]:
    model = candidate.model
    row = candidate.fact_row
    decision = candidate.decision
    attrs: dict[str, object] = {
        "ID": candidate.audit_id,
        "moselecta_alt_gene": candidate.public_gene_id,
        "moselecta_alt_primary": candidate.public_primary_transcript_id,
        "moselecta_alt_candidate_order_rank": decision.alternative_candidate_order_rank,
        "moselecta_alt_preset": candidate.preset,
        "moselecta_alt_support_policy": candidate.support_policy,
        "moselecta_alt_evidence_requirement": candidate.evidence_requirement,
        "moselecta_alt_decision": decision.stage2b_decision,
        "moselecta_alt_rejection_reason": decision.stage2b_rejection_reason,
        "moselecta_alt_support_routes": (
            ",".join(decision.factual_support_routes)
            if decision.factual_support_routes
            else "none"
        ),
        "moselecta_alt_support_policy_eligible": _boolean_text(
            decision.support_policy_eligible
        ),
        "moselecta_alt_external_evidence_requirement_met": _boolean_text(
            decision.external_evidence_requirement_met
        ),
        "moselecta_alt_external_evidence_requirement_routes": (
            ",".join(decision.external_evidence_requirement_routes)
            if decision.external_evidence_requirement_routes
            else "none"
        ),
        "moselecta_alt_prethreshold_eligible": _boolean_text(
            decision.prethreshold_eligible
        ),
        "moselecta_alt_mosaic_policy_required": _boolean_text(
            decision.mosaic_policy_required
        ),
        "moselecta_alt_frame_concordance_status": (
            decision.stage2b_frame_concordance_status
        ),
        "moselecta_alt_primary_score_ratio_status": (
            decision.stage2b_primary_score_ratio_status
        ),
        "moselecta_alt_cds_base_union_containment": _containment_value(
            decision.candidate_cds_bases_contained_in_primary,
            decision.primary_cds_bases_contained_in_candidate,
        ),
        "moselecta_alt_relaxed_chain_containment": _containment_value(
            decision.candidate_relaxed_cds_chain_contained_in_primary,
            decision.primary_relaxed_cds_chain_contained_in_candidate,
        ),
        "moselecta_alt_full_structural_mosaic": _boolean_text(
            _truth(row.get("full_structural_mosaic"))
        ),
        "moselecta_alt_mosaic_components_complete": _boolean_text(
            decision.mosaic_components_complete
        ),
        "moselecta_alt_mosaic_additional_corroboration": _boolean_text(
            decision.mosaic_additional_corroboration
        ),
        "original_source": model.source_name,
        "original_gene_id": model.input_gene_id,
        "original_transcript_id": model.input_transcript_id,
        "moselecta_internal_gene_id": model.output_gene_id,
        "moselecta_internal_transcript_id": model.output_transcript_id,
    }
    for key, value in (
        (
            "moselecta_alt_reciprocal_cds_overlap",
            decision.stage2b_reciprocal_cds_overlap,
        ),
        (
            "moselecta_alt_frame_concordance",
            decision.stage2b_frame_concordance,
        ),
        (
            "moselecta_alt_model_selection_score_scaled",
            _optional_float(row.get("candidate_model_selection_score_scaled")),
        ),
        (
            "moselecta_alt_primary_score_ratio",
            decision.stage2b_candidate_primary_model_selection_score_ratio,
        ),
    ):
        if value is not None:
            attrs[key] = fmt_optional_float(value, 6)

    sequence_status = str(row.get("sequence_integrity_status", "."))
    if sequence_status != ".":
        attrs["moselecta_alt_sequence_status"] = sequence_status
    sequence_concerns = _concern_names(row, (
        "sequence_internal_stop_concern",
        "sequence_missing_start_concern",
        "sequence_missing_stop_concern",
        "sequence_phase_inconsistency_concern",
        "sequence_noncanonical_splice_concern",
    ))
    if sequence_concerns:
        attrs["moselecta_alt_sequence_concerns"] = ",".join(sequence_concerns)
    intrinsic_concerns = tuple(
        name
        for name in (
            "intrinsic_cds_length_not_multiple_of_3",
            "intrinsic_short_cds",
            "intrinsic_short_cds_intron",
            "intrinsic_short_internal_cds_exon",
            "intrinsic_short_terminal_cds_exon",
            "intrinsic_long_cds_intron",
        )
        if (_optional_float(row.get(name)) or 0.0) > 0.0
    )
    if intrinsic_concerns:
        attrs["moselecta_alt_intrinsic_concerns"] = ",".join(intrinsic_concerns)
    if decision.accepted_public_transcript_id is not None:
        attrs["moselecta_alt_public_transcript_id"] = (
            decision.accepted_public_transcript_id
        )
    if decision.greedy_conflict_public_transcript_id is not None:
        attrs["moselecta_alt_conflict_transcript"] = (
            decision.greedy_conflict_public_transcript_id
        )
    if decision.greedy_conflict_type is not None:
        attrs["moselecta_alt_conflict_type"] = decision.greedy_conflict_type
    return attrs


def write_alternative_splice_candidate_audit_gff3(
    *,
    candidate_rows: Iterable[Mapping[str, object]],
    selection_result: AlternativeSpliceSelectionResult,
    all_scored_models: Iterable[CodingModel],
    primary_annotation: FinalAnnotation,
    path: Path,
) -> tuple[AlternativeSpliceAuditCandidate, ...]:
    """Write the diagnostic browser track for uniquely owned structural candidates."""
    candidates = build_alternative_splice_audit_candidates(
        candidate_rows=candidate_rows,
        selection_result=selection_result,
        all_scored_models=all_scored_models,
        primary_annotation=primary_annotation,
    )
    with path.open("w") as out:
        out.write("##gff-version 3\n")
        for comment in alternative_splice_configuration_comments(
            selection_result.config,
            output_description=(
                "diagnostic candidate browser track; not retained gene annotation"
            ),
        ):
            out.write(comment + "\n")
        for candidate in candidates:
            model = candidate.model
            out.write(feature_line(
                seqid=model.chrom,
                source=f"{GFF3_SOURCE}:Audit:{model.source_name}",
                feature_type="mRNA",
                start0=model.model_start0,
                end0=model.model_end0,
                score=model.transcript_record.score if model.transcript_record else ".",
                strand=model.strand,
                phase=".",
                attrs=_audit_attrs(candidate),
            ) + "\n")
            write_source_model_children(out, model, candidate.audit_id)
            out.write("###\n")
    return candidates
