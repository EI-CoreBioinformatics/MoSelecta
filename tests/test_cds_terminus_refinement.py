from __future__ import annotations

from collections import Counter
import csv
import json
from pathlib import Path
import os
import subprocess
import sys

import pytest
import moselecta.pipeline as pipeline_module

from moselecta.annotation_outputs import (
    FinalAnnotationOutputMetadata,
    write_final_annotation_gff3,
)
from moselecta.cds_refinement import (
    CDSTerminusRefinementSettings,
    _FrozenRelationshipIndex,
    _SplicedPath,
    _simultaneous_conflicting_refinements,
    refine_cds_termini,
    write_cds_refinement_summary_tsv,
)
from moselecta.final_annotation import (
    AnnotationSummary,
    FinalAnnotation,
    FinalGene,
    FinalTranscript,
    public_gene_id,
    public_transcript_id,
)
from moselecta.genome_sequence import IndexedFasta
from moselecta.models import CodingModel, FeatureRecord, SequenceIntegrityMetrics
from moselecta.utr_augmentation import (
    UTREvidence,
    UTRLearnedLimits,
    UTRSettings,
    augment_final_annotation,
    learn_utr_limits_for_annotation,
)


def _record(
    feature: str, start0: int, end0: int, phase: str = ".", strand: str = "+"
) -> FeatureRecord:
    return FeatureRecord(
        seqid="chr1",
        source="prediction",
        feature_type=feature,
        start=start0 + 1,
        end=end0,
        score=".",
        strand=strand,
        phase=phase,
        attrs={"Parent": "tx"},
        raw_attrs=".",
        input_file=Path("input.gff3"),
        source_name="prediction",
        line_number=1,
        record_index=1,
        raw_line="",
    )


def _metrics(*, missing_start: bool = False, missing_stop: bool = False) -> SequenceIntegrityMetrics:
    flags = []
    if missing_start:
        flags.append("missing_start_codon")
    if missing_stop:
        flags.append("missing_stop_codon")
    return SequenceIntegrityMetrics(
        status="evaluated_with_concern" if flags else "evaluated_no_concern",
        start_codon="CCC" if missing_start else "ATG",
        start_codon_status="missing" if missing_start else "present",
        stop_codon="CCC" if missing_stop else "TAA",
        stop_codon_status="missing" if missing_stop else "included_in_cds",
        internal_stop_count=0,
        phase_status="consistent",
        phase_mismatch_count=0,
        noncanonical_splice_intron_count=0,
        noncanonical_splice_intron_fraction=0.0,
        noncanonical_splice_details=".",
        has_ambiguous_bases=False,
        penalty_score=0.3 * missing_start + 0.5 * missing_stop,
        penalty_components=".",
        flags=",".join(flags) or ".",
    )


def _model(
    parts: list[tuple[int, int]],
    *,
    missing_start: bool = False,
    missing_stop: bool = False,
    name: str = "selected",
    strand: str = "+",
    phases: list[int] | None = None,
) -> CodingModel:
    phase_values = phases or [0] * len(parts)
    children = [_record("exon", start, end, strand=strand) for start, end in parts]
    children += [
        _record("CDS", start, end, str(phase), strand=strand)
        for (start, end), phase in zip(parts, phase_values)
    ]
    return CodingModel(
        source_name="prediction",
        input_file=Path("input.gff3"),
        input_gene_id=name,
        input_transcript_id=f"{name}.t1",
        output_gene_id=f"{name}.gene",
        output_transcript_id=f"{name}.mrna",
        chrom="chr1",
        strand=strand,
        cds_parts=parts,
        gene_record=None,
        transcript_record=None,
        child_records=children,
        sequence_integrity=_metrics(
            missing_start=missing_start, missing_stop=missing_stop
        ),
    )


def _annotation(*models: CodingModel) -> FinalAnnotation:
    genes = []
    for index, model in enumerate(models, start=1):
        gene_id = public_gene_id(index)
        transcript = FinalTranscript(
            public_transcript_id(gene_id, 1),
            model,
            annotation_summary=AnnotationSummary(
                model_selection_score_scaled=0.77,
                sequence_integrity_status=model.sequence_integrity.status,
            ),
        )
        genes.append(FinalGene(gene_id, (transcript,)))
    return FinalAnnotation(tuple(genes))


def _evidence(
    evidence_id: str,
    *parts: tuple[int, int],
    sources: int = 2,
    strand: str = "+",
) -> UTREvidence:
    return UTREvidence(
        evidence_id=evidence_id,
        kind="transcript",
        chrom="chr1",
        strand=strand,
        exon_parts=tuple(parts),
        source_member_counts=tuple((f"source{i}", 1) for i in range(sources)),
    )


def _sequence(length: int, *values: tuple[int, str]) -> str:
    bases = ["C"] * length
    for start, value in values:
        bases[start : start + len(value)] = value
    return "".join(bases)


def _write_fasta(tmp_path: Path, sequence: str) -> Path:
    path = tmp_path / "genome.fa"
    path.write_text(">chr1\n" + sequence + "\n")
    return path


def _run(
    tmp_path: Path,
    model: CodingModel,
    sequence: str,
    evidence: list[UTREvidence],
    *,
    settings: CDSTerminusRefinementSettings | None = None,
    annotation: FinalAnnotation | None = None,
    utr_settings: UTRSettings | None = None,
    limits: UTRLearnedLimits | None = None,
):
    final = annotation or _annotation(model)
    with IndexedFasta(_write_fasta(tmp_path, sequence)) as fasta:
        return refine_cds_termini(
            final,
            evidence,
            [],
            fasta,
            settings or CDSTerminusRefinementSettings(min_cds_length_bp=50),
            utr_settings=utr_settings,
            learned_utr_limits=limits,
        )


def _summary_counts(result, field: str) -> dict[str, int]:
    return dict(getattr(result.summary, field))


def _assert_summary_reconciles_with_audit(result) -> None:
    rows = result.audit_rows
    summary = result.summary
    assert len({(row["transcript_id"], row["side"]) for row in rows}) == len(rows)
    assert summary.audit_sides == len(rows)
    assert summary.final_refined_sides + summary.final_failed_sides == len(rows)
    assert sum(_summary_counts(result, "final_outcome_side_counts").values()) == len(rows)
    assert sum(_summary_counts(result, "failure_counts").values()) == (
        summary.final_failed_sides
    )
    assert _summary_counts(result, "failure_counts") == dict(
        Counter(
            row["failure_reason"]
            for row in rows
            if row["result"] != "refined"
        )
    )
    for reason, transcript_count in _summary_counts(
        result, "failure_transcript_counts"
    ).items():
        assert transcript_count == len({
            row["transcript_id"]
            for row in rows
            if row["result"] != "refined" and row["failure_reason"] == reason
        })
    encountered = Counter()
    for row in rows:
        encountered.update(json.loads(row["encountered_rejection_counts"]))
        assert row["encountered_rejection_event_count"] == sum(
            json.loads(row["encountered_rejection_counts"]).values()
        )
    assert _summary_counts(result, "encountered_rejection_event_counts") == dict(
        encountered
    )


def test_simple_upstream_completion_and_immutable_selection_provenance(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True)
    annotation = _annotation(model)
    result = _run(
        tmp_path,
        model,
        _sequence(250, (70, "ATG"), (157, "TAA")),
        [_evidence("tx", (50, 160))],
        annotation=annotation,
    )
    refined = result.annotation.genes[0].transcripts[0]
    assert refined.source_model.cds_parts == [(70, 160)]
    assert model.cds_parts == [(100, 160)]
    assert refined.annotation_summary is annotation.genes[0].transcripts[0].annotation_summary
    assert result.summary.refined_sides_upstream_start == 1


def test_most_upstream_atg_after_preceding_stop_wins(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(
            250,
            (49, "ATG"),
            (55, "TGA"),
            (61, "ATG"),
            (91, "ATG"),
            (157, "TAA"),
        ),
        [_evidence("tx", (40, 160))],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_start0 == 61
    assert refined.cds_start0 != 49


def test_preceding_stop_clears_an_earlier_atg_on_minus_strand(
    tmp_path: Path,
) -> None:
    model = _model([(100, 160)], missing_start=True, strand="-")
    result = _run(
        tmp_path,
        model,
        _sequence(
            240,
            (100, "TTA"),  # reverse-complemented terminal TAA
            (196, "CAT"),  # first transcript-oriented ATG after the stop
            (202, "TCA"),  # reverse-complemented transcript-oriented TGA
            (208, "CAT"),  # earlier ATG that must not survive the stop
        ),
        [_evidence("minus_path", (100, 220), strand="-")],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(100, 199)]
    assert refined.cds_end0 != 211


def test_open_search_boundary_uses_most_upstream_visible_atg_and_audits_open_orf(tmp_path: Path) -> None:
    model = _model([(150, 210)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(260, (51, "ATG"), (120, "ATG"), (207, "TAA")),
        [_evidence("tx", (40, 210))],
        settings=CDSTerminusRefinementSettings(max_relative_extension=2.0),
    )
    assert result.annotation.genes[0].transcripts[0].source_model.cds_start0 == 51
    assert result.audit_rows[0]["orf_open_at_search_boundary"] == "true"


def test_farther_atg_beyond_former_100_bp_window_is_the_only_upstream_choice(
    tmp_path: Path,
) -> None:
    model = _model([(200, 380)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(
            430,
            (50, "ATG"),
            (140, "ATG"),
            (377, "TAA"),
        ),
        [_evidence("full_path", (20, 380))],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_start0 == 50
    assert result.audit_rows[0]["transcript_bases_searched"] > 100
    assert result.audit_rows[0]["configured_search_limit_bp"] == "."


def test_million_base_supported_structure_uses_part_scaled_coordinate_state(
    tmp_path: Path,
) -> None:
    structure_length = 1_000_001
    model = _model([(500_000, structure_length)], missing_start=True)
    sequence = "C" * structure_length
    sequence = sequence[:101] + "ATG" + sequence[104:999_998] + "TAA"
    result = _run(
        tmp_path,
        model,
        sequence,
        [_evidence("million_base_path", (0, structure_length))],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_start0 == 101
    assert refined.cds_end0 == structure_length
    assert result.audit_rows[0]["transcript_bases_searched"] == 499_998

    path = _SplicedPath.from_parts(
        ((0, 500_000), (600_000, 1_100_001)), "+"
    )
    assert path.length == structure_length
    assert len(path.parts) == 2
    assert len(path.offsets) == 3
    assert not hasattr(path, "positions")


def test_farther_multiexon_start_crosses_splice_on_plus_strand(
    tmp_path: Path,
) -> None:
    model = _model([(300, 360), (400, 520)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(
            560,
            (130, "ATG"),
            (220, "ATG"),
            (250, "GT"),
            (298, "AG"),
            (360, "GT"),
            (398, "AG"),
            (517, "TAA"),
        ),
        [_evidence("spliced_plus", (100, 250), (300, 360), (400, 520))],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(130, 250), (300, 360), (400, 520)]
    assert result.audit_rows[0]["splice_junctions_crossed"] == 1


def test_farther_start_search_is_symmetric_on_minus_strand(
    tmp_path: Path,
) -> None:
    model = _model(
        [(100, 220), (260, 320)], missing_start=True, strand="-"
    )
    result = _run(
        tmp_path,
        model,
        _sequence(
            540,
            (100, "TTA"),
            (220, "CT"),
            (258, "AC"),
            (350, "CAT"),
            (440, "CAT"),
        ),
        [_evidence("spliced_minus", (100, 220), (260, 500), strand="-")],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(100, 220), (260, 443)]
    assert result.audit_rows[0]["transcript_bases_searched"] > 100


def test_upstream_splice_completion_and_codon_crossing_junction(tmp_path: Path) -> None:
    model = _model([(154, 160)], missing_start=True)
    sequence = _sequence(
        220,
        (99, "A"),
        (100, "GT"),
        (150, "AG"),
        (152, "T"),
        (153, "G"),
        (157, "TAA"),
    )
    result = _run(tmp_path, model, sequence, [_evidence("tx", (90, 100), (152, 160))])
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(99, 100), (152, 160)]


def test_minus_strand_start_codon_can_cross_a_splice_junction(
    tmp_path: Path,
) -> None:
    model = _model([(90, 96)], missing_start=True, strand="-")
    result = _run(
        tmp_path,
        model,
        _sequence(
            220,
            (90, "TTA"),
            (96, "CA"),
            (98, "CT"),
            (148, "AC"),
            (150, "T"),
        ),
        [_evidence("minus_spliced", (90, 98), (150, 160), strand="-")],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(90, 98), (150, 151)]
    assert result.audit_rows[0]["codon_genomic_positions"] == "151,98,97"


def test_terminal_compatibility_can_complete_without_evidence_covering_full_cds(
    tmp_path: Path,
) -> None:
    model = _model([(100, 130), (200, 230)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(
            360,
            (91, "ATG"),
            (130, "GT"),
            (198, "AG"),
            (227, "TAA"),
        ),
        [_evidence("terminal", (70, 130), (300, 330))],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(91, 130), (200, 230)]


@pytest.mark.parametrize("has_structure", [False, True])
def test_downstream_fallback_without_start_or_without_upstream_atg(
    tmp_path: Path, has_structure: bool
) -> None:
    model = _model([(100, 160)], missing_start=True)
    evidence = [_evidence("tx", (70, 160))] if has_structure else []
    result = _run(
        tmp_path,
        model,
        _sequence(220, (109, "ATG"), (157, "TAA")),
        evidence,
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_start0 == 109
    assert result.summary.refined_sides_downstream_start_fallback == 1


def test_upstream_success_precedes_downstream_fallback(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(220, (70, "ATG"), (109, "ATG"), (157, "TAA")),
        [_evidence("tx", (50, 160))],
    )
    assert result.annotation.genes[0].transcripts[0].source_model.cds_start0 == 70


def test_relative_extension_and_truncation_guards(tmp_path: Path) -> None:
    upstream = _model([(100, 130)], missing_start=True, name="up")
    upstream_result = _run(
        tmp_path,
        upstream,
        _sequence(200, (40, "ATG"), (127, "TAA")),
        [_evidence("tx", (30, 130))],
        settings=CDSTerminusRefinementSettings(max_relative_extension=0.5),
    )
    assert upstream_result.annotation.genes[0].transcripts[0].source_model is upstream
    assert _summary_counts(upstream_result, "failure_counts") == {
        "cds_extension_relative_change_exceeded": 1
    }

    downstream = _model([(100, 130)], missing_start=True, name="down")
    downstream_result = _run(
        tmp_path,
        downstream,
        _sequence(200, (112, "ATG"), (127, "TAA")),
        [],
        settings=CDSTerminusRefinementSettings(
            downstream_start_search_bp=30, max_relative_truncation=0.2
        ),
    )
    assert downstream_result.annotation.genes[0].transcripts[0].source_model is downstream
    assert _summary_counts(downstream_result, "failure_counts") == {
        "cds_truncation_relative_change_exceeded": 1
    }


def test_rejected_true_upstream_start_does_not_use_internal_or_downstream_fallback(
    tmp_path: Path,
) -> None:
    model = _model([(100, 130)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(200, (40, "ATG"), (109, "ATG"), (127, "TAA")),
        [_evidence("tx", (30, 130))],
        settings=CDSTerminusRefinementSettings(
            max_relative_extension=0.5,
            max_relative_truncation=0.33,
        ),
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined is model
    assert _summary_counts(result, "failure_counts") == {
        "cds_extension_relative_change_exceeded": 1
    }
    assert result.summary.refined_sides_downstream_start_fallback == 0
    assert result.audit_rows[0]["upstream_completion_failure_reason"] == "cds_extension_relative_change_exceeded"
    assert result.audit_rows[0]["completion_direction"] == "upstream_start"

    summary_path = tmp_path / "summary.tsv"
    write_cds_refinement_summary_tsv(result, summary_path)
    with summary_path.open() as handle:
        summary_row = next(csv.DictReader(handle, delimiter="\t"))
    assert json.loads(summary_row["failure_counts"]) == {
        "cds_extension_relative_change_exceeded": 1
    }
    assert json.loads(summary_row["failure_transcript_counts"]) == {
        "cds_extension_relative_change_exceeded": 1
    }
    assert json.loads(summary_row["final_outcome_side_counts"]) == {
        "cds_extension_relative_change_exceeded": 1
    }
    assert json.loads(summary_row["encountered_rejection_event_counts"]) == {
        "cds_extension_relative_change_exceeded": 1
    }
    assert summary_row["audit_sides"] == "1"
    assert summary_row["final_failed_sides"] == "1"
    assert "rejected_structural_conflict" not in summary_row
    assert "rejected_integrity" not in summary_row
    assert "structural_rejections_encountered" not in summary_row


def test_rejected_upstream_splice_proposal_does_not_run_downstream_fallback(
    tmp_path: Path,
) -> None:
    model = _model([(154, 220)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(
            260,
            (93, "ATG"),
            (163, "ATG"),
            (217, "TAA"),
        ),
        [_evidence("noncanonical_path", (90, 100), (152, 220))],
    )
    assert result.annotation.genes[0].transcripts[0].source_model is model
    assert result.summary.refined_sides_downstream_start_fallback == 0
    assert result.audit_rows[0]["failure_reason"] == "new_noncanonical_cds_splice"
    assert result.audit_rows[0]["completion_direction"] == "upstream_start"


def test_same_exon_first_stop_completion(tmp_path: Path) -> None:
    model = _model([(100, 130)], missing_stop=True)
    result = _run(
        tmp_path,
        model,
        _sequence(220, (100, "ATG"), (133, "TAG"), (142, "TGA")),
        [_evidence("tx", (100, 180))],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(100, 136)]
    assert result.summary.refined_sides_downstream_stop == 1


def test_stop_requires_supported_structure(tmp_path: Path) -> None:
    model = _model([(100, 130)], missing_stop=True)
    result = _run(tmp_path, model, _sequence(180, (100, "ATG"), (133, "TAA")), [])
    assert result.annotation.genes[0].transcripts[0].source_model is model
    assert result.audit_rows[0]["failure_reason"] == "no_compatible_completion_structure"


def test_stop_across_canonical_splice(tmp_path: Path) -> None:
    model = _model([(100, 130)], missing_stop=True)
    sequence = _sequence(
        240,
        (100, "ATG"),
        (130, "T"),
        (131, "GT"),
        (178, "AG"),
        (180, "GA"),
    )
    result = _run(tmp_path, model, sequence, [_evidence("tx", (100, 131), (180, 210))])
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(100, 131), (180, 182)]


def test_ambiguity_rule_allows_single_n_and_nn_is_an_upstream_search_barrier(
    tmp_path: Path,
) -> None:
    model = _model([(100, 160)], missing_start=True)
    single = _run(
        tmp_path,
        model,
        _sequence(220, (70, "ATG"), (80, "N"), (157, "TAA")),
        [_evidence("tx", (50, 160))],
    )
    assert single.summary.refined_transcripts == 1

    blocked = _run(
        tmp_path,
        model,
        _sequence(
            220,
            (70, "ATG"),
            (80, "NN"),
            (109, "ATG"),
            (157, "TAA"),
        ),
        [_evidence("tx", (50, 160))],
    )
    blocked_refined = blocked.annotation.genes[0].transcripts[0].source_model
    assert blocked_refined.cds_start0 == 109
    assert blocked.summary.refined_sides_downstream_start_fallback == 1
    assert _summary_counts(blocked, "encountered_rejection_event_counts").get(
        "coding_extension_crosses_ambiguous_run", 0
    ) == 0

    spliced_model = _model([(154, 160)], missing_start=True, name="spliced")
    intronic = _run(
        tmp_path,
        spliced_model,
        _sequence(
            220,
                (99, "A"),
                (100, "GT"),
                (110, "NN"),
                (150, "AG"),
                (152, "T"),
                (153, "G"),
                (157, "TAA"),
            ),
            [_evidence("tx2", (90, 100), (152, 160))],
    )
    assert intronic.summary.refined_transcripts == 1


def test_new_frozen_gene_splice_relationship_rolls_back(tmp_path: Path) -> None:
    selected = _model([(154, 160)], missing_start=True, name="selected")
    other = _model([(90, 100), (152, 180)], name="other")
    annotation = _annotation(selected, other)
    result = _run(
        tmp_path,
        selected,
        _sequence(
            240,
            (99, "A"),
            (100, "GT"),
            (150, "AG"),
            (152, "T"),
            (153, "G"),
            (157, "TAA"),
        ),
        [_evidence("tx", (90, 100), (152, 160))],
        annotation=annotation,
    )
    assert result.annotation.genes[0].transcripts[0].source_model is selected
    assert result.summary.final_failed_sides == 1


def test_new_canonical_single_exon_overlap_relationship_rolls_back(tmp_path: Path) -> None:
    selected = _model([(100, 130)], missing_stop=True, name="selected")
    other = _model([(148, 165)], name="other")
    annotation = _annotation(selected, other)
    result = _run(
        tmp_path,
        selected,
        _sequence(220, (100, "ATG"), (157, "TAA")),
        [_evidence("tx", (100, 180))],
        annotation=annotation,
    )
    assert result.annotation.genes[0].transcripts[0].source_model is selected
    assert _summary_counts(result, "failure_counts") == {
        "frozen_final_gene_structural_conflict": 1
    }
    assert result.audit_rows[0]["failure_reason"] == "frozen_final_gene_structural_conflict"


def test_upstream_extension_new_same_strand_cross_gene_cds_overlap_is_rejected(
    tmp_path: Path,
) -> None:
    selected = _model([(100, 160)], missing_start=True, name="selected")
    other = _model([(50, 80), (200, 230)], name="other")
    result = _run(
        tmp_path,
        selected,
        _sequence(260, (70, "ATG"), (109, "ATG"), (157, "TAA")),
        [_evidence("target_path", (40, 160))],
        annotation=_annotation(selected, other),
    )
    assert result.annotation.genes[0].transcripts[0].source_model is selected
    assert result.summary.refined_sides_downstream_start_fallback == 0
    assert _summary_counts(result, "failure_counts") == {
        "new_same_strand_cross_gene_cds_overlap": 1
    }
    assert result.audit_rows[0]["failure_reason"] == (
        "new_same_strand_cross_gene_cds_overlap"
    )


def test_stop_extension_new_same_strand_cross_gene_cds_overlap_is_rejected(
    tmp_path: Path,
) -> None:
    selected = _model([(100, 130)], missing_stop=True, name="selected")
    other = _model([(148, 165), (200, 230)], name="other")
    result = _run(
        tmp_path,
        selected,
        _sequence(260, (100, "ATG"), (157, "TAA")),
        [_evidence("target_path", (100, 180))],
        annotation=_annotation(selected, other),
    )
    assert result.annotation.genes[0].transcripts[0].source_model is selected
    assert _summary_counts(result, "failure_counts") == {
        "new_same_strand_cross_gene_cds_overlap": 1
    }
    assert _summary_counts(result, "failure_transcript_counts") == {
        "new_same_strand_cross_gene_cds_overlap": 1
    }
    assert _summary_counts(result, "encountered_rejection_event_counts") == {
        "new_same_strand_cross_gene_cds_overlap": 1
    }
    _assert_summary_reconciles_with_audit(result)
    assert result.audit_rows[0]["failure_reason"] == (
        "new_same_strand_cross_gene_cds_overlap"
    )


def test_preexisting_cross_gene_cds_overlap_is_not_retroactively_rejected(
    tmp_path: Path,
) -> None:
    selected = _model([(100, 160)], missing_stop=True, name="selected")
    other = _model([(120, 130), (200, 230)], name="other")
    result = _run(
        tmp_path,
        selected,
        _sequence(260, (100, "ATG"), (169, "TAA")),
        [_evidence("target_path", (100, 190))],
        annotation=_annotation(selected, other),
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_end0 == 172
    assert _summary_counts(result, "failure_counts").get(
        "new_same_strand_cross_gene_cds_overlap", 0
    ) == 0


def test_other_gene_utr_overlap_alone_does_not_block_cds_extension(
    tmp_path: Path,
) -> None:
    selected = _model([(100, 130)], missing_stop=True, name="selected")
    other = _model([(200, 230)], name="other")
    other.child_records.append(_record("five_prime_UTR", 140, 180))
    result = _run(
        tmp_path,
        selected,
        _sequence(260, (100, "ATG"), (157, "TAA"), (200, "ATG"), (227, "TAA")),
        [_evidence("target_path", (100, 180))],
        annotation=_annotation(selected, other),
    )
    assert result.annotation.genes[0].transcripts[0].source_model.cds_end0 == 160
    assert _summary_counts(result, "failure_counts").get(
        "new_same_strand_cross_gene_cds_overlap", 0
    ) == 0


def test_opposite_strand_cds_overlap_does_not_block_extension(
    tmp_path: Path,
) -> None:
    selected = _model([(100, 160)], missing_start=True, name="selected")
    other = _model(
        [(50, 80), (200, 230)], name="other", strand="-"
    )
    result = _run(
        tmp_path,
        selected,
        _sequence(260, (70, "ATG"), (157, "TAA")),
        [_evidence("target_path", (40, 160))],
        annotation=_annotation(selected, other),
    )
    assert result.annotation.genes[0].transcripts[0].source_model.cds_start0 == 70


def test_same_gene_alternative_cds_overlap_does_not_block_extension(
    tmp_path: Path,
) -> None:
    selected = _model([(100, 160)], missing_start=True, name="selected")
    alternative = _model([(50, 80), (200, 230)], name="alternative")
    gene_id = "MoSelecta_000001"
    annotation = FinalAnnotation((FinalGene(gene_id, (
        FinalTranscript(public_transcript_id(gene_id, 1), selected),
        FinalTranscript(
            public_transcript_id(gene_id, 2), alternative, transcript_role="alternative"
        ),
    )),))
    result = _run(
        tmp_path,
        selected,
        _sequence(260, (70, "ATG"), (157, "TAA")),
        [_evidence("target_path", (40, 160))],
        annotation=annotation,
    )
    refined_gene = result.annotation.genes[0]
    assert refined_gene.transcripts[0].source_model.cds_start0 == 70
    assert refined_gene.transcripts[1].source_model is alternative


@pytest.mark.parametrize("reverse_order", [False, True])
def test_refined_vs_refined_new_cds_overlap_rolls_back_both_deterministically(
    tmp_path: Path,
    reverse_order: bool,
) -> None:
    left = _model(
        [(100, 161)], missing_start=True, name="left", phases=[1]
    )
    right = _model([(40, 68)], missing_stop=True, name="right")
    models = (right, left) if reverse_order else (left, right)
    result = _run(
        tmp_path,
        models[0],
        _sequence(
            220,
            (40, "ATG"),
            (68, "ATG"),
            (97, "TAA"),
            (158, "TAA"),
        ),
        [
            _evidence("left_path", (50, 161)),
            _evidence("right_path", (40, 120)),
        ],
        annotation=_annotation(*models),
        settings=CDSTerminusRefinementSettings(max_relative_extension=2.0),
    )
    assert result.summary.refined_transcripts == 0
    assert _summary_counts(result, "failure_counts") == {
        "new_same_strand_cross_gene_cds_overlap": 2
    }
    assert {row["failure_reason"] for row in result.audit_rows} == {
        "new_same_strand_cross_gene_cds_overlap"
    }


def test_residual_utr_is_audited_without_affecting_cds_commit(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True)
    limits = UTRLearnedLimits(
        five_prime_exonic_bp=50,
        three_prime_exonic_bp=50,
        five_prime_genomic_bp=None,
        three_prime_genomic_bp=None,
        five_prime_exonic_source="test",
        three_prime_exonic_source="test",
        five_prime_genomic_source="unavailable",
        three_prime_genomic_source="unavailable",
    )
    result = _run(
        tmp_path,
        model,
        _sequence(220, (61, "ATG"), (157, "TAA")),
        [_evidence("tx", (20, 160))],
        utr_settings=UTRSettings(),
        limits=limits,
    )
    assert result.summary.refined_transcripts == 1
    assert result.audit_rows[0]["residual_utr_parts"] == "21-61"


def test_refined_geometry_and_recomputed_phase_are_emitted(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(220, (70, "ATG"), (157, "TAA")),
        [_evidence("tx", (50, 160))],
    )
    output = tmp_path / "final.gff3"
    write_final_annotation_gff3(
        result.annotation,
        output,
        FinalAnnotationOutputMetadata(False, False, False),
        cds_only_structure=False,
    )
    text = output.read_text()
    assert "\tCDS\t71\t160\t.\t+\t0\t" in text
    assert "\tmRNA\t71\t160\t" in text
    assert "\tgene\t71\t160\t" in text
    assert "\tCDS\t101\t160\t" not in text
    assert "cds_terminus_refined=upstream_start" in text


def test_candidate_tie_is_deterministic_across_insertion_order(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True)
    sequence = _sequence(220, (70, "ATG"), (80, "ATG"), (157, "TAA"))
    left = _evidence("a", (50, 160), sources=2)
    right = _evidence("b", (60, 160), sources=2)
    first = _run(tmp_path, model, sequence, [right, left])
    second = _run(tmp_path, model, sequence, [left, right])
    assert first.audit_rows[0]["candidate_evidence_id"] == second.audit_rows[0]["candidate_evidence_id"]
    assert first.annotation.genes[0].transcripts[0].source_model.cds_parts == second.annotation.genes[0].transcripts[0].source_model.cds_parts


def test_failed_preferred_structure_does_not_try_lower_ranked_codon_path(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True)
    sequence = _sequence(
        220,
        (32, "ATG"),
        (35, "GT"),
        (98, "AG"),
        (157, "TAA"),
    )
    preferred = _evidence("preferred", (40, 160), sources=3)
    lower = _evidence("lower", (30, 35), (100, 160), sources=1)
    result = _run(tmp_path, model, sequence, [lower, preferred])
    assert result.summary.refined_transcripts == 0
    assert result.audit_rows[0]["candidate_evidence_id"] == "preferred"

    lower_only = _run(tmp_path, model, sequence, [lower])
    assert lower_only.summary.refined_transcripts == 1
    assert lower_only.annotation.genes[0].transcripts[0].source_model.cds_start0 == 32


@pytest.mark.parametrize(
    ("start_at", "stop_at", "expected_start", "expected_stop", "refined"),
    [
        (70, 163, 70, 166, 1),
        (70, None, 70, 160, 1),
        (None, 163, 100, 166, 1),
        (None, None, 100, 160, 0),
    ],
)
def test_missing_both_sides_are_independent(
    tmp_path: Path,
    start_at: int | None,
    stop_at: int | None,
    expected_start: int,
    expected_stop: int,
    refined: int,
) -> None:
    model = _model([(100, 160)], missing_start=True, missing_stop=True)
    values = []
    if start_at is not None:
        values.append((start_at, "ATG"))
    if stop_at is not None:
        values.append((stop_at, "TAA"))
    result = _run(
        tmp_path,
        model,
        _sequence(240, *values),
        [_evidence("tx", (50, 200))],
    )
    final = result.annotation.genes[0].transcripts[0].source_model
    assert (final.cds_start0, final.cds_end0) == (expected_start, expected_stop)
    assert result.summary.refined_transcripts == refined
    assert result.summary.transcripts_with_two_incomplete_sides == 1
    outcomes = (
        result.summary.two_incomplete_sides_both_refined_transcripts,
        result.summary.two_incomplete_sides_start_only_refined_transcripts,
        result.summary.two_incomplete_sides_stop_only_refined_transcripts,
        result.summary.two_incomplete_sides_neither_refined_transcripts,
    )
    expected_outcome = (
        (1, 0, 0, 0)
        if start_at is not None and stop_at is not None
        else (0, 1, 0, 0)
        if start_at is not None
        else (0, 0, 1, 0)
        if stop_at is not None
        else (0, 0, 0, 1)
    )
    assert outcomes == expected_outcome
    assert result.summary.final_refined_sides == int(start_at is not None) + int(
        stop_at is not None
    )
    _assert_summary_reconciles_with_audit(result)


def test_combined_integrity_failure_rolls_back_both_sides(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True, missing_stop=True)
    result = _run(
        tmp_path,
        model,
        _sequence(240, (70, "ATG"), (130, "TGA"), (163, "TAA")),
        [_evidence("tx", (50, 200))],
    )
    assert result.annotation.genes[0].transcripts[0].source_model is model
    assert _summary_counts(result, "failure_counts") == {
        "refined_cds_integrity_failed": 2
    }
    assert _summary_counts(result, "failure_transcript_counts") == {
        "refined_cds_integrity_failed": 1
    }
    assert _summary_counts(result, "encountered_rejection_event_counts") == {
        "refined_cds_integrity_failed": 4
    }
    assert _summary_counts(result, "encountered_rejection_side_counts") == {
        "refined_cds_integrity_failed": 2
    }
    assert _summary_counts(result, "encountered_rejection_transcript_counts") == {
        "refined_cds_integrity_failed": 1
    }
    assert {row["result"] for row in result.audit_rows} == {"rolled_back"}
    _assert_summary_reconciles_with_audit(result)


def test_stop_relative_change_rejection_preserves_discovered_codon_in_audit(tmp_path: Path) -> None:
    model = _model([(100, 106)], missing_stop=True)
    result = _run(
        tmp_path,
        model,
        _sequence(200, (100, "ATG"), (112, "TAA")),
        [_evidence("tx", (100, 150))],
        settings=CDSTerminusRefinementSettings(max_relative_extension=0.5),
    )
    assert result.summary.refined_transcripts == 0
    assert result.audit_rows[0]["failure_reason"] == "cds_extension_relative_change_exceeded"
    assert result.audit_rows[0]["codon_found"] == "TAA"
    assert result.audit_rows[0]["coding_bp_added"] == 9


def test_start_or_stop_codon_containing_n_is_not_accepted(tmp_path: Path) -> None:
    start_model = _model([(100, 160)], missing_start=True)
    start = _run(
        tmp_path,
        start_model,
        _sequence(220, (70, "ATN"), (157, "TAA")),
        [_evidence("tx", (50, 160))],
    )
    assert start.summary.refined_transcripts == 0

    stop_model = _model([(100, 130)], missing_stop=True, name="stop")
    stop = _run(
        tmp_path,
        stop_model,
        _sequence(220, (100, "ATG"), (130, "TAN")),
        [_evidence("tx2", (100, 180))],
    )
    assert stop.summary.refined_transcripts == 0


def test_residual_utr_failure_does_not_roll_back_valid_cds_refinement(
    tmp_path: Path,
) -> None:
    model = _model([(100, 160)], missing_start=True)
    limits = UTRLearnedLimits(
        five_prime_exonic_bp=20,
        three_prime_exonic_bp=50,
        five_prime_genomic_bp=None,
        three_prime_genomic_bp=None,
        five_prime_exonic_source="test",
        three_prime_exonic_source="test",
        five_prime_genomic_source="unavailable",
        three_prime_genomic_source="unavailable",
    )
    without_utr_limits = _run(
        tmp_path,
        model,
        _sequence(220, (61, "ATG"), (157, "TAA")),
        [_evidence("tx", (20, 160))],
    )
    result = _run(
        tmp_path,
        model,
        _sequence(220, (61, "ATG"), (157, "TAA")),
        [_evidence("tx", (20, 160))],
        utr_settings=UTRSettings(),
        limits=limits,
    )
    assert result.summary.refined_transcripts == 1
    assert (
        result.annotation.genes[0].transcripts[0].source_model.cds_parts
        == without_utr_limits.annotation.genes[0].transcripts[0].source_model.cds_parts
        == [(61, 160)]
    )
    augmented = augment_final_annotation(
        result.annotation,
        [_evidence("tx", (20, 160))],
        [],
        UTRSettings(),
        learned_limits=limits,
        excluded_training_sides=set(result.originally_incomplete_sides),
    )
    assert augmented.annotation.genes[0].transcripts[0].five_prime_utr_parts == ()


def test_no_utr_writer_emits_refined_cds_without_utr_features(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(220, (70, "ATG"), (157, "TAA")),
        [_evidence("tx", (50, 160))],
    )
    output = tmp_path / "refined.gff3"
    write_final_annotation_gff3(
        result.annotation,
        output,
        FinalAnnotationOutputMetadata(False, False, False),
        cds_only_structure=True,
    )
    text = output.read_text()
    assert "\tCDS\t71\t160\t" in text
    assert "five_prime_UTR" not in text and "three_prime_UTR" not in text


def test_incomplete_side_does_not_train_limits_but_refined_side_can_consume_them(
    tmp_path: Path,
) -> None:
    incomplete = _model([(100, 160)], missing_start=True)
    complete = _model([(300, 360)], name="complete")
    annotation = _annotation(incomplete, complete)
    evidence = [
        _evidence("incomplete_tx", (20, 180)),
        _evidence("complete_tx", (280, 380)),
    ]
    settings = UTRSettings(limit_min_observations=1)
    limits = learn_utr_limits_for_annotation(annotation, evidence, [], settings)
    complete_only_limits = learn_utr_limits_for_annotation(
        annotation, [evidence[1]], [], settings
    )
    assert limits.five_prime_exonic_bp == complete_only_limits.five_prime_exonic_bp

    strict_limits = UTRLearnedLimits(
        five_prime_exonic_bp=20,
        three_prime_exonic_bp=20,
        five_prime_genomic_bp=None,
        three_prime_genomic_bp=None,
        five_prime_exonic_source="test",
        three_prime_exonic_source="test",
        five_prime_genomic_source="unavailable",
        three_prime_genomic_source="unavailable",
    )

    result = _run(
        tmp_path,
        incomplete,
        _sequence(420, (61, "ATG"), (157, "TAA"), (300, "ATG"), (357, "TAA")),
        evidence,
        annotation=annotation,
        utr_settings=settings,
        limits=strict_limits,
    )
    assert result.summary.refined_transcripts == 1
    assert result.annotation.genes[0].transcripts[0].source_model.cds_start0 == 61
    augmented = augment_final_annotation(
        result.annotation,
        evidence,
        [],
        settings,
        learned_limits=strict_limits,
        excluded_training_sides=set(result.originally_incomplete_sides),
    )
    assert augmented.annotation.genes[0].transcripts[0].five_prime_utr_parts == ()


def test_refined_annotation_can_feed_ordinary_residual_utr_augmentation(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True)
    evidence = [_evidence("tx", (50, 180))]
    refined = _run(
        tmp_path,
        model,
        _sequence(220, (70, "ATG"), (157, "TAA")),
        evidence,
    )
    augmented = augment_final_annotation(
        refined.annotation,
        evidence,
        [],
        UTRSettings(
            five_prime_exonic_override_bp=100,
            three_prime_exonic_override_bp=100,
        ),
        excluded_training_sides=set(refined.originally_incomplete_sides),
    )
    transcript = augmented.annotation.genes[0].transcripts[0]
    assert transcript.model.cds_start0 == 70
    assert transcript.five_prime_utr_parts == ((50, 70),)


def test_minus_strand_upstream_completion_is_symmetric(tmp_path: Path) -> None:
    model = _model([(100, 160)], missing_start=True, strand="-")
    result = _run(
        tmp_path,
        model,
        _sequence(240, (100, "TTA"), (187, "CAT")),
        [_evidence("tx", (100, 190), strand="-")],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(100, 190)]
    assert refined.sequence_integrity.start_codon_status == "present"


def test_joint_refinements_creating_new_shared_splice_are_both_rolled_back(
    tmp_path: Path,
) -> None:
    left = _model([(154, 160)], missing_start=True, name="left")
    right = _model([(180, 186)], missing_start=True, name="right")
    annotation = _annotation(left, right)
    result = _run(
        tmp_path,
        left,
        _sequence(
            260,
            (99, "A"),
            (100, "GT"),
            (150, "AG"),
            (152, "TG"),
            (157, "TAA"),
            (176, "AG"),
            (178, "TG"),
            (183, "TAA"),
        ),
        [
            _evidence("left_path", (90, 100), (152, 160)),
            _evidence("right_path", (90, 100), (178, 186)),
        ],
        annotation=annotation,
    )
    assert result.annotation.genes[0].transcripts[0].source_model is left
    assert result.annotation.genes[1].transcripts[0].source_model is right
    assert result.summary.refined_transcripts == 0
    assert result.summary.final_failed_sides == 2
    assert {row["failure_reason"] for row in result.audit_rows} == {
        "new_same_strand_cross_gene_cds_overlap"
    }
    assert _summary_counts(result, "failure_counts") == {
        "new_same_strand_cross_gene_cds_overlap": 2
    }
    assert {row["final_gene_structural_check"] for row in result.audit_rows} == {
        "fail"
    }


def test_population_validation_queries_only_spatial_neighbours(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    models = [
        _model([(0, 10)], name="target"),
        _model([(20, 30)], name="nearby"),
        *[
            _model([(index * 100, index * 100 + 10)], name=f"far{index}")
            for index in range(1, 1001)
        ],
    ]
    annotation = _annotation(*models)
    replacement = _model([(0, 25)], name="target_refined")
    calls = 0
    import moselecta.cds_refinement as refinement_module

    real_predicate = refinement_module.models_share_sublocus_relationship

    def counted_predicate(*args, **kwargs):
        nonlocal calls
        calls += 1
        return real_predicate(*args, **kwargs)

    monkeypatch.setattr(
        refinement_module, "models_share_sublocus_relationship", counted_predicate
    )
    frozen = _FrozenRelationshipIndex(annotation, 0.2)
    conflicting = _simultaneous_conflicting_refinements(
        annotation,
        {annotation.genes[0].transcripts[0].transcript_id: replacement},
        frozen,
    )
    assert conflicting == {annotation.genes[0].transcripts[0].transcript_id}
    assert calls < 20


@pytest.mark.parametrize(
    "order",
    [
        ("left", "middle", "right"),
        ("right", "left", "middle"),
        ("middle", "right", "left"),
    ],
)
def test_three_way_indirect_refinement_bridge_rolls_back_all_independent_of_order(
    tmp_path: Path, order: tuple[str, str, str]
) -> None:
    models = {
        "left": _model([(154, 160)], missing_start=True, name="left"),
        "middle": _model([(298, 307)], missing_start=True, name="middle"),
        "right": _model([(204, 207)], missing_start=True, name="right"),
    }
    annotation = _annotation(*(models[name] for name in order))
    result = _run(
        tmp_path,
        models[order[0]],
        _sequence(
            360,
            (94, "ATG"), (99, "A"), (100, "GT"),
            (126, "ATG"), (130, "GT"),
            (150, "AG"), (152, "TG"), (157, "TAA"),
            (176, "AG"), (178, "TG"),
            (199, "GT"), (296, "AG"), (304, "TAA"),
            (204, "TAA"),
        ),
        [
            _evidence("left_path", (90, 100), (152, 160)),
            _evidence("middle_path", (90, 100), (178, 199), (298, 307)),
            _evidence("right_path", (120, 130), (178, 207)),
        ],
        annotation=annotation,
        settings=CDSTerminusRefinementSettings(max_relative_extension=20.0),
    )
    assert result.summary.refined_transcripts == 0
    assert _summary_counts(result, "failure_counts") == {
        "new_same_strand_cross_gene_cds_overlap": 3
    }
    assert all(
        transcript.source_model is models[name]
        for gene, name in zip(result.annotation.genes, order)
        for transcript in gene.transcripts
    )
    assert {row["failure_reason"] for row in result.audit_rows} == {
        "new_same_strand_cross_gene_cds_overlap"
    }


def test_bridge_detection_includes_complete_final_genes(tmp_path: Path) -> None:
    incomplete = _model([(151, 160)], missing_start=True, name="incomplete")
    complete = _model([(90, 120)], name="complete")
    result = _run(
        tmp_path,
        incomplete,
        _sequence(240, (120, "ATG"), (157, "TAA")),
        [_evidence("bridging", (80, 180))],
        annotation=_annotation(incomplete, complete),
        settings=CDSTerminusRefinementSettings(max_relative_extension=5.0),
    )
    assert result.annotation.genes[0].transcripts[0].source_model is incomplete
    assert result.audit_rows[0]["failure_reason"] == (
        "frozen_final_gene_structural_conflict"
    )
    assert result.audit_rows[0]["final_gene_structural_check"] == "not_run"


def test_intermediate_structural_rejection_is_not_reported_as_final_failure(
    tmp_path: Path,
) -> None:
    incomplete = _model([(151, 211)], missing_start=True, name="incomplete")
    complete = _model([(90, 120)], name="complete")
    result = _run(
        tmp_path,
        incomplete,
        _sequence(260, (160, "ATG"), (208, "TAA")),
        [_evidence("bridging", (80, 230))],
        annotation=_annotation(incomplete, complete),
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_start0 == 160
    assert result.summary.refined_sides_downstream_start_fallback == 1
    assert _summary_counts(result, "encountered_rejection_event_counts") == {
        "frozen_final_gene_structural_conflict": 1
    }
    assert result.summary.final_failed_sides == 0
    assert result.audit_rows[0]["result"] == "refined"
    assert result.audit_rows[0]["failure_reason"] == "."
    assert json.loads(result.audit_rows[0]["encountered_rejection_counts"]) == {
        "frozen_final_gene_structural_conflict": 1
    }
    assert result.summary.refined_sides_with_encountered_rejection_events == 1
    _assert_summary_reconciles_with_audit(result)


def test_multiple_rejected_structural_paths_are_counted_as_events_not_failures(
    tmp_path: Path,
) -> None:
    incomplete = _model([(151, 211)], missing_start=True, name="incomplete")
    complete = _model([(90, 120)], name="complete")
    result = _run(
        tmp_path,
        incomplete,
        _sequence(260, (160, "ATG"), (208, "TAA")),
        [
            _evidence("bridging_a", (80, 230)),
            _evidence("bridging_b", (70, 220)),
        ],
        annotation=_annotation(incomplete, complete),
    )
    assert result.annotation.genes[0].transcripts[0].source_model.cds_start0 == 160
    assert result.summary.final_refined_sides == 1
    assert result.summary.final_failed_sides == 0
    assert _summary_counts(result, "encountered_rejection_event_counts") == {
        "frozen_final_gene_structural_conflict": 2
    }
    assert _summary_counts(result, "encountered_rejection_side_counts") == {
        "frozen_final_gene_structural_conflict": 1
    }
    assert _summary_counts(result, "encountered_rejection_transcript_counts") == {
        "frozen_final_gene_structural_conflict": 1
    }
    _assert_summary_reconciles_with_audit(result)


def test_integrity_rejection_can_be_encountered_before_other_side_succeeds(
    tmp_path: Path,
) -> None:
    model = _model(
        [(154, 160)], missing_start=True, missing_stop=True, name="both"
    )
    result = _run(
        tmp_path,
        model,
        _sequence(240, (99, "A"), (152, "TG"), (160, "TAA")),
        [_evidence("spliced_both", (90, 100), (152, 200))],
    )
    rows = {row["side"]: row for row in result.audit_rows}
    assert rows["five_prime"]["failure_reason"] == "new_noncanonical_cds_splice"
    assert rows["three_prime"]["result"] == "refined"
    assert json.loads(rows["three_prime"]["encountered_rejection_counts"])[
        "new_noncanonical_cds_splice"
    ] == 1
    assert _summary_counts(result, "failure_counts") == {
        "new_noncanonical_cds_splice": 1
    }
    assert _summary_counts(result, "encountered_rejection_event_counts") == {
        "new_noncanonical_cds_splice": 3
    }
    assert _summary_counts(result, "encountered_rejection_side_counts") == {
        "new_noncanonical_cds_splice": 2
    }
    assert _summary_counts(result, "encountered_rejection_transcript_counts") == {
        "new_noncanonical_cds_splice": 1
    }
    assert result.summary.refined_sides_with_encountered_rejection_events == 1
    _assert_summary_reconciles_with_audit(result)


def test_one_sided_single_exon_evidence_is_not_a_complete_gene_utr_bridge(
    tmp_path: Path,
) -> None:
    incomplete = _model([(151, 160)], missing_start=True, name="incomplete")
    complete = _model([(120, 140)], name="complete")
    result = _run(
        tmp_path,
        incomplete,
        _sequence(220, (148, "ATG"), (157, "TAA")),
        [_evidence("one_sided", (120, 160))],
        annotation=_annotation(incomplete, complete),
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(148, 160)]
    assert result.summary.refined_transcripts == 1


@pytest.mark.parametrize(("five_limit", "three_limit"), [(20, 50), (50, 20)])
def test_missing_both_refinement_is_independent_of_residual_utr_limits(
    tmp_path: Path,
    five_limit: int,
    three_limit: int,
) -> None:
    model = _model([(100, 160)], missing_start=True, missing_stop=True)
    model.child_records.extend([
        _record("start_codon", 100, 103, "0"),
        _record("stop_codon", 157, 160, "0"),
    ])
    limits = UTRLearnedLimits(
        five_prime_exonic_bp=five_limit,
        three_prime_exonic_bp=three_limit,
        five_prime_genomic_bp=None,
        three_prime_genomic_bp=None,
        five_prime_exonic_source="test",
        three_prime_exonic_source="test",
        five_prime_genomic_source="unavailable",
        three_prime_genomic_source="unavailable",
    )
    result = _run(
        tmp_path,
        model,
        _sequence(240, (70, "ATG"), (163, "TAA")),
        [_evidence("tx", (20, 200))],
        utr_settings=UTRSettings(),
        limits=limits,
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(70, 166)]
    assert {row["side"] for row in result.audit_rows if row["result"] == "refined"} == {
        "five_prime",
        "three_prime",
    }
    child_types = {
        record.feature_type
        for record in refined.child_records
        if record.feature_type in {"start_codon", "stop_codon"}
    }
    assert child_types == {"start_codon", "stop_codon"}


def test_refined_children_are_rebuilt_without_stale_codon_or_overlapping_utr(
    tmp_path: Path,
) -> None:
    model = _model(
        [(100, 160)], missing_start=True, missing_stop=True
    )
    model.child_records.extend([
        _record("five_prime_UTR", 50, 100),
        _record("three_prime_UTR", 160, 190),
        _record("start_codon", 100, 103, "0"),
        _record("stop_codon", 157, 160, "0"),
        _record("selenocysteine", 130, 133, "0"),
        _record("intron", 80, 90),
        _record("transcription_start_site", 50, 51),
        _record("transcription_end_site", 189, 190),
        _record("source_translation_hint", 120, 123),
    ])
    annotation = _annotation(model)
    result = _run(
        tmp_path,
        model,
        _sequence(240, (70, "ATG"), (163, "TAA")),
        [_evidence("tx", (50, 190))],
        annotation=annotation,
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    by_type = {}
    for record in refined.child_records:
        by_type.setdefault(record.feature_type, []).append(
            (record.start0, record.end0)
        )
    assert by_type["five_prime_UTR"] == [(50, 70)]
    assert by_type["three_prime_UTR"] == [(166, 190)]
    assert by_type["start_codon"] == [(70, 73)]
    assert by_type["stop_codon"] == [(163, 166)]
    assert "selenocysteine" not in by_type
    assert "intron" not in by_type
    assert "transcription_start_site" not in by_type
    assert "transcription_end_site" not in by_type
    assert "source_translation_hint" not in by_type
    output = tmp_path / "rebuilt.gff3"
    write_final_annotation_gff3(
        result.annotation,
        output,
        FinalAnnotationOutputMetadata(False, False, False),
    )
    text = output.read_text()
    assert "\tfive_prime_UTR\t51\t70\t" in text
    assert "\tthree_prime_UTR\t167\t190\t" in text
    assert "\tstart_codon\t71\t73\t" in text
    assert "\tstop_codon\t164\t166\t" in text
    assert "\tstart_codon\t101\t103\t" not in text
    assert "\tintron\t" not in text
    assert "\ttranscription_start_site\t" not in text
    assert "\ttranscription_end_site\t" not in text


def test_emitted_integrity_is_current_and_selected_status_remains_provenance(
    tmp_path: Path,
) -> None:
    model = _model([(100, 160)], missing_start=True)
    result = _run(
        tmp_path,
        model,
        _sequence(220, (70, "ATG"), (157, "TAA")),
        [_evidence("tx", (50, 160))],
    )
    output = tmp_path / "integrity.gff3"
    write_final_annotation_gff3(
        result.annotation,
        output,
        FinalAnnotationOutputMetadata(False, False, False),
    )
    text = output.read_text()
    assert "sequence_integrity_status=evaluated_no_concern" in text
    assert "selected_sequence_integrity_status=evaluated_with_concern" in text
    assert result.audit_rows[0]["original_integrity_status"] == "evaluated_with_concern"
    assert result.audit_rows[0]["refined_integrity_status"] == "evaluated_no_concern"


@pytest.mark.parametrize("phase", [1, 2])
def test_upstream_completion_respects_nonzero_initial_phase_and_rebuilds_phase_zero(
    tmp_path: Path, phase: int
) -> None:
    end0 = 161 if phase == 1 else 162
    start_at = 71 if phase == 1 else 72
    stop_at = end0 - 3
    model = _model(
        [(100, end0)], missing_start=True, phases=[phase]
    )
    result = _run(
        tmp_path,
        model,
        _sequence(240, (start_at, "ATG"), (stop_at, "TAA")),
        [_evidence("tx", (50, end0))],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_start0 == start_at
    cds = [record for record in refined.child_records if record.feature_type == "CDS"]
    assert cds[0].phase == "0"


def test_minus_strand_upstream_completion_respects_nonzero_initial_phase(
    tmp_path: Path,
) -> None:
    model = _model(
        [(100, 161)], missing_start=True, strand="-", phases=[1]
    )
    result = _run(
        tmp_path,
        model,
        _sequence(240, (100, "TTA"), (187, "CAT")),
        [_evidence("tx", (100, 190), strand="-")],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(100, 190)]
    assert _first_cds_phases(refined) == ["0"]
    assert refined.sequence_integrity.start_codon_status == "present"


def test_minus_strand_source_utrs_are_clipped_away_from_refined_cds(
    tmp_path: Path,
) -> None:
    model = _model([(100, 160)], missing_start=True, strand="-")
    model.child_records.extend([
        _record("five_prime_UTR", 160, 200, strand="-"),
        _record("three_prime_UTR", 70, 100, strand="-"),
    ])
    result = _run(
        tmp_path,
        model,
        _sequence(240, (100, "TTA"), (187, "CAT")),
        [_evidence("tx", (70, 200), strand="-")],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    utrs = {
        record.feature_type: (record.start0, record.end0)
        for record in refined.child_records
        if record.feature_type in {"five_prime_UTR", "three_prime_UTR"}
    }
    assert utrs == {
        "five_prime_UTR": (190, 200),
        "three_prime_UTR": (70, 100),
    }


def test_nonzero_phase_splice_start_and_terminal_stop_continuation(
    tmp_path: Path,
) -> None:
    start_model = _model(
        [(154, 161)], missing_start=True, phases=[1], name="phase_start"
    )
    start_result = _run(
        tmp_path,
        start_model,
        _sequence(
            260,
            (99, "AT"),
            (101, "G"),
            (102, "GT"),
            (150, "AG"),
            (152, "G"),
            (158, "TAA"),
        ),
        [_evidence("start_path", (90, 102), (152, 161))],
    )
    start_refined = start_result.annotation.genes[0].transcripts[0].source_model
    assert start_refined.sequence_integrity.phase_status == "consistent"
    assert start_refined.cds_parts == [(99, 102), (152, 161)]
    assert _first_cds_phases(start_refined) == ["0", "0"]

    stop_model = _model(
        [(100, 111), (150, 158)],
        missing_stop=True,
        phases=[1, 2],
        name="phase_stop",
    )
    stop_result = _run(
        tmp_path,
        stop_model,
        _sequence(
            260,
            (101, "ATG"),
            (111, "GT"),
            (148, "AG"),
            (161, "TAA"),
        ),
        [_evidence("stop_path", (100, 111), (150, 200))],
    )
    stop_refined = stop_result.annotation.genes[0].transcripts[0].source_model
    assert stop_refined.cds_end0 == 164
    assert _first_cds_phases(stop_refined) == ["1", "2"]
    output = tmp_path / "nonzero_phase.gff3"
    write_final_annotation_gff3(
        stop_result.annotation,
        output,
        FinalAnnotationOutputMetadata(False, False, False),
    )
    cds_lines = [
        line.split("\t")
        for line in output.read_text().splitlines()
        if "\tCDS\t" in line
    ]
    assert [fields[7] for fields in cds_lines] == ["1", "2"]
    cds_only_output = tmp_path / "nonzero_phase.cds_only.gff3"
    write_final_annotation_gff3(
        stop_result.annotation,
        cds_only_output,
        FinalAnnotationOutputMetadata(False, False, False),
        cds_only_structure=True,
    )
    cds_only_lines = [
        line.split("\t")
        for line in cds_only_output.read_text().splitlines()
        if "\tCDS\t" in line
    ]
    assert [fields[7] for fields in cds_only_lines] == ["1", "2"]


def test_missing_stop_continues_a_non_codon_aligned_terminal_frame(
    tmp_path: Path,
) -> None:
    model = _model(
        [(100, 159)], missing_stop=True, phases=[1], name="partial_terminal"
    )
    result = _run(
        tmp_path,
        model,
        _sequence(240, (101, "ATG"), (158, "TAA")),
        [_evidence("stop_path", (100, 200))],
    )
    refined = result.annotation.genes[0].transcripts[0].source_model
    assert refined.cds_parts == [(100, 161)]
    assert refined.sequence_integrity.stop_codon_status == "included_in_cds"
    assert _first_cds_phases(refined) == ["1"]


def _first_cds_phases(model: CodingModel) -> list[str]:
    phases = {
        (record.start0, record.end0): record.phase
        for record in model.child_records
        if record.feature_type == "CDS"
    }
    return [
        phases[part]
        for part in sorted(model.cds_parts, reverse=model.strand == "-")
    ]


def test_endpoint_choice_is_stable_across_python_hash_seeds() -> None:
    script = r'''
from moselecta.utr_augmentation import UTRCandidateSide, UTREvidence, UTRSettings, UTRStructure, _select_endpoint

structure = UTRStructure("same")
for endpoint, evidence_id, members in (
    (10, "a", {"alpha", "beta", "gamma"}),
    (20, "b", {"alpha", "beta", "delta"}),
):
    evidence = UTREvidence(
        evidence_id=evidence_id,
        kind="transcript",
        chrom="chr1",
        strand="+",
        exon_parts=((endpoint, 100),),
        source_member_counts=tuple((source, 1) for source in members),
    )
    structure.add(UTRCandidateSide(
        gene_id="g", transcript_id="t", target=None, evidence=evidence,
        compatibility_class="test", side="five_prime", endpoint0=endpoint,
        extension_interval=(endpoint, 100), added_exons=((endpoint, 100),),
        added_introns=(), complete_contributing_exon_lengths=(100-endpoint,),
        extension_exonic_bp=100-endpoint, extension_genomic_bp=100-endpoint,
        structure_key="same", source_cds_class=".", source_cds_overlap_bp=0,
    ))
print(_select_endpoint(structure, UTRSettings()).endpoint0)
'''
    outputs = []
    for seed in ("1", "2", "17", "101"):
        environment = dict(os.environ)
        environment["PYTHONHASHSEED"] = seed
        environment["PYTHONPATH"] = str(Path(__file__).parents[1] / "src")
        outputs.append(
            subprocess.check_output(
                [sys.executable, "-c", script], env=environment, text=True
            ).strip()
        )
    assert len(set(outputs)) == 1


def test_cli_refines_without_add_utrs_and_opt_out_preserves_original(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    model_gff = tmp_path / "model.gff3"
    model_gff.write_text(
        "##gff-version 3\n"
        "chr1\tmodel\tgene\t101\t160\t.\t+\t.\tID=g1\n"
        "chr1\tmodel\tmRNA\t101\t160\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tmodel\texon\t101\t160\t.\t+\t.\tParent=t1\n"
        "chr1\tmodel\tCDS\t101\t160\t.\t+\t0\tParent=t1\n"
    )
    transcript_gff = tmp_path / "transcripts.gff3"
    transcript_gff.write_text(
        "##gff-version 3\n"
        "chr1\trna\ttranscript\t51\t180\t.\t+\t.\tID=rna1\n"
        "chr1\trna\texon\t51\t180\t.\t+\t.\tParent=rna1\n"
    )
    fasta = _write_fasta(
        tmp_path,
        _sequence(220, (70, "ATG"), (157, "TAA")),
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(root / "src")

    def run(name: str, *extra: str) -> Path:
        prefix = tmp_path / name
        subprocess.run(
            [
                sys.executable,
                "-m",
                "moselecta",
                "--inputs",
                str(model_gff),
                "--source-name",
                "Model",
                "--evidence-transcripts",
                str(transcript_gff),
                "--evidence-transcript-source-name",
                "RNA",
                "--genome-fasta",
                str(fasta),
                "--out-prefix",
                str(prefix),
                "--pruning-level",
                "none",
                "--no-alternative-splice",
                "--output-profile",
                "standard",
                "--progress-log-level",
                "none",
                *extra,
            ],
            cwd=root,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        )
        return prefix

    refined_prefix = run("refined")
    ordinary = Path(str(refined_prefix) + ".final_annotation.gff3").read_text()
    assert "\tCDS\t71\t160\t.\t+\t0\t" in ordinary
    assert "five_prime_UTR" not in ordinary and "three_prime_UTR" not in ordinary
    assert Path(str(refined_prefix) + ".cds_terminus_refinement.audit.tsv").exists()
    performed_summary = Path(
        str(refined_prefix) + ".cds_terminus_refinement.summary.tsv"
    ).read_text()
    assert "performed" in performed_summary

    disabled_prefix = run("disabled", "--no-refine-cds-termini")
    disabled = Path(str(disabled_prefix) + ".final_annotation.gff3").read_text()
    assert "\tCDS\t101\t160\t.\t+\t0\t" in disabled
    disabled_summary = Path(
        str(disabled_prefix) + ".cds_terminus_refinement.summary.tsv"
    ).read_text()
    assert "disabled" in disabled_summary

    utr_prefix = run("with_utrs", "--add-utrs")
    augmented = Path(str(utr_prefix) + ".final_annotation.with_utrs.gff3").read_text()
    assert "\tCDS\t71\t160\t.\t+\t0\t" in augmented
    assert "\tfive_prime_UTR\t51\t70\t" in augmented

    rejected_utr_prefix = run(
        "with_rejected_utr",
        "--add-utrs",
        "--utr-five-prime-exonic-override-bp",
        "10",
    )
    ordinary_with_rejected_utr = Path(
        str(rejected_utr_prefix) + ".final_annotation.gff3"
    ).read_text()
    augmented_with_rejected_utr = Path(
        str(rejected_utr_prefix) + ".final_annotation.with_utrs.gff3"
    ).read_text()
    def cds_lines(text: str) -> list[str]:
        return [line for line in text.splitlines() if "\tCDS\t" in line]
    assert cds_lines(ordinary_with_rejected_utr) == cds_lines(ordinary)
    assert "\tCDS\t71\t160\t.\t+\t0\t" in augmented_with_rejected_utr
    assert "\tfive_prime_UTR\t" not in augmented_with_rejected_utr


def test_zero_incomplete_run_skips_deferred_catalogue_and_writes_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    model_gff = tmp_path / "complete.gff3"
    model_gff.write_text(
        "##gff-version 3\n"
        "chr1\tmodel\tgene\t71\t160\t.\t+\t.\tID=g1\n"
        "chr1\tmodel\tmRNA\t71\t160\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tmodel\texon\t71\t160\t.\t+\t.\tParent=t1\n"
        "chr1\tmodel\tCDS\t71\t160\t.\t+\t0\tParent=t1\n"
    )
    transcript_gff = tmp_path / "transcripts.gff3"
    transcript_gff.write_text(
        "##gff-version 3\n"
        "chr1\trna\ttranscript\t51\t180\t.\t+\t.\tID=rna1\n"
        "chr1\trna\texon\t51\t180\t.\t+\t.\tParent=rna1\n"
    )
    fasta = _write_fasta(
        tmp_path, _sequence(220, (70, "ATG"), (157, "TAA"))
    )

    def unexpected_catalogue(*_args, **_kwargs):
        raise AssertionError("completion catalogue must not be built")

    monkeypatch.setattr(
        pipeline_module, "build_transcript_utr_catalogue", unexpected_catalogue
    )
    prefix = tmp_path / "complete"
    exit_code = pipeline_module.main([
        "--inputs", str(model_gff),
        "--source-name", "Model",
        "--evidence-transcripts", str(transcript_gff),
        "--evidence-transcript-source-name", "RNA",
        "--genome-fasta", str(fasta),
        "--out-prefix", str(prefix),
        "--pruning-level", "none",
        "--no-alternative-splice",
        "--output-profile", "standard",
        "--progress-log-level", "none",
    ])
    assert exit_code == 0
    summary = Path(
        str(prefix) + ".cds_terminus_refinement.summary.tsv"
    ).read_text()
    assert "execution_status" in summary
    assert "zero_incomplete_models" in summary


def test_summary_is_written_when_genome_prerequisite_is_unavailable(
    tmp_path: Path,
) -> None:
    model_gff = tmp_path / "model.gff3"
    model_gff.write_text(
        "##gff-version 3\n"
        "chr1\tmodel\tgene\t101\t160\t.\t+\t.\tID=g1\n"
        "chr1\tmodel\tmRNA\t101\t160\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tmodel\texon\t101\t160\t.\t+\t.\tParent=t1\n"
        "chr1\tmodel\tCDS\t101\t160\t.\t+\t0\tParent=t1\n"
    )
    prefix = tmp_path / "no_genome"
    exit_code = pipeline_module.main([
        "--inputs", str(model_gff),
        "--source-name", "Model",
        "--out-prefix", str(prefix),
        "--pruning-level", "none",
        "--no-alternative-splice",
        "--output-profile", "standard",
        "--progress-log-level", "none",
    ])
    assert exit_code == 0
    summary = Path(
        str(prefix) + ".cds_terminus_refinement.summary.tsv"
    ).read_text()
    assert "prerequisite_unavailable:no_genome_fasta" in summary


def test_summary_reports_not_applicable_without_final_annotation(
    tmp_path: Path,
) -> None:
    model_gff = tmp_path / "model.gff3"
    model_gff.write_text(
        "##gff-version 3\n"
        "chr1\tmodel\tgene\t71\t160\t.\t+\t.\tID=g1\n"
        "chr1\tmodel\tmRNA\t71\t160\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tmodel\texon\t71\t160\t.\t+\t.\tParent=t1\n"
        "chr1\tmodel\tCDS\t71\t160\t.\t+\t0\tParent=t1\n"
    )
    fasta = _write_fasta(
        tmp_path, _sequence(220, (70, "ATG"), (157, "TAA"))
    )
    prefix = tmp_path / "no_final"
    exit_code = pipeline_module.main([
        "--inputs", str(model_gff),
        "--source-name", "Model",
        "--genome-fasta", str(fasta),
        "--out-prefix", str(prefix),
        "--stop-after", "scored-candidates",
        "--pruning-level", "none",
        "--no-alternative-splice",
        "--output-profile", "standard",
        "--progress-log-level", "none",
    ])
    assert exit_code == 0
    summary = Path(
        str(prefix) + ".cds_terminus_refinement.summary.tsv"
    ).read_text()
    assert "not_applicable:no_final_annotation" in summary
