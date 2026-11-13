# Output files and interpretation

> Authoritative output reference for MoSelecta `0.2.16`.
> Paths append the shown suffix to `--out-prefix` unless overridden.

## Which annotation should I use?

`<prefix>.final_annotation.gff3` is the normal integrated final annotation. It
contains frozen primary transcripts and, unless `--no-alternative-splice` is
used, accepted alternatives. `<prefix>.final_annotation.with_utrs.gff3` is a
separate optional product written by `--add-utrs`; it never replaces the
ordinary final GFF3.

`<prefix>.standardised_primary_models.gff3` contains only frozen `.1`
primaries. Files named `selected`, `final_selected`, `pruned`, or `review` are
stage/audit products, not preferred substitutes for the final annotation.

## Output profiles

| Profile | Intended use | Membership |
|---|---|---|
| `minimal` | Compact final run | Final annotation, configuration, input/validation/exclusion and score/source/evidence/sequence summaries, plus enabled correction, refinement, UTR and automatic-pruning summaries. |
| `standard` | Default auditable run | Minimal plus primary-only annotation, selected/suppressed/backfilled/pruned ledgers, removed/shadow ledgers, and action/audit files for enabled corrections. |
| `full` | Scientific investigation | All 85 output keys except two feature-gated strict-replacement reports; those are added when active. |
| `debug` | Complete debugging | Same output-key membership as full with debug execution context. |

Feature-specific products may be omitted or automatically enabled according to
inputs and stages. `--write-output` and `--skip-output` amend a profile.

## Principal products and ledgers

| Suffix | Purpose |
|---|---|
| `.final_annotation.gff3` | Normal final integrated annotation. |
| `.standardised_primary_models.gff3` | Frozen retained primaries only. |
| `.final_annotation.with_utrs.gff3` | Separate optional UTR-augmented annotation. |
| `.utr_junction_validation.summary.tsv` | Bounded aggregate UTR validation counts by evidence origin; present in every profile. |
| `.utr_junction_validation.audit.tsv` | Exception-only inherited and candidate UTR validation rows; present in full/debug or when explicitly requested. |
| `.run_configuration.json` | Requested/resolved options, inputs and run metadata. |
| `.input_summary.tsv` | Input source/file counts and parsing summary. |
| `.excluded_models.tsv` | Models excluded during parsing or eligibility. |
| `.validation_failures.tsv` | Model/superlocus validation findings. |
| `.model_selection_scores.tsv` | Per-model components, aggregate, scaled score, ranks and ties. |
| `.model_selection_score_summary.tsv` | Source-level descriptive score summaries. |
| `.selected_models.tsv` | Monosublocus-selected models before backfill. |
| `.suppressed_monosublocus_models.tsv` | Winners suppressed by overlap conflict. |
| `.backfilled_models.tsv` | Models restored by ordinary or strict backfill. |
| `.backfill_rejected_models.tsv` | Backfill candidates and blockers/reasons. |
| `.r1_backfill_replacements.tsv` | Atomic strict score-improving replacements. |
| `.final_selected_models.pruned.tsv` | Selected models retained after pruning. |
| `.final_selected_models.pruned.removed.tsv` | Models removed by final pruning with reasons/protection fields. |
| `.pruning_scenario_summary.tsv` | Simulated fixed pruning scenarios. |
| `.auto_pruning_diagnostics.tsv` | Automatic-pruning inputs and decision. |
| `.alternative_splice_candidates.tsv` | Candidate-to-primary eligibility and decision ledger. |
| `.alternative_splice_genes.tsv` | Per-gene alternative outcome. |
| `.alternative_splice_run_summary.tsv` | Run-level selector counts and settings. |

## Scientific diagnostics

Structural products include `.model_structural_support.tsv`,
`.structural_corroboration.tsv`, `.structural_coherence_regions.tsv`,
`.structural_features.tsv`, and the superlocus/sublocus/relationship tables.
Protein and transcript evidence each provide a summary, raw-to-collapsed map,
collapsed records, loci/members, and per-model support. Genome assessment
provides `.sequence_integrity_summary.tsv` and
`.sequence_integrity_model_metrics.tsv`.

Early correction/exclusion families use `.summary.tsv`, `.actions.tsv` and/or
`.audit.tsv` for terminal-stop canonicalisation, 5-prime CDS boundary
canonicalisation, phase recovery, definite internal-stop exclusion, and extreme
short-CDS exclusion. CDS refinement and UTR augmentation have summary/audit
products. These ledgers are not retained annotations.

The compact UTR validation summary reports assessed, accepted, rejected and
unevaluable sides, junction and motif totals, side-specific rejection counts,
and deterministic reason counts. The detailed UTR junction audit contains
exceptions only: definite rejection rows and unusual unevaluable cases. It
does not repeat accepted sides, passing junctions or routine no-FASTA motif
states. The audit separates side-level reasons from reasons attached to a
specific intron and reports the final side action.

<!-- GENERATED OUTPUT INDEX START -->

## Complete output-key and default-suffix index

This index is generated from `ALL_OUTPUT_KEYS`, `output_keys_for_profile()` and the current `with_default()` paths in `pipeline.py`. A check mark means membership in the base profile before feature/stage forcing or explicit write/skip controls.

| Output key | Default suffix | minimal | standard | full | debug |
|---|---|:---:|:---:|:---:|:---:|
| `alternative_splice_candidate_audit_gff3` | `.alternative_splice_candidate_audit.gff3` | — | — | ✓ | ✓ |
| `alternative_splice_candidates_tsv` | `.alternative_splice_candidates.tsv` | — | — | ✓ | ✓ |
| `alternative_splice_genes_tsv` | `.alternative_splice_genes.tsv` | — | — | ✓ | ✓ |
| `alternative_splice_run_summary_tsv` | `.alternative_splice_run_summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `auto_pruning_diagnostics_tsv` | `.auto_pruning_diagnostics.tsv` | ✓ | ✓ | ✓ | ✓ |
| `backfill_candidates_tsv` | `.backfill_candidates.tsv` | — | — | ✓ | ✓ |
| `backfill_rejected_models_tsv` | `.backfill_rejected_models.tsv` | — | — | ✓ | ✓ |
| `backfilled_models_tsv` | `.backfilled_models.tsv` | — | ✓ | ✓ | ✓ |
| `cds_5prime_boundary_canonicalisation_actions_tsv` | `.cds_5prime_boundary_canonicalisation.actions.tsv` | — | ✓ | ✓ | ✓ |
| `cds_5prime_boundary_canonicalisation_audit_tsv` | `.cds_5prime_boundary_canonicalisation.audit.tsv` | — | — | ✓ | ✓ |
| `cds_5prime_boundary_canonicalisation_summary_tsv` | `.cds_5prime_boundary_canonicalisation.summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `cds_phase_recovery_actions_tsv` | `.cds_phase_recovery.actions.tsv` | — | ✓ | ✓ | ✓ |
| `cds_phase_recovery_audit_tsv` | `.cds_phase_recovery.audit.tsv` | — | — | ✓ | ✓ |
| `cds_phase_recovery_summary_tsv` | `.cds_phase_recovery.summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `cds_terminus_refinement_audit_tsv` | `.cds_terminus_refinement.audit.tsv` | — | ✓ | ✓ | ✓ |
| `cds_terminus_refinement_summary_tsv` | `.cds_terminus_refinement.summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `coding_models_tsv` | `.coding_models.tsv` | — | — | ✓ | ✓ |
| `combined_gff3` | `.coding_models_with_superloci.gff3` | — | — | ✓ | ✓ |
| `definite_internal_stop_exclusion_summary_tsv` | `.definite_internal_stop_exclusion.summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `definite_internal_stop_exclusions_tsv` | `.definite_internal_stop_exclusions.tsv` | — | ✓ | ✓ | ✓ |
| `excluded_models_tsv` | `.excluded_models.tsv` | ✓ | ✓ | ✓ | ✓ |
| `extreme_short_cds_exclusion_summary_tsv` | `.extreme_short_cds_exclusion.summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `extreme_short_cds_exclusions_tsv` | `.extreme_short_cds_exclusions.tsv` | — | ✓ | ✓ | ✓ |
| `final_annotation_gff3` | `.final_annotation.gff3` | ✓ | ✓ | ✓ | ✓ |
| `final_pruned_selected_models_gff3` | `.final_selected_models.pruned.gff3` | — | — | ✓ | ✓ |
| `final_pruned_selected_models_tsv` | `.final_selected_models.pruned.tsv` | — | ✓ | ✓ | ✓ |
| `final_pruned_selected_removed_tsv` | `.final_selected_models.pruned.removed.tsv` | — | ✓ | ✓ | ✓ |
| `final_selected_models_gff3` | `.final_selected_models.gff3` | — | — | ✓ | ✓ |
| `final_selected_models_tsv` | `.final_selected_models.tsv` | — | — | ✓ | ✓ |
| `final_shadow_gene_pairs_tsv` | `.final_shadow_gene_pairs.tsv` | — | ✓ | ✓ | ✓ |
| `input_summary_tsv` | `.input_summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `model_intrinsic_metrics_tsv` | `.model_intrinsic_metrics.tsv` | — | — | ✓ | ✓ |
| `model_selection_score_summary_tsv` | `.model_selection_score_summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `model_selection_scores_tsv` | `.model_selection_scores.tsv` | — | — | ✓ | ✓ |
| `model_structural_support_tsv` | `.model_structural_support.tsv` | — | — | ✓ | ✓ |
| `monosubloci_tsv` | `.monosubloci.tsv` | — | — | ✓ | ✓ |
| `monosublocus_edges_tsv` | `.monosublocus_edges.tsv` | — | — | ✓ | ✓ |
| `monosublocus_members_tsv` | `.monosublocus_members.tsv` | — | — | ✓ | ✓ |
| `monosublocus_review_gff3` | `.monosublocus_review.gff3` | — | — | ✓ | ✓ |
| `other_features_gff3` | `.other_or_no_cds_features.gff3` | — | — | ✓ | ✓ |
| `protein_evidence_collapsed_alignments_tsv` | `.protein_evidence_collapsed_alignments.tsv` | — | — | ✓ | ✓ |
| `protein_evidence_model_support_tsv` | `.protein_evidence_model_support.tsv` | — | — | ✓ | ✓ |
| `protein_evidence_raw_to_collapsed_tsv` | `.protein_evidence_raw_to_collapsed.tsv` | — | — | ✓ | ✓ |
| `protein_evidence_subloci_tsv` | `.protein_evidence_subloci.tsv` | — | — | ✓ | ✓ |
| `protein_evidence_sublocus_members_tsv` | `.protein_evidence_sublocus_members.tsv` | — | — | ✓ | ✓ |
| `protein_evidence_summary_tsv` | `.protein_evidence_summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `protein_evidence_superloci_tsv` | `.protein_evidence_superloci.tsv` | — | — | ✓ | ✓ |
| `pruning_scenario_summary_tsv` | `.pruning_scenario_summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `r1_backfill_replacements_tsv` | `.r1_backfill_replacements.tsv` | — | — | — | — |
| `r1_backfill_summary_tsv` | `.r1_backfill_summary.tsv` | — | — | — | — |
| `relationship_edges_tsv` | `.relationship_edges.tsv` | — | — | ✓ | ✓ |
| `run_configuration_json` | `.run_configuration.json` | ✓ | ✓ | ✓ | ✓ |
| `scored_pruned_gff3` | `.coding_models_with_superloci.no_singleton_subloci.gff3` | — | — | ✓ | ✓ |
| `selected_models_gff3` | `.selected_models.gff3` | — | — | ✓ | ✓ |
| `selected_models_tsv` | `.selected_models.tsv` | — | ✓ | ✓ | ✓ |
| `sequence_integrity_model_metrics_tsv` | `.sequence_integrity_model_metrics.tsv` | — | — | ✓ | ✓ |
| `sequence_integrity_summary_tsv` | `.sequence_integrity_summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `source_auto_weights_tsv` | `.source_auto_weights.tsv` | ✓ | ✓ | ✓ | ✓ |
| `standardised_primary_models_gff3` | `.standardised_primary_models.gff3` | — | ✓ | ✓ | ✓ |
| `structural_coherence_regions_tsv` | `.structural_coherence_regions.tsv` | — | — | ✓ | ✓ |
| `structural_corroboration_tsv` | `.structural_corroboration.tsv` | — | — | ✓ | ✓ |
| `structural_features_tsv` | `.structural_features.tsv` | — | — | ✓ | ✓ |
| `subloci_tsv` | `.subloci.tsv` | — | — | ✓ | ✓ |
| `sublocus_members_tsv` | `.sublocus_members.tsv` | — | — | ✓ | ✓ |
| `sublocus_winners_gff3` | `.sublocus_winners.gff3` | — | — | ✓ | ✓ |
| `sublocus_winners_tsv` | `.sublocus_winners.tsv` | — | — | ✓ | ✓ |
| `superloci_tsv` | `.superloci.tsv` | — | — | ✓ | ✓ |
| `superlocus_members_tsv` | `.superlocus_members.tsv` | — | — | ✓ | ✓ |
| `suppressed_monosublocus_models_tsv` | `.suppressed_monosublocus_models.tsv` | — | ✓ | ✓ | ✓ |
| `terminal_stop_canonicalisation_actions_tsv` | `.terminal_stop_canonicalisation.actions.tsv` | — | ✓ | ✓ | ✓ |
| `terminal_stop_canonicalisation_audit_tsv` | `.terminal_stop_canonicalisation.audit.tsv` | — | — | ✓ | ✓ |
| `terminal_stop_canonicalisation_summary_tsv` | `.terminal_stop_canonicalisation.summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `transcript_evidence_collapsed_transcripts_tsv` | `.transcript_evidence_collapsed_transcripts.tsv` | — | — | ✓ | ✓ |
| `transcript_evidence_model_support_tsv` | `.transcript_evidence_model_support.tsv` | — | — | ✓ | ✓ |
| `transcript_evidence_raw_to_collapsed_tsv` | `.transcript_evidence_raw_to_collapsed.tsv` | — | — | ✓ | ✓ |
| `transcript_evidence_subloci_tsv` | `.transcript_evidence_subloci.tsv` | — | — | ✓ | ✓ |
| `transcript_evidence_sublocus_members_tsv` | `.transcript_evidence_sublocus_members.tsv` | — | — | ✓ | ✓ |
| `transcript_evidence_summary_tsv` | `.transcript_evidence_summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `transcript_evidence_superloci_tsv` | `.transcript_evidence_superloci.tsv` | — | — | ✓ | ✓ |
| `utr_augmentation_audit_tsv` | `.utr_augmentation.audit.tsv` | — | ✓ | ✓ | ✓ |
| `utr_augmentation_summary_tsv` | `.utr_augmentation.summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `utr_augmented_final_annotation_gff3` | `.final_annotation.with_utrs.gff3` | ✓ | ✓ | ✓ | ✓ |
| `utr_junction_validation_audit_tsv` | `.utr_junction_validation.audit.tsv` | — | — | ✓ | ✓ |
| `utr_junction_validation_summary_tsv` | `.utr_junction_validation.summary.tsv` | ✓ | ✓ | ✓ | ✓ |
| `validation_failures_tsv` | `.validation_failures.tsv` | ✓ | ✓ | ✓ | ✓ |

The strict-replacement keys `r1_backfill_replacements_tsv` and `r1_backfill_summary_tsv` are deliberately absent from the base full/debug sets and are enabled when strict replacement is active. Other feature outputs can be forced by their stage or input even when absent from a base profile.

<!-- GENERATED OUTPUT INDEX END -->

## Important columns and conventions

TSVs have a header row. `.` means not applicable, unavailable, not evaluated,
or no value as defined by the column; it is not numeric zero. Boolean fields
normally use `true`/`false`. Common provenance includes source and input IDs,
internal/output model ID, locus IDs, action/reason/status, blocker/replaced IDs,
selection stage, protection/rescue route, and warnings.

Per-model score fields are:

- `model_structural_support_score`: signed model-only aggregate;
- `model_selection_score`: integrated signed aggregate used for ranking;
- `model_selection_score_scaled`: clipped positive reporting/pruning transform;
- `model_selection_score_rank`: epsilon-grouped competition rank within the scored canonical sublocus (`1, 1, 3` after a two-way top tie);
- `model_selection_order_rank`: unique `1..N` rank using the exact shared primary-selector total order within that same pre-selection population;
- tie-group and tie-break fields: the epsilon group and the exact criterion that first determines total order within it.

Order rank is assigned before later suppression, backfill and pruning. Rank 1 matches the first-pass sublocus winner but does not guarantee final retention, and ranks are not globally comparable across subloci.

`mean_scaled_score_all` and `mean_scaled_score_non_singletons` are source-level
means. The latter excludes canonical single-candidate subloci where
`sublocus_model_count == 1`. These are descriptive diagnostics, not validated
estimates of source accuracy.

## GFF3 attributes

Final records use stable public gene/transcript IDs, standard `ID`/`Parent`
relationships, source/provenance attributes, and retained selection/evidence
diagnostics. The normal annotation is CDS-based; the optional augmented file
adds selected UTR/exon structure. Diagnostic GFF3s may contain synthetic
top-level audit transcripts and must not be treated as retained genes.

## Model-state glossary

| State | Meaning and responsible stage |
|---|---|
| Candidate | Parsed CDS model that passed early eligibility. |
| Excluded | Removed before competitive selection by parsing, geometry, stop or validation rules. |
| Sublocus winner | Top-ranked candidate in one first-pass sublocus. |
| Selected | Winner retained after monosublocus conflict resolution. |
| Suppressed | Sublocus winner blocked by a selected conflicting winner. |
| Backfilled | Eligible model added after conflict resolution. |
| Replaced | Holder removed atomically by a strictly higher-scoring eligible candidate. |
| Pruned | Selected model removed by final pruning. |
| Retained | Selected model surviving final pruning. |
| Primary | Frozen retained model assigned `.1`. |
| Alternative transcript | Post-freeze accepted isoform assigned `.2+`. |

See [Scoring](technical/scoring.md),
[Selection and backfill](technical/selection-and-backfill.md), and
[Alternative transcripts](technical/alternative-transcripts.md).

