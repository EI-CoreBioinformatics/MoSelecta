# Changelog

## v0.3.0

### Initial public release

- Integrates multiple CDS-containing prediction annotations into one deterministic final annotation.
- Uses explicit structural, completeness, source, integrity, and optional protein/transcript evidence components in auditable model-selection scores.
- Supports source-collapsed run-aware structural corroboration and candidate structural-coherence penalties.
- Applies sublocus selection, monosublocus conflict resolution, ordinary backfill, and safeguarded score-improving replacement of eligible blockers.
- Provides fixed and reference-free automatic final-pruning modes with rescue/protection rules and detailed diagnostics.
- Freezes primary transcript identities before optional alternative-transcript selection.
- Supports optional genome-based sequence assessment and post-selection CDS terminus refinement.
- Supports separate evidence-based UTR augmentation without changing the ordinary CDS-based final annotation.
- Reports model-level `model_selection_score` fields and source-level descriptive means, including `mean_scaled_score_non_singletons`.
- Provides Basic, Full, and Expert command help through `--help`, `--full-help`, and `--expert-help`.
- Reorganises the documentation into practical user guidance, output/configuration references, focused scientific methods, and limited maintainer material.
