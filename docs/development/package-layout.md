# MoSelecta package layout

> Maintainer reference for the v0.3.0 release. Validation and
> distribution checks are in [Testing and regression](testing-and-regression.md).

MoSelecta uses a `src/` package layout. The major modules are:

```text
src/moselecta/
  auto_pruning.py reference-free pruning scenarios, diagnostics, safeguards, and decisions
  cli.py          CLI entry point; calls the pipeline main function
  compare.py      tolerance-aware output comparison helper
  config.py       argparse construction and runtime option resolution
  console_reporting.py shared deterministic command-line report formatting
  constants.py    version, presets, output-key constants, score weights
  core.py         compatibility wrapper for older imports
  gffio.py        GFF/GTF parsing, attribute handling, model extraction
  genome_sequence.py indexed FASTA access and genome-derived model integrity metrics
  intervals.py    generic interval helpers plus the evidence-locus overlap index
  loci.py         superlocus/sublocus construction and source-agreement helpers
  log_utils.py    progress logging and timing utilities
  models.py       dataclasses for records, models, loci, metrics, and reports
  outputs.py      GFF3/TSV output writers and formatting helpers
  pipeline.py     command-line pipeline orchestration
  pipeline_stages.py typed final-stage and alternative-splice orchestration boundaries
  protein_evidence.py    protein evidence parsing, collapse, scoring, and support
  run_context.py  run-level state shared across pipeline stages
  scoring.py      generic score and integrated model-selection-ranking helpers
  selection.py    monosublocus de-overlap, backfill, evidence rescue, and final pruning
  structural.py   structural support and pairwise/consensus metrics
  structural_coherence.py candidate bridge-region detection, source-state voting, and coherence audits
  transcript_evidence.py transcript parsing, exact-chain collapse, scoring, and cDNA context
  tsvio.py        TSV table writing helpers
  utr_augmentation.py post-selection UTR catalogues, selection, reconstruction, and audit writers
  validation.py   intrinsic model metrics and validation helpers
```

Current design aim:

- `pipeline.py` coordinates stages.
- `pipeline_stages.py` resolves optional selector execution, freezes retained
  primaries, and coordinates candidate facts through selection and immutable
  annotation augmentation without implementing biological rules or writing
  files.
- `console_reporting.py` owns reusable console formatting and source-score
  summary calculations; `pipeline.py` retains report ordering and stage timing.
- Biological/scoring logic lives in focused modules.
- `core.py` remains only for compatibility.
- New evidence/genome functionality should be added as new modules rather than folded into `pipeline.py`.

The dependency direction is core models/configuration, then biological stage
modules, then annotation/report writers, then orchestration. Lower-level
biological and writer modules must not import `pipeline.py`.

UTR augmentation is isolated in `utr_augmentation.py`. `pipeline.py` builds
the compact endpoint-paired transcript catalogue before raw transcript objects
are released, then calls the augmentation module only after alternative-splice
selection has finalized the public transcript population. The module consumes
the immutable `FinalAnnotation` and all normalized prediction models, and
returns a separate augmented representation; it never mutates `CodingModel` or
the ordinary public annotation.

Possible future modules include shared evidence abstractions, BED12/junction loaders, automatic-pruning diagnostics, and deterministic parallel workers. Protein and transcript evidence remain separate focused modules until their common behaviour is stable enough to factor safely.

The current design keeps cDNA-context classification in `scoring.py`, transcript construction/context measurement in `transcript_evidence.py`, pruning decisions in `selection.py`, and orchestration/option reporting in `pipeline.py`. This separation keeps missing context neutral and avoids coupling transcript construction to final-pruning policy.

The current design keeps unsupported-short-exon feature comparison in `structural.py`, its aggregate negative term in `scoring.py`, option validation/reporting in `pipeline.py`, and GFF/TSV presentation in `outputs.py`. The feature remains separate from intrinsic metric flags so it cannot independently trigger intrinsic singleton pruning.

Candidate-level structural coherence is isolated in `structural_coherence.py`. It resolves source-gene separation and continuity routes, applies same-gene and bridge-hypothesis grouping, and emits the auditable region rows consumed by `scoring.py`. Candidate-anchor overlap values are cached per superlocus; this cache is a performance optimization and does not alter source votes, penalties, or output fields.


The current design keeps the generic `IntervalOverlapIndex` in `intervals.py` because it is shared by protein and transcript evidence. Evidence-specific association and scoring remain in their existing modules. High-count records use fixed slotted layouts; candidate and collapsed-evidence records retain only compact repeatedly used geometry, per-model evidence feature caches are temporary, and raw evidence is released after collapse. No parallel-worker or serialization layer has been introduced yet.


The current genome-sequence design is isolated in `genome_sequence.py`: it builds/reads a standard `.fai`, fetches one candidate-superlocus window at a time, reconstructs strand-aware CDS sequence, validates phase/start/stop/splice properties, attaches compact slotted metrics to each candidate model, and builds all-source/per-source summaries without retaining a second model-sized row table. `pipeline.py` controls optional activation and output routing, `scoring.py` applies the negative term to the existing model-only aggregate, and `selection.py` applies the explicit sequence-integrity pruning whitelist and rescue/protection rules. This keeps absence of a FASTA neutral and leaves the superlocus as a natural future parallel work unit. The `very-high` preset is resolved in `config.py`; it does not duplicate pruning logic. `selection.py` implements sequence-integrity `scope=all` as a superset of backfill by retaining the backfill-standard route and adding the global-severe route, with route-specific audit fields returned to `pipeline.py` for reporting.
