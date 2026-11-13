# MoSelecta v0.3.0 configuration reference

> Documentation for MoSelecta `0.2.16`. CLI names and defaults below are derived from the accepted parser. `moselecta --expert-help` is the complete machine-local authority.

## How to use this reference

The documentation chain is:

- the [user guide](user-guide/MoSelecta-v0.3.0-user-guide.pdf) explains what happens and why;
- this page explains what users can control and which stage is affected;
- specialist references define exact formulas and decision rules;
- [Outputs](outputs.md) identifies the resulting files and fields.

MoSelecta has no general configuration-file input. Values resolve from built-in defaults, then preset-derived settings, then explicit command-line overrides. Requested and effective values are recorded in `<prefix>.run_configuration.json`.

Help levels change presentation, not parseability:

- `moselecta --help`: routine inputs and workflow choices;
- `moselecta --full-help`: advanced, biologically interpretable controls;
- `moselecta --expert-help`: every supported option, including low-level weights, thresholds and output filenames.

Changing an option can alter an implemented calculation or decision without guaranteeing a particular biological model will be retained. Score weights, conflict geometry and pruning predicates act at different stages.

## Inputs, labels and output identity

`--inputs` accepts repeated, comma-separated or `@FILE` GFF/GTF paths. Each candidate annotation needs usable CDS features. `--source-name` supplies one label per expanded input; labels identify prediction sources throughout grouping, scoring and reporting.

Optional inputs are independent:

| Input | Main option | Stages enabled or informed | Absent-input behaviour |
|---|---|---|---|
| Genome sequence | `--genome-fasta` | Sequence integrity, hard definite-stop check, CDS terminus refinement and sequence-related pruning | Coordinate-only candidate processing continues; sequence fields are inactive/neutral |
| Protein alignments | `--evidence-proteins` | Protein scoring, evidence rescue/protection and evidence-aware pruning | Effective protein selection contribution is zero |
| Transcript alignments | `--evidence-transcripts` | Transcript scoring/context, rescue/protection, pruning and UTR evidence catalogue | Effective transcript contribution is zero; UTR workflow has no transcript catalogue |

Corresponding `--evidence-*-source-name` options label evidence files. Manifest labels are resolved after expansion; explicit label lists must cover every expanded path and override manifest labels. Duplicate resolved paths and nested manifests are rejected.

`--strict` turns applicable malformed-input warnings/skips into errors. It changes parsing policy, not score weights.

The required `--out-prefix` is a filename prefix. Its parent must already exist and be writable. A genome directory may also need to be writable so a missing/stale `.fai` can be created or rebuilt.

## Normal workflow defaults

| Function | Option | Default | Scientific/operational effect |
|---|---|---:|---|
| Pipeline extent | `--stop-after` | `final` | Run selection, backfill, pruning, ID freeze, alternatives/refinement and final outputs where applicable |
| Reports | `--output-profile` | `standard` | Normal products and audit ledgers |
| Final pruning | `--pruning-level` | `auto` | Simulate canonical fixed scenarios and select an effective level |
| Alternative preset | `--alternative-splice-preset` | `intermediate` | Numerical eligibility/redundancy thresholds after primaries freeze |
| Alternative support | `--alternative-splice-support-policy` | `exact` | Require an exact qualifying support route |
| Alternatives emitted | `--no-alternative-splice` | false | Alternatives are evaluated by default at the final stage |
| UTR augmentation | `--add-utrs` | false | Write a separate post-selection UTR-augmented annotation |
| CDS refinement | `--refine-cds-termini` | true | Attempt bounded repair only when genome context makes it actionable |
| Strict backfill replacement | `--strict-r1-backfill-replacement` | true | Permit safeguarded atomic one-/two-blocker replacement |
| Multi-exon blocker extension | `--strict-r1-allow-multiexon-blockers` | true | Permit eligible multi-exon holder blockers as well as single-exon blockers |
| Progress | `--progress-log-level` | `basic` | Stage and periodic stderr progress |

`--stop-after scored-candidates`, `monosublocus`, or `backfill` is intended for investigation. Feature-dependent outputs are enabled when needed unless explicitly skipped.

## Structural comparison and grouping controls

These controls change relationships or feature matching, not aggregate weights:

| Option | Default | Direct mechanism | Downstream consequences |
|---|---:|---|---|
| `--single-exon-min-reciprocal-overlap` | 0.50 | Minimum overlap fraction on both single-exon CDS spans for a first-pass sublocus edge | Changes sublocus membership, score context and provisional winners |
| `--single-exon-support-min-reciprocal-overlap` | inherited | Single-exon structural-support relationship | Changes support components; does not connect single and multi-exon subloci |
| `--terminal-cds-3prime-wobble-bp` | 3 | Permitted strand-aware outer-boundary difference in explicitly wobble-aware terminal/chain matching | Can change terminal/exact support and relevant evidence matching |
| `--monosublocus-min-cds-exon-overlap-fraction` | 0.10 | Non-reciprocal maximum directional CDS-exon overlap for winner edges and backfill conflict | Changes suppression/backfill membership without changing scores |

Multi-exon first-pass sublocus edges always require at least one exact shared CDS splice-site coordinate. Single- and multi-exon models do not connect in that first pass. See [Selection and backfill](technical/selection-and-backfill.md).

## Model validation and hard eligibility

Model validation defaults to `exclude` above a 2,000,000-bp CDS span or a 1,000,000-bp CDS intron. Superlocus validation defaults to `stop` above a 5,000,000-bp span, 250 input genes or 100 genes from one source. These are operational safeguards, not biological plausibility claims.

Two early rules are distinct from scoring/pruning:

- definite internal-stop exclusion is enabled by default when sequence can be evaluated; `--no-definite-internal-stop-exclusion` is a global historical/diagnostic opt-out;
- `--minimum-coding-length-bp` defaults to 30 and is a hard pre-locus floor.

The separate `--min-cds-length-bp 150` only flags an intrinsic score concern for eligible candidates. Protein/transcript evidence cannot rescue a hard eligibility exclusion. See [Scoring](technical/scoring.md#genome-sequence-integrity-penalty).

## Scoring weights

Positive/negative outer aggregate defaults are:

| Component | Weight | Stage/evidence dependency |
|---|---:|---|
| Structural support | 1.00 | Model-only |
| Consensus similarity | 0.75 | Model-only |
| Pairwise similarity | 1.00 | Model-only |
| Completeness | 0.50 | Model-only |
| Source prior | 0.25 | Model-only; prior mode determines component value |
| Terminal extension | 0.10 | Model-only |
| Exact structure support | 0.75 | Model-only |
| Run-aware structural corroboration | 0.50 | Model-only; alpha 16 |
| Protein evidence | 3.00 | Effective only when protein evidence is supplied |
| Transcript evidence | 2.00 | Effective only when transcript evidence is supplied; uses context-adjusted score |
| Unsupported short CDS-exon penalty | 0.50 | Model-only |
| Structural-coherence penalty | 0.30 | Model-only |
| Sequence-integrity penalty | 1.00 | Effective only with genome sequence unless explicitly overridden in another context |
| Intrinsic penalty | 1.00 | Model-only |
| Partial-structure penalty | 0.50 | Model-only |
| Internal-extra penalty | 0.50 | Model-only |
| Terminal-extra penalty | 0.00 | Model-only; terminal extension is not penalized by default |

Positive options begin `--selection-weight-`; subtractive outer terms begin `--selection-penalty-`. Increasing a positive weight directly increases that component's contribution and the theoretical positive scaling denominator. Increasing a penalty weight makes the corresponding active concern more negative but does not itself create a hard exclusion.

The default model-only positive denominator is 4.85; active protein and transcript weights raise it to 7.85, 6.85 or 9.85. Changing positive weights therefore affects both signed scores and scaled threshold values. Read [Scoring](technical/scoring.md#model-only-aggregate) before tuning.

`--model-selection-tie-epsilon` defaults to `1e-6` and controls reported score-rank groups, which use competition ranking. It does not make the primary selector epsilon-based: exact raw scores define the total order, and `model_selection_order_rank` reports that same deterministic selector order.

## Intrinsic and sequence controls

Intrinsic metric thresholds default to 150-bp short CDS, 3-bp terminal CDS exon, 6-bp internal CDS exon, 20-bp short CDS intron, and a 99.999th-percentile run-derived long-intron cutoff unless an explicit long-intron value is supplied. The `--intrinsic-penalty-*` options construct `intrinsic_penalty_score`; `--selection-penalty-intrinsic` multiplies that completed score in the model aggregate.

Unsupported-moderately-short soft thresholds default to 15 bp for biological 5-prime terminal, 3-prime terminal and internal CDS exons. Their component weights are 0.25, 0.25 and 0.50. They require lack of other-source support and do not duplicate exons already below the intrinsic thresholds.

Sequence defaults are internal-stop 1.00 per stop capped at 2, phase inconsistency 0.75, missing stop 0.75, missing start 0.50 and noncanonical-intron fraction 1.00. `--no-sequence-integrity-scoring` preserves sequence metrics but makes their selection contribution zero and disables preset-driven sequence pruning. It does not disable the separate hard definite-internal-stop exclusion.

## Source prior

`--source-prior-mode auto` is the default. With fewer than three prediction sources, automatic weights are neutral at 1.0. With at least three, unique exact-coordinate agreement is combined from splice sites (0.50), introns (0.30) and internal CDS exons (0.20), renormalized over available classes, then mapped to the default 0.5–1.0 weight range.

`manual` uses repeated `--source-weight SOURCE=VALUE`, with unspecified sources at 1.0. `none` sets the effective prior component to zero. `--source-prior-scale` (default 1.0) multiplies the effective source weight before the outer source-prior weight is applied.

This changes candidate scores; it does not establish or report source accuracy. Inspect `source_auto_weights.tsv` and [Scoring](technical/scoring.md#source-prior).

## External-evidence scoring

| Family | Important controls | Direct effect | Inactive/missing behaviour |
|---|---|---|---|
| Protein | `--selection-weight-protein-evidence`; support-score threshold | Scales bounded source-collapsed protein structural support; threshold classifies supporting sources and affects evidence policies | No protein files gives effective selection weight zero |
| Transcript | `--selection-weight-transcript-evidence`; association/scoring controls | Scales context-adjusted transcript score; only scoring-eligible associations contribute | No transcript files gives effective selection weight zero |
| cDNA context | neutral fraction 0.30, severe fraction 0.15, maximum scoring adjustment 0.15 | Reduces the transcript score used in ranking for calculated low CDS/min-cDNA fractions; informs conservative pruning concern | Not calculable is neutral, not zero-quality evidence |

Setting an external-evidence selection weight to zero keeps parsing and diagnostics but makes that class ranking-neutral. Evidence can also participate in pruning rescue/protection under separate thresholds; changing a ranking weight does not automatically rewrite those predicates. Exact formulas are in [Scoring](technical/scoring.md#protein-evidence) and pruning use is in [Automatic pruning](technical/automatic-pruning.md).

## Selection and backfill controls

Strict replacement is enabled by default. It can replace exactly one or two uncorroborated monosublocus-selected holder models with a corroborated multi-exon candidate only when the candidate's **unscaled** score is strictly greater than every blocker. Scores are not summed. The multi-exon-blocker option extends eligible blockers; it does not loosen corroboration, provenance, blocker-count, strict-score or atomicity requirements.

Disabling strict replacement makes the multi-exon extension ineffective; an explicitly contradictory request is rejected. Exact conditions and rejection codes are in [Selection and backfill](technical/selection-and-backfill.md#strict-score-improving-replacement).

## Fixed and automatic pruning

Levels are `none`, `low`, `medium`, `high`, `very-high`, `auto` and `auto-conservative`; default is `auto`.

- Fixed levels apply named rules/thresholds directly. Use a fixed level when overriding per-threshold values.
- Automatic modes simulate canonical `medium`, `high` and `very-high` against the same pre-pruning population, assess evidence adequacy and the model-source consensus backbone, apply flood-aware safeguards where eligible, then choose an effective fixed result.
- `auto-conservative` uses the same simulations but a more conservative decision policy and default maximum.
- Recommendation and effective level can differ after user min/max constraints.
- Advisory source-count warnings never alter scoring, source inclusion or scenario contents.

Per-threshold singleton/nonsingleton overrides cannot be combined with automatic levels because automatic assessment is defined over canonical scenarios. Choose a fixed level for customized thresholds. See [Automatic pruning](technical/automatic-pruning.md) for the exact matrix, evidence adequacy, source-count-dependent consensus rules, rescue/protection ordering, monotonic closure and audit fields.

## Alternative-transcript controls

Alternative selection runs only at `--stop-after final`, after retained primaries and public `.1` identities are frozen. Main categories are:

- preset/threshold controls for numerical quality and redundancy;
- `--alternative-splice-support-policy` for qualifying prediction-model support;
- external-evidence requirements and thresholds;
- maximum alternatives/cap and deterministic greedy redundancy settings;
- `--alternative-splice-report` for detailed diagnostics;
- `--no-alternative-splice` for primary-only final emission.

A non-`none` external-evidence requirement needs compatible supplied evidence. Alternative controls cannot replace a primary or alter public gene/`.1` identity. See [Alternative transcripts](technical/alternative-transcripts.md).

## CDS refinement and UTR augmentation

`--refine-cds-termini` is enabled by default but is actionable only with suitable genome sequence and explicit incomplete-boundary status. Refinement is post-selection: it does not rescore candidates or change membership. Expert search-distance and safety controls affect bounded proposals/rollback. See [CDS terminus refinement](technical/cds-terminus-refinement.md).

`--add-utrs` requests a separate `<prefix>.final_annotation.with_utrs.gff3` after public transcript population/IDs are fixed. Transcript evidence association, source pooling, structure/endpoints, extension caps and sequence bounds are controlled by the UTR option family. UTR choices do not feed CDS scoring/selection and do not replace the ordinary final annotation. See [UTR augmentation](technical/utr-augmentation.md).

UTR splice-junction validation is active for the ordinary final annotation as
well as `--add-utrs`. It uses the established splice-motif definitions and the
existing effective UTR minimum-intron setting (20 bp by default). Inherited
source UTR sides are filtered only after selection; prediction- and
transcript-derived candidate sides are filtered before pooling and automatic
limit learning. Missing genome sequence leaves motif checks explicitly
unevaluated but does not by itself reject a side. The CDS-specific
`--sequence-penalty-noncanonical-splice` weight retains its existing meaning
and is not a UTR-policy control.

All profiles include `utr_junction_validation_summary_tsv`, which contains
bounded aggregate counts by evidence origin. Full/debug profiles include the
exception-only `utr_junction_validation_audit_tsv`; it can also be requested
explicitly in another profile. The detailed audit omits accepted sides,
passing junctions and routine `no_genome_fasta` rows. Reporting choices do not
alter validation decisions, support pooling, limit learning or UTR selection.

## Output and reporting controls

Profiles are `minimal`, `standard`, `full` and `debug`. Full/debug enable all normal output keys; strict-replacement reports remain feature-gated. `--write-output KEY` adds products and `--skip-output KEY` removes them; both accept repeated/comma-separated keys and `all`.

`--gff3-include-model-selection-score-components` can make GFF3 much larger and is intended for explicit audit needs. `--show-zero-weight-score-components` exposes disabled terms in component strings; by default they are omitted to avoid suggesting participation.

Reporting controls do not change scientific behavior unless an option explicitly requests/omits a post-selection product. File membership, fields and state meanings are in [Outputs](outputs.md).

## Important interactions

- Protein/transcript selection weights act only when the class is supplied; absence is neutral.
- `--add-utrs`, alternative selection and their detailed final reports require the complete final stage.
- Alternative selection, CDS refinement and UTR augmentation are post-selection and cannot revise selection membership.
- Pruning occurs after backfill and does not trigger a new backfill pass.
- Fixed pruning threshold overrides require a fixed level; automatic modes keep canonical scenario definitions.
- A stale/missing FASTA index can require write permission beside the genome.
- Source-level mean scaled scores are descriptive within-run summaries, not source accuracy or source-subset recommendations.

For every exact option spelling, accepted range and parser validation, use `moselecta --expert-help`; use `run_configuration.json` to verify what the run actually applied.
