# Scoring and model ranking

> Technical reference for MoSelecta `0.2.16`. The implementation in `structural.py`, `scoring.py`, `protein_evidence.py`, `transcript_evidence.py`, `loci.py` and `config.py` is authoritative.

## Scope and terminology

MoSelecta scores CDS-containing candidate models within a run. Scores describe implemented structural agreement, candidate properties and optional evidence. They are not probabilities that a gene model is biologically correct.

Four score levels must remain distinct:

1. **Component metrics** describe one aspect of one candidate, such as `structural_support_score` or `protein_evidence_score`.
2. **`model_structural_support_score`** is the signed, weighted model-only aggregate. Despite its name, it includes every active model-only positive and negative component, not only `structural_support_score`.
3. **`model_selection_score`** is the signed integrated aggregate. It equals the model-only aggregate plus active protein and transcript contributions and is the principal score used by candidate selection and strict replacement.
4. **Scaled and source-summary fields** transform or summarize scores; they do not replace the model-level aggregate.

`model_selection_score` is not a raw input metric. It already combines weighted positive and negative components, so calling it a “raw score” would be misleading. “Provisional” describes a first-pass winner that can later be suppressed or replaced; it does not describe the score.

See [Selection, conflict resolution and backfill](selection-and-backfill.md) for decision rules, [Configuration](../configuration.md) for controls, and [Outputs](../outputs.md) for schemas.

## Coordinate and feature conventions

Candidate comparison is CDS-centred and strand-aware. The main feature families are CDS introns, CDS splice-site coordinates, internal CDS exon intervals, biological 5-prime and 3-prime terminal CDS exon intervals, CDS-base intervals, and complete intron or CDS-exon chains.

Exact features are source-collapsed where a score asks for independent prediction-source support: repeated models from one source do not create multiple source votes. A biological 3-prime terminal CDS outer boundary may differ by `--terminal-cds-3prime-wobble-bp` (default 3 bp) where the matching routine explicitly permits it. Internal boundaries and all other exact coordinates remain exact.

For precision/sensitivity/F1 metrics:

```text
precision   = matched query features / query features
sensitivity = matched target features / target features
F1          = harmonic mean of precision and sensitivity
```

An inapplicable component is omitted and remaining weights are renormalized. This matters particularly for internal CDS exons: a model with fewer than three CDS exons has no internal CDS exon, so the internal component is not applicable rather than a failed score of zero. A genuinely applicable but unsupported feature receives zero.

## Prediction-model components

### Structural support

`structural_support_score` measures direct support for a candidate's features from **other prediction sources in the same first-pass sublocus**. For each relevant feature, support means at least one model from a different source supplies the matching feature. Complete-chain fractions are scaled by the possible other-source count.

For a multi-exon model the weighted mean is:

| Component | Weight | Meaning |
|---|---:|---|
| CDS intron support fraction | 0.25 | Fraction of query CDS introns supported by another source |
| CDS splice-site support fraction | 0.20 | Fraction of query CDS splice sites supported by another source |
| Internal CDS exon support fraction | 0.20 | Fraction of applicable internal CDS exons supported by another source |
| Terminal CDS exon support fraction | 0.15 | Fraction of terminal CDS exons supported, with configured 3-prime wobble |
| Exact intron-chain support fraction | 0.20 | Other-source support for the complete intron chain divided by possible other sources |

For a single-exon model:

| Component | Weight | Meaning |
|---|---:|---|
| Single-exon reciprocal-overlap support fraction | 0.60 | Other-source support at the configured overlap relationship |
| Exact single-exon interval support fraction | 0.40 | Exact other-source interval support divided by possible other sources |

When there is no possible other source, source-scaled support is zero rather than not applicable: the run provides no independent prediction-source support. Diagnostic fields include `structural_support_score_type`, `structural_support_score_components` and feature-level counts/fractions.

Increasing `--selection-weight-structural-support` increases this aggregate contribution. It does not change how features match, and it does not guarantee retention because other scores and conflict rules still apply.

### Consensus similarity

`consensus_similarity_score` compares the candidate with features supported by at least two distinct prediction sources in its sublocus. It asks how closely the candidate resembles the source-collapsed sublocus consensus; it is not the same as counting support for the candidate's own features.

Multi-exon weights are:

| Consensus F1 component | Weight |
|---|---:|
| CDS introns | 0.30 |
| CDS splice sites | 0.25 |
| Internal CDS exons | 0.20 |
| CDS bases | 0.15 |
| Terminal CDS exons | 0.10 |

Single-exon weights are 0.75 for CDS-base consensus F1 and 0.25 for the single-exon overlap-support fraction. The implementation also records weighted consensus precision and sensitivity. Those fields drive:

```text
partial_structure_penalty = max(0, consensus_precision - consensus_sensitivity)
extra_structure_penalty   = max(0, consensus_sensitivity - consensus_precision)
```

Only `partial_structure_penalty` is a separately weighted aggregate penalty. Extra CDS is represented more directly by `internal_extra_fraction` and `terminal_extra_fraction`.

### Pairwise similarity

`pairwise_similarity_score` compares a candidate with individual models from other prediction sources. MoSelecta calculates every pairwise score, keeps the best-scoring model independently for each other source, then takes the arithmetic mean of those source-best scores. Each other source contributes at most once.

Multi-exon pairwise weights are:

| Pairwise F1 component | Weight |
|---|---:|
| CDS introns | 0.30 |
| CDS splice sites | 0.25 |
| Internal CDS exons | 0.20 |
| CDS bases | 0.15 |
| Terminal CDS exons | 0.10 |

A single-exon pairwise score is CDS-base F1 with weight 1.0. The score can recognize a strong near match even when no whole-chain exact match exists. `pairwise_similarity_score` is distinct from exact support and corroboration.

`pairwise_other_source_count`, the best model/source, source-best feature means and `pairwise_similarity_score_components` support audit. A candidate with no other-source comparison has no available pairwise mean; aggregate conversion treats the unavailable numeric field as zero.

### Completeness

`model_completeness_score` asks whether the candidate recovers source-supported CDS structure without looking like a terminally contained fragment. Consensus-supported bases are split relative to the candidate span into missing terminal and missing internal structure. A directional containment term is used only when another-source intron chain contains the query, the reverse containment is absent, and supported terminal structure is missing.

```text
model_completeness_penalty = min(
    1,
    0.50 × missing_supported_terminal_structure_fraction
  + 0.35 × missing_supported_internal_structure_fraction
  + 0.15 × contained_model_source_fraction
)

model_completeness_score = max(0, 1 - model_completeness_penalty)
```

If no consensus CDS bases exist, the completeness score is not applicable and its aggregate numeric fallback is zero. This score rewards recovery of supported structure; it does not assert that longer is always better.

### Exact structure support

`exact_structure_support_score` counts other sources in the current sublocus whose complete ordered CDS-exon chain matches the candidate, allowing only the configured biological 3-prime terminal wobble. It divides that count by the **run-wide possible other-source count**.

This differs from the local exact-chain fraction, which uses the sublocus's possible other sources, and from exact intron-chain support, which ignores terminal CDS boundaries. It also differs from run-aware corroboration: exact support is derived from local sublocus comparisons, whereas corroboration is a separately constructed whole-run fact with saturating normalization.

### Run-aware structural corroboration

`structural_corroboration_score` asks whether distinct other prediction sources supply the same complete CDS-exon chain anywhere in the run. Duplicate models/isoforms from one source count once. Let `k` be the corroborating-other-source count, `N` the run-wide possible other-source count, and `a` `--structural-corroboration-alpha` (default 16):

```text
structural_corroboration_score = 0                         if k = 0 or N = 0
structural_corroboration_score = (k / (k + a)) / (N / (N + a)) otherwise
```

The score is bounded and gives full score when all possible other sources corroborate. A smaller alpha changes curvature; increasing the aggregate weight changes influence but not source facts. Corroboration is within one execution and must not be interpreted as independent evidence when inputs are correlated. `structural_corroboration.tsv` reports the support facts independently of whether its aggregate weight is zero.

### Source prior

`source_prior_score` is a source-level component copied to every candidate from that prediction source:

```text
source_prior_score = effective_source_weight × source_prior_scale
```

`--source-prior-mode` selects:

- `auto` (default): derive agreement from unique exact-coordinate features. Components are splice-site agreement 0.50, intron agreement 0.30 and internal-CDS-exon agreement 0.20, renormalized over classes that exist. If none applies, exact CDS-exon agreement is the fallback. With at least three prediction sources, agreement maps linearly to `--auto-source-weight-min`..`--auto-source-weight-max` (defaults 0.5..1.0). With fewer than three sources, automatic weight is neutral at 1.0.
- `manual`: use repeated `--source-weight SOURCE=VALUE`; unspecified sources use 1.0.
- `none`: effective weight is zero.

The default `--source-prior-scale` is 1.0 and the default aggregate weight is 0.25. The automatic prior describes agreement within this run; it is not a validated predictor-accuracy estimate. `source_auto_weights.tsv` records its inputs and effective values.

### Terminal extension and extra structure

Consensus-relative extra candidate CDS bases are classified as terminal or internal. `terminal_extension_score` is `terminal_extra_bp / total_extra_cds_bp` when extra CDS exists; otherwise it is not applicable. It provides a small positive term (default weight 0.10) when extra sequence is terminal rather than internal. `internal_extra_fraction` is penalized at 0.50 by default. `terminal_extra_fraction` has default penalty weight 0, so useful terminal extension is not simultaneously punished.

### Intrinsic penalty

`intrinsic_penalty_score` represents candidate properties that do not require genome sequence. It is an unbounded weighted sum of counts or fractions:

| Flag or fraction | Default intrinsic weight | Trigger/default threshold |
|---|---:|---|
| Short CDS | 1.00 | CDS length below `--min-cds-length-bp 150` |
| CDS length not divisible by 3 | 1.00 | Modulo-three failure |
| Short terminal CDS exon | 0.25 per exon | Below 3 bp by default |
| Short internal CDS exon | 0.75 per exon | Below 6 bp by default |
| Short CDS intron | 1.00 per intron | Below 20 bp |
| Long CDS intron | 0.75 per intron | Above explicit cutoff or default run percentile 99.999 |
| Likely contained model | 0.50 × source fraction | Directional containment with missing supported terminal structure |
| Unsupported intron inside another CDS exon | 0.25 × affected-intron fraction | No other-source exact intron support and full containment in another sublocus model's CDS exon |

The aggregate subtracts `--selection-penalty-intrinsic` (default 1.0) times this score. The 150-bp short-CDS threshold is not the hard minimum coding length; eligible candidates between the hard floor and 149 bp remain scoreable.

### Unsupported moderately short CDS-exon penalty

This source-support-aware penalty is separate from `intrinsic_penalty_score`, so it cannot by itself trigger intrinsic-singleton pruning. It applies only to moderately short exons at or above the active intrinsic minimum but below the soft threshold (15 bp by default) and lacking other-source support.

For exon length `L`, intrinsic threshold `I` and soft upper threshold `S`:

```text
severity = (S - L) / (S - I)   for I <= L < S
severity = 0                   otherwise
```

Exact other-source internal-exon support neutralizes internal severity. Terminal support is biological-side aware and applies configured 3-prime wobble. The normalized penalty is:

```text
0.25 × unsupported 5-prime terminal severity
+ 0.25 × unsupported 3-prime terminal severity
+ 0.50 × min(1, sum of unsupported internal severities)
```

It is bounded at 1 and multiplied in the aggregate by `--selection-penalty-unsupported-short-cds-exon` (default 0.50). Exons below the intrinsic threshold are handled only by the intrinsic penalty. If an intrinsic threshold equals or exceeds its soft threshold, that interval is disabled.

### Structural-coherence penalty

The structural-coherence subsystem asks whether a candidate connects two supported coding anchors while source-collapsed evidence favours separate gene hypotheses. It produces `structural_coherence_raw_score`; the model aggregate subtracts that score times `--selection-penalty-structural-coherence` (default 0.30).

A region is eligible only when separation votes `S >= 2`, `S > C` continuity votes and each anchor has at least two source supports. Directional evidence is `(S-C)/(S+C)` and broader separation support can increase the raw value, bounded at 1.5. A candidate uses its maximum eligible region value, not a sum. See [Structural coherence](structural-coherence.md) for source votes, bridge groups and exact gates.

### Genome sequence-integrity penalty

When `--genome-fasta` is supplied, MoSelecta assesses translated CDS and splice properties under genetic code 1. The default score is:

```text
sequence_integrity_penalty_score =
    1.00 × min(internal_stop_count, 2)
  + 0.75 × phase_inconsistent
  + 0.75 × missing_terminal_stop
  + 0.50 × missing_start
  + 1.00 × noncanonical_splice_intron_fraction
```

GT-AG, GC-AG and AT-AC are canonical transcript-oriented introns. A terminal stop may be included in the CDS or be the immediately adjacent in-frame codon. Ambiguous start/stop contexts and an absent start/stop at the exact corresponding contig boundary are neutral. Motifs containing ambiguous bases are noncanonical for the fraction.

The aggregate subtracts the penalty times `--selection-penalty-sequence-integrity` (default 1.0). `--no-sequence-integrity-scoring` keeps metrics but sets the selection contribution to zero and disables preset-driven sequence pruning; explicit pruning overrides remain possible. No FASTA means the component is inactive/neutral.

Sequence-integrity scoring is distinct from two hard eligibility rules:

- the default definite-internal-stop stage excludes a resolved frame with a literal internal TAA/TAG/TGA, or an unresolved model only if all three frames contain a literal internal stop; invalid overlapping CDS geometry is excluded at the same stage;
- the default minimum coding length is 30 coding nucleotides after accepted early corrections and included-stop accounting.

These hard pre-locus exclusions cannot be rescued by protein/transcript evidence and are not score penalties. `--min-cds-length-bp 150` remains only an intrinsic scoring threshold.

## Model-only aggregate

With `w(x)` denoting configured weights and `f` the low-source factor below:

```text
model_structural_support_score =
  + w(structural_support)        × structural_support_score × f
  + w(consensus_similarity)     × consensus_similarity_score × f
  + w(pairwise_similarity)      × pairwise_similarity_score
  + w(completeness)             × model_completeness_score
  + w(source_prior)             × source_prior_score
  + w(terminal_extension)       × terminal_extension_score
  + w(exact_structure_support)  × exact_structure_support_score
  + w(structural_corroboration) × structural_corroboration_score
  - w(unsupported_short)        × unsupported_short_cds_exon_penalty_score
  - w(structural_coherence)     × structural_coherence_raw_score
  - w(sequence_integrity)       × sequence_integrity_penalty_score
  - w(intrinsic)                × intrinsic_penalty_score
  - w(partial_structure)        × partial_structure_penalty
  - w(internal_extra)           × internal_extra_fraction
  - w(terminal_extra)           × terminal_extra_fraction
```

Default weights are:

| Positive term | Weight | Negative term | Weight |
|---|---:|---|---:|
| Structural support | 1.00 | Unsupported short CDS exon | 0.50 |
| Consensus similarity | 0.75 | Structural coherence | 0.30 |
| Pairwise similarity | 1.00 | Sequence integrity | 1.00 |
| Completeness | 0.50 | Intrinsic | 1.00 |
| Source prior | 0.25 | Partial structure | 0.50 |
| Terminal extension | 0.10 | Internal extra | 0.50 |
| Exact structure support | 0.75 | Terminal extra | 0.00 |
| Structural corroboration | 0.50 |  |  |

The default model-only positive theoretical maximum is 4.85. Penalty weights do not enter that positive denominator. Missing numerical components fall back to zero in the aggregate. `model_structural_support_score_components` records active signed contributions; zero-weight terms are omitted unless requested.

### Low-source conflict dampener

The agreement factor `f` is normally 1.0. By default it becomes 0.25 only when all are true:

```text
sublocus_source_count == 2
sublocus_model_count == 2
pairwise_other_source_count > 0
exact_structure_support_score == 0
pairwise_similarity_score < 0.80
```

Only structural-support and consensus-similarity contributions are damped. Component values and other terms are unchanged. `--no-low-source-conflict-dampen-agreement-scores` disables multiplication, while diagnostic flag/factor/reason fields preserve the detected case.

## Protein evidence

Each protein-evidence file is a separate source. CDS alignments are grouped in same-strand evidence superloci, and exact CDS-exon chains are collapsed before comparison. For each evidence source, the best collapsed alignment supplies pairwise metrics.

For single-exon models:

```text
protein_evidence_score =
    0.70 × mean source-best CDS-base F1
  + 0.15 × exact-CDS-chain source fraction
  + 0.15 × supporting-source fraction
```

For multi-exon models the renormalized weighted mean is:

```text
0.25 × CDS-intron F1
+ 0.20 × CDS-splice-site F1
+ 0.20 × internal-CDS-exon F1 (when applicable)
+ 0.15 × CDS-base F1
+ 0.10 × exact-CDS-chain source fraction
+ 0.10 × supporting-source fraction
```

A source supports a model when its best score reaches `--protein-evidence-support-score-threshold` (default 0.5). With no configured sources the status is `protein_evidence_not_configured` and effective aggregate weight is zero. Overlap without sufficient similarity receives an explicit low-similarity status. The score is already bounded 0–1; its scaled field is the clamped form of the same value.

The default integrated weight is 3.0. Setting `--selection-weight-protein-evidence 0` retains reports but removes ranking contribution.

## Transcript evidence

Transcript evidence is built from exons. Each file is a source; multi-exon transcripts are collapsed by exact intron chain and single-exon transcripts by exact exon interval. Broad spatial associations remain diagnostic, but only scoring-eligible associations contribute to `transcript_evidence_score` or evidence rescue.

The pairwise score for a scoring-eligible multi-exon association is a renormalized weighted mean:

| Pairwise transcript component | Weight |
|---|---:|
| Model CDS-intron support | 0.30 |
| Splice-site F1 | 0.20 |
| Full intron-chain F1 | 0.20 |
| Internal CDS-exon compatibility | 0.15 |
| Terminal CDS-exon compatibility | 0.10 |
| Model CDS exonic coverage | 0.05 |

For a single-exon model it is `0.85 × model CDS exonic coverage + 0.15 × terminal compatibility`.

The model-level multi-exon score is:

```text
0.40 × mean best-source pairwise score
+ 0.25 × source-weighted mean scoring-eligible pairwise score
+ 0.25 × chain-source score
+ 0.10 × support confidence
```

The single-exon score is:

```text
0.45 × mean best-source pairwise score
+ 0.35 × source-weighted mean scoring-eligible pairwise score
+ 0.15 × supporting-source fraction
+ 0.05 × log-scaled member support
```

`support confidence = 0.80 × supporting-source fraction + 0.20 × log-scaled member support`. The source support threshold is 0.50. Detailed association, chain-relation and source-best fields provide audit; source-best internal-exon values preserve not-applicable state.

### cDNA-context adjustment

The cDNA-context diagnostic selects a compatible transcript independently of the general best-pairwise transcript. A single-exon CDS must be contained by an outermost-envelope exon. A multi-exon CDS chain must be a contiguous subchain of the transcript chain with every CDS exon contained by its corresponding transcript exon.

For calculated context, let `q = CDS length / minimum compatible cDNA length`. With default neutral fraction 0.30:

```text
concern = max(0, min(1, (0.30 - q) / 0.30))
requested penalty = 0.15 × concern
transcript_evidence_score_context_adjusted =
    max(0, transcript_evidence_score - requested penalty)
```

Context that cannot be calculated is neutral: concern and penalty are zero and undefined fields remain `.`. A separate severe fraction 0.15 contributes to conservative pruning classification; a fraction from 0.15 to below 0.30 requires one of the documented structural-context flags to become a pruning concern.

The default integrated transcript weight is 2.0. Setting it to zero keeps diagnostics but removes the score and context adjustment from ranking.

## Integrated model-selection score

When evidence classes are active:

```text
model_selection_score =
    model_structural_support_score
  + protein_weight × protein_evidence_score
  + transcript_weight × transcript_evidence_score_context_adjusted
```

An evidence class is active only when supplied and its effective weight is positive. Supplying no files gives effective weight zero; it does not subtract a missing-evidence penalty. The unadjusted transcript score and applied context penalty remain auditable. Component, score-delta and rank-change fields show how integration changed ordering.

Increasing an evidence weight directly increases a positive evidence contribution and the positive scaling denominator. It may change ordering and downstream decisions, but does not guarantee retention because conflicts, tie fields, backfill eligibility and pruning operate separately.

## Scaling

Scaled model aggregates are review-friendly positive transforms, not probabilities or alternative ranking scores:

```text
positive_max = sum of configured positive model weights
             + active protein weight
             + active transcript weight

scaled = 0                            if score <= 0 or positive_max <= 0
scaled = min(score / positive_max, 1) otherwise
```

Negative penalty bounds are not in the denominator. With default weights:

| Active score family | Positive denominator |
|---|---:|
| Model-only | 4.85 |
| Model + protein | 7.85 |
| Model + transcript | 6.85 |
| Model + protein + transcript | 9.85 |

`model_structural_support_score_scaled` always uses the model-only denominator. `model_selection_score_scaled` includes only active positive external-evidence weights. Zero and negative signed scores collapse to 0, so scaled values cannot reconstruct ordering among non-positive candidates. Unscaled scores drive ordering and strict replacement. Scaled values are used by configured pruning and alternative-transcript thresholds and source summaries.

Do not compare scaled scores across runs with materially different positive weights/evidence participation without accounting for the denominator.

## Ranks and deterministic ordering

MoSelecta reports competition score-rank and unique order-rank fields within each scored canonical sublocus. Scores whose absolute difference is at most `--model-selection-tie-epsilon` (default `1e-6`) share a score rank. Competition numbering advances by the tied group size, for example `1, 1, 3`; it is not dense ranking. Epsilon grouping does not replace exact-value total ordering.

`model_selection_order_rank` applies the shared primary-selector order:

1. higher aggregate score;
2. higher completeness;
3. lower intrinsic penalty;
4. higher consensus similarity;
5. higher structural support;
6. higher source prior;
7. longer CDS;
8. lexical output transcript ID.

The aggregate and all secondary numeric values are compared exactly in this total order. Secondary fields are reached only for exact aggregate-score equality. Tie-group fields still describe the epsilon-based score-rank group. Within such a group, the tie-break-reason field records `model_selection_score` when near-but-unequal raw scores determine the order; for exactly equal scores it records the first differentiating shared secondary field.

The same authoritative helper supplies first-pass primary selection, monosublocus ordering, ordinary backfill ordering and this reported order rank. The rank is assigned before those later stateful stages and reports their shared candidate order; the rank field itself is not read back as a selection input and does not guarantee final retention. `model_structural_support_order_rank` retains its separate model-only diagnostic ordering. See [deterministic comparators](selection-and-backfill.md#deterministic-comparators).

## Source-level summaries

`mean_scaled_score_all` is the arithmetic mean of `model_selection_score_scaled` for every scored candidate from a source. `mean_scaled_score_non_singletons` excludes candidates whose canonical `sublocus_model_count == 1`; when none qualify, the mean is `.` and `n_scored_non_singletons` is 0.

A single-candidate sublocus is not a single-exon model and is not simply “a model from one source.” These means are within-run descriptive diagnostics, not model-level selection scores, calibrated source accuracy, or source-subset recommendations.

## Main controls and direct effects

| Control family | Direct effect | Does not directly do |
|---|---|---|
| Positive `--selection-weight-*` | Multiplies named positive contribution and enters positive scaling maximum | Change feature matching |
| Aggregate `--selection-penalty-*` | Multiplies named subtractive contribution | Create a hard exclusion unless a separate rule says so |
| `--intrinsic-penalty-*` | Changes how intrinsic facts form the intrinsic score | Change the outer aggregate intrinsic weight |
| `--terminal-cds-3prime-wobble-bp` | Changes permitted outer-boundary matching in wobble-aware comparisons | Relax internal boundaries |
| `--single-exon-*-overlap` | Changes single-exon clustering or support relationships | Connect single- and multi-exon first-pass subloci |
| Protein/transcript weights | Activate/scale evidence contribution when supplied | Penalize an absent evidence class |
| cDNA-context controls | Change context concern, adjusted transcript contribution and pruning classification | Alter the unadjusted transcript score |
| Tie epsilon | Changes score-rank tie grouping | Replace deterministic tie fields |

Exact defaults and interactions are in [Configuration](../configuration.md). `--expert-help` is the complete machine-current CLI authority.

## Audit trace and interpretation limits

For one model, inspect component facts, optional evidence model-support rows, `model_selection_scores.tsv`, winner/selected/suppressed/backfill ledgers, pruning ledgers, then alternative-splice reports. Exact output membership and fields are in [Outputs](../outputs.md).

- Agreement among correlated prediction sources is not independent validation.
- A score is meaningful only with its configured weights, active evidence and comparison context.
- Scaling clips negative values and is not a calibrated percentage of quality.
- Hard exclusion, conflict suppression and pruning are different mechanisms from score subtraction.
- A higher score does not override structural conflict rules or every pruning rule.
