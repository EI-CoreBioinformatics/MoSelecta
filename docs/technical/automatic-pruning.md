# Automatic final pruning

> Technical reference for MoSelecta `0.2.16`.
> Automatic pruning is a post-selection assessment; it does not change scoring,
> conflict resolution, or backfill. See [Outputs](../outputs.md) for its reports.

For the workflow-level rationale, start with the user guide. User controls and
defaults are summarized in [Configuration](../configuration.md); the signed and
scaled score fields consumed by pruning are defined in [Scoring](scoring.md),
and the population entering this stage is defined in
[Selection and backfill](selection-and-backfill.md).

Automatic modes select and report an effective final-pruning level. The fixed
medium/high/very-high scenarios retain their own pruning predicates; automatic
assessment does not alter model selection, scoring, or backfill.

## Nonsingleton model-confidence fallback

The medium, high, and very-high simulations enable
low-score nonsingleton pruning even when external evidence is absent or no
supplied evidence class passes adequacy. Those cases use the canonical
model-confidence fallback defined below. Adequate
evidence continues to select the unchanged evidence-aware regime. The
evidence-driven automatic pathway and the no-adequate-evidence maximum of high
are unchanged.

## Per-threshold overrides

`auto` and `auto-conservative` always simulate the canonical `medium`, `high`,
and `very-high` configurations. The eight per-threshold singleton/nonsingleton
model, protein, transcript, and model-selection-score low-score options therefore cannot
be combined with either automatic pruning level; MoSelecta stops at startup
with an error that names the conflicting option(s).

Choose a fixed `none`, `low`, `medium`, `high`, or `very-high` level to use
custom threshold values. Explicit values override that fixed preset. With a
fixed level plus `--auto-pruning-report-only`, the automatic recommendations
remain based on canonical scenarios and the requested fixed output uses the
explicit threshold values.

## Monotonicity and rescue semantics

Each medium, high, and very-high scenario is still evaluated natively against
the same immutable pre-pruning selected-model population.  A monotonic closure
is then applied:

```text
effective removed(high) = native removed(high) union effective removed(medium)
effective removed(very-high) = native removed(very-high) union effective removed(high)
```

This preserves:

```text
retained(very-high) subset retained(high) subset retained(medium)
```

When a lower-level removal is inherited, the current level records
`inherited_lower_level_removal` as its effective primary reason while separately
retaining the original level, original primary reason, and original prune
classes.  Thus a model removed as an opposite-strand shadow at medium remains
removed at high even if its partner is removed natively at high, but the audit
does not falsely claim that high independently made a shadow decision.

Strong local protein and transcript singleton protection is now evaluated as a
route-independent eligibility state.  It protects an otherwise removable
singleton during low-score pruning even when the earlier intrinsic/sequence
singleton route did not need to invoke rescue.  The `*_rescue_applied` fields
mean that protection actually prevented removal in the current scenario; they
are not used as the underlying eligibility test.

Automatic pruning is a reference-free choice among the independently usable
fixed `medium`, `high`, and `very-high` presets. It is performed after final
model selection and before the pruned outputs are written. All scenarios use
the same selected-model population and cached structural, sequence, protein,
and transcript metrics.

`auto` and `auto-conservative` simulate identical scenarios. They differ only
in how a scenario is selected.

## Evidence adequacy

Protein adequacy requires:

```text
model overlap fraction >= 0.75
supported-model fraction >= 0.50
```

Transcript adequacy requires:

```text
scoring-eligible-model fraction >= 0.75
supported-model fraction >= 0.50
```

Each evidence class is reported as:

- `absent`: not supplied;
- `pass`: both thresholds pass;
- `marginal`: the limiting margin is below zero but at least `-0.05`;
- `fail`: the limiting margin is below `-0.05`.

Only `pass` counts as adequacy-passing. Marginal or failed global evidence
cannot select the evidence-aware nonsingleton regime; it selects the guarded
model-confidence fallback instead. Supplied high local protein or transcript
scores veto its conventional route.

## Model-source consensus metrics

For `L_k`, the number of model-input subloci represented by exactly `k`
prediction sources:

```text
L = sum(L_k)
I = sum(k * L_k)
```

MoSelecta reports locus fractions for one source, at least two sources, at
least three sources, and all available sources. It also reports:

```text
singleton-incidence fraction = L_1 / I
consensus-incidence fraction = sum(k*L_k, k>=2) / I
deep-consensus incidence fraction = sum(k*L_k, k>=3) / I
corroboration fraction = (I-L) / I
```

Incidence-weighted metrics distinguish a large consensus backbone plus a tail
of isolated predictions from widespread disagreement between predictors.

The current implementation treats every source label as independent:

```text
effective source count = raw source count
```

Source-family grouping is deferred. Where input labels are related annotation
tiers or share most of a prediction pipeline, confidence should be interpreted
accordingly.

## Backbone classes

Deep consensus means at least three sources.

Strong:

```text
source count >= 3
multisource subloci >= 100
singleton-incidence fraction <= 0.35
corroboration fraction >= 0.45
deep-consensus incidence fraction >= 0.55 with exactly 3 sources
deep-consensus incidence fraction >= 0.60 with 4 or more sources
deep-consensus sublocus fraction >= 0.25
```

Intermediate:

```text
source count >= 3
multisource subloci >= 100
singleton-incidence fraction <= 0.55
corroboration fraction >= 0.30
deep-consensus incidence fraction >= 0.35
deep-consensus sublocus fraction >= 0.15
```

Anything failing the intermediate definition is `weak_or_discordant`. With
three sources, deep consensus is unanimity. With four or more sources it still
means support from at least three sources. All other strong and intermediate
thresholds are the canonical current values shown above.

For each source, MoSelecta reports its share of all source-singleton subloci,
its source-specific singleton rate, and the difference and ratio relative to
the median rate of other sources. A strong source excess requires:

```text
singleton share >= 0.60
source-specific singleton rate >= 0.20
rate >= 2 * median rate of other sources
```

Exactly one source must satisfy all three criteria to qualify as strongly
excessive. Detection always uses the original, unadjusted distribution.

## Flood-aware backbone assessment

When source-excess handling is enabled and exactly one source qualifies,
MoSelecta discounts all subloci belonging exclusively to that source from a
second statistical view of the source-consensus distribution. If `d` such
subloci exist:

```text
adjusted L1 = L1 - d
adjusted L = L - d
adjusted I = I - d
```

Every multisource sublocus is preserved, including all corroborated
contributions from the excessive source. The original number of input sources
is also preserved. Metrics and denominators are then recalculated and the
adjusted distribution is classified using the ordinary thresholds above.

This is a backbone-assessment adjustment only. It does not delete a sublocus,
remove a prediction model, change the scenario populations, or exclude a
source. The detector is not rerun on adjusted metrics.

## Decision rules without adequacy-passing evidence

- Fewer than three model sources: both modes select `medium`.
- Weak/discordant class used by a mode: it selects `medium`.
- Strong class used by a mode: it selects `high`.
- Intermediate backbone: `auto-conservative` selects `medium`. Standard
  `auto` selects `high` only when exactly one source has a strong singleton
  excess, the high-versus-medium increment is at least 10 models and at least
  0.0005 of the original selected population, and at least 70% of those
  additional removals come from the same excess source.
- Neither mode selects `very-high` without an adequacy-passing evidence class.

Standard `auto` uses the flood-adjusted class when the adjustment is eligible
and enabled. Otherwise it uses the unadjusted class. `auto-conservative` always
uses the unadjusted class (Option B), although adjusted metrics may still be
reported. Both modes apply the revised source-count-dependent deep-consensus
threshold. This keeps balanced or distributed disagreement at `medium` while
allowing one flooding source to stop masking an otherwise strong backbone.

When an adjustment is applied, the normal run log displays raw and adjusted
singleton incidence, corroboration, deep-consensus incidence, and
deep-consensus sublocus fraction beside their strong thresholds. It also names
the excessive source, discount count, class transition, and the raw or adjusted
class that determined the recommendation. Automatic-decision explanations
separately state the unconstrained recommendation and any user minimum or
maximum that changes the effective level. These are reporting changes only.

## Decision rules with adequacy-passing evidence

Standard `auto` selects at least `high` when:

```text
selected singleton-model fraction >= 0.10
OR
high increment / original selected models >= 0.02 and increment >= 10
```

From `high`, it selects `very-high` only when the increment is meaningful and:

```text
very-high increment / high-retained models <= 0.05
```

`auto-conservative` selects `high` when:

```text
selected singleton-model fraction >= 0.15
OR
high increment / original selected models >= 0.03 and increment >= 10
```

Its default maximum is `high`. It considers `very-high` only when the user
explicitly supplies `--auto-pruning-max-level very-high`; the additional
removal limit is then 0.02 of high-retained models.

## Consistency safeguards

Every simulation must satisfy:

```text
starting = removed + retained
retained(very-high) subset retained(high) subset retained(medium)
```

If reconciliation or monotonicity fails, automatic selection falls back to
`medium`, confidence is `low`, and a warning is written to stdout and the
diagnostics table.

## Controls and reporting

```text
--auto-pruning-min-level {low,medium,high}
--auto-pruning-max-level {medium,high,very-high}
--auto-pruning-max-extra-removal-fraction 0.05
--auto-pruning-high-min-singleton-fraction 0.10
--auto-pruning-high-min-prospective-removal-fraction 0.02
--auto-pruning-use-source-excess
--no-auto-pruning-use-source-excess
--auto-pruning-report-only
```

The diagnostics report both the unconstrained recommendation and the effective
level after user minimum/maximum constraints. `--auto-pruning-report-only`
preserves a fixed requested output level while writing the same scenario table
and both automatic recommendations.

The stable output files are:

```text
<prefix>.pruning_scenario_summary.tsv
<prefix>.auto_pruning_diagnostics.tsv
```

Established `source_consensus` row names retain their raw meanings. New
append-only rows report the unadjusted class, any adjusted metrics and
class, the qualifying source and discount count, thresholds, criterion
pass/fail states and margins, and the class actually used by each mode.

Each input source also receives an advisory coding-model count comparison
against one common median calculated across all input sources in that run.
Ratios at or above 2.0 receive a high-count warning; ratios at or below 0.5
receive a low-count warning. The established leave-one-out median and ratio
remain available as explicitly legacy/reference diagnostics, but they no
longer determine warning codes. An undefined ratio caused by a zero run median
produces no warning. These warnings may indicate over-prediction,
under-prediction, different gene definitions, or an incomplete prediction
job. They never alter automatic pruning, model selection, scoring, source
inclusion, or scenario contents. This calculation is separate from the
biological source-singleton excess detector, which retains its existing
other-source rate comparison.

Reason and warning codes are machine-readable. Scenario rows continue to
include cumulative and incremental route/source/category counts,
source-normalised fractions, protections, and effective scenario settings.
Automatic diagnostics additionally expose the number and fraction of starting
selected models removed by high beyond medium and attribute that increment by
source.
