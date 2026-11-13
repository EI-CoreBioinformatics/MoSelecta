# Alternative-splice selection

> Technical reference for MoSelecta `0.2.16`.
> The normal product and diagnostic files are mapped in [Outputs](../outputs.md).

The [Configuration reference](../configuration.md) explains the option
families; [Scoring](scoring.md) defines the model-level and scaled score fields
used here. This document owns the exact post-freeze eligibility and redundancy
rules.

Alternative-splice selection runs by default when MoSelecta reaches the final
annotation stage. It consumes the frozen retained-primary set and cannot
change primary selection, scoring, pruning, automatic pruning, public gene
IDs, or primary `.1` identities. Use `--no-alternative-splice` to emit a
primary-only final annotation.

## Eligibility sequence

Diagnostic accounting keeps the following populations distinct:

```text
structural candidate
→ uniquely owned
→ fixed-rule eligible
→ support-policy eligible
→ external-evidence requirement met
→ prethreshold eligible
→ numerical-threshold survivor
→ greedy filtering
→ cap
→ accepted
```

`fixed_rule_potential_candidate` is retained only as a historical fixed-rule
diagnostic. It is not the current prethreshold population.

Ownership is structural: zero owners are rejected, one owner is unique, and
more than one owner is ambiguous. A single-exon primary can own only a
candidate with exactly two CDS exons under the 0.50/0.50 discovery-overlap
rule.

Structural rejection precedence is:

1. pairwise wobble-aware exact CDS equivalence;
2. same CDS intron chain with terminal-only CDS differences;
3. relaxed contiguous CDS-chain containment, applied bidirectionally.

For relaxed chain containment, the contained model must have at least two CDS
exons that map in transcript order to a contiguous slice of the other
CDS-exon chain. Every mapped splice junction and internal mapped exon boundary
must match. Only the outer CDS boundary of the first and last mapped exon may
be trimmed inside the corresponding containing exon. The comparison is
strand-aware. Exon skipping, intron retention, an added internal intron,
altered internal donor/acceptor usage, and same-span but internally
splice-distinct structures therefore remain distinct rather than being
rejected as containment.

The single-exon-primary/two-exon-candidate ownership route receives no
relaxed-chain containment exclusion. Normalized union-of-CDS-base subset
relationships remain diagnostics only and are calculated without double
counting duplicate or overlapping intervals. Historical exact exon-chain
containment also remains a separate unchanged diagnostic.

Frame must be fully calculated. Configured genome sequence is optional and
unevaluable states are neutral unless an actual defined concern is present.
Evaluated ambiguous splice motifs such as `NN` remain noncanonical. These
alternative-specific exclusions do not alter ordinary pruning policy.

## Presets and support policy

| Preset | Reciprocal CDS overlap | Frame concordance | Scaled model selection score | Unscaled primary-score ratio |
|---|---:|---:|---:|---:|
| permissive | 0.80 | 0.95 | 0.50 | 0.70 |
| intermediate | 0.85 | 0.95 | 0.60 | 0.70 |
| stringent | 0.85 | 0.95 | 0.60 | 0.80 |

All comparisons are inclusive. Explicit threshold options override only the
active result. The three canonical comparisons always use the table values,
the active support policy, external-evidence requirement, and configured
maximum-alternative cap.

The default `exact` support policy requires at least one of exact multisource
model corroboration, exact protein CDS-chain support, or one coherent exact
transcript association. Strong-partial evidence remains diagnostic only. The
`exact-plus-mosaic` policy additionally permits a qualifying full structural
mosaic route.

## Prediction-model source consensus

`--alternative-splice-min-model-support-sources` sets the minimum number of
distinct prediction sources for model-based support and defaults to 2. The
candidate source counts once; duplicate models from that source do not add
support. `--alternative-splice-min-model-support-source-fraction` optionally
sets a fraction greater than zero and at most one. Its denominator is the
number of distinct prediction sources in the valid/scored population for the
whole run. Protein and transcript evidence sources are excluded.

The effective requirement is the maximum of two, the configured absolute
minimum, and `ceil(prediction_source_count × configured_fraction)`. The same
effective count applies to exact whole-model support and mosaic support. Exact
support counts distinct sources containing the qualifying whole model. Mosaic
support requires every applicable structural component to meet the count;
different components may be supported by different source combinations.

`full_structural_mosaic` remains the historical K=2 diagnostic. Active mosaic
eligibility uses `alternative_splice_mosaic_source_requirement_met`, so an
exact whole-model count below the active requirement does not block a deeper
mosaic consensus. Raw exact and mosaic source counts and global-denominator
fractions remain visible even when the active threshold is not met.

## External-evidence requirements

`--alternative-splice-evidence-requirement` accepts:

- `none` (default);
- `protein-or-transcript`;
- `protein`;
- `transcript`;
- `protein-and-transcript`.

The required evidence inputs must be configured. Inputs that are configured
but yield no qualifying support are a valid biological outcome, not a
configuration error. Only exact protein and coherent exact transcript support
satisfy these requirements; multisource-model and mosaic routes do not.
Consequently, for every non-`none` requirement, `exact` and
`exact-plus-mosaic` produce the same accepted set under otherwise identical
configuration.

## Score-ratio availability

The unscaled candidate/primary ratio is eligible only when the retained primary
model selection score is greater than zero. Zero and negative primary scores
receive distinct unavailable statuses and conservative rejection. No
absolute-value, scaled-ratio, offset, or score-difference substitute is used.
With a positive primary, zero and negative candidate scores retain their
ordinary numeric ratio and fail the configured threshold naturally.

## Greedy processing and cap

Candidates are processed by the established model-selection comparator.
Every candidate continues to use `.1` as its ownership, evidence, overlap,
frame, and score anchor. Accepted alternatives participate only in redundancy
checks. Before acceptance, a candidate is compared with `.1` and all earlier
accepted alternatives using the same exact, intron-chain, and relaxed-chain
containment precedence.

Only accepted alternatives consume the cap. The default cap is 50 additional
transcripts and does not count `.1`. Processing continues after the cap so all
candidates retain diagnostics. Active accepted representatives receive `.2+`
ordinals; canonical comparisons report internal model identities only.

## Outputs

`<prefix>.final_annotation.gff3` is the final annotation. It
contains every frozen public gene and primary `.1`, plus accepted `.2+`
models when alternatives are enabled. `--no-alternative-splice` keeps it
primary-only while preserving all gene/`.1` identities and features.

`--alternative-splice-report` writes the candidate, gene, and run-summary TSV
files. It may be combined with `--no-alternative-splice`: the selector is then
evaluated diagnostically but its result is not applied to the final
annotation. The run metadata distinguishes `selector_evaluated` from
`alternatives_emitted`; emission means application of the result, even when
zero candidates happened to be accepted.

The optional `<prefix>.alternative_splice_candidate_audit.gff3` is requested
with `--write-output alternative_splice_candidate_audit_gff3` or configured
with `--alternative-splice-candidate-audit-gff3`. It contains every uniquely
owned structural candidate, including fixed-rule, support, evidence,
threshold, greedy, and cap rejections. Zero-owner and ambiguous-owner
assessments remain TSV-only.

Audit candidates are top-level synthetic mRNAs without a public-gene
`Parent`. Their exon/CDS/UTR children use the synthetic audit transcript as
`Parent`. Public gene and `.1` associations are attributes only. The audit is
a diagnostic browser track, not a retained gene annotation. Its IDs remain
stable across threshold presets and support policies whenever the frozen
primary and uniquely owned candidate populations are unchanged.

The GFF support-route attribute records factual routes present on a model,
independently of policy. Policy eligibility and the active decision state
whether those routes qualified.

Run-summary ratio-unavailability fields ending in `_count` count unique
candidate models with at least one matching assessment. Fields ending in
`_assessment_count` count candidate-primary assessment rows; an ambiguous
candidate can contribute several such rows.

Native greedy nonmonotonicity is compared once for each adjacent preset pair:
permissive-to-intermediate and intermediate-to-stringent. The diagnostic does
not impose monotonic closure or change selector decisions.
