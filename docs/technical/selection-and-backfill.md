# Selection, conflict resolution and backfill

> Technical reference for MoSelecta `0.2.16`. The accepted implementation in `loci.py`, `selection.py`, `scoring.py` and the strict-replacement tests is authoritative.

## Purpose and stage boundaries

Selection converts scored candidates into a structurally compatible set. It is not one “take the highest score” operation. MoSelecta first chooses one candidate in each conservative first-pass competition group, resolves overlaps among those provisional winners, then revisits other compatible candidates through backfill. A strict route can replace a narrowly defined set of blockers when a corroborated multi-exon candidate is better than every blocker.

The principal states are:

| State | Meaning |
|---|---|
| Eligible/scored candidate | Passed pre-locus exclusion and has a model-selection score |
| Sublocus winner | Highest candidate in one first-pass sublocus; genuinely provisional |
| Monosublocus selected | Sublocus winner retained by greedy winner de-overlap |
| Suppressed | Sublocus winner blocked by a higher-ordered selected winner |
| Backfilled | Non-winner candidate accepted because it does not conflict with current selected models |
| Replaced holder | Monosublocus-selected blocker removed atomically by strict replacement |
| Pruned | Selected/backfilled model removed later by final pruning |
| Retained primary | Model surviving pruning and frozen as a public `.1` transcript |

Exclusion, suppression and pruning are distinct. Exclusion precedes scoring/selection, suppression resolves winner conflicts, and pruning is a later final filter. Pruning does not trigger another backfill pass.

## Locus hierarchy

### Superloci

Candidates are first grouped by chromosome and strand. Sorted CDS spans belong to the same superlocus when they overlap by at least one base; connectivity is transitive. A superlocus therefore bounds potentially related or conflicting candidates and is the scope of backfill. It is not asserted to be a final biological gene.

### First-pass subloci

Within each superlocus, an undirected graph defines conservative competition groups:

- two multi-exon models connect when they share at least one exact CDS splice-site coordinate;
- two single-exon models connect when the reciprocal CDS-span overlap on **both** models reaches `--single-exon-min-reciprocal-overlap` (default 0.50);
- single-exon and multi-exon models do not connect at this stage.

Connected components form subloci. Because connectivity is transitive, a component can contain endpoints that are not directly pairwise related; `splice_site_chain` and `single_exon_overlap_chain` flags expose those cases. A sublocus is the scoring/ranking context and supplies one provisional winner.

The first-pass split deliberately avoids forcing every spatial overlap into one competition. Later monosublocus and backfill conflict rules use CDS-exon overlap and can still connect or block candidates from different subloci.

## Deterministic comparators

### Comparator used by selection

Sublocus winner choice, monosublocus member ordering and backfill candidate ordering call the same `selection.py` comparator:

1. higher unscaled `model_selection_score`;
2. higher `model_completeness_score`;
3. lower `intrinsic_penalty_score`;
4. higher `consensus_similarity_score`;
5. higher `structural_support_score`;
6. higher `source_prior_score`;
7. longer CDS;
8. lexical `output_transcript_id`.

The comparator is fully deterministic. The score is not rounded, and the score-tie epsilon is not applied inside this selector key.

### Reported ranks and the shared comparator

`model_selection_score_rank` groups scores within the configured epsilon (`1e-6` by default) and uses competition numbering (`1, 1, 3`). A near-but-unequal pair can therefore share a score rank.

`model_selection_order_rank` is the unique `1..N` order over the same scored canonical-sublocus population. It calls the same authoritative helper as the selector comparator above, compares raw scores exactly and reaches secondary fields only for exact score equality. Thus sorting candidates by this rank reproduces the selector order, and rank 1 is the first-pass sublocus winner. The rank is assigned before suppression, backfill and pruning; it does not guarantee final retention and is not comparable across unrelated subloci. The rank reports the shared order rather than serving as an input to selection.

## Sublocus winner selection

Subloci are traversed in stable genomic/ID order. All scored members of one sublocus are sorted with the selection comparator; the first becomes its sublocus winner. A sublocus with no scored member contributes no winner.

The winner is provisional because:

- a winner from another sublocus may conflict with it in monosublocus de-overlap;
- it can be suppressed by greedy overlap resolution;
- a retained winner may later be atomically replaced during strict backfill;
- a selected winner can subsequently be removed by final pruning.

The score itself is not provisional and is not recalculated after each membership decision.

## Monosublocus winner de-overlap

### Overlap geometry

For two same-strand sublocus winners, overlap is measured across the union of their CDS exon intervals, not their outer CDS spans:

```text
overlap_bp = total CDS-exon intersection
fraction_A = overlap_bp / CDS length A
fraction_B = overlap_bp / CDS length B
max_fraction = max(fraction_A, fraction_B)
```

An edge forms when `overlap_bp > 0` and `max_fraction` reaches `--monosublocus-min-cds-exon-overlap-fraction` (default 0.10). The threshold is non-reciprocal: substantial coverage of either model is sufficient. Only same-chromosome, same-strand winners are compared.

Connected components of this edge graph are called monosublocus holders. A component need not be a clique; the `not_all_pairwise_overlapping` diagnostic records transitive components.

### Greedy retention and suppression

Members of each component are sorted with the selection comparator. They are visited in that order:

1. keep a model if it has no edge to any already kept model;
2. otherwise suppress it;
3. identify the suppressor as the conflicting retained model with the largest `max_overlap_fraction`, using the selection comparator to break equal overlap.

This is greedy set construction, not a global optimization of score sums. More than one mutually compatible model can remain in a transitive component. Suppression records the suppressor, overlap bases, both directional fractions and the threshold-bearing maximum fraction.

Suppressed winners are permanently excluded from ordinary backfill consideration. This prevents immediate reintroduction of a conflict that the de-overlap stage just resolved.

## Ordinary backfill

### Candidate population and ordering

Backfill is evaluated independently within each superlocus. It begins from the models retained by monosublocus de-overlap. Candidate population includes every other scored model in that superlocus except:

- a model already selected; and
- a sublocus winner suppressed by monosublocus de-overlap.

Candidates are sorted by the selection comparator. Acceptance is greedy: every accepted candidate is immediately added to the current selected set and can block a later candidate.

Backfill is useful because first-pass selection intentionally chooses only one candidate per conservative sublocus. A non-winner can still represent a compatible additional coding locus or structure after higher-priority conflicts have been resolved.

### Conflict rules

For a candidate and every currently selected model in the same superlocus, conflict exists when either condition is true:

1. both are multi-exon and share at least one exact CDS splice-site coordinate; or
2. CDS-exon overlap is positive and `max(overlap / candidate CDS length, overlap / selected CDS length)` reaches the same default 0.10 monosublocus threshold.

The shared-splice rule can reject models even when CDS-exon overlap fraction alone is below threshold. Conversely, the overlap rule applies to any exon-count combination. Conflicts are reported in deterministic priority order: combined shared-splice-and-overlap, shared splice only, then overlap only; larger splice counts/overlap and blocker ID break ties.

If there is no conflict, the candidate is accepted as `final_selection_stage=backfill`. If there is a conflict, ordinary backfill rejects it and records the strongest blocker in the main ledger. All conflicts remain available to the strict replacement test.

## Strict score-improving replacement

This route is implemented as strict R1 and its multi-exon-blocker extension R1-ME. Functional requirements, rather than shorthand, define the rule.

### Complete eligibility rule

A blocked backfill candidate is eligible only when all of the following hold:

1. `--strict-r1-backfill-replacement` is enabled (default true).
2. The candidate is multi-exon.
3. `corroborating_other_source_count >= 1` for the candidate.
4. The complete conflict set contains exactly one or two blockers.
5. Every blocker is an existing model selected at the monosublocus stage. A model accepted by earlier backfill is not an eligible blocker.
6. Every blocker has `corroborating_other_source_count < 1`.
7. Under original R1, every blocker is single-exon. With `--strict-r1-allow-multiexon-blockers` enabled (default true when strict R1 is enabled), qualifying multi-exon holder winners are also allowed.
8. The candidate's unscaled `model_selection_score` is **strictly greater than the score of each blocker individually**.

The last rule is:

```text
candidate_score > blocker_1_score
candidate_score > blocker_2_score  # when a second blocker exists
```

MoSelecta does not compare the candidate against the sum, mean or scaled scores of blockers. Equality fails the strict test. Corroboration eligibility uses the source-count field, not the weighted corroboration contribution, so setting the corroboration aggregate weight to zero does not invent or erase the factual count.

### Atomicity and state change

Eligibility is evaluated against the full conflict set before state mutation. If any requirement fails, no blocker is removed and the candidate remains rejected. If all requirements pass:

- all one or two blockers are removed from the current selected set;
- the candidate is inserted once;
- final membership and per-model provenance are updated together;
- the replacement report records blocker classes, scores, corroboration counts and every conflict geometry.

This atomic behavior prevents partial replacement when one blocker qualifies and another does not. Replaced holder IDs are removed from the final selected rows after the superlocus scan.

### Rejection reasons

The implementation records precise reasons, including:

- `candidate_not_multi_exon`;
- `candidate_corroboration_count_lt_1`;
- blocker count not one or two / greater than two;
- blocker not found;
- blocker not single-exon when the extension is disabled;
- blocker not selected at the monosublocus stage;
- blocker corroboration count at least one;
- candidate score not strictly greater than every blocker.

These are algorithmic reasons, not biological verdicts.

### Control consistency

`--strict-r1-allow-multiexon-blockers` requires strict replacement. A contradictory effective configuration is rejected. The configured/effective states are captured in `run_configuration.json` and replacement summaries.

## Conceptual trace

This example is conceptual and intentionally uses order labels rather than invented MoSelecta scores.

1. Models A and B share a CDS splice site and enter one multi-exon sublocus; A is ordered first and becomes its provisional winner.
2. A second sublocus supplies winner C. A and C overlap across at least 10% of one model's CDS exons, so they share a monosublocus edge. If A is ordered first, A is retained and C is suppressed.
3. B was not a sublocus winner and was not suppressed, so it is a backfill candidate. If B has neither a shared-splice nor threshold overlap conflict with the current selected set, it is accepted; otherwise ordinary backfill rejects it.
4. If B is a corroborated multi-exon model blocked only by one or two uncorroborated monosublocus holder winners and its unscaled score is strictly greater than each, strict replacement may remove those blockers and accept B atomically.
5. The selected set then enters pruning. No new backfill occurs if pruning removes a member.

The maintained two-file smoke fixture is a simpler validated trace: two multi-exon candidates occupy one sublocus; one becomes the winner and the other is later rejected by backfill because it shares a CDS splice site and exceeds the CDS-exon overlap threshold.

## Interaction with final pruning

Selection/backfill establishes the unpruned selected population written to `final_selected_models.*`. Final pruning then applies singleton, low-score, sequence-integrity, opposite-strand-shadow, rescue and protection rules to that immutable starting population.

Important consequences:

- a high selection score does not make a model immune to every pruning rule;
- scaled thresholds used by pruning do not replace unscaled score ordering or strict comparison;
- backfilled models can be in a pruning scope distinct from monosublocus-selected models;
- pruning removal does not cause a suppressed/rejected model to be reconsidered;
- retained primaries and public identifiers are frozen only after pruning.

Exact fixed and automatic rules are in [Automatic pruning](automatic-pruning.md).

## Configuration controls

| Option | Default | Direct effect |
|---|---:|---|
| `--single-exon-min-reciprocal-overlap` | 0.50 | Single-single first-pass sublocus edges |
| `--monosublocus-min-cds-exon-overlap-fraction` | 0.10 | Non-reciprocal winner edges and ordinary backfill overlap conflicts |
| `--strict-r1-backfill-replacement` | enabled | Enables strict one-/two-blocker replacement testing |
| `--strict-r1-allow-multiexon-blockers` | enabled | Extends eligible holder blockers beyond single-exon models |
| Score/evidence weights | See scoring reference | Change model-selection scores before this stage |
| `--model-selection-tie-epsilon` | `1e-6` | Changes reported score-rank tie groups, not the selector comparator key |
| `--stop-after` | `final` | Can stop at scored candidates, monosublocus or backfill for diagnostics |

Changing a scoring weight can change ordering indirectly; it does not alter the geometric conflict predicates. Changing the overlap threshold directly changes monosublocus edges and backfill overlap conflicts, and can therefore change membership even if scores are unchanged.

## Audit outputs

Inspect in stage order:

| Output | Question answered |
|---|---|
| `sublocus_members.tsv` / relationship ledgers | Why were candidates grouped together? |
| `model_selection_scores.tsv` | What components, score and report ranks were assigned? |
| `sublocus_winners.tsv` | Which first-pass candidate won each sublocus? |
| `monosubloci.tsv`, member and edge TSVs | Which winner-overlap component formed and what geometry linked it? |
| `selected_models.tsv` | Which winners survived greedy de-overlap? |
| `suppressed_monosublocus_models.tsv` | Which winner was suppressed by which retained winner? |
| `backfill_candidates.tsv` | In what order was each eligible non-winner considered? |
| `backfilled_models.tsv` | Which candidates were accepted normally or by replacement? |
| `backfill_rejected_models.tsv` | Which blocker and conflict rejected the candidate? |
| strict replacement detail/summary TSVs | Which one-/two-blocker cases were accepted or rejected, and why? |
| `final_selected_models.*` | What unpruned membership left selection/backfill? |
| pruned retained/removed ledgers | What membership survived final pruning? |

Output-profile membership and exact fields are maintained in [Outputs](../outputs.md). `run_configuration.json` preserves effective thresholds, weights and strict-replacement states.

## Determinism and limitations

- Genomic traversal, component formation, candidate ordering, conflict prioritization and final row sorting use stable keys.
- Greedy algorithms are order-dependent by design; deterministic ordering makes the result reproducible but does not make it a global optimum.
- CDS-exon overlap is a computational conflict predicate, not a claim that every same-strand biological overlap represents one gene.
- First-pass subloci and monosubloci serve different purposes and must not be treated as interchangeable gene definitions.
- Structural corroboration counts prediction sources; correlated methods may not be independent biological evidence.
- Strict replacement is intentionally narrow. A candidate that seems preferable by another informal criterion remains rejected if any implemented safeguard fails.
