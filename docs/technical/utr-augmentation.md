# Post-selection UTR augmentation

> Technical reference for MoSelecta `0.2.16`.
> UTR reconstruction writes a separate annotation after CDS-based integration;
> it does not replace the ordinary final annotation. See [Outputs](../outputs.md).

The [Configuration reference](../configuration.md) summarizes activation and
option families. This document owns evidence eligibility, extension limits,
proposal selection and reconstruction semantics.

UTR splice-junction validation is a default post-selection quality stage. It
validates inherited source-model UTRs in the ordinary final annotation and
candidate UTR sides in the optional augmentation pathway. UTR augmentation
itself is an optional, default-off final stage enabled with
`--add-utrs`. It operates only after primary pruning, stable public identifier
assignment, and alternative-splice acceptance. It does not participate in CDS
scoring, ranking, de-overlap, backfill, pruning, or alternative-splice
selection. Independent CDS terminus refinement may run immediately before it;
see `cds_terminus_refinement.md`.

## Command and products

```bash
moselecta ... --add-utrs --out-prefix integrated
```

The stage and shared validator write these additional products according to the
selected output profile:

```text
integrated.final_annotation.with_utrs.gff3
integrated.utr_augmentation.audit.tsv
integrated.utr_augmentation.summary.tsv
integrated.utr_junction_validation.summary.tsv
integrated.utr_junction_validation.audit.tsv
```

The ordinary `integrated.final_annotation.gff3` remains a separate product.
When default terminus refinement succeeds it contains refined CDS geometry but
no newly selected UTR features. Named UTR output keys and path options are:

| Product | Named output key | Path option |
|---|---|---|
| UTR-augmented GFF3 | `utr_augmented_final_annotation_gff3` | `--utr-augmented-final-annotation-gff3` |
| Per-transcript audit | `utr_augmentation_audit_tsv` | `--utr-augmentation-audit-tsv` |
| Run summary | `utr_augmentation_summary_tsv` | `--utr-augmentation-summary-tsv` |
| Junction validation summary | `utr_junction_validation_summary_tsv` | `--utr-junction-validation-summary-tsv` |
| Junction validation audit | `utr_junction_validation_audit_tsv` | `--utr-junction-validation-audit-tsv` |

## Invariants

Model selection is authoritative. Enabling `--add-utrs` cannot alter:

- retained genes or transcripts;
- public gene or transcript IDs and `.1/.2+` numbering;
- transcript membership within genes;
- any score, rank, pruning decision, or alternative-splice decision.

CDS coordinates and phases may already have been changed by the independent,
default-enabled terminus-refinement stage. The UTR switch itself does not
enable, disable, score, or select that refinement.

Exactly one augmented form is emitted for each final public transcript. It can
receive both UTR sides, only the biological 5-prime side, only the biological
3-prime side, or neither. No UTR-only alternative isoforms are created.

An explicitly incomplete CDS terminus is not eligible as an ordinary UTR
anchor. `start_codon_status == missing` excludes the biological 5-prime side;
`stop_codon_status == missing` excludes the biological 3-prime side. This
ordinary filter remains in place during limit learning, so incomplete sides
never train the distributions. A separate completion-only candidate mode may
use the same supported structure to repair the CDS. After success, residual
sequence is re-evaluated as ordinary UTR against the already frozen limits.

## Evidence catalogues

During transcript-evidence processing, MoSelecta retains a compact exact
endpoint-paired exon-chain catalogue before releasing the much larger raw
transcript objects. It preserves strand, exon and intron chains, paired
terminal endpoints, source-specific member counts, representative evidence
IDs, and complete source-exon geometry.

All valid prediction models are also eligible UTR evidence, whether selected
or unselected. Their normalized exon structure is used; when exon records are
absent, the existing normalized CDS/UTR child representation is used. A source
prediction can contribute sequence that it labelled CDS, but the selected
MoSelecta CDS remains authoritative. The audit distinguishes selected from
unselected prediction evidence and source-CDS from source-non-CDS extension.

## Compatibility and hard eligibility

Evidence must satisfy full CDS-chain containment, strand-aware terminal
compatibility, or strict single-exon containment. Coordinate overlap alone is
not sufficient. Biological 5-prime and 3-prime logic is strand-aware.

Before support is counted, each candidate side must pass:

1. an explicitly complete/non-missing CDS anchor on that biological side;
2. multi-final-gene exclusion;
3. minimum complete source/evidence exon length;
4. minimum and maximum newly introduced UTR intron lengths;
5. the shared strand-aware splice-motif policy when sequence is evaluable;
6. the side-specific UTR exon-count limit;
7. the learned or explicitly configured exonic extension limit;
8. the learned or explicitly configured genomic extension limit.

The shared validator applies to every materialized candidate side regardless
of evidence origin. Prediction-derived and transcript-derived candidates are
therefore assessed at the same point, before support pooling, endpoint
consensus, automatic extension-limit learning, or winner selection. Rejecting
a transcript-derived UTR side does not remove that transcript alignment from
its separate CDS-scoring, rescue, pruning, or alternative-transcript roles.

The accepted motif definitions are the same `GT-AG`, `GC-AG`, and `AT-AC`
definitions used by genome-derived CDS integrity assessment. MoSelecta does
not maintain a separate UTR motif list. The
`--sequence-penalty-noncanonical-splice` option remains a CDS model-selection
penalty weight and does not change UTR validation policy.

At the default, any UTR-associated intron shorter than 20 bp rejects the whole
biological side. The threshold is inclusive: 20 bp passes the length check.
Introns of 1-5 bp use the same `short_utr_intron` reason as every other
sub-20 bp intron. The audit records the measured length and effective minimum.

The per-transcript audit reports blocked anchors as
`five_prime_anchor_block_reason=missing_start_codon` or
`three_prime_anchor_block_reason=missing_stop_codon`, with quality status
`blocked_incomplete_cds_anchor`. Run-level rejection counts use
`incomplete_cds_anchor_missing_start_codon` and
`incomplete_cds_anchor_missing_stop_codon`.

The 5-bp exon rule applies to the complete source exon, not only the clipped
piece that becomes UTR. Thus a 300-bp evidence exon leaving a 4-bp UTR segment
passes, while a genuinely new 4-bp evidence exon fails.

Evidence structurally compatible with transcripts from more than one final
MoSelecta gene is excluded before it can vote. Sharing among `.1/.2+`
transcripts in the same final gene is allowed. Geometric overlap with another
gene is not itself a reason for exclusion.

The final-gene separation is also protected at the *structural-feature* level.
Before a side proposal may vote or enter automatic-limit learning, MoSelecta
materializes that side against its authoritative CDS and rejects any newly
introduced exact exon, intron, or splice-site coordinate already belonging to
another same-strand final gene's frozen CDS or ordinary source-exon structure.
The generic splice-site rule matches MoSelecta's multi-exon locus grouping:
one exact shared site is a structural connection even when full introns differ.
Exact exons additionally protect terminal and single-exon cases. Mere exon or
gene-span overlap, partial containment, and opposite-strand coincidence do not
trigger this safeguard.

After provisional selections, a global pass detects exact links created only
when two different final genes independently choose optional UTRs. Both
conflicting side proposals are invalidated simultaneously; each side falls
through to its next ranked eligible proposal, and validation repeats until no
new cross-final-gene link remains. If no safe proposal remains, that side is
left without a UTR. This never merges or renumbers frozen genes, never changes
their CDSs, and never assigns ownership of an ambiguous UTR to the first gene
processed or merely to the side with the larger support score. Same-final-gene
alternatives may still share UTR structure.

## Automatic extension limits

Limits are learned once per run from candidate sides that passed all
non-extension checks. Biological sides and geometries are independent:

- 5-prime total newly added exonic sequence;
- 3-prime total newly added exonic sequence;
- 5-prime genomic boundary-to-endpoint distance;
- 3-prime genomic boundary-to-endpoint distance.

For each side and geometry, the automatic limit is the next whole kilobase at
or above the larger of the all-candidate P99.5 and the multi-source P99.9. The
multi-source tail is used only with at least 1,000 observations; otherwise the
all-candidate P99.9 is used. A sparse side pools both biological sides for the
same geometry. If the pooled population remains sparse, exonic extension uses
the 20-kb fail-safe and genomic extension remains uncapped apart from the
individual-intron safeguard. The summary records the effective value, source,
and any sparse genomic warning.

## Structure and endpoint selection

Candidate sides are first grouped by a canonical side-specific splice
`structure_key`; terminal-coordinate variation is normalized within a
structure, while incompatible splice topology remains separate. Observations
are pooled by biological source across all endpoints.

For a transcript source with `N` members, its structure/window weight is:

```text
1 + min(3, log2(N))
```

A prediction source contributes at most one vote regardless of the number of
models or endpoints. Structures rank by effective support, independent source
breadth, multi-member transcript-source count, and deterministic lexical key.

Endpoint selection occurs only inside the winning structure. MoSelecta builds
bounded windows from observed coordinates using three times the median positive
adjacent gap, clamped to 20-250 bp. A window is local consensus only when its
source-pooled support is strictly greater than half the structure support. A
single endpoint is also local consensus. Otherwise all observed endpoints are
used in the diffuse fallback.

Source influence is normalized across its endpoints and the selected endpoint
is an actually observed weighted medoid. Strong disjoint alternative endpoint
modes are audited but do not create extra transcripts. The two biological
sides are selected independently. Shared-evidence coherence can resolve only
an exact numerical tie on both sides; it cannot replace a stronger side-specific
choice.

## Reconstruction and GFF3 behavior

The augmented writer reconstructs exon, `five_prime_UTR`, `three_prime_UTR`,
and CDS children from the authoritative CDS plus the selected side proposals.
Transcript spans become the union of these parts; gene spans become the union
of all augmented transcript spans in that public gene. A transcript with no
selected UTR is emitted with a CDS-derived exon chain.

Compact selection, protein/transcript evidence and sequence-integrity
attributes on each augmented mRNA are copied from the immutable public
annotation metadata for its authoritative selected CDS model. They are kept
for interpretation and provenance; they are not scores of the selected UTR
extension and are never recomputed from the augmented transcript span.

`--final-gff3-cds-only-structure` still makes the ordinary final annotation
CDS-only. The separately requested `with_utrs` annotation deliberately retains
its selected UTR structure while matching the active ordinary-writer CDS phase
semantics.

The GFF3 contains only concise augmentation attributes. The detailed audit TSV
has one row per final transcript with side-specific structure, endpoint method,
support, source breadth, compact evidence provenance, prediction provenance,
competing-structure and multimodality flags, effective limits, extension
lengths, final-recheck status, cross-final-gene rejection reasons, the initially
selected structure/endpoint, and any conflict fallback. The run summary counts
hard candidate rejections, selected-UTR conflicts, fallback/no-alternative
sides, and affected genes/transcripts. Filtering against frozen other-gene
structure occurs before extension-limit learning, so learned limits may differ
from earlier runs if invalid candidate tails had influenced them.

## Inherited UTR validation and missing sequence

The ordinary final writer validates complete biological 5-prime and 3-prime
UTR sides independently after CDS selection is frozen. A definite geometry,
minimum-length, or resolved motif failure suppresses only the affected side.
The opposite valid side, authoritative CDS coordinates and phases, public IDs,
scores, and selection membership are retained. Exon and parent spans are
reconstructed from the surviving structure; MoSelecta neither joins separated
exons nor shifts a splice boundary.

Without genome FASTA, length and geometry checks still run. Motif status is
`not_fully_evaluated` with reason `no_genome_fasta`, and the side is retained
unless another definite violation exists. Missing contigs, out-of-range
contexts, and ambiguous bases are handled in the same conservative way. An
unspliced side has no motif to test and is recorded as not applicable rather
than as motif-passed.

The validator runs before every eligible transcript- or prediction-derived
candidate can contribute to support pooling, consensus, automatic extension
limits or side selection. A rejected side therefore cannot influence any UTR
decision. The same transcript evidence remains available to the separate CDS
evidence and model-selection stages.

Every profile receives a compact summary with assessed, accepted, rejected and
unevaluable side totals plus junction, motif, side and reason counters. The
detailed junction audit is enabled by full/debug profiles or explicit output
selection. It contains exceptions only and is streamed while validation runs;
accepted sides, passing junctions and routine `no_genome_fasta` outcomes are
not retained as detail rows. Side-level and intron-level reason columns are
separate. Rejected inherited sides use `suppress_source_side`; rejected
augmentation candidates use `exclude_candidate_side`.

## Advanced controls

These settings are shown by `moselecta --full-help`:

| Option | Default |
|---|---:|
| `--utr-min-complete-exon-bp` | 5; 0 disables |
| `--utr-min-intron-bp` | inherit active CDS minimum; 0 disables |
| `--utr-max-intron-bp` | inherit resolved long-CDS-intron cutoff; 0 disables |
| `--utr-max-five-prime-exons` | 4; 0 disables |
| `--utr-max-three-prime-exons` | 2; 0 disables |
| `--utr-exonic-failsafe-bp` | 20000; 0 disables |
| `--utr-five-prime-exonic-override-bp` | automatic; 0 disables |
| `--utr-three-prime-exonic-override-bp` | automatic; 0 disables |
| `--utr-five-prime-genomic-override-bp` | automatic; 0 disables |
| `--utr-three-prime-genomic-override-bp` | automatic; 0 disables |
| `--utr-member-depth-cap` | 3 |
| `--utr-endpoint-window-min-bp` | 20 |
| `--utr-endpoint-window-max-bp` | 250 |
| `--utr-local-consensus-min-fraction` | 0.5, applied strictly greater than |
| `--utr-limit-min-observations` | 1000 |

These controls are intended for explicit project policies and validation, not
routine tuning. An observed/supported UTR endpoint is not necessarily a known
biological transcription start site or polyadenylation site.
