# Post-selection CDS terminus refinement

> Technical reference for MoSelecta `0.2.16`.
> Refinement is applied only after selection and does not rescore models or
> change membership. See [Outputs](../outputs.md) for the audit products.

The [Configuration reference](../configuration.md) summarizes activation and
bounded-search controls. This document owns the exact proposal, safety and
rollback rules.

MoSelecta automatically attempts conservative repair of final retained CDS
models whose genome-derived start or stop status is explicitly `missing`.
Refinement requires `--genome-fasta`, runs after selection/pruning and public
transcript construction, and is independent of `--add-utrs`.

Before locus construction, a separate bounded representation-normalisation
stage may already include an immediately implied terminal stop by one unique
`+3/+2/+1/-1/-2` local boundary correction. That early stage uses only the
model's established frame and simple terminal exon geometry; it is not this
post-selection evidence-supported refinement. Successfully canonicalised stops
subsequently evaluate as `included_in_cds` and therefore do not enter the
missing-stop refinement route. Structurally nonlocal cases remain unchanged for
the existing downstream policy.

The original selected model remains the immutable selection provenance. A
successful repair creates a consistent cloned final-annotation model with new
CDS parts, geometry caches, phases, sequence-integrity values, exon/transcript
bounds, and emitted children. It does not rescore, re-prune, reselect, renumber,
or alter O5a policy.

## Structure and search rules

Compatible structures may come from other input models or supplied transcript
evidence. They are grouped and ranked by the existing source-normalized support
rules without inspecting codons. Only the preferred eligible structure is
searched; codon failure does not fall through to a weaker structure.

For a missing start, MoSelecta first searches the preferred upstream spliced
path. It walks the current frame to the nearest preceding in-frame stop, the
selected structure boundary, or a contiguous `NN+` barrier. Within that finite
open segment it chooses the first `ATG` in 5-prime-to-3-prime order: the
most-upstream valid ATG on the preferred path. Structure end is a valid
boundary and an actual preceding stop is not required.

The smaller fallback searches inside the existing CDS for the nearest
downstream in-frame `ATG` and does not require external structure. It runs only
when no usable upstream structure exists or the preferred structure has no
in-frame upstream ATG. Once a true upstream ATG is identified, rejection by a
safety rule fails that completion; MoSelecta does not substitute a nearer
internal ATG or trim farther downstream.

For a missing stop, a compatible downstream structure is required. MoSelecta
chooses the first in-frame `TAA`, `TAG`, or `TGA` along that spliced path and
includes the codon in the CDS. Both strands, multi-exon paths, and codons that
cross exon junctions are supported. Models missing both ends are searched on
both sides. Start and stop proposals are assembled and the complete CDS is
validated before any optional residual UTR decision.

## Safety and rollback

The preferred upstream structure, rather than an absolute base window, bounds
start discovery. The former bounded upstream-start distance option has been
removed; discovery follows the supported ORF/structure. Default search windows
remain 30 CDS bases for the downstream-start fallback and 100 spliced bases
downstream for a stop. Extensions may add at most 100% of the original selected
CDS length by default; downstream-start correction may remove at most 33%.
Discovery reaches the natural structure boundary and the selected natural
codon is then explicitly rejected if its relative guard fails.

A single isolated `N` in newly coding sequence is tolerated. `NN` or a longer
run bounds upstream traversal; other unsupported ambiguity or ambiguity in the
proposed start/stop codon blocks that proposal. Intronic sequence is irrelevant. New CDS introns
must remain canonical under the existing genetic-code-1 integrity machinery.

Every proposed final geometry is compared with the frozen final-gene
population using MoSelecta's canonical sublocus relationship predicate. Before
commit, the complete population of refined and unchanged models is checked a
second time. If refinements jointly create a new cross-gene relationship, all
refinements participating in that relationship are rolled back and the
population is revalidated. Candidate bridge detection includes complete final
genes as well as incomplete targets. Relationship checks use indexed
same-chromosome/strand spatial neighbours of replacement transcripts and cache
the corresponding frozen gene-pair result; they do not scan all gene pairs. A
new shared splice-site relationship or
qualifying single-exon reciprocal overlap causes rollback, as do internal
stops, inconsistent phase, invalid start/stop status, and other failed
integrity checks. Evidence cannot override an internal stop.

An additional refinement-only guard considers newly introduced coding bases.
They must not overlap another frozen final MoSelecta gene's CDS on the same
strand. Existing CDS overlap is not reconsidered, opposite-strand overlap and
same-gene alternatives are allowed, and UTR-only overlap is irrelevant. A
second deterministic population pass applies the same rule to refined-vs-
refined genes and rolls back every contributing refinement.

## UTR interaction

Without `--add-utrs`, only incomplete final transcripts are inspected and no
ordinary UTR candidate-selection pass is run. The completion catalogue is
deferred until the retained final population is known and is skipped entirely
when no retained CDS is incomplete. Raw transcript evidence is released after
the ordinary evidence stage as before. If completion is required, MoSelecta
re-reads it once with the same canonical parser to build the exact
endpoint-paired catalogue, then immediately releases those temporary raw
records. This exchanges bounded additional input I/O for avoiding long-lived
raw-transcript memory. Successful refined CDS geometry is written to the
ordinary final annotation, with no newly inferred UTR features.

With `--add-utrs`, limits are learned first from complete ordinary anchors.
Originally incomplete sides never contribute observations. A CDS repair is
committed solely from CDS-specific checks. The downstream UTR stage then
evaluates residual sequence independently; a residual UTR that fails its
length, exon, or structural checks is omitted without rolling back the valid
CDS repair. It may choose another valid UTR structure under the ordinary UTR
policy.

## Controls and products

`--refine-cds-termini` is enabled by default when a genome FASTA is supplied;
`--no-refine-cds-termini` disables it. Advanced controls shown by
`--full-help` are:

| Option | Default |
|---|---:|
| `--cds-refinement-downstream-start-search-bp` | 30 |
| `--cds-refinement-downstream-stop-search-bp` | 100 |
| `--cds-refinement-max-relative-extension` | 1.0 |
| `--cds-refinement-max-relative-truncation` | 0.33 |

Standard runs write
`<prefix>.cds_terminus_refinement.audit.tsv`; the compact summary
`<prefix>.cds_terminus_refinement.summary.tsv` is included in minimal output
and is written with an explicit execution status even when refinement is
disabled, prerequisites are unavailable, or no retained CDS is incomplete.
The audit records original/refined geometry and length, candidate provenance,
direction, search limit, ORF-boundary state, codon positions, relative change,
junctions, ambiguity, short-CDS threshold state, structural/integrity checks,
residual UTR geometry, encountered rejection events, final result, and rollback
reason.

### Audit and summary accounting

A **side** is one explicitly incomplete biological CDS terminus on one retained
transcript: `five_prime` for a missing start or `three_prime` for a missing
stop. Each side appears exactly once in the audit and has exactly one
`final_outcome`. A refined side has outcome `refined`; otherwise the outcome is
its final `failure_reason` after all permitted fallback and population rollback
processing.

An **affected transcript** for a reason is a transcript with at least one side
whose final failure has that reason. Transcript counts are deduplicated by
public transcript ID, so a transcript with two failed sides contributes two to
`failure_counts` but one to `failure_transcript_counts` for the same reason.

An **encountered rejection event** is a rejected candidate path, proposal guard,
or proposal-validation attempt observed while evaluating a side. It is not a
terminal outcome: a side can encounter one or more events and later refine via
a permitted fallback. The audit columns `encountered_rejection_event_count` and
`encountered_rejection_counts` preserve these per-side events. Summary event,
side, and transcript maps report respectively the number of events, the number
of sides encountering at least one event of that reason, and the number of
deduplicated transcripts encountering it.

Summary units are explicit:

| Field | Unit and meaning |
|---|---|
| `audit_sides` | Sides represented in the audit. |
| `incomplete_five_prime_sides`, `incomplete_three_prime_sides` | Input incomplete sides by biological terminus. |
| `transcripts_with_two_incomplete_sides` | Transcripts initially missing both termini. |
| `final_refined_sides`, `final_failed_sides` | Final terminal side outcomes. Their sum equals `audit_sides`. |
| `refined_sides_upstream_start`, `refined_sides_downstream_start_fallback`, `refined_sides_downstream_stop` | Final refined sides by completion route. |
| `refined_transcripts`, `unchanged_transcripts` | Deduplicated transcripts with at least one final refined side, or with none. |
| `two_incomplete_sides_*_transcripts` | Deduplicated two-incomplete-side transcripts by final combination. |
| `refined_sides_with_encountered_rejection_events` | Refined sides that encountered at least one earlier rejected route. |
| `final_outcome_side_counts` | JSON map of every terminal outcome to side count; sums to `audit_sides`. |
| `failure_counts` | JSON map of final failure reason to failed-side count only. |
| `failure_transcript_counts` | JSON map of final failure reason to deduplicated affected-transcript count. |
| `encountered_rejection_event_counts` | JSON map of reason to rejected evaluation-event count. |
| `encountered_rejection_side_counts` | JSON map of reason to sides encountering at least one such event. |
| `encountered_rejection_transcript_counts` | JSON map of reason to deduplicated transcripts encountering at least one such event. |

All maps are serialized as deterministic, key-sorted JSON. The former ambiguous
`rejected_*`, `structural_rejections_encountered`, and `unchanged` fields have
been replaced by these unit-specific fields. The former completion counters
were likewise renamed with `*_sides` or `*_transcripts` units.

Refined transcripts carry a concise `cds_terminus_refined` attribute. Refined
ordinary GFF3 children are rebuilt from the final geometry. CDS,
exons, phases, and any retained start/stop-codon children therefore match the
emitted CDS; source UTR fragments are retained only where they remain outside
the refined coding sequence. A source start/stop-codon child is suppressed when
the corresponding final calculated integrity status remains explicitly
`missing`; otherwise supported codon children are retained or rebuilt only
where they are consistent with the emitted final geometry. CDS-only output
preserves a non-zero first CDS phase. Other recognized translation-geometry
children that cannot be safely regenerated are omitted rather than copied at
stale coordinates. The ordinary `sequence_integrity_status` reports
the emitted model, while the selected model's original integrity remains
available as provenance and in the refinement audit. More generally, final
annotation writers use an explicit supported child set: exon, CDS, valid
start/stop codons, and intentionally retained or constructed UTR features.
Source intron, TSS/TES, and other arbitrary positional children are not
propagated into the final annotation. Standardized transcript bounds follow
the outer union of the exon geometry actually emitted, including recognized
CDS/UTR extensions, and standardized gene bounds follow the union of their
final transcript spans. Discarded source children cannot widen either parent.
