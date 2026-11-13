# Structural-coherence evidence and audit

> Technical reference for MoSelecta `0.2.16`.

The [Configuration reference](../configuration.md) identifies the aggregate
penalty control, [Scoring](scoring.md) defines its contribution, and
[Outputs](../outputs.md) identifies the candidate and region diagnostics. This
document owns the source-vote and bridge-group definitions.

The candidate structural-coherence term asks whether a candidate connects two
coding blocks that independent annotation sources instead represent as separate
source-gene hypotheses. Region discovery and connector diagnostics remain as in
the accepted structural-coherence implementation. The revised second pass
collects factual evidence before awarding independent source votes.

This term is part of the accepted model-only aggregate at default weight
`0.30`; see [Scoring](scoring.md). Exact model-anchor overlap values are cached within each superlocus for
performance; the cache does not change region discovery, voting, audit values,
or scoring semantics.

For each region, one source contributes at most one `C` vote if it supplies an
exact complete CDS chain, the exact simple intron or complete contiguous
multi-feature intron connector, >=50% coverage of both candidate anchors, a
credible spanning model in a `mixed` source, or continuity propagated across a
compatible bridge-hypothesis group. The candidate's own source is excluded from
independent `C` only after the factual routes are collected.

One source contributes `S` only if it has left and right local models from
resolvably distinct source genes and no continuity route. A shared gene ID on
either side suppresses pure separation; unresolved IDs also suppress `S`.
Identity resolution prefers `input_gene_id`, then source gene/transcript record
gene metadata, then a recognizable transcript-root suffix normalization.
Unequal transcript IDs alone are not evidence of separate genes.

The scoring formula is unchanged: `S >= 2`, `S > C`, at least two independent
sources covering each anchor, direction `(S-C)/(S+C)`, breadth
`min(1.5, 1+0.10*(S-2))`, maximum eligible region raw score, and the default
weight `0.30`. The >=50% anchor threshold is unchanged.

## Region TSV serialization

The existing `bilateral_separated_source_count` / `sources` columns now report
sources actually contributing `S`. The existing
`independent_spanning_source_count` / `sources` columns report sources actually
contributing `C`. `mixed_source_count` / `sources` and
`source_anchor_states` retain the original human-readable state labels;
`mixed` can therefore appear with `source_contributes_C=1` and
`source_contributes_S=0`.

Each new `source_*` column is a semicolon-separated `source:value` map in
alphabetical source order. Boolean values are `1` or `0`; no model/gene IDs
are serialized as `.`. The source-evidence model-ID field uses, for each source,
`left=...|right=...|crossing=...|chain=...|connector=...`. IDs within one field
are comma-separated. It identifies the actual models behind each flag.

The following columns make each vote reconstructable:

| Columns | Meaning |
|---|---|
| `source_has_exact_chain`, `source_has_exact_connector_continuity`, `source_has_anchor_crossing` | Direct factual continuity routes. Multi-feature exact connector requires one model with the complete contiguous intron chain. |
| `source_has_bridge_group_continuity`, `source_bridge_group_continuity_direct`, `source_bridge_group_continuity_propagated` | Continuity available through group membership. `direct=1` means direct local evidence or direct compatibility with an evidence-bearing member; `propagated=1` means the evidence is available only through the connected component. |
| `source_has_left_support`, `source_has_right_support`, `source_left_gene_ids`, `source_right_gene_ids` | Local support and resolved gene identities on both sides. |
| `source_gene_id_resolution_method`, `source_gene_relationship`, `source_same_gene_alternative` | Gene-ID method (`input_gene_id`, `source_record`, `normalized_transcript_id`, `unresolved`, or `.`) and why bilateral support can or cannot become `S`. |
| `source_state_label`, `source_contributes_S`, `source_contributes_C`, `source_S_reason`, `source_C_reason` | Human-readable anchor state, final independent votes, and primary reason. |
| `bridge_hypothesis_group_id`, `bridge_hypothesis_member_candidate_ids` | Group identity and member candidates; `.` means no multi-candidate group. |

`source_C_reason` uses the priority order `exact_complete_chain`,
`exact_connector`, `anchor_based_non_exact_continuity`,
`mixed_source_spanning_interpretation`, `bridge_hypothesis_grouping`. All
Boolean route flags remain visible even when a stronger primary reason wins.
`source_S_reason` is `distinct_source_genes` only for an awarded `S` vote;
the gene-relationship column explains suppressed votes. A source can never
contribute both `S` and `C` to one region.

Bridge groups are formed only within one sublocus and strand. Regions from
different candidate models must have overlapping connector intervals that are
nested or share a boundary, >=50% overlap of both corresponding shorter
anchors, and preliminary separation sets sharing at least two sources and
>=2/3 of the smaller set. Connected-component propagation remains auditable
through the direct/propagated distinction. Large or non-pairwise-equivalent
components should be inspected, not presumed biologically identical merely
because they share a sublocus.

The model-selection score table retains its existing candidate-level coherence
columns and is not widened with these per-source fields. An identical CDS chain
within one sublocus must have identical coherence raw scores and penalties;
source priors remain separate and may distinguish the final model identity.
