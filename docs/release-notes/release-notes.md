# Release notes

## MoSelecta v0.3.0 release notes

MoSelecta integrates CDS-containing gene-prediction annotations into a
deterministic, auditable annotation. This first public release provides
structural and evidence-aware scoring, conflict resolution and backfill, fixed
or automatic pruning, stable primaries, optional alternatives, bounded CDS
terminus refinement, and separate optional UTR augmentation.

Inputs are prediction GFF/GTF annotations. Genome FASTA, protein alignment
annotations, and transcript exon annotations are optional. The normal product
is `<prefix>.final_annotation.gff3`; the standard profile also writes the main
decision ledgers. The UTR-augmented annotation is separate.

The integrated unscaled `model_selection_score` drives ranking. Its scaled form
is a clipped positive transform used for reporting and configured pruning
thresholds. Source-level means are descriptive run diagnostics, not validated
predictor-accuracy measures. Corroboration reflects distinct prediction sources
within the current run and does not make correlated sources independent.

MoSelecta supports Python 3.10–3.15 and is released under the MIT licence.

