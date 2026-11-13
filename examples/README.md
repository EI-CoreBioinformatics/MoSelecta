# MoSelecta Bombus 4 Mb functional example

This archive demonstrates MoSelecta on the first 4,000,000 bp
of the Bombus reference sequence `NC_063270.1` (1-based interval
`NC_063270.1:1-4000000`). Coordinates and the sequence identifier are unchanged.

The example combines seven gene-prediction sources with protein and transcript
evidence, runs strict input parsing and UTR augmentation, and supplies a compact
set of validated reference outputs. It is a functional software demonstration,
not a genome-wide biological performance benchmark.

## Archive contents

- `inputs/genome/`: 4 Mb FASTA and its `.fai` index.
- `inputs/predictions/`: seven prediction GFF3 files.
- `inputs/protein_evidence/`: seven Miniprot GFF files.
- `inputs/transcript_evidence/`: 33 StringTie GTF files.
- `inputs/*_inputs.txt` and `inputs/*_evidence.txt`: the three ordered MoSelecta
  manifests, containing 47 prediction/evidence entries. Together with the FASTA
  and index, the archive has 49 validated data input files.
- `inputs/PROVENANCE.tsv`: original/extracted filenames, SHA-256 values, sizes
  and extraction counts retained from release preparation.
- `outputs/`: 11 representative products from the verified run.
- `SHA256SUMS`: checksums for every distributed file except the checksum file
  itself.

The seven predictors are Annevo, BRAKER, GALBA, Helixer, OrionGeno, Tiberius and
Vipsania. Their extracted transcript counts are 83, 132, 189, 102, 92, 88 and
89 respectively (775 total).

## Provenance and known gaps

The source records available for release preparation establish the following.
They do not contain enough information to infer missing versions, download URLs
or licences, so those fields remain explicit publication checks.

| Material | Established identifier/generator | Original record | Interval and selection | Provider, version, licence or download record |
|---|---|---|---|---|
| Reference genome | `NC_063270.1` | supplied `genome.fa` | bases 1–4,000,000 copied without offset | Upstream provider/download and redistribution terms not recorded in the supplied provenance |
| Predictions | Annevo, BRAKER, GALBA, Helixer, OrionGeno, Tiberius, Vipsania | `bombus_terrestris_<predictor>.gff3` | complete CDS-containing transcript models wholly inside the interval; boundary-crossing models were not clipped | Predictor versions/settings, original dataset location and redistribution terms not recorded |
| Protein evidence | Miniprot; `GCF_000188095.3`, `GCF_011952205.1`, `GCF_011952255.1`, `GCF_011952275.1`, `GCF_014825855.1`, `GCF_024516045.1`, `GCF_024542735.1` | `<assembly>.miniprot.gff` | complete alignment chains inside the interval; no artificial clipping | Miniprot version/settings, protein source/download and redistribution terms not recorded |
| Transcript evidence | StringTie; 33 run-labelled files | `<run>.stringtie.gtf` | complete exon-bearing transcript chains inside the interval; no artificial clipping | StringTie version/settings, original study/download records and redistribution terms not recorded |

The 33 transcript run identifiers are listed in
`inputs/transcript_evidence.txt`. Exact source-file and extracted-file hashes are
in `inputs/PROVENANCE.tsv`.

The project owner intends to publish this example alongside MoSelecta. That
intention is not evidence of permission from external data providers. Before
public release, confirm the upstream genome, predictor-output, protein-evidence
and transcript-evidence redistribution terms and add the confirmed source URLs,
versions and licences to this table.

## Verify the archive

The ZIP extracts into this single directory:

```bash
unzip MoSelecta-bombus-4mb-example.zip
cd bombus_4mb
```

On Linux:

```bash
sha256sum --check SHA256SUMS
```

On macOS:

```bash
shasum -a 256 --check SHA256SUMS
```

## Run MoSelecta

```bash
mkdir -p run
moselecta \
  --inputs @inputs/model_inputs.txt \
  --genome-fasta inputs/genome/NC_063270.1_1-4000000.fa \
  --evidence-proteins @inputs/protein_evidence.txt \
  --evidence-transcripts @inputs/transcript_evidence.txt \
  --strict \
  --add-utrs \
  --output-profile standard \
  --progress-log-level basic \
  --out-prefix run/bombus_4mb
```

The parent directory of `--out-prefix` must exist. The manifests preserve the
validated file order and source labels; do not change their membership when
checking reproducibility.

## Expected result

The verified run reports:

- 775 input prediction models;
- 144 superloci and 171 subloci;
- 95 final genes;
- 135 final transcripts: 95 primary and 40 alternative;
- 1,120 CDS features in each final annotation;
- 173 five-prime UTR and 141 three-prime UTR features in the UTR-augmented
  annotation; and
- zero structural GFF3 errors in the standardised and final products.

The principal reference products are:

- `outputs/bombus_4mb.final_annotation.gff3` — ordinary final annotation;
- `outputs/bombus_4mb.final_annotation.with_utrs.gff3` — separate
  UTR-augmented annotation;
- `outputs/bombus_4mb.standardised_primary_models.gff3` — retained primaries;
- input, protein, transcript, scoring, pruning, source-weight and UTR summaries;
  and
- `outputs/bombus_4mb.run_configuration.json` — portable record of the verified
  configuration.

To compare the two core regenerated annotations with the supplied references:

```bash
cmp run/bombus_4mb.final_annotation.gff3 \
  outputs/bombus_4mb.final_annotation.gff3
cmp run/bombus_4mb.final_annotation.with_utrs.gff3 \
  outputs/bombus_4mb.final_annotation.with_utrs.gff3
```

The four supplied metadata/report files that record input filenames use portable
paths relative to `bombus_4mb/`. Runtime and output-prefix fields in a new run
will naturally reflect that run's location and timing; these do not affect the
biological annotations.

## Scope and limitations

Models or evidence chains crossing either interval boundary were excluded rather
than truncated. Counts, evidence density, runtime and memory are specific to this
4 Mb region and the supplied source collection. They must not be interpreted as
whole-genome accuracy, sensitivity or predictor-performance estimates.
