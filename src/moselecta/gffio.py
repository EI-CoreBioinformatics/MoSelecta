"""GFF/GTF parsing and identifier/attribute helpers."""

from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional
from urllib.parse import quote, unquote

from .constants import CHILD_FEATURE_ORDER, TRANSCRIPT_FEATURES
from .models import CodingModel, ExcludedModel, FeatureRecord, ParsedInput

def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def sanitize_id_part(text: str) -> str:
    text = str(text)
    text = re.sub(r"[^A-Za-z0-9_.:-]+", "_", text)
    text = text.strip("_")
    return text or "unknown"


def gff3_escape(value: object) -> str:
    return quote(str(value), safe="A-Za-z0-9_.:-")


def format_attrs(attrs: dict[str, object]) -> str:
    return ";".join(f"{gff3_escape(k)}={gff3_escape(v)}" for k, v in attrs.items())


def parse_attributes(attr_text: str) -> dict[str, str]:
    """Parse GFF3 key=value attributes, with a permissive GTF fallback."""
    attrs: dict[str, str] = {}
    attr_text = attr_text.strip()
    if not attr_text or attr_text == ".":
        return attrs

    # GFF3 style: ID=x;Parent=y
    if "=" in attr_text:
        for part in attr_text.split(";"):
            part = part.strip()
            if not part:
                continue
            if "=" in part:
                key, value = part.split("=", 1)
                attrs[unquote(key.strip())] = unquote(value.strip())
            else:
                attrs[unquote(part)] = ""
        return attrs

    # GTF style: gene_id "g1"; transcript_id "t1";
    for match in re.finditer(r'(\S+)\s+"([^"]*)"', attr_text):
        attrs[match.group(1)] = match.group(2)

    # Very permissive fallback for key value; entries.
    if not attrs:
        for part in attr_text.split(";"):
            part = part.strip()
            if not part:
                continue
            bits = part.split(None, 1)
            if len(bits) == 2:
                attrs[bits[0]] = bits[1].strip('"')
            else:
                attrs[bits[0]] = ""
    return attrs


def first_parent(record: FeatureRecord) -> Optional[str]:
    parents = record.parents
    return parents[0] if parents else None


def make_unique_source_names(inputs: list[Path], supplied: Optional[list[str]]) -> list[str]:
    if supplied:
        if len(supplied) != len(inputs):
            raise SystemExit("--source-name must be provided once per --input file")
        raw_names = supplied
    else:
        raw_names = [p.stem for p in inputs]

    seen: Counter[str] = Counter()
    unique: list[str] = []
    for i, name in enumerate(raw_names, start=1):
        base = sanitize_id_part(name)
        seen[base] += 1
        if seen[base] == 1:
            unique.append(base)
        else:
            unique.append(f"{base}_{seen[base]}")
    return unique


def read_gff_records(path: Path, source_name: str, strict: bool = False) -> list[FeatureRecord]:
    records: list[FeatureRecord] = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                message = f"{path}:{line_number}: expected 9 GFF/GFF3 columns, got {len(fields)}"
                if strict:
                    raise ValueError(message)
                warn(message + "; skipping line")
                continue
            seqid, source, feature_type, start_s, end_s, score, strand, phase, raw_attrs = fields
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError as exc:
                message = f"{path}:{line_number}: non-integer coordinates"
                if strict:
                    raise ValueError(message) from exc
                warn(message + "; skipping line")
                continue
            if start > end:
                message = f"{path}:{line_number}: start > end"
                if strict:
                    raise ValueError(message)
                warn(message + "; skipping line")
                continue
            records.append(
                FeatureRecord(
                    seqid=seqid,
                    source=source,
                    feature_type=feature_type,
                    start=start,
                    end=end,
                    score=score,
                    strand=strand,
                    phase=phase,
                    attrs=parse_attributes(raw_attrs),
                    raw_attrs=raw_attrs,
                    input_file=path,
                    source_name=source_name,
                    line_number=line_number,
                    record_index=len(records),
                    raw_line=line,
                )
            )
    return records


def infer_gene_id(
    tx_id: str,
    tx_gene: dict[str, str],
    synthetic_tx_to_gene: dict[str, str],
    gene_records: dict[str, FeatureRecord],
) -> str:
    if tx_id in tx_gene:
        return tx_gene[tx_id]
    if tx_id in synthetic_tx_to_gene:
        return synthetic_tx_to_gene[tx_id]
    # If a transcript_id is also a gene ID, use it as the gene ID.
    if tx_id in gene_records:
        return tx_id
    return f"gene_for_{tx_id}"


def parse_input(path: Path, source_name: str, strict: bool = False) -> ParsedInput:
    records = read_gff_records(path, source_name, strict=strict)

    gene_records: dict[str, FeatureRecord] = {}
    tx_records: dict[str, FeatureRecord] = {}
    tx_gene: dict[str, str] = {}

    for rec in records:
        if rec.feature_type == "gene":
            gene_id = rec.id or rec.attrs.get("gene_id")
            if gene_id:
                gene_records[gene_id] = rec
        elif rec.feature_type in TRANSCRIPT_FEATURES:
            tx_id = rec.id or rec.attrs.get("transcript_id")
            if tx_id:
                tx_records[tx_id] = rec
                gene_id = first_parent(rec) or rec.attrs.get("gene_id")
                if gene_id:
                    tx_gene[tx_id] = gene_id

    cds_by_tx: dict[str, list[FeatureRecord]] = defaultdict(list)
    synthetic_tx_to_gene: dict[str, str] = {}

    for rec in records:
        if rec.feature_type != "CDS":
            continue
        parents = rec.parents
        if not parents:
            # GTF fallback if attributes were parsed but Parent was absent.
            tx_id = rec.attrs.get("transcript_id")
            if tx_id:
                parents = [tx_id]
        if not parents:
            message = f"{path}:{rec.line_number}: CDS feature has no Parent/transcript_id"
            if strict:
                raise ValueError(message)
            warn(message + "; skipping CDS")
            continue
        for parent in parents:
            if parent in tx_records or parent in tx_gene:
                tx_id = parent
            elif parent in gene_records:
                tx_id = f"{parent}.t1"
                synthetic_tx_to_gene[tx_id] = parent
            else:
                tx_id = parent
            cds_by_tx[tx_id].append(rec)

    # Child features are normally parented by transcript ID. If a CDS was directly
    # parented by a gene and we generated a synthetic transcript, child records
    # parented by that gene are included under the synthetic transcript.
    children_by_tx: dict[str, list[FeatureRecord]] = defaultdict(list)
    tx_ids = set(cds_by_tx) | set(tx_records)
    for rec in records:
        if rec.feature_type in {"gene"} or rec.feature_type in TRANSCRIPT_FEATURES:
            continue
        for parent in rec.parents:
            if parent in tx_ids:
                children_by_tx[parent].append(rec)
            for synthetic_tx, gene_id in synthetic_tx_to_gene.items():
                if parent == gene_id:
                    children_by_tx[synthetic_tx].append(rec)

    coding_models: list[CodingModel] = []
    excluded_models: list[ExcludedModel] = []
    used_record_indices: set[int] = set()

    for tx_id, cds_records in sorted(cds_by_tx.items()):
        gene_id = infer_gene_id(tx_id, tx_gene, synthetic_tx_to_gene, gene_records)
        cds_chroms = sorted({r.seqid for r in cds_records})
        cds_strands = sorted({r.strand for r in cds_records})
        valid_strands = {s for s in cds_strands if s in {"+", "-"}}

        if len(cds_chroms) != 1:
            excluded_models.append(
                ExcludedModel(
                    source_name=source_name,
                    input_file=path,
                    input_gene_id=gene_id,
                    input_transcript_id=tx_id,
                    reason="mixed_cds_chromosomes",
                    cds_count=len(cds_records),
                    cds_chroms=cds_chroms,
                    cds_strands=cds_strands,
                )
            )
            continue
        if len(valid_strands) != 1 or len(cds_strands) != 1:
            excluded_models.append(
                ExcludedModel(
                    source_name=source_name,
                    input_file=path,
                    input_gene_id=gene_id,
                    input_transcript_id=tx_id,
                    reason="mixed_or_unknown_cds_strands",
                    cds_count=len(cds_records),
                    cds_chroms=cds_chroms,
                    cds_strands=cds_strands,
                )
            )
            continue

        chrom = cds_chroms[0]
        strand = cds_strands[0]
        gene_record = gene_records.get(gene_id)
        transcript_record = tx_records.get(tx_id)
        child_records = sorted(
            children_by_tx.get(tx_id, []),
            key=lambda r: (r.start, CHILD_FEATURE_ORDER.get(r.feature_type, 99), r.end, r.line_number),
        )
        if not child_records:
            # At least include CDS records if child grouping failed.
            child_records = sorted(cds_records, key=lambda r: (r.start, r.end, r.line_number))

        cds_parts = sorted((r.start0, r.end0) for r in cds_records)

        output_base = f"{source_name}__{sanitize_id_part(gene_id)}__{sanitize_id_part(tx_id)}"
        output_gene_id = f"{output_base}__gene"
        output_transcript_id = f"{output_base}__mRNA"

        model = CodingModel(
            source_name=source_name,
            input_file=path,
            input_gene_id=gene_id,
            input_transcript_id=tx_id,
            output_gene_id=output_gene_id,
            output_transcript_id=output_transcript_id,
            chrom=chrom,
            strand=strand,
            cds_parts=cds_parts,
            gene_record=gene_record,
            transcript_record=transcript_record,
            child_records=child_records,
        )

        if gene_record:
            model.used_record_indices.add(gene_record.record_index)
        if transcript_record:
            model.used_record_indices.add(transcript_record.record_index)
        for rec in child_records:
            model.used_record_indices.add(rec.record_index)
        used_record_indices.update(model.used_record_indices)
        coding_models.append(model)

    return ParsedInput(
        input_file=path,
        source_name=source_name,
        records=records,
        coding_models=sorted(coding_models, key=lambda m: (m.chrom, m.cds_start0, m.cds_end0, m.input_transcript_id)),
        excluded_models=excluded_models,
        used_record_indices=used_record_indices,
    )
