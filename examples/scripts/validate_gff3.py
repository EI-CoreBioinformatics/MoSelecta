#!/usr/bin/env python3
"""Perform independent structural checks on one or more GFF3 files."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path


def attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def validate(path: Path) -> tuple[Counter[str], list[str]]:
    ids: dict[str, tuple[str, int, int, str]] = {}
    parents: list[tuple[int, str, str, int, int, str]] = []
    counts: Counter[str] = Counter()
    errors: list[str] = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, start=1):
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                errors.append(f"line {line_number}: expected 9 columns")
                continue
            seqid, _, feature_type, start_text, end_text, _, strand, phase, attr_text = fields
            try:
                start, end = int(start_text), int(end_text)
            except ValueError:
                errors.append(f"line {line_number}: non-integer coordinate")
                continue
            if start < 1 or end < start:
                errors.append(f"line {line_number}: invalid coordinates {start}-{end}")
            if strand not in {"+", "-", ".", "?"}:
                errors.append(f"line {line_number}: invalid strand {strand}")
            if feature_type == "CDS" and phase not in {"0", "1", "2"}:
                errors.append(f"line {line_number}: invalid CDS phase {phase}")
            attrs = attributes(attr_text)
            feature_id = attrs.get("ID")
            if feature_id:
                if feature_id in ids:
                    errors.append(f"line {line_number}: duplicate ID {feature_id}")
                ids[feature_id] = (seqid, start, end, strand)
            for parent in attrs.get("Parent", "").split(","):
                if parent:
                    parents.append((line_number, parent, seqid, start, end, strand))
            counts[feature_type] += 1
    for line_number, parent, seqid, start, end, strand in parents:
        if parent not in ids:
            errors.append(f"line {line_number}: missing Parent {parent}")
            continue
        parent_seqid, parent_start, parent_end, parent_strand = ids[parent]
        if seqid != parent_seqid or start < parent_start or end > parent_end:
            errors.append(f"line {line_number}: child outside Parent {parent}")
        if strand in {"+", "-"} and parent_strand in {"+", "-"} and strand != parent_strand:
            errors.append(f"line {line_number}: strand differs from Parent {parent}")
    return counts, errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", type=Path)
    args = parser.parse_args()
    failed = False
    for path in args.paths:
        counts, errors = validate(path)
        print(f"{path}: records={sum(counts.values())}; types={dict(sorted(counts.items()))}; errors={len(errors)}")
        for error in errors[:20]:
            print(f"  {error}")
        failed = failed or bool(errors)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
