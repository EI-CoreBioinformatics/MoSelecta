"""Small TSV writing helpers shared across MoSelecta outputs."""

from __future__ import annotations

from pathlib import Path

from .scoring import fmt_optional_float


def format_tsv_value(value: object) -> str:
    if value is None:
        return "."
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return fmt_optional_float(value, 6)
    return str(value)


def write_dict_tsv(rows: list[dict[str, object]], path: Path, columns: list[str]) -> None:
    with path.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(format_tsv_value(row.get(col, ".")) for col in columns) + "\n")
