"""Tolerance-aware output comparison for MoSelecta regression checks."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

_NUMBER_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
_NUMBER_TOKEN_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")


@dataclass
class FileComparison:
    left: Path
    right: Path
    status: str
    message: str = ""
    differing_cells: int = 0
    max_abs_delta: float = 0.0
    max_rel_delta: float = 0.0

    @property
    def equivalent(self) -> bool:
        return self.status in {"identical", "equivalent"}


def _is_number(value: str) -> bool:
    if value in {"", ".", "NA", "NaN", "nan"}:
        return False
    return bool(_NUMBER_RE.match(value))


def _as_float(value: str) -> Optional[float]:
    if not _is_number(value):
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    if math.isnan(number):
        return None
    return number


def _numbers_close(left: float, right: float, *, abs_tol: float, rel_tol: float) -> tuple[bool, float, float]:
    abs_delta = abs(left - right)
    denom = max(abs(left), abs(right), 1.0)
    rel_delta = abs_delta / denom
    return abs_delta <= (max(abs_tol, rel_tol * denom) + 1e-15), abs_delta, rel_delta


def _compare_numeric_or_component_values(
    left_value: str,
    right_value: str,
    *,
    abs_tol: float,
    rel_tol: float,
) -> tuple[bool, float, float]:
    """Return whether two cells differ only by tolerated numeric values.

    In addition to plain numeric TSV cells, MoSelecta reports some score
    component columns as deterministic strings such as
    ``Annevo:0.701416;Helixer:0.480938``.  Refactors can produce harmless
    last-decimal float differences inside these component strings.  This
    helper treats them as equivalent only when the surrounding text/template is
    identical and all paired numeric tokens are within tolerance.
    """
    left_float = _as_float(left_value)
    right_float = _as_float(right_value)
    if left_float is not None and right_float is not None:
        return _numbers_close(left_float, right_float, abs_tol=abs_tol, rel_tol=rel_tol)

    left_matches = list(_NUMBER_TOKEN_RE.finditer(left_value))
    right_matches = list(_NUMBER_TOKEN_RE.finditer(right_value))
    if not left_matches or len(left_matches) != len(right_matches):
        return False, 0.0, 0.0

    left_template = _NUMBER_TOKEN_RE.sub("{#}", left_value)
    right_template = _NUMBER_TOKEN_RE.sub("{#}", right_value)
    if left_template != right_template:
        return False, 0.0, 0.0

    max_abs_delta = 0.0
    max_rel_delta = 0.0
    for left_match, right_match in zip(left_matches, right_matches):
        try:
            left_number = float(left_match.group(0))
            right_number = float(right_match.group(0))
        except ValueError:
            return False, 0.0, 0.0
        if math.isnan(left_number) or math.isnan(right_number):
            return False, 0.0, 0.0
        close, abs_delta, rel_delta = _numbers_close(left_number, right_number, abs_tol=abs_tol, rel_tol=rel_tol)
        max_abs_delta = max(max_abs_delta, abs_delta)
        max_rel_delta = max(max_rel_delta, rel_delta)
        if not close:
            return False, max_abs_delta, max_rel_delta
    return True, max_abs_delta, max_rel_delta


def _read_text(path: Path) -> str:
    return path.read_text(errors="replace")


def compare_tsv_files(left: Path, right: Path, *, abs_tol: float, rel_tol: float) -> FileComparison:
    with left.open(newline="") as left_handle, right.open(newline="") as right_handle:
        left_rows = list(csv.reader(left_handle, delimiter="\t"))
        right_rows = list(csv.reader(right_handle, delimiter="\t"))

    if len(left_rows) != len(right_rows):
        return FileComparison(left, right, "different", f"row count differs: {len(left_rows)} vs {len(right_rows)}")

    differing_cells = 0
    tolerated_numeric_cells = 0
    max_abs_delta = 0.0
    max_rel_delta = 0.0
    first_difference = ""

    for row_index, (left_row, right_row) in enumerate(zip(left_rows, right_rows), start=1):
        if len(left_row) != len(right_row):
            return FileComparison(left, right, "different", f"column count differs at row {row_index}: {len(left_row)} vs {len(right_row)}")
        for col_index, (left_value, right_value) in enumerate(zip(left_row, right_row), start=1):
            if left_value == right_value:
                continue
            close, abs_delta, rel_delta = _compare_numeric_or_component_values(
                left_value,
                right_value,
                abs_tol=abs_tol,
                rel_tol=rel_tol,
            )
            max_abs_delta = max(max_abs_delta, abs_delta)
            max_rel_delta = max(max_rel_delta, rel_delta)
            if close:
                tolerated_numeric_cells += 1
                continue
            differing_cells += 1
            if not first_difference:
                column = left_rows[0][col_index - 1] if left_rows and len(left_rows[0]) >= col_index else f"column_{col_index}"
                first_difference = (
                    f"first differing cell row={row_index} col={col_index} ({column!r}): "
                    f"{left_value!r} != {right_value!r}"
                )

    if differing_cells:
        return FileComparison(left, right, "different", first_difference, differing_cells, max_abs_delta, max_rel_delta)
    if tolerated_numeric_cells:
        return FileComparison(
            left,
            right,
            "equivalent",
            f"only tolerated numeric/component differences: cells={tolerated_numeric_cells}; max_abs_delta={max_abs_delta:.6g}; max_rel_delta={max_rel_delta:.6g}",
            0,
            max_abs_delta,
            max_rel_delta,
        )
    return FileComparison(left, right, "identical")


def compare_text_files(left: Path, right: Path, *, abs_tol: float, rel_tol: float) -> FileComparison:
    left_text = _read_text(left)
    right_text = _read_text(right)
    if left_text == right_text:
        return FileComparison(left, right, "identical")
    if left.suffix.lower() in {".tsv", ".tab"} or right.suffix.lower() in {".tsv", ".tab"}:
        return compare_tsv_files(left, right, abs_tol=abs_tol, rel_tol=rel_tol)
    return FileComparison(left, right, "different", "text differs")


def compare_file_pair(left: Path, right: Path, *, abs_tol: float = 1e-6, rel_tol: float = 1e-9) -> FileComparison:
    if not left.exists():
        return FileComparison(left, right, "missing-left", "left file is missing")
    if not right.exists():
        return FileComparison(left, right, "missing-right", "right file is missing")
    if left.is_dir() or right.is_dir():
        return FileComparison(left, right, "different", "expected two files, got at least one directory")
    return compare_text_files(left, right, abs_tol=abs_tol, rel_tol=rel_tol)


def iter_directory_pairs(left_dir: Path, right_dir: Path, pattern: str) -> Iterable[tuple[Path, Path]]:
    for left_file in sorted(p for p in left_dir.glob(pattern) if p.is_file()):
        yield left_file, right_dir / left_file.name


def iter_prefix_pairs(left_prefix: Path, right_prefix: Path) -> Iterable[tuple[Path, Path]]:
    left_parent = left_prefix.parent if str(left_prefix.parent) else Path(".")
    prefix_name = left_prefix.name
    for left_file in sorted(p for p in left_parent.glob(prefix_name + "*") if p.is_file()):
        suffix = left_file.name[len(prefix_name):]
        yield left_file, right_prefix.parent / (right_prefix.name + suffix)


def compare_pairs(pairs: Iterable[tuple[Path, Path]], *, abs_tol: float, rel_tol: float) -> list[FileComparison]:
    return [compare_file_pair(left, right, abs_tol=abs_tol, rel_tol=rel_tol) for left, right in pairs]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare MoSelecta output files with tolerance for last-decimal floating-point differences.",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--files", nargs=2, metavar=("LEFT", "RIGHT"), help="Compare two files.")
    mode.add_argument("--dirs", nargs=2, metavar=("LEFT_DIR", "RIGHT_DIR"), help="Compare files with the same names in two directories.")
    mode.add_argument("--prefixes", nargs=2, metavar=("LEFT_PREFIX", "RIGHT_PREFIX"), help="Compare all files matching LEFT_PREFIX* to RIGHT_PREFIX* with matching suffixes.")
    parser.add_argument("--glob", default="*", help="Glob used with --dirs. Default: *")
    parser.add_argument("--abs-tol", type=float, default=1e-6, help="Absolute tolerance for numeric TSV cells. Default: 1e-6")
    parser.add_argument("--rel-tol", type=float, default=1e-9, help="Relative tolerance for numeric TSV cells. Default: 1e-9")
    parser.add_argument("--quiet", action="store_true", help="Only print failures and the final summary.")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.files:
        pairs = [(Path(args.files[0]), Path(args.files[1]))]
    elif args.dirs:
        pairs = list(iter_directory_pairs(Path(args.dirs[0]), Path(args.dirs[1]), args.glob))
    else:
        pairs = list(iter_prefix_pairs(Path(args.prefixes[0]), Path(args.prefixes[1])))

    if not pairs:
        print("No files matched for comparison.", file=sys.stderr)
        return 2

    results = compare_pairs(pairs, abs_tol=args.abs_tol, rel_tol=args.rel_tol)
    counts: dict[str, int] = {}
    for result in results:
        counts[result.status] = counts.get(result.status, 0) + 1
        if not args.quiet or not result.equivalent:
            relative_label = result.left.name
            if result.status == "identical":
                print(f"IDENTICAL\t{relative_label}")
            elif result.status == "equivalent":
                print(f"EQUIVALENT\t{relative_label}\t{result.message}")
            else:
                print(f"DIFFERENT\t{relative_label}\t{result.message}")

    total = len(results)
    summary = ", ".join(f"{status}={count}" for status, count in sorted(counts.items()))
    print(f"Compared {total} file(s): {summary}")
    return 0 if all(result.equivalent for result in results) else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
