"""Shared console-report formatting for the command pipeline."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Callable, Collection, Sequence
from typing import TypeVar

from .scoring import value_as_float


T = TypeVar("T")


def log_section(title: str) -> None:
    """Print a top-level console report heading."""
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def log_subsection(title: str) -> None:
    """Print a console report subsection heading."""
    print("")
    print(title)
    print("-" * len(title))


def format_counter(counter: Counter[object]) -> str:
    """Format a deterministic key/count summary."""
    if not counter:
        return "."
    return ", ".join(
        f"{key}={value:,}" for key, value in sorted(counter.items())
    )


def log_distribution(
    title: str,
    counter: Counter[object],
    key_label: str,
    value_label: str,
) -> None:
    """Print a deterministic multiline distribution."""
    print(f"{title}:")
    if not counter:
        print("  .")
        return
    for key, value in sorted(counter.items()):
        print(f"  {key_label}={key}: {value_label}={value:,}")


def log_source_presence_breakdown(
    title: str,
    unit_label: str,
    items: Sequence[T],
    get_sources: Callable[[T], Collection[str]],
    all_sources: Sequence[str],
) -> None:
    """Print source presence grouped by the number of sources per item."""
    print(f"{title}:")
    source_count_values = sorted(Counter(len(get_sources(item)) for item in items))
    if not source_count_values:
        print("  .")
        return
    for source_count in source_count_values:
        subset = [item for item in items if len(get_sources(item)) == source_count]
        presence_counts = {
            source: sum(1 for item in subset if source in get_sources(item))
            for source in all_sources
        }
        presence_text = "; ".join(
            f"{source}={presence_counts[source]:,}" for source in all_sources
        )
        print(
            f"  source_count={source_count}: {unit_label}={len(subset):,}; "
            f"{presence_text}"
        )


def log_source_combination_breakdown(
    title: str,
    items: Sequence[T],
    get_sources: Callable[[T], Collection[str]],
) -> None:
    """Print exact source combinations grouped by source count."""
    print(f"{title}:")
    source_count_values = sorted(Counter(len(get_sources(item)) for item in items))
    if not source_count_values:
        print("  .")
        return
    for source_count in source_count_values:
        combo_counts = Counter(
            "+".join(sorted(get_sources(item)))
            for item in items
            if len(get_sources(item)) == source_count
        )
        combo_text = "; ".join(
            f"{combo}={count:,}" for combo, count in sorted(combo_counts.items())
        )
        print(f"  source_count={source_count}: {combo_text}")


def summarise_score_rows_by_source(
    rows: Sequence[dict[str, object]],
) -> dict[str, dict[str, float]]:
    """Calculate source-level counts and raw/scaled score sums."""
    stats: dict[str, dict[str, float]] = defaultdict(
        lambda: {
            "count": 0.0,
            "scaled_sum": 0.0,
            "raw_sum": 0.0,
        }
    )
    for row in rows:
        source = str(row.get("source_name", "."))
        stats[source]["count"] += 1.0
        stats[source]["scaled_sum"] += value_as_float(
            row.get("model_selection_score_scaled"), 0.0
        )
        stats[source]["raw_sum"] += value_as_float(
            row.get("model_selection_score"), 0.0
        )
    return stats


def summarise_source_scaled_scores(
    rows: Sequence[dict[str, object]],
) -> dict[str, dict[str, float]]:
    """Aggregate all-model and non-singleton scaled scores by source.

    ``sublocus_model_count`` is populated by the canonical sublocus builder on
    every scored candidate row.  A singleton sublocus contains exactly one
    candidate model, irrespective of exon count, corroboration, source labels,
    later selection, or pruning status.
    """
    stats: dict[str, dict[str, float]] = defaultdict(
        lambda: {
            "count": 0.0,
            "scaled_sum": 0.0,
            "raw_sum": 0.0,
            "non_singleton_count": 0.0,
            "non_singleton_scaled_sum": 0.0,
            "singleton_count": 0.0,
        }
    )
    for row in rows:
        source = str(row.get("source_name", "."))
        sublocus_model_count = value_as_float(
            row.get("sublocus_model_count"), -1.0
        )
        if sublocus_model_count < 1 or not sublocus_model_count.is_integer():
            raise ValueError(
                "Scored candidate row is missing a valid integer "
                "sublocus_model_count"
            )
        scaled_score = value_as_float(
            row.get("model_selection_score_scaled"), 0.0
        )
        stats[source]["count"] += 1.0
        stats[source]["scaled_sum"] += scaled_score
        stats[source]["raw_sum"] += value_as_float(
            row.get("model_selection_score"), 0.0
        )
        if sublocus_model_count == 1:
            stats[source]["singleton_count"] += 1.0
        else:
            stats[source]["non_singleton_count"] += 1.0
            stats[source]["non_singleton_scaled_sum"] += scaled_score
    return stats
