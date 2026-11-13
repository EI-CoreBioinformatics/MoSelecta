"""Shared normalization for direct and ``@manifest`` input paths."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


@dataclass(frozen=True, slots=True)
class ResolvedInput:
    path: Path
    manifest_label: str | None
    origin: str


@dataclass(frozen=True, slots=True)
class ResolvedInputFamily:
    entries: tuple[ResolvedInput, ...]
    labels: tuple[str, ...] | None
    expansion_messages: tuple[str, ...]

    @property
    def paths(self) -> list[Path]:
        return [entry.path for entry in self.entries]


def split_input_tokens(values: Iterable[object] | None, *, comma_separated: bool) -> list[str]:
    """Flatten the CLI forms already accepted by one input family."""
    tokens: list[str] = []
    for value in values or ():
        parts = str(value).split(",") if comma_separated else (str(value),)
        tokens.extend(part.strip() for part in parts if part.strip())
    return tokens


def split_label_tokens(values: Iterable[object] | None, *, comma_separated: bool) -> list[str]:
    labels: list[str] = []
    for value in values or ():
        parts = str(value).split(",") if comma_separated else (str(value),)
        labels.extend(part.strip() for part in parts if part.strip())
    return labels


def _manifest_entries(manifest: Path) -> list[ResolvedInput]:
    if not manifest.is_file():
        raise SystemExit(f"Input manifest does not exist or is not a file: {manifest}")
    entries: list[ResolvedInput] = []
    for line_number, raw_line in enumerate(manifest.read_text().splitlines(), start=1):
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        fields = raw_line.split("\t")
        if len(fields) > 2:
            raise SystemExit(
                f"{manifest}:{line_number}: expected one or two tab-separated columns, "
                f"got {len(fields)}"
            )
        path_text = fields[0].strip()
        if not path_text:
            raise SystemExit(f"{manifest}:{line_number}: empty input path")
        if path_text.startswith("@"):
            raise SystemExit(
                f"{manifest}:{line_number}: nested input manifests are not supported: "
                f"{path_text}"
            )
        label = fields[1].strip() if len(fields) == 2 else ""
        path = Path(path_text)
        if not path.is_absolute():
            path = manifest.parent / path
        path = path.resolve(strict=False)
        if not path.exists():
            raise SystemExit(
                f"{manifest}:{line_number}: referenced input does not exist: {path}"
            )
        entries.append(ResolvedInput(
            path=path,
            manifest_label=label or None,
            origin=f"manifest:{manifest}:{line_number}",
        ))
    return entries


def resolve_input_family(
    raw_values: Iterable[object] | None,
    explicit_labels: Iterable[object] | None,
    *,
    family_name: str,
    input_option: str,
    label_option: str,
    comma_separated_inputs: bool,
    comma_separated_labels: bool,
    fallback_labels: Sequence[str] | None = None,
    fallback_prefix: str | None = None,
) -> ResolvedInputFamily:
    """Expand manifests, reject duplicates and resolve label precedence.

    ``fallback_labels`` must correspond to the final expanded entries. When it
    is omitted, path stems are used (the established model-input fallback).
    """
    entries: list[ResolvedInput] = []
    messages: list[str] = []
    for token in split_input_tokens(raw_values, comma_separated=comma_separated_inputs):
        if token.startswith("@"):
            manifest_text = token[1:].strip()
            if not manifest_text:
                raise SystemExit(f"{input_option} contains an empty @manifest token")
            manifest = Path(manifest_text).resolve(strict=False)
            expanded = _manifest_entries(manifest)
            entries.extend(expanded)
            messages.append(
                f"Expanded @{manifest}: {len(expanded):,} {family_name} inputs"
            )
        else:
            entries.append(ResolvedInput(Path(token), None, "direct"))

    first_by_path: dict[Path, ResolvedInput] = {}
    for entry in entries:
        normalized = entry.path.resolve(strict=False)
        previous = first_by_path.get(normalized)
        if previous is not None:
            raise SystemExit(
                f"Duplicate {family_name} input after manifest expansion: {normalized} "
                f"(first: {previous.origin}; again: {entry.origin})"
            )
        first_by_path[normalized] = entry

    supplied = split_label_tokens(
        explicit_labels, comma_separated=comma_separated_labels
    )
    if supplied:
        if len(supplied) != len(entries):
            raise SystemExit(
                f"{label_option} supplied {len(supplied)} labels for "
                f"{len(entries)} expanded {family_name} inputs"
            )
        effective_labels: tuple[str, ...] | None = tuple(supplied)
    elif any(entry.manifest_label is not None for entry in entries):
        if fallback_labels is None:
            fallbacks = (
                [f"{fallback_prefix}_{index:03d}" for index in range(1, len(entries) + 1)]
                if fallback_prefix is not None
                else [entry.path.stem for entry in entries]
            )
        else:
            if len(fallback_labels) != len(entries):
                raise ValueError("fallback label count must match expanded input count")
            fallbacks = list(fallback_labels)
        effective_labels = tuple(
            entry.manifest_label or fallback
            for entry, fallback in zip(entries, fallbacks)
        )
    else:
        effective_labels = None

    return ResolvedInputFamily(tuple(entries), effective_labels, tuple(messages))
