"""Indexed genomic FASTA access and model sequence-integrity metrics.

The implementation is deliberately superlocus-oriented: one bounded genomic
window is fetched for each candidate superlocus, all models in that window are
assessed, and the sequence is then released.  This keeps memory bounded for
large genomes and leaves a natural unit for future parallel execution.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import os
import statistics
import tempfile
from typing import BinaryIO, Iterable, Optional

from .cds_phase import (
    CANONICALISED_FIVE_PRIME_BOUNDARY,
    PRESERVED_NUMERIC_SOURCE,
    RECOVERED_MISSING_SOURCE_PHASE,
    STANDARDISED_SOURCE_MISSING_TO_ZERO,
    STANDARD_START_CODONS,
    STANDARD_STOP_CODONS,
    UNRESOLVED_SOURCE_MISSING,
    analyse_cds_frame,
    biological_cds_parts,
    source_first_cds_phase,
)
from .constants import PROGRAM_NAME, PROGRAM_VERSION
from .models import CDSPhaseOutput, CodingModel, SequenceIntegrityMetrics, Superlocus
from .tsvio import format_tsv_value, write_dict_tsv


CANONICAL_SPLICE_MOTIFS = frozenset({"GT-AG", "GC-AG", "AT-AC"})
@dataclass(frozen=True, slots=True)
class FastaIndexEntry:
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFasta:
    """Small dependency-free reader for ordinary FASTA plus ``.fai`` index."""

    def __init__(self, path: Path):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"Genome FASTA does not exist: {self.path}")
        if self.path.suffix.lower() in {".gz", ".bgz", ".bgzf"}:
            raise ValueError(
                f"{PROGRAM_NAME} {PROGRAM_VERSION} sequence-integrity analysis currently requires an uncompressed FASTA: "
                f"{self.path}"
            )
        self.index_path = Path(str(self.path) + ".fai")
        self.index_action = "reused"
        self.index_rebuild_reason = "."
        if not self.index_path.exists():
            build_fasta_index(self.path, self.index_path)
            self.index_action = "created"
        else:
            rebuild_reason = fasta_index_rebuild_reason(self.path, self.index_path)
            if rebuild_reason is not None:
                build_fasta_index(self.path, self.index_path)
                self.index_action = "rebuilt"
                self.index_rebuild_reason = rebuild_reason
        # Backwards-compatible indicator used by existing summaries/tests. It is
        # true whenever this process created or replaced the index.
        self.index_built = self.index_action in {"created", "rebuilt"}
        self.index_rebuilt = self.index_action == "rebuilt"
        self.entries = read_fasta_index(self.index_path)
        self._handle: Optional[BinaryIO] = None

    def __enter__(self) -> "IndexedFasta":
        self._handle = self.path.open("rb")
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def contig_length(self, name: str) -> Optional[int]:
        entry = self.entries.get(name)
        return entry.length if entry else None

    def fetch(self, name: str, start0: int, end0: int) -> str:
        """Fetch a 0-based half-open interval without loading a full contig."""
        if self._handle is None:
            raise RuntimeError("IndexedFasta must be used as a context manager")
        entry = self.entries.get(name)
        if entry is None:
            raise KeyError(name)
        if start0 < 0 or end0 < start0 or end0 > entry.length:
            raise ValueError(
                f"FASTA interval out of range for {name}: {start0}-{end0}; contig_length={entry.length}"
            )
        if start0 == end0:
            return ""
        if entry.line_bases <= 0 or entry.line_width <= 0:
            raise ValueError(f"Cannot fetch from empty or invalid FASTA contig: {name}")

        chunks: list[bytes] = []
        position = start0
        while position < end0:
            line_index, column = divmod(position, entry.line_bases)
            take = min(end0 - position, entry.line_bases - column)
            byte_offset = entry.offset + line_index * entry.line_width + column
            self._handle.seek(byte_offset)
            chunk = self._handle.read(take)
            if len(chunk) != take:
                raise ValueError(
                    f"Unexpected end of FASTA while reading {name}:{position + 1}-{position + take}"
                )
            chunks.append(chunk)
            position += take
        return b"".join(chunks).decode("ascii").upper()


def build_fasta_index(fasta_path: Path, index_path: Optional[Path] = None) -> Path:
    """Build a samtools-compatible ``.fai`` for an uncompressed FASTA.

    Sequence lines must have a consistent width within each contig, except for
    the final line.  Both LF and CRLF files are supported.
    """
    fasta_path = Path(fasta_path)
    index_path = Path(index_path) if index_path is not None else Path(str(fasta_path) + ".fai")
    if fasta_path.suffix.lower() in {".gz", ".bgz", ".bgzf"}:
        raise ValueError(
            "Automatic FASTA indexing currently supports ordinary uncompressed FASTA only. "
            f"Provide an uncompressed FASTA: {fasta_path}"
        )

    entries: list[FastaIndexEntry] = []
    with fasta_path.open("rb") as handle:
        name: Optional[str] = None
        length = 0
        seq_offset = 0
        line_bases = 0
        line_width = 0
        previous_line_was_short = False

        def finish_contig() -> None:
            nonlocal name, length, seq_offset, line_bases, line_width, previous_line_was_short
            if name is None:
                return
            entries.append(
                FastaIndexEntry(
                    name=name,
                    length=length,
                    offset=seq_offset,
                    line_bases=line_bases,
                    line_width=line_width,
                )
            )
            name = None
            length = 0
            seq_offset = 0
            line_bases = 0
            line_width = 0
            previous_line_was_short = False

        while True:
            line_start = handle.tell()
            raw = handle.readline()
            if not raw:
                finish_contig()
                break
            if raw.startswith(b">"):
                finish_contig()
                header = raw[1:].strip().split(None, 1)
                if not header:
                    raise ValueError(f"Empty FASTA header at byte offset {line_start}: {fasta_path}")
                name = header[0].decode("utf-8")
                seq_offset = handle.tell()
                continue
            if name is None:
                if raw.strip():
                    raise ValueError(f"FASTA sequence encountered before first header: {fasta_path}")
                continue
            seq = raw.rstrip(b"\r\n")
            if not seq:
                raise ValueError(f"Blank FASTA sequence line within contig {name}")
            if line_bases == 0:
                line_bases = len(seq)
                line_width = len(raw)
            else:
                if previous_line_was_short:
                    raise ValueError(
                        f"Inconsistent FASTA wrapping for {name}: a short sequence line was followed by more sequence"
                    )
                if len(seq) > line_bases or (len(seq) == line_bases and len(raw) != line_width):
                    raise ValueError(f"Inconsistent FASTA line width for contig {name}")
                previous_line_was_short = len(seq) < line_bases
            length += len(seq)

    seen: set[str] = set()
    for entry in entries:
        if entry.name in seen:
            raise ValueError(f"Duplicate FASTA contig name: {entry.name}")
        seen.add(entry.name)

    index_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=index_path.parent,
            prefix=f".{index_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as out:
            temp_path = Path(out.name)
            for entry in entries:
                out.write(
                    f"{entry.name}\t{entry.length}\t{entry.offset}\t{entry.line_bases}\t{entry.line_width}\n"
                )
            out.flush()
            os.fsync(out.fileno())
        os.replace(temp_path, index_path)
        # Ensure a rebuilt index is not immediately considered stale when the
        # FASTA carries a future or preserved modification timestamp.
        fasta_mtime_ns = fasta_path.stat().st_mtime_ns
        index_stat = index_path.stat()
        if index_stat.st_mtime_ns < fasta_mtime_ns:
            os.utime(index_path, ns=(index_stat.st_atime_ns, fasta_mtime_ns))
    finally:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink()
    return index_path


def read_fasta_index(index_path: Path) -> dict[str, FastaIndexEntry]:
    entries: dict[str, FastaIndexEntry] = {}
    with Path(index_path).open() as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Malformed FASTA index {index_path}:{line_number}")
            name = fields[0]
            try:
                entry = FastaIndexEntry(
                    name=name,
                    length=int(fields[1]),
                    offset=int(fields[2]),
                    line_bases=int(fields[3]),
                    line_width=int(fields[4]),
                )
            except ValueError as exc:
                raise ValueError(f"Malformed numeric FASTA index field {index_path}:{line_number}") from exc
            if name in entries:
                raise ValueError(f"Duplicate FASTA index contig name: {name}")
            entries[name] = entry
    if not entries:
        raise ValueError(f"FASTA index contains no contigs: {index_path}")
    return entries


def validate_fasta_index_entries(
    entries: dict[str, FastaIndexEntry],
    *,
    fasta_size: int,
) -> Optional[str]:
    """Return a rebuild reason when index geometry cannot fit the FASTA file."""
    if fasta_size <= 0:
        return "genome_fasta_is_empty"
    for entry in entries.values():
        if entry.length < 0 or entry.offset < 0:
            return f"invalid_negative_index_value:{entry.name}"
        if entry.length == 0:
            # Empty contigs have no fetchable sequence, but their offset must
            # still lie within or immediately after the FASTA file.
            if entry.offset > fasta_size:
                return f"index_offset_out_of_range:{entry.name}"
            continue
        if entry.line_bases <= 0 or entry.line_width < entry.line_bases:
            return f"invalid_line_geometry:{entry.name}"
        if entry.offset >= fasta_size:
            return f"index_offset_out_of_range:{entry.name}"
        final_base_offset = (
            entry.offset
            + ((entry.length - 1) // entry.line_bases) * entry.line_width
            + ((entry.length - 1) % entry.line_bases)
        )
        if final_base_offset >= fasta_size:
            return f"indexed_sequence_exceeds_fasta_size:{entry.name}"
    return None


def fasta_index_rebuild_reason(fasta_path: Path, index_path: Path) -> Optional[str]:
    """Return why an existing index should be rebuilt, or ``None`` if reusable.

    Standard ``.fai`` files do not store a FASTA checksum. The practical stale
    index checks are therefore modification time plus validation that every
    indexed offset and line geometry can fit inside the current FASTA file.
    """
    fasta_path = Path(fasta_path)
    index_path = Path(index_path)
    try:
        fasta_stat = fasta_path.stat()
        index_stat = index_path.stat()
    except OSError as exc:
        return f"index_stat_error:{type(exc).__name__}"
    if index_stat.st_mtime_ns < fasta_stat.st_mtime_ns:
        return "genome_fasta_newer_than_index"
    try:
        entries = read_fasta_index(index_path)
    except (OSError, ValueError) as exc:
        return f"invalid_index:{type(exc).__name__}"
    return validate_fasta_index_entries(entries, fasta_size=fasta_stat.st_size)


def reverse_complement(sequence: str) -> str:
    table = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn",
    )
    return sequence.translate(table)[::-1].upper()


def parse_cds_phase(value: object) -> Optional[int]:
    try:
        phase = int(str(value))
    except (TypeError, ValueError):
        return None
    return phase if phase in {0, 1, 2} else None


def model_cds_phases(model: CodingModel) -> dict[tuple[int, int], Optional[int]]:
    phases: dict[tuple[int, int], Optional[int]] = {}
    for record in model.child_records:
        if record.feature_type != "CDS":
            continue
        key = (record.start0, record.end0)
        parsed = parse_cds_phase(record.phase)
        if key not in phases or (phases[key] is None and parsed is not None):
            phases[key] = parsed
    return phases


# Private aliases retained for existing internal callers and downstream users.
_parse_phase = parse_cds_phase
_model_cds_phases = model_cds_phases


def _slice_window(window: str, window_start0: int, start0: int, end0: int) -> str:
    local_start = start0 - window_start0
    local_end = end0 - window_start0
    if local_start < 0 or local_end > len(window) or local_end < local_start:
        raise ValueError(
            f"Requested interval {start0}-{end0} lies outside fetched window "
            f"{window_start0}-{window_start0 + len(window)}"
        )
    return window[local_start:local_end]


def _transcript_oriented_interval(
    window: str,
    window_start0: int,
    start0: int,
    end0: int,
    strand: str,
) -> str:
    sequence = _slice_window(window, window_start0, start0, end0)
    return reverse_complement(sequence) if strand == "-" else sequence


def infer_cds_phase_output(
    model: CodingModel,
    *,
    window: str,
    window_start0: int,
    contig_length: int,
) -> CDSPhaseOutput:
    """Infer an output-only first phase without changing biological framing.

    Numeric source phase remains authoritative.  Missing source phase can only
    be standardized to zero, the same effective frame used by the existing
    integrity calculation, and only when sequence termini support that clean
    frame without an equally clean competing supported frame.
    """
    boundary = model.five_prime_boundary_canonicalisation
    if boundary is not None and boundary.status == "canonicalised":
        return CDSPhaseOutput(
            status=CANONICALISED_FIVE_PRIME_BOUNDARY,
            first_phase=0,
            reason="five_prime_boundary_canonicalised",
            source_first_phase=None,
        )
    recovery = model.cds_phase_recovery
    if (
        recovery is not None
        and recovery.status == RECOVERED_MISSING_SOURCE_PHASE
        and recovery.recovered_first_phase is not None
    ):
        return CDSPhaseOutput(
            status=RECOVERED_MISSING_SOURCE_PHASE,
            first_phase=recovery.recovered_first_phase,
            reason=recovery.reason,
            source_first_phase=None,
        )

    source_first_phase = source_first_cds_phase(model)
    if source_first_phase is not None:
        return CDSPhaseOutput(
            status=PRESERVED_NUMERIC_SOURCE,
            first_phase=source_first_phase,
            reason="numeric_source_phase",
            source_first_phase=source_first_phase,
        )

    parts = biological_cds_parts(model)
    if not parts:
        return CDSPhaseOutput(
            status=UNRESOLVED_SOURCE_MISSING,
            reason="no_cds_parts",
        )
    try:
        cds_sequence = "".join(
            _transcript_oriented_interval(
                window, window_start0, start0, end0, model.strand
            )
            for start0, end0 in parts
        )
    except ValueError:
        return CDSPhaseOutput(
            status=UNRESOLVED_SOURCE_MISSING,
            reason="sequence_context_unavailable",
        )
    if not cds_sequence:
        return CDSPhaseOutput(
            status=UNRESOLVED_SOURCE_MISSING,
            reason="coding_sequence_not_evaluable",
        )

    if model.strand == "+":
        adjacent_start0, adjacent_end0 = model.cds_end0, model.cds_end0 + 3
    else:
        adjacent_start0, adjacent_end0 = model.cds_start0 - 3, model.cds_start0
    adjacent_stop = ""
    if adjacent_start0 >= 0 and adjacent_end0 <= contig_length:
        try:
            adjacent_stop = _transcript_oriented_interval(
                window,
                window_start0,
                adjacent_start0,
                adjacent_end0,
                model.strand,
            )
        except ValueError:
            adjacent_stop = ""

    frames = {
        phase: analyse_cds_frame(cds_sequence, adjacent_stop, phase)
        for phase in (0, 1, 2)
    }
    phase_zero = frames[0]
    if not phase_zero.definitely_supported:
        supported_nonzero = [
            phase for phase in (1, 2)
            if frames[phase].supported_for_output_competition
        ]
        reason = (
            "only_nonzero_frame_supported"
            if supported_nonzero
            else "no_informative_terminus"
        )
        return CDSPhaseOutput(
            status=UNRESOLVED_SOURCE_MISSING,
            reason=reason,
        )
    if phase_zero.definite_internal_stop_count:
        return CDSPhaseOutput(
            status=UNRESOLVED_SOURCE_MISSING,
            reason="phase_zero_internal_stop",
        )

    competing_supported = [
        phase for phase in (1, 2)
        if frames[phase].supported_for_output_competition
    ]
    if any(
        frames[phase].definite_internal_stop_count == 0
        for phase in competing_supported
    ):
        return CDSPhaseOutput(
            status=UNRESOLVED_SOURCE_MISSING,
            reason="competing_supported_clean_frame",
        )

    start_supported = phase_zero.definite_start
    stop_supported = phase_zero.definite_stop
    if start_supported and stop_supported:
        reason = "phase_zero_start_and_stop_clean"
    elif start_supported:
        reason = "phase_zero_start_only_clean"
    else:
        reason = "phase_zero_stop_only_clean"
    if competing_supported:
        reason = "phase_zero_uniquely_clean_among_supported_frames"
    return CDSPhaseOutput(
        status=STANDARDISED_SOURCE_MISSING_TO_ZERO,
        first_phase=0,
        reason=reason,
    )


def unresolved_cds_phase_output(model: CodingModel, reason: str) -> CDSPhaseOutput:
    """Preserve numeric source phase or record why missing phase was unresolved."""
    boundary = model.five_prime_boundary_canonicalisation
    if boundary is not None and boundary.status == "canonicalised":
        return CDSPhaseOutput(
            status=CANONICALISED_FIVE_PRIME_BOUNDARY,
            first_phase=0,
            reason="five_prime_boundary_canonicalised",
            source_first_phase=None,
        )
    recovery = model.cds_phase_recovery
    if (
        recovery is not None
        and recovery.status == RECOVERED_MISSING_SOURCE_PHASE
        and recovery.recovered_first_phase is not None
    ):
        return CDSPhaseOutput(
            status=RECOVERED_MISSING_SOURCE_PHASE,
            first_phase=recovery.recovered_first_phase,
            reason=recovery.reason,
            source_first_phase=None,
        )
    source_first_phase = source_first_cds_phase(model)
    if source_first_phase is not None:
        return CDSPhaseOutput(
            status=PRESERVED_NUMERIC_SOURCE,
            first_phase=source_first_phase,
            reason="numeric_source_phase",
            source_first_phase=source_first_phase,
        )
    return CDSPhaseOutput(
        status=UNRESOLVED_SOURCE_MISSING,
        reason=reason,
    )


def _default_identity_row(model: CodingModel) -> dict[str, object]:
    return {
        "superlocus_id": model.superlocus_id or ".",
        "sublocus_id": model.sublocus_id or ".",
        "source_name": model.source_name,
        "input_file": str(model.input_file),
        "input_gene_id": model.input_gene_id,
        "input_transcript_id": model.input_transcript_id,
        "output_gene_id": model.output_gene_id,
        "output_transcript_id": model.output_transcript_id,
        "chrom": model.chrom,
        "strand": model.strand,
        "cds_start": model.cds_start0 + 1,
        "cds_end": model.cds_end0,
        "cds_length": model.cds_length,
        "cds_exon_count": model.cds_exon_count,
        "cds_intron_count": len(model.cds_introns),
    }


def unavailable_sequence_integrity_row(model: CodingModel, status: str, flags: str) -> dict[str, object]:
    row = _default_identity_row(model)
    phase_output = unresolved_cds_phase_output(model, flags or "sequence_unavailable")
    row.update(
        {
            "sequence_integrity_status": status,
            "sequence_start_codon": ".",
            "sequence_start_codon_status": "not_evaluable",
            "sequence_stop_codon": ".",
            "sequence_stop_codon_status": "not_evaluable",
            "sequence_internal_stop_count": 0,
            "sequence_phase_status": "not_evaluable",
            "sequence_phase_mismatch_count": 0,
            "cds_phase_output_status": phase_output.status,
            "cds_phase_output_first_phase": (
                phase_output.first_phase
                if phase_output.first_phase is not None
                else "."
            ),
            "cds_phase_output_reason": phase_output.reason,
            "sequence_noncanonical_splice_intron_count": 0,
            "sequence_noncanonical_splice_intron_fraction": 0.0,
            "sequence_noncanonical_splice_details": ".",
            "sequence_has_ambiguous_bases": "false",
            "sequence_integrity_penalty_score": 0.0,
            "sequence_integrity_penalty_components": ".",
            "sequence_integrity_flags": flags,
        }
    )
    return row


def calculate_model_sequence_integrity(
    model: CodingModel,
    *,
    window: str,
    window_start0: int,
    contig_length: int,
    genetic_code: int = 1,
    internal_stop_penalty: float = 1.0,
    internal_stop_count_cap: int = 2,
    phase_inconsistent_penalty: float = 0.75,
    missing_stop_penalty: float = 0.50,
    missing_start_penalty: float = 0.30,
    noncanonical_splice_penalty: float = 1.00,
) -> dict[str, object]:
    if genetic_code != 1:
        raise ValueError(f"{PROGRAM_NAME} {PROGRAM_VERSION} currently supports --genetic-code 1 only")

    row = _default_identity_row(model)
    phases_by_part = _model_cds_phases(model)
    parts = list(model.cds_parts)
    transcript_parts = list(reversed(parts)) if model.strand == "-" else parts
    segment_sequences = [
        _transcript_oriented_interval(window, window_start0, start0, end0, model.strand)
        for start0, end0 in transcript_parts
    ]
    transcript_phases = [phases_by_part.get(part) for part in transcript_parts]

    first_phase = transcript_phases[0] if transcript_phases else None
    translation_offset = first_phase if first_phase is not None else 0
    if not segment_sequences or translation_offset > len(segment_sequences[0]):
        return unavailable_sequence_integrity_row(model, "not_evaluable_invalid_cds", "invalid_cds_sequence")

    coding_sequence = segment_sequences[0][translation_offset:] + "".join(segment_sequences[1:])
    has_ambiguous = any(base not in "ACGT" for base in coding_sequence)

    present_phase_count = sum(phase is not None for phase in transcript_phases)
    phase_mismatch_count = 0
    cumulative_length = len(segment_sequences[0]) - translation_offset if segment_sequences else 0
    for phase, segment in zip(transcript_phases[1:], segment_sequences[1:]):
        expected = (3 - (cumulative_length % 3)) % 3
        if phase is not None and phase != expected:
            phase_mismatch_count += 1
        cumulative_length += len(segment)
    if phase_mismatch_count:
        phase_status = "inconsistent"
    elif present_phase_count == 0:
        phase_status = "not_available"
    elif present_phase_count < len(transcript_phases):
        phase_status = "partially_missing"
    else:
        phase_status = "consistent"

    five_prime_boundary = (
        model.cds_start0 == 0 if model.strand == "+" else model.cds_end0 == contig_length
    )
    three_prime_boundary = (
        model.cds_end0 == contig_length if model.strand == "+" else model.cds_start0 == 0
    )

    first_codon = coding_sequence[:3] if len(coding_sequence) >= 3 else coding_sequence
    if len(first_codon) < 3:
        start_status = "not_evaluable_at_contig_boundary" if five_prime_boundary else "ambiguous"
    elif any(base not in "ACGT" for base in first_codon):
        start_status = "not_evaluable_at_contig_boundary" if five_prime_boundary else "ambiguous"
        has_ambiguous = True
    elif first_codon in STANDARD_START_CODONS:
        start_status = "present"
    elif five_prime_boundary:
        start_status = "not_evaluable_at_contig_boundary"
    else:
        start_status = "missing"

    complete_codon_length = (len(coding_sequence) // 3) * 3
    codons = [coding_sequence[i : i + 3] for i in range(0, complete_codon_length, 3)]
    terminal_stop_included = bool(
        len(coding_sequence) % 3 == 0 and codons and codons[-1] in STANDARD_STOP_CODONS
    )

    adjacent_stop = ""
    stop_status: str
    stop_codon = codons[-1] if terminal_stop_included else "."
    if terminal_stop_included:
        stop_status = "included_in_cds"
    elif three_prime_boundary:
        stop_status = "not_evaluable_at_contig_boundary"
    elif len(coding_sequence) % 3 != 0:
        stop_status = "missing"
    else:
        if model.strand == "+":
            adjacent_start0, adjacent_end0 = model.cds_end0, model.cds_end0 + 3
        else:
            adjacent_start0, adjacent_end0 = model.cds_start0 - 3, model.cds_start0
        if adjacent_start0 < 0 or adjacent_end0 > contig_length:
            stop_status = "not_evaluable_at_contig_boundary"
        else:
            adjacent_stop = _transcript_oriented_interval(
                window,
                window_start0,
                adjacent_start0,
                adjacent_end0,
                model.strand,
            )
            if any(base not in "ACGT" for base in adjacent_stop):
                stop_status = "ambiguous"
                stop_codon = adjacent_stop or "."
                has_ambiguous = True
            elif adjacent_stop in STANDARD_STOP_CODONS:
                stop_status = "adjacent_to_cds"
                stop_codon = adjacent_stop
            else:
                stop_status = "missing"
                stop_codon = adjacent_stop or "."

    internal_codons = codons[:-1] if terminal_stop_included else codons
    internal_stop_count = sum(codon in STANDARD_STOP_CODONS for codon in internal_codons)

    noncanonical_details: list[str] = []
    for intron_start0, intron_end0 in model.cds_introns:
        if intron_end0 - intron_start0 < 2:
            donor = _transcript_oriented_interval(
                window,
                window_start0,
                intron_start0,
                intron_end0,
                model.strand,
            )
            acceptor = donor
        elif model.strand == "+":
            donor = _slice_window(window, window_start0, intron_start0, intron_start0 + 2)
            acceptor = _slice_window(window, window_start0, intron_end0 - 2, intron_end0)
        else:
            donor = reverse_complement(
                _slice_window(window, window_start0, intron_end0 - 2, intron_end0)
            )
            acceptor = reverse_complement(
                _slice_window(window, window_start0, intron_start0, intron_start0 + 2)
            )
        motif = f"{donor}-{acceptor}"
        if any(base not in "ACGT" for base in donor + acceptor):
            has_ambiguous = True
        if motif not in CANONICAL_SPLICE_MOTIFS:
            noncanonical_details.append(f"{intron_start0 + 1}-{intron_end0}:{motif}")

    intron_count = len(model.cds_introns)
    noncanonical_count = len(noncanonical_details)
    noncanonical_fraction = noncanonical_count / intron_count if intron_count else 0.0

    penalty_components: list[str] = []
    penalty_score = 0.0
    internal_stop_component = min(internal_stop_count, internal_stop_count_cap) * internal_stop_penalty
    if internal_stop_component > 0:
        penalty_score += internal_stop_component
        penalty_components.append(
            f"internal_stop_count:{internal_stop_count}*{internal_stop_penalty:.3f}_cap{internal_stop_count_cap}"
        )
    if phase_status == "inconsistent":
        penalty_score += phase_inconsistent_penalty
        penalty_components.append(f"phase_inconsistent:1*{phase_inconsistent_penalty:.3f}")
    if stop_status == "missing":
        penalty_score += missing_stop_penalty
        penalty_components.append(f"missing_stop:1*{missing_stop_penalty:.3f}")
    if start_status == "missing":
        penalty_score += missing_start_penalty
        penalty_components.append(f"missing_start:1*{missing_start_penalty:.3f}")
    if noncanonical_fraction > 0:
        contribution = noncanonical_fraction * noncanonical_splice_penalty
        penalty_score += contribution
        penalty_components.append(
            f"noncanonical_splice_fraction:{noncanonical_fraction:.6f}*{noncanonical_splice_penalty:.3f}"
        )

    flags: list[str] = []
    if start_status == "missing":
        flags.append("missing_start_codon")
    if stop_status == "missing":
        flags.append("missing_stop_codon")
    if internal_stop_count:
        flags.append("internal_stop_codon")
    if phase_status == "inconsistent":
        flags.append("phase_inconsistent")
    if noncanonical_count:
        flags.append("noncanonical_splice")
        if intron_count and noncanonical_count == intron_count:
            flags.append("all_splices_noncanonical")
    if has_ambiguous:
        flags.append("ambiguous_sequence")
    if five_prime_boundary:
        flags.append("five_prime_contig_boundary")
    if three_prime_boundary:
        flags.append("three_prime_contig_boundary")
    if stop_status == "adjacent_to_cds":
        flags.append("terminal_stop_adjacent_to_cds")

    row.update(
        {
            "sequence_integrity_status": "evaluated_with_concern" if penalty_score > 0 else "evaluated_no_concern",
            "sequence_start_codon": first_codon or ".",
            "sequence_start_codon_status": start_status,
            "sequence_stop_codon": stop_codon,
            "sequence_stop_codon_status": stop_status,
            "sequence_internal_stop_count": internal_stop_count,
            "sequence_phase_status": phase_status,
            "sequence_phase_mismatch_count": phase_mismatch_count,
            "sequence_noncanonical_splice_intron_count": noncanonical_count,
            "sequence_noncanonical_splice_intron_fraction": noncanonical_fraction,
            "sequence_noncanonical_splice_details": ",".join(noncanonical_details) if noncanonical_details else ".",
            "sequence_has_ambiguous_bases": "true" if has_ambiguous else "false",
            "sequence_integrity_penalty_score": penalty_score,
            "sequence_integrity_penalty_components": ",".join(penalty_components) if penalty_components else ".",
            "sequence_integrity_flags": ",".join(flags) if flags else ".",
        }
    )
    return row


def calculate_sequence_integrity_rows(
    superloci: Iterable[Superlocus],
    fasta: IndexedFasta,
    *,
    strict: bool = False,
    genetic_code: int = 1,
    internal_stop_penalty: float = 1.0,
    internal_stop_count_cap: int = 2,
    phase_inconsistent_penalty: float = 0.75,
    missing_stop_penalty: float = 0.50,
    missing_start_penalty: float = 0.30,
    noncanonical_splice_penalty: float = 1.00,
    progress=None,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    superlocus_list = list(superloci)
    for index, superlocus in enumerate(superlocus_list, start=1):
        contig_length = fasta.contig_length(superlocus.chrom)
        if contig_length is None:
            message = f"Genome FASTA has no contig matching annotation seqid {superlocus.chrom}"
            if strict:
                raise ValueError(message)
            for model in superlocus.models:
                rows.append(
                    unavailable_sequence_integrity_row(
                        model,
                        "not_evaluable_missing_contig",
                        "genome_contig_missing",
                    )
                )
            continue
        if superlocus.start0 < 0 or superlocus.end0 > contig_length:
            message = (
                f"Annotation superlocus {superlocus.superlocus_id} lies outside FASTA contig "
                f"{superlocus.chrom}: {superlocus.start0}-{superlocus.end0}; contig_length={contig_length}"
            )
            if strict:
                raise ValueError(message)
            for model in superlocus.models:
                rows.append(
                    unavailable_sequence_integrity_row(
                        model,
                        "not_evaluable_out_of_range",
                        "genome_coordinates_out_of_range",
                    )
                )
            continue

        window_start0 = max(0, superlocus.start0 - 3)
        window_end0 = min(contig_length, superlocus.end0 + 3)
        window = fasta.fetch(superlocus.chrom, window_start0, window_end0)
        for model in superlocus.models:
            try:
                rows.append(
                    calculate_model_sequence_integrity(
                        model,
                        window=window,
                        window_start0=window_start0,
                        contig_length=contig_length,
                        genetic_code=genetic_code,
                        internal_stop_penalty=internal_stop_penalty,
                        internal_stop_count_cap=internal_stop_count_cap,
                        phase_inconsistent_penalty=phase_inconsistent_penalty,
                        missing_stop_penalty=missing_stop_penalty,
                        missing_start_penalty=missing_start_penalty,
                        noncanonical_splice_penalty=noncanonical_splice_penalty,
                    )
                )
            except (ValueError, IndexError) as exc:
                if strict:
                    raise
                rows.append(
                    unavailable_sequence_integrity_row(
                        model,
                        "not_evaluable_invalid_sequence_context",
                        f"sequence_context_error:{type(exc).__name__}",
                    )
                )
        if progress is not None:
            progress.every(
                index,
                len(superlocus_list),
                "Genome sequence integrity superloci",
                level="verbose",
            )
    return rows


SEQUENCE_INTEGRITY_MODEL_COLUMNS = [
    "superlocus_id",
    "sublocus_id",
    "source_name",
    "input_file",
    "input_gene_id",
    "input_transcript_id",
    "output_gene_id",
    "output_transcript_id",
    "chrom",
    "strand",
    "cds_start",
    "cds_end",
    "cds_length",
    "cds_exon_count",
    "cds_intron_count",
    "sequence_integrity_status",
    "sequence_start_codon",
    "sequence_start_codon_status",
    "sequence_stop_codon",
    "sequence_stop_codon_status",
    "sequence_internal_stop_count",
    "sequence_phase_status",
    "sequence_phase_mismatch_count",
    "cds_phase_output_status",
    "cds_phase_output_first_phase",
    "cds_phase_output_reason",
    "sequence_noncanonical_splice_intron_count",
    "sequence_noncanonical_splice_intron_fraction",
    "sequence_noncanonical_splice_details",
    "sequence_has_ambiguous_bases",
    "sequence_integrity_penalty_score",
    "sequence_integrity_penalty_components",
    "sequence_integrity_flags",
]


def write_sequence_integrity_model_metrics_tsv(rows: list[dict[str, object]], path: Path) -> None:
    ordered = sorted(
        rows,
        key=lambda row: (
            str(row.get("chrom", ".")),
            int(row.get("cds_start", 0)),
            int(row.get("cds_end", 0)),
            str(row.get("source_name", ".")),
            str(row.get("input_transcript_id", ".")),
        ),
    )
    write_dict_tsv(ordered, path, SEQUENCE_INTEGRITY_MODEL_COLUMNS)


def _sequence_metric_fraction(count: int | float, denominator: int | float) -> object:
    return (float(count) / float(denominator)) if denominator else "."


def _sequence_summary_rows_from_stats(
    stats_by_source: dict[str, Counter[str]],
    penalty_values_by_source: dict[str, list[float]],
    *,
    fasta_path: Path,
    index_path: Path,
    index_built: bool,
    index_action: Optional[str] = None,
    index_rebuild_reason: str = ".",
    genetic_code: int,
    penalty_active: bool,
) -> list[dict[str, object]]:
    """Build compact all-source and per-source sequence-integrity summaries."""
    summary: list[dict[str, object]] = []

    def add(
        scope: str,
        source_name: str,
        metric: str,
        value: object,
        fraction: object = ".",
        details: object = ".",
    ) -> None:
        summary.append(
            {
                "scope": scope,
                "source_name": source_name,
                "metric": metric,
                "value": value,
                "fraction": fraction,
                "details": details,
            }
        )

    add("configuration", ".", "genome_fasta", str(fasta_path))
    add("configuration", ".", "genome_fasta_index", str(index_path))
    add(
        "configuration",
        ".",
        "genome_fasta_index_built_automatically",
        "true" if index_built else "false",
    )
    add(
        "configuration",
        ".",
        "genome_fasta_index_action",
        index_action or ("created_or_rebuilt" if index_built else "reused"),
        details=index_rebuild_reason,
    )
    add("configuration", ".", "genetic_code", genetic_code)
    add("configuration", ".", "sequence_integrity_penalty_active", "true" if penalty_active else "false")
    add(
        "notes",
        ".",
        "flag_counts_nonexclusive",
        "true",
        details="A model may contribute to more than one individual sequence-integrity flag count.",
    )

    ordered_sources = ["."] + sorted(source for source in stats_by_source if source != ".")
    for source_name in ordered_sources:
        counters = stats_by_source.get(source_name, Counter())
        scope = "all_sources" if source_name == "." else "source"
        total = counters["models_total"]
        evaluated = counters["models_evaluated"]
        metric_counts = [
            ("models_total", total),
            ("models_evaluated", evaluated),
            ("models_evaluated_clean", counters["models_evaluated_clean"]),
            ("models_evaluated_with_concern", counters["models_evaluated_with_concern"]),
            ("models_not_evaluable", counters["models_not_evaluable"]),
            ("models_with_penalty_gt_0", counters["models_with_penalty_gt_0"]),
            ("models_with_penalty_applied", counters["models_with_penalty_applied"]),
            ("models_with_missing_start", counters["models_with_missing_start"]),
            ("models_with_missing_stop", counters["models_with_missing_stop"]),
            ("models_with_internal_stop", counters["models_with_internal_stop"]),
            ("models_with_phase_inconsistent", counters["models_with_phase_inconsistent"]),
            ("models_with_noncanonical_splice", counters["models_with_noncanonical_splice"]),
            ("models_with_all_splices_noncanonical", counters["models_with_all_splices_noncanonical"]),
            ("models_with_ambiguous_sequence", counters["models_with_ambiguous_sequence"]),
            ("models_at_5prime_contig_boundary", counters["models_at_5prime_contig_boundary"]),
            ("models_at_3prime_contig_boundary", counters["models_at_3prime_contig_boundary"]),
        ]
        for metric, count in metric_counts:
            denominator = total
            add(scope, source_name, metric, count, _sequence_metric_fraction(count, denominator))

        assessable_introns = counters["total_assessable_splice_intron_count"]
        noncanonical_introns = counters["noncanonical_splice_intron_count"]
        add(
            scope,
            source_name,
            "total_assessable_splice_intron_count",
            assessable_introns,
            ".",
        )
        add(
            scope,
            source_name,
            "noncanonical_splice_intron_count",
            noncanonical_introns,
            _sequence_metric_fraction(noncanonical_introns, assessable_introns),
        )
        add(
            scope,
            source_name,
            "noncanonical_splice_intron_fraction",
            _sequence_metric_fraction(noncanonical_introns, assessable_introns),
        )

        penalties = penalty_values_by_source.get(source_name, [])
        if penalties:
            add(scope, source_name, "mean_sequence_integrity_penalty_score", statistics.fmean(penalties))
            add(scope, source_name, "median_sequence_integrity_penalty_score", statistics.median(penalties))
            add(scope, source_name, "max_sequence_integrity_penalty_score", max(penalties))
        else:
            add(scope, source_name, "mean_sequence_integrity_penalty_score", 0.0)
            add(scope, source_name, "median_sequence_integrity_penalty_score", 0.0)
            add(scope, source_name, "max_sequence_integrity_penalty_score", 0.0)

        for key in sorted(k for k in counters if k.startswith("status::")):
            count = counters[key]
            add(
                scope,
                source_name,
                f"status_{key.split('::', 1)[1]}",
                count,
                _sequence_metric_fraction(count, total),
            )
    return summary


def sequence_integrity_summary_rows(
    rows: list[dict[str, object]],
    *,
    fasta_path: Path,
    index_path: Path,
    index_built: bool,
    genetic_code: int,
    penalty_active: bool = True,
    index_action: Optional[str] = None,
    index_rebuild_reason: str = ".",
) -> list[dict[str, object]]:
    """Build summary rows from an in-memory set of model metric rows."""
    stats_by_source: dict[str, Counter[str]] = defaultdict(Counter)
    penalty_values_by_source: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        source = str(row.get("source_name", "."))
        for key in (".", source):
            _update_sequence_counters(stats_by_source[key], row, penalty_active=penalty_active)
            penalty_values_by_source[key].append(float(row.get("sequence_integrity_penalty_score", 0.0) or 0.0))
    return _sequence_summary_rows_from_stats(
        stats_by_source,
        penalty_values_by_source,
        fasta_path=fasta_path,
        index_path=index_path,
        index_built=index_built,
        index_action=index_action,
        index_rebuild_reason=index_rebuild_reason,
        genetic_code=genetic_code,
        penalty_active=penalty_active,
    )


SEQUENCE_INTEGRITY_SUMMARY_COLUMNS = [
    "scope", "source_name", "metric", "value", "fraction", "details"
]


def write_sequence_integrity_summary_tsv(rows: list[dict[str, object]], path: Path) -> None:
    write_dict_tsv(rows, path, SEQUENCE_INTEGRITY_SUMMARY_COLUMNS)


def _format_count_fraction(count: int, total: int) -> str:
    fraction = (count / total) if total else 0.0
    return f"{count:,}({fraction:.2%})"


def sequence_integrity_log_lines_from_stats(
    stats_by_source: dict[str, Counter[str]],
) -> list[str]:
    overall = stats_by_source.get(".", Counter())
    lines = [
        (
            "Genome sequence integrity: "
            f"evaluated={overall['models_evaluated']:,}; "
            f"penalized={overall['models_with_penalty_gt_0']:,}; "
            f"internal_stop={overall['models_with_internal_stop']:,}; "
            f"missing_start={overall['models_with_missing_start']:,}; "
            f"missing_stop={overall['models_with_missing_stop']:,}; "
            f"phase_inconsistent={overall['models_with_phase_inconsistent']:,}; "
            f"noncanonical_splice={overall['models_with_noncanonical_splice']:,}; "
            f"all_splices_noncanonical={overall['models_with_all_splices_noncanonical']:,}"
        ),
        "Sequence integrity by source (model fractions use source model totals; individual flag counts are non-exclusive):",
        (
            "source\tmodels\tevaluated\tclean\twith_concern\tpenalty_applied\t"
            "missing_start\tmissing_stop\tinternal_stop\tnoncanonical_splice\t"
            "all_splices_noncanonical\tambiguous_sequence\tnoncanonical_introns/assessable"
        ),
    ]
    for source in sorted(key for key in stats_by_source if key != "."):
        counters = stats_by_source[source]
        total = counters["models_total"]
        assessable = counters["total_assessable_splice_intron_count"]
        noncanonical_introns = counters["noncanonical_splice_intron_count"]
        intron_text = (
            f"{noncanonical_introns:,}/{assessable:,}"
            f"({(noncanonical_introns / assessable if assessable else 0.0):.2%})"
        )
        lines.append(
            "\t".join(
                [
                    source,
                    f"{total:,}",
                    _format_count_fraction(counters["models_evaluated"], total),
                    _format_count_fraction(counters["models_evaluated_clean"], total),
                    _format_count_fraction(counters["models_evaluated_with_concern"], total),
                    _format_count_fraction(counters["models_with_penalty_applied"], total),
                    _format_count_fraction(counters["models_with_missing_start"], total),
                    _format_count_fraction(counters["models_with_missing_stop"], total),
                    _format_count_fraction(counters["models_with_internal_stop"], total),
                    _format_count_fraction(counters["models_with_noncanonical_splice"], total),
                    _format_count_fraction(counters["models_with_all_splices_noncanonical"], total),
                    _format_count_fraction(counters["models_with_ambiguous_sequence"], total),
                    intron_text,
                ]
            )
        )
    return lines


def sequence_integrity_log_lines(rows: list[dict[str, object]]) -> list[str]:
    stats_by_source: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        source = str(row.get("source_name", "."))
        for key in (".", source):
            _update_sequence_counters(stats_by_source[key], row, penalty_active=True)
    return sequence_integrity_log_lines_from_stats(stats_by_source)

def sequence_integrity_metrics_from_row(row: dict[str, object]) -> SequenceIntegrityMetrics:
    """Convert one transient output row to compact model-attached metrics."""
    return SequenceIntegrityMetrics(
        status=str(row.get("sequence_integrity_status", "not_evaluable")),
        start_codon=str(row.get("sequence_start_codon", ".")),
        start_codon_status=str(row.get("sequence_start_codon_status", "not_evaluable")),
        stop_codon=str(row.get("sequence_stop_codon", ".")),
        stop_codon_status=str(row.get("sequence_stop_codon_status", "not_evaluable")),
        internal_stop_count=int(row.get("sequence_internal_stop_count", 0) or 0),
        phase_status=str(row.get("sequence_phase_status", "not_evaluable")),
        phase_mismatch_count=int(row.get("sequence_phase_mismatch_count", 0) or 0),
        noncanonical_splice_intron_count=int(row.get("sequence_noncanonical_splice_intron_count", 0) or 0),
        noncanonical_splice_intron_fraction=float(row.get("sequence_noncanonical_splice_intron_fraction", 0.0) or 0.0),
        noncanonical_splice_details=str(row.get("sequence_noncanonical_splice_details", ".")),
        has_ambiguous_bases=str(row.get("sequence_has_ambiguous_bases", "false")).lower() == "true",
        penalty_score=float(row.get("sequence_integrity_penalty_score", 0.0) or 0.0),
        penalty_components=str(row.get("sequence_integrity_penalty_components", ".")),
        flags=str(row.get("sequence_integrity_flags", ".")),
    )


def sequence_integrity_row_from_model(model: CodingModel) -> dict[str, object]:
    """Reconstruct a full dedicated-output row from compact attached metrics."""
    metrics = model.sequence_integrity
    if metrics is None:
        return unavailable_sequence_integrity_row(
            model,
            "genome_not_configured",
            ".",
        )
    row = _default_identity_row(model)
    row.update(
        {
            "sequence_integrity_status": metrics.status,
            "sequence_start_codon": metrics.start_codon,
            "sequence_start_codon_status": metrics.start_codon_status,
            "sequence_stop_codon": metrics.stop_codon,
            "sequence_stop_codon_status": metrics.stop_codon_status,
            "sequence_internal_stop_count": metrics.internal_stop_count,
            "sequence_phase_status": metrics.phase_status,
            "sequence_phase_mismatch_count": metrics.phase_mismatch_count,
            "cds_phase_output_status": (
                model.cds_phase_output.status
                if model.cds_phase_output is not None
                else UNRESOLVED_SOURCE_MISSING
            ),
            "cds_phase_output_first_phase": (
                model.cds_phase_output.first_phase
                if model.cds_phase_output is not None
                and model.cds_phase_output.first_phase is not None
                else "."
            ),
            "cds_phase_output_reason": (
                model.cds_phase_output.reason
                if model.cds_phase_output is not None
                else "not_evaluated"
            ),
            "sequence_noncanonical_splice_intron_count": metrics.noncanonical_splice_intron_count,
            "sequence_noncanonical_splice_intron_fraction": metrics.noncanonical_splice_intron_fraction,
            "sequence_noncanonical_splice_details": metrics.noncanonical_splice_details,
            "sequence_has_ambiguous_bases": "true" if metrics.has_ambiguous_bases else "false",
            "sequence_integrity_penalty_score": metrics.penalty_score,
            "sequence_integrity_penalty_components": metrics.penalty_components,
            "sequence_integrity_flags": metrics.flags,
        }
    )
    return row


def _update_sequence_counters(
    counters: Counter[str],
    row: dict[str, object],
    *,
    penalty_active: bool,
) -> None:
    status = str(row.get("sequence_integrity_status", "."))
    counters["models_total"] += 1
    counters[f"status::{status}"] += 1
    if status.startswith("evaluated_"):
        counters["models_evaluated"] += 1
    else:
        counters["models_not_evaluable"] += 1
    if status == "evaluated_no_concern":
        counters["models_evaluated_clean"] += 1
    if status == "evaluated_with_concern":
        counters["models_evaluated_with_concern"] += 1
    penalty_score = float(row.get("sequence_integrity_penalty_score", 0.0) or 0.0)
    if penalty_score > 0:
        counters["models_with_penalty_gt_0"] += 1
        if penalty_active:
            counters["models_with_penalty_applied"] += 1
    if int(row.get("sequence_internal_stop_count", 0) or 0) > 0:
        counters["models_with_internal_stop"] += 1
    if row.get("sequence_start_codon_status") == "missing":
        counters["models_with_missing_start"] += 1
    if row.get("sequence_stop_codon_status") == "missing":
        counters["models_with_missing_stop"] += 1
    if row.get("sequence_phase_status") == "inconsistent":
        counters["models_with_phase_inconsistent"] += 1
    noncanonical_count = int(row.get("sequence_noncanonical_splice_intron_count", 0) or 0)
    if noncanonical_count > 0:
        counters["models_with_noncanonical_splice"] += 1
    counters["noncanonical_splice_intron_count"] += noncanonical_count
    counters["total_assessable_splice_intron_count"] += int(row.get("cds_intron_count", 0) or 0)
    raw_flags = str(row.get("sequence_integrity_flags", "."))
    flags = {flag for flag in raw_flags.split(",") if flag and flag != "."}
    for flag in flags:
        counters[f"flag::{flag}"] += 1
    if "all_splices_noncanonical" in flags:
        counters["models_with_all_splices_noncanonical"] += 1
    if "ambiguous_sequence" in flags:
        counters["models_with_ambiguous_sequence"] += 1
    if "five_prime_contig_boundary" in flags:
        counters["models_at_5prime_contig_boundary"] += 1
    if "three_prime_contig_boundary" in flags:
        counters["models_at_3prime_contig_boundary"] += 1

def calculate_and_attach_sequence_integrity(
    superloci: Iterable[Superlocus],
    fasta: IndexedFasta,
    *,
    strict: bool = False,
    genetic_code: int = 1,
    internal_stop_penalty: float = 1.0,
    internal_stop_count_cap: int = 2,
    phase_inconsistent_penalty: float = 0.75,
    missing_stop_penalty: float = 0.50,
    missing_start_penalty: float = 0.30,
    noncanonical_splice_penalty: float = 1.00,
    penalty_active: bool = True,
    progress=None,
) -> tuple[list[dict[str, object]], list[str]]:
    """Calculate metrics by superlocus and attach compact state to models.

    Only one transient row dictionary exists at a time.  The model-level output
    can later be streamed from the attached slotted metrics, avoiding a second
    full model-sized row collection on million-model genomes.
    """
    stats_by_source: dict[str, Counter[str]] = defaultdict(Counter)
    penalty_values_by_source: dict[str, list[float]] = defaultdict(list)
    superlocus_list = superloci if isinstance(superloci, list) else list(superloci)

    def record(row: dict[str, object]) -> None:
        source = str(row.get("source_name", "."))
        penalty_score = float(row.get("sequence_integrity_penalty_score", 0.0) or 0.0)
        for key in (".", source):
            _update_sequence_counters(stats_by_source[key], row, penalty_active=penalty_active)
            penalty_values_by_source[key].append(penalty_score)
    for index, superlocus in enumerate(superlocus_list, start=1):
        contig_length = fasta.contig_length(superlocus.chrom)
        if contig_length is None:
            message = f"Genome FASTA has no contig matching annotation seqid {superlocus.chrom}"
            if strict:
                raise ValueError(message)
            for model in superlocus.models:
                row = unavailable_sequence_integrity_row(
                    model, "not_evaluable_missing_contig", "genome_contig_missing"
                )
                model.sequence_integrity = sequence_integrity_metrics_from_row(row)
                model.cds_phase_output = unresolved_cds_phase_output(
                    model, "genome_contig_missing"
                )
                record(row)
            continue
        if superlocus.start0 < 0 or superlocus.end0 > contig_length:
            message = (
                f"Annotation superlocus {superlocus.superlocus_id} lies outside FASTA contig "
                f"{superlocus.chrom}: {superlocus.start0}-{superlocus.end0}; contig_length={contig_length}"
            )
            if strict:
                raise ValueError(message)
            for model in superlocus.models:
                row = unavailable_sequence_integrity_row(
                    model, "not_evaluable_out_of_range", "genome_coordinates_out_of_range"
                )
                model.sequence_integrity = sequence_integrity_metrics_from_row(row)
                model.cds_phase_output = unresolved_cds_phase_output(
                    model, "genome_coordinates_out_of_range"
                )
                record(row)
            continue

        window_start0 = max(0, superlocus.start0 - 3)
        window_end0 = min(contig_length, superlocus.end0 + 3)
        window = fasta.fetch(superlocus.chrom, window_start0, window_end0)
        for model in superlocus.models:
            try:
                row = calculate_model_sequence_integrity(
                    model,
                    window=window,
                    window_start0=window_start0,
                    contig_length=contig_length,
                    genetic_code=genetic_code,
                    internal_stop_penalty=internal_stop_penalty,
                    internal_stop_count_cap=internal_stop_count_cap,
                    phase_inconsistent_penalty=phase_inconsistent_penalty,
                    missing_stop_penalty=missing_stop_penalty,
                    missing_start_penalty=missing_start_penalty,
                    noncanonical_splice_penalty=noncanonical_splice_penalty,
                )
            except (ValueError, IndexError) as exc:
                if strict:
                    raise
                row = unavailable_sequence_integrity_row(
                    model,
                    "not_evaluable_invalid_sequence_context",
                    f"sequence_context_error:{type(exc).__name__}",
                )
            model.sequence_integrity = sequence_integrity_metrics_from_row(row)
            try:
                model.cds_phase_output = infer_cds_phase_output(
                    model,
                    window=window,
                    window_start0=window_start0,
                    contig_length=contig_length,
                )
            except (ValueError, IndexError):
                model.cds_phase_output = unresolved_cds_phase_output(
                    model, "sequence_context_unavailable"
                )
            record(row)
        if progress is not None:
            progress.every(
                index,
                len(superlocus_list),
                "Genome sequence integrity superloci",
                level="verbose",
            )

    summary = _sequence_summary_rows_from_stats(
        stats_by_source,
        penalty_values_by_source,
        fasta_path=fasta.path,
        index_path=fasta.index_path,
        index_built=fasta.index_built,
        index_action=fasta.index_action,
        index_rebuild_reason=fasta.index_rebuild_reason,
        genetic_code=genetic_code,
        penalty_active=penalty_active,
    )
    return summary, sequence_integrity_log_lines_from_stats(stats_by_source)


def write_sequence_integrity_model_metrics_from_models(
    models: Iterable[CodingModel],
    path: Path,
) -> None:
    """Stream the dedicated model table from compact attached metrics.

    Sorting creates only a list of model references, rather than a second list
    of large row dictionaries, and preserves the established deterministic TSV
    order.
    """
    ordered_models = sorted(
        models,
        key=lambda model: (
            model.chrom,
            model.cds_start0,
            model.cds_end0,
            model.source_name,
            model.input_transcript_id,
        ),
    )
    with Path(path).open("w") as out:
        out.write("\t".join(SEQUENCE_INTEGRITY_MODEL_COLUMNS) + "\n")
        for model in ordered_models:
            row = sequence_integrity_row_from_model(model)
            out.write(
                "\t".join(
                    format_tsv_value(row.get(column, "."))
                    for column in SEQUENCE_INTEGRITY_MODEL_COLUMNS
                )
                + "\n"
            )
