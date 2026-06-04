from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, Iterable, List, Optional, TextIO, Tuple

def open_text_auto(path: Path) -> TextIO:
    """Open plain text or gzip-compressed text."""

    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def open_text_auto_write(path: Path) -> TextIO:
    """Open plain text or gzip-compressed text for writing."""

    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return path.open("wt", encoding="utf-8")


def wrap_fasta(seq: str, width: int = 80) -> Iterable[str]:
    """Yield wrapped FASTA sequence lines."""

    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def write_fasta_record(handle: TextIO, seq_id: str, seq: str) -> None:
    """Write one FASTA record safely."""

    print(f">{seq_id}", file=handle)
    for part in wrap_fasta(seq):
        print(part, file=handle)


def read_fasta(path: Path) -> Dict[str, str]:
    """Read FASTA or FASTA.gz into {sequence_id: sequence}."""

    seqs: Dict[str, List[str]] = {}
    current_id: Optional[str] = None
    with open_text_auto(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                seqs[current_id] = []
            elif current_id is not None:
                seqs[current_id].append(line)
    return {seq_id: "".join(parts).upper() for seq_id, parts in seqs.items()}


def parse_silva_header(header: str) -> Tuple[str, str]:
    """Parse SILVA FASTA header into subject ID and taxonomy."""

    text = header[1:] if header.startswith(">") else header
    parts = text.strip().split(maxsplit=1)
    seq_id = parts[0]
    taxonomy = parts[1].strip() if len(parts) > 1 else ""
    return seq_id, taxonomy


def taxonomy_domain(taxonomy: str) -> str:
    """Return first taxonomy field."""

    parts = [p.strip() for p in taxonomy.split(";") if p.strip()]
    return parts[0] if parts else "NA"


def species_key_from_taxonomy(taxonomy: str) -> Optional[str]:
    """Return a species-level key when taxonomy contains a clear binomial."""

    parts = [p.strip() for p in taxonomy.split(";") if p.strip()]
    if len(parts) < 7:
        return None
    species_name = parts[6]
    species_words = species_name.split()
    if len(species_words) < 2:
        return None
    species_epithet = species_words[1]
    if any(char.isdigit() for char in species_epithet):
        return None
    return ";".join(parts[:6] + [" ".join(species_words[:2])])


def clean_dna_sequence(seq_parts: List[str]) -> str:
    """Join sequence lines, normalize U to T, and remove whitespace."""

    seq = "".join(seq_parts).upper().replace("U", "T")
    return "".join(seq.split())


def is_strict_atgc(seq: str) -> bool:
    """Return True only if sequence is non-empty and contains only A/T/G/C."""

    return bool(seq) and set(seq) <= {"A", "T", "G", "C"}
