from __future__ import annotations

from pathlib import Path
from typing import List

from .fasta import iter_fasta_records, write_fasta_record_with_description
from .logging import LOG


def split_fasta(query: Path, chunks_dir: Path, chunk_size: int, prefix: str = "query") -> List[Path]:
    """Split a FASTA/FASTA.gz file into fixed-record-count chunk FASTA files."""

    chunks_dir.mkdir(parents=True, exist_ok=True)
    chunk_paths: List[Path] = []
    handle = None
    try:
        for record_index, (_record_id, description, seq) in enumerate(iter_fasta_records(query), start=1):
            chunk_index = (record_index - 1) // chunk_size + 1
            if (record_index - 1) % chunk_size == 0:
                if handle is not None:
                    handle.close()
                chunk_path = chunks_dir / f"{prefix}.{chunk_index:06d}.fa"
                chunk_paths.append(chunk_path)
                handle = chunk_path.open("wt", encoding="utf-8")
            if handle is None:
                raise RuntimeError("Failed to open chunk FASTA")
            write_fasta_record_with_description(handle, description, seq)
    finally:
        if handle is not None:
            handle.close()

    manifest = chunks_dir / "chunks.tsv"
    with manifest.open("wt", encoding="utf-8") as out:
        print("chunk_index\tquery_fasta", file=out)
        for index, path in enumerate(chunk_paths, start=1):
            print(f"{index}\t{path.resolve()}", file=out)
    LOG.info("Split %s into %s chunks under %s", query, len(chunk_paths), chunks_dir)
    return chunk_paths
