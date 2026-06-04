from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

from .fasta import open_text_auto
from .logging import LOG
from .models import BlastQueryStats, HSP


def run_external_command(cmd: List[str], label: str) -> None:
    """Run an external command and include stderr in failures."""

    completed = subprocess.run(cmd, text=True, capture_output=True)
    if completed.stdout:
        LOG.debug("%s stdout:\n%s", label, completed.stdout.rstrip())
    if completed.stderr:
        LOG.debug("%s stderr:\n%s", label, completed.stderr.rstrip())
    if completed.returncode != 0:
        details = completed.stderr.strip() or completed.stdout.strip() or "no stderr/stdout captured"
        raise RuntimeError(f"{label} failed with exit code {completed.returncode}:\n{details}")


def make_blast_db(fasta: Path, db_prefix: Path, makeblastdb_bin: str) -> None:
    """Build nucleotide BLAST database."""

    cmd = [makeblastdb_bin, "-in", str(fasta), "-dbtype", "nucl", "-out", str(db_prefix)]
    LOG.info("Running makeblastdb: %s", " ".join(cmd))
    run_external_command(cmd, "makeblastdb")


def run_blastn_to_file(query: Path, db: Path, out_file: Path, args: argparse.Namespace, max_targets: int, max_hsps: int, label: str) -> Path:
    """Run BLASTN and return output table path."""

    cmd = [
        args.blastn_bin,
        "-query",
        str(query),
        "-db",
        str(db),
        "-outfmt",
        "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
        "-max_target_seqs",
        str(max_targets),
        "-max_hsps",
        str(max_hsps),
        "-num_threads",
        str(args.threads),
        "-out",
        str(out_file),
    ]
    if args.blast_task:
        cmd.extend(["-task", args.blast_task])
    if args.blast_evalue:
        cmd.extend(["-evalue", str(args.blast_evalue)])
    LOG.info("Running %s: %s", label, " ".join(cmd))
    run_external_command(cmd, label)
    return out_file


def run_query_blastn(args: argparse.Namespace) -> Path:
    """Run query BLASTN and return table path."""

    return run_blastn_to_file(
        query=args.query,
        db=args.db,
        out_file=args.query_blast,
        args=args,
        max_targets=args.top_subjects,
        max_hsps=args.blast_max_hsps,
        label="BLASTN",
    )


def parse_blast_row(parts: List[str]) -> Optional[HSP]:
    """Parse one BLAST outfmt 6 row, returning None for malformed rows."""

    if len(parts) < 10:
        LOG.debug("Skipping BLAST row with <10 columns: %s", parts)
        return None
    try:
        return HSP(
            qseqid=parts[0],
            sseqid=parts[1],
            pident=float(parts[2]),
            length=int(parts[3]),
            qstart=int(parts[4]),
            qend=int(parts[5]),
            sstart=int(parts[6]),
            send=int(parts[7]),
            evalue=parts[8],
            bitscore=float(parts[9]),
        )
    except ValueError:
        LOG.debug("Skipping malformed BLAST row: %s", parts)
        return None


def iter_blast_hsps(path: Path) -> Iterator[HSP]:
    """Yield parsed HSPs from plain or gzip-compressed BLAST outfmt 6."""

    with open_text_auto(path) as handle:
        reader = csv.reader(handle, delimiter=chr(9))
        for parts in reader:
            if not parts or parts[0].startswith("#"):
                continue
            hsp = parse_blast_row(parts)
            if hsp is not None:
                yield hsp


def parse_blast_with_stats(path: Path, min_pident: float, min_hsp_len: int) -> Tuple[Dict[str, Dict[str, List[HSP]]], Dict[str, BlastQueryStats]]:
    """Parse BLAST outfmt 6 and return retained HSPs plus per-query row counts."""

    grouped: Dict[str, Dict[str, List[HSP]]] = defaultdict(lambda: defaultdict(list))
    stats: Dict[str, BlastQueryStats] = defaultdict(BlastQueryStats)
    raw_subjects: Dict[str, Set[str]] = defaultdict(set)
    retained_subjects: Dict[str, Set[str]] = defaultdict(set)
    n_lines = 0
    n_kept = 0
    n_bad = 0
    with open_text_auto(path) as handle:
        reader = csv.reader(handle, delimiter=chr(9))
        for parts in reader:
            if not parts or parts[0].startswith("#"):
                continue
            n_lines += 1
            if len(parts) >= 2:
                stats[parts[0]].raw_hsps += 1
                raw_subjects[parts[0]].add(parts[1])
            hsp = parse_blast_row(parts)
            if hsp is None:
                n_bad += 1
                continue
            if hsp.pident < min_pident or hsp.length < min_hsp_len:
                continue
            grouped[hsp.qseqid][hsp.sseqid].append(hsp)
            stats[hsp.qseqid].retained_hsps += 1
            retained_subjects[hsp.qseqid].add(hsp.sseqid)
            n_kept += 1
    for query_id, query_stats in stats.items():
        query_stats.raw_subjects = len(raw_subjects.get(query_id, set()))
        query_stats.retained_subjects = len(retained_subjects.get(query_id, set()))
    LOG.info("Parsed BLAST HSPs: %s rows, %s kept after filters, %s malformed skipped", n_lines, n_kept, n_bad)
    return grouped, dict(stats)


def parse_blast(path: Path, min_pident: float, min_hsp_len: int) -> Dict[str, Dict[str, List[HSP]]]:
    """Parse BLAST outfmt 6 and group retained HSPs by query and subject."""

    grouped, _stats = parse_blast_with_stats(path, min_pident, min_hsp_len)
    return grouped
