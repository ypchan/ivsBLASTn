from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

from .logging import LOG
from .models import HSP

def make_blast_db(fasta: Path, db_prefix: Path, makeblastdb_bin: str) -> None:
    """Build nucleotide BLAST database."""

    cmd = [makeblastdb_bin, "-in", str(fasta), "-dbtype", "nucl", "-out", str(db_prefix)]
    LOG.info("Running makeblastdb: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


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
    subprocess.run(cmd, check=True)
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


def parse_blast(path: Path, min_pident: float, min_hsp_len: int) -> Dict[str, Dict[str, List[HSP]]]:
    """Parse BLAST outfmt 6 and group HSPs by query and subject."""

    grouped: Dict[str, Dict[str, List[HSP]]] = defaultdict(lambda: defaultdict(list))
    n_lines = 0
    n_kept = 0
    with path.open("rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=chr(9))
        for parts in reader:
            if not parts or parts[0].startswith("#"):
                continue
            n_lines += 1
            if len(parts) < 10:
                LOG.debug("Skipping BLAST row with <10 columns: %s", parts)
                continue
            pident = float(parts[2])
            length = int(parts[3])
            if pident < min_pident or length < min_hsp_len:
                continue
            hsp = HSP(
                qseqid=parts[0],
                sseqid=parts[1],
                pident=pident,
                length=length,
                qstart=int(parts[4]),
                qend=int(parts[5]),
                sstart=int(parts[6]),
                send=int(parts[7]),
                evalue=parts[8],
                bitscore=float(parts[9]),
            )
            grouped[hsp.qseqid][hsp.sseqid].append(hsp)
            n_kept += 1
    LOG.info("Parsed BLAST HSPs: %s rows, %s kept after filters", n_lines, n_kept)
    return grouped
