from __future__ import annotations

import argparse
from pathlib import Path

def output_path(prefix: Path, suffix: str) -> Path:
    """Return PREFIX + suffix without replacing dotted suffixes."""

    return Path(str(prefix) + suffix)


def fasta_label(path: Path) -> str:
    """Create a compact label from FASTA/FASTA.gz filename."""

    name = path.name
    for suffix in [".fasta.gz", ".fa.gz", ".fna.gz", ".faa.gz", ".fasta", ".fa", ".fna", ".faa", ".gz"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    chars = [c if (c.isalnum() or c in {".", "_", "-"}) else "_" for c in name]
    return "".join(chars).strip("._-") or "query"


def setup_output_paths(args: argparse.Namespace) -> None:
    """Create output directories and deterministic output file paths."""

    args.blast_dir = args.outdir / "blast"
    args.results_dir = args.outdir / "results"
    for directory in [args.outdir, args.blast_dir, args.results_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    args.query_label = fasta_label(args.query)
    args.query_blast = args.blast_dir / f"{args.query_label}.vs_reference.blastn.tsv"
    args.result_prefix = args.results_dir / args.query_label
    args.summary_tsv = output_path(args.result_prefix, ".summary.tsv")
    args.supporting_tsv = output_path(args.result_prefix, ".supporting_hsps.tsv")
    fasta_suffix = ".fa.gz" if args.gzip_fasta_output else ".fa"
    args.ivs_free_fa = output_path(args.result_prefix, f".ivs_free{fasta_suffix}")
    args.ivs_fa = output_path(args.result_prefix, f".ivs{fasta_suffix}")
    args.exons_bed = output_path(args.result_prefix, ".exons.bed")
    args.ivs_bed = output_path(args.result_prefix, ".ivs.bed")
    args.report_md = output_path(args.result_prefix, ".report.md")
