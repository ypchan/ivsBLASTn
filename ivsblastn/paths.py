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

    args.reference_dir = args.outdir / "reference"
    args.blast_dir = args.outdir / "blast"
    args.results_dir = args.outdir / "results"
    for directory in [args.outdir, args.reference_dir, args.blast_dir, args.results_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    args.query_label = fasta_label(args.query)
    args.raw_ref_fa = args.reference_dir / "raw_reference.fa"
    args.raw_ref_tax = args.reference_dir / "raw_reference.tax.tsv"
    args.raw_ref_db = args.reference_dir / "raw_reference_db"
    args.ref_self_blast = args.reference_dir / "reference_self.blastn.tsv"
    args.ref_self_clean_prefix = args.reference_dir / "reference_self_clean"
    args.ref_self_clean_introns_fa = output_path(args.ref_self_clean_prefix, ".ivs.fa")
    args.cleaned_ref_fa = args.reference_dir / "cleaned_reference.fa"
    args.cleaned_ref_tax = args.reference_dir / "cleaned_reference.tax.tsv"
    args.cleaned_ref_db = args.reference_dir / "cleaned_reference_db"

    args.query_blast = args.blast_dir / f"{args.query_label}.vs_reference.blastn.tsv"
    args.result_prefix = args.results_dir / args.query_label
    args.summary_tsv = output_path(args.result_prefix, ".summary.tsv")
    args.supporting_tsv = output_path(args.result_prefix, ".supporting_hsps.tsv")
    fasta_suffix = ".fa.gz" if args.gzip_fasta_output else ".fa"
    args.intron_free_fa = output_path(args.result_prefix, f".ivs_free{fasta_suffix}")
    args.introns_fa = output_path(args.result_prefix, f".ivs{fasta_suffix}")
    args.exons_bed = output_path(args.result_prefix, ".exons.bed")
    args.introns_bed = output_path(args.result_prefix, ".ivs.bed")
    args.report_md = output_path(args.result_prefix, ".report.md")
