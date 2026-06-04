from __future__ import annotations

import argparse
import concurrent.futures as futures
from pathlib import Path
from typing import List

from rich.progress import BarColumn, Progress, SpinnerColumn, TaskProgressColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn
from rich_argparse import RichHelpFormatter

from .algorithm import analyze_query, confidence_rank
from .blast import parse_blast, run_query_blastn
from .fasta import read_fasta
from .logging import CONSOLE, LOG, setup_logging
from .models import QueryResult
from .outputs import print_summary, write_bed_outputs, write_fasta_outputs, write_report, write_summary, write_supporting_hsps
from .paths import setup_output_paths
from .reference import clean_reference_introns, preprocess_reference
from .taxonomy import parse_taxonomy

def positive_int(value: str) -> int:
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("Value must be positive")
    return ivalue


def nonnegative_int(value: str) -> int:
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("Value must be non-negative")
    return ivalue


def probability_percent(value: str) -> float:
    fvalue = float(value)
    if fvalue < 0 or fvalue > 100:
        raise argparse.ArgumentTypeError("Value must be between 0 and 100")
    return fvalue


def build_parser() -> argparse.ArgumentParser:
    """Create argument parser."""

    parser = argparse.ArgumentParser(
        prog="ivsBLASTn",
        description="ivsBLASTn -- detect intervening sequences by blastn search.",
        formatter_class=RichHelpFormatter,
    )

    required = parser.add_argument_group("Required inputs")
    required.add_argument("--query", required=True, type=Path, help="Query 16S FASTA/FASTA.gz file. No default.")
    required.add_argument("--blast", default=None, type=Path, help="Existing BLASTN outfmt 6 table. Default: not used; provide exactly one of --blast, --db, or --ref-fasta.")
    required.add_argument("--db", default=None, type=Path, help="Existing BLAST database prefix. Default: not used; provide exactly one of --blast, --db, or --ref-fasta.")
    required.add_argument("--ref-fasta", default=None, type=Path, help="SILVA-style reference FASTA/FASTA.gz with taxonomy in headers. Default: not used; provide exactly one of --blast, --db, or --ref-fasta.")
    required.add_argument("--outdir", required=True, type=Path, help="Output directory. No default.")

    algorithm = parser.add_argument_group("Algorithm")
    algorithm.add_argument("--algorithm", default="hsp-gap-support", choices=["hsp-gap-support"], help="Detection algorithm. Default: hsp-gap-support.")

    filters = parser.add_argument_group("Reference preprocessing and HSP filters")
    filters.add_argument("--ref-domains", default="Archaea,Bacteria", help="Comma-separated SILVA domains retained with --ref-fasta. Default: Archaea,Bacteria.")
    filters.add_argument("--ref-per-species", default=1, type=nonnegative_int, help="Maximum sequences per clear species in reference preprocessing, keeping the longest sequences first. Use 0 to disable. Default: 1.")
    filters.add_argument("--clean-ref-introns", action="store_true", help="Self-BLAST reference and remove candidate introns before query BLAST. Default: disabled.")
    filters.add_argument("--ref-clean-min-confidence", default="LOW", choices=["LOW", "MEDIUM", "HIGH"], help="Minimum confidence required to remove a reference intron. Default: LOW.")
    filters.add_argument("--ref-self-blast-max-target-seqs", default=100, type=positive_int, help="Reference self-BLAST -max_target_seqs. Default: 100.")
    filters.add_argument("--ref-self-blast-max-hsps", default=20, type=positive_int, help="Reference self-BLAST -max_hsps. Default: 20.")
    filters.add_argument("--min-pident", default=75.0, type=probability_percent, help="Minimum HSP percent identity. Default: 75.0.")
    filters.add_argument("--min-hsp-len", default=100, type=positive_int, help="Minimum HSP length in bp. Default: 100.")
    filters.add_argument("--top-subjects", default=100, type=positive_int, help="Top subjects retained per query after parsing BLAST. Default: 100.")
    filters.add_argument("--makeblastdb-bin", default="makeblastdb", help="makeblastdb executable. Default: makeblastdb.")
    filters.add_argument("--blastn-bin", default="blastn", help="blastn executable. Default: blastn.")
    filters.add_argument("--blast-max-target-seqs", default=100, type=positive_int, help="Query BLASTN -max_target_seqs. Default: 100.")
    filters.add_argument("--blast-max-hsps", default=20, type=positive_int, help="Query BLASTN -max_hsps. Default: 20.")
    filters.add_argument("--blast-task", default="blastn", choices=["blastn", "megablast", "dc-megablast", "blastn-short"], help="BLASTN task. Default: blastn.")
    filters.add_argument("--blast-evalue", default="1e-20", help="BLASTN e-value. Default: 1e-20.")

    intron = parser.add_argument_group("Candidate intron geometry")
    intron.add_argument("--min-intron-len", default=25, type=nonnegative_int, help="Minimum query gap size. Default: 25 bp.")
    intron.add_argument("--max-intron-len", default=2000, type=positive_int, help="Maximum query gap size. Default: 2000 bp.")
    intron.add_argument("--max-ref-gap", default=30, type=nonnegative_int, help="Maximum absolute reference gap/overlap. Default: 30 bp.")
    intron.add_argument("--max-query-overlap", default=20, type=nonnegative_int, help="Maximum allowed query HSP overlap. Default: 20 bp.")
    intron.add_argument("--breakpoint-window", default=30, type=nonnegative_int, help="Breakpoint clustering window. Default: 30 bp.")

    taxonomy = parser.add_argument_group("Taxonomy support")
    taxonomy.add_argument("--taxonomy", default=None, type=Path, help="Optional subject taxonomy TSV. If --ref-fasta is used, generated automatically.")
    taxonomy.add_argument("--tax-rank", default="genus", choices=["domain", "phylum", "class", "order", "family", "genus", "species"], help="Taxonomic rank used for confidence support. Default: genus.")

    confidence = parser.add_argument_group("Confidence thresholds")
    confidence.add_argument("--min-support-subjects", default=1, type=positive_int, help="Minimum subjects for LOW confidence. Default: 1.")
    confidence.add_argument("--medium-support-subjects", default=3, type=positive_int, help="Minimum subjects for MEDIUM confidence. Default: 3.")
    confidence.add_argument("--medium-support-taxa", default=3, type=positive_int, help="Minimum taxa for MEDIUM confidence. Default: 3.")
    confidence.add_argument("--high-support-subjects", default=10, type=positive_int, help="Minimum subjects for HIGH confidence. Default: 10.")
    confidence.add_argument("--high-support-taxa", default=3, type=positive_int, help="Minimum taxa for HIGH confidence. Default: 3.")
    confidence.add_argument("--min-output-confidence", default="LOW", choices=["LOW", "MEDIUM", "HIGH"], help="Minimum confidence removed in intron-free FASTA and written to intron FASTA/BED outputs. Default: LOW.")
    confidence.add_argument("--gzip-fasta-output", action="store_true", help="Write intron-free and intron FASTA outputs as .fa.gz. Default: disabled.")

    runtime = parser.add_argument_group("Runtime and logging")
    runtime.add_argument("--threads", default=4, type=positive_int, help="Worker threads. BLASTN also uses this value. Default: 4.")
    runtime.add_argument("--verbose", action="store_true", help="Print debug logs. Default: disabled.")

    return parser


def validate_args(args: argparse.Namespace) -> None:
    """Validate input arguments."""

    if not args.query.exists():
        raise FileNotFoundError(f"Query FASTA not found: {args.query}")
    supplied_modes = sum(x is not None for x in [args.blast, args.db, args.ref_fasta])
    if supplied_modes != 1:
        raise ValueError("Provide exactly one of --blast, --db, or --ref-fasta")
    if args.blast is not None and not args.blast.exists():
        raise FileNotFoundError(f"BLAST table not found: {args.blast}")
    if args.ref_fasta is not None and not args.ref_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {args.ref_fasta}")
    if args.taxonomy is not None and not args.taxonomy.exists():
        raise FileNotFoundError(f"Taxonomy table not found: {args.taxonomy}")
    if args.min_intron_len > args.max_intron_len:
        raise ValueError("--min-intron-len must be <= --max-intron-len")


def main() -> int:
    """CLI entry point."""

    parser = build_parser()
    args = parser.parse_args()
    setup_logging(args.verbose)
    validate_args(args)
    setup_output_paths(args)

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=CONSOLE,
    )

    with progress:
        task = progress.add_task("Reading query FASTA", total=None)
        seqs = read_fasta(args.query)
        progress.update(task, completed=1, total=1)
        LOG.info("Loaded query sequences: %s", len(seqs))

        if args.ref_fasta is not None:
            task = progress.add_task("Preprocessing reference and building BLAST DB", total=None)
            generated_ref, generated_taxonomy, generated_db = preprocess_reference(args)
            args.taxonomy = generated_taxonomy
            args.db = generated_db
            progress.update(task, completed=1, total=1)

            if args.clean_ref_introns:
                task = progress.add_task("Cleaning reference introns by self-BLAST", total=None)
                _, cleaned_taxonomy, cleaned_db = clean_reference_introns(args, generated_ref, generated_taxonomy, generated_db)
                args.taxonomy = cleaned_taxonomy
                args.db = cleaned_db
                progress.update(task, completed=1, total=1)

        task = progress.add_task("Reading taxonomy", total=None)
        taxonomy = parse_taxonomy(args.taxonomy)
        progress.update(task, completed=1, total=1)

        task = progress.add_task("Preparing BLAST HSPs", total=None)
        blast_path = args.blast if args.blast is not None else run_query_blastn(args)
        args.blast = blast_path
        blast_by_query = parse_blast(blast_path, args.min_pident, args.min_hsp_len)
        progress.update(task, completed=1, total=1)

        query_ids = sorted(seqs.keys())
        task = progress.add_task("Detecting introns", total=len(query_ids))

        def worker(query_id: str) -> QueryResult:
            return analyze_query(query_id, blast_by_query.get(query_id, {}), len(seqs.get(query_id, "")), taxonomy, args)

        results: List[QueryResult] = []
        if args.threads == 1:
            for query_id in query_ids:
                results.append(worker(query_id))
                progress.advance(task)
        else:
            with futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
                for result in executor.map(worker, query_ids, chunksize=128):
                    results.append(result)
                    progress.advance(task)

        results.sort(key=lambda r: (confidence_rank(r.confidence), r.support_subjects, r.support_taxa, r.query_id), reverse=True)

        task = progress.add_task("Writing outputs", total=6)
        write_summary(args.summary_tsv, results, args.tax_rank)
        progress.advance(task)
        write_supporting_hsps(args.supporting_tsv, results)
        progress.advance(task)
        write_fasta_outputs(args, results, seqs)
        progress.advance(task)
        write_bed_outputs(args, results)
        progress.advance(task)
        write_report(args.report_md, results, args)
        progress.advance(task)
        print_summary(results, args)
        progress.advance(task)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
