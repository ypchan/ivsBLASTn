from __future__ import annotations

import argparse
import concurrent.futures as futures
import shlex
import sys
from pathlib import Path
from typing import List, Optional, Sequence

from rich.progress import BarColumn, Progress, SpinnerColumn, TaskProgressColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn
from rich_argparse import RichHelpFormatter

from . import __version__
from .algorithm import analyze_query, confidence_rank
from .blast import parse_blast, run_query_blastn
from .fasta import read_fasta
from .logging import CONSOLE, LOG, setup_logging
from .merge import merge_chunk_outputs
from .models import QueryResult
from .outputs import print_summary, write_bed_outputs, write_fasta_outputs, write_report, write_summary, write_supporting_hsps
from .paths import setup_output_paths
from .reference import clean_reference_introns, preprocess_reference
from .slurm import render_slurm_array_script, submit_sbatch, write_slurm_array_script
from .split import split_fasta
from .taxonomy import parse_taxonomy


SUBCOMMANDS = {"run", "init-reference", "split", "submit-slurm", "merge"}
BLAST_DB_EXTENSIONS = (".nhr", ".nin", ".nsq", ".nal", ".ndb", ".njs", ".nog", ".nos", ".not", ".ntf", ".nto")
DEFAULT_REF_CLEAN_MIN_CONFIDENCE = "MEDIUM"
DEFAULT_MIN_PIDENT = 70.0
DEFAULT_MIN_HSP_LEN = 100
DEFAULT_TOP_SUBJECTS = 100
DEFAULT_REF_SELF_BLAST_MAX_TARGET_SEQS = 100
DEFAULT_REF_SELF_BLAST_MAX_HSPS = 20
DEFAULT_BLAST_MAX_HSPS = 5
DEFAULT_MIN_INTRON_LEN = 25
DEFAULT_MAX_INTRON_LEN = 2000
DEFAULT_MAX_REF_GAP = 15
DEFAULT_MAX_QUERY_OVERLAP = 20
DEFAULT_BREAKPOINT_WINDOW = 20
DEFAULT_MIN_OUTPUT_CONFIDENCE = "MEDIUM"


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


def blast_db_prefix_exists(prefix: Path) -> bool:
    """Return True when files for a nucleotide BLAST DB prefix are present."""

    return any(Path(str(prefix) + suffix).exists() for suffix in BLAST_DB_EXTENSIONS)


def quoted_option(option: str, value: object) -> str:
    """Render one CLI option/value pair for a shell script."""

    return f"{option} {shlex.quote(str(value))}"


def slurm_run_args(args: argparse.Namespace) -> List[str]:
    """Return detection arguments forwarded from submit-slurm to each run task."""

    forwarded = [
        ("--algorithm", args.algorithm),
        ("--min-pident", args.min_pident),
        ("--min-hsp-len", args.min_hsp_len),
        ("--min-intron-len", args.min_intron_len),
        ("--max-intron-len", args.max_intron_len),
        ("--max-ref-gap", args.max_ref_gap),
        ("--max-query-overlap", args.max_query_overlap),
        ("--breakpoint-window", args.breakpoint_window),
        ("--blastn-bin", args.blastn_bin),
        ("--blast-task", args.blast_task),
        ("--blast-evalue", args.blast_evalue),
        ("--tax-rank", args.tax_rank),
        ("--min-support-subjects", args.min_support_subjects),
        ("--medium-support-subjects", args.medium_support_subjects),
        ("--medium-support-taxa", args.medium_support_taxa),
        ("--high-support-subjects", args.high_support_subjects),
        ("--high-support-taxa", args.high_support_taxa),
        ("--min-output-confidence", args.min_output_confidence),
    ]
    rendered = [quoted_option(option, value) for option, value in forwarded]
    if args.gzip_fasta_output:
        rendered.append("--gzip-fasta-output")
    return rendered


def add_common_logging_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--verbose", action="store_true", help="Print debug logs. Default: disabled.")


def add_run_args(parser: argparse.ArgumentParser) -> None:
    required = parser.add_argument_group("Required inputs")
    required.add_argument("--query", required=True, type=Path, help="Query 16S FASTA/FASTA.gz file. No default.")
    required.add_argument("--blast", default=None, type=Path, help="Existing BLASTN outfmt 6 table. Provide exactly one of --blast, --db, or --ref-fasta.")
    required.add_argument("--db", default=None, type=Path, help="Existing BLAST database prefix. Provide exactly one of --blast, --db, or --ref-fasta.")
    required.add_argument("--ref-fasta", default=None, type=Path, help="SILVA-style reference FASTA/FASTA.gz with taxonomy in headers. Provide exactly one of --blast, --db, or --ref-fasta.")
    required.add_argument("--outdir", required=True, type=Path, help="Output directory. No default.")

    algorithm = parser.add_argument_group("Algorithm")
    algorithm.add_argument("--algorithm", default="hsp-gap-support", choices=["hsp-gap-support"], help="Detection algorithm. Default: hsp-gap-support.")

    filters = parser.add_argument_group("Reference preprocessing and HSP filters")
    filters.add_argument("--ref-domains", default="Archaea,Bacteria", help="Comma-separated SILVA domains retained with --ref-fasta. Default: Archaea,Bacteria.")
    filters.add_argument("--ref-per-species", default=1, type=nonnegative_int, help="Maximum sequences per clear species in reference preprocessing, keeping the longest sequences first. Use 0 to disable. Default: 1.")
    filters.add_argument("--ref-unclear-per-genus", default=5, type=nonnegative_int, help="Maximum unclear-species records retained per genus, keeping the longest sequences first. Use 0 to skip all unclear species. Default: 5.")
    filters.add_argument("--clean-ref-introns", action="store_true", help="Self-BLAST reference and remove candidate introns before query BLAST. Default: disabled.")
    filters.add_argument("--ref-clean-min-confidence", default=DEFAULT_REF_CLEAN_MIN_CONFIDENCE, choices=["LOW", "MEDIUM", "HIGH"], help=f"Minimum confidence required to remove a reference intron. Default: {DEFAULT_REF_CLEAN_MIN_CONFIDENCE}.")
    filters.add_argument("--ref-self-blast-max-target-seqs", default=DEFAULT_REF_SELF_BLAST_MAX_TARGET_SEQS, type=positive_int, help=f"Reference self-BLAST -max_target_seqs. Default: {DEFAULT_REF_SELF_BLAST_MAX_TARGET_SEQS}.")
    filters.add_argument("--ref-self-blast-max-hsps", default=DEFAULT_REF_SELF_BLAST_MAX_HSPS, type=positive_int, help=f"Reference self-BLAST -max_hsps. Default: {DEFAULT_REF_SELF_BLAST_MAX_HSPS}.")
    filters.add_argument("--min-pident", default=DEFAULT_MIN_PIDENT, type=probability_percent, help=f"Minimum HSP percent identity. Default: {DEFAULT_MIN_PIDENT}.")
    filters.add_argument("--min-hsp-len", default=DEFAULT_MIN_HSP_LEN, type=positive_int, help=f"Minimum HSP length in bp. Default: {DEFAULT_MIN_HSP_LEN}.")
    filters.add_argument("--top-subjects", default=DEFAULT_TOP_SUBJECTS, type=positive_int, help=f"Subjects requested from query BLASTN and retained per query after parsing. Default: {DEFAULT_TOP_SUBJECTS}.")
    filters.add_argument("--makeblastdb-bin", default="makeblastdb", help="makeblastdb executable. Default: makeblastdb.")
    filters.add_argument("--blastn-bin", default="blastn", help="blastn executable. Default: blastn.")
    filters.add_argument("--blast-max-hsps", default=DEFAULT_BLAST_MAX_HSPS, type=positive_int, help=f"Query BLASTN -max_hsps per subject. Default: {DEFAULT_BLAST_MAX_HSPS}.")
    filters.add_argument("--blast-task", default="blastn", choices=["blastn", "megablast", "dc-megablast", "blastn-short"], help="BLASTN task. Default: blastn.")
    filters.add_argument("--blast-evalue", default="1e-20", help="BLASTN e-value. Default: 1e-20.")

    intron = parser.add_argument_group("Candidate IVS geometry")
    intron.add_argument("--min-intron-len", default=DEFAULT_MIN_INTRON_LEN, type=nonnegative_int, help=f"Minimum query gap size. Default: {DEFAULT_MIN_INTRON_LEN} bp.")
    intron.add_argument("--max-intron-len", default=DEFAULT_MAX_INTRON_LEN, type=positive_int, help=f"Maximum query gap size. Default: {DEFAULT_MAX_INTRON_LEN} bp.")
    intron.add_argument("--max-ref-gap", default=DEFAULT_MAX_REF_GAP, type=nonnegative_int, help=f"Maximum absolute reference gap/overlap. Default: {DEFAULT_MAX_REF_GAP} bp.")
    intron.add_argument("--max-query-overlap", default=DEFAULT_MAX_QUERY_OVERLAP, type=nonnegative_int, help=f"Maximum allowed query HSP overlap. Default: {DEFAULT_MAX_QUERY_OVERLAP} bp.")
    intron.add_argument("--breakpoint-window", default=DEFAULT_BREAKPOINT_WINDOW, type=nonnegative_int, help=f"Breakpoint clustering window. Default: {DEFAULT_BREAKPOINT_WINDOW} bp.")

    taxonomy = parser.add_argument_group("Taxonomy support")
    taxonomy.add_argument("--taxonomy", default=None, type=Path, help="Optional subject taxonomy TSV. If --ref-fasta is used, generated automatically.")
    taxonomy.add_argument("--tax-rank", default="genus", choices=["domain", "phylum", "class", "order", "family", "genus", "species"], help="Taxonomic rank used for confidence support. Default: genus.")

    confidence = parser.add_argument_group("Confidence thresholds")
    confidence.add_argument("--min-support-subjects", default=1, type=positive_int, help="Minimum subjects for LOW confidence. Default: 1.")
    confidence.add_argument("--medium-support-subjects", default=3, type=positive_int, help="Minimum subjects for MEDIUM confidence. Default: 3.")
    confidence.add_argument("--medium-support-taxa", default=3, type=positive_int, help="Minimum taxa for MEDIUM confidence. Default: 3.")
    confidence.add_argument("--high-support-subjects", default=10, type=positive_int, help="Minimum subjects for HIGH confidence. Default: 10.")
    confidence.add_argument("--high-support-taxa", default=3, type=positive_int, help="Minimum taxa for HIGH confidence. Default: 3.")
    confidence.add_argument("--min-output-confidence", default=DEFAULT_MIN_OUTPUT_CONFIDENCE, choices=["LOW", "MEDIUM", "HIGH"], help=f"Minimum confidence removed in intron-free FASTA and written to intron FASTA/BED outputs. Default: {DEFAULT_MIN_OUTPUT_CONFIDENCE}.")
    confidence.add_argument("--gzip-fasta-output", action="store_true", help="Write intron-free and intron FASTA outputs as .fa.gz. Default: disabled.")

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("--threads", default=4, type=positive_int, help="Worker threads. BLASTN also uses this value. Default: 4.")
    add_common_logging_args(parser)


def build_parser() -> argparse.ArgumentParser:
    """Create argument parser."""

    parser = argparse.ArgumentParser(
        prog="ivsBLASTn",
        description="Detect 16S/SSU intervening sequences by BLASTN HSP geometry.",
        formatter_class=RichHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"ivsBLASTn {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    run_parser = subparsers.add_parser("run", help="Run IVS detection on one query FASTA or one chunk.", formatter_class=RichHelpFormatter)
    add_run_args(run_parser)
    run_parser.set_defaults(func=run_command)

    init_ref_parser = subparsers.add_parser("init-reference", help="Prepare a reusable reference FASTA, taxonomy TSV, and BLAST DB.", formatter_class=RichHelpFormatter)
    init_inputs = init_ref_parser.add_argument_group("Required inputs and outputs")
    init_inputs.add_argument("--ref-fasta", required=True, type=Path, help="SILVA-style reference FASTA/FASTA.gz with taxonomy in headers.")
    init_inputs.add_argument("--outdir", required=True, type=Path, help="Reference initialization output directory.")

    init_selection = init_ref_parser.add_argument_group("Reference selection")
    init_selection.add_argument("--ref-domains", default="Archaea,Bacteria", help="Comma-separated SILVA domains retained. Default: Archaea,Bacteria.")
    init_selection.add_argument("--ref-per-species", default=1, type=nonnegative_int, help="Maximum sequences per clear species, longest first. Use 0 to disable. Default: 1.")
    init_selection.add_argument("--ref-unclear-per-genus", default=5, type=nonnegative_int, help="Maximum unclear-species records retained per genus, longest first. Use 0 to skip all unclear species. Default: 5.")

    init_self_clean = init_ref_parser.add_argument_group("Optional reference self-cleaning")
    init_self_clean.add_argument("--clean-ref-introns", action="store_true", help="Self-BLAST reference and remove candidate IVSs before final DB creation.")
    init_self_clean.add_argument("--ref-clean-min-confidence", default=DEFAULT_REF_CLEAN_MIN_CONFIDENCE, choices=["LOW", "MEDIUM", "HIGH"], help=f"Minimum confidence required to remove a reference IVS. Default: {DEFAULT_REF_CLEAN_MIN_CONFIDENCE}.")
    init_self_clean.add_argument("--ref-self-blast-max-target-seqs", default=DEFAULT_REF_SELF_BLAST_MAX_TARGET_SEQS, type=positive_int, help=f"Reference self-BLAST -max_target_seqs. Default: {DEFAULT_REF_SELF_BLAST_MAX_TARGET_SEQS}.")
    init_self_clean.add_argument("--ref-self-blast-max-hsps", default=DEFAULT_REF_SELF_BLAST_MAX_HSPS, type=positive_int, help=f"Reference self-BLAST -max_hsps. Default: {DEFAULT_REF_SELF_BLAST_MAX_HSPS}.")
    init_self_clean.add_argument("--min-pident", default=DEFAULT_MIN_PIDENT, type=probability_percent, help=f"Minimum HSP percent identity for optional self-cleaning. Default: {DEFAULT_MIN_PIDENT}.")
    init_self_clean.add_argument("--min-hsp-len", default=DEFAULT_MIN_HSP_LEN, type=positive_int, help=f"Minimum HSP length for optional self-cleaning. Default: {DEFAULT_MIN_HSP_LEN}.")
    init_self_clean.add_argument("--top-subjects", default=DEFAULT_TOP_SUBJECTS, type=positive_int, help=f"Top subjects retained per reference during optional self-cleaning. Default: {DEFAULT_TOP_SUBJECTS}.")
    init_self_clean.add_argument("--algorithm", default="hsp-gap-support", choices=["hsp-gap-support"], help="Detection algorithm for optional self-cleaning. Default: hsp-gap-support.")
    init_self_clean.add_argument("--tax-rank", default="genus", choices=["domain", "phylum", "class", "order", "family", "genus", "species"], help="Taxonomic rank used during optional self-cleaning. Default: genus.")

    init_geometry = init_ref_parser.add_argument_group("Self-cleaning IVS geometry")
    init_geometry.add_argument("--min-intron-len", default=DEFAULT_MIN_INTRON_LEN, type=nonnegative_int, help=f"Minimum query gap size for optional self-cleaning. Default: {DEFAULT_MIN_INTRON_LEN} bp.")
    init_geometry.add_argument("--max-intron-len", default=DEFAULT_MAX_INTRON_LEN, type=positive_int, help=f"Maximum query gap size for optional self-cleaning. Default: {DEFAULT_MAX_INTRON_LEN} bp.")
    init_geometry.add_argument("--max-ref-gap", default=DEFAULT_MAX_REF_GAP, type=nonnegative_int, help=f"Maximum absolute reference gap/overlap for optional self-cleaning. Default: {DEFAULT_MAX_REF_GAP} bp.")
    init_geometry.add_argument("--max-query-overlap", default=DEFAULT_MAX_QUERY_OVERLAP, type=nonnegative_int, help=f"Maximum allowed query HSP overlap. Default: {DEFAULT_MAX_QUERY_OVERLAP} bp.")
    init_geometry.add_argument("--breakpoint-window", default=DEFAULT_BREAKPOINT_WINDOW, type=nonnegative_int, help=f"Breakpoint clustering window. Default: {DEFAULT_BREAKPOINT_WINDOW} bp.")

    init_confidence = init_ref_parser.add_argument_group("Self-cleaning confidence thresholds")
    init_confidence.add_argument("--min-support-subjects", default=1, type=positive_int, help="Minimum subjects for LOW confidence. Default: 1.")
    init_confidence.add_argument("--medium-support-subjects", default=3, type=positive_int, help="Minimum subjects for MEDIUM confidence. Default: 3.")
    init_confidence.add_argument("--medium-support-taxa", default=3, type=positive_int, help="Minimum taxa for MEDIUM confidence. Default: 3.")
    init_confidence.add_argument("--high-support-subjects", default=10, type=positive_int, help="Minimum subjects for HIGH confidence. Default: 10.")
    init_confidence.add_argument("--high-support-taxa", default=3, type=positive_int, help="Minimum taxa for HIGH confidence. Default: 3.")
    init_confidence.add_argument("--min-output-confidence", default=DEFAULT_MIN_OUTPUT_CONFIDENCE, choices=["LOW", "MEDIUM", "HIGH"], help=argparse.SUPPRESS)

    init_runtime = init_ref_parser.add_argument_group("Runtime and external tools")
    init_runtime.add_argument("--threads", default=4, type=positive_int, help="Worker and BLASTN threads for optional self-cleaning. Default: 4.")
    init_runtime.add_argument("--makeblastdb-bin", default="makeblastdb", help="makeblastdb executable. Default: makeblastdb.")
    init_runtime.add_argument("--blastn-bin", default="blastn", help="blastn executable. Default: blastn.")
    init_runtime.add_argument("--blast-task", default="blastn", choices=["blastn", "megablast", "dc-megablast", "blastn-short"], help="BLASTN task for optional self-cleaning. Default: blastn.")
    init_runtime.add_argument("--blast-evalue", default="1e-20", help="BLASTN e-value for optional self-cleaning. Default: 1e-20.")
    add_common_logging_args(init_ref_parser)
    init_ref_parser.set_defaults(func=init_reference_command)

    split_parser = subparsers.add_parser("split", help="Split a large query FASTA into chunk FASTA files.", formatter_class=RichHelpFormatter)
    split_inputs = split_parser.add_argument_group("Inputs and outputs")
    split_inputs.add_argument("--query", required=True, type=Path, help="Input query FASTA/FASTA.gz.")
    split_inputs.add_argument("--chunks-dir", required=True, type=Path, help="Output directory for chunk FASTA files.")
    split_settings = split_parser.add_argument_group("Chunk settings")
    split_settings.add_argument("--chunk-size", default=5000, type=positive_int, help="Records per chunk. Default: 5000.")
    split_settings.add_argument("--chunk-prefix", default="query", help="Chunk filename prefix. Default: query.")
    add_common_logging_args(split_parser)
    split_parser.set_defaults(func=split_command)

    submit_parser = subparsers.add_parser("submit-slurm", help="Generate and optionally submit a Slurm array over query chunks.", formatter_class=RichHelpFormatter)
    submit_inputs = submit_parser.add_argument_group("Required inputs and outputs")
    submit_inputs.add_argument("--chunks-dir", required=True, type=Path, help="Directory created by `ivsBLASTn split`.")
    submit_inputs.add_argument("--db", required=True, type=Path, help="BLAST database prefix shared by all chunks.")
    submit_inputs.add_argument("--taxonomy", default=None, type=Path, help="Optional taxonomy TSV shared by all chunks.")
    submit_inputs.add_argument("--outdir", required=True, type=Path, help="Batch run output directory.")
    submit_inputs.add_argument("--chunk-runs-dir", default=None, type=Path, help="Directory for per-chunk run output directories. Default: OUTDIR.")

    submit_detection = submit_parser.add_argument_group("Detection parameters forwarded to each chunk")
    submit_detection.add_argument("--threads", default=8, type=positive_int, help="Threads passed to each ivsBLASTn run. Default: 8.")
    submit_detection.add_argument("--top-subjects", default=DEFAULT_TOP_SUBJECTS, type=positive_int, help=f"Subjects requested and analyzed per query. Default: {DEFAULT_TOP_SUBJECTS}.")
    submit_detection.add_argument("--blast-max-hsps", default=DEFAULT_BLAST_MAX_HSPS, type=positive_int, help=f"HSPs requested per query-subject pair. Default: {DEFAULT_BLAST_MAX_HSPS}.")
    submit_detection.add_argument("--algorithm", default="hsp-gap-support", choices=["hsp-gap-support"], help="Detection algorithm forwarded to each run. Default: hsp-gap-support.")
    submit_detection.add_argument("--min-pident", default=DEFAULT_MIN_PIDENT, type=probability_percent, help=f"Minimum HSP percent identity forwarded to each run. Default: {DEFAULT_MIN_PIDENT}.")
    submit_detection.add_argument("--min-hsp-len", default=DEFAULT_MIN_HSP_LEN, type=positive_int, help=f"Minimum HSP length forwarded to each run. Default: {DEFAULT_MIN_HSP_LEN}.")
    submit_detection.add_argument("--min-intron-len", default=DEFAULT_MIN_INTRON_LEN, type=nonnegative_int, help=f"Minimum query gap size forwarded to each run. Default: {DEFAULT_MIN_INTRON_LEN} bp.")
    submit_detection.add_argument("--max-intron-len", default=DEFAULT_MAX_INTRON_LEN, type=positive_int, help=f"Maximum query gap size forwarded to each run. Default: {DEFAULT_MAX_INTRON_LEN} bp.")
    submit_detection.add_argument("--max-ref-gap", default=DEFAULT_MAX_REF_GAP, type=nonnegative_int, help=f"Maximum absolute reference gap/overlap forwarded to each run. Default: {DEFAULT_MAX_REF_GAP} bp.")
    submit_detection.add_argument("--max-query-overlap", default=DEFAULT_MAX_QUERY_OVERLAP, type=nonnegative_int, help=f"Maximum allowed query HSP overlap forwarded to each run. Default: {DEFAULT_MAX_QUERY_OVERLAP} bp.")
    submit_detection.add_argument("--breakpoint-window", default=DEFAULT_BREAKPOINT_WINDOW, type=nonnegative_int, help=f"Breakpoint clustering window forwarded to each run. Default: {DEFAULT_BREAKPOINT_WINDOW} bp.")
    submit_detection.add_argument("--blastn-bin", default="blastn", help="blastn executable forwarded to each run. Default: blastn.")
    submit_detection.add_argument("--blast-task", default="blastn", choices=["blastn", "megablast", "dc-megablast", "blastn-short"], help="BLASTN task forwarded to each run. Default: blastn.")
    submit_detection.add_argument("--blast-evalue", default="1e-20", help="BLASTN e-value forwarded to each run. Default: 1e-20.")
    submit_detection.add_argument("--tax-rank", default="genus", choices=["domain", "phylum", "class", "order", "family", "genus", "species"], help="Taxonomic rank forwarded to each run. Default: genus.")
    submit_detection.add_argument("--min-support-subjects", default=1, type=positive_int, help="Minimum subjects for LOW confidence forwarded to each run. Default: 1.")
    submit_detection.add_argument("--medium-support-subjects", default=3, type=positive_int, help="Minimum subjects for MEDIUM confidence forwarded to each run. Default: 3.")
    submit_detection.add_argument("--medium-support-taxa", default=3, type=positive_int, help="Minimum taxa for MEDIUM confidence forwarded to each run. Default: 3.")
    submit_detection.add_argument("--high-support-subjects", default=10, type=positive_int, help="Minimum subjects for HIGH confidence forwarded to each run. Default: 10.")
    submit_detection.add_argument("--high-support-taxa", default=3, type=positive_int, help="Minimum taxa for HIGH confidence forwarded to each run. Default: 3.")
    submit_detection.add_argument("--min-output-confidence", default=DEFAULT_MIN_OUTPUT_CONFIDENCE, choices=["LOW", "MEDIUM", "HIGH"], help=f"Minimum confidence for FASTA/BED outputs forwarded to each run. Default: {DEFAULT_MIN_OUTPUT_CONFIDENCE}.")
    submit_detection.add_argument("--gzip-fasta-output", action="store_true", help="Write per-chunk FASTA outputs as .fa.gz and merge them as gzip streams. Default: disabled.")

    submit_slurm = submit_parser.add_argument_group("Slurm resources and platform directives")
    submit_slurm.add_argument("--cpus-per-task", default=8, type=positive_int, help="Slurm CPUs per array task. Default: 8.")
    submit_slurm.add_argument("--mem", default="16G", help="Slurm memory per task. Default: 16G.")
    submit_slurm.add_argument("--time", default="12:00:00", help="Slurm time limit. Default: 12:00:00.")
    submit_slurm.add_argument("--partition", default=None, help="Slurm partition/queue, e.g. normal_fcp1.")
    submit_slurm.add_argument("--account", default=None, help="Slurm account/project, e.g. prj_219_3. Writes #SBATCH --account.")
    submit_slurm.add_argument("--qos", default=None, help="Slurm QoS, e.g. qos_prj_219_3. Writes #SBATCH --qos.")
    submit_slurm.add_argument("--nodes", default=None, type=positive_int, help="Optional Slurm node count, written as #SBATCH --nodes.")
    submit_slurm.add_argument("--ntasks", default=None, type=positive_int, help="Optional Slurm task count, written as #SBATCH --ntasks.")
    submit_slurm.add_argument("--constraint", default=None, help="Optional Slurm node constraint, written as #SBATCH --constraint.")
    submit_slurm.add_argument("--gres", default=None, help="Optional Slurm generic resources, written as #SBATCH --gres.")
    submit_slurm.add_argument("--exclude", default=None, help="Optional Slurm excluded node list, written as #SBATCH --exclude.")
    submit_slurm.add_argument("--nodelist", default=None, help="Optional Slurm node list, written as #SBATCH --nodelist.")
    submit_slurm.add_argument("--sbatch-option", action="append", default=[], help="Extra SBATCH directive, for example '--mail-type=END'. Can be repeated.")

    submit_control = submit_parser.add_argument_group("Submission control")
    submit_control.add_argument("--array-concurrency", default=None, type=positive_int, help="Optional Slurm array concurrency limit.")
    submit_control.add_argument("--extra-run-args", default="", help="Extra arguments appended to each `ivsBLASTn run` command.")
    submit_control.add_argument("--submit", action="store_true", help="Submit with sbatch after writing the script. Default: write only.")
    add_common_logging_args(submit_parser)
    submit_parser.set_defaults(func=submit_slurm_command)

    merge_parser = subparsers.add_parser("merge", help="Merge per-chunk ivsBLASTn outputs.", formatter_class=RichHelpFormatter)
    merge_inputs = merge_parser.add_argument_group("Inputs and outputs")
    merge_inputs.add_argument("--chunk-results-dir", required=True, type=Path, help="Directory containing per-chunk run output directories, usually the submit-slurm OUTDIR.")
    merge_inputs.add_argument("--outdir", required=True, type=Path, help="Final merged output directory.")
    merge_settings = merge_parser.add_argument_group("Merge settings")
    merge_settings.add_argument("--label", default="merged", help="Output file prefix. Default: merged.")
    add_common_logging_args(merge_parser)
    merge_parser.set_defaults(func=merge_command)

    return parser


def normalize_legacy_argv(argv: Sequence[str]) -> List[str]:
    """Map old `ivsBLASTn --query ...` usage to `ivsBLASTn run --query ...`."""

    args = list(argv)
    if not args:
        return args
    if args[0] in SUBCOMMANDS or args[0] in {"-h", "--help", "--version"}:
        return args
    if args[0].startswith("-"):
        return ["run", *args]
    return args


def validate_run_args(args: argparse.Namespace) -> None:
    """Validate run command arguments."""

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


def run_pipeline(args: argparse.Namespace) -> int:
    """Run the original single-query-FASTA ivsBLASTn workflow."""

    validate_run_args(args)
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

        query_ids = list(seqs.keys())
        task = progress.add_task("Detecting IVSs", total=len(query_ids))

        def worker(query_id: str) -> QueryResult:
            return analyze_query(query_id, blast_by_query.get(query_id, {}), len(seqs.get(query_id, "")), taxonomy, args)

        results: List[QueryResult] = []
        if args.threads == 1:
            for query_id in query_ids:
                results.append(worker(query_id))
                progress.advance(task)
        else:
            with futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
                for result in executor.map(worker, query_ids, chunksize=256):
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


def run_command(args: argparse.Namespace) -> int:
    return run_pipeline(args)


def setup_reference_output_paths(args: argparse.Namespace) -> None:
    """Populate paths expected by reference preprocessing for init-reference."""

    args.reference_dir = args.outdir
    args.blast_dir = args.outdir / "blast"
    args.results_dir = args.outdir / "results"
    for directory in [args.outdir, args.blast_dir, args.results_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    args.raw_ref_fa = args.outdir / "raw_reference.fa"
    args.raw_ref_tax = args.outdir / "raw_reference.tax.tsv"
    args.raw_ref_db = args.outdir / "raw_reference_db"
    args.ref_self_blast = args.blast_dir / "reference_self.blastn.tsv"
    args.ref_self_clean_prefix = args.results_dir / "reference_self_clean"
    args.ref_self_clean_introns_fa = args.results_dir / "reference_self_clean.introns.fa"
    args.cleaned_ref_fa = args.outdir / "cleaned_reference.fa"
    args.cleaned_ref_tax = args.outdir / "cleaned_reference.tax.tsv"
    args.cleaned_ref_db = args.outdir / "cleaned_reference_db"
    args.query = args.ref_fasta
    args.blast = None
    args.db = None
    args.taxonomy = None
    args.blast_max_hsps = getattr(args, "ref_self_blast_max_hsps", 20)
    args.min_output_confidence = getattr(args, "ref_clean_min_confidence", DEFAULT_REF_CLEAN_MIN_CONFIDENCE)


def init_reference_command(args: argparse.Namespace) -> int:
    if not args.ref_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {args.ref_fasta}")
    if args.min_intron_len > args.max_intron_len:
        raise ValueError("--min-intron-len must be <= --max-intron-len")
    setup_reference_output_paths(args)
    ref_fa, taxonomy_tsv, db_prefix = preprocess_reference(args)
    if args.clean_ref_introns:
        ref_fa, taxonomy_tsv, db_prefix = clean_reference_introns(args, ref_fa, taxonomy_tsv, db_prefix)

    manifest = args.outdir / "reference_manifest.tsv"
    with manifest.open("wt", encoding="utf-8") as handle:
        print("key\tpath", file=handle)
        print(f"reference_fasta\t{ref_fa.resolve()}", file=handle)
        print(f"taxonomy_tsv\t{taxonomy_tsv.resolve()}", file=handle)
        print(f"blast_db_prefix\t{db_prefix.resolve()}", file=handle)
        if args.clean_ref_introns:
            print(f"reference_self_clean_report\t{(args.results_dir / 'reference_self_clean.report.md').resolve()}", file=handle)
            print(f"reference_self_clean_summary\t{(args.results_dir / 'reference_self_clean.summary.tsv').resolve()}", file=handle)
            print(f"reference_introns_fasta\t{args.ref_self_clean_introns_fa.resolve()}", file=handle)

    CONSOLE.print("Reference initialized")
    CONSOLE.print(f"  DB prefix: {db_prefix}")
    CONSOLE.print(f"  Taxonomy:  {taxonomy_tsv}")
    if args.clean_ref_introns:
        CONSOLE.print(f"  Self-clean report: {args.results_dir / 'reference_self_clean.report.md'}")
        CONSOLE.print(f"  Reference IVSs:    {args.ref_self_clean_introns_fa}")
    CONSOLE.print("Use with:")
    CONSOLE.print(f"  ivsBLASTn run --query query.fa --db {db_prefix} --taxonomy {taxonomy_tsv} --outdir ivs_run")
    return 0


def split_command(args: argparse.Namespace) -> int:
    if not args.query.exists():
        raise FileNotFoundError(f"Query FASTA not found: {args.query}")
    split_fasta(args.query, args.chunks_dir, args.chunk_size, args.chunk_prefix)
    return 0


def submit_slurm_command(args: argparse.Namespace) -> int:
    manifest = args.chunks_dir / "chunks.tsv"
    if not manifest.exists():
        raise FileNotFoundError(f"Chunk manifest not found: {manifest}")
    if args.min_intron_len > args.max_intron_len:
        raise ValueError("--min-intron-len must be <= --max-intron-len")
    if not blast_db_prefix_exists(args.db):
        LOG.warning("No BLAST DB files found for prefix on this filesystem: %s", args.db)
    if not args.taxonomy and "--taxonomy" not in args.extra_run_args:
        LOG.warning("No taxonomy TSV supplied; MEDIUM/HIGH confidence will be harder to reach")
    script = render_slurm_array_script(
        chunks_dir=args.chunks_dir,
        outdir=args.outdir,
        chunk_runs_dir=args.chunk_runs_dir,
        db=args.db,
        taxonomy=args.taxonomy,
        threads=args.threads,
        top_subjects=args.top_subjects,
        blast_max_hsps=args.blast_max_hsps,
        cpus_per_task=args.cpus_per_task,
        mem=args.mem,
        time=args.time,
        partition=args.partition,
        account=args.account,
        qos=args.qos,
        nodes=args.nodes,
        ntasks=args.ntasks,
        constraint=args.constraint,
        gres=args.gres,
        exclude=args.exclude,
        nodelist=args.nodelist,
        extra_sbatch_options=args.sbatch_option,
        array_concurrency=args.array_concurrency,
        extra_run_args=args.extra_run_args,
        run_args=slurm_run_args(args),
    )
    script_path = write_slurm_array_script(script, args.outdir)
    if args.submit:
        job_id = submit_sbatch(script_path)
        (args.outdir / "slurm" / "job_id.txt").write_text(f"{job_id}\n", encoding="utf-8")
    else:
        CONSOLE.print(f"Slurm script written: {script_path}")
        CONSOLE.print(f"Chunk run outputs: {args.chunk_runs_dir or args.outdir}")
        CONSOLE.print("Review it, then submit with:")
        CONSOLE.print(f"  sbatch {shlex.quote(str(script_path))}")
    return 0


def merge_command(args: argparse.Namespace) -> int:
    merge_chunk_outputs(args.chunk_results_dir, args.outdir, args.label)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry point."""

    parser = build_parser()
    parsed_argv = normalize_legacy_argv(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(parsed_argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    setup_logging(getattr(args, "verbose", False))
    try:
        return args.func(args)
    except Exception:
        if getattr(args, "verbose", False):
            LOG.exception("ivsBLASTn failed")
        else:
            LOG.error("ivsBLASTn failed: %s", sys.exc_info()[1])
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
