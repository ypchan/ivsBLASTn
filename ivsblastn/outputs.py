from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from statistics import median
from typing import Dict, List

from rich.table import Table

from .algorithm import confidence_rank, is_intron_result
from .fasta import open_text_auto_write, write_fasta_record
from .logging import CONSOLE
from .models import QueryResult, SupportPair

def summary_row(result: QueryResult, tax_rank: str) -> Dict[str, str]:
    """Convert one query result to a TSV row."""

    return {
        "query_id": result.query_id,
        "query_len": str(result.query_len),
        "classification": result.classification,
        "confidence": result.confidence,
        "intron_start": str(result.intron_start) if result.intron_start else "",
        "intron_end": str(result.intron_end) if result.intron_end else "",
        "intron_len": str(result.intron_len) if result.intron_len else "",
        "exon1": f"{result.exon1_start}-{result.exon1_end}" if result.exon1_start else "",
        "exon2": f"{result.exon2_start}-{result.exon2_end}" if result.exon2_start else "",
        "intron_free_len": str(result.intron_free_len) if result.intron_free_len else "",
        "support_subjects": str(result.support_subjects),
        f"support_taxa_at_{tax_rank}": str(result.support_taxa),
        "support_species": str(result.support_species),
        "support_genera": str(result.support_genera),
        "median_subject_gap": f"{result.median_subject_gap:.2f}",
        "median_pident": f"{result.median_pident:.2f}",
        "mean_bitscore": f"{result.mean_bitscore:.2f}",
        "best_subject": result.best_subject,
        "best_subject_taxonomy": result.best_subject_taxonomy,
        "reasons": "|".join(result.reasons or []),
    }


def support_row(pair: SupportPair) -> Dict[str, str]:
    """Convert one support pair to a TSV row."""

    h1, h2 = sorted([pair.hsp1, pair.hsp2], key=lambda h: (h.qlo, h.qhi))
    return {
        "query_id": pair.query_id,
        "subject_id": pair.subject_id,
        "taxonomy": pair.taxonomy,
        "taxon_at_rank": pair.taxon_at_rank,
        "intron_start": str(pair.intron_start),
        "intron_end": str(pair.intron_end),
        "intron_len": str(pair.intron_len),
        "query_gap": str(pair.query_gap),
        "subject_gap": str(pair.subject_gap),
        "hsp1_q": f"{h1.qlo}-{h1.qhi}",
        "hsp2_q": f"{h2.qlo}-{h2.qhi}",
        "hsp1_s": f"{h1.slo}-{h1.shi}",
        "hsp2_s": f"{h2.slo}-{h2.shi}",
        "hsp1_pident": f"{h1.pident:.2f}",
        "hsp2_pident": f"{h2.pident:.2f}",
        "hsp1_len": str(h1.length),
        "hsp2_len": str(h2.length),
        "bitscore_sum": f"{h1.bitscore + h2.bitscore:.2f}",
        "pair_score": f"{pair.pair_score:.2f}",
    }


def write_summary(path: Path, results: List[QueryResult], tax_rank: str) -> None:
    """Write query-level summary TSV."""

    fields = [
        "query_id",
        "query_len",
        "classification",
        "confidence",
        "intron_start",
        "intron_end",
        "intron_len",
        "exon1",
        "exon2",
        "intron_free_len",
        "support_subjects",
        f"support_taxa_at_{tax_rank}",
        "support_species",
        "support_genera",
        "median_subject_gap",
        "median_pident",
        "mean_bitscore",
        "best_subject",
        "best_subject_taxonomy",
        "reasons",
    ]
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=chr(9), lineterminator=chr(10))
        writer.writeheader()
        for result in results:
            writer.writerow(summary_row(result, tax_rank))


def write_supporting_hsps(path: Path, results: List[QueryResult]) -> None:
    """Write subject-level supporting HSP pairs."""

    fields = [
        "query_id",
        "subject_id",
        "taxonomy",
        "taxon_at_rank",
        "intron_start",
        "intron_end",
        "intron_len",
        "query_gap",
        "subject_gap",
        "hsp1_q",
        "hsp2_q",
        "hsp1_s",
        "hsp2_s",
        "hsp1_pident",
        "hsp2_pident",
        "hsp1_len",
        "hsp2_len",
        "bitscore_sum",
        "pair_score",
    ]
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=chr(9), lineterminator=chr(10))
        writer.writeheader()
        for result in results:
            for pair in result.support_pairs or []:
                writer.writerow(support_row(pair))


def write_fasta_outputs(args: argparse.Namespace, results: List[QueryResult], seqs: Dict[str, str]) -> None:
    """Write intron-free 16S and intron FASTA outputs."""

    result_by_query = {result.query_id: result for result in results}
    with open_text_auto_write(args.intron_free_fa) as free_out, open_text_auto_write(args.introns_fa) as intron_out:
        for query_id, seq in seqs.items():
            if not seq:
                continue
            result = result_by_query.get(query_id)
            if result is None:
                write_fasta_record(free_out, f"{query_id}|intron_free|intron=none|confidence=NONE|action=unchanged", seq)
                continue
            if is_intron_result(result, args.min_output_confidence):
                exon1 = seq[result.exon1_start - 1 : result.exon1_end]
                exon2 = seq[result.exon2_start - 1 : result.exon2_end]
                intron = seq[result.intron_start - 1 : result.intron_end]
                write_fasta_record(free_out, f"{result.query_id}|intron_free|intron={result.intron_start}-{result.intron_end}|confidence={result.confidence}|action=removed", exon1 + exon2)
                write_fasta_record(intron_out, f"{result.query_id}|intron|{result.intron_start}-{result.intron_end}|len={result.intron_len}|confidence={result.confidence}", intron)
            elif result.intron_start and result.intron_end:
                write_fasta_record(free_out, f"{result.query_id}|intron_free|intron={result.intron_start}-{result.intron_end}|confidence={result.confidence}|action=unchanged_below_min_output_confidence", seq)
            else:
                write_fasta_record(free_out, f"{result.query_id}|intron_free|intron=none|confidence={result.confidence}|action=unchanged", seq)


def write_bed_outputs(args: argparse.Namespace, results: List[QueryResult]) -> None:
    """Write BED files."""

    with args.introns_bed.open("wt", encoding="utf-8", newline="") as intron_file, args.exons_bed.open("wt", encoding="utf-8", newline="") as exon_file:
        intron_writer = csv.writer(intron_file, delimiter=chr(9), lineterminator=chr(10))
        exon_writer = csv.writer(exon_file, delimiter=chr(9), lineterminator=chr(10))
        for result in results:
            if not is_intron_result(result, args.min_output_confidence):
                continue
            intron_writer.writerow([result.query_id, result.intron_start - 1, result.intron_end, f"{result.query_id}|intron|{result.intron_start}-{result.intron_end}|confidence={result.confidence}", ".", "+"])
            exon_writer.writerow([result.query_id, result.exon1_start - 1, result.exon1_end, f"{result.query_id}|exon1|{result.exon1_start}-{result.exon1_end}|confidence={result.confidence}", ".", "+"])
            exon_writer.writerow([result.query_id, result.exon2_start - 1, result.exon2_end, f"{result.query_id}|exon2|{result.exon2_start}-{result.exon2_end}|confidence={result.confidence}", ".", "+"])


def write_report(path: Path, results: List[QueryResult], args: argparse.Namespace) -> None:
    """Write Markdown report."""

    counts = Counter(r.confidence for r in results)
    class_counts = Counter(r.classification for r in results)
    intron_lens = [r.intron_len for r in results if r.intron_len > 0]

    with path.open("wt", encoding="utf-8") as handle:
        print("# ivsBLASTn report", file=handle)
        print("", file=handle)
        print("## Inputs", file=handle)
        print("", file=handle)
        print(f"- Query FASTA: `{args.query}`", file=handle)
        print(f"- BLAST table: `{args.blast}`", file=handle)
        print(f"- BLAST database: `{args.db}`", file=handle)
        print(f"- Reference FASTA: `{args.ref_fasta}`", file=handle)
        print(f"- Taxonomy table: `{args.taxonomy}`", file=handle)
        print(f"- Output directory: `{args.outdir}`", file=handle)
        print("", file=handle)
        print("## Algorithm", file=handle)
        print("", file=handle)
        print(f"- Algorithm: `{args.algorithm}`", file=handle)
        print("- Detects HSP1--query insertion--HSP2 patterns where reference subject coordinates are nearly continuous.", file=handle)
        print("", file=handle)
        print("## Parameters", file=handle)
        print("", file=handle)
        for name in [
            "min_pident",
            "min_hsp_len",
            "min_intron_len",
            "max_intron_len",
            "max_ref_gap",
            "max_query_overlap",
            "breakpoint_window",
            "top_subjects",
            "ref_domains",
            "ref_per_species",
            "ref_unclear_per_genus",
            "clean_ref_introns",
            "ref_clean_min_confidence",
            "ref_self_blast_max_target_seqs",
            "ref_self_blast_max_hsps",
            "blast_max_hsps",
            "blast_task",
            "blast_evalue",
            "tax_rank",
            "min_support_subjects",
            "medium_support_subjects",
            "medium_support_taxa",
            "high_support_subjects",
            "high_support_taxa",
            "min_output_confidence",
            "threads",
        ]:
            print(f"- `{name}`: `{getattr(args, name, 'NA')}`", file=handle)
        print("", file=handle)
        print("## Summary", file=handle)
        print("", file=handle)
        print(f"- Total query results: `{len(results)}`", file=handle)
        for label in ["HIGH", "MEDIUM", "LOW", "NONE"]:
            print(f"- {label}: `{counts.get(label, 0)}`", file=handle)
        print("", file=handle)
        print("## Classification counts", file=handle)
        print("", file=handle)
        for cls, count in class_counts.most_common():
            print(f"- `{cls}`: `{count}`", file=handle)
        if intron_lens:
            print("", file=handle)
            print("## Candidate intron length distribution", file=handle)
            print("", file=handle)
            print(f"- Min: `{min(intron_lens)}` bp", file=handle)
            print(f"- Median: `{int(round(median(intron_lens)))}` bp", file=handle)
            print(f"- Max: `{max(intron_lens)}` bp", file=handle)
        print("", file=handle)
        print("## Top candidates", file=handle)
        print("", file=handle)
        for result in sorted(results, key=lambda r: (confidence_rank(r.confidence), r.support_subjects, r.support_taxa), reverse=True)[:30]:
            if result.confidence == "NONE":
                continue
            print(f"- `{result.query_id}`: {result.confidence}, intron={result.intron_start}-{result.intron_end} ({result.intron_len} bp), support_subjects={result.support_subjects}, support_taxa={result.support_taxa}, best_subject=`{result.best_subject}`", file=handle)

def print_summary(results: List[QueryResult], args: argparse.Namespace) -> None:
    """Print terminal summary."""

    counts = Counter(r.confidence for r in results)
    table = Table(title="ivsBLASTn summary")
    table.add_column("Metric", style="bold")
    table.add_column("Value", justify="right")
    table.add_row("Queries analyzed", f"{len(results):,}")
    table.add_row("High confidence", f"{counts.get('HIGH', 0):,}")
    table.add_row("Medium confidence", f"{counts.get('MEDIUM', 0):,}")
    table.add_row("Low confidence", f"{counts.get('LOW', 0):,}")
    table.add_row("No signal", f"{counts.get('NONE', 0):,}")
    table.add_row("Summary TSV", str(args.summary_tsv))
    table.add_row("Intron-free FASTA", str(args.intron_free_fa))
    table.add_row("Introns FASTA", str(args.introns_fa))
    CONSOLE.print(table)
