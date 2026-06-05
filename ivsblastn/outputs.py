from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, List, Tuple

from rich.table import Table

from .algorithm import confidence_rank, is_ivs_result
from .fasta import open_text_auto_write, write_fasta_record
from .logging import CONSOLE
from .models import QueryResult, SupportPair

def summary_row(result: QueryResult, tax_rank: str) -> Dict[str, str]:
    """Convert one query result to a TSV row."""

    return {
        "query_id": result.query_id,
        "query_len": str(result.query_len),
        "ivs_index": str(result.ivs_index) if result.ivs_index else "",
        "ivs_count": str(result.ivs_count) if result.ivs_count else "",
        "blast_status": result.blast_status,
        "blast_raw_hsps": str(result.blast_raw_hsps),
        "blast_raw_subjects": str(result.blast_raw_subjects),
        "blast_retained_hsps": str(result.blast_retained_hsps),
        "blast_retained_subjects": str(result.blast_retained_subjects),
        "blast_subjects_analyzed": str(result.blast_subjects_analyzed),
        "best_blast_subject": result.best_blast_subject,
        "best_blast_pident": f"{result.best_blast_pident:.2f}",
        "best_blast_bitscore": f"{result.best_blast_bitscore:.2f}",
        "best_blast_taxonomy": result.best_blast_taxonomy,
        "classification": result.classification,
        "confidence": result.confidence,
        "ivs_start": str(result.ivs_start) if result.ivs_start else "",
        "ivs_end": str(result.ivs_end) if result.ivs_end else "",
        "ivs_len": str(result.ivs_len) if result.ivs_len else "",
        "exon1": f"{result.exon1_start}-{result.exon1_end}" if result.exon1_start else "",
        "exon2": f"{result.exon2_start}-{result.exon2_end}" if result.exon2_start else "",
        "ivs_free_len": str(result.ivs_free_len) if result.ivs_free_len else "",
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


def support_row(pair: SupportPair, ivs_index: int) -> Dict[str, str]:
    """Convert one support pair to a TSV row."""

    h1, h2 = sorted([pair.hsp1, pair.hsp2], key=lambda h: (h.qlo, h.qhi))
    return {
        "query_id": pair.query_id,
        "ivs_index": str(ivs_index),
        "subject_id": pair.subject_id,
        "taxonomy": pair.taxonomy,
        "taxon_at_rank": pair.taxon_at_rank,
        "ivs_start": str(pair.ivs_start),
        "ivs_end": str(pair.ivs_end),
        "ivs_len": str(pair.ivs_len),
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
        "ivs_index",
        "ivs_count",
        "blast_status",
        "blast_raw_hsps",
        "blast_raw_subjects",
        "blast_retained_hsps",
        "blast_retained_subjects",
        "blast_subjects_analyzed",
        "best_blast_subject",
        "best_blast_pident",
        "best_blast_bitscore",
        "best_blast_taxonomy",
        "classification",
        "confidence",
        "ivs_start",
        "ivs_end",
        "ivs_len",
        "exon1",
        "exon2",
        "ivs_free_len",
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
        "ivs_index",
        "subject_id",
        "taxonomy",
        "taxon_at_rank",
        "ivs_start",
        "ivs_end",
        "ivs_len",
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
                writer.writerow(support_row(pair, result.ivs_index))


def write_fasta_outputs(args: argparse.Namespace, results: List[QueryResult], seqs: Dict[str, str]) -> None:
    """Write IVS-free 16S and IVS FASTA outputs."""

    results_by_query: Dict[str, List[QueryResult]] = defaultdict(list)
    for result in results:
        results_by_query[result.query_id].append(result)
    with open_text_auto_write(args.ivs_free_fa) as free_out, open_text_auto_write(args.ivs_fa) as ivs_out:
        for query_id, seq in seqs.items():
            if not seq:
                continue
            query_results = sorted(results_by_query.get(query_id, []), key=lambda r: (r.ivs_start, r.ivs_end))
            if not query_results:
                write_fasta_record(free_out, f"{query_id}|ivs_free|ivs=none|confidence=NONE|action=unchanged", seq)
                continue
            passing_results = [result for result in query_results if is_ivs_result(result, args.min_output_confidence)]
            candidate_results = [result for result in query_results if result.ivs_start and result.ivs_end]
            if passing_results:
                seq_parts: List[str] = []
                cursor = 1
                for result in passing_results:
                    seq_parts.append(seq[cursor - 1 : result.ivs_start - 1])
                    cursor = result.ivs_end + 1
                seq_parts.append(seq[cursor - 1 :])
                intervals = ",".join(f"{result.ivs_start}-{result.ivs_end}" for result in passing_results)
                confidence = max((result.confidence for result in passing_results), key=confidence_rank)
                write_fasta_record(free_out, f"{query_id}|ivs_free|ivs={intervals}|ivs_count={len(passing_results)}|confidence={confidence}|action=removed", "".join(seq_parts))
                for result in passing_results:
                    ivs = seq[result.ivs_start - 1 : result.ivs_end]
                    write_fasta_record(ivs_out, f"{result.query_id}|ivs_{result.ivs_index or 1}|{result.ivs_start}-{result.ivs_end}|len={result.ivs_len}|confidence={result.confidence}", ivs)
            elif candidate_results:
                intervals = ",".join(f"{result.ivs_start}-{result.ivs_end}" for result in candidate_results)
                confidence = max((result.confidence for result in candidate_results), key=confidence_rank)
                write_fasta_record(free_out, f"{query_id}|ivs_free|ivs={intervals}|ivs_count={len(candidate_results)}|confidence={confidence}|action=unchanged_below_min_output_confidence", seq)
            else:
                confidence = max((result.confidence for result in query_results), key=confidence_rank)
                write_fasta_record(free_out, f"{query_id}|ivs_free|ivs=none|confidence={confidence}|action=unchanged", seq)


def write_bed_outputs(args: argparse.Namespace, results: List[QueryResult]) -> None:
    """Write BED files."""

    passing_by_query: Dict[str, List[QueryResult]] = defaultdict(list)
    for result in results:
        if is_ivs_result(result, args.min_output_confidence):
            passing_by_query[result.query_id].append(result)

    with args.ivs_bed.open("wt", encoding="utf-8", newline="") as ivs_file, args.exons_bed.open("wt", encoding="utf-8", newline="") as exon_file:
        ivs_writer = csv.writer(ivs_file, delimiter=chr(9), lineterminator=chr(10))
        exon_writer = csv.writer(exon_file, delimiter=chr(9), lineterminator=chr(10))
        for query_id, query_results in passing_by_query.items():
            sorted_ivs = sorted(query_results, key=lambda r: (r.ivs_start, r.ivs_end))
            for result in sorted_ivs:
                ivs_label = f"ivs_{result.ivs_index or 1}"
                ivs_writer.writerow([result.query_id, result.ivs_start - 1, result.ivs_end, f"{result.query_id}|{ivs_label}|{result.ivs_start}-{result.ivs_end}|confidence={result.confidence}", ".", "+"])

            exon_intervals = retained_exon_intervals(sorted_ivs[0].query_len, [(r.ivs_start, r.ivs_end) for r in sorted_ivs])
            confidence = max((result.confidence for result in sorted_ivs), key=confidence_rank)
            for exon_index, (start, end) in enumerate(exon_intervals, start=1):
                exon_writer.writerow([query_id, start - 1, end, f"{query_id}|exon_{exon_index}|{start}-{end}|ivs_removed={len(sorted_ivs)}|confidence={confidence}", ".", "+"])


def retained_exon_intervals(query_len: int, ivs_intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Return 1-based closed intervals left after removing sorted IVS intervals."""

    intervals: List[Tuple[int, int]] = []
    cursor = 1
    for start, end in sorted(ivs_intervals):
        if cursor <= start - 1:
            intervals.append((cursor, start - 1))
        cursor = max(cursor, end + 1)
    if cursor <= query_len:
        intervals.append((cursor, query_len))
    return intervals


def write_report(path: Path, results: List[QueryResult], args: argparse.Namespace) -> None:
    """Write Markdown report."""

    counts = Counter(r.confidence for r in results)
    class_counts = Counter(r.classification for r in results)
    blast_status_counts = Counter(r.blast_status for r in results)
    ivs_lens = [r.ivs_len for r in results if r.ivs_len > 0]
    candidates = [r for r in results if r.ivs_start and r.ivs_end]
    output_candidates = [r for r in results if is_ivs_result(r, args.min_output_confidence)]
    query_count = len({r.query_id for r in results})

    with path.open("wt", encoding="utf-8") as handle:
        print("# ivsBLASTn report", file=handle)
        print("", file=handle)
        print("## Inputs", file=handle)
        print("", file=handle)
        print(f"- Query FASTA: `{args.query}`", file=handle)
        print(f"- BLAST table: `{getattr(args, 'blast', None)}`", file=handle)
        print(f"- BLAST database: `{getattr(args, 'db', None)}`", file=handle)
        if hasattr(args, "ref_fasta"):
            print(f"- Reference FASTA: `{args.ref_fasta}`", file=handle)
        print(f"- Taxonomy table: `{getattr(args, 'taxonomy', None)}`", file=handle)
        print(f"- Output directory: `{args.outdir}`", file=handle)
        print("", file=handle)
        print("## Algorithm", file=handle)
        print("", file=handle)
        print(f"- Algorithm: `{args.algorithm}`", file=handle)
        print("- Detects HSP1--query insertion--HSP2 patterns where reference subject coordinates are nearly continuous.", file=handle)
        print("", file=handle)
        print("## Parameters", file=handle)
        print("", file=handle)
        parameter_names = [
            "min_pident",
            "min_hsp_len",
            "min_ivs_len",
            "max_ivs_len",
            "max_ref_gap",
            "max_query_overlap",
            "breakpoint_window",
            "top_subjects",
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
        ]
        if hasattr(args, "ref_domains"):
            parameter_names.extend(
                [
                    "ref_domains",
                    "ref_per_species",
                    "ref_unclear_per_genus",
                    "clean_ref_ivs",
                    "ref_clean_min_confidence",
                    "ref_self_blast_max_target_seqs",
                    "ref_self_blast_max_hsps",
                ]
            )
        for name in parameter_names:
            print(f"- `{name}`: `{getattr(args, name, 'NA')}`", file=handle)
        print("", file=handle)
        print("## Summary", file=handle)
        print("", file=handle)
        print(f"- Query sequences analyzed: `{query_count}`", file=handle)
        print(f"- Result rows: `{len(results)}`", file=handle)
        print(f"- Candidate IVS events detected: `{len(candidates)}`", file=handle)
        print(f"- Candidate IVSs passing `{args.min_output_confidence}` output threshold: `{len(output_candidates)}`", file=handle)
        for label in ["HIGH", "MEDIUM", "LOW", "NONE"]:
            print(f"- {label}: `{counts.get(label, 0)}`", file=handle)
        print("", file=handle)
        print("## Classification counts", file=handle)
        print("", file=handle)
        for cls, count in class_counts.most_common():
            print(f"- `{cls}`: `{count}`", file=handle)
        print("", file=handle)
        print("## BLAST status counts", file=handle)
        print("", file=handle)
        for status, count in blast_status_counts.most_common():
            print(f"- `{status}`: `{count}`", file=handle)
        if ivs_lens:
            print("", file=handle)
            print("## Candidate IVS length distribution", file=handle)
            print("", file=handle)
            print(f"- Min: `{min(ivs_lens)}` bp", file=handle)
            print(f"- Median: `{int(round(median(ivs_lens)))}` bp", file=handle)
            print(f"- Max: `{max(ivs_lens)}` bp", file=handle)
        print("", file=handle)
        print("## Top candidates", file=handle)
        print("", file=handle)
        for result in sorted(results, key=lambda r: (confidence_rank(r.confidence), r.support_subjects, r.support_taxa), reverse=True)[:30]:
            if result.confidence == "NONE":
                continue
            print(f"- `{result.query_id}` ivs_{result.ivs_index or 1}: {result.confidence}, IVS={result.ivs_start}-{result.ivs_end} ({result.ivs_len} bp), support_subjects={result.support_subjects}, support_taxa={result.support_taxa}, best_subject=`{result.best_subject}`", file=handle)

def print_summary(results: List[QueryResult], args: argparse.Namespace) -> None:
    """Print terminal summary."""

    counts = Counter(r.confidence for r in results)
    blast_status_counts = Counter(r.blast_status for r in results)
    query_count = len({r.query_id for r in results})
    candidates = sum(1 for r in results if r.ivs_start and r.ivs_end)
    output_candidates = sum(1 for r in results if is_ivs_result(r, args.min_output_confidence))
    table = Table(title="ivsBLASTn summary")
    table.add_column("Metric", style="bold")
    table.add_column("Value", justify="right")
    table.add_row("Queries analyzed", f"{query_count:,}")
    table.add_row("Result rows", f"{len(results):,}")
    table.add_row("High confidence IVS events", f"{counts.get('HIGH', 0):,}")
    table.add_row("Medium confidence IVS events", f"{counts.get('MEDIUM', 0):,}")
    table.add_row("Low confidence IVS events", f"{counts.get('LOW', 0):,}")
    table.add_row("No signal", f"{counts.get('NONE', 0):,}")
    table.add_row("No BLAST hit", f"{blast_status_counts.get('NO_BLAST_HIT', 0):,}")
    table.add_row("BLAST hit, no IVS pattern", f"{blast_status_counts.get('BLAST_HIT_NO_IVS_PATTERN', 0):,}")
    table.add_row("Candidate IVS events detected", f"{candidates:,}")
    table.add_row(f"IVSs written/removed (>= {args.min_output_confidence})", f"{output_candidates:,}")
    table.add_row("Summary TSV", str(args.summary_tsv))
    table.add_row("IVS-free FASTA", str(args.ivs_free_fa))
    table.add_row("IVSs FASTA", str(args.ivs_fa))
    CONSOLE.print(table)
