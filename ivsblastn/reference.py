from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from .algorithm import analyze_query_all, confidence_rank, is_intron_result
from .blast import make_blast_db, parse_blast_with_stats, run_blastn_to_file
from .fasta import (
    clean_dna_sequence,
    genus_key_from_taxonomy,
    iter_fasta_records,
    is_strict_atgc,
    parse_silva_header,
    read_fasta,
    species_key_from_taxonomy,
    taxonomy_domain,
    write_fasta_record,
)
from .logging import LOG
from .models import QueryResult, ReferenceRecord
from .outputs import write_report, write_summary, write_supporting_hsps
from .paths import output_path
from .taxonomy import parse_taxonomy


def write_reference_introns_fasta(path: Path, ref_seqs: Dict[str, str], results: List[QueryResult], min_confidence: str) -> int:
    """Write reference IVS sequences removed during self-cleaning."""

    written = 0
    with path.open("wt", encoding="utf-8") as handle:
        for result in results:
            if not is_intron_result(result, min_confidence):
                continue
            seq = ref_seqs.get(result.query_id, "")
            if not seq:
                continue
            intron = seq[result.intron_start - 1 : result.intron_end]
            write_fasta_record(
                handle,
                f"{result.query_id}|reference_ivs_{result.ivs_index or 1}|{result.intron_start}-{result.intron_end}|len={result.intron_len}|confidence={result.confidence}|action=removed",
                intron,
            )
            written += 1
    return written


def preprocess_reference(args: argparse.Namespace) -> Tuple[Path, Path, Path]:
    """Filter SILVA reference, write clean FASTA/taxonomy, and build BLAST DB."""

    allowed_domains = {x.strip() for x in args.ref_domains.split(",") if x.strip()}
    selected_by_species: Dict[str, List[ReferenceRecord]] = defaultdict(list)
    selected_unclear_by_genus: Dict[str, List[ReferenceRecord]] = defaultdict(list)
    eligible_clear = 0
    eligible_unclear = 0
    skipped_domain = 0
    skipped_unclear_species = 0
    skipped_species_cap = 0
    skipped_unclear_genus_cap = 0
    skipped_non_atgc = 0
    skipped_empty = 0

    def keep_capped_record(selection: Dict[str, List[ReferenceRecord]], key: str, record: ReferenceRecord, limit: int) -> None:
        if limit == 0:
            return
        selected = selection[key]
        selected.append(record)
        selected.sort(key=lambda r: (-len(r.seq), r.order))
        del selected[limit:]

    def write_record(record: ReferenceRecord) -> None:
        write_fasta_record(fa_out, record.seq_id, record.seq)
        tax_writer.writerow([record.seq_id, record.taxonomy])

    def keep_selected_species_record(record: ReferenceRecord, species_key: str) -> None:
        selected = selected_by_species[species_key]
        selected.append(record)
        selected.sort(key=lambda r: (-len(r.seq), r.order))
        del selected[args.ref_per_species :]

    with args.raw_ref_fa.open("wt", encoding="utf-8") as fa_out, args.raw_ref_tax.open("wt", encoding="utf-8", newline="") as tax_file:
        tax_writer = csv.writer(tax_file, delimiter=chr(9), lineterminator=chr(10))
        for record_order, (_record_id, description, raw_seq) in enumerate(iter_fasta_records(args.ref_fasta), start=1):
            seq_id, taxonomy = parse_silva_header(description)
            domain = taxonomy_domain(taxonomy)
            if domain not in allowed_domains:
                skipped_domain += 1
                continue
            seq = clean_dna_sequence([raw_seq])
            if not seq:
                skipped_empty += 1
                continue
            if not is_strict_atgc(seq):
                skipped_non_atgc += 1
                continue
            species_key = species_key_from_taxonomy(taxonomy)
            record = ReferenceRecord(order=record_order, seq_id=seq_id, taxonomy=taxonomy, seq=seq)
            if species_key is not None:
                eligible_clear += 1
                if args.ref_per_species == 0:
                    write_record(record)
                else:
                    keep_selected_species_record(record, species_key)
                continue

            genus_key = genus_key_from_taxonomy(taxonomy)
            if genus_key is None or args.ref_unclear_per_genus == 0:
                skipped_unclear_species += 1
                continue
            eligible_unclear += 1
            keep_capped_record(selected_unclear_by_genus, genus_key, record, args.ref_unclear_per_genus)

        selected_records: List[ReferenceRecord] = []
        if args.ref_per_species > 0:
            selected_species_records = [record for records in selected_by_species.values() for record in records]
            skipped_species_cap = eligible_clear - len(selected_species_records)
            selected_records.extend(selected_species_records)
        if args.ref_unclear_per_genus > 0:
            selected_unclear_records = [record for records in selected_unclear_by_genus.values() for record in records]
            skipped_unclear_genus_cap = eligible_unclear - len(selected_unclear_records)
            selected_records.extend(selected_unclear_records)
        for record in sorted(selected_records, key=lambda r: r.order):
            if args.ref_per_species == 0 and species_key_from_taxonomy(record.taxonomy) is not None:
                continue
            write_record(record)

    kept_clear = eligible_clear - skipped_species_cap
    kept_unclear = eligible_unclear - skipped_unclear_genus_cap
    kept = kept_clear + kept_unclear
    LOG.info("Reference preprocessing kept %s sequences", kept)
    LOG.info("Kept clear species records: %s", kept_clear)
    LOG.info("Kept unclear-species genus fallback records: %s", kept_unclear)
    LOG.info("Skipped by domain filter: %s", skipped_domain)
    LOG.info("Skipped unclear species without genus fallback: %s", skipped_unclear_species)
    LOG.info("Skipped by per-species cap: %s", skipped_species_cap)
    LOG.info("Skipped by unclear-genus cap: %s", skipped_unclear_genus_cap)
    LOG.info("Skipped empty sequences: %s", skipped_empty)
    LOG.info("Skipped sequences containing non-ATGC characters: %s", skipped_non_atgc)
    make_blast_db(args.raw_ref_fa, args.raw_ref_db, args.makeblastdb_bin)
    return args.raw_ref_fa, args.raw_ref_tax, args.raw_ref_db

def clean_reference_introns(args: argparse.Namespace, ref_fa: Path, tax_tsv: Path, db_prefix: Path) -> Tuple[Path, Path, Path]:
    """Self-BLAST reference, remove candidate IVSs, and build cleaned DB."""

    self_blast = run_blastn_to_file(
        query=ref_fa,
        db=db_prefix,
        out_file=args.ref_self_blast,
        args=args,
        max_targets=args.ref_self_blast_max_target_seqs,
        max_hsps=args.ref_self_blast_max_hsps,
        label="reference self-BLASTN",
    )
    ref_seqs = read_fasta(ref_fa)
    ref_taxonomy = parse_taxonomy(tax_tsv)
    ref_blast_by_query, ref_blast_stats_by_query = parse_blast_with_stats(self_blast, args.min_pident, args.min_hsp_len)
    ref_results: List[QueryResult] = []
    for query_id in sorted(ref_seqs):
        ref_results.extend(
            analyze_query_all(
                query_id,
                ref_blast_by_query.get(query_id, {}),
                len(ref_seqs.get(query_id, "")),
                ref_taxonomy,
                args,
                ref_blast_stats_by_query.get(query_id),
            )
        )
    ref_results.sort(key=lambda r: (confidence_rank(r.confidence), r.support_subjects, r.support_taxa, r.query_id), reverse=True)
    write_summary(output_path(args.ref_self_clean_prefix, ".summary.tsv"), ref_results, args.tax_rank)
    write_supporting_hsps(output_path(args.ref_self_clean_prefix, ".supporting_hsps.tsv"), ref_results)
    write_report(output_path(args.ref_self_clean_prefix, ".report.md"), ref_results, args)

    remove_by_id: Dict[str, List[QueryResult]] = defaultdict(list)
    for result in ref_results:
        if is_intron_result(result, args.ref_clean_min_confidence):
            remove_by_id[result.query_id].append(result)
    removed_introns_fa = getattr(args, "ref_self_clean_introns_fa", output_path(args.ref_self_clean_prefix, ".introns.fa"))
    removed_introns = write_reference_introns_fasta(removed_introns_fa, ref_seqs, ref_results, args.ref_clean_min_confidence)
    LOG.info("Reference sequences with candidate IVSs to remove: %s", len(remove_by_id))
    LOG.info("Reference IVS sequences written: %s", removed_introns)
    with args.cleaned_ref_fa.open("wt", encoding="utf-8") as fa_out, args.cleaned_ref_tax.open("wt", encoding="utf-8", newline="") as tax_file:
        tax_writer = csv.writer(tax_file, delimiter=chr(9), lineterminator=chr(10))
        for seq_id in sorted(ref_seqs):
            seq = ref_seqs[seq_id]
            for result in sorted(remove_by_id.get(seq_id, []), key=lambda r: r.intron_start, reverse=True):
                seq = seq[: result.intron_start - 1] + seq[result.intron_end :]
            write_fasta_record(fa_out, seq_id, seq)
            tax_writer.writerow([seq_id, ref_taxonomy.get(seq_id, "")])
    make_blast_db(args.cleaned_ref_fa, args.cleaned_ref_db, args.makeblastdb_bin)
    return args.cleaned_ref_fa, args.cleaned_ref_tax, args.cleaned_ref_db
