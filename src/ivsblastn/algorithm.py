from __future__ import annotations

import argparse
from statistics import median
from typing import Dict, List, Optional, Tuple

from .models import HSP, QueryResult, SupportPair
from .taxonomy import get_taxon_at_rank

def subject_gap_by_query_order(left: HSP, right: HSP) -> Optional[int]:
    """Compute subject gap for two HSPs ordered by query coordinate."""

    if left.orientation != right.orientation:
        return None
    if left.orientation == 1:
        if left.slo > right.slo:
            return None
        return right.slo - left.shi - 1
    if left.orientation == -1:
        if left.slo < right.slo:
            return None
        return left.slo - right.shi - 1
    return None


def hsp_pair_supports_intron(h1: HSP, h2: HSP, args: argparse.Namespace) -> Optional[Tuple[int, int, int, int, int]]:
    """Return intron geometry if two HSPs support an intron-like insertion."""

    left, right = sorted([h1, h2], key=lambda h: (h.qlo, h.qhi))
    query_gap = right.qlo - left.qhi - 1
    if query_gap < 0 and abs(query_gap) > args.max_query_overlap:
        return None
    if query_gap < args.min_intron_len or query_gap > args.max_intron_len:
        return None
    subject_gap = subject_gap_by_query_order(left, right)
    if subject_gap is None:
        return None
    if abs(subject_gap) > args.max_ref_gap:
        return None
    intron_start = left.qhi + 1
    intron_end = right.qlo - 1
    intron_len = intron_end - intron_start + 1
    if intron_len < args.min_intron_len or intron_len > args.max_intron_len:
        return None
    return query_gap, subject_gap, intron_start, intron_end, intron_len


def best_support_pair_for_subject(query_id: str, subject_id: str, hsps: List[HSP], taxonomy: Dict[str, str], args: argparse.Namespace) -> Optional[SupportPair]:
    """Find best intron-like HSP pair for one query-subject pair."""

    if len(hsps) < 2:
        return None
    best_pair: Optional[SupportPair] = None
    tax = taxonomy.get(subject_id, "")
    taxon = get_taxon_at_rank(tax, args.tax_rank)
    hsps_sorted = sorted(hsps, key=lambda h: (-h.bitscore, h.qlo, h.qhi))
    for i in range(len(hsps_sorted)):
        for j in range(i + 1, len(hsps_sorted)):
            support = hsp_pair_supports_intron(hsps_sorted[i], hsps_sorted[j], args)
            if support is None:
                continue
            query_gap, subject_gap, intron_start, intron_end, intron_len = support
            pair_score = hsps_sorted[i].bitscore + hsps_sorted[j].bitscore - 2.0 * abs(subject_gap)
            candidate = SupportPair(
                query_id=query_id,
                subject_id=subject_id,
                hsp1=hsps_sorted[i],
                hsp2=hsps_sorted[j],
                query_gap=query_gap,
                subject_gap=subject_gap,
                intron_start=intron_start,
                intron_end=intron_end,
                intron_len=intron_len,
                pair_score=pair_score,
                taxonomy=tax,
                taxon_at_rank=taxon,
            )
            if best_pair is None or candidate.pair_score > best_pair.pair_score:
                best_pair = candidate
    return best_pair


def top_subjects_by_bitscore(subject_hsps: Dict[str, List[HSP]], top_subjects: int) -> Dict[str, List[HSP]]:
    """Keep top subjects by summed HSP bitscore."""

    ranked = sorted(subject_hsps.items(), key=lambda item: sum(h.bitscore for h in item[1]), reverse=True)
    return dict(ranked[:top_subjects])


def cluster_support_pairs(pairs: List[SupportPair], breakpoint_window: int) -> List[List[SupportPair]]:
    """Cluster support pairs by similar query intron coordinates."""

    clusters: List[List[SupportPair]] = []
    for pair in sorted(pairs, key=lambda p: (p.intron_start, p.intron_end, -p.pair_score)):
        assigned = False
        for cluster in clusters:
            med_start = int(median([p.intron_start for p in cluster]))
            med_end = int(median([p.intron_end for p in cluster]))
            if abs(pair.intron_start - med_start) <= breakpoint_window and abs(pair.intron_end - med_end) <= breakpoint_window:
                cluster.append(pair)
                assigned = True
                break
        if not assigned:
            clusters.append([pair])
    return clusters


def unique_taxa_at_rank(pairs: List[SupportPair], rank: str) -> int:
    """Count unique taxa at a rank."""

    taxa = set()
    for pair in pairs:
        taxon = get_taxon_at_rank(pair.taxonomy, rank)
        if taxon != "NA":
            taxa.add(taxon)
    return len(taxa)


def classify_confidence(support_subjects: int, support_taxa: int, args: argparse.Namespace) -> Tuple[str, str, List[str]]:
    """Assign confidence label from support counts."""

    if support_subjects >= args.high_support_subjects and support_taxa >= args.high_support_taxa:
        return "HIGH_CONFIDENCE_16S_INTRON", "HIGH", [f"high_support_subjects>={args.high_support_subjects}", f"high_support_taxa>={args.high_support_taxa}"]
    if support_subjects >= args.medium_support_subjects and support_taxa >= args.medium_support_taxa:
        return "MEDIUM_CONFIDENCE_16S_INTRON", "MEDIUM", [f"medium_support_subjects>={args.medium_support_subjects}", f"medium_support_taxa>={args.medium_support_taxa}"]
    if support_subjects >= args.min_support_subjects:
        return "LOW_CONFIDENCE_16S_INTRON", "LOW", [f"min_support_subjects>={args.min_support_subjects}"]
    return "NO_INTRON_SIGNAL", "NONE", ["no_supported_intron_cluster"]


def analyze_query(query_id: str, subject_hsps: Dict[str, List[HSP]], query_len: int, taxonomy: Dict[str, str], args: argparse.Namespace) -> QueryResult:
    """Analyze one query sequence."""

    top_subjects = top_subjects_by_bitscore(subject_hsps, args.top_subjects)
    support_pairs: List[SupportPair] = []
    for subject_id, hsps in top_subjects.items():
        pair = best_support_pair_for_subject(query_id, subject_id, hsps, taxonomy, args)
        if pair is not None:
            support_pairs.append(pair)
    if not support_pairs:
        return QueryResult(query_id=query_id, query_len=query_len, classification="NO_INTRON_SIGNAL", confidence="NONE", reasons=["no_subject_supported_hsp_gap_pattern"])

    clusters = cluster_support_pairs(support_pairs, args.breakpoint_window)
    clusters.sort(key=lambda c: (len(c), len({p.taxon_at_rank for p in c if p.taxon_at_rank != "NA"}), sum(p.pair_score for p in c)), reverse=True)
    best_cluster = clusters[0]
    intron_start = int(round(median([p.intron_start for p in best_cluster])))
    intron_end = int(round(median([p.intron_end for p in best_cluster])))
    intron_len = intron_end - intron_start + 1
    support_subjects = len({p.subject_id for p in best_cluster})
    support_taxa = len({p.taxon_at_rank for p in best_cluster if p.taxon_at_rank != "NA"})
    support_species = unique_taxa_at_rank(best_cluster, "species")
    support_genera = unique_taxa_at_rank(best_cluster, "genus")
    classification, confidence, reasons = classify_confidence(support_subjects, support_taxa, args)
    reasons.extend([f"algorithm={args.algorithm}", f"support_subjects={support_subjects}", f"support_taxa_at_{args.tax_rank}={support_taxa}", f"breakpoint_window={args.breakpoint_window}bp"])
    best_pair = sorted(best_cluster, key=lambda p: p.pair_score, reverse=True)[0]
    median_subject_gap = float(median([p.subject_gap for p in best_cluster]))
    median_pident = float(median([(p.hsp1.pident + p.hsp2.pident) / 2.0 for p in best_cluster]))
    mean_bitscore = sum(p.hsp1.bitscore + p.hsp2.bitscore for p in best_cluster) / len(best_cluster)
    exon1_start = 1
    exon1_end = intron_start - 1
    exon2_start = intron_end + 1
    exon2_end = query_len
    intron_free_len = max(0, exon1_end - exon1_start + 1) + max(0, exon2_end - exon2_start + 1)
    return QueryResult(
        query_id=query_id,
        query_len=query_len,
        classification=classification,
        confidence=confidence,
        intron_start=intron_start,
        intron_end=intron_end,
        intron_len=intron_len,
        exon1_start=exon1_start,
        exon1_end=exon1_end,
        exon2_start=exon2_start,
        exon2_end=exon2_end,
        intron_free_len=intron_free_len,
        support_subjects=support_subjects,
        support_taxa=support_taxa,
        support_species=support_species,
        support_genera=support_genera,
        median_subject_gap=median_subject_gap,
        median_pident=median_pident,
        mean_bitscore=mean_bitscore,
        best_subject=best_pair.subject_id,
        best_subject_taxonomy=best_pair.taxonomy,
        reasons=reasons,
        support_pairs=best_cluster,
    )


def confidence_rank(label: str) -> int:
    """Rank confidence labels."""

    return {"NONE": 0, "LOW": 1, "MEDIUM": 2, "HIGH": 3}.get(label, 0)


def is_intron_result(result: QueryResult, min_confidence: str) -> bool:
    """Return True if a result passes output confidence threshold."""

