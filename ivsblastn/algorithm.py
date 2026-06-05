from __future__ import annotations

import argparse
import heapq
from statistics import median
from typing import Dict, List, Optional, Tuple

from .models import BlastQueryStats, HSP, QueryResult, SupportPair
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
    """Return IVS geometry if two HSPs support an insertion-like gap."""

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


def support_pairs_for_subject(query_id: str, subject_id: str, hsps: List[HSP], taxonomy: Dict[str, str], args: argparse.Namespace) -> List[SupportPair]:
    """Find IVS-like HSP pairs for one query-subject comparison."""

    if len(hsps) < 2:
        return []
    pairs: List[SupportPair] = []
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
            pairs.append(candidate)
    return pairs


def best_support_pair_for_subject(query_id: str, subject_id: str, hsps: List[HSP], taxonomy: Dict[str, str], args: argparse.Namespace) -> Optional[SupportPair]:
    """Find the best IVS-like HSP pair for one query-subject pair."""

    pairs = support_pairs_for_subject(query_id, subject_id, hsps, taxonomy, args)
    if not pairs:
        return None
    return max(pairs, key=lambda p: p.pair_score)


def top_subjects_by_bitscore(subject_hsps: Dict[str, List[HSP]], top_subjects: int) -> Dict[str, List[HSP]]:
    """Keep top subjects by summed HSP bitscore."""

    if top_subjects >= len(subject_hsps):
        return dict(subject_hsps)
    indexed_items = list(enumerate(subject_hsps.items()))
    score = lambda item: (sum(h.bitscore for h in item[1][1]), -item[0])
    ranked = heapq.nlargest(top_subjects, indexed_items, key=score)
    return dict(item for _index, item in ranked)


def blast_stats_from_retained_hsps(subject_hsps: Dict[str, List[HSP]]) -> BlastQueryStats:
    """Build fallback BLAST stats when only retained HSPs are available."""

    retained_hsps = sum(len(hsps) for hsps in subject_hsps.values())
    retained_subjects = len(subject_hsps)
    return BlastQueryStats(
        raw_hsps=retained_hsps,
        raw_subjects=retained_subjects,
        retained_hsps=retained_hsps,
        retained_subjects=retained_subjects,
    )


def best_blast_subject_fields(subject_hsps: Dict[str, List[HSP]], taxonomy: Dict[str, str]) -> Dict[str, object]:
    """Return summary fields for the strongest retained BLAST subject."""

    if not subject_hsps:
        return {
            "best_blast_subject": "",
            "best_blast_pident": 0.0,
            "best_blast_bitscore": 0.0,
            "best_blast_taxonomy": "",
        }
    best_subject, best_hsps = max(
        subject_hsps.items(),
        key=lambda item: (sum(h.bitscore for h in item[1]), max(h.bitscore for h in item[1])),
    )
    best_hsp = max(best_hsps, key=lambda h: h.bitscore)
    return {
        "best_blast_subject": best_subject,
        "best_blast_pident": best_hsp.pident,
        "best_blast_bitscore": sum(h.bitscore for h in best_hsps),
        "best_blast_taxonomy": taxonomy.get(best_subject, ""),
    }


def blast_result_fields(blast_status: str, stats: BlastQueryStats, top_subjects: Dict[str, List[HSP]], taxonomy: Dict[str, str]) -> Dict[str, object]:
    """Return QueryResult BLAST summary fields."""

    fields = best_blast_subject_fields(top_subjects, taxonomy)
    fields.update(
        {
            "blast_status": blast_status,
            "blast_raw_hsps": stats.raw_hsps,
            "blast_raw_subjects": stats.raw_subjects,
            "blast_retained_hsps": stats.retained_hsps,
            "blast_retained_subjects": stats.retained_subjects,
            "blast_subjects_analyzed": len(top_subjects),
        }
    )
    return fields


def deduplicate_cluster_subjects(cluster: List[SupportPair]) -> List[SupportPair]:
    """Keep the best support pair per subject within one IVS cluster."""

    best_by_subject: Dict[str, SupportPair] = {}
    for pair in sorted(cluster, key=lambda p: p.pair_score, reverse=True):
        best_by_subject.setdefault(pair.subject_id, pair)
    return sorted(best_by_subject.values(), key=lambda p: p.pair_score, reverse=True)


def cluster_support_pairs(pairs: List[SupportPair], breakpoint_window: int) -> List[List[SupportPair]]:
    """Cluster support pairs by similar query IVS coordinates."""

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
        return "HIGH_CONFIDENCE_16S_IVS", "HIGH", [f"high_support_subjects>={args.high_support_subjects}", f"high_support_taxa>={args.high_support_taxa}"]
    if support_subjects >= args.medium_support_subjects and support_taxa >= args.medium_support_taxa:
        return "MEDIUM_CONFIDENCE_16S_IVS", "MEDIUM", [f"medium_support_subjects>={args.medium_support_subjects}", f"medium_support_taxa>={args.medium_support_taxa}"]
    if support_subjects >= args.min_support_subjects:
        return "LOW_CONFIDENCE_16S_IVS", "LOW", [f"min_support_subjects>={args.min_support_subjects}"]
    return "NO_IVS_SIGNAL", "NONE", ["no_supported_ivs_cluster"]


def no_signal_result(
    query_id: str,
    query_len: int,
    stats: BlastQueryStats,
    top_subjects: Dict[str, List[HSP]],
    taxonomy: Dict[str, str],
    args: argparse.Namespace,
) -> QueryResult:
    """Build a no-IVS result with BLAST diagnostic context."""

    if stats.raw_hsps == 0:
        blast_status = "NO_BLAST_HIT"
        reasons = ["no_blast_hsp_reported"]
    elif stats.retained_hsps == 0:
        blast_status = "BLAST_HITS_FILTERED"
        reasons = ["blast_hsps_failed_filters", f"min_pident={args.min_pident}", f"min_hsp_len={args.min_hsp_len}"]
    else:
        blast_status = "BLAST_HIT_NO_IVS_PATTERN"
        reasons = ["blast_hsps_present_but_no_supported_ivs_gap_pattern"]
    reasons.extend(
        [
            f"blast_raw_hsps={stats.raw_hsps}",
            f"blast_raw_subjects={stats.raw_subjects}",
            f"blast_retained_hsps={stats.retained_hsps}",
            f"blast_retained_subjects={stats.retained_subjects}",
            f"blast_subjects_analyzed={len(top_subjects)}",
        ]
    )
    return QueryResult(
        query_id=query_id,
        query_len=query_len,
        classification="NO_IVS_SIGNAL",
        confidence="NONE",
        reasons=reasons,
        **blast_result_fields(blast_status, stats, top_subjects, taxonomy),
    )


def cluster_to_result(
    query_id: str,
    query_len: int,
    cluster: List[SupportPair],
    stats: BlastQueryStats,
    top_subjects: Dict[str, List[HSP]],
    taxonomy: Dict[str, str],
    args: argparse.Namespace,
    ivs_index: int,
    ivs_count: int,
) -> QueryResult:
    """Convert one support-pair cluster into one IVS result row."""

    best_cluster = deduplicate_cluster_subjects(cluster)
    intron_start = int(round(median([p.intron_start for p in best_cluster])))
    intron_end = int(round(median([p.intron_end for p in best_cluster])))
    intron_len = intron_end - intron_start + 1
    support_subjects = len({p.subject_id for p in best_cluster})
    support_taxa = len({p.taxon_at_rank for p in best_cluster if p.taxon_at_rank != "NA"})
    support_species = unique_taxa_at_rank(best_cluster, "species")
    support_genera = unique_taxa_at_rank(best_cluster, "genus")
    classification, confidence, reasons = classify_confidence(support_subjects, support_taxa, args)
    reasons.extend(
        [
            f"algorithm={args.algorithm}",
            "blast_status=IVS_PATTERN_DETECTED",
            f"ivs_index={ivs_index}",
            f"ivs_count={ivs_count}",
            f"support_subjects={support_subjects}",
            f"support_taxa_at_{args.tax_rank}={support_taxa}",
            f"breakpoint_window={args.breakpoint_window}bp",
        ]
    )
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
        ivs_index=ivs_index,
        ivs_count=ivs_count,
        **blast_result_fields("IVS_PATTERN_DETECTED", stats, top_subjects, taxonomy),
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


def intervals_overlap(left: QueryResult, right: QueryResult) -> bool:
    """Return True when two IVS intervals overlap on the query."""

    return left.intron_start <= right.intron_end and right.intron_start <= left.intron_end


def analyze_query_all(
    query_id: str,
    subject_hsps: Dict[str, List[HSP]],
    query_len: int,
    taxonomy: Dict[str, str],
    args: argparse.Namespace,
    blast_stats: Optional[BlastQueryStats] = None,
) -> List[QueryResult]:
    """Analyze one query sequence and return all supported non-overlapping IVSs."""

    stats = blast_stats or blast_stats_from_retained_hsps(subject_hsps)
    top_subjects = top_subjects_by_bitscore(subject_hsps, args.top_subjects)
    support_pairs: List[SupportPair] = []
    for subject_id, hsps in top_subjects.items():
        support_pairs.extend(support_pairs_for_subject(query_id, subject_id, hsps, taxonomy, args))
    if not support_pairs:
        return [no_signal_result(query_id, query_len, stats, top_subjects, taxonomy, args)]

    clusters = cluster_support_pairs(support_pairs, args.breakpoint_window)
    candidate_results = [
        cluster_to_result(query_id, query_len, cluster, stats, top_subjects, taxonomy, args, 0, 0)
        for cluster in clusters
    ]
    candidate_results = [result for result in candidate_results if confidence_rank(result.confidence) > 0]
    if not candidate_results:
        result = no_signal_result(query_id, query_len, stats, top_subjects, taxonomy, args)
        result.reasons.append("ivs_gap_clusters_below_support_threshold")
        return [result]

    ranked_candidates = sorted(
        candidate_results,
        key=lambda r: (confidence_rank(r.confidence), r.support_subjects, r.support_taxa, r.mean_bitscore),
        reverse=True,
    )
    selected: List[QueryResult] = []
    for candidate in ranked_candidates:
        if any(intervals_overlap(candidate, existing) for existing in selected):
            continue
        selected.append(candidate)

    selected.sort(key=lambda r: (r.intron_start, r.intron_end))
    ivs_count = len(selected)
    ivs_free_len = max(0, query_len - sum(result.intron_len for result in selected))
    results: List[QueryResult] = []
    for ivs_index, result in enumerate(selected, start=1):
        cluster = result.support_pairs or []
        final_result = cluster_to_result(query_id, query_len, cluster, stats, top_subjects, taxonomy, args, ivs_index, ivs_count)
        final_result.intron_free_len = ivs_free_len
        results.append(final_result)
    return results


def analyze_query(
    query_id: str,
    subject_hsps: Dict[str, List[HSP]],
    query_len: int,
    taxonomy: Dict[str, str],
    args: argparse.Namespace,
    blast_stats: Optional[BlastQueryStats] = None,
) -> QueryResult:
    """Analyze one query sequence and return the first result for compatibility."""

    return analyze_query_all(query_id, subject_hsps, query_len, taxonomy, args, blast_stats)[0]


def confidence_rank(label: str) -> int:
    """Rank confidence labels."""

    return {"NONE": 0, "LOW": 1, "MEDIUM": 2, "HIGH": 3}.get(label, 0)


def is_intron_result(result: QueryResult, min_confidence: str) -> bool:
    """Return True if a result passes output confidence threshold."""

    return bool(result.intron_start and result.intron_end) and confidence_rank(result.confidence) >= confidence_rank(min_confidence)
