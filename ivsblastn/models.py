from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

@dataclass(frozen=True)
class HSP:
    """One BLASTN high-scoring pair."""

    qseqid: str
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: str
    bitscore: float

    @property
    def qlo(self) -> int:
        return min(self.qstart, self.qend)

    @property
    def qhi(self) -> int:
        return max(self.qstart, self.qend)

    @property
    def slo(self) -> int:
        return min(self.sstart, self.send)

    @property
    def shi(self) -> int:
        return max(self.sstart, self.send)

    @property
    def qdir(self) -> int:
        return 1 if self.qend >= self.qstart else -1

    @property
    def sdir(self) -> int:
        return 1 if self.send >= self.sstart else -1

    @property
    def orientation(self) -> int:
        return self.qdir * self.sdir


@dataclass
class BlastQueryStats:
    """Per-query BLAST row counts before and after HSP filtering."""

    raw_hsps: int = 0
    raw_subjects: int = 0
    retained_hsps: int = 0
    retained_subjects: int = 0


@dataclass(frozen=True)
class SupportPair:
    """Best IVS-like HSP pair for one query-subject comparison."""

    query_id: str
    subject_id: str
    hsp1: HSP
    hsp2: HSP
    query_gap: int
    subject_gap: int
    ivs_start: int
    ivs_end: int
    ivs_len: int
    pair_score: float
    taxonomy: str
    taxon_at_rank: str


@dataclass(frozen=True)
class ReferenceRecord:
    """One selected reference sequence candidate."""

    order: int
    seq_id: str
    taxonomy: str
    seq: str


@dataclass
class QueryResult:
    """Final IVS detection result for one query or one IVS event."""

    query_id: str
    query_len: int
    classification: str
    confidence: str
    ivs_index: int = 0
    ivs_count: int = 0
    blast_status: str = "NO_BLAST_HIT"
    blast_raw_hsps: int = 0
    blast_raw_subjects: int = 0
    blast_retained_hsps: int = 0
    blast_retained_subjects: int = 0
    blast_subjects_analyzed: int = 0
    best_blast_subject: str = ""
    best_blast_pident: float = 0.0
    best_blast_bitscore: float = 0.0
    best_blast_taxonomy: str = ""
    ivs_start: int = 0
    ivs_end: int = 0
    ivs_len: int = 0
    exon1_start: int = 0
    exon1_end: int = 0
    exon2_start: int = 0
    exon2_end: int = 0
    ivs_free_len: int = 0
    support_subjects: int = 0
    support_taxa: int = 0
    support_species: int = 0
    support_genera: int = 0
    median_subject_gap: float = 0.0
    median_pident: float = 0.0
    mean_bitscore: float = 0.0
    best_subject: str = ""
    best_subject_taxonomy: str = ""
    reasons: Optional[List[str]] = None
    support_pairs: Optional[List[SupportPair]] = None

    def __post_init__(self) -> None:
        if self.reasons is None:
            self.reasons = []
        if self.support_pairs is None:
            self.support_pairs = []
