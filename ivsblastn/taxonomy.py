from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Optional

from .logging import LOG

def parse_taxonomy(path: Optional[Path]) -> Dict[str, str]:
    """Read subject taxonomy mapping."""

    if path is None:
        return {}
    taxonomy: Dict[str, str] = {}
    with path.open("rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=chr(9))
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 2:
                continue
            taxonomy[row[0]] = row[1]
    LOG.info("Parsed taxonomy mappings: %s subjects", len(taxonomy))
    return taxonomy

def get_taxon_at_rank(taxonomy: str, rank: str) -> str:
    """Extract taxon at a rank from semicolon-delimited taxonomy."""

    if not taxonomy:
        return "NA"
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    parts = [p.strip() for p in taxonomy.split(";") if p.strip()]
    index = {name: i for i, name in enumerate(ranks)}.get(rank, 6)
    if index < len(parts):
        return parts[index]
    return parts[-1] if parts else "NA"
