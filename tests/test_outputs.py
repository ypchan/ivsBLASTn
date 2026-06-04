from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from ivsblastn.models import QueryResult
from ivsblastn.outputs import write_report


class OutputTests(unittest.TestCase):
    def test_write_report_tolerates_missing_optional_parameter_fields(self) -> None:
        with TemporaryDirectory() as tmpdir:
            report = Path(tmpdir) / "report.md"
            args = SimpleNamespace(
                query="query.fa",
                blast="hits.tsv",
                db="ref_db",
                ref_fasta=None,
                taxonomy="ref.tax.tsv",
                outdir="out",
                algorithm="hsp-gap-support",
                min_pident=75.0,
                min_hsp_len=100,
                min_intron_len=25,
                max_intron_len=2000,
                max_ref_gap=30,
                max_query_overlap=20,
                breakpoint_window=30,
                top_subjects=100,
                ref_domains="Archaea,Bacteria",
                ref_per_species=1,
                clean_ref_introns=True,
                ref_clean_min_confidence="LOW",
                ref_self_blast_max_target_seqs=100,
                ref_self_blast_max_hsps=20,
                blast_task="blastn",
                blast_evalue="1e-20",
                tax_rank="genus",
                min_support_subjects=1,
                medium_support_subjects=3,
                medium_support_taxa=3,
                high_support_subjects=10,
                high_support_taxa=3,
                min_output_confidence="LOW",
                threads=4,
            )

            write_report(report, [QueryResult(query_id="q1", query_len=100, classification="NO_INTRON_SIGNAL", confidence="NONE")], args)

            text = report.read_text(encoding="utf-8")
            self.assertIn("`blast_max_hsps`: `NA`", text)


if __name__ == "__main__":
    unittest.main()
