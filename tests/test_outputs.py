from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import unittest

from ivsblastn.algorithm import is_intron_result
from ivsblastn.models import QueryResult
from ivsblastn.outputs import print_summary, write_bed_outputs, write_fasta_outputs, write_report, write_summary


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

            write_report(report, [QueryResult(query_id="q1", query_len=100, classification="NO_IVS_SIGNAL", confidence="NONE")], args)

            text = report.read_text(encoding="utf-8")
            self.assertIn("`blast_max_hsps`: `NA`", text)
            self.assertIn("Candidate IVS events detected: `0`", text)
            self.assertIn("Candidate IVSs passing `LOW` output threshold: `0`", text)

    def test_intron_result_threshold_controls_fasta_removal(self) -> None:
        low = QueryResult(
            query_id="q1",
            query_len=10,
            classification="LOW_CONFIDENCE_16S_IVS",
            confidence="LOW",
            intron_start=4,
            intron_end=6,
            intron_len=3,
            exon1_start=1,
            exon1_end=3,
            exon2_start=7,
            exon2_end=10,
        )
        high = QueryResult(
            query_id="q2",
            query_len=10,
            classification="HIGH_CONFIDENCE_16S_IVS",
            confidence="HIGH",
            intron_start=4,
            intron_end=6,
            intron_len=3,
            exon1_start=1,
            exon1_end=3,
            exon2_start=7,
            exon2_end=10,
        )

        self.assertTrue(is_intron_result(low, "LOW"))
        self.assertFalse(is_intron_result(low, "MEDIUM"))
        self.assertTrue(is_intron_result(high, "MEDIUM"))

    def test_write_fasta_outputs_removes_passing_introns(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            args = SimpleNamespace(
                intron_free_fa=tmp / "intron_free.fa",
                introns_fa=tmp / "introns.fa",
                min_output_confidence="LOW",
            )
            result = QueryResult(
                query_id="q1",
                query_len=10,
                classification="HIGH_CONFIDENCE_16S_IVS",
                ivs_index=1,
                ivs_count=1,
                confidence="HIGH",
                intron_start=4,
                intron_end=6,
                intron_len=3,
                exon1_start=1,
                exon1_end=3,
                exon2_start=7,
                exon2_end=10,
            )

            write_fasta_outputs(args, [result], {"q1": "AAACCCGGGG"})

            free_text = args.intron_free_fa.read_text(encoding="utf-8")
            intron_text = args.introns_fa.read_text(encoding="utf-8")
            self.assertIn("action=removed", free_text)
            self.assertIn("AAAGGGG", free_text)
            self.assertIn(">q1|ivs_1|4-6|len=3|confidence=HIGH", intron_text)
            self.assertIn("CCC", intron_text)

    def test_write_fasta_outputs_removes_multiple_passing_ivs(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            args = SimpleNamespace(
                intron_free_fa=tmp / "ivs_free.fa",
                introns_fa=tmp / "ivs.fa",
                min_output_confidence="MEDIUM",
            )
            results = [
                QueryResult(
                    query_id="q1",
                    query_len=12,
                    classification="MEDIUM_CONFIDENCE_16S_IVS",
                    confidence="MEDIUM",
                    ivs_index=1,
                    ivs_count=2,
                    intron_start=4,
                    intron_end=5,
                    intron_len=2,
                    exon1_start=1,
                    exon1_end=3,
                    exon2_start=6,
                    exon2_end=12,
                ),
                QueryResult(
                    query_id="q1",
                    query_len=12,
                    classification="MEDIUM_CONFIDENCE_16S_IVS",
                    confidence="MEDIUM",
                    ivs_index=2,
                    ivs_count=2,
                    intron_start=9,
                    intron_end=10,
                    intron_len=2,
                    exon1_start=1,
                    exon1_end=8,
                    exon2_start=11,
                    exon2_end=12,
                ),
            ]

            write_fasta_outputs(args, results, {"q1": "AAACCCTTGGGG"})

            free_text = args.intron_free_fa.read_text(encoding="utf-8")
            ivs_text = args.introns_fa.read_text(encoding="utf-8")
            self.assertIn("ivs=4-5,9-10", free_text)
            self.assertIn("AAACTTGG", free_text)
            self.assertIn(">q1|ivs_1|4-5|len=2|confidence=MEDIUM", ivs_text)
            self.assertIn(">q1|ivs_2|9-10|len=2|confidence=MEDIUM", ivs_text)

    def test_write_bed_outputs_reports_final_exon_segments_for_multiple_ivs(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            args = SimpleNamespace(
                introns_bed=tmp / "q1.ivs.bed",
                exons_bed=tmp / "q1.exons.bed",
                min_output_confidence="MEDIUM",
            )
            results = [
                QueryResult(
                    query_id="q1",
                    query_len=12,
                    classification="MEDIUM_CONFIDENCE_16S_IVS",
                    confidence="MEDIUM",
                    ivs_index=1,
                    ivs_count=2,
                    intron_start=4,
                    intron_end=5,
                    intron_len=2,
                ),
                QueryResult(
                    query_id="q1",
                    query_len=12,
                    classification="MEDIUM_CONFIDENCE_16S_IVS",
                    confidence="MEDIUM",
                    ivs_index=2,
                    ivs_count=2,
                    intron_start=9,
                    intron_end=10,
                    intron_len=2,
                ),
            ]

            write_bed_outputs(args, results)

            ivs_bed = args.introns_bed.read_text(encoding="utf-8")
            exon_bed = args.exons_bed.read_text(encoding="utf-8")
            self.assertEqual(ivs_bed.count("\n"), 2)
            self.assertIn("q1\t3\t5\tq1|ivs_1|4-5|confidence=MEDIUM", ivs_bed)
            self.assertIn("q1\t8\t10\tq1|ivs_2|9-10|confidence=MEDIUM", ivs_bed)
            self.assertIn("q1\t0\t3\tq1|exon_1|1-3|ivs_removed=2|confidence=MEDIUM", exon_bed)
            self.assertIn("q1\t5\t8\tq1|exon_2|6-8|ivs_removed=2|confidence=MEDIUM", exon_bed)
            self.assertIn("q1\t10\t12\tq1|exon_3|11-12|ivs_removed=2|confidence=MEDIUM", exon_bed)

    def test_write_summary_includes_blast_diagnostic_columns(self) -> None:
        with TemporaryDirectory() as tmpdir:
            summary = Path(tmpdir) / "summary.tsv"
            result = QueryResult(
                query_id="q1",
                query_len=2417,
                classification="NO_IVS_SIGNAL",
                confidence="NONE",
                blast_status="BLAST_HIT_NO_IVS_PATTERN",
                blast_raw_hsps=5,
                blast_raw_subjects=2,
                blast_retained_hsps=3,
                blast_retained_subjects=1,
                blast_subjects_analyzed=1,
                best_blast_subject="s1",
                best_blast_pident=88.0,
                best_blast_bitscore=900.0,
                best_blast_taxonomy="Archaea;P;C;O;F;Genus;Genus species",
                reasons=["blast_hsps_present_but_no_supported_hsp_gap_pattern"],
            )

            write_summary(summary, [result], "genus")

            text = summary.read_text(encoding="utf-8")
            self.assertIn("blast_status", text)
            self.assertIn("BLAST_HIT_NO_IVS_PATTERN", text)
            self.assertIn("best_blast_subject", text)
            self.assertIn("ivs_index", text)
            self.assertIn("ivs_start", text)
            self.assertNotIn("intron_start", text)

    def test_print_summary_handles_blast_status_counts(self) -> None:
        args = SimpleNamespace(
            min_output_confidence="MEDIUM",
            summary_tsv="summary.tsv",
            intron_free_fa="intron_free.fa",
            introns_fa="introns.fa",
        )
        result = QueryResult(
            query_id="q1",
            query_len=2417,
            classification="NO_IVS_SIGNAL",
            confidence="NONE",
            blast_status="NO_BLAST_HIT",
        )

        with patch("ivsblastn.outputs.CONSOLE.print") as print_mock:
            print_summary([result], args)

        self.assertTrue(print_mock.called)


if __name__ == "__main__":
    unittest.main()
