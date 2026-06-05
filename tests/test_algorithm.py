from types import SimpleNamespace
import unittest

from ivsblastn.algorithm import analyze_query, analyze_query_all
from ivsblastn.models import BlastQueryStats, HSP


def algorithm_args() -> SimpleNamespace:
    return SimpleNamespace(
        algorithm="hsp-gap-support",
        top_subjects=100,
        min_pident=70.0,
        min_hsp_len=100,
        min_intron_len=25,
        max_intron_len=2000,
        max_ref_gap=15,
        max_query_overlap=20,
        breakpoint_window=20,
        tax_rank="genus",
        min_support_subjects=1,
        medium_support_subjects=3,
        medium_support_taxa=3,
        high_support_subjects=10,
        high_support_taxa=3,
    )


class AlgorithmStatusTests(unittest.TestCase):
    def test_no_raw_blast_rows_are_reported_as_no_blast_hit(self) -> None:
        result = analyze_query(
            "q1",
            {},
            2417,
            {},
            algorithm_args(),
            BlastQueryStats(raw_hsps=0, raw_subjects=0, retained_hsps=0, retained_subjects=0),
        )

        self.assertEqual(result.blast_status, "NO_BLAST_HIT")
        self.assertIn("no_blast_hsp_reported", result.reasons)

    def test_filtered_blast_rows_are_reported_separately(self) -> None:
        result = analyze_query(
            "q1",
            {},
            2417,
            {},
            algorithm_args(),
            BlastQueryStats(raw_hsps=4, raw_subjects=2, retained_hsps=0, retained_subjects=0),
        )

        self.assertEqual(result.blast_status, "BLAST_HITS_FILTERED")
        self.assertIn("blast_hsps_failed_filters", result.reasons)
        self.assertEqual(result.blast_raw_subjects, 2)

    def test_retained_hit_without_gap_is_reported_as_no_ivs_pattern(self) -> None:
        hsp = HSP(
            qseqid="q1",
            sseqid="s1",
            pident=88.0,
            length=500,
            qstart=1,
            qend=500,
            sstart=1,
            send=500,
            evalue="1e-50",
            bitscore=900.0,
        )
        result = analyze_query(
            "q1",
            {"s1": [hsp]},
            2417,
            {"s1": "Archaea;P;C;O;F;Genus;Genus species"},
            algorithm_args(),
            BlastQueryStats(raw_hsps=1, raw_subjects=1, retained_hsps=1, retained_subjects=1),
        )

        self.assertEqual(result.blast_status, "BLAST_HIT_NO_IVS_PATTERN")
        self.assertEqual(result.best_blast_subject, "s1")
        self.assertEqual(result.best_blast_pident, 88.0)
        self.assertIn("blast_hsps_present_but_no_supported_ivs_gap_pattern", result.reasons)

    def test_two_non_overlapping_ivs_are_reported_for_one_query(self) -> None:
        hsps = [
            HSP("q1", "s1", 99.0, 100, 1, 100, 1, 100, "1e-30", 200.0),
            HSP("q1", "s1", 99.0, 100, 151, 250, 101, 200, "1e-30", 200.0),
            HSP("q1", "s1", 99.0, 100, 251, 350, 201, 300, "1e-30", 200.0),
            HSP("q1", "s1", 99.0, 100, 401, 500, 301, 400, "1e-30", 200.0),
            HSP("q1", "s2", 98.0, 100, 1, 100, 1, 100, "1e-30", 190.0),
            HSP("q1", "s2", 98.0, 100, 151, 250, 101, 200, "1e-30", 190.0),
            HSP("q1", "s2", 98.0, 100, 251, 350, 201, 300, "1e-30", 190.0),
            HSP("q1", "s2", 98.0, 100, 401, 500, 301, 400, "1e-30", 190.0),
        ]
        results = analyze_query_all(
            "q1",
            {"s1": hsps[:4], "s2": hsps[4:]},
            600,
            {
                "s1": "Archaea;P;C;O;F;GenusA;GenusA species",
                "s2": "Archaea;P;C;O;F;GenusB;GenusB species",
            },
            algorithm_args(),
            BlastQueryStats(raw_hsps=8, raw_subjects=2, retained_hsps=8, retained_subjects=2),
        )

        self.assertEqual(len(results), 2)
        self.assertEqual([(r.intron_start, r.intron_end) for r in results], [(101, 150), (351, 400)])
        self.assertEqual([r.ivs_index for r in results], [1, 2])
        self.assertEqual([r.ivs_count for r in results], [2, 2])
        self.assertEqual([r.intron_free_len for r in results], [500, 500])
        self.assertTrue(all(r.classification == "LOW_CONFIDENCE_16S_IVS" for r in results))


if __name__ == "__main__":
    unittest.main()
