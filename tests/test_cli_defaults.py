import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ivsblastn.cli import build_parser, normalize_legacy_argv, slurm_run_args, validate_run_args


class CliDefaultTests(unittest.TestCase):
    def test_query_blast_limits_are_top_subjects_times_hsps(self) -> None:
        args = build_parser().parse_args(
            [
                "run",
                "--query",
                "query.fa",
                "--db",
                "reference_db",
                "--outdir",
                "out",
            ]
        )

        self.assertEqual(args.top_subjects, 100)
        self.assertEqual(args.blast_max_hsps, 5)
        self.assertEqual(args.ref_unclear_per_genus, 5)
        self.assertFalse(hasattr(args, "blast_max_target_seqs"))

    def test_legacy_run_arguments_are_mapped_to_run_subcommand(self) -> None:
        self.assertEqual(normalize_legacy_argv(["--query", "query.fa"]), ["run", "--query", "query.fa"])
        self.assertEqual(normalize_legacy_argv(["run", "--query", "query.fa"]), ["run", "--query", "query.fa"])

    def test_init_reference_subcommand_exists(self) -> None:
        args = build_parser().parse_args(
            [
                "init-reference",
                "--ref-fasta",
                "silva.fa",
                "--outdir",
                "reference",
            ]
        )

        self.assertEqual(args.ref_per_species, 1)
        self.assertEqual(args.ref_unclear_per_genus, 5)
        self.assertEqual(args.ref_clean_min_confidence, "MEDIUM")
        self.assertFalse(args.clean_ref_introns)

    def test_ivs_named_options_keep_legacy_destinations(self) -> None:
        args = build_parser().parse_args(
            [
                "run",
                "--query",
                "query.fa",
                "--db",
                "reference_db",
                "--outdir",
                "out",
                "--clean-ref-ivs",
                "--min-ivs-len",
                "30",
                "--max-ivs-len",
                "1500",
            ]
        )

        self.assertTrue(args.clean_ref_introns)
        self.assertEqual(args.min_intron_len, 30)
        self.assertEqual(args.max_intron_len, 1500)

    def test_submit_slurm_forwards_detection_defaults(self) -> None:
        args = build_parser().parse_args(
            [
                "submit-slurm",
                "--chunks-dir",
                "chunks",
                "--db",
                "reference_db",
                "--outdir",
                "batch",
            ]
        )

        self.assertEqual(args.min_pident, 70.0)
        self.assertEqual(args.min_hsp_len, 100)
        self.assertEqual(args.max_ref_gap, 15)
        self.assertEqual(args.breakpoint_window, 20)
        self.assertEqual(args.min_output_confidence, "MEDIUM")
        self.assertFalse(args.gzip_fasta_output)

        forwarded = "\n".join(slurm_run_args(args))
        self.assertIn("--min-ivs-len 25", forwarded)
        self.assertIn("--max-ivs-len 2000", forwarded)
        self.assertIn("--min-output-confidence MEDIUM", forwarded)

    def test_run_validation_reports_missing_blast_db_prefix(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query = tmp / "query.fa"
            taxonomy = tmp / "reference.tax.tsv"
            query.write_text(">q1\nACGT\n", encoding="utf-8")
            taxonomy.write_text("s1\tArchaea;Genus species\n", encoding="utf-8")

            args = build_parser().parse_args(
                [
                    "run",
                    "--query",
                    str(query),
                    "--db",
                    str(tmp / "missing_db"),
                    "--taxonomy",
                    str(taxonomy),
                    "--outdir",
                    str(tmp / "out"),
                ]
            )

            with self.assertRaisesRegex(FileNotFoundError, "BLAST DB files not found"):
                validate_run_args(args)


if __name__ == "__main__":
    unittest.main()
