import unittest

from ivsblastn.cli import build_parser, normalize_legacy_argv


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
        self.assertFalse(args.clean_ref_introns)


if __name__ == "__main__":
    unittest.main()
