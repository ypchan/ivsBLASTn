import unittest

from ivsblastn.cli import build_parser


class CliDefaultTests(unittest.TestCase):
    def test_query_blast_limits_are_top_subjects_times_hsps(self) -> None:
        args = build_parser().parse_args(
            [
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


if __name__ == "__main__":
    unittest.main()
