from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import sys
import unittest

from ivsblastn.blast import parse_blast, run_external_command


class BlastParsingTests(unittest.TestCase):
    def test_parse_blast_reads_gzip_and_skips_malformed_rows(self) -> None:
        with TemporaryDirectory() as tmpdir:
            blast_path = Path(tmpdir) / "hits.tsv.gz"
            with gzip.open(blast_path, "wt", encoding="utf-8") as handle:
                handle.write("q1\ts1\t99.0\t120\t1\t120\t5\t124\t1e-50\t200\n")
                handle.write("q1\ts2\tbad\t120\t1\t120\t5\t124\t1e-50\t200\n")
                handle.write("q1\ts3\t99.0\t50\t1\t50\t5\t54\t1e-10\t80\n")

            grouped = parse_blast(blast_path, min_pident=75.0, min_hsp_len=100)

            self.assertIn("q1", grouped)
            self.assertIn("s1", grouped["q1"])
            self.assertNotIn("s2", grouped["q1"])
            self.assertNotIn("s3", grouped["q1"])

    def test_external_command_failure_reports_stderr(self) -> None:
        cmd = [
            sys.executable,
            "-c",
            "import sys; sys.stderr.write('blast database missing\\n'); raise SystemExit(2)",
        ]

        with self.assertRaisesRegex(RuntimeError, "blast database missing"):
            run_external_command(cmd, "BLASTN")


if __name__ == "__main__":
    unittest.main()
