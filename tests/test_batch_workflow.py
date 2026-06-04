from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import unittest

from ivsblastn.merge import merge_chunk_outputs
from ivsblastn.slurm import render_slurm_array_script
from ivsblastn.split import split_fasta


class BatchWorkflowTests(unittest.TestCase):
    def test_split_fasta_writes_chunks_and_manifest(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query = tmp / "query.fa"
            query.write_text(">q1 desc\nAAAA\n>q2\nCCCC\n>q3\nGGGG\n", encoding="utf-8")

            chunks = split_fasta(query, tmp / "chunks", chunk_size=2, prefix="part")

            self.assertEqual(len(chunks), 2)
            self.assertTrue((tmp / "chunks" / "chunks.tsv").exists())
            self.assertIn(">q1 desc", chunks[0].read_text(encoding="utf-8"))
            self.assertIn(">q3", chunks[1].read_text(encoding="utf-8"))

    def test_merge_keeps_one_tsv_header(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            for chunk_name, query_id in [("query.000001", "q1"), ("query.000002", "q2")]:
                results = tmp / "chunks" / chunk_name / "results"
                results.mkdir(parents=True)
                (results / f"{chunk_name}.summary.tsv").write_text("query_id\tconfidence\n" f"{query_id}\tNONE\n", encoding="utf-8")

            merge_chunk_outputs(tmp / "chunks", tmp / "final", label="all")
            merged = (tmp / "final" / "results" / "all.summary.tsv").read_text(encoding="utf-8")

            self.assertEqual(merged.count("query_id\tconfidence"), 1)
            self.assertIn("q1\tNONE", merged)
            self.assertIn("q2\tNONE", merged)

    def test_merge_concatenates_gzip_fasta_outputs(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            for chunk_name, query_id in [("query.000001", "q1"), ("query.000002", "q2")]:
                results = tmp / "chunks" / chunk_name / "results"
                results.mkdir(parents=True)
                with gzip.open(results / f"{chunk_name}.intron_free.fa.gz", "wt", encoding="utf-8") as handle:
                    handle.write(f">{query_id}\nAAAA\n")

            merge_chunk_outputs(tmp / "chunks", tmp / "final", label="all")
            merged_path = tmp / "final" / "results" / "all.intron_free.fa.gz"

            with gzip.open(merged_path, "rt", encoding="utf-8") as handle:
                merged = handle.read()
            self.assertIn(">q1", merged)
            self.assertIn(">q2", merged)

    def test_slurm_script_uses_manifest_chunk_count(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            chunks = tmp / "chunks"
            chunks.mkdir()
            (chunks / "chunks.tsv").write_text("chunk_index\tquery_fasta\n1\t/a.fa\n2\t/b.fa\n", encoding="utf-8")

            script = render_slurm_array_script(
                chunks_dir=chunks,
                outdir=tmp / "run",
                db=Path("ref_db"),
                taxonomy=Path("ref.tax.tsv"),
                threads=8,
                top_subjects=100,
                blast_max_hsps=5,
                cpus_per_task=8,
                mem="16G",
                time="12:00:00",
                partition=None,
                account=None,
                qos=None,
                array_concurrency=10,
                extra_run_args="--min-pident 80",
                run_args=["--max-ref-gap 12", "--min-output-confidence MEDIUM"],
            )

            self.assertIn("#SBATCH --array=1-2%10", script)
            self.assertIn("ivsBLASTn run", script)
            self.assertIn("--min-pident 80", script)
            self.assertIn("--max-ref-gap 12", script)
            self.assertIn("--min-output-confidence MEDIUM", script)
            self.assertIn("awk -F '\\t'", script)

    def test_slurm_script_supports_commercial_resource_directives(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            chunks = tmp / "chunks"
            chunks.mkdir()
            (chunks / "chunks.tsv").write_text("chunk_index\tquery_fasta\n1\t/a.fa\n", encoding="utf-8")

            script = render_slurm_array_script(
                chunks_dir=chunks,
                outdir=tmp / "run",
                db=Path("ref_db"),
                taxonomy=None,
                threads=20,
                top_subjects=100,
                blast_max_hsps=5,
                cpus_per_task=20,
                mem="32G",
                time="24:00:00",
                partition="normal_fcp1",
                account="prj_219_3",
                qos="qos_prj_219_3",
                nodes=1,
                ntasks=1,
                constraint="cpu",
                gres=None,
                exclude="node001",
                nodelist=None,
                extra_sbatch_options=["--mail-type=END", "#SBATCH --comment=ivsBLASTn"],
                array_concurrency=None,
                extra_run_args="",
            )

            self.assertIn("#SBATCH --partition=normal_fcp1", script)
            self.assertIn("#SBATCH --qos=qos_prj_219_3", script)
            self.assertIn("#SBATCH --account=prj_219_3", script)
            self.assertIn("#SBATCH --cpus-per-task=20", script)
            self.assertIn("#SBATCH --nodes=1", script)
            self.assertIn("#SBATCH --ntasks=1", script)
            self.assertIn("#SBATCH --constraint=cpu", script)
            self.assertIn("#SBATCH --exclude=node001", script)
            self.assertIn("#SBATCH --mail-type=END", script)
            self.assertIn("#SBATCH --comment=ivsBLASTn", script)


if __name__ == "__main__":
    unittest.main()
