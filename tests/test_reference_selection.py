from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import unittest

from ivsblastn.fasta import genus_key_from_taxonomy, read_fasta, species_key_from_taxonomy
from ivsblastn.cli import init_reference_command, setup_reference_output_paths
from ivsblastn.models import QueryResult
from ivsblastn.reference import preprocess_reference, write_reference_introns_fasta


class ReferenceSelectionTests(unittest.TestCase):
    def test_species_key_requires_clear_binomial_without_digits(self) -> None:
        self.assertEqual(
            species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio halioticoli"),
            "Bacteria;P;C;O;F;Vibrio;Vibrio halioticoli",
        )
        self.assertIsNone(species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio"))
        self.assertIsNone(species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio 1234"))
        self.assertIsNone(species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio sp001"))
        self.assertIsNone(species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio sp."))
        self.assertIsNone(species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio cf."))

    def test_genus_key_uses_available_genus_for_unclear_species(self) -> None:
        self.assertEqual(
            genus_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio sp001"),
            "Bacteria;P;C;O;F;Vibrio",
        )

    def test_read_fasta_rejects_duplicate_ids(self) -> None:
        with TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "query.fa"
            fasta.write_text(">q1\nAAAA\n>q1 duplicate\nCCCC\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Duplicate FASTA sequence ID"):
                read_fasta(fasta)

    def test_preprocess_reference_keeps_longest_record_per_species(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ref_fasta = tmp / "reference.fa"
            raw_ref_fa = tmp / "raw_reference.fa"
            raw_ref_tax = tmp / "raw_reference.tax.tsv"
            raw_ref_db = tmp / "raw_reference_db"
            ref_fasta.write_text(
                "\n".join(
                    [
                        ">short Bacteria;P;C;O;F;Vibrio;Vibrio halioticoli",
                        "ATGC",
                        ">long Bacteria;P;C;O;F;Vibrio;Vibrio halioticoli",
                        "ATGCATGC",
                        ">unclear Bacteria;P;C;O;F;Vibrio;Vibrio sp001",
                        "ATGCATGCATGC",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            args = SimpleNamespace(
                ref_fasta=ref_fasta,
                ref_domains="Archaea,Bacteria",
                ref_per_species=1,
                ref_unclear_per_genus=0,
                raw_ref_fa=raw_ref_fa,
                raw_ref_tax=raw_ref_tax,
                raw_ref_db=raw_ref_db,
                makeblastdb_bin="true",
            )

            preprocess_reference(args)

            selected_fasta = raw_ref_fa.read_text(encoding="utf-8")
            self.assertIn(">long", selected_fasta)
            self.assertNotIn(">short", selected_fasta)
            self.assertNotIn(">unclear", selected_fasta)

    def test_preprocess_reference_keeps_unclear_species_by_genus_fallback(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ref_fasta = tmp / "reference.fa"
            raw_ref_fa = tmp / "raw_reference.fa"
            raw_ref_tax = tmp / "raw_reference.tax.tsv"
            raw_ref_db = tmp / "raw_reference_db"
            ref_fasta.write_text(
                "\n".join(
                    [
                        ">short_unclear Bacteria;P;C;O;F;Vibrio;Vibrio sp001",
                        "ATGC",
                        ">long_unclear Bacteria;P;C;O;F;Vibrio;Vibrio sp002",
                        "ATGCATGC",
                        ">other_genus Bacteria;P;C;O;F;Photobacterium;Photobacterium",
                        "ATGCATGCATGC",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            args = SimpleNamespace(
                ref_fasta=ref_fasta,
                ref_domains="Archaea,Bacteria",
                ref_per_species=1,
                ref_unclear_per_genus=1,
                raw_ref_fa=raw_ref_fa,
                raw_ref_tax=raw_ref_tax,
                raw_ref_db=raw_ref_db,
                makeblastdb_bin="true",
            )

            preprocess_reference(args)

            selected_fasta = raw_ref_fa.read_text(encoding="utf-8")
            self.assertIn(">long_unclear", selected_fasta)
            self.assertIn(">other_genus", selected_fasta)
            self.assertNotIn(">short_unclear", selected_fasta)

    def test_init_reference_command_writes_reusable_manifest(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ref_fasta = tmp / "reference.fa"
            outdir = tmp / "prepared"
            ref_fasta.write_text(
                ">ref1 Bacteria;P;C;O;F;Vibrio;Vibrio halioticoli\nATGCATGC\n",
                encoding="utf-8",
            )
            args = SimpleNamespace(
                ref_fasta=ref_fasta,
                outdir=outdir,
                ref_domains="Archaea,Bacteria",
                ref_per_species=1,
                ref_unclear_per_genus=5,
                clean_ref_introns=False,
                min_intron_len=25,
                max_intron_len=2000,
                makeblastdb_bin="true",
            )

            with patch("ivsblastn.cli.CONSOLE.print"):
                self.assertEqual(init_reference_command(args), 0)

            manifest = (outdir / "reference_manifest.tsv").read_text(encoding="utf-8")
            self.assertIn("blast_db_prefix", manifest)
            self.assertTrue((outdir / "raw_reference.fa").exists())
            self.assertTrue((outdir / "raw_reference.tax.tsv").exists())

    def test_init_reference_manifest_includes_self_clean_outputs(self) -> None:
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            ref_fasta = tmp / "reference.fa"
            outdir = tmp / "prepared"
            ref_fasta.write_text(">ref1 Bacteria;P;C;O;F;Vibrio;Vibrio halioticoli\nATGCATGC\n", encoding="utf-8")
            raw_ref_fa = outdir / "raw_reference.fa"
            raw_ref_tax = outdir / "raw_reference.tax.tsv"
            raw_ref_db = outdir / "raw_reference_db"
            cleaned_ref_fa = outdir / "cleaned_reference.fa"
            cleaned_ref_tax = outdir / "cleaned_reference.tax.tsv"
            cleaned_ref_db = outdir / "cleaned_reference_db"
            args = SimpleNamespace(
                ref_fasta=ref_fasta,
                outdir=outdir,
                clean_ref_introns=True,
                min_intron_len=25,
                max_intron_len=2000,
                ref_self_blast_max_hsps=20,
            )

            with patch("ivsblastn.cli.preprocess_reference", return_value=(raw_ref_fa, raw_ref_tax, raw_ref_db)), patch(
                "ivsblastn.cli.clean_reference_introns",
                return_value=(cleaned_ref_fa, cleaned_ref_tax, cleaned_ref_db),
            ), patch("ivsblastn.cli.CONSOLE.print"):
                self.assertEqual(init_reference_command(args), 0)

            manifest = (outdir / "reference_manifest.tsv").read_text(encoding="utf-8")
            self.assertIn("reference_self_clean_report", manifest)
            self.assertIn("reference_self_clean_summary", manifest)
            self.assertIn("reference_introns_fasta", manifest)

    def test_reference_output_paths_add_self_clean_report_compatibility_fields(self) -> None:
        with TemporaryDirectory() as tmpdir:
            args = SimpleNamespace(
                ref_fasta=Path(tmpdir) / "reference.fa",
                outdir=Path(tmpdir) / "prepared",
                ref_self_blast_max_hsps=20,
                ref_clean_min_confidence="MEDIUM",
            )

            setup_reference_output_paths(args)

            self.assertEqual(args.blast_max_hsps, 20)
            self.assertEqual(args.min_output_confidence, "MEDIUM")
            self.assertEqual(args.ref_self_clean_introns_fa, Path(tmpdir) / "prepared" / "results" / "reference_self_clean.introns.fa")

    def test_write_reference_introns_fasta_outputs_removed_reference_ivs(self) -> None:
        with TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "reference_self_clean.introns.fa"
            result = QueryResult(
                query_id="ref1",
                query_len=10,
                classification="HIGH_CONFIDENCE_16S_INTRON",
                confidence="HIGH",
                intron_start=4,
                intron_end=6,
                intron_len=3,
                exon1_start=1,
                exon1_end=3,
                exon2_start=7,
                exon2_end=10,
            )

            written = write_reference_introns_fasta(output, {"ref1": "AAACCCGGGG"}, [result], "LOW")

            text = output.read_text(encoding="utf-8")
            self.assertEqual(written, 1)
            self.assertIn(">ref1|reference_intron|4-6|len=3|confidence=HIGH|action=removed", text)
            self.assertIn("CCC", text)


if __name__ == "__main__":
    unittest.main()
