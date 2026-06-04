from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from ivsblastn.fasta import species_key_from_taxonomy
from ivsblastn.reference import preprocess_reference


class ReferenceSelectionTests(unittest.TestCase):
    def test_species_key_requires_clear_binomial_without_digits(self) -> None:
        self.assertEqual(
            species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio halioticoli"),
            "Bacteria;P;C;O;F;Vibrio;Vibrio halioticoli",
        )
        self.assertIsNone(species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio"))
        self.assertIsNone(species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio 1234"))
        self.assertIsNone(species_key_from_taxonomy("Bacteria;P;C;O;F;Vibrio;Vibrio sp001"))

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


if __name__ == "__main__":
    unittest.main()
