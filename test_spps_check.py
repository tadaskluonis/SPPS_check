import io
import json
import tempfile
import unittest
from pathlib import Path

import spps_check


class AnalysisTests(unittest.TestCase):
    def test_complete_record_and_ranked_recommendations(self):
        result = spps_check.analyze_sequence("ACVVVVDGCMW", "hard")

        self.assertEqual(result["schema_version"], "1.0")
        self.assertEqual(result["features"]["length"], 11)
        self.assertEqual(result["features"]["max_beta_branched_run"], 4)
        self.assertEqual(result["features"]["aspartimide_motif_count"], 1)
        self.assertEqual(result["features"]["cysteine_count"], 2)
        self.assertEqual(result["sites"]["aspartimide_motifs"]["DG"], [7])
        self.assertLessEqual(len(result["recommendations"]), 3)
        self.assertEqual(result["risk_level"], "high")
        self.assertEqual(set(result["features"]["residue_counts"]), set(spps_check.AA_ORDER))

    def test_invalid_residue_is_not_silently_removed(self):
        with self.assertRaisesRegex(ValueError, "unsupported residue"):
            spps_check.analyze_sequence("ACDX")

    def test_low_risk_sequence_has_default_recommendation(self):
        result = spps_check.analyze_sequence("AEKPGQ")
        self.assertEqual(result["risk_level"], "none")
        self.assertEqual(len(result["recommendations"]), 1)

    def test_fasta_batch_and_jsonl(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "batch.fasta"
            path.write_text(">one\nACDE\n>two\nVVVVV\n")
            records = spps_check.collect_records([str(path)])

        output = io.StringIO()
        spps_check.render_records(records, "jsonl", output)
        decoded = [json.loads(line) for line in output.getvalue().splitlines()]
        self.assertEqual([record["id"] for record in decoded], ["one", "two"])

    def test_directory_kind_filters_mixed_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "same.fasta").write_text(">same\nACDE\n")
            (root / "same.pdb").write_text(
                "ATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00  0.00           C  \n"
            )
            records = spps_check.collect_records([directory], kind="structure")

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["chain"], "A")

    def test_pdb_chains_use_ca_residues_once(self):
        pdb = (
            "ATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00  0.00           C  \n"
            "ATOM      2  CA  CYS A   2      0.000   0.000   0.000  1.00  0.00           C  \n"
            "ATOM      3  CA  GLY B   1      0.000   0.000   0.000  1.00  0.00           C  \n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "small.pdb"
            path.write_text(pdb)
            records = spps_check.collect_records([str(path)], selected_chain="A")

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["sequence"], "AC")
        self.assertEqual(records[0]["chain"], "A")

    def test_common_mmcif_atom_loop(self):
        cif = (
            "data_small\n"
            "loop_\n"
            "_atom_site.group_PDB\n"
            "_atom_site.id\n"
            "_atom_site.label_atom_id\n"
            "_atom_site.label_alt_id\n"
            "_atom_site.label_comp_id\n"
            "_atom_site.label_asym_id\n"
            "_atom_site.auth_seq_id\n"
            "_atom_site.pdbx_PDB_model_num\n"
            "ATOM 1 CA . ALA A 1 1\n"
            "ATOM 2 CA . GLY A 2 1\n"
            "#\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "small.cif"
            path.write_text(cif)
            records = spps_check.collect_records([str(path)])

        self.assertEqual(records[0]["sequence"], "AG")
        self.assertEqual(records[0]["structure"]["format"], "mmcif")

    def test_legacy_structure_path_argument(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sheet.pdb"
            path.write_text("SHEET    1   A 1 THR A   1  VAL A   2  0\n")
            record = spps_check.analyze_sequence("TV", "legacy", str(path))

        self.assertEqual(record["structure"]["beta_sheet_records"], 1)

    def test_pipeline_failure_threshold(self):
        high = spps_check.analyze_sequence("VVVVV")
        low = spps_check.analyze_sequence("AEKPGQ")
        self.assertTrue(spps_check._should_fail([high], "high"))
        self.assertFalse(spps_check._should_fail([low], "any"))


if __name__ == "__main__":
    unittest.main()
