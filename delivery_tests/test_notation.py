"""Check that the new notation is reversible and cannot hide changed witnesses."""

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from notation import current_catalogue, legacy_catalogue, resource_parameters


class NotationTests(unittest.TestCase):
    def setUp(self):
        self.old = json.loads((ROOT / "provenance/pareto_frontier_before_notation.json").read_bytes())
        self.new = current_catalogue(self.old)

    def test_exact_round_trip(self):
        self.assertEqual(legacy_catalogue(self.new), self.old)
        self.assertEqual(current_catalogue(self.new), self.new)

    def test_current_file(self):
        self.assertEqual(self.new, json.loads((ROOT / "data/protocols/pareto_frontier.json").read_bytes()))

    def test_resource_identity(self):
        for row in self.new["protocols"]:
            params = resource_parameters(row["generator_matrix_rows"], row["k"])
            self.assertTrue(all(row[key] == value for key, value in params.items()))

    def test_reject_inconsistent_rows(self):
        self.new["protocols"][0]["r"] += 1
        with self.assertRaises(ValueError):
            legacy_catalogue(self.new)

    def test_reject_conflicting_legacy_key(self):
        self.new["protocols"][0]["S"] = 99
        with self.assertRaises(ValueError):
            legacy_catalogue(self.new)

    def test_gate_s_unchanged(self):
        for old, new in zip(self.old["outputs"], self.new["outputs"]):
            self.assertEqual(old["representative_gate"], new["representative_gate"])
            self.assertEqual(old["factorisation_certificate"], new["factorisation_certificate"])

    def test_matrices_and_ids_unchanged(self):
        for old, new in zip(self.old["protocols"], self.new["protocols"]):
            for field in ("generator_matrix_rows", "output_id", "n", "d_Z", "error_coefficient"):
                self.assertEqual(old[field], new[field])


if __name__ == "__main__":
    unittest.main()
