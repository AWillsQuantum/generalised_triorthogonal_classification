import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_length52_large_quotients import check_class


class LargeQuotientBindingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.domain = json.loads((ROOT / "data/protocol_sectors/length52_large_quotient_domain.json").read_bytes())
        cls.primitive = json.loads((ROOT / cls.domain["primitive_domain"]).read_bytes())
        cls.profiles = json.loads((ROOT / cls.domain["profile_dataset"]).read_bytes())["profiles"]

    def test_positive_distance_four_quotient_below_logical_dimension(self):
        result = check_class(self.domain["classes"][385], 385, self.primitive, self.profiles[385])
        self.assertEqual(result, dict(index=385, d3=15, d4=1, q5_positive=True))

    def test_additional_primitive_zero_result(self):
        result = check_class(self.domain["classes"][55], 55, self.primitive, self.profiles[55])
        self.assertEqual(result["d3"], 15)

    def test_reject_profile_assigned_to_another_geometry(self):
        with self.assertRaisesRegex(ValueError, "profile is not bound"):
            check_class(self.domain["classes"][55], 55, self.primitive, self.profiles[56])

    def test_reject_nonzero_primitive_result(self):
        row = copy.deepcopy(self.domain["classes"][55])
        row["primitive_target_census"]["weighted_subspace_count"] = 1
        with self.assertRaisesRegex(ValueError, "zero result"):
            check_class(row, 55, self.primitive, self.profiles[55])

    def test_reject_ambiguous_primitive_binding(self):
        row = copy.deepcopy(self.domain["classes"][55])
        row["primitive_union_case_index"] = 0
        with self.assertRaisesRegex(ValueError, "exactly one"):
            check_class(row, 55, self.primitive, self.profiles[55])

    def test_reject_profile_mass_contradiction(self):
        profile = copy.deepcopy(self.profiles[385])
        profile["nondegenerate_subspaces"] = 0
        with self.assertRaisesRegex(ValueError, "profile is not bound"):
            check_class(self.domain["classes"][385], 385, self.primitive, profile)


if __name__ == "__main__":
    unittest.main()
