import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_cubic_parent_sector import check_raw_census


class CubicParentCensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = json.loads((ROOT / "data/protocol_sectors/length52_separate_cubic_censuses.json").read_bytes())["cases"]

    def test_complete_positive_profile(self):
        self.assertEqual(check_raw_census(self.cases[0]["censuses"][0], 17), {(0xe400,)})

    def test_complete_hyperplane_zero_census(self):
        self.assertEqual(check_raw_census(self.cases[58]["censuses"][1], 15), set())

    def test_reject_lost_radical_mass(self):
        row = copy.deepcopy(self.cases[58]["censuses"][1])
        row["statistics"]["radical_dimension_counts"]["5"] -= 1
        with self.assertRaisesRegex(ValueError, "mass"):
            check_raw_census(row, 15)

    def test_reject_truncated_output_profile(self):
        row = copy.deepcopy(self.cases[0]["censuses"][0])
        row["output_keys"] = []
        with self.assertRaisesRegex(ValueError, "profile"):
            check_raw_census(row, 17)

    def test_reject_filtered_input(self):
        row = copy.deepcopy(self.cases[0]["censuses"][0])
        row["statistics"]["quotient_dimension_filtered_supports"] = 1
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            check_raw_census(row, 17)


if __name__ == "__main__":
    unittest.main()
