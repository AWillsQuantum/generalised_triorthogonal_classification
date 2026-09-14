import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_space_m08_cover import check_primary_partition, check_quotient


class SpaceCoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.domain = json.loads((ROOT / "data/space_contractions/c54_m08/domain.json").read_bytes())
        cls.census = json.loads((ROOT / "data/space_contractions/c54_m08/primary_census.json").read_bytes())
        cls.alternative = json.loads((ROOT / "certificates/alternative_contractions.json").read_bytes())

    def test_complete_primary_partition_includes_empty_family(self):
        result = check_primary_partition(self.domain, self.census, self.alternative)
        self.assertEqual(result["complete_raw_marked_pairs"], 154850689024)
        self.assertEqual(result["primary_cores"], 317)
        self.assertEqual(sum(r["marked_representatives"] == 0 for r in self.census["primary_cores"]), 1)

    def test_deleted_repeated_or_misbound_core_rejected(self):
        for mode in ("delete", "repeat", "misbind", "mass"):
            bad = copy.deepcopy(self.census)
            rows = bad["primary_cores"]
            if mode == "delete":
                rows.pop()
            elif mode == "repeat":
                rows[-1] = copy.deepcopy(rows[0])
            elif mode == "misbind":
                rows[0]["input_sha256"] = "0"*64
            else:
                rows[0]["raw_marked_pairs"] += 1
            with self.assertRaises(ValueError):
                check_primary_partition(self.domain, bad, self.alternative)

    def test_exact_quotient_accounting(self):
        result = check_quotient(self.census)
        self.assertEqual(result["affine_classes"], 17777766)
        self.assertEqual(result["quotient_slices"], 5262)

    def test_quotient_gaps_or_incomplete_negative_tests_rejected(self):
        for field, value in (("first_candidate", 1), ("first_bucket", 1),
                             ("first_representative", 1), ("exact_transporter_tests_complete", False),
                             ("equivalent_comparisons", 0)):
            bad = copy.deepcopy(self.census)
            bad["exact_quotient_slices"][0][field] = value
            with self.assertRaises(ValueError):
                check_quotient(bad)


if __name__ == "__main__":
    unittest.main()
