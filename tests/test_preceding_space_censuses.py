import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_preceding_space_censuses import check_core_partition, check_generation, check_quotient


class PrecedingSpaceCensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = {(c, m): json.loads((ROOT / f"data/space_contractions/preceding_sectors/c{c}_m{m:02d}.json").read_bytes())
                    for c, last in ((50, 14), (52, 16)) for m in range(8, last+1)}

    def test_all_aggregate_censuses_close(self):
        for key, data in self.data.items():
            with self.subTest(sector=key):
                self.assertEqual(check_quotient(data)["affine_classes"], data["class_count"])
                check_generation(data)

    def test_quotient_gap_rejected(self):
        bad = copy.deepcopy(self.data[52, 9])
        bad["exact_quotient_intervals"][1]["first_bucket"] += 1
        with self.assertRaises(ValueError):
            check_quotient(bad)

    def test_partition_mass_and_class_counts_are_independent(self):
        for field in ("candidate_count", "class_count"):
            bad = copy.deepcopy(self.data[52, 8])
            bad[field] += 1
            with self.assertRaises(ValueError):
                check_quotient(bad)

    def test_missing_negative_verification_rejected(self):
        bad = copy.deepcopy(self.data[52, 10])
        bad["independent_negative_replay"]["negative_comparison_count"] -= 1
        with self.assertRaises(ValueError):
            check_generation(bad)

    def test_filtered_domain_cannot_drop_source(self):
        bad = copy.deepcopy(self.data[52, 9])
        del bad["generation"]["source_counts_by_multiplicity"]["2"]
        with self.assertRaises(ValueError):
            check_generation(bad)

    def test_alternative_core_cannot_be_omitted(self):
        data = dict(c=50, primary_cores=[], alternative_cores=[])
        with self.assertRaises(ValueError):
            check_core_partition(data, {(48, 7, 96): {}}, {})

    def test_empty_predecessor_cannot_generate_children(self):
        bad = copy.deepcopy(self.data[52, 16])
        bad["candidate_count"] = bad["generation"]["candidate_count"] = 1
        with self.assertRaises(ValueError):
            check_generation(bad)


if __name__ == "__main__":
    unittest.main()
