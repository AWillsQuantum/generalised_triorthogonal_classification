import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_length52_small_blocks import check_branch, check_quotient_cover, check_raw_pointings


class Length52BlockCoverTests(unittest.TestCase):
    def branch(self):
        return dict(q=5, minimum_distance=3, enumeration_mode="raw_isotropic_subspaces",
            quotient_dimension_interval=[5, 14], support_range=dict(start=0, count=3, end_exclusive=3),
            input_quotient_profile=[[5, 2], [15, 1]], positive_supports=[], witnesses=[],
            statistics=dict(source_spaces=1, supports_processed=3, eligible_supports=2,
                raw_enumerated_supports=2, marked_orbit_enumerated_supports=0,
                quotient_dimension_filtered_supports=1, isotropic_subspaces=9, nondegenerate_subspaces=0,
                isotropic_subspace_statistics_exact=True, radical_statistics_exact=True,
                radical_dimension_counts={"5": 9}, marked_code_orbits=0, canonical_output_orbits=0))

    def test_raw_pointing_translation_period(self):
        self.assertEqual(check_raw_pointings(7, 1, 181, [128, 52, 0, 1, 0]), 0)
        self.assertEqual(check_raw_pointings(7, 1, 91, [64, 26, 0, 1, 0]), 1)

    def test_invalid_pointing_count(self):
        with self.assertRaises(ValueError):
            check_raw_pointings(7, 1, 180, [127, 52, 0, 1, 0])

    def test_complete_small_quotient_with_deferred_input(self):
        branch = self.branch()
        check_branch(branch, {0: 9, 5: 2, 15: 1}, 1)
        check_quotient_cover([branch], {0: 9, 5: 2, 15: 1}, 5)

    def test_missing_eligible_quotient(self):
        branch = self.branch()
        branch["quotient_dimension_interval"] = [6, 14]
        with self.assertRaises(ValueError):
            check_quotient_cover([branch], {5: 2, 15: 1}, 5)

    def test_overlapping_intervals(self):
        branch = self.branch()
        with self.assertRaises(ValueError):
            check_quotient_cover([branch, copy.deepcopy(branch)], {5: 2, 15: 1}, 5)

    def test_missing_subspace_mass(self):
        branch = self.branch()
        branch["statistics"]["isotropic_subspaces"] += 1
        with self.assertRaises(ValueError):
            check_branch(branch, {5: 2, 15: 1}, 1)

    def test_q5_orbit_census_cannot_replace_raw_domain(self):
        branch = self.branch()
        branch["enumeration_mode"] = "hybrid_raw_and_marked_code_orbits"
        branch["statistics"].update(raw_enumerated_supports=0, marked_orbit_enumerated_supports=2,
                                    isotropic_subspaces=None, isotropic_subspace_statistics_exact=False,
                                    radical_statistics_exact=False)
        with self.assertRaises(ValueError):
            check_branch(branch, {5: 2, 15: 1}, 1)

    def test_alternating_tensor_exclusion_is_not_valid_for_q3(self):
        branch = self.branch()
        branch.update(q=3, minimum_distance=4, enumeration_mode="alternating_tensor_exclusion", statistics=None)
        with self.assertRaises(ValueError):
            check_branch(branch, {5: 2, 15: 1}, 1)


if __name__ == "__main__":
    unittest.main()
