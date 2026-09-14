from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_selected_cubic_sector import check_profile


class SelectedCubicProfileTests(unittest.TestCase):
    def profile(self):
        return dict(pointing_index=0, enumeration_mode="hybrid_raw_and_marked_code_orbits",
            output_keys=[[0x604120]], statistics=dict(source_spaces=1, supports_processed=1,
                quotient_dimension_filtered_supports=0, eligible_supports=1,
                raw_enumerated_supports=0, marked_orbit_enumerated_supports=1,
                isotropic_subspace_statistics_exact=False, radical_statistics_exact=False,
                isotropic_subspaces=None, radical_dimension_counts={},
                canonical_output_orbits=1, marked_code_orbits=2, nondegenerate_subspaces=72))

    def test_complete_weighted_profile(self):
        self.assertEqual(check_profile(self.profile(), 5, 15), {(0x604120,)})

    def test_duplicate_output(self):
        value = self.profile()
        value["output_keys"] *= 2
        with self.assertRaises(ValueError):
            check_profile(value, 5, 15)

    def test_marked_mass_not_raw_mass(self):
        value = self.profile()
        value["statistics"]["isotropic_subspaces"] = 72
        with self.assertRaises(ValueError):
            check_profile(value, 5, 15)

    def test_no_missing_output(self):
        value = self.profile()
        value["output_keys"] = []
        with self.assertRaises(ValueError):
            check_profile(value, 5, 15)

    def test_ineligible_quotient(self):
        with self.assertRaises(ValueError):
            check_profile(self.profile(), 5, 4)


if __name__ == "__main__":
    unittest.main()
