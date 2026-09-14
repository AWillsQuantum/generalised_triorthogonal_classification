from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from enumerate_source_block import check_raw_census


class SourceBlockCensusTests(unittest.TestCase):
    def empty(self):
        return dict(status="complete", logical_qubits=8, minimum_distance=3,
            matrix_scope="full_projective", enumeration_mode="raw_isotropic_subspaces",
            support_range=dict(start=0, count=2, end_exclusive=2), pareto_protocols=[],
            statistics=dict(source_spaces=1, supports_processed=2,
                quotient_dimension_filtered_supports=0, eligible_supports=1,
                raw_enumerated_supports=1, marked_orbit_enumerated_supports=0,
                isotropic_subspace_statistics_exact=True, isotropic_subspaces=0,
                radical_statistics_exact=True, nondegenerate_subspaces=0,
                marked_code_orbits=0, canonical_output_orbits=0, radical_dimension_counts={}))

    def test_valid_zero_census(self):
        stats, witnesses = check_raw_census(self.empty(), 8, 3, 2)
        self.assertEqual(stats["isotropic_subspaces"], 0)
        self.assertEqual(witnesses, [])

    def test_no_unaccounted_subspaces(self):
        value = self.empty()
        value["statistics"]["isotropic_subspaces"] = 1
        with self.assertRaises(ValueError):
            check_raw_census(value, 8, 3, 2)

    def test_no_missing_output_witness(self):
        value = self.empty()
        value["statistics"]["canonical_output_orbits"] = 1
        with self.assertRaises(ValueError):
            check_raw_census(value, 8, 3, 2)

    def test_orbit_count_cannot_replace_raw_count(self):
        value = self.empty()
        value["statistics"]["isotropic_subspace_statistics_exact"] = False
        with self.assertRaises(ValueError):
            check_raw_census(value, 8, 3, 2)


if __name__ == "__main__":
    unittest.main()
