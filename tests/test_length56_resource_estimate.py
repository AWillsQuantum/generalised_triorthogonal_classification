"""Check the resource forecast arithmetic, not classification completeness."""

import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from estimate_length56_core_hours import action_orbits, binary_basis, estimate, in_span


class ResourceEstimateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = json.loads((ROOT / "resource_estimates/length56_inputs.json").read_text())

    def test_reproduce_recorded_result(self):
        expected = json.loads((ROOT / "resource_estimates/length56_estimate.json").read_text())
        self.assertEqual(estimate(self.inputs), expected)
        self.assertEqual(expected["estimate_core_hours"], 100_000_000)

    def test_pilot_inventory(self):
        self.assertEqual(len(self.inputs["sampled_n4_cores"]), 16)
        self.assertEqual(sum(len(c["samples"]) for c in self.inputs["sampled_n4_cores"]), 256)
        for core in self.inputs["sampled_n4_cores"]:
            for row in core["samples"]:
                self.assertEqual(len(set(row["support"])), 56)
                self.assertLessEqual(max(map(int, row["profile"]["quotient_dimension_histogram"])), 16)
        self.assertEqual(sum(r["marked_orbit_count"] for r in self.inputs["n4_population_strata"]),
                         14_453_126_733)

    def test_population_weights_are_not_equal_core_weights(self):
        data = copy.deepcopy(self.inputs)
        profile = dict(marked_automorphism_order=1, admissible_direction_orbits=1,
                       quotient_dimension_histogram={"15": 1})
        empty = dict(profile, quotient_dimension_histogram={})
        data["sampled_n4_cores"] = [
            dict(core_id="small", quadratic_nullity=2, marked_orbit_count=1,
                 samples=[dict(profile=profile)]),
            dict(core_id="large", quadratic_nullity=2, marked_orbit_count=99,
                 samples=[dict(profile=empty)]),
        ]
        data["n4_population_strata"] = [dict(quadratic_nullity=2, core_count=2,
                                             marked_orbit_count=100)]
        result = estimate(data)["arithmetic"]["projected_hard_supports"]
        self.assertEqual(result["15"], 1)
        self.assertEqual(result["16"], 0)

    def test_marked_orbit_weight_and_direction_multiplicity(self):
        data = copy.deepcopy(self.inputs)
        data["sampled_n4_cores"] = [dict(
            core_id="example", quadratic_nullity=2, marked_orbit_count=100,
            samples=[
                dict(profile=dict(marked_automorphism_order=3, admissible_direction_orbits=2,
                                  quotient_dimension_histogram={"15": 4})),
                dict(profile=dict(marked_automorphism_order=1, admissible_direction_orbits=1,
                                  quotient_dimension_histogram={})),
            ])]
        data["n4_population_strata"] = [dict(quadratic_nullity=2, core_count=1,
                                             marked_orbit_count=100)]
        self.assertEqual(estimate(data)["arithmetic"]["projected_hard_supports"]["15"], 150)

    def test_small_algebra_helpers(self):
        basis = binary_basis([3, 5, 6])
        self.assertEqual(len(basis), 2)
        self.assertTrue(in_span(6, basis))
        self.assertFalse(in_span(1, basis))
        self.assertEqual(action_orbits(range(4), [[1, 0, 3, 2]]), [{0, 1}, {2, 3}])
        with self.assertRaises(ValueError):
            action_orbits([0], [[1, 0]])


if __name__ == "__main__":
    unittest.main()
