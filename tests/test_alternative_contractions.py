import copy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_alternative_contractions import check_counts, check_cover, contract, difference_signature, valid_cover


class AlternativeContractionTests(unittest.TestCase):
    def test_odd_fibres(self):
        self.assertEqual(contract((0, 1, 2, 4, 5, 7), 1), (1, 3))
        self.assertEqual(contract((0, 1, 2, 4, 5, 7), 2), (1, 2))
        for direction in (0, 256, -1):
            with self.assertRaises(ValueError):
                contract((0, 1), direction)

    def test_invariant_under_affine_coordinates(self):
        points = (0, 2, 5, 6, 9, 12, 16, 23)
        images = tuple(sorted((x ^ (2 if x & 1 else 0)) ^ 37 for x in points))
        self.assertEqual(difference_signature(points), difference_signature(images))
        self.assertFalse(valid_cover(points, 1, set()))

    def test_no_false_complete_cover(self):
        case = dict(target_length=54, profile=dict(raw_marked_pairs=100))
        summary = dict(status="complete_alternative_contraction_cover", target_length=54,
                       candidate_count=7, covered_count=7, unresolved_count=0,
                       raw_marked_mass=100, covered_raw_marked_mass=100, direction_counts={"3": 7})
        check_cover(summary, case, 7)
        for field, value in (("covered_count", 6), ("unresolved_count", 1),
                             ("covered_raw_marked_mass", 99), ("direction_counts", {"0": 7}),
                             ("target_length", 52)):
            bad = copy.deepcopy(summary)
            bad[field] = value
            with self.assertRaises(ValueError):
                check_cover(bad, case, 7)

    def test_rank_deficient_accounting(self):
        profile = dict(core_length=52, full_fibres=0, compatible_fibre_sets=1,
                       lift_quotient_dimension=15, stabiliser_order=8,
                       raw_marked_pairs=32767, rank_deficient_pairs=1)
        case = dict(id="source", profile=profile)
        summary = dict(all_checks_pass=True, count_only=False, parent_id="source", parent_weight=52,
                       full_fiber_multiplicity=0, compatible_fiber_set_count=1, lift_quotient_dimension=15,
                       core_stabilizer_order=8, raw_marked_pair_count=32767, excluded_rank_deficient_pair_count=1)
        check_counts(summary, case)
        summary["excluded_rank_deficient_pair_count"] = 0
        with self.assertRaises(ValueError):
            check_counts(summary, case)


if __name__ == "__main__":
    unittest.main()
