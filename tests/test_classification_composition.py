import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_classification import check_frontier_scope, check_source_partition


class ClassificationCompositionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frontier = json.loads((ROOT / "data/protocols/pareto_frontier.json").read_bytes())
        cls.cover = json.loads((ROOT / "certificates/through54_protocol_cover.json").read_bytes())

    def test_matching_frontier(self):
        check_frontier_scope(self.frontier, self.cover)

    def test_exact_distance_not_replaced_by_a_floor(self):
        data = copy.deepcopy(self.cover)
        data["exact_distance_is_part_of_key"] = False
        with self.assertRaises(ValueError):
            check_frontier_scope(self.frontier, data)

    def test_missing_metric(self):
        data = copy.deepcopy(self.cover)
        data["frontier_metrics"].pop()
        with self.assertRaises(ValueError):
            check_frontier_scope(self.frontier, data)

    def test_full_clifford_is_not_claimed(self):
        data = copy.deepcopy(self.cover)
        data["output_equivalence"] = "Clifford"
        with self.assertRaises(ValueError):
            check_frontier_scope(self.frontier, data)

    def test_disjoint_complete_parent_partition(self):
        index = dict(lengths=[dict(c=16, sectors=[dict(m=4, included_classes=3)])])
        self.assertEqual(check_source_partition(index, [(16, 4, 0, 1), (16, 4, 1, 3)])["16"], 3)
        for intervals in ([(16, 4, 0, 2)], [(16, 4, 0, 2), (16, 4, 1, 3)], [(16, 5, 0, 3)]):
            with self.assertRaises(ValueError):
                check_source_partition(index, intervals)


if __name__ == "__main__":
    unittest.main()
