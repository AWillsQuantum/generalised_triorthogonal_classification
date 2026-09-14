from collections import Counter
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_space_classification_cover import compose_sectors, dimension_bound


class SpaceInductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT / "certificates/space_classification_cover.json").read_bytes())
        cls.counts = {(r["c"], r["m"]): r["classes"] for r in cls.data["sectors"]}
        cls.base = Counter({key: count for key, count in cls.counts.items() if key[1] <= 7})
        cls.routes = {(r["c"], r["m"]): dict(classes=r["classes"], proof=r["proof"])
                      for r in cls.data["sectors"] if r["m"] > 7 and r["proof"] != "empty_zero_fibre_successor"}

    def test_complete_induction_and_conservative_bounds(self):
        self.assertEqual([dimension_bound(c) for c in (48, 50, 52, 54)], [16, 16, 17, 18])
        rows = compose_sectors(self.counts, self.base, self.routes)
        self.assertEqual(sum(r["classes"] for r in rows), 301029259)

    def test_nonempty_sector_cannot_be_omitted(self):
        routes = copy.deepcopy(self.routes)
        del routes[44, 9]
        with self.assertRaises(ValueError):
            compose_sectors(self.counts, self.base, routes)

    def test_empty_successor_requires_empty_predecessor(self):
        routes = copy.deepcopy(self.routes)
        del routes[54, 16]
        with self.assertRaises(ValueError):
            compose_sectors(self.counts, self.base, routes)

    def test_no_nonempty_class_outside_dimension_bound(self):
        counts = dict(self.counts)
        counts[54, 19] = 1
        with self.assertRaises(ValueError):
            compose_sectors(counts, self.base, self.routes)


if __name__ == "__main__":
    unittest.main()
