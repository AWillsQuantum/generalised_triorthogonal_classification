import copy
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from affine_graph import transform
from prepare_c48_marked_inputs import stabiliser_generators
from enumerate_c48_marked_spaces import check_contractions


class ProperSpanMarkedInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = json.loads((ROOT / "certificates/low_dimensional_spaces.json").read_bytes())
        cls.form, = [r for r in data["canonical_maps"] if (r["c"], r["m"], r["index"]) == (40, 6, 0)]
        mask = transform(int(cls.form["canonical_mask"], 16), cls.form["frame"], 6)
        cls.support = tuple(x for x in range(64) if mask >> x & 1)

    def test_full_ambient_stabiliser_kernel(self):
        generators, order = stabiliser_generators(self.support, 6, self.form)
        self.assertEqual(order, 64*self.form["stabilizer_order"])
        for images in generators:
            self.assertEqual(sorted(images), list(range(128)))
            self.assertEqual({images[x] for x in self.support}, set(self.support))
        kernel = [g for g in generators if tuple(g[:64]) == tuple(range(64))]
        reached = {64}
        while True:
            larger = reached | {g[x] for g in kernel for x in reached}
            if larger == reached:
                break
            reached = larger
        self.assertEqual(reached, set(range(64, 128)))

    def test_bad_group_and_wrong_core_rejected(self):
        for field, value in (("generators", []), ("stabilizer_order", 11), ("m", 7),
                             ("canonical_mask", "0")):
            form = copy.deepcopy(self.form)
            form[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                stabiliser_generators(self.support, 6, form)

    def test_child_contraction_is_checked(self):
        payload = sum(1 << x for x in (0, 6, 14, 15)).to_bytes(32, "little")
        bitmaps = np.frombuffer(payload, dtype=np.uint8).reshape(1, 32)
        check_contractions(bitmaps, (0, 3), 1)
        with self.assertRaises(ValueError):
            check_contractions(bitmaps, (0, 4), 1)
        with self.assertRaises(ValueError):
            check_contractions(bitmaps, (0, 3), 2)


if __name__ == "__main__":
    unittest.main()
