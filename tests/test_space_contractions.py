import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from affine_graph import transform
from prepare_space_contractions import task_text


class ContractionInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = json.loads((ROOT / "certificates/low_dimensional_spaces.json").read_bytes())
        cls.form, = [r for r in data["canonical_maps"] if (r["c"], r["m"], r["index"]) == (52, 7, 169)]
        mask = transform(int(cls.form["canonical_mask"], 16), cls.form["frame"], 7)
        cls.support = tuple(x for x in range(128) if mask >> x & 1)

    def test_exact_task_reproduction(self):
        text, profile = task_text("c52_m07_000000169", self.support, self.form)
        path = ROOT / "data/space_contractions/c54_m08/c52_m07_000000169.task"
        self.assertEqual(path.read_text(), text)
        self.assertEqual(profile["raw_marked_pairs"], 1835008)
        self.assertEqual(profile["compatible_fibre_sets"], 28)
        self.assertEqual(profile["lift_quotient_dimension"], 16)

    def test_bad_group_and_wrong_domain_rejected(self):
        for field, value in (("stabilizer_order", 11), ("generators", []),
                              ("canonical_mask", "0"), ("m", 6)):
            form = copy.deepcopy(self.form)
            form[field] = value
            with self.assertRaises(ValueError):
                task_text("case", self.support, form)
        with self.assertRaises(ValueError):
            task_text("case two", self.support, self.form)

    def test_zero_fibre_removes_the_affine_graph(self):
        text, profile = task_text("case", self.support, self.form, 52)
        self.assertTrue(text.startswith("L52_M08_NATIVE_TASK_V1\n"))
        self.assertEqual(profile["full_fibres"], 0)
        self.assertEqual(profile["compatible_fibre_sets"], 1)
        self.assertEqual(profile["rank_deficient_pairs"], 1)
        self.assertEqual(profile["raw_marked_pairs"], (1 << 16)-1)
        for length in (48, 50, 53, 56):
            with self.assertRaises(ValueError):
                task_text("case", self.support, self.form, length)


if __name__ == "__main__":
    unittest.main()
