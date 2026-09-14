from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from space_codec import validate_unital_support
from space_lifts import evaluation_rows, lift_family, reconstruct_lift


class LightSectionMomentTests(unittest.TestCase):
    def test_cubic_section_moments_need_not_vanish(self):
        core = (12, 14, 15, 28, 30, 31, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
                54, 55, 56, 57, 58, 59, 77, 93, 141, 157, 205, 221, 261, 265,
                269, 273, 325, 329, 333, 337, 389, 393, 397, 401, 453, 457, 461, 465)
        section = tuple(range(48, 56))
        validate_unital_support(core, 9)
        self.assertTrue(all(row.bit_count() % 2 == 0 for row in evaluation_rows(section, 9, 2)))
        self.assertTrue(any(row.bit_count() % 2 for row in evaluation_rows(section, 9, 3)))
        labels = sum(1 << i for i, x in enumerate(core) if x in section)
        self.assertEqual(labels, 261120)
        self.assertIsNotNone(lift_family(core, (), 9))
        child = reconstruct_lift(core, (), labels, 9)
        validate_unital_support(child, 10)
        self.assertEqual(len(child), 44)


if __name__ == "__main__":
    unittest.main()
