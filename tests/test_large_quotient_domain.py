from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_large_quotient_domain import pointed_support


class PointingTests(unittest.TestCase):
    def test_all_five_pointing_cases(self):
        support = [0, 1, 2, 3]
        self.assertEqual(pointed_support(support, 3, 0, 4), ([4, 5, 6, 7], 3))
        self.assertEqual(pointed_support(support, 3, 1, 2), ([1, 2, 3], 3))
        self.assertEqual(pointed_support(support, 3, 2, 4), ([0, 4, 5, 6, 7], 3))
        self.assertEqual(pointed_support(support, 3, 3, -1), ([8, 9, 10, 11], 4))
        self.assertEqual(pointed_support(support, 3, 4, -1), ([0, 8, 9, 10, 11], 4))

    def test_invalid_pointing_rejected(self):
        for parity, origin in ((0, -1), (0, 8), (1, 4), (2, 0), (5, 0)):
            with self.subTest(parity=parity, origin=origin), self.assertRaises(ValueError):
                pointed_support([0, 1, 2, 3], 3, parity, origin)


if __name__ == "__main__":
    unittest.main()
