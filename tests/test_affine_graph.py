from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from affine_graph import affine_images, group_order, transform
from verify_low_dimensional_spaces import lower_degree_seeds
from verify_rm37 import polynomial_support


class AffineGraphInputTests(unittest.TestCase):
    def test_frames(self):
        self.assertEqual(affine_images((3, 1, 2, 4), 3), [x ^ 3 for x in range(8)])
        self.assertEqual(transform(0b11, (2, 1, 2, 4), 3), 0b1100)
        for frame in ((0, 1, 1, 4), (0, 1, 2), (8, 1, 2, 4), (0, 1, 2, 8)):
            with self.assertRaises(ValueError):
                affine_images(frame, 3)
        self.assertEqual(group_order(7), 20972799094947840)

    def test_polynomial_domain(self):
        self.assertEqual(polynomial_support([], 3, 2), 0)
        self.assertEqual(polynomial_support([0], 3, 2), 255)
        self.assertEqual(polynomial_support([1, 2], 2, 1), 0b0110)
        for monomials in ([3, 3], [7], [-1], [8]):
            with self.assertRaises(ValueError):
                polynomial_support(monomials, 3, 2)
        self.assertEqual([len(list(lower_degree_seeds(m))) for m in (4, 5, 6)], [2, 3, 11])


if __name__ == "__main__":
    unittest.main()
