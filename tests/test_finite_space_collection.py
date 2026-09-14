import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_finite_space_censuses import check_frame


class FiniteSpaceCollectionTests(unittest.TestCase):
    def test_valid_affine_frame(self):
        check_frame((0, 1, 2), (0, 1, 2), (0, 1, 2), 2)
        check_frame((1, 2, 3), (0, 1, 2), (3, 1, 2), 2)

    def test_invalid_frames_and_point_domains(self):
        for points, canonical, frame in (
                ((0, 1, 2), (0, 1, 2), (0, 1, 1)),
                ((0, 1, 2), (0, 1, 6), (0, 1, 2)),
                ((0, 1, 2), (0, 1, 2), (3, 1, 2)),
                ((0, 1, 2), (0, 1, 1), (0, 1, 2))):
            with self.assertRaises(ValueError):
                check_frame(points, canonical, frame, 2)


if __name__ == "__main__":
    unittest.main()
