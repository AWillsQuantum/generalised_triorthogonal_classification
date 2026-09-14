from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_low_dimensional_blocks import check_raw_pointing_count, histogram


class BlockCoverageTests(unittest.TestCase):
    def test_histogram(self):
        self.assertEqual(histogram([[0, 8], [4, 2]], 10), {0: 8, 4: 2})
        for values in ([[0, 8], [0, 2]], [[0, -1], [4, 11]], [[0, 9]]):
            with self.assertRaises(ValueError):
                histogram(values, 10)

    def test_translation_period_counts(self):
        self.assertEqual(check_raw_pointing_count(52, 7, 1, 39), 1)
        self.assertEqual(check_raw_pointing_count(52, 7, 1, 77), 0)
        self.assertEqual(check_raw_pointing_count(52, 7, 2, 116), 1)
        self.assertEqual(check_raw_pointing_count(54, 8, 2, 622), 0)
        for arguments in ((52, 7, 1, 38), (52, 7, 1, 78), (54, 8, 2, 621)):
            with self.assertRaises(ValueError):
                check_raw_pointing_count(*arguments)
