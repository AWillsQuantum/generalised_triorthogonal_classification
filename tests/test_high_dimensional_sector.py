from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_high_dimensional_sector import check_intervals


class SectorIntervalTests(unittest.TestCase):
    def test_exact_cover(self):
        check_intervals([(4, 9), (0, 4), (9, 12)], 12)
        check_intervals([], 0)

    def test_missing_or_overlapping(self):
        for intervals, count in (([(0, 3), (4, 9)], 9), ([(0, 4), (3, 9)], 9),
                                 ([(0, 3)], 4), ([(0, 0)], 0), ([(0, 5)], 4)):
            with self.assertRaises(ValueError):
                check_intervals(intervals, count)
