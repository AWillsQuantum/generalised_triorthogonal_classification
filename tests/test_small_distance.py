import random
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from protocol_checks import distance_up_to, exact_distance_and_coefficient


class SmallDistanceTests(unittest.TestCase):
    def test_matches_syndrome_dp(self):
        rng = random.Random(1291)
        for q in range(1, 4):
            for _ in range(60):
                columns = [rng.randrange(1 << (q+3)) for _ in range(rng.randrange(1, 13))]
                distance, coefficient = exact_distance_and_coefficient(columns, q)
                expected = (distance, coefficient) if distance is not None and distance <= 5 else (None, 0)
                self.assertEqual(distance_up_to(columns, q), expected)

    def test_weight_five(self):
        # Every proper subset has a nonzero stabiliser syndrome.
        columns = [2, 4, 8, 16, 31]
        self.assertEqual(distance_up_to(columns, 1, 4), (None, 0))
        self.assertEqual(distance_up_to(columns, 1, 5), (5, 1))

    def test_repeated_zero_and_nonzero_columns(self):
        for columns in ([0, 0, 0], [1, 1, 0], [2, 3, 2, 3, 0], []):
            self.assertEqual(distance_up_to(columns, 1), exact_distance_and_coefficient(columns, 1))

    def test_invalid(self):
        for arguments in (([1], 0, 5), ([1], 1, 0), ([-1], 1, 5)):
            with self.assertRaises(ValueError):
                distance_up_to(*arguments)

    def test_zero_column_deletion(self):
        rng = random.Random(53054)
        for q in range(1, 4):
            for _ in range(10):
                columns = rng.sample(range(1, 1 << (q+3)), 12)
                self.assertEqual(exact_distance_and_coefficient(columns, q),
                                 exact_distance_and_coefficient(columns+[0], q))
                self.assertEqual(distance_up_to(columns, q), distance_up_to(columns+[0], q))
