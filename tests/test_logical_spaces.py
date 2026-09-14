from itertools import combinations
from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from logical_spaces import logical_label_space, zero_syndrome_constraints


class LogicalSpaceTests(unittest.TestCase):
    def test_small_zero_sum_subsets_against_exhaustion(self):
        rng = random.Random(743)
        for trial in range(40):
            points = rng.sample(range(32), rng.randrange(1, 11))
            for maximum in range(6):
                expected = set()
                for weight in range(1, maximum+1):
                    for subset in combinations(range(len(points)), weight):
                        value, mask = 0, 0
                        for i in subset:
                            value ^= points[i]
                            mask |= 1 << i
                        if not value:
                            expected.add(mask)
                actual = list(zero_syndrome_constraints(points, maximum))
                self.assertEqual(set(actual), expected)
                self.assertEqual(len(actual), len(expected))

    def test_one_t_label_space(self):
        points = tuple(range(1, 16))
        d3 = logical_label_space(points, 4, 3)
        self.assertEqual(d3["quotient_dimension"], 1)
        self.assertEqual(d3["common_isotropy_forms"], ((0,),)*4)
        self.assertEqual(logical_label_space(points, 4, 4)["quotient_dimension"], 0)
        with self.assertRaises(ValueError):
            logical_label_space(points[:-1], 4, 3)


if __name__ == "__main__":
    unittest.main()
