from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from space_codec import binary_rank, validate_unital_support
from space_lifts import evaluation_rows, iter_lifts, lift_family, linear_solutions, reconstruct_lift


class LiftTests(unittest.TestCase):
    def test_binary_systems_against_exhaustive_solutions(self):
        rng = random.Random(5455)
        for width in range(1, 7):
            for _ in range(25):
                rows = [rng.randrange(1 << width) for _ in range(7)]
                rhs = [rng.randrange(2) for _ in rows]
                expected = {x for x in range(1 << width)
                            if all((x & r).bit_count() % 2 == b for r, b in zip(rows, rhs))}
                result = linear_solutions(rows, rhs, width)
                if result is None:
                    self.assertFalse(expected)
                    continue
                particular, kernel = result
                actual = {particular}
                for row in kernel:
                    actual |= {x ^ row for x in actual}
                self.assertEqual(actual, expected)

    def test_graph_lifts(self):
        self.assertEqual(lift_family(tuple(range(16)), (), 4), (0, ()))
        self.assertEqual(list(iter_lifts(tuple(range(16)), (), 4)), [])
        core = tuple(range(32))
        _, basis = lift_family(core, (), 5)
        self.assertEqual(len(basis), 10)
        lifts = list(iter_lifts(core, (), 5))
        self.assertEqual(len(lifts), 1023)
        self.assertEqual(len(set(lifts)), 1023)
        for support in lifts:
            self.assertTrue(validate_unital_support(support, 6))

    def test_fibres_with_rank_deficient_core(self):
        # A rank-four quadratic on six variables has a 24-point support.
        support = tuple(x for x in range(64) if ((x >> 5) & 1) * ((x >> 4) & 1)
                        ^ ((x >> 3) & 1) * ((x >> 2) & 1))
        self.assertEqual(len(support), 24)
        self.assertTrue(validate_unital_support(support, 6))
        # Use a change of basis so the contracted direction has both
        # singleton and complete fibres.
        support = tuple(sorted(x ^ ((x & 1) << 5) for x in support))
        fibre_points = {}
        for x in support:
            fibre_points.setdefault(x >> 1, []).append(x & 1)
        core = tuple(sorted(x for x, values in fibre_points.items() if len(values) == 1))
        fibres = tuple(sorted(x for x, values in fibre_points.items() if len(values) == 2))
        self.assertEqual((len(core), len(fibres)), (16, 4))
        self.assertEqual(binary_rank(x ^ core[0] for x in core), 4)
        labels = sum(fibre_points[x][0] << i for i, x in enumerate(core))
        self.assertEqual(reconstruct_lift(core, fibres, labels, 5), support)
        family = lift_family(core, fibres, 5)
        self.assertIsNotNone(family)
        particular, quotient = family
        self.assertTrue(all((r & particular).bit_count() % 2 == b.bit_count() % 2
                            for r, b in zip(evaluation_rows(core, 5, 2), evaluation_rows(fibres, 5, 2))))
        for lifted in iter_lifts(core, fibres, 5):
            self.assertTrue(validate_unital_support(lifted, 6))


if __name__ == "__main__":
    unittest.main()
