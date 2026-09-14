from functools import reduce
from itertools import combinations
from operator import xor
from pathlib import Path
import random
import sys
import unittest

CODE = Path(__file__).resolve().parents[1] / "code"
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(CODE / "space_algorithms"))
from decomposition_audit import decomposition_profile
from inverse_contraction_extensions import (CoreExtensionContext, count_zero_xor_index_subsets,
    count_zero_xor_index_subsets_walsh, iter_zero_xor_index_subsets, xor_quadratic_evaluations)
from space_codec import generator_rows
from space_lifts import lift_family


class SpaceAlgorithmTests(unittest.TestCase):
    def test_xor_iterators_and_independent_counts(self):
        rng = random.Random(174)
        for length in range(11):
            for trial in range(8):
                labels = [rng.randrange(16) for _ in range(length)]
                for size in range(6):
                    expected = {indices for indices in combinations(range(length), size)
                                if reduce(xor, (labels[i] for i in indices), 0) == 0}
                    actual = list(iter_zero_xor_index_subsets(labels, size))
                    self.assertEqual(set(actual), expected)
                    self.assertEqual(len(actual), len(expected))
                    self.assertEqual(count_zero_xor_index_subsets(labels, size), len(expected))
                    self.assertEqual(count_zero_xor_index_subsets_walsh(labels, size), len(expected))

    def test_lift_quotient_against_reference(self):
        for dimension, core, fibres in ((4, tuple(range(16)), ()),
                                        (5, tuple(range(32)), ()),
                                        (5, tuple(range(16)), (16, 17, 18, 19))):
            context = CoreExtensionContext.build(dimension, core)
            rhs = xor_quadratic_evaluations(fibres, dimension)
            reference = lift_family(core, fibres, dimension)
            if reference is None:
                self.assertIsNone(context.solver.solve(rhs))
                continue
            particular, kernel = reference
            affine = context.affine_rows
            def canonical(value):
                return min(value ^ reduce(xor, (row for i, row in enumerate(affine)
                                                if mask >> i & 1), 0)
                           for mask in range(1 << len(affine)))
            expected = {canonical(particular ^ reduce(xor, (row for i, row in enumerate(kernel)
                                                            if mask >> i & 1), 0))
                        for mask in range(1 << len(kernel))}
            self.assertEqual({canonical(x) for x in context.lift_representatives(rhs)}, expected)

    def test_two_decomposition_methods(self):
        rows = generator_rows(tuple(range(16)), 4)
        self.assertTrue(decomposition_profile(rows)["indecomposable"])
        block = tuple(row + "0" * 16 for row in rows) + tuple("0" * 16 + row for row in rows)
        result = decomposition_profile(block)
        self.assertTrue(result["methods_agree"])
        self.assertEqual(result["matroid_component_lengths"], [16, 16])


if __name__ == "__main__":
    unittest.main()
