from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from space_codec import validate_unital_support
from space_recursion import iter_space_children, multiplicity_bound


class SpaceRecursionTests(unittest.TestCase):
    def test_rank_deficient_core_family(self):
        children = list(iter_space_children(tuple(range(16)), 4, 24, 6))
        self.assertEqual(len(children), 140)
        self.assertEqual(len(set(children)), 140)
        for child in children:
            self.assertTrue(validate_unital_support(child, 6))

    def test_full_rank_graph_family(self):
        children = list(iter_space_children(tuple(range(32)), 5, 32, 6))
        self.assertEqual(len(children), 1023)
        self.assertEqual(len(set(children)), 1023)
        self.assertEqual(list(iter_space_children(tuple(range(16)), 4, 16, 5)), [])
        with self.assertRaises(ValueError):
            list(iter_space_children(tuple(range(16)), 4, 32, 5))

    def test_direction_bounds(self):
        self.assertEqual([multiplicity_bound(54, m) for m in range(8, 12)], [5, 2, 1, 0])


if __name__ == "__main__":
    unittest.main()
