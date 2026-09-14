from itertools import combinations
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_tensor_authorities import standard_tensor, subspace_bases
from space_codec import binary_rank


class TensorAuthorityTests(unittest.TestCase):
    def test_standard_signature_coordinates(self):
        for q in range(1, 6):
            coordinates = [indices for d in (1, 2, 3) for indices in combinations(range(q), d)]
            for position, selected in enumerate(coordinates):
                tensor = standard_tensor(q, 1 << position)
                for indices in coordinates:
                    args = tuple(1 << i for i in indices)
                    args += (args[-1],) * (3 - len(args))
                    self.assertEqual(tensor.value(*args), indices == selected)

    def test_echelon_spaces(self):
        for n, k, count in ((3, 1, 7), (3, 2, 7), (4, 3, 15), (5, 3, 155), (5, 4, 31)):
            bases = list(subspace_bases(n, k))
            spans = set()
            for basis in bases:
                self.assertEqual(binary_rank(basis), k)
                span = {0}
                for row in basis:
                    span |= {v ^ row for v in span}
                spans.add(tuple(sorted(span)))
            self.assertEqual(len(bases), count)
            self.assertEqual(len(spans), count)


if __name__ == "__main__":
    unittest.main()
