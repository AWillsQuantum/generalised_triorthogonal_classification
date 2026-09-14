import itertools
from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from protocol_checks import (LogicalTensor, check_matrix, exact_distance_and_coefficient,
                             find_equivalence_basis, gate_tensor, verify_basis)


class ProtocolChecksTests(unittest.TestCase):
    def test_distance_against_all_subsets(self):
        rng = random.Random(541)
        for _ in range(30):
            columns = [rng.randrange(16) for _ in range(9)]
            best, count = 10, 0
            for subset in range(1 << len(columns)):
                value = 0
                for i, column in enumerate(columns):
                    if subset & (1 << i):
                        value ^= column
                if 0 < value < 4:
                    weight = subset.bit_count()
                    if weight < best:
                        best, count = weight, 1
                    elif weight == best:
                        count += 1
            self.assertEqual(exact_distance_and_coefficient(columns, 2),
                             (None, 0) if not count else (best, count))

    def test_gate_tensor_matches_parity_identities(self):
        for q, gate, term in ((1, "T1", (0,)), (2, "CS12", (0, 1)),
                               (3, "CCZ123", (0, 1, 2))):
            tensor = gate_tensor(q, gate)
            self.assertTrue(tensor.intrinsic())
            for degree in (1, 2, 3):
                for indices in itertools.combinations(range(q), degree):
                    args = list(indices) + [indices[-1]] * (3 - degree)
                    self.assertEqual(tensor.value(*(1 << i for i in args)), int(indices == term))

    def test_basis_recovery(self):
        for q, gate in ((3, "T1CS12CCZ123"), (4, "CS12CCZ134"),
                        (5, "CCZ145CCZ235")):
            source = gate_tensor(q, gate)
            basis = [1 << i for i in range(q)]
            for i in range(q - 1):
                basis[i] ^= basis[i+1]
            target = LogicalTensor([source.words[v] for v in basis])
            recovered = find_equivalence_basis(source, target)
            self.assertIsNotNone(recovered)
            self.assertTrue(verify_basis(source, target, recovered))

    def test_inequivalent_tensors(self):
        self.assertIsNone(find_equivalence_basis(gate_tensor(2, "T1T2"), gate_tensor(2, "CS12")))

    def test_two_qubit_isomorphism_against_all_bases(self):
        tensors = [LogicalTensor((0, 0))]
        factors = ("T1", "T2", "CS12")
        tensors.extend(gate_tensor(2, "".join(gate for i, gate in enumerate(factors) if mask & (1 << i)))
                       for mask in range(1, 8))
        bases = list(itertools.permutations((1, 2, 3), 2))
        for source in tensors:
            for target in tensors:
                expected = False
                for basis in bases:
                    if all(source.value(basis[i], basis[j], basis[k]) == target.value(1 << i, 1 << j, 1 << k)
                           for i, j, k in itertools.product(range(2), repeat=3)):
                        expected = True
                self.assertEqual(find_equivalence_basis(source, target) is not None, expected)

    def test_matrix_rejections(self):
        for rows in (("10", "10"), ("111", "110"), ("101", "00"), ("12", "10")):
            with self.assertRaises(ValueError):
                check_matrix(rows, 1)


if __name__ == "__main__":
    unittest.main()
