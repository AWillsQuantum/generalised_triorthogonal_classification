from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from prepare_zero_fibre_input import source_record, SOURCE_SUFFIX
from verify_zero_fibre_spaces import graph_labels, verify


class ZeroFibreTests(unittest.TestCase):
    def test_affine_cube_graph_quotient(self):
        self.assertEqual(graph_labels(tuple(range(16)), 4), (11, ()))
        rank, basis = graph_labels(tuple(range(32)), 5)
        self.assertEqual((rank, len(basis)), (16, 10))

    def test_reject_invalid_core(self):
        with self.assertRaises(ValueError):
            graph_labels(tuple(range(15)), 4)

    def test_requires_a_proved_zero_fibre_cover(self):
        with self.assertRaisesRegex(ValueError, "not a zero-fibre"):
            verify(ROOT, 54, 10)

    def test_dimension_generic_native_layout(self):
        for m in range(10, 15):
            row = source_record(tuple(range(54)), m, 19, 53-m, 0)
            self.assertEqual(len(row), (1 << m)//8+32)
            self.assertEqual(int.from_bytes(row[:-32], "little"), (1 << 54)-1)
            self.assertEqual(SOURCE_SUFFIX.unpack(row[-32:]), (0, 0, 53-m, 0, 19, 1, 1))

    def test_reject_wrong_rank_nullity(self):
        with self.assertRaisesRegex(ValueError, "rank-nullity"):
            source_record(tuple(range(54)), 12, 0, 40, 2)


if __name__ == "__main__":
    unittest.main()
