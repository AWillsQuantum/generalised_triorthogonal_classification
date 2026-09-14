from pathlib import Path
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from finite_space_census import replay_quotient
from space_codec import validate_unital_support


class FiniteQuotientReplayTests(unittest.TestCase):
    def fixture(self, directory, mutation=None):
        points = (*range(16), *(128+(a << 3) for a in range(16)))
        validate_unital_support(points, 8)
        first = sum(1 << x for x in points).to_bytes(32, "little")
        last = sum(1 << (x ^ 1) for x in points).to_bytes(32, "little")
        ledger, buckets, classes, assignments = [Path(directory)/name for name in ("ledger", "buckets", "classes", "assignments")]
        ledger.write_bytes(bytes(32)+first+bytes(32)+last)
        buckets.write_bytes(bytes(32)+struct.pack("<QQ", 0, 2)+first+last)
        classes.write_bytes(struct.pack("<QII", 0, 0, 2 if mutation == "position" else 0)+first+struct.pack("<Q", 1 if mutation == "members" else 2))
        rows = [struct.pack("<I9H2x", 0, origin, *(1 << bit for bit in range(8)))
                for origin in (0, 3 if mutation == "map" else 1)]
        assignments.write_bytes(b"".join(rows)[:-1] if mutation == "truncation" else b"".join(rows))
        return ledger, buckets, classes, assignments

    def test_every_affine_map_and_member_count(self):
        with tempfile.TemporaryDirectory() as directory:
            result = replay_quotient(*self.fixture(directory), 32, 8)
            self.assertEqual(result["candidate_count"], 2)
            self.assertEqual(result["affine_classes"], 1)
            self.assertEqual(result["member_counts"], [2])

    def test_false_maps_members_positions_and_truncation_fail(self):
        for mutation in ("map", "members", "position", "truncation"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError):
                    replay_quotient(*self.fixture(directory, mutation), 32, 8)


if __name__ == "__main__":
    unittest.main()
