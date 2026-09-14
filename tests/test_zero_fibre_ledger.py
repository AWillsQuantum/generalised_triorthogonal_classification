from pathlib import Path
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_zero_fibre_ledger import build_buckets


def record(signature, mask):
    return struct.pack("<4Q", signature, 0, 0, 0)+mask.to_bytes(32, "little")


class SignatureIndexTests(unittest.TestCase):
    def check(self, raw):
        with tempfile.TemporaryDirectory() as directory:
            ledger, index = Path(directory)/"ledger", Path(directory)/"index"
            ledger.write_bytes(raw)
            result = build_buckets(ledger, index, 8)
            return result, index.read_bytes()

    def test_exact_group_boundaries(self):
        first = (1 << 54)-1
        second = first << 1
        third = first << 2
        result, raw = self.check(record(1, first)+record(1, second)+record(2, third))
        self.assertEqual(result, dict(candidate_count=3, signature_bucket_count=2,
                                     bucket_size_distribution={1: 1, 2: 1}))
        self.assertEqual(len(raw), 2*(32+16+64))
        self.assertEqual(struct.unpack_from("<2Q", raw, 32), (0, 2))
        self.assertEqual(struct.unpack_from("<2Q", raw, 112+32), (2, 1))

    def test_repeats_order_truncation_and_bad_weight(self):
        mask = (1 << 54)-1
        bad = (record(1, mask)*2, record(2, mask)+record(1, mask),
               record(1, mask << 1)+record(1, mask), record(1, mask)[:-1], record(1, mask-1))
        for raw in bad:
            with self.subTest(raw_size=len(raw)), self.assertRaises(ValueError):
                self.check(raw)

    def test_empty_ledger(self):
        result, raw = self.check(b"")
        self.assertEqual(result["candidate_count"], 0)
        self.assertEqual(raw, b"")


if __name__ == "__main__":
    unittest.main()
