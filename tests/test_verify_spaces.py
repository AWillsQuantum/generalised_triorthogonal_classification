import io
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from space_codec import BITMAP, COLEX, ShardHeader, encode_support
from verify_spaces import checked_content


class ShortReads(io.BytesIO):
    def read(self, count=-1):
        return super().read(min(count, 3) if count >= 0 else 3)


class StreamingAuditTests(unittest.TestCase):
    def collect(self, stream, header, validity=True, reference=False):
        chunks = list(checked_content(stream, header, validity, reference, batch_size=2))
        return sum(n for n, _ in chunks), b"".join(data for _, data in chunks)

    def test_reference_and_batched_digest_bytes(self):
        for m in (4, 5, 6):
            support = tuple(range(1 << m))
            for codec in (BITMAP, COLEX):
                header = ShardHeader(codec, m, len(support), 17, 5)
                payload = encode_support(support, m, codec) * header.count
                expected = self.collect(io.BytesIO(payload), header, reference=True)
                self.assertEqual(expected, self.collect(ShortReads(payload), header))
                self.assertEqual(expected[0], header.count)

    def test_truncated_trailing_and_invalid_data(self):
        header = ShardHeader(BITMAP, 4, 16, 0, 3)
        payload = encode_support(tuple(range(16)), 4, BITMAP) * 3
        for reference in (False, True):
            for bad in (payload[:-1], payload + b"x", b"\x00" + payload[1:]):
                with self.assertRaises(ValueError):
                    self.collect(io.BytesIO(bad), header, reference=reference)
        with self.assertRaises(ValueError):
            list(checked_content(io.BytesIO(payload), header, True, batch_size=0))

    def test_empty_sector_and_rank_deficient_support(self):
        header = ShardHeader(BITMAP, 5, 16, 0, 0)
        self.assertEqual(self.collect(io.BytesIO(), header), (0, b""))
        header = ShardHeader(BITMAP, 5, 16, 0, 1)
        payload = encode_support(tuple(range(16)), 5, BITMAP)
        for reference in (False, True):
            with self.assertRaises(ValueError):
                self.collect(io.BytesIO(payload), header, reference=reference)


if __name__ == "__main__":
    unittest.main()
