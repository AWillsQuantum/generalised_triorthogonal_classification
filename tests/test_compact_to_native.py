import io
from pathlib import Path
import struct
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from compact_to_native import BASE_RECORD, MANIFEST_HEADER, write_base_manifest


class NativeInputTests(unittest.TestCase):
    def test_native_layout(self):
        buffer = io.BytesIO()
        supports = [tuple(range(16)), tuple(range(16, 32))]
        write_base_manifest(buffer, supports, c=16, m=5, first_index=27, count=2)
        buffer.seek(0)
        self.assertEqual(MANIFEST_HEADER.unpack(buffer.read(MANIFEST_HEADER.size)),
                         (b"UTSPTS1\0", 1, 2, 2, 54, 3))
        for name in (b"c16_m5_i27", b"c16_m5", b"c16_m5_i28", b"c16_m5"):
            size, = struct.unpack("<H", buffer.read(2))
            self.assertEqual(buffer.read(size), name)
        for index, support in enumerate(supports):
            self.assertEqual(BASE_RECORD.unpack(buffer.read(BASE_RECORD.size)), (index, 16, 5, 5, 0, 16, 0))
            self.assertEqual(struct.unpack("<16I", buffer.read(64)), support)
        self.assertEqual(buffer.read(), b"")

    def test_incorrect_count_and_support(self):
        for supports, count in (([], 1), ([tuple(range(16))], 0), ([tuple([0] * 16)], 1)):
            with self.assertRaises(ValueError):
                write_base_manifest(io.BytesIO(), supports, c=16, m=4, first_index=0, count=count)


if __name__ == "__main__":
    unittest.main()
