from pathlib import Path
import random
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from space_batch import (decode_bitmaps, decode_text_supports, encode_bitmaps,
                         ordered_point_bytes, validate_batch)
from space_codec import BITMAP, encode_support, generator_rows, validate_unital_support


class BatchSpaceTests(unittest.TestCase):
    def test_codecs_against_reference(self):
        rng = random.Random(54002)
        for m in range(2, 14):
            c = min(54, 1 << m)
            supports = [tuple(sorted(rng.sample(range(1 << m), c))) for _ in range(7)]
            points = np.array(supports, dtype=np.uint32)
            encoded = encode_bitmaps(points, m, c)
            self.assertEqual(encoded, b"".join(encode_support(s, m, BITMAP) for s in supports))
            self.assertTrue(np.array_equal(decode_bitmaps(encoded, m, c), points))
            expected = b"".join(x.to_bytes((m+7)//8, "little") for s in supports for x in s)
            self.assertEqual(ordered_point_bytes(points, m), expected)
            strings = [[format(x, f"0{m}b") for x in s] for s in supports]
            self.assertTrue(np.array_equal(decode_text_supports(strings, [generator_rows(s, m) for s in supports], m, c), points))

    def test_validity_against_reference(self):
        rng = random.Random(54003)
        for m in range(4, 10):
            for _ in range(12):
                points = sorted(rng.sample(range(1 << m), min(54, 1 << m)))
                try:
                    validate_unital_support(points, m)
                    valid = True
                except ValueError:
                    valid = False
                try:
                    validate_batch(np.array([points], dtype=np.uint32), m, len(points))
                    observed = True
                except ValueError:
                    observed = False
                self.assertEqual(observed, valid)
        for m in (4, 5, 6):
            self.assertTrue(validate_batch(np.array([list(range(1 << m))], dtype=np.uint32), m, 1 << m))

    def test_corruptions(self):
        support = list(range(16))
        strings = [[format(x, "04b") for x in support]]
        matrix = generator_rows(support, 4)
        matrix[2] = "x" + matrix[2][1:]
        with self.assertRaises(ValueError):
            decode_text_supports(strings, [matrix], 4, 16)
        with self.assertRaises(ValueError):
            decode_bitmaps(b"\xff", 4, 16)
        with self.assertRaises(ValueError):
            encode_bitmaps(np.array([[1, 1]], dtype=np.uint32), 4, 2)
        with self.assertRaises(ValueError):
            encode_bitmaps(np.array([[-1, 1]], dtype=np.int32), 4, 2)
        with self.assertRaises(ValueError):
            validate_batch(np.array([list(range(16))], dtype=np.uint32), 5, 16)


if __name__ == "__main__":
    unittest.main()
