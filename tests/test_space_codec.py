import io
from itertools import combinations
from math import comb
from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from space_codec import (BITMAP, COLEX, HEADER, ShardHeader, decode_support,
                         encode_support, generator_rows, indicator_anf,
                         polynomial_text, rank_support, read_header,
                         read_records, record_width, unrank_support,
                         validate_unital_support, write_shard)


class SupportCodecTests(unittest.TestCase):
    def test_exhaustive_fixed_weight_bijection(self):
        for m in range(4):
            for c in range((1 << m) + 1):
                seen = set()
                for points in combinations(range(1 << m), c):
                    rank = rank_support(points, m)
                    seen.add(rank)
                    self.assertEqual(points, unrank_support(rank, m, c))
                    for codec in (COLEX, BITMAP):
                        self.assertEqual(points, decode_support(encode_support(points, m, codec), m, c, codec))
                self.assertEqual(seen, set(range(comb(1 << m, c))))

    def test_random_large_domains(self):
        rng = random.Random(5481)
        for m in range(6, 18):
            for _ in range(4):
                points = tuple(sorted(rng.sample(range(1 << m), 54)))
                for codec in (COLEX, BITMAP):
                    self.assertEqual(points, decode_support(encode_support(points, m, codec), m, 54, codec))

    def test_invalid_records(self):
        for points in ((1, 1), (2, 1), (-1,), (8,)):
            with self.assertRaises(ValueError):
                encode_support(points, 3)
        for rank in (-1, comb(16, 4)):
            with self.assertRaises(ValueError):
                unrank_support(rank, 4, 4)
        with self.assertRaises(ValueError):
            decode_support(b"\xff", 2, 2, BITMAP)
        with self.assertRaises(ValueError):
            record_width(3, 4, 3)
        with self.assertRaises(ValueError):
            decode_support(b"", 3, 2)

    def test_matrix_variable_order(self):
        self.assertEqual(generator_rows((1, 2, 4), 3), ["111", "001", "010", "100"])

    def test_indicator_polynomial_truth_table(self):
        rng = random.Random(361)
        for m in range(1, 8):
            points = tuple(sorted(rng.sample(range(1 << m), (1 << m) // 2)))
            monomials = indicator_anf(points, m)
            recovered = tuple(x for x in range(1 << m)
                              if sum((x & mask) == mask for mask in monomials) % 2)
            self.assertEqual(points, recovered)

    def test_unused_variable_keeps_ambient_dimension(self):
        points = tuple(range(8, 16))
        self.assertEqual(indicator_anf(points, 4), (8,))
        self.assertEqual(polynomial_text((8,), 4), "x_1")
        self.assertEqual(len(generator_rows(points, 4)), 5)

    def test_unital_validity(self):
        self.assertTrue(validate_unital_support(tuple(range(16)), 4))
        with self.assertRaises(ValueError):
            validate_unital_support(tuple(range(15)), 4)
        with self.assertRaises(ValueError):
            validate_unital_support(tuple(range(16)), 5)

    def test_shard_round_trip(self):
        points = [(0, 1, 3), (2, 4, 5)]
        for codec in (COLEX, BITMAP):
            header = ShardHeader(codec, 3, 3, 200, 2)
            stream = io.BytesIO()
            write_shard(stream, header, points)
            self.assertEqual(len(stream.getvalue()), HEADER.size + 2 * header.width)
            stream.seek(0)
            self.assertEqual(read_header(stream), header)
            self.assertEqual(list(read_records(stream, header)), points)

    def test_shard_rejects_truncation_and_extra_data(self):
        header = ShardHeader(COLEX, 3, 2, 0, 1)
        good = header.pack() + encode_support((1, 3), 3)
        for data in (good[:-1], good + b"x", b"bad"):
            with self.assertRaises(ValueError):
                stream = io.BytesIO(data)
                list(read_records(stream, read_header(stream)))
        with self.assertRaises(ValueError):
            write_shard(io.BytesIO(), header, [])
        with self.assertRaises(ValueError):
            write_shard(io.BytesIO(), header, [(1, 3), (1, 4)])


if __name__ == "__main__":
    unittest.main()
