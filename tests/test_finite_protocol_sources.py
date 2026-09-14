import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from space_codec import validate_unital_support
from verify_finite_protocol_sources import census_metric, raw_pointing_counts, source_keys
from verify_small_source_closure import all_pointings


class FiniteProtocolSourcesTests(unittest.TestCase):
    def test_translation_period_count_matches_all_pointings(self):
        a = set(range(16))
        b = {i + (j << 4) for i in range(4) for j in range(4)}
        supports = [(4, tuple(range(16))), (6, tuple(sorted(a ^ b)))]
        for m, points in supports:
            validate_unital_support(points, m)
            actual = [0] * 5
            for row in all_pointings(points, m, 54):
                actual[row["parity_case"]] += 1
            self.assertEqual(raw_pointing_counts(points, m), actual)

    def test_source_partition_rejects_overlap_and_count_gaps(self):
        block = dict(c=44, m=9, source_interval=[3, 5], source_spaces=2)
        data = dict(source_spaces=2, blocks=[block])
        self.assertEqual(source_keys(data), {(44, 9, 3), (44, 9, 4)})
        with self.assertRaises(ValueError):
            source_keys(dict(source_spaces=4, blocks=[block, block]))
        with self.assertRaises(ValueError):
            source_keys(dict(source_spaces=3, blocks=[block]))

    def test_metric_formats_agree(self):
        compact = dict(output_key_words=["0x01"], d_Z=3, n=15, S=5)
        native = dict(output=dict(canonical_key_words=["0x01"]), d_Z=3,
                      protocol_length_n=15, space_footprint_S=5)
        self.assertEqual(census_metric(compact), census_metric(native))


if __name__ == "__main__":
    unittest.main()
