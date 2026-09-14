from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_profile_exclusions import verify


class ProfileExclusionTests(unittest.TestCase):
    def test_positive_subset_with_mixed_word_encoding(self):
        profiles = [{(4,), (9,)}, {(6,)}]
        data = dict(id="example", predecessor_logical_dimension=5, minimum_distance=3,
                    profile_scope="individual_support", profile_count=3, expected_q6_candidates=2,
                    profiles=[dict(index=0, output_keys=[["0x4"], [9], [13]]),
                              dict(index=1, output_keys=[[4]]),
                              dict(index=2, output_keys=[["0x6"]])])
        result = verify(data, profiles)
        self.assertEqual([r["index"] for r in result["candidates"]], [0, 2])
        self.assertEqual(result["excluded"], 1)
        data["expected_q6_candidates"] = 0
        with self.assertRaises(ValueError):
            verify(data, profiles)

    def test_union_is_not_an_individual_positive_claim(self):
        data = dict(id="union", predecessor_logical_dimension=5, minimum_distance=3,
                    profile_scope="block_union", profile_count=1, expected_q6_candidates=1,
                    profiles=[dict(index=0, output_keys=[[4], [9]])])
        result = verify(data, [{(4,), (9,)}])
        self.assertEqual(result["profile_scope"], "block_union")
        self.assertFalse(result["all_q6_absent_on_listed_profiles"])


if __name__ == "__main__":
    unittest.main()
