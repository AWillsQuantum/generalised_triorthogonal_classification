import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_middle_space_censuses import check9, check10, check_upper


class MiddleCensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = {m: json.loads((ROOT / f"data/space_contractions/middle_sectors/c54_m{m:02d}.json").read_bytes()) for m in (9, 10, 12)}
        cls.profile12 = json.loads((ROOT / "certificates/zero_fibre_native_c54_m12_input.json").read_bytes())

    def test_closed_censuses(self):
        self.assertEqual(check9(self.data[9])["affine_classes"], 259202914)
        self.assertEqual(check10(self.data[10])["affine_classes"], 15617101)
        self.assertEqual(check_upper(self.data[12], self.profile12)["affine_classes"], 21585)

    def test_two_fibre_route_cannot_be_silently_dropped(self):
        bad = copy.deepcopy(self.data[9])
        bad["generation_parts"].pop()
        with self.assertRaises(ValueError):
            check9(bad)

    def test_weighted_member_loss_is_rejected(self):
        bad = copy.deepcopy(self.data[9])
        bad["weighted_merge"]["output_member_count"] -= 1
        with self.assertRaises(ValueError):
            check9(bad)
        bad = copy.deepcopy(self.data[9])
        bad["global_exact_quotient"]["independently_replayed_negative_comparison_count"] -= 1
        with self.assertRaises(ValueError):
            check9(bad)

    def test_incomplete_source_profile_or_quotient_rejected(self):
        for mode in ("interval", "rank", "negative"):
            bad = copy.deepcopy(self.data[10])
            if mode == "interval":
                bad["source_profile_intervals"][1]["first_source_index"] += 1
            elif mode == "rank":
                bad["source_profiles"]["0"]["quadratic_profile_distribution"]["rank31_lift14"] = 1
            else:
                bad["exact_quotient"]["independent_negative_replay_count"] = 0
            with self.assertRaises(ValueError):
                check10(bad)

    def test_graph_domain_and_exact_quotient_must_agree(self):
        profile = copy.deepcopy(self.profile12)
        profile["nonaffine_graph_lifts"] -= 1
        with self.assertRaises(ValueError):
            check_upper(self.data[12], profile)
        bad = copy.deepcopy(self.data[12])
        bad["exact_quotient_intervals"][1]["first_candidate"] -= 1
        with self.assertRaises(ValueError):
            check_upper(bad, self.profile12)


if __name__ == "__main__":
    unittest.main()
