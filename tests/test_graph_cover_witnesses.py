import copy
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_separated_graph_cover import replay


def image(points, frame):
    values = []
    for point in points:
        value = frame[0]
        for bit, column in enumerate(frame[1:]):
            if point >> bit & 1:
                value ^= column
        values.append(value)
    return tuple(sorted(values))


class GraphCoverWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data = json.loads((ROOT / "certificates/zero_fibre_c54_m14.json").read_bytes())
        orbit, lift = next((orbit, lift) for orbit in data["orbits"] for lift in orbit["lifts"]
                           if data["source_profile"][lift["source_index"]]["lift_dimension"] == 1)
        child = image(orbit["representative_points"], lift["representative_to_lift"])
        cls.parent = tuple(x >> 1 for x in child)
        cls.target = image(orbit["representative_points"], orbit["representative_to_target"])
        profile = dict(data["source_profile"][lift["source_index"]], source_index=0)
        single_orbit = dict(orbit, target_index=0, lifts=[dict(lift, source_index=0)])
        cls.fixture = dict(schema=data["schema"], status="pass", c=54, m=14, source_classes=1, target_classes=1,
                           nonaffine_graph_lifts=1, source_profile=[profile], orbits=[single_orbit], dependencies={})

    def run_proof(self, data):
        with patch("verify_separated_graph_cover.sector", side_effect=lambda e, idx, c, m:
                   [self.parent] if m == 13 else [self.target]):
            return replay(ROOT, data)

    def test_complete_single_parent_proof(self):
        self.assertTrue(self.run_proof(self.fixture)["complete_given_predecessor_catalogue"])

    def test_reject_repeated_lift_coset(self):
        data = copy.deepcopy(self.fixture)
        data["orbits"][0]["lifts"] *= 2
        with self.assertRaisesRegex(ValueError, "coset"):
            self.run_proof(data)

    def test_reject_missing_orbit(self):
        data = copy.deepcopy(self.fixture)
        data["orbits"] = []
        with self.assertRaisesRegex(ValueError, "omits"):
            self.run_proof(data)

    def test_reject_affine_graph(self):
        data = copy.deepcopy(self.fixture)
        data["orbits"][0]["lifts"][0]["labels"] = 0
        with self.assertRaisesRegex(ValueError, "coset"):
            self.run_proof(data)

    def test_reject_false_affine_map(self):
        data = copy.deepcopy(self.fixture)
        data["orbits"][0]["lifts"][0]["representative_to_lift"][0] ^= 1
        with self.assertRaisesRegex(ValueError, "does not reproduce"):
            self.run_proof(data)

    def test_reject_incorrect_profile(self):
        data = copy.deepcopy(self.fixture)
        data["source_profile"][0]["nonaffine_lifts"] = 2
        with self.assertRaisesRegex(ValueError, "profile"):
            self.run_proof(data)


if __name__ == "__main__":
    unittest.main()
