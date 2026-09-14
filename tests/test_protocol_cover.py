import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from verify_protocol_cover import Evidence, check52_partition
from verify_pointed_sources import check_raw_block_position


class EvidenceTests(unittest.TestCase):
    def test_certificate_requires_exact_scope_and_dependencies(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = b'{"count": 3}'
            (root / "input.json").write_bytes(raw)
            certificate = dict(schema="example-v1", status="pass", complete=True, count=3,
                               dependencies={"input.json": hashlib.sha256(raw).hexdigest()})
            (root / "certificate.json").write_text(json.dumps(certificate))
            evidence = Evidence(root)
            self.assertEqual(evidence.certificate("certificate.json", "example-v1", ("complete",),
                                                  dict(count=3))["count"], 3)
            with self.assertRaisesRegex(ValueError, "scope"):
                evidence.certificate("certificate.json", "example-v1", expected=dict(count=4))
            with self.assertRaisesRegex(ValueError, "implication"):
                evidence.certificate("certificate.json", "example-v1", ("inequivalent",))
            (root / "input.json").write_bytes(b'{}')
            with self.assertRaisesRegex(ValueError, "binding"):
                Evidence(root).certificate("certificate.json", "example-v1")

    def test_path_cannot_leave_release(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "escapes"):
                Evidence(directory).path("../external.json")


class SourcePartitionTests(unittest.TestCase):
    def setUp(self):
        self.domain = dict(source_spaces=13, main_blocks=list(range(5)), separate_sources=["a", "b"])
        self.small = dict(totals=dict(source_spaces=8, larger_quotient_supports=7, pointed_supports=80),
                          source_blocks=2, selected_cubic_source_blocks=3, separately_certified_sources=["a", "b"])
        self.selected = dict(source_spaces=3, pointings=30, all_q_at_least7_excluded=True)
        self.large = dict(pointings=7, all_q_at_least7_excluded=True)
        self.separate = dict(source_spaces=2, raw_pointings=20, all_q_at_least7_excluded=True)

    def check(self):
        return check52_partition(self.domain, self.small, self.selected, self.large, self.separate)

    def test_deferred_inputs_are_not_counted_twice(self):
        self.assertEqual(self.check()["source_pointings"], 130)

    def test_reject_missing_parent(self):
        self.domain["source_spaces"] += 1
        with self.assertRaisesRegex(ValueError, "exhaust"):
            self.check()

    def test_reject_unclosed_logical_dimension(self):
        self.large["all_q_at_least7_excluded"] = False
        with self.assertRaisesRegex(ValueError, "remains open"):
            self.check()


class PointedBlockPositionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.row = json.loads((ROOT / "data/protocol_sectors/length54_q6_candidates.json").read_bytes())["cases"][0]
        cls.parent = sorted(x ^ cls.row["space"]["origin"] for x in cls.row["points"])
        cls.block = dict(index=2290, c=54, m=8, source_interval=[14072000, 14080000],
                         pointing_mode="all_origins", protocol_length_interval=[53, 54],
                         distance_three_quotient_profile=[[5, 2488000]], pointing_supports=2488000)

    def test_exact_native_input_order(self):
        self.assertEqual(check_raw_block_position(self.row, self.block, self.parent), 87)

    def test_reject_adjacent_pointing(self):
        row = copy.deepcopy(self.row)
        row["block_support_index"] += 1
        with self.assertRaisesRegex(ValueError, "exact source-block position"):
            check_raw_block_position(row, self.block, self.parent)

    def test_reject_different_parent(self):
        row = copy.deepcopy(self.row)
        row["block_support_index"] += 311
        with self.assertRaisesRegex(ValueError, "different parent"):
            check_raw_block_position(row, self.block, self.parent)

    def test_reject_filtered_block_index(self):
        block = copy.deepcopy(self.block)
        block["distance_three_quotient_profile"] = [[4, 1], [5, 2487999]]
        with self.assertRaisesRegex(ValueError, "full raw source-block index"):
            check_raw_block_position(self.row, block, self.parent)


if __name__ == "__main__":
    unittest.main()
