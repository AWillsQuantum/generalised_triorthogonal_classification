import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from selected_spaces import iter_interval, read_selected
from space_codec import BITMAP, ShardHeader, write_shard


class SparseCatalogueTests(unittest.TestCase):
    def test_sparse_indices_and_missing_records(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            directory = root / "certificates/support_encoding"
            directory.mkdir(parents=True)
            path = root / "example.utspace"
            points = [tuple(range(i, i+16)) for i in range(5)]
            with path.open("wb") as stream:
                write_shard(stream, ShardHeader(BITMAP, 5, 16, 23, 5), points)
            certificate = dict(status="pass", c=16, m=5, shards=[dict(path=path.name,
                first_index=23, count=5, sha256=hashlib.sha256(path.read_bytes()).hexdigest())])
            (directory / "example.json").write_text(json.dumps(certificate))
            selected, bindings = read_selected(root, {(16, 5, 24), (16, 5, 27)})
            self.assertEqual(selected, {(16, 5, 24): points[1], (16, 5, 27): points[4]})
            self.assertEqual(len(bindings), 1)
            with self.assertRaises(ValueError):
                read_selected(root, {(16, 5, 28)})
            path.write_bytes(path.read_bytes()[:-1])
            with self.assertRaises(ValueError):
                read_selected(root, {(16, 5, 24)})

    def interval_fixture(self, root, intervals=((0, 3), (3, 2))):
        directory = root / "certificates/support_encoding"
        directory.mkdir(parents=True)
        (root / "data/spaces").mkdir(parents=True)
        shards = []
        for number, (first, count) in enumerate(intervals):
            path = root / f"part{number}.utspace"
            with path.open("wb") as stream:
                write_shard(stream, ShardHeader(BITMAP, 5, 16, first, count),
                            [tuple(range(i, i+16)) for i in range(first, first+count)])
            shards.append(dict(path=path.name, first_index=first, count=count,
                               sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        certificate = directory / "sector.json"
        certificate.write_text(json.dumps(dict(status="pass", c=16, m=5, shards=shards)))
        sector = dict(m=5, data_complete=True, included_classes=5, expected_classes=5,
                      encoding_certificates=[dict(path=certificate.relative_to(root).as_posix(),
                        sha256=hashlib.sha256(certificate.read_bytes()).hexdigest())])
        (root / "data/spaces/index.json").write_text(json.dumps(dict(lengths=[dict(c=16, sectors=[sector])])))
        return shards

    def test_streaming_interval_across_shards(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.interval_fixture(root)
            bindings = {}
            result = list(iter_interval(root, 16, 5, 1, 4, bindings))
            self.assertEqual(result, [(i, tuple(range(i, i+16))) for i in range(1, 5)])
            self.assertEqual(len(bindings), 4)
            self.assertEqual(list(iter_interval(root, 16, 5, 5, 0)), [])

    def test_streaming_gaps_overlaps_and_out_of_bounds(self):
        for intervals in (((0, 2), (3, 2)), ((0, 4), (3, 2))):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                self.interval_fixture(root, intervals)
                with self.assertRaises(ValueError):
                    list(iter_interval(root, 16, 5, 0, 5))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.interval_fixture(root)
            for first, count in ((-1, 1), (0, -1), (4, 2)):
                with self.assertRaises(ValueError):
                    list(iter_interval(root, 16, 5, first, count))

    def test_streaming_changed_payload_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            shards = self.interval_fixture(root)
            path = root / shards[1]["path"]
            path.write_bytes(path.read_bytes()[:-1])
            with self.assertRaises(ValueError):
                list(iter_interval(root, 16, 5, 0, 5))


if __name__ == "__main__":
    unittest.main()
