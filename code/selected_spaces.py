"""Read sparse catalogue selections without decoding intervening records."""

import hashlib
import json
from pathlib import Path

from read_spaces import open_shard
from space_codec import decode_support, read_exact, read_header


def iter_interval(root, c, m, first, count, bindings=None):
    """Yield (index, support) with bounded memory and checked interval coverage."""
    root = Path(root).resolve()
    if any(type(x) is not int for x in (c, m, first, count)) or first < 0 or count < 0:
        raise ValueError("Invalid compact catalogue interval")
    if bindings is None:
        bindings = {}

    def bound(name, expected=None):
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError("External catalogue path")
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if expected is not None and digest != expected:
            raise ValueError("Changed catalogue binding")
        if name in bindings and bindings[name] != digest:
            raise ValueError("Inconsistent catalogue binding")
        bindings[name] = digest
        return path

    index = json.loads(bound("data/spaces/index.json").read_bytes())
    sectors = [s for row in index["lengths"] if row["c"] == c for s in row["sectors"] if s["m"] == m]
    if len(sectors) != 1:
        raise ValueError("The requested sector is not uniquely declared")
    sector, = sectors
    if not sector["data_complete"] or sector["included_classes"] != sector["expected_classes"] or first+count > sector["included_classes"]:
        raise ValueError("The interval is outside the complete compact sector")
    if count == 0:
        return
    shards = []
    for item in sector["encoding_certificates"]:
        certificate = json.loads(bound(item["path"], item["sha256"]).read_bytes())
        if certificate["status"] != "pass" or (certificate["c"], certificate["m"]) != (c, m):
            raise ValueError("Wrong encoding certificate scope")
        for shard in certificate["shards"]:
            start, size = shard["first_index"], shard["count"]
            if max(first, start) < min(first+count, start+size):
                shards.append(shard)
    shards.sort(key=lambda row: row["first_index"])
    position = first
    for shard in shards:
        start, size = shard["first_index"], shard["count"]
        left, right = max(first, start), min(first+count, start+size)
        if left != position:
            raise ValueError("The selected shard cover has a gap or overlap")
        position = right
    if position != first+count:
        raise ValueError("The selected shard cover is incomplete")
    for shard in shards:
        start, size = shard["first_index"], shard["count"]
        left, right = max(first, start), min(first+count, start+size)
        path = bound(shard["path"], shard["sha256"])
        with open_shard(path) as stream:
            header = read_header(stream)
            if (header.c, header.m, header.first_index, header.count) != (c, m, start, size):
                raise ValueError("Wrong compact source interval")
            remaining = (left-start)*header.width
            while remaining:
                step = min(remaining, 1 << 20)
                read_exact(stream, step)
                remaining -= step
            for i in range(left, right):
                yield i, decode_support(read_exact(stream, header.width), m, c, header.codec)
            if right == start+size and stream.read(1):
                raise ValueError("Unexpected trailing compact catalogue data")


def read_selected(root, keys):
    root = Path(root).resolve()
    keys = set(keys)
    values, bindings = {}, {}
    for path in sorted((root / "certificates/support_encoding").glob("*.json")):
        certificate = json.loads(path.read_bytes())
        c, m = certificate["c"], certificate["m"]
        indices = sorted(i for cc, mm, i in keys if (cc, mm) == (c, m))
        if not indices:
            continue
        if certificate["status"] != "pass":
            raise ValueError("Unverified compact encoding")
        for shard in certificate["shards"]:
            start, count = shard["first_index"], shard["count"]
            selection = [i for i in indices if start <= i < start+count]
            if not selection:
                continue
            source = (root / shard["path"]).resolve()
            if not source.is_relative_to(root):
                raise ValueError("External catalogue path")
            with source.open("rb") as raw:
                digest = hashlib.file_digest(raw, "sha256").hexdigest()
            if digest != shard["sha256"]:
                raise ValueError("Changed source catalogue")
            bindings[shard["path"]] = digest
            with open_shard(source) as stream:
                header = read_header(stream)
                if (header.c, header.m, header.first_index, header.count) != (c, m, start, count):
                    raise ValueError("Wrong compact source interval")
                position = start
                for index in selection:
                    remaining = (index-position)*header.width
                    while remaining:
                        size = min(remaining, 1 << 20)
                        read_exact(stream, size)
                        remaining -= size
                    key = c, m, index
                    if key in values:
                        raise ValueError("Overlapping source catalogue intervals")
                    values[key] = decode_support(read_exact(stream, header.width), m, c, header.codec)
                    position = index+1
    if set(values) != keys:
        raise ValueError("A requested catalogue record is absent")
    return values, bindings
