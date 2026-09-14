"""Translate a compact support interval into the native mathematical input format."""

import argparse
import hashlib
from itertools import islice
import json
from pathlib import Path
import struct

from read_spaces import open_shard
from space_codec import checked_support, read_header, read_records


MANIFEST_HEADER = struct.Struct("<8sIIQII")
BASE_RECORD = struct.Struct("<IHBBBBi")


def write_base_manifest(stream, points, *, c, m, first_index, count, length_limit=54):
    if not 0 < c <= 64 or c & 1 or not 0 <= m <= 26:
        raise ValueError("Invalid native support length or affine dimension")
    if not 0 < length_limit <= 64 or not 0 <= count < 1 << 32 or first_index < 0:
        raise ValueError("Invalid native interval or length limit")
    stream.write(MANIFEST_HEADER.pack(b"UTSPTS1\0", 1, count, count, length_limit, 3))
    family = f"c{c}_m{m}".encode("ascii")
    for index in range(first_index, first_index + count):
        name = f"c{c}_m{m}_i{index}".encode("ascii")
        for value in (name, family):
            stream.write(struct.pack("<H", len(value)))
            stream.write(value)
    seen = 0
    for seen, support in enumerate(points, 1):
        if seen > count:
            raise ValueError("More supports than declared")
        checked_support(support, m, c)
        stream.write(BASE_RECORD.pack(seen - 1, c, m, m, 0, c, 0))
        stream.write(struct.pack(f"<{c}I", *support))
    if seen != count:
        raise ValueError("Fewer supports than declared")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start", type=int, default=0, help="Zero-based offset within the compact shard")
    parser.add_argument("--count", type=int)
    parser.add_argument("--length-limit", type=int, default=54)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output already exists")
    with args.input.open("rb") as raw:
        source_hash = hashlib.file_digest(raw, "sha256").hexdigest()
    temporary = args.output.with_name(args.output.name + ".partial")
    if temporary.exists():
        parser.error("Partial output already exists")
    with open_shard(args.input) as decoded:
        header = read_header(decoded)
        count = header.count - args.start if args.count is None else args.count
        if args.start < 0 or count < 0 or args.start + count > header.count:
            parser.error("Requested interval is outside the compact shard")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        try:
            with temporary.open("xb") as target:
                write_base_manifest(target,
                                    islice(read_records(decoded, header), args.start, args.start + count),
                                    c=header.c, m=header.m, first_index=header.first_index + args.start,
                                    count=count, length_limit=args.length_limit)
            temporary.replace(args.output)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    with args.output.open("rb") as raw:
        output_hash = hashlib.file_digest(raw, "sha256").hexdigest()
    print(json.dumps(dict(schema="triorthogonal-native-input-v1", status="pass",
                          c=header.c, m=header.m, first_index=header.first_index + args.start,
                          count=count, source_sha256=source_hash, manifest_sha256=output_hash,
                          maximum_protocol_length=args.length_limit)))


if __name__ == "__main__":
    main()
