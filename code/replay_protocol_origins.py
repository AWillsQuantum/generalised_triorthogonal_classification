"""Recompute a finite origin quotient and compare its complete pointed domain."""

import argparse
import hashlib
from itertools import groupby
import json
from pathlib import Path
import subprocess

from protocol_domain import (open_domain, read_header, read_native_header,
    read_native_records, read_pointings, reconstruct, source_key, write_base_native)
from selected_spaces import read_selected


def compare_domains(native_path, compact_path, sources, supports, expected_count):
    total = 0
    with native_path.open("rb") as stream, open_domain(compact_path) as original:
        native = read_native_header(stream)
        compact = read_header(original)
        if (native["source_count"] != compact["source_count"] or native["count"] != compact["count"]
                or native["count"] != expected_count or native["maximum_protocol_length"] != compact["maximum_protocol_length"]
                or native["flags"] != 1):
            raise ValueError("Recomputed finite origin domain has different cardinality or scope")
        left = groupby(read_native_records(stream, native), key=lambda row: row["source_index"])
        right = groupby(read_pointings(original, compact), key=lambda row: row[0])
        for index, ((a, actual), (b, expected)) in enumerate(zip(left, right, strict=True)):
            if a != b or a != index:
                raise ValueError("Pointing sources are not in consecutive blocks")
            actual_keys = []
            for row in actual:
                points, h = reconstruct(sources, supports, (a, row["parity_case"], row["origin"]))
                if (tuple(points), h) != (row["points"], row["ambient_dimension"]):
                    raise ValueError("Native pointing is inconsistent with its compact source")
                actual_keys.append((h, tuple(points)))
            expected_keys = [(h, tuple(points)) for points, h in (reconstruct(sources, supports, row) for row in expected)]
            if len(set(actual_keys)) != len(actual_keys) or sorted(actual_keys) != sorted(expected_keys):
                raise ValueError("Recomputed origin quotient does not equal the recorded finite support set")
            total += len(actual_keys)
        if index+1 != len(sources) or total != expected_count:
            raise ValueError("Missing source or pointing in the complete origin quotient")
    return total


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/through48_source_domain.json")
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("Workers must be positive")
    raw = args.input.read_bytes()
    data = json.loads(raw)
    if data["schema"] != "finite-protocol-source-domain-v1":
        raise ValueError("Incorrect source domain")
    compact_path = root / data["pointing_data"]
    with compact_path.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != data["pointing_sha256"]:
            raise ValueError("Changed finite input domain")
    sources = data["sources"]
    supports, bindings = read_selected(root, {source_key(row) for row in sources})
    work = args.work_directory.resolve()
    work.mkdir(parents=True, exist_ok=True)
    base, pointed = work / "base.utsp", work / "pointed.utsp"
    with base.open("wb") as stream:
        write_base_native(stream, sources, supports, data["maximum_protocol_length"])
    process = subprocess.run([str(args.native), "manifest-origin-reduce", "--input", str(base),
        "--output", str(pointed), "--maximum-protocol-length", str(data["maximum_protocol_length"]),
        "--workers", str(args.workers)], capture_output=True, text=True, check=True)
    summary = json.loads(process.stdout)
    total = compare_domains(pointed, compact_path, sources, supports, data["pointing_count"])
    result = dict(schema="finite-origin-quotient-replay-v1", status="pass",
        source_domain_sha256=hashlib.sha256(raw).hexdigest(),
        pointing_sha256=data["pointing_sha256"], sources=len(sources), pointed_supports=total,
        all_origin_orbits_freshly_recomputed=True, every_finite_point_set_agrees=True,
        exact_search_nodes=summary["exact_search_nodes"],
        raw_origin_fallback_sources=summary["raw_origin_fallback_sources"],
        source_shard_bindings=bindings, is_global_completeness_certificate=False)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "source_shard_bindings"}, indent=2))


if __name__ == "__main__":
    main()
