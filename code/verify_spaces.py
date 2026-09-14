"""Check included compact shards and their ordered mathematical contents."""

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

from read_spaces import open_shard
from space_codec import BITMAP, read_exact, read_header, read_records, validate_unital_support


def checked_content(stream, header, validity, reference=False, batch_size=4096):
    """Yield canonical point bytes after decoding with bounded memory."""
    if batch_size < 1:
        raise ValueError("Batch size must be positive")
    if reference or header.codec != BITMAP or (validity and not 0 < header.c <= 64):
        width = (header.m + 7) // 8
        for points in read_records(stream, header):
            if validity:
                validate_unital_support(points, header.m)
            yield 1, b"".join(x.to_bytes(width, "little") for x in points)
        return
    from space_batch import decode_bitmaps, ordered_point_bytes, validate_batch
    batch_size = min(batch_size, max(1, (32 * 1024 * 1024) // (1 << header.m)))
    remaining = header.count
    while remaining:
        number = min(batch_size, remaining)
        block = read_exact(stream, number * header.width)
        points = decode_bitmaps(block, header.m, header.c)
        if validity:
            validate_batch(points, header.m, header.c)
        yield number, ordered_point_bytes(points, header.m)
        remaining -= number
    if stream.read(1):
        raise ValueError("Trailing support shard data")


def verify(root, validity=False, reference=False, batch_size=4096, maximum_length=None,
           require_complete=False):
    root = root.resolve()
    sectors = defaultdict(list)
    files = set()
    count = 0
    total_bytes = 0
    if maximum_length is not None and maximum_length < 1:
        raise ValueError("Maximum length must be positive")
    for certificate_path in sorted((root / "certificates/support_encoding").glob("*.json")):
        certificate = json.loads(certificate_path.read_bytes())
        if maximum_length is not None and certificate["c"] > maximum_length:
            continue
        if certificate["status"] != "pass":
            raise ValueError("Unsuccessful encoding certificate")
        m, c = certificate["m"], certificate["c"]
        first = certificate["first_index"]
        seen = 0
        overall = hashlib.sha256()
        for row in certificate["shards"]:
            path = (root / row["path"]).resolve()
            if not path.is_relative_to(root) or path in files:
                raise ValueError("External or repeated shard path")
            files.add(path)
            if path.stat().st_size != row["bytes"]:
                raise ValueError("Incorrect compressed byte count")
            with path.open("rb") as raw:
                if hashlib.file_digest(raw, "sha256").hexdigest() != row["sha256"]:
                    raise ValueError("Incorrect compressed digest")
            content = hashlib.sha256()
            with open_shard(path) as stream:
                header = read_header(stream)
                if (header.m, header.c, header.first_index, header.count) != (
                        m, c, first + seen, row["count"]):
                    raise ValueError("Inconsistent shard header")
                if header.first_index != row["first_index"]:
                    raise ValueError("Incorrect shard interval")
                for number, encoded in checked_content(stream, header, validity, reference, batch_size):
                    content.update(encoded)
                    overall.update(encoded)
                    seen += number
            if content.hexdigest() != row["support_sha256"]:
                raise ValueError("Ordered support digest mismatch")
            total_bytes += row["bytes"]
        if seen != certificate["count"] or overall.hexdigest() != certificate["ordered_support_sha256"]:
            raise ValueError("Certificate count or combined support digest mismatch")
        sectors[c, m].append((first, first + seen))
        count += seen
    actual_files = set()
    for path in (root / "data/spaces").rglob("*.utspace.zst"):
        with open_shard(path) as stream:
            header = read_header(stream)
        if maximum_length is None or header.c <= maximum_length:
            actual_files.add(path.resolve())
    if not files or actual_files != files:
        raise ValueError("Missing or uncertified included shards")
    coverage = []
    for (c, m), intervals in sorted(sectors.items()):
        intervals.sort()
        if any(a[1] > b[0] for a, b in zip(intervals, intervals[1:])):
            raise ValueError("Overlapping record intervals")
        coverage.append(dict(c=c, m=m, included_intervals=intervals,
                             count=sum(b-a for a, b in intervals)))
    index_digest = None
    if require_complete:
        if maximum_length is not None:
            raise ValueError("Complete-catalogue audit cannot be length-filtered")
        raw_index = (root / "data/spaces/index.json").read_bytes()
        index = json.loads(raw_index)
        if index.get("data_complete") is not True or index["expected_classes"] != count:
            raise ValueError("The full declared catalogue is not present")
        index_digest = hashlib.sha256(raw_index).hexdigest()
        expected = {}
        indexed_certificates = set()
        for length in index["lengths"]:
            for sector in length["sectors"]:
                key = length["c"], sector["m"]
                if key in expected or sector["data_complete"] is not True:
                    raise ValueError("Repeated or incomplete catalogue sector")
                expected[key] = sector["expected_classes"]
                cursor = 0
                for first, end in sorted(sectors.get(key, [])):
                    if first != cursor:
                        raise ValueError("A full catalogue sector has an interval gap")
                    cursor = end
                if cursor != expected[key]:
                    raise ValueError("Catalogue sector count does not match the declared index")
                for binding in sector["encoding_certificates"]:
                    path = (root / binding["path"]).resolve()
                    if not path.is_relative_to(root) or path in indexed_certificates:
                        raise ValueError("External or repeated indexed certificate")
                    if hashlib.sha256(path.read_bytes()).hexdigest() != binding["sha256"]:
                        raise ValueError("The indexed encoding certificate changed")
                    indexed_certificates.add(path)
        actual_certificates = {p.resolve() for p in (root / "certificates/support_encoding").glob("*.json")}
        if (indexed_certificates != actual_certificates or set(sectors)-set(expected)
                or sum(expected.values()) != count or index["compressed_bytes"] != total_bytes):
            raise ValueError("The complete data index and audited files disagree")
    return dict(schema="triorthogonal-included-space-verification-v1", status="pass",
                shards=len(files), representatives=count, compressed_bytes=total_bytes,
                validity_checked_for_all_included_representatives=validity,
                decoder="reference" if reference else "batched_bitmap_with_reference_colex",
                maximum_length=maximum_length,
                all_declared_catalogue_records_verified=require_complete,
                catalogue_index_sha256=index_digest,
                classification_completeness_established=False,
                coverage=coverage)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--validity", action="store_true")
    parser.add_argument("--reference", action="store_true",
                        help="Use the independent scalar decoder and validity checker")
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--maximum-length", type=int,
                        help="Audit only completed space lengths at or below this bound")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-complete", action="store_true",
                        help="Require every indexed sector, record interval and certificate")
    args = parser.parse_args()
    report = verify(args.root, args.validity, args.reference, args.batch_size,
                    args.maximum_length, args.require_complete)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
