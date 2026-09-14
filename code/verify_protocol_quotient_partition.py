"""Recompute the whole quotient-dimension distribution of a finite domain."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from protocol_domain import open_domain, read_header, read_pointings, source_key, write_native
from selected_spaces import read_selected


def check_distribution(domain, census, result):
    if (result["manifest_supports"] != domain["pointing_count"]
            or result["records_processed"] != domain["pointing_count"]
            or result["record_start"] != 0 or result["logical_dimension"] != 64
            or result["minimum_distance"] != 3 or result["manifest_length_filter"] != "protocol"
            or result["length_limit"] != domain["maximum_protocol_length"]
            or result["eligible_supports"] or result["isotropic_subspace_count"]):
        raise ValueError("Incomplete quotient-dimension count")
    rows = result["label_spaces_by_dimension"]
    counts = dict(rows)
    if len(rows) != len(counts) or sum(counts.values()) != domain["pointing_count"]:
        raise ValueError("Invalid complete quotient histogram")
    low, high = census["small_quotient_interval"]
    if (sum(n for d, n in rows if low <= d <= high) != census["small_quotient_supports"]
            or [[d, n] for d, n in rows if d > high] != census["quotient_selected_counts"]
            or any(d < low for d in counts)):
        raise ValueError("The complete histogram disagrees with the declared quotient partition")
    return rows


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/through48_logical_censuses.json")
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    census = json.loads(raw)
    domain_bytes = (root / census["source_domain"]).read_bytes()
    domain = json.loads(domain_bytes)
    if domain["schema"] != "finite-protocol-source-domain-v1" or args.workers < 1:
        raise ValueError("Invalid input domain or worker count")
    path = root / domain["pointing_data"]
    with path.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != domain["pointing_sha256"]:
            raise ValueError("Changed compact pointing domain")
    supports, bindings = read_selected(root, {source_key(row) for row in domain["sources"]})
    work = args.work_directory.resolve()
    work.mkdir(parents=True, exist_ok=True)
    manifest = work / "inputs.utsp"
    with open_domain(path) as stream, manifest.open("wb") as target:
        header = read_header(stream)
        write_native(target, header, domain["sources"], supports, read_pointings(stream, header))
    process = subprocess.run([str(args.native), "manifest-count", "--input", str(manifest),
        "--q", "64", "--minimum-distance", "3", "--workers", str(args.workers)],
        capture_output=True, text=True, check=True)
    histogram = check_distribution(domain, census, json.loads(process.stdout))
    result = dict(schema="complete-quotient-domain-verification-v1", status="pass",
        source_domain_sha256=hashlib.sha256(domain_bytes).hexdigest(),
        logical_census_sha256=hashlib.sha256(raw).hexdigest(),
        pointing_supports=domain["pointing_count"], quotient_dimension_counts=histogram,
        entire_quotient_partition_freshly_recomputed=True,
        omitted_larger_quotient_inputs=0, source_shard_bindings=bindings,
        is_global_completeness_certificate=False)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "source_shard_bindings"}, indent=2))


if __name__ == "__main__":
    main()
