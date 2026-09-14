"""Reconstruct a length-52 source interval and enumerate a finite logical branch."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from compact_to_native import write_base_manifest
from enumerate_source_block import mathematical_statistics
from selected_spaces import read_selected
from verify_length52_small_blocks import check_branch


def replay(root, domain_path, block_index, q, native, work, workers, maximum_quotient=None):
    raw = domain_path.read_bytes()
    domain = json.loads(raw)
    if domain["schema"] != "finite-length52-protocol-source-domain-v1" or q not in range(1, 6) or workers < 1:
        raise ValueError("Invalid source domain or logical branch")
    source = domain["main_blocks"][block_index]
    if source["index"] != block_index or block_index < 0:
        raise ValueError("Incorrect source index")
    c, m = source["c"], source["m"]
    first, end = source["source_interval"]
    keys = [(c, m, i) for i in range(first, end)]
    points, bindings = read_selected(root, keys)
    work.mkdir(parents=True, exist_ok=True)
    base, pointed = work / "base.utsp", work / "pointed.utsp"
    with base.open("wb") as stream:
        write_base_manifest(stream, (points[key] for key in keys), c=c, m=m,
            first_index=first, count=end-first, length_limit=52)
    def call(*arguments):
        result = subprocess.run([str(native), *map(str, arguments)], capture_output=True, text=True)
        if result.returncode:
            raise ValueError("Native mathematical calculation failed: "+result.stderr.strip())
        return json.loads(result.stdout)
    origin = call("manifest-origin-reduce", "--input", base, "--output", pointed,
        "--minimum-protocol-length", 51, "--maximum-protocol-length", 52, "--workers", workers,
        *(["--raw-origins"] if m <= 10 else ["--node-limit", 50000000]))
    if origin["source_spaces"] != end-first or origin["source_start"] != 0:
        raise ValueError("Incomplete pointing source interval")
    manifests = {d: work / f"d{d}.utsp" for d in (3, 4)}
    caches = {d: work / f"d{d}.utslc" for d in (3, 4)}
    labels = call("manifest-initial-label-caches", "--input", pointed,
        "--distance-three-manifest", manifests[3], "--distance-three-cache", caches[3],
        "--distance-four-manifest", manifests[4], "--distance-four-cache", caches[4],
        "--workers", workers, *(["--source-schur"] if m <= 10 else []))
    if labels["status"] != "complete" or labels["input_supports"] != origin["reduced_supports"]:
        raise ValueError("Incomplete logical quotient construction")
    floor = 3 if q == 5 else 4
    selected = labels["distance_three" if q == 5 else "distance_four"]
    minimum = 5 if q == 5 else 1
    profile = selected["quotient_dimension_histogram"]
    if (selected["minimum_selected_quotient_dimension"] != minimum
            or sum(n for _, n in profile) != origin["reduced_supports"]
            or sum(n for d, n in profile if d >= minimum) != selected["supports"]):
        raise ValueError("The selected quotient does not cover the logical domain")
    high = (14 if q == 5 else 64) if maximum_quotient is None else maximum_quotient
    if not q <= high <= 64:
        raise ValueError("Invalid maximum quotient dimension")
    output = work / "frontier.json"
    call("manifest-catalogue", "--input", manifests[floor], "--label-space-cache", caches[floor],
        "--output", output, "--q", q, "--minimum-distance", floor,
        "--minimum-quotient-dimension", minimum, "--maximum-quotient-dimension", high,
        "--support-workers", workers, "--emit-positive-supports")
    result = json.loads(output.read_bytes())
    if (result["status"] != "complete" or result["matrix_scope"] != "full_projective"
            or result["enumeration_mode"] != "raw_isotropic_subspaces" or result["length_limit"] != 52
            or result["manifest_length_filter"] != "protocol"):
        raise ValueError("Incorrect direct finite enumeration mode")
    branch = dict(q=q, minimum_distance=floor, enumeration_mode=result["enumeration_mode"],
        quotient_dimension_interval=result["quotient_dimension_interval"], support_range=result["support_range"],
        input_quotient_profile=[[d, n] for d, n in profile if d >= minimum],
        statistics=mathematical_statistics(result["statistics"]),
        positive_supports=[dict(selected_support_index=row["manifest_support_index"],
                               output_keys=[key["canonical_key_words"] for key in row["output_keys"]])
                           for row in result["positive_supports"]],
        witnesses=[dict(q=q, n=row["protocol_length_n"], S=row["space_footprint_S"], d_Z=row["d_Z"],
                        error_coefficient=row["error_coefficient"], output_key_words=row["output"]["canonical_key_words"],
                        generator_matrix_rows=row["generator_rows"])
                   for row in result["pareto_protocols"]])
    check_branch(branch, dict(profile), end-first)
    return dict(schema="finite-length52-branch-replay-v1", status="pass", c=c, m=m, index=block_index,
        source_interval=[first, end], source_domain_sha256=hashlib.sha256(raw).hexdigest(),
        source_shard_bindings=bindings, protocol_length_interval=[51, 52],
        pointing_supports=origin["reduced_supports"], parity_counts=origin["parity_counts"],
        distance_three_quotient_profile=labels["distance_three"]["quotient_dimension_histogram"],
        distance_four_quotient_profile=labels["distance_four"]["quotient_dimension_histogram"],
        branch=branch, all_censuses_freshly_recomputed=True, is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/length52_source_domain.json")
    parser.add_argument("--block", type=int, required=True)
    parser.add_argument("--q", type=int, required=True)
    parser.add_argument("--maximum-quotient-dimension", type=int)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = replay(root, args.input, args.block, args.q, args.native.resolve(), args.work_directory.resolve(),
                    args.workers, args.maximum_quotient_dimension)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps(dict(status="pass", block=args.block, q=args.q, statistics=result["branch"]["statistics"]), indent=2))


if __name__ == "__main__":
    main()
