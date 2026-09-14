"""Enumerate a finite protocol sector using shared distance-filtered quotients."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from compact_to_native import write_base_manifest
from protocol_checks import LogicalTensor, check_matrix, distance_up_to
from selected_spaces import read_selected


def mathematical_statistics(stats):
    fields = ("source_spaces", "supports_processed", "quotient_dimension_filtered_supports",
              "eligible_supports", "raw_enumerated_supports", "marked_orbit_enumerated_supports",
              "isotropic_subspace_statistics_exact", "isotropic_subspaces", "radical_statistics_exact",
              "nondegenerate_subspaces", "marked_code_orbits", "canonical_output_orbits",
              "radical_dimension_counts")
    return {key: stats[key] for key in fields}


def check_raw_census(result, q, distance, count):
    stats = mathematical_statistics(result["statistics"])
    if (result["status"] != "complete" or result["logical_qubits"] != q
            or result["minimum_distance"] != distance or result["matrix_scope"] != "full_projective"
            or result["enumeration_mode"] != "raw_isotropic_subspaces"
            or result["support_range"] != dict(start=0, count=count, end_exclusive=count)
            or stats["supports_processed"] != count or stats["quotient_dimension_filtered_supports"]
            or stats["marked_orbit_enumerated_supports"] or not stats["isotropic_subspace_statistics_exact"]
            or not stats["radical_statistics_exact"]
            or stats["raw_enumerated_supports"] != stats["eligible_supports"]
            or sum(stats["radical_dimension_counts"].values())+stats["nondegenerate_subspaces"] != stats["isotropic_subspaces"]):
        raise ValueError("Incomplete or inconsistent direct finite census")
    witnesses = []
    for row in result["pareto_protocols"]:
        masks, columns = check_matrix(row["generator_rows"], q)
        tensor = LogicalTensor(masks[:q])
        actual = distance_up_to(columns, q, maximum=row["d_Z"])
        if (not tensor.intrinsic() or actual != (row["d_Z"], row["error_coefficient"])
                or row["d_Z"] < distance or row["protocol_length_n"] != len(columns)
                or row["space_footprint_S"] != len(masks)):
            raise ValueError("Native witness failed independent matrix or distance checks")
        witnesses.append(dict(q=q, n=len(columns), S=len(masks), d_Z=row["d_Z"],
            error_coefficient=row["error_coefficient"], output_key_words=row["output"]["canonical_key_words"],
            generator_matrix_rows=row["generator_rows"]))
    if len({tuple(row["output_key_words"]) for row in witnesses}) != stats["canonical_output_orbits"]:
        raise ValueError("Not every observed output has a finite witness")
    return stats, witnesses


def enumerate_block(root, domain_path, block_index, native, work, workers):
    raw = domain_path.read_bytes()
    data = json.loads(raw)
    if data["schema"] != "finite-protocol-source-block-domain-v1" or workers < 1:
        raise ValueError("Invalid source domain or worker count")
    block = data["blocks"][block_index]
    if block["index"] != block_index or block_index < 0:
        raise ValueError("Invalid source block")
    c, m = block["c"], block["m"]
    first, end = block["source_interval"]
    keys = [(c, m, i) for i in range(first, end)]
    supports, bindings = read_selected(root, keys)
    work.mkdir(parents=True, exist_ok=True)
    base, pointed = work / "base.utsp", work / "pointed.utsp"
    with base.open("wb") as stream:
        write_base_manifest(stream, (supports[key] for key in keys), c=c, m=m,
            first_index=first, count=end-first, length_limit=data["maximum_protocol_length"])
    def call(*arguments):
        result = subprocess.run([str(native), *map(str, arguments)], capture_output=True, text=True)
        if result.returncode:
            raise ValueError("Native mathematical calculation failed: "+result.stderr.strip())
        return json.loads(result.stdout)
    origin = call("manifest-origin-reduce", "--input", base, "--output", pointed,
        "--minimum-protocol-length", data["minimum_protocol_length"],
        "--maximum-protocol-length", data["maximum_protocol_length"],
        "--node-limit", 50000000, "--workers", workers)
    if (origin["source_spaces"] != end-first
            or origin["reduced_supports"] != block["expected_pointed_supports"]):
        raise ValueError("Fresh origin quotient differs from the bound source-domain count")
    caches = {d: work / f"d{d}.utslc" for d in (3, 4)}
    manifests = {d: work / f"d{d}.utsp" for d in (3, 4)}
    # The word-packed source-Schur equations support affine dimensions <=10.
    # Larger sources still share the directly constructed label cache.
    schur = ["--source-schur"] if m <= 10 else []
    labels = call("manifest-initial-label-caches", "--input", pointed,
        "--distance-three-cache", caches[3], "--distance-three-manifest", manifests[3],
        "--distance-four-cache", caches[4], "--distance-four-manifest", manifests[4],
        "--workers", workers, *schur)
    if labels["status"] != "complete" or labels["input_supports"] != origin["reduced_supports"]:
        raise ValueError("Incomplete shared logical-label quotient computation")
    label_data = {d: labels["distance_three" if d == 3 else "distance_four"] for d in (3, 4)}
    for d, minimum in ((3, 5), (4, 1)):
        row = label_data[d]
        if (row["minimum_distance"] != d or row["minimum_selected_quotient_dimension"] != minimum
                or sum(n for _, n in row["quotient_dimension_histogram"]) != origin["reduced_supports"]
                or sum(n for dim, n in row["quotient_dimension_histogram"] if dim >= minimum) != row["supports"]):
            raise ValueError("Shared quotient filter does not cover the declared distance domain")
    branches = []
    for specification in data["aggregate_censuses"]:
        q, distance = specification["q"], specification["minimum_distance"]
        minimum = label_data[distance]["minimum_selected_quotient_dimension"]
        if q < minimum:
            raise ValueError("The shared quotient filter excludes eligible logical subspaces")
        output = work / f"q{q}.json"
        subprocess.run([str(native), "manifest-catalogue", "--input", str(manifests[distance]),
            "--label-space-cache", str(caches[distance]), "--output", str(output),
            "--q", str(q), "--minimum-distance", str(distance), "--support-workers", str(workers)],
            check=True, capture_output=True, text=True)
        stats, witnesses = check_raw_census(json.loads(output.read_bytes()), q, distance, label_data[distance]["supports"])
        eligible = sum(n for dim, n in label_data[distance]["quotient_dimension_histogram"] if dim >= q)
        if stats["eligible_supports"] != eligible or stats["source_spaces"] != end-first:
            raise ValueError("Logical census eligibility disagrees with the complete quotient histogram")
        branches.append(dict(q=q, minimum_distance=distance, minimum_input_quotient_dimension=minimum,
                             statistics=stats, witnesses=witnesses))
    return dict(schema="finite-protocol-source-block-census-v1", status="pass",
        source_domain_sha256=hashlib.sha256(raw).hexdigest(), index=block_index, c=c, m=m,
        source_interval=[first, end], minimum_protocol_length=data["minimum_protocol_length"],
        maximum_protocol_length=data["maximum_protocol_length"],
        pointed_supports=origin["reduced_supports"], parity_counts=origin["parity_counts"],
        raw_origins=origin["raw_origins"], origin_orbits=origin["origin_orbits"],
        source_shard_bindings=bindings,
        distance_three_quotient_profile=label_data[3]["quotient_dimension_histogram"],
        distance_four_quotient_profile=label_data[4]["quotient_dimension_histogram"],
        branches=branches, all_censuses_freshly_recomputed=True,
        is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/length50_source_domain.json")
    parser.add_argument("--block", type=int, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = enumerate_block(root, args.input, args.block, args.native.resolve(), args.work_directory.resolve(), args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix+".partial")
    temporary.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    temporary.replace(args.output)
    print(json.dumps(dict(status="pass", block=args.block, sources=result["source_interval"][1]-result["source_interval"][0],
        pointings=result["pointed_supports"], witnesses=sum(len(b["witnesses"]) for b in result["branches"])), indent=2))


if __name__ == "__main__":
    main()
