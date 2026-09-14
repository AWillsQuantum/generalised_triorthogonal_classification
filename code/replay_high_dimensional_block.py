"""Recompute one finite sector block directly from the compact space catalogue."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from compact_to_native import write_base_manifest
from protocol_checks import check_matrix, distance_up_to
from selected_spaces import read_selected


def replay(root, data_path, block_index, binary, work, workers):
    raw = data_path.read_bytes()
    data = json.loads(raw)
    block = data["blocks"][block_index]
    if block["index"] != block_index or workers < 1:
        raise ValueError("Invalid block selection")
    c, m = block["c"], block["m"]
    first, end = block["source_interval"]
    keys = [(c, m, index) for index in range(first, end)]
    supports, bindings = read_selected(root, keys)
    work.mkdir(parents=True, exist_ok=True)
    source = work / "base.utsp"
    with source.open("wb") as stream:
        write_base_manifest(stream, (supports[key] for key in keys), c=c, m=m,
                            first_index=first, count=end-first, length_limit=54)
    reduced = work / "pointed.utsp"
    process = subprocess.run([str(binary), "manifest-origin-reduce", "--input", str(source),
                              "--output", str(reduced), "--node-limit", "50000000",
                              "--workers", str(workers)], check=True, capture_output=True, text=True)
    origins = json.loads(process.stdout)
    if (origins["parity_counts"] != block["pointing_counts"]
            or origins["reduced_supports"] != block["reduced_supports"]):
        raise ValueError("Recomputed pointing counts differ")
    censuses = []
    for branch in block["branches"]:
        q, floor = branch["q"], branch["minimum_distance"]
        output = work / f"q{q}.json"
        subprocess.run([str(binary), "manifest-catalogue", "--input", str(reduced),
                        "--output", str(output), "--q", str(q), "--minimum-distance", str(floor),
                        "--support-workers", str(workers)], check=True, capture_output=True, text=True)
        result = json.loads(output.read_bytes())
        if (result["status"] != "complete" or result["logical_qubits"] != q
                or result["minimum_distance"] != floor or result["matrix_scope"] != "full_projective"
                or result["enumeration_mode"] != "raw_isotropic_subspaces"
                or result["support_range"] != dict(start=0, count=block["reduced_supports"],
                                                   end_exclusive=block["reduced_supports"])):
            raise ValueError("Incomplete native replay")
        stats = result["statistics"]
        if (not stats["isotropic_subspace_statistics_exact"] or not stats["radical_statistics_exact"]
                or any(stats[key] != value for key, value in branch["statistics"].items())):
            raise ValueError("Recomputed finite census differs")
        actual_metrics = set()
        for witness in result["pareto_protocols"]:
            rows = witness["generator_rows"]
            _, columns = check_matrix(rows, q)
            if distance_up_to(columns, q) != (witness["d_Z"], witness["error_coefficient"]):
                raise ValueError("Fresh native witness failed the independent distance test")
            actual_metrics.add((witness["protocol_length_n"], witness["space_footprint_S"], witness["d_Z"]))
        expected_metrics = {(data["witnesses"][i]["n"], data["witnesses"][i]["S"], data["witnesses"][i]["d_Z"])
                            for i in branch["witness_indices"]}
        if actual_metrics != expected_metrics:
            raise ValueError("Recomputed finite Pareto metric set differs")
        censuses.append(dict(q=q, minimum_distance=floor, statistics=branch["statistics"]))
    return dict(schema="finite-high-dimensional-block-replay-v1", status="pass",
                sector_data_sha256=hashlib.sha256(raw).hexdigest(), block=block_index,
                c=c, m=m, source_interval=[first, end], source_shards=bindings,
                pointing_counts=origins["parity_counts"], censuses=censuses,
                is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/length54_high_dimensional.json")
    parser.add_argument("--block", type=int, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = replay(root, args.input, args.block, args.native.resolve(), args.work_directory.resolve(), args.workers)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps({key: value for key, value in result.items() if key not in ("source_shards", "censuses")}, indent=2))


if __name__ == "__main__":
    main()
