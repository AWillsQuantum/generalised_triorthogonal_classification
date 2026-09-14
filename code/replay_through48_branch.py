"""Reproduce one finite logical census in the through-48 source domain."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from protocol_domain import open_domain, read_header, read_pointings, source_key, write_native
from selected_spaces import read_selected
from verify_protocol_case_census import write_pointed_cases


def metric(row):
    if "output" in row:
        return (tuple(row["output"]["canonical_key_words"]), row["d_Z"],
                row["protocol_length_n"], row["space_footprint_S"])
    return tuple(row["output_key_words"]), row["d_Z"], row["n"], row["S"]


def check_census(branch, actual):
    if (actual["status"] != "complete" or actual["matrix_scope"] != "full_projective"
            or actual["logical_qubits"] != branch["q"] or actual["minimum_distance"] != 3
            or actual["length_limit"] != 48 or actual["enumeration_mode"] != branch["enumeration_mode"]):
        raise ValueError("Wrong logical census scope")
    counts = ("eligible_supports", "raw_enumerated_supports", "marked_orbit_enumerated_supports",
              "isotropic_subspace_statistics_exact", "isotropic_subspaces", "radical_statistics_exact",
              "nondegenerate_subspaces", "marked_code_orbits", "canonical_output_orbits", "radical_dimension_counts")
    for name in counts:
        if actual["statistics"][name] != branch["statistics"][name]:
            raise ValueError("Recomputed census differs on "+name)
    if (actual["statistics"]["supports_processed"]-actual["statistics"]["quotient_dimension_filtered_supports"]
            != branch["selected_supports"]):
        raise ValueError("Wrong filtered input count")
    if {metric(row) for row in actual["pareto_protocols"]} != {metric(row) for row in branch["witnesses"]}:
        raise ValueError("Recomputed frontier metrics differ")
    return {name: actual["statistics"][name] for name in counts}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--branch", type=int, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    path = root / "data/protocol_sectors/through48_logical_censuses.json"
    raw = path.read_bytes()
    data = json.loads(raw)
    if not 0 <= args.branch < len(data["branches"]) or args.workers < 1:
        parser.error("Invalid branch index or worker count")
    branch = data["branches"][args.branch]
    domain_path = root / data["source_domain"]
    domain_bytes = domain_path.read_bytes()
    domain = json.loads(domain_bytes)
    work = args.work_directory.resolve()
    work.mkdir(parents=True, exist_ok=True)
    manifest, output = work / "inputs.utsp", work / "census.json"
    low, high = branch["quotient_dimension_interval"]
    options = []
    if low == 0:
        supports, _ = read_selected(root, {source_key(row) for row in domain["sources"]})
        with open_domain(root / domain["pointing_data"]) as stream, manifest.open("wb") as target:
            header = read_header(stream)
            write_native(target, header, domain["sources"], supports, read_pointings(stream, header))
        options += ["--maximum-quotient-dimension", str(high)]
    else:
        selected = [row for row in data["selected_quotients"] if low <= row["quotient_dimension"] <= high]
        cases = [dict(index=i, points=row["points"], ambient_dimension=row["ambient_dimension"])
                 for i, row in enumerate(selected)]
        with manifest.open("wb") as stream:
            write_pointed_cases(stream, dict(cases=cases, maximum_protocol_length=48))
    if branch["enumeration_mode"] == "complete_marked_code_orbits":
        options += ["--marked-orbits", "--workers", str(args.workers)]
    elif branch["enumeration_mode"] == "hybrid_raw_and_marked_code_orbits":
        options += ["--hybrid-orbits", "--raw-max-quotient-dimension", "14", "--workers", str(args.workers)]
    subprocess.run([str(args.native), "manifest-catalogue", "--input", str(manifest), "--output", str(output),
        "--q", str(branch["q"]), "--minimum-distance", "3", "--support-workers", str(args.workers), *options],
        check=True, capture_output=True, text=True)
    counts = check_census(branch, json.loads(output.read_bytes()))
    result = dict(schema="finite-through48-logical-branch-replay-v1", status="pass",
        branch_index=args.branch, q=branch["q"], quotient_dimension_interval=[low, high],
        selected_supports=branch["selected_supports"], logical_census_sha256=hashlib.sha256(raw).hexdigest(),
        source_domain_sha256=hashlib.sha256(domain_bytes).hexdigest(),
        complete_census_counts=counts, complete_frontier_metrics_agree=True,
        logical_census_freshly_recomputed=True, is_global_completeness_certificate=False)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
