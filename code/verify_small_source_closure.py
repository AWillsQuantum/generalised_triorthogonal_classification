"""Certify small source sectors using quotient, tensor and Pareto exclusions."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from logical_spaces import logical_label_space
from selected_spaces import read_selected
from space_codec import binary_rank
from space_lifts import evaluation_rows
from verify_large_quotient_domain import pointed_support
from verify_low_q_pruning import verify as verify_small_outputs
from verify_protocol_case_census import write_pointed_cases


def all_pointings(points, m, maximum_length):
    seen = set()
    records = []
    for parity in range(5):
        origins = range(1 << m) if parity < 3 else (-1,)
        for origin in origins:
            if parity == 1 and origin not in points or parity == 2 and origin in points:
                continue
            values, h = pointed_support(points, m, parity, origin)
            key = h, tuple(values)
            if len(values) > maximum_length or key in seen:
                continue
            seen.add(key)
            records.append(dict(index=len(records), points=values, ambient_dimension=h,
                                parity_case=parity, origin=origin))
    return records


def check_zero_census(result, q, count, distance):
    if (result["status"] != "complete" or result["matrix_scope"] != "full_projective"
            or result["logical_qubits"] != q or result["minimum_distance"] != distance
            or result["enumeration_mode"] != "raw_isotropic_subspaces"
            or result["support_range"] != dict(start=0, count=count, end_exclusive=count)):
        raise ValueError("Wrong finite exclusion census")
    stats = result["statistics"]
    if (stats["supports_processed"] != count or stats["quotient_dimension_filtered_supports"]
            or stats["marked_orbit_enumerated_supports"] or not stats["isotropic_subspace_statistics_exact"]
            or not stats["radical_statistics_exact"] or stats["nondegenerate_subspaces"]
            or result["pareto_protocols"]
            or sum(stats["radical_dimension_counts"].values()) != stats["isotropic_subspaces"]):
        raise ValueError("The full finite census does not exclude the claimed logical dimension")
    return dict(q=q, minimum_distance=distance, supports=count,
                isotropic_subspaces=stats["isotropic_subspaces"], nondegenerate_subspaces=0)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/through48_source_domain.json")
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    domain = json.loads(raw)
    if domain["schema"] != "finite-protocol-source-domain-v1" or args.workers < 1:
        raise ValueError("Invalid source domain or worker count")
    keys = {(row["c"], row["m"], row["index"]) for row in domain["separate_sources"]}
    supports, bindings = read_selected(root, keys)
    small_outputs = verify_small_outputs(root)
    work = args.work_directory.resolve()
    work.mkdir(parents=True, exist_ok=True)
    records = []
    for index, (key, points) in enumerate(sorted(supports.items())):
        c, m, _ = key
        cases = all_pointings(points, m, domain["maximum_protocol_length"])
        if not cases:
            raise ValueError("Empty finite source domain")
        minimum_n = min(len(row["points"]) for row in cases)
        minimum_h = min(row["ambient_dimension"] for row in cases)
        for output in small_outputs["witnesses"]:
            if not (output["incumbent_n"] <= minimum_n and output["incumbent_S"] < output["q"]+minimum_h):
                raise ValueError("The source's exact-distance-three small outputs are not all dominated")
        quotients, constants = [], []
        for row in cases:
            pts, h = row["points"], row["ambient_dimension"]
            d3 = logical_label_space(pts, h, 3)["quotient_dimension"]
            d4 = logical_label_space(pts, h, 4)["quotient_dimension"]
            quotients.append(dict(index=row["index"], d3=d3, d4=d4))
            if d4:
                stabilisers = evaluation_rows(pts, h, 1)[1:]
                if binary_rank([*stabilisers, (1 << len(pts))-1]) != h:
                    raise ValueError("A positive higher-distance quotient lacks the constant stabiliser")
                constants.append(row["index"])
        manifest = work / f"source_{index}.utsp"
        with manifest.open("wb") as stream:
            write_pointed_cases(stream, dict(cases=cases, maximum_protocol_length=domain["maximum_protocol_length"]))
        results = []
        for q, distance in ((5, 3), (3, 4)):
            path = work / f"source_{index}_q{q}_d{distance}.json"
            subprocess.run([str(args.native), "manifest-catalogue", "--input", str(manifest),
                "--output", str(path), "--q", str(q), "--minimum-distance", str(distance),
                "--support-workers", str(args.workers)], check=True, capture_output=True, text=True)
            results.append(check_zero_census(json.loads(path.read_bytes()), q, len(cases), distance))
        records.append(dict(space=dict(c=c, m=m, index=key[2]), pointings=len(cases),
            minimum_protocol_length=minimum_n, minimum_stabiliser_dimension=minimum_h,
            quotient_dimensions=quotients, positive_d4_constant_stabiliser_cases=constants,
            complete_raw_censuses=results))
    result = dict(schema="finite-small-source-closure-v1", status="pass",
        source_domain_sha256=hashlib.sha256(raw).hexdigest(),
        maximum_protocol_length=domain["maximum_protocol_length"], sources=records,
        all_origins_explicitly_enumerated=True, all_quotient_dimensions_recomputed=True,
        all_required_logical_censuses_freshly_recomputed=True,
        exact_d3_q_at_most4_dominated=True, all_q_at_least5_excluded_at_distance3=True,
        all_outputs_excluded_at_distance_at_least4=True,
        zero_column_extensions_pareto_redundant=True, source_shard_bindings=bindings,
        low_q_frontier_sha256=small_outputs["frontier_sha256"],
        low_q_census_sha256=small_outputs["q5_census_sha256"],
        is_global_completeness_certificate=False)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps(dict(status="pass", sources=len(records),
        pointings=sum(row["pointings"] for row in records), logical_censuses=2*len(records)), indent=2))


if __name__ == "__main__":
    main()
