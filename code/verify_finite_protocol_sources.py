"""Verify or reproduce complete pointed-source blocks and their Pareto exclusions."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess

from compact_to_native import write_base_manifest
from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from selected_spaces import read_selected
from space_codec import validate_unital_support


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source_keys(data):
    keys = set()
    for block in data["blocks"]:
        first, end = block["source_interval"]
        require(type(first) is int and type(end) is int and 0 <= first < end,
                "Invalid source interval")
        selected = {(block["c"], block["m"], i) for i in range(first, end)}
        require(not keys & selected and len(selected) == block["source_spaces"], "Overlapping source blocks")
        keys.update(selected)
    require(len(keys) == data["source_spaces"], "Incomplete source count")
    return keys


def raw_pointing_counts(points, m):
    """Count all five cases, quotienting only identical translated point sets."""
    points = frozenset(points)
    origin = min(points)
    periods = [point ^ origin for point in points
               if all((value ^ point ^ origin) in points for value in points)]
    t = len(periods)
    require(t > 0 and not t & (t - 1), "Translation periods do not form a binary subspace")
    require(len(points) % t == 0, "The support is not a union of translation-period cosets")
    return [(1 << m) // t, len(points) // t, ((1 << m) - len(points)) // t, 1, 1]


def census_metric(row):
    if "output" in row:
        return (tuple(int(w, 0) for w in row["output"]["canonical_key_words"]),
                row["d_Z"], row["protocol_length_n"], row["space_footprint_S"])
    return tuple(int(w, 0) for w in row["output_key_words"]), row["d_Z"], row["n"], row["S"]


def native_json(native, arguments, receipt):
    result = subprocess.run([str(native), *map(str, arguments)], capture_output=True, text=True, check=True)
    data = json.loads(result.stdout)
    receipt.write_text(result.stdout, encoding="ascii")
    return data


def replay_block(native, work, block, supports, workers, independent_labels):
    work.mkdir(parents=True, exist_ok=True)
    base, pointed, cache = [work / name for name in ("base.utsp", "pointed.utsp", "labels.bin")]
    c, m = block["c"], block["m"]
    first, end = block["source_interval"]
    with base.open("wb") as stream:
        write_base_manifest(stream, [supports[c, m, i] for i in range(first, end)],
                            c=c, m=m, first_index=first, count=end-first, length_limit=54)
    origin = native_json(native, ["manifest-origin-reduce", "--input", base, "--output", pointed,
        "--raw-origins", "--maximum-protocol-length", 54, "--workers", workers], work / "origins.json")
    require(origin["reduced_supports"] == block["pointing_supports"]
            and origin["parity_counts"] == block["parity_counts"], "The complete pointing domain differs")
    profile = native_json(native, ["manifest-label-cache", "--input", pointed, "--output", cache,
        "--minimum-distance", 3, "--workers", workers, *(["--source-schur"] if m <= 10 else [])],
        work / "label_profile.json")
    require(profile["quotient_dimension_histogram"] == block["quotient_dimension_histogram"],
            "The complete logical quotient profile differs")
    if independent_labels and m <= 10:
        reference = work / "reference_labels.bin"
        native_json(native, ["manifest-label-cache", "--input", pointed, "--output", reference,
            "--minimum-distance", 3, "--workers", workers], work / "reference_profile.json")
        native_json(native, ["manifest-label-cache-differential", "--input", pointed,
            "--reference-cache", reference, "--candidate-cache", cache,
            "--minimum-distance", 3, "--workers", workers], work / "label_differential.json")
    for branch in block["branches"]:
        q = branch["q"]
        output = work / f"q{q}.json"
        native_json(native, ["manifest-catalogue", "--input", pointed, "--output", output,
            "--label-space-cache", cache, "--q", q, "--minimum-distance", 3,
            "--support-workers", workers], work / f"q{q}_summary.json")
        result = json.loads(output.read_bytes())
        require(result["status"] == "complete" and result["logical_qubits"] == q
                and result["minimum_distance"] == 3 and result["length_limit"] == 54
                and result["matrix_scope"] == "full_projective"
                and result["enumeration_mode"] == "raw_isotropic_subspaces", "Incorrect reproduced census scope")
        require(all(result["statistics"][k] == v for k, v in branch["statistics"].items()),
                "Recomputed logical census statistics disagree")
        require({census_metric(r) for r in result["pareto_protocols"]}
                == {census_metric(r) for r in branch["witnesses"]}, "Recomputed finite frontier differs")


def verify(root, native=None, work=None, workers=1, independent_labels=False):
    root = Path(root).resolve()
    path = root / "data/protocol_sectors/finite_shorter_source_blocks.json"
    raw = path.read_bytes()
    data = json.loads(raw)
    require(data["schema"] == "finite-pointed-source-block-censuses-v1"
            and data["maximum_protocol_length"] == 54 and data["minimum_distance"] == 3
            and data["matrix_scope"] == "full_projective" and data["output_equivalence"] == "CNOT+S"
            and data["pareto_objectives"] == ["n", "S"] and data["exact_distance_is_part_of_key"],
            "Incorrect finite-source classification scope")
    supports, bindings = read_selected(root, source_keys(data))
    require(bindings == data["source_shard_bindings"], "Changed finite-source catalogue bindings")
    frontier_path = root / "data/protocols/pareto_frontier.json"
    frontier = json.loads(frontier_path.read_bytes())
    outputs = {(r["q"], tuple(int(w, 0) for w in r["tensor_key_words"])): r["output_id"] for r in frontier["outputs"]}
    tensors = {r["output_id"]: LogicalTensor(tuple(int(s, 2) for s in r["generator_matrix_rows"][:r["q"]]))
               for r in frontier["protocols"]}
    total, witnesses, metrics = Counter(), [], set()
    intervals = []
    for index, block in enumerate(data["blocks"]):
        require(block["index"] == index and block["pointing_mode"] == "all_origins", "Invalid finite block identity")
        c, m = block["c"], block["m"]
        require(c + 1 <= 48, "This finite shorter-source cover exceeds protocol length 48")
        first, end = block["source_interval"]
        counts = [0] * 5
        for i in range(first, end):
            points = supports[c, m, i]
            validate_unital_support(points, m)
            counts = [a + b for a, b in zip(counts, raw_pointing_counts(points, m))]
        require(counts == block["parity_counts"] and sum(counts) == block["pointing_supports"],
                "The five pointing cases do not exhaust the source block")
        histogram = dict(block["quotient_dimension_histogram"])
        maximum = block["maximum_quotient_dimension"]
        require(len(histogram) == len(block["quotient_dimension_histogram"])
                and sum(histogram.values()) == sum(counts) and max(histogram) == maximum,
                "Invalid logical quotient profile")
        require([branch["q"] for branch in block["branches"]] == list(range(1, min(maximum, 8) + 1)),
                "Incomplete logical-dimension census")
        for branch in block["branches"]:
            q, stats = branch["q"], branch["statistics"]
            require(branch["enumeration_mode"] == "raw_isotropic_subspaces"
                    and stats["isotropic_subspace_statistics_exact"] and stats["radical_statistics_exact"]
                    and not stats["marked_orbit_enumerated_supports"] and not stats["quotient_dimension_filtered_supports"],
                    "The finite logical census is not exhaustive raw enumeration")
            require(stats["supports_processed"] == sum(counts)
                    and stats["eligible_supports"] == sum(n for d, n in histogram.items() if d >= q)
                    and stats["raw_enumerated_supports"] == stats["eligible_supports"]
                    and stats["nondegenerate_subspaces"] + sum(stats["radical_dimension_counts"].values()) == stats["isotropic_subspaces"],
                    "Inconsistent complete logical census")
            if q == 8 and maximum > 8:
                require(stats["isotropic_subspaces"] == 0, "The higher-dimensional isotropic containment obligation is open")
            total["logical_censuses"] += 1
            total["raw_isotropic_subspaces"] += stats["isotropic_subspaces"]
            total["raw_nondegenerate_subspaces"] += stats["nondegenerate_subspaces"]
            for position, row in enumerate(branch["witnesses"]):
                masks, columns = check_matrix(row["generator_matrix_rows"], q)
                tensor = LogicalTensor(masks[:q])
                identifier = outputs.get((q, tuple(int(w, 0) for w in row["output_key_words"])))
                require(identifier is not None and tensor.intrinsic(), "Unknown or degenerate finite output")
                basis = find_equivalence_basis(tensor, tensors[identifier])
                distance, coefficient = distance_up_to(columns, q, max(3, row["d_Z"]))
                require(basis is not None and (row["q"], row["n"], row["S"], row["d_Z"], row["error_coefficient"])
                        == (q, len(columns), len(masks), distance, coefficient), "Invalid finite protocol witness")
                dominators = [r["index"] for r in frontier["protocols"] if r["output_id"] == identifier
                    and r["d_Z"] == distance and r["n"] <= row["n"] and r["S"] <= row["S"]
                    and (r["n"] < row["n"] or r["S"] < row["S"])]
                require(bool(dominators), "A finite witness is not strictly dominated")
                metrics.add((identifier, distance, row["n"], row["S"]))
                witnesses.append(dict(block=index, q=q, witness=position, output_id=identifier,
                                      dominator=min(dominators), output_basis=list(basis)))
        if native is not None:
            require(work is not None and workers >= 1, "A native replay needs a work directory and positive worker count")
            replay_block(Path(native).resolve(), Path(work).resolve() / f"block_{index:03d}", block,
                         supports, workers, independent_labels)
            print(json.dumps(dict(replayed_block=index, source_spaces=end-first)), flush=True)
        total["pointing_supports"] += sum(counts)
        total["source_spaces"] += end-first
        intervals.append(dict(c=c, m=m, first_index=first, count=end-first))
    require(all(total[k] == data[k] for k in total), "Finite block totals disagree")
    return dict(schema="finite-pointed-source-closure-v1", status="pass", **total,
        maximum_protocol_length=54, minimum_distance=3, exact_distance_is_part_of_key=True,
        all_intrinsic_logical_dimensions_exhausted=True, all_five_pointing_cases_counted_independently=True,
        finite_censuses_freshly_recomputed=native is not None,
        source_schur_and_direct_label_spaces_compared=native is not None and independent_labels,
        all_emitted_witnesses_strictly_dominated=True, witness_checks=witnesses,
        frontier_metrics=[dict(output_id=o, d_Z=d, n=n, S=s) for o, d, n, s in sorted(metrics)],
        source_intervals=intervals, source_shard_bindings=bindings,
        data_sha256=hashlib.sha256(raw).hexdigest(), frontier_sha256=hashlib.sha256(frontier_path.read_bytes()).hexdigest(),
        is_global_completeness_certificate=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--native", type=Path)
    parser.add_argument("--work-directory", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--independent-labels", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.native, args.work_directory, args.workers, args.independent_labels)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in
                     ("witness_checks", "frontier_metrics", "source_intervals", "source_shard_bindings")}, indent=2))


if __name__ == "__main__":
    main()
