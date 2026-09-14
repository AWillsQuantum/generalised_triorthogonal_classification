"""Check a complete finite source-block sector and its cumulative Pareto frontier."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from verify_high_dimensional_sector import check_intervals
from verify_low_q_pruning import verify as verify_low_q
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_frontier
from verify_through48_sector import nondominated_metrics


EXACT_COUNTS = ("eligible_supports", "raw_enumerated_supports", "marked_orbit_enumerated_supports",
                "isotropic_subspaces", "nondegenerate_subspaces", "marked_code_orbits", "canonical_output_orbits")


def verify(root, path, result_directory, allow_incomplete=False):
    raw = path.read_bytes()
    data = json.loads(raw)
    source_hash = hashlib.sha256(raw).hexdigest()
    if (data["schema"] != "finite-protocol-source-block-domain-v1" or data["matrix_scope"] != "full_projective"
            or data["minimum_protocol_length"] < 49 or data["maximum_protocol_length"] > 54):
        raise ValueError("Source sector does not satisfy the stated small-output Pareto bound")
    catalogue = json.loads((root / "data/spaces/index.json").read_bytes())
    length, = [row for row in catalogue["lengths"] if row["c"] == data["source_length"]]
    expected = {row["m"]: row["expected_classes"] for row in length["sectors"] if row["expected_classes"]}
    intervals = defaultdict(list)
    for index, block in enumerate(data["blocks"]):
        if block["index"] != index or block["c"] != data["source_length"]:
            raise ValueError("Incorrect source block identity")
        intervals[block["m"]].append(tuple(block["source_interval"]))
    if set(intervals) != set(expected) or sum(expected.values()) != data["source_spaces"]:
        raise ValueError("Source block dimension partition is incomplete")
    for m, count in expected.items():
        check_intervals(intervals[m], count)
    if sum(block["expected_pointed_supports"] for block in data["blocks"]) != data["expected_pointed_supports"]:
        raise ValueError("Wrong complete origin-domain cardinality")
    frontier_path = root / "data/protocols/pareto_frontier.json"
    verify_frontier(frontier_path)
    low_q = verify_low_q(root)
    frontier = json.loads(frontier_path.read_bytes())
    outputs = {(row["q"], tuple(int(word, 0) for word in row["tensor_key_words"])): row["output_id"]
               for row in frontier["outputs"]}
    tensors = {row["output_id"]: LogicalTensor(tuple(int(s, 2) for s in row["generator_matrix_rows"][:row["q"]]))
               for row in frontier["protocols"]}
    branch_specs = {row["q"]: row for row in data["aggregate_censuses"]}
    if set(branch_specs) != set(range(1, 9)) or any(
            row["minimum_distance"] != (4 if q <= 4 else 3) for q, row in branch_specs.items()):
        raise ValueError("Incomplete or incorrect logical branch scope")
    missing, records, witnesses, results = [], {}, [], {}
    aggregates, radicals = defaultdict(Counter), defaultdict(Counter)
    candidates = set()
    for index, block in enumerate(data["blocks"]):
        file = result_directory / f"block{index:03d}.json"
        if not file.exists():
            missing.append(index)
            continue
        content = file.read_bytes()
        result = json.loads(content)
        if (result["schema"] != "finite-protocol-source-block-census-v1" or result["status"] != "pass"
                or result["source_domain_sha256"] != source_hash or result["index"] != index
                or (result["c"], result["m"], result["source_interval"]) != (block["c"], block["m"], block["source_interval"])
                or result["minimum_protocol_length"] != data["minimum_protocol_length"]
                or result["maximum_protocol_length"] != data["maximum_protocol_length"]
                or result["pointed_supports"] != block["expected_pointed_supports"]
                or not result["all_censuses_freshly_recomputed"]):
            raise ValueError("Finite result is not bound to its complete source block")
        records[file.relative_to(root).as_posix()] = hashlib.sha256(content).hexdigest()
        results[index] = result
        profiles = {3: result["distance_three_quotient_profile"], 4: result["distance_four_quotient_profile"]}
        for pairs in profiles.values():
            if (len(pairs) != len(dict(pairs)) or any(type(d) is not int or type(n) is not int or d < 0 or n < 0 for d, n in pairs)
                    or sum(n for _, n in pairs) != result["pointed_supports"]):
                raise ValueError("Incomplete quotient histogram")
        if len(result["branches"]) != 8:
            raise ValueError("Missing logical branch")
        for q, branch in enumerate(result["branches"], 1):
            distance = 4 if q <= 4 else 3
            minimum = 1 if q <= 4 else 5
            stats = branch["statistics"]
            if (branch["q"] != q or branch["minimum_distance"] != distance
                    or branch["minimum_input_quotient_dimension"] != minimum
                    or stats["source_spaces"] != block["source_interval"][1]-block["source_interval"][0]
                    or stats["supports_processed"] != sum(n for d, n in profiles[distance] if d >= minimum)
                    or stats["eligible_supports"] != sum(n for d, n in profiles[distance] if d >= q)
                    or stats["raw_enumerated_supports"] != stats["eligible_supports"]
                    or stats["quotient_dimension_filtered_supports"] or stats["marked_orbit_enumerated_supports"]
                    or not stats["isotropic_subspace_statistics_exact"] or not stats["radical_statistics_exact"]
                    or sum(stats["radical_dimension_counts"].values())+stats["nondegenerate_subspaces"] != stats["isotropic_subspaces"]):
                raise ValueError("Incomplete direct logical census or omitted eligible quotient")
            aggregates[q].update({k: stats[k] for k in EXACT_COUNTS})
            radicals[q].update(stats["radical_dimension_counts"])
            for position, witness in enumerate(branch["witnesses"]):
                masks, columns = check_matrix(witness["generator_matrix_rows"], q)
                tensor = LogicalTensor(masks[:q])
                key, = normalized_keys([witness["output_key_words"]])
                identifier = outputs.get((q, key))
                if identifier is None or not tensor.intrinsic():
                    raise ValueError("Unrecognised or degenerate finite witness output")
                basis = find_equivalence_basis(tensor, tensors[identifier])
                actual = distance_up_to(columns, q, maximum=witness["d_Z"])
                if (basis is None or actual != (witness["d_Z"], witness["error_coefficient"])
                        or witness["d_Z"] < distance or witness["n"] != len(columns) or witness["S"] != len(masks)
                        or not data["minimum_protocol_length"] <= witness["n"] <= data["maximum_protocol_length"]):
                    raise ValueError("Incorrect finite witness parameters or output equivalence")
                candidates.add((identifier, witness["d_Z"], witness["n"], witness["S"]))
                dominators = [row for row in frontier["protocols"] if row["output_id"] == identifier
                              and row["d_Z"] == witness["d_Z"] and row["n"] <= witness["n"] and row["S"] <= witness["S"]]
                if not dominators:
                    raise ValueError("Finite witness is not covered by the final frontier")
                witnesses.append(dict(block=index, q=q, witness=position, output_basis=list(basis),
                                      dominator=min(row["index"] for row in dominators)))
    if missing and not allow_incomplete:
        raise ValueError(f"Missing {len(missing)} finite source blocks")
    complete = not missing
    for q, branch in branch_specs.items():
        for name in EXACT_COUNTS:
            count, expected_count = aggregates[q][name], branch["statistics"][name]
            if count > expected_count or complete and count != expected_count:
                raise ValueError(f"Fresh aggregate q{q} differs on {name}")
        for dimension in set(radicals[q]) | set(branch["statistics"]["radical_dimension_counts"]):
            count = radicals[q][dimension]
            expected_count = branch["statistics"]["radical_dimension_counts"].get(dimension, 0)
            if count > expected_count or complete and count != expected_count:
                raise ValueError("Fresh aggregate radical profile differs")
    if aggregates[8]["isotropic_subspaces"]:
        raise ValueError("The zero-isotropic higher-logical-dimension bound does not hold")
    previous = {(row["output_id"], row["d_Z"], row["n"], row["S"]) for row in frontier["protocols"]
                if row["n"] < data["minimum_protocol_length"]}
    merged = nondominated_metrics(previous | candidates)
    released = {(row["output_id"], row["d_Z"], row["n"], row["S"]) for row in frontier["protocols"]
                if row["n"] <= data["maximum_protocol_length"]}
    if complete and merged != released:
        raise ValueError("Cumulative finite metrics do not equal the released frontier")
    return dict(schema="finite-source-block-sector-verification-v1", status="pass" if complete else "partial",
        source_domain_sha256=source_hash, maximum_protocol_length=data["maximum_protocol_length"],
        source_spaces=data["source_spaces"], complete_source_blocks=len(results), missing_source_blocks=missing,
        source_intervals_exactly_cover_compact_catalogue=True,
        pointed_supports_checked=sum(row["pointed_supports"] for row in results.values()),
        witnesses_independently_checked=len(witnesses), all_finite_censuses_freshly_recomputed=complete,
        exact_aggregate_counts_agree=complete, all_q_at_least8_excluded_by_isotropic_subspace_containment=complete,
        exact_d3_small_output_exclusions_verified=True,
        lower_parent_zero_column_extensions_pareto_redundant=True,
        cumulative_frontier_points=len(merged), cumulative_frontier_equals_release=complete,
        frontier_sha256=hashlib.sha256(frontier_path.read_bytes()).hexdigest(),
        low_q_census_sha256=low_q["q5_census_sha256"], result_bindings=records,
        witnesses=witnesses, is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/length50_source_domain.json")
    parser.add_argument("--results", type=Path, default=root / "data/protocol_sectors/length50_blocks")
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(root, args.input, args.results, args.allow_incomplete)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("result_bindings", "witnesses", "missing_source_blocks")}, indent=2))


if __name__ == "__main__":
    main()
