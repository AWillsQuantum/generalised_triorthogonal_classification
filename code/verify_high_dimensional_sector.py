"""Check a complete finite block cover and its exact-distance Pareto exclusions."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from protocol_checks import LogicalTensor, check_matrix, distance_up_to
from verify_low_q_pruning import verify as verify_small_outputs


def check_intervals(intervals, count):
    cursor = 0
    for first, end in sorted(intervals):
        if type(first) is not int or type(end) is not int or first != cursor or end <= first:
            raise ValueError("Source intervals overlap or leave a gap")
        cursor = end
    if cursor != count:
        raise ValueError("Source intervals do not exhaust the catalogue sector")


def verify(root, data_path):
    root = Path(root).resolve()
    raw = data_path.read_bytes()
    data = json.loads(raw)
    if (data["schema"] != "finite-high-dimensional-protocol-sector-v1"
            or data["maximum_protocol_length"] != 54
            or data["pareto_application_length_interval"] != [53, 54]
            or data["source_lengths"] != [52, 54]
            or data["minimum_source_affine_dimension"] != 11
            or data["matrix_scope"] != "full_projective"
            or data["all_five_pointing_cases"] is not True
            or data["census_mode"] != "all_common_isotropic_subspaces"):
        raise ValueError("Incorrect finite-sector scope")
    catalogue = json.loads((root / "data/spaces/index.json").read_bytes())
    expected = {(length["c"], sector["m"]): sector["expected_classes"]
                for length in catalogue["lengths"] if length["c"] in (52, 54)
                for sector in length["sectors"] if sector["m"] >= 11 and sector["expected_classes"]}
    binding_intervals = defaultdict(list)
    for row in data["source_bindings"]:
        c, m, first, count = row["c"], row["m"], row["first_index"], row["count"]
        if (c, m) not in expected:
            raise ValueError("Unexpected source binding sector")
        encoding = json.loads((root / f"certificates/support_encoding/c{c}_m{m:02d}_{first:09d}.json").read_bytes())
        if (encoding["status"] != "pass" or encoding["first_index"] != first or encoding["count"] != count
                or row["ordered_support_sha256"] != encoding["ordered_support_sha256"]):
            raise ValueError("Unbound ordered source support catalogue")
        binding_intervals[c, m].append((first, first+count))
    if set(binding_intervals) != set(expected):
        raise ValueError("Incomplete source catalogue domain")
    for key, count in expected.items():
        check_intervals(binding_intervals[key], count)
    intervals = defaultdict(list)
    totals = defaultdict(Counter)
    witness_ownership = {}
    for index, block in enumerate(data["blocks"]):
        if block["index"] != index:
            raise ValueError("Noncontiguous finite block indices")
        key = block["c"], block["m"]
        if key not in expected:
            raise ValueError("Unexpected source sector")
        first, end = block["source_interval"]
        if end-first != block["source_count"]:
            raise ValueError("Incorrect block cardinality")
        intervals[key].append((first, end))
        pointing = block["pointing_counts"]
        if len(pointing) != 5 or any(type(n) is not int or n < 0 for n in pointing):
            raise ValueError("Invalid pointing counts")
        if sum(pointing) != block["reduced_supports"] or pointing[3] != block["source_count"]:
            raise ValueError("Incomplete five-case pointing cover")
        profiles = {int(distance): dict(profile) for distance, profile in block["quotient_profiles"].items()}
        if set(profiles) != {3, 4} or any(sum(profile.values()) != block["reduced_supports"]
                                         for profile in profiles.values()):
            raise ValueError("Incomplete filtered label-space profiles")
        if len(block["branches"]) != 8:
            raise ValueError("Missing logical dimension")
        for q, branch in enumerate(block["branches"], 1):
            distance = 4 if q <= 4 else 3
            if (branch["q"], branch["minimum_distance"]) != (q, distance):
                raise ValueError("Incorrect logical branch threshold")
            stats = branch["statistics"]
            eligible = sum(n for dimension, n in profiles[distance].items() if dimension >= q)
            if (stats["source_spaces"] != block["source_count"]
                    or stats["supports_processed"] != block["reduced_supports"]
                    or stats["eligible_supports"] != eligible
                    or stats["raw_enumerated_supports"] != eligible
                    or stats["marked_orbit_enumerated_supports"] != 0
                    or stats["quotient_dimension_filtered_supports"] != 0):
                raise ValueError("Raw enumeration fails to cover its filtered domain")
            if (any(not 1 <= int(r) <= q for r in stats["radical_dimension_counts"])
                    or sum(stats["radical_dimension_counts"].values()) + stats["nondegenerate_subspaces"]
                    != stats["isotropic_subspaces"]):
                raise ValueError("Inconsistent radical census")
            if q >= 2 and stats["nondegenerate_subspaces"]:
                raise ValueError("Unexpected positive higher-dimensional output")
            if q >= 7 and stats["isotropic_subspaces"]:
                raise ValueError("The common-isotropic dimension bound is not established")
            for witness in branch["witness_indices"]:
                if witness in witness_ownership:
                    raise ValueError("Repeated finite witness ownership")
                witness_ownership[witness] = index, q
            for name in ("source_spaces", "supports_processed", "eligible_supports",
                         "isotropic_subspaces", "nondegenerate_subspaces"):
                totals[q][name] += stats[name]
    if set(intervals) != set(expected):
        raise ValueError("Missing source sector")
    for key, count in expected.items():
        check_intervals(intervals[key], count)
    if set(witness_ownership) != set(range(len(data["witnesses"]))):
        raise ValueError("Missing finite-sector witnesses")
    frontier_path = root / "data/protocols/pareto_frontier.json"
    frontier = json.loads(frontier_path.read_bytes())
    small_output_certificate = verify_small_outputs(root)
    if not small_output_certificate["all_exact_distance_three_candidates_strictly_dominated"]:
        raise ValueError("Missing small-output distance-three pruning")
    dominators, checked_incumbents = [], {}
    for index, witness in enumerate(data["witnesses"]):
        if (witness["index"] != index or witness_ownership[index] != (witness["block"], witness["q"])
                or witness["q"] != 1):
            raise ValueError("Incorrect finite witness identity")
        rows = witness["generator_matrix_rows"]
        masks, columns = check_matrix(rows, 1)
        if LogicalTensor(masks[:1]).key_words() != ["0x0000000000000001"]:
            raise ValueError("The finite witness is not a T output")
        distance, coefficient = distance_up_to(columns, 1)
        if (witness["n"], witness["S"], witness["d_Z"], witness["error_coefficient"]) != (
                len(rows[0]), len(rows), distance, coefficient) or distance is None:
            raise ValueError("Incorrect finite witness metrics")
        matches = [r for r in frontier["protocols"] if r["q"] == 1 and r["d_Z"] == distance
                   and r["n"] <= witness["n"] and r["S"] <= witness["S"]]
        if not matches:
            raise ValueError("A finite witness is not covered by the frontier")
        incumbent = min(matches, key=lambda row: (row["n"], row["S"], row["index"]))
        if incumbent["index"] not in checked_incumbents:
            imasks, icolumns = check_matrix(incumbent["generator_matrix_rows"], 1)
            if (LogicalTensor(imasks[:1]).key_words() != ["0x0000000000000001"]
                    or distance_up_to(icolumns, 1) != (incumbent["d_Z"], incumbent["error_coefficient"])):
                raise ValueError("Invalid incumbent witness")
            checked_incumbents[incumbent["index"]] = True
        dominators.append(incumbent["index"])
    return dict(schema="finite-high-dimensional-sector-verification-v1", status="pass",
                data_sha256=hashlib.sha256(raw).hexdigest(),
                frontier_sha256=hashlib.sha256(frontier_path.read_bytes()).hexdigest(),
                source_spaces=sum(expected.values()), blocks=len(data["blocks"]),
                witnesses=len(data["witnesses"]), source_interval_cover_exact=True,
                recorded_census_domains_checked=True, all_finite_witnesses_independently_verified=True,
                exact_distance_dominator_indices=dominators,
                statistics={str(q): dict(values) for q, values in totals.items()},
                logical_dimensions_at_least_seven_excluded_by_common_isotropic_containment=True,
                all_block_enumerations_freshly_recomputed=False,
                is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/length54_high_dimensional.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(root, args.input)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("statistics", "exact_distance_dominator_indices")}, indent=2))


if __name__ == "__main__":
    main()
