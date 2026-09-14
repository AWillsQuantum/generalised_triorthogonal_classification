"""Check the complete lower-dimensional source intervals and finite branch evidence."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from protocol_checks import (LogicalTensor, check_matrix, distance_up_to,
                             find_equivalence_basis)
from verify_high_dimensional_sector import check_intervals
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_frontier


def histogram(values, expected_count):
    result = dict(values)
    if (len(values) != len(result) or any(type(d) is not int or d < 0 or type(n) is not int or n < 0
                                        for d, n in result.items()) or sum(result.values()) != expected_count):
        raise ValueError("Incomplete filtered quotient histogram")
    return result


def check_raw_pointing_count(c, m, source_count, points):
    if c == 54:
        if points != source_count*((1 << m)+c+1):
            raise ValueError("Incomplete raw length-54 pointing domain")
        return 0
    if c != 52:
        raise ValueError("Unexpected parent length")
    # A length-52 unital support can have a one-dimensional translation
    # period. Exact duplicate point sets are removed even in raw-origin mode.
    half = ((1 << m)-c)//2
    difference = source_count*(2*half+1)-points
    if difference < 0 or difference % half or difference//half > source_count:
        raise ValueError("Raw pointing count is inconsistent with translation periods")
    return difference//half


def verify(root, path):
    root = Path(root).resolve()
    raw = path.read_bytes()
    data = json.loads(raw)
    if (data["schema"] != "finite-low-dimensional-protocol-blocks-v1"
            or data["maximum_protocol_length"] != 54
            or data["pareto_application_length_interval"] != [53, 54]
            or data["matrix_scope"] != "full_projective"
            or data["q5_small_quotient_interval"] != [5, 14]
            or data["q5_large_quotient_minimum"] != 15):
        raise ValueError("Incorrect finite block scope")
    profile_path = root / data["profile_dataset"]
    profiles = json.loads(profile_path.read_bytes())["profiles"]
    catalogue = json.loads((root / "data/spaces/index.json").read_bytes())
    expected = {(length["c"], sector["m"]): sector["expected_classes"]
                for length in catalogue["lengths"] if length["c"] in (52, 54)
                for sector in length["sectors"] if sector["m"] <= 10 and sector["expected_classes"]}
    frontier_path = root / "data/protocols/pareto_frontier.json"
    verify_frontier(frontier_path)
    frontier = json.loads(frontier_path.read_bytes())
    outputs = {(row["q"], tuple(int(word, 0) for word in row["tensor_key_words"])): row for row in frontier["outputs"]}
    tensors = {}
    for row in frontier["protocols"]:
        if row["output_id"] not in tensors:
            tensors[row["output_id"]] = LogicalTensor(tuple(int(s, 2) for s in row["generator_matrix_rows"][:row["q"]]))
    intervals, totals = defaultdict(list), Counter()
    witnesses = []
    if len(data["blocks"]) != len(profiles):
        raise ValueError("Q5 profiles do not match the source block domain")
    for index, (block, profile) in enumerate(zip(data["blocks"], profiles, strict=True)):
        if block["index"] != index or block["q5_profile_index"] != index or profile["index"] != index:
            raise ValueError("Conflicting finite block and profile identities")
        c, m = key = block["c"], block["m"]
        if key not in expected:
            raise ValueError("Unexpected source sector")
        first, end = block["source_interval"]
        source_count = end-first
        intervals[key].append((first, end))
        if block["protocol_length_interval"] != ([53, 53] if c == 52 else [53, 54]):
            raise ValueError("Wrong exact protocol length filter")
        if block["pointing_mode"] == "all_origins":
            totals["translation_period_sources"] += check_raw_pointing_count(c, m, source_count, block["pointing_supports"])
            totals["raw_pointing_blocks"] += 1
        elif block["pointing_mode"] == "affine_automorphism_orbits":
            totals["origin_orbit_blocks"] += 1
        else:
            raise ValueError("Unknown pointing quotient")
        d3 = histogram(block["distance_three_quotient_profile"], block["pointing_supports"])
        d4 = histogram(block["distance_four_quotient_profile"], block["pointing_supports"])
        large = sum(n for d, n in d3.items() if d >= 15)
        if large != block["q5_large_quotient_supports"]:
            raise ValueError("Incomplete larger-quotient successor domain")
        totals["larger_quotient_supports"] += large
        totals["pointed_supports"] += block["pointing_supports"]
        if len(block["branches"]) != 5:
            raise ValueError("A logical branch is missing")
        for q, branch in enumerate(block["branches"], 1):
            floor = 3 if q == 5 else 4
            if branch["q"] != q or branch["minimum_distance"] != floor:
                raise ValueError("Incorrect branch threshold")
            stats = branch["statistics"]
            distribution = d3 if q == 5 else d4
            input_count = sum(n for d, n in distribution.items() if d >= (5 if q == 5 else 1))
            eligible = sum(n for d, n in distribution.items() if (5 <= d <= 14 if q == 5 else d >= q))
            if (branch["support_range"] != dict(start=0, count=input_count, end_exclusive=input_count)
                    or stats["supports_processed"] != input_count or stats["source_spaces"] != source_count
                    or stats["eligible_supports"] != eligible
                    or stats["raw_enumerated_supports"]+stats["marked_orbit_enumerated_supports"] != eligible
                    or branch["quotient_dimension_interval"] != ([5, 14] if q == 5 else [0, None])):
                raise ValueError("Incomplete logical-subspace domain")
            if q == 5:
                if (stats["marked_orbit_enumerated_supports"] != 0
                        or stats["quotient_dimension_filtered_supports"] != large):
                    raise ValueError("The small-quotient census is not wholly raw")
                for name, value in profile["counts"].items():
                    if stats[name] != value:
                        raise ValueError("The complete output profile is bound to a different census")
            elif stats["quotient_dimension_filtered_supports"]:
                raise ValueError("A higher-distance branch has an unaccounted quotient exclusion")
            if stats["isotropic_subspace_statistics_exact"]:
                if (not stats["radical_statistics_exact"] or stats["marked_orbit_enumerated_supports"]
                        or sum(stats["radical_dimension_counts"].values())+stats["nondegenerate_subspaces"]
                        != stats["isotropic_subspaces"]):
                    raise ValueError("Inconsistent exact raw subspace counts")
            elif stats["isotropic_subspaces"] is not None or not stats["marked_orbit_enumerated_supports"]:
                raise ValueError("An inexact orbit count is incorrectly presented as a raw count")
            totals[f"q{q}_raw_supports"] += stats["raw_enumerated_supports"]
            totals[f"q{q}_marked_supports"] += stats["marked_orbit_enumerated_supports"]
            for position, witness in enumerate(branch["witnesses"]):
                rows = witness["generator_matrix_rows"]
                masks, columns = check_matrix(rows, q)
                tensor = LogicalTensor(masks[:q])
                output_key, = normalized_keys([witness["output_key_words"]])
                output = outputs.get((q, output_key))
                if output is None or output["q"] != q or not tensor.intrinsic():
                    raise ValueError("Unknown or degenerate finite witness output")
                basis = find_equivalence_basis(tensor, tensors[output["output_id"]])
                if basis is None:
                    raise ValueError("Incorrect canonical output key")
                distance, coefficient = distance_up_to(columns, q)
                if (witness["q"], witness["n"], witness["S"], witness["d_Z"], witness["error_coefficient"]) != (
                        q, len(rows[0]), len(rows), distance, coefficient) or distance is None or distance < floor:
                    raise ValueError("Incorrect exact finite witness parameters")
                dominators = [row for row in frontier["protocols"] if row["output_id"] == output["output_id"]
                              and row["d_Z"] == distance and row["n"] <= witness["n"] and row["S"] <= witness["S"]]
                if not dominators:
                    raise ValueError("A recorded finite witness is not covered by the frontier")
                incumbent = min(dominators, key=lambda row: (row["n"], row["S"], row["index"]))
                witnesses.append(dict(block=index, q=q, position=position, dominator=incumbent["index"],
                                      output_basis_in_witness_coordinates=list(basis)))
    if set(intervals) != set(expected):
        raise ValueError("A source sector is absent")
    for key, count in expected.items():
        # All length-52 parents yield a zero column at protocol length 53,
        # so their entire sector is Pareto-redundant by the zero-column lemma.
        # The included auxiliary census covers the initial 186 cubic entries.
        checked_count = 186 if key == (52, 7) else count
        check_intervals(intervals[key], checked_count)
    return dict(schema="finite-low-dimensional-block-verification-v1", status="pass",
                data_sha256=hashlib.sha256(raw).hexdigest(),
                output_profile_sha256=hashlib.sha256(profile_path.read_bytes()).hexdigest(),
                frontier_sha256=hashlib.sha256(frontier_path.read_bytes()).hexdigest(),
                blocks=len(data["blocks"]), witnesses=len(witnesses),
                all_length54_source_intervals_through_dimension10_covered=True,
                length52_parent_extensions_excluded_by_zero_column_domination=True,
                larger_quotient_enumeration_is_separate_obligation=True,
                all_enumerations_freshly_recomputed=False,
                totals=dict(totals), witness_checks=witnesses,
                is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/length54_low_dimensional_blocks.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(root, args.input)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "witness_checks"}, indent=2))


if __name__ == "__main__":
    main()
