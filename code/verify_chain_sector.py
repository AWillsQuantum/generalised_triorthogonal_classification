"""Check the finite chain-sector partition, aggregate results and witnesses."""

import argparse
import hashlib
import json
from pathlib import Path

from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from space_lifts import linear_solutions
from verify_profile_exclusions import necessary_profiles
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_frontier


def verify(root, path):
    bindings = {}
    def load(path):
        raw = path.read_bytes()
        if not path.resolve().is_relative_to(root):
            raise ValueError("External chain-sector dependency")
        bindings[path.relative_to(root).as_posix()] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    data = load(path)
    if (data["schema"] != "finite-nonprimitive-chain-sector-v1"
            or data["maximum_protocol_length"] != 54 or data["minimum_distance"] != 3
            or data["evidence_granularity"] != "aggregate_censuses_over_explicit_finite_input_domains"):
        raise ValueError("Wrong finite chain-sector scope")
    domain = load(root / data["domain_dataset"])
    finite = load(root / data["finite_subsector_dataset"])
    primitive = load(root / "data/protocol_sectors/primitive_support_union.json")
    primitive_certificate = load(root / "certificates/primitive_support_sector.json")
    if (primitive_certificate["status"] != "pass" or not primitive_certificate["normalisation_equivalences_freshly_recomputed"]
            or primitive_certificate["input_sha256"] != bindings["data/protocol_sectors/primitive_support_union.json"]):
        raise ValueError("Unverified primitive input cover")
    sets = [data["q5_input_support_class_indices"], data["zero_column_dominated_support_class_indices"],
            [row["support_class_index"] for row in finite["cases"]]]
    if (any(len(values) != len(set(values)) for values in sets)
            or sum(map(len, sets)) != len(domain["classes"])
            or set().union(*map(set, sets)) != set(range(len(domain["classes"])))):
        raise ValueError("The larger-quotient domain is not partitioned exactly once")
    for index in sets[1]:
        cls = domain["classes"][index]
        normalisation = primitive["normalisation_inputs"][cls["normalisation_input_index"]]
        if (not normalisation["deleted_zero_column"] or cls["support_length"] != 53
                or len(normalisation["points"]) != 52):
            raise ValueError("Zero-column exclusion is not justified")
    census = data["q5_census"]
    stats5, stats6 = data["statistics"]["5"], data["statistics"]["6"]
    if (census["complete_input_classes"] != len(sets[0])
            or census["positive_input_classes"]+census["zero_input_classes"] != len(sets[0])
            or census["positive_input_classes"] != stats5["supports_processed"]
            or census["output_candidates"] != stats5["canonical_output_orbits"]):
        raise ValueError("Incomplete aggregate chain census")
    profiles = necessary_profiles(load(root / "certificates/output_profiles/q6_necessary_profiles.json"))
    candidates = set()
    matched_profiles = []
    for index, row in enumerate(data["q6_candidates"]):
        available = normalized_keys(row["predecessor_output_keys"])
        matched = [i for i, profile in enumerate(profiles) if profile.issubset(available)]
        cls = row["support_class_index"]
        if (row["index"] != index or cls not in sets[0] or cls in candidates or not matched
                or row["predecessor_weighted_subspaces"] < 1):
            raise ValueError("Invalid recorded successor obligation")
        candidates.add(cls)
        matched_profiles.append(dict(index=index, necessary_profiles=matched))
    if (data["q6_complete_input_classes"] != len(candidates)
            or stats6["supports_processed"] != len(candidates)
            or stats6["canonical_output_orbits"] != len(candidates)):
        raise ValueError("Successor aggregate does not cover its input set")
    gate_certificate = load(root / "certificates/output_profiles/certificate.json")
    gate = load(root / "certificates/output_profiles/q7_anchored_exclusions.json")
    if (gate_certificate["status"] != "pass"
            or gate_certificate["files"]["q7_anchored_exclusions.json"]
            != bindings["certificates/output_profiles/q7_anchored_exclusions.json"]):
        raise ValueError("Unbound finite higher-dimensional profile gate")
    q6_keys = normalized_keys(data["q6_output_key_union"])
    excluded = [normalized_keys([[key] for key in row["q6_keys"]]) for row in gate["support_profiles"]
                if row["surviving_anchored_extensions"] == 0]
    if not any(q6_keys.issubset(profile) for profile in excluded):
        raise ValueError("Nonprimitive q7 outputs not excluded by the aggregate profile")
    frontier_path = root / "data/protocols/pareto_frontier.json"
    verify_frontier(frontier_path)
    frontier = load(frontier_path)
    outputs = {(row["q"], tuple(int(word, 0) for word in row["tensor_key_words"])): row for row in frontier["outputs"]}
    tensors = {row["output_id"]: LogicalTensor(tuple(int(s, 2) for s in row["generator_matrix_rows"][:row["q"]]))
               for row in frontier["protocols"]}
    witness_checks, q5_keys, q6_witness_keys = [], set(), set()
    for index, witness in enumerate(data["witnesses"]):
        q = witness["q"]
        if witness["index"] != index or q not in (5, 6):
            raise ValueError("Incorrect aggregate witness identity")
        rows = witness["generator_matrix_rows"]
        masks, columns = check_matrix(rows, q)
        tensor = LogicalTensor(masks[:q])
        key, = normalized_keys([witness["output_key_words"]])
        output = outputs.get((q, key))
        if output is None or not tensor.intrinsic():
            raise ValueError("Unknown or degenerate aggregate output")
        basis = find_equivalence_basis(tensor, tensors[output["output_id"]])
        distance, coefficient = distance_up_to(columns, q)
        if (basis is None or distance is None or distance < 3 or
                (witness["n"], witness["S"], witness["d_Z"], witness["error_coefficient"])
                != (len(rows[0]), len(rows), distance, coefficient)):
            raise ValueError("Incorrect aggregate protocol witness")
        choices = [row for row in frontier["protocols"] if row["output_id"] == output["output_id"]
                   and row["d_Z"] == distance and row["n"] <= witness["n"] and row["S"] <= witness["S"]]
        if not choices:
            raise ValueError("Aggregate witness is not covered by the frontier")
        incumbent = min(choices, key=lambda row: (row["n"], row["S"], row["index"]))
        witness_checks.append(dict(index=index, dominator=incumbent["index"], output_basis=list(basis)))
        (q5_keys if q == 5 else q6_witness_keys).add(key)
    if (len(q5_keys) != census["frontier_points"] or q6_witness_keys != q6_keys):
        raise ValueError("Incomplete aggregate output witnesses")
    tensor_census = load(root / "certificates/tensors/q5_orbits.json")
    alternating = [row for row in tensor_census["orbits"] if row["radical_dimension"] == 0
                   and int(row["canonical_standard_signature"], 0) & ((1 << 15)-1) == 0]
    if (tensor_census["status"] != "pass" or tensor_census["tensor_count"] != 1 << 25
            or len(alternating) != 1 or normalized_keys([alternating[0]["canonical_key_words"]]) != {(0x600000,)}):
        raise ValueError("The alternating q5 tensor exclusion is not established")
    constant_rows = []
    for index in sets[0]:
        cls = domain["classes"][index]
        case = primitive["cases"][cls["primitive_case_index"]]
        solution = linear_solutions(case["points"], [1]*len(case["points"]), case["ambient_dimension"])
        if solution is not None:
            constant_rows.append(dict(support_class_index=index, primitive_case_index=case["index"],
                                      constant_row_coefficients=solution[0]))
    return dict(schema="finite-chain-sector-consistency-v1", status="pass", dependencies=bindings,
        q5_input_classes=len(sets[0]), zero_column_dominated_classes=len(sets[1]), separate_finite_classes=len(sets[2]),
        q6_input_classes=len(candidates), aggregate_witnesses=len(witness_checks),
        exact_input_partition_checked=True, aggregate_count_consistency_checked=True,
        every_recorded_successor_satisfies_necessary_profile=True,
        recorded_successor_exclusions_on_other_inputs_freshly_recomputed=False,
        complete_aggregate_enumerations_freshly_recomputed=False,
        higher_dimension_profile_implication_checked=True,
        constant_stabiliser_q5_exclusions=constant_rows,
        witness_checks=witness_checks, candidate_profile_checks=matched_profiles,
        premise="The supplied aggregate censuses are complete on their explicit finite input domains",
        is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/nonprimitive_chain_sector.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(root, args.input)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("dependencies", "witness_checks", "candidate_profile_checks")}, indent=2))


if __name__ == "__main__":
    main()
