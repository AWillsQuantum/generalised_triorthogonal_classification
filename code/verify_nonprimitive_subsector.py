"""Verify and reproduce a complete finite q5 marked-orbit subsector."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from verify_primitive_support_sector import support_key
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_frontier


def orbit_set(rows):
    result = {tuple(int(word, 16) for word in row["marked_key_words"]): row["orbit_size"] for row in rows}
    if len(result) != len(rows) or any(type(mass) is not int or mass <= 0 for mass in result.values()):
        raise ValueError("Invalid or repeated marked orbit")
    return result


def verify(root, data):
    if (data["schema"] != "finite-nonprimitive-marked-sector-v1" or data["logical_qubits"] != 5
            or data["minimum_distance"] != 3 or data["maximum_protocol_length"] != 54
            or data["enumeration_family"] != "q3_chain_cover"):
        raise ValueError("Wrong finite enumeration family")
    paths = [root / data[name] for name in ("domain_dataset", "profile_dataset")]
    if any(not path.resolve().is_relative_to(root) for path in paths):
        raise ValueError("External finite-sector data")
    raw = [path.read_bytes() for path in paths]
    domain, profiles = map(json.loads, raw)
    frontier_path = root / "data/protocols/pareto_frontier.json"
    verify_frontier(frontier_path)
    frontier = json.loads(frontier_path.read_bytes())
    outputs = {(row["q"], tuple(int(word, 0) for word in row["tensor_key_words"])): row for row in frontier["outputs"]}
    tensors = {row["output_id"]: LogicalTensor(tuple(int(s, 2) for s in row["generator_matrix_rows"][:row["q"]]))
               for row in frontier["protocols"]}
    checks = []
    for index, (case, profile) in enumerate(zip(data["cases"], profiles["profiles"], strict=True)):
        cls = domain["classes"][case["support_class_index"]]
        if (case["index"] != index or profile["index"] != index
                or support_key(case["points"], case["ambient_dimension"], cls["key_words"]) != cls["canonical_key_sha256"]
                or profile["canonical_geometry_sha256"] != cls["canonical_key_sha256"]
                or case["support_automorphism_group_order"] != cls["automorphism_group_order"]):
            raise ValueError("Finite census is not bound to its source domain")
        orbits = orbit_set(case["orbits"])
        if (len(orbits) != case["marked_orbits"] or sum(orbits.values()) != case["weighted_subspaces"]
                or normalized_keys(case["output_keys"]) != normalized_keys(profile["output_keys"])
                or profile["nondegenerate_subspaces"] != case["weighted_subspaces"]
                or profile["quotient_dimension"] != case["quotient_dimension"]):
            raise ValueError("Inconsistent finite orbit census or complete output profile")
        seen_outputs = set()
        for position, witness in enumerate(case["orbits"]):
            if witness["index"] != position or case["support_automorphism_group_order"] % witness["orbit_size"]:
                raise ValueError("Invalid orbit index or orbit-stabiliser mass")
            rows = witness["generator_matrix_rows"]
            masks, columns = check_matrix(rows, 5)
            tensor = LogicalTensor(masks[:5])
            key, = normalized_keys([witness["output_key_words"]])
            output = outputs.get((5, key))
            if output is None or not tensor.intrinsic():
                raise ValueError("Unknown or degenerate finite witness output")
            basis = find_equivalence_basis(tensor, tensors[output["output_id"]])
            if basis is None:
                raise ValueError("Finite witness output key is incorrect")
            distance, coefficient = distance_up_to(columns, 5)
            if (witness["q"], witness["n"], witness["S"], witness["d_Z"], witness["error_coefficient"]) != (
                    5, len(rows[0]), len(rows), distance, coefficient) or distance is None or distance < 3:
                raise ValueError("Incorrect finite protocol parameters")
            candidates = [row for row in frontier["protocols"] if row["output_id"] == output["output_id"]
                          and row["d_Z"] == distance and row["n"] <= witness["n"] and row["S"] <= witness["S"]]
            if not candidates:
                raise ValueError("Finite witness not covered by the frontier")
            incumbent = min(candidates, key=lambda row: (row["n"], row["S"], row["index"]))
            seen_outputs.add(key)
            checks.append(dict(case_index=index, orbit_index=position, dominator=incumbent["index"],
                               output_basis_in_witness_coordinates=list(basis)))
        if seen_outputs != normalized_keys(case["output_keys"]):
            raise ValueError("Output profile lacks a marked-orbit witness")
    return dict(schema="finite-nonprimitive-marked-sector-verification-v1", status="pass",
        cases=len(data["cases"]), orbits=len(checks), all_witness_metrics_independently_verified=True,
        all_finite_candidates_covered_by_frontier=True, marked_enumerations_freshly_recomputed=False,
        is_global_completeness_certificate=False, witness_checks=checks,
        supporting_data_sha256={path.relative_to(root).as_posix(): hashlib.sha256(content).hexdigest()
                               for path, content in zip(paths, raw, strict=True)},
        frontier_sha256=hashlib.sha256(frontier_path.read_bytes()).hexdigest())


def replay(data, native, work, workers):
    checked = []
    for case in data["cases"]:
        result = subprocess.run([str(native), "canonical-isotropic-orbits", "--ambient",
            str(case["ambient_dimension"]), "--points", ",".join(map(str, case["points"])), "--q", "5",
            "--q5-q3-chain-cover", "--workers", str(workers), "--emit", "1000000", "--emit-orbit-data"],
            capture_output=True, text=True, check=True)
        result = json.loads(result.stdout)
        actual = {tuple(int(word, 16) for word in key): mass for key, mass in
                  zip(result["canonical_keys"], result["orbit_sizes"], strict=True)}
        if (not result["q5_q3_chain_cover"] or result["logical_dimension"] != 5
                or result["quotient_dimension"] != case["quotient_dimension"]
                or result["support_automorphism_group_order"] != case["support_automorphism_group_order"]
                or result["subspace_orbits"] != len(actual)
                or result["weighted_subspace_count"] != case["weighted_subspaces"]
                or len(result["canonical_keys"]) != len(actual) or actual != orbit_set(case["orbits"])):
            raise ValueError("Fresh complete marked-orbit census differs")
        checked.append(case["index"])
        (work / f"case_{case['index']:03d}.json").write_text(json.dumps(dict(index=case["index"],
            status="pass", exact_marked_keys_and_orbit_masses_match=True,
            orbits=len(actual), weighted_subspaces=result["weighted_subspace_count"]), indent=2)+"\n", encoding="ascii")
        print(json.dumps(dict(replayed_cases=len(checked), total_cases=len(data["cases"]))), flush=True)
    return dict(marked_enumerations_freshly_recomputed=True, replayed_case_indices=checked)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/nonprimitive_marked_subsector.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--work-directory", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1 or (args.native is None) != (args.work_directory is None):
        parser.error("A native replay needs an executable, work directory and positive workers")
    raw = args.input.read_bytes()
    data = json.loads(raw)
    result = verify(root, data)
    if args.native:
        args.work_directory.mkdir(parents=True, exist_ok=True)
        result.update(replay(data, args.native, args.work_directory, args.workers))
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "witness_checks"}, indent=2))


if __name__ == "__main__":
    main()
