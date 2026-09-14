"""Reproduce either q5 target sector on one canonical length-52 stabiliser support."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from verify_primitive_support_sector import support_key, validate_zero_result
from verify_protocol_case_census import normalized_keys, write_pointed_cases


def replay(root, index, sector, native, work, workers):
    bindings = {}
    def load(name):
        raw = (root / name).read_bytes()
        bindings[name] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    data = load("data/protocol_sectors/length52_large_quotient_domain.json")
    profiles = load(data["profile_dataset"])
    if (data["schema"] != "finite-length52-large-quotient-domain-v1"
            or not 0 <= index < len(data["classes"]) or sector not in ("primitive", "nonprimitive") or workers < 1):
        raise ValueError("Invalid finite replay case")
    row = data["classes"][index]
    points, h = row["points"], row["ambient_dimension"]
    if support_key(points, h, row["key_words"]) != row["canonical_key_sha256"]:
        raise ValueError("Canonical support does not reconstruct its key")
    work.mkdir(parents=True, exist_ok=True)
    if sector == "primitive":
        if "primitive_target_census" in row:
            expected = row["primitive_target_census"]
            validate_zero_result(expected)
        else:
            union = load(data["primitive_domain"])
            if bindings[data["primitive_domain"]] != data["primitive_domain_sha256"]:
                raise ValueError("Changed primitive input domain")
            shared = union["cases"][row["primitive_union_case_index"]]
            if shared["canonical_key_sha256"] != row["canonical_key_sha256"]:
                raise ValueError("Primitive census refers to a different support")
            expected = shared["result"]
        result = subprocess.run([str(native), "canonical-isotropic-orbits", "--ambient", str(h),
            "--points", ",".join(map(str, points)), "--q", "5", "--target-signatures", "0x60000,0x60001",
            "--target-nondegenerate-seed", "3", "--workers", str(workers), "--emit", "100000", "--emit-orbit-data"],
            capture_output=True, text=True, check=True)
        value = json.loads(result.stdout)
        validate_zero_result(value)
        if (value["quotient_dimension"] != expected["quotient_dimension"]
                or value["support_automorphism_group_order"] != row["automorphism_group_order"]):
            raise ValueError("Fresh primitive-target geometry differs")
        checks = dict(quotient_dimension=value["quotient_dimension"], subspace_orbits=0, weighted_subspace_count=0)
    else:
        manifest, output = work / "case.utsp", work / "census.json"
        with manifest.open("wb") as stream:
            write_pointed_cases(stream, dict(maximum_protocol_length=52,
                cases=[dict(index=0, points=sorted(points), ambient_dimension=h)]))
        subprocess.run([str(native), "manifest-catalogue", "--input", str(manifest), "--output", str(output),
            "--q", "5", "--minimum-distance", "3", "--q5-q3-chain-cover-orbits", "--workers", str(workers),
            "--emit-positive-supports"], check=True)
        value = json.loads(output.read_bytes())
        expected = profiles["profiles"][index]
        if (value["status"] != "complete" or value["matrix_scope"] != "full_projective"
                or value["logical_qubits"] != 5 or value["minimum_distance"] != 3
                or value["enumeration_mode"] != "q5_q3_chain_cover_marked_code_orbits"
                or value["support_range"] != dict(start=0, count=1, end_exclusive=1)
                or value["statistics"]["supports_processed"] != 1
                or expected["canonical_geometry_sha256"] != row["canonical_key_sha256"]):
            raise ValueError("Incomplete or unbound finite chain census")
        keys = normalized_keys(key["canonical_key_words"] for item in value["positive_supports"] for key in item["output_keys"])
        if (keys != normalized_keys(expected["output_keys"])
                or value["statistics"]["nondegenerate_subspaces"] != expected["nondegenerate_subspaces"]):
            raise ValueError("Fresh complete profile or weighted subspace count differs")
        frontier = load("data/protocols/pareto_frontier.json")
        outputs = {(item["q"], tuple(int(word, 0) for word in item["tensor_key_words"])): item["output_id"]
                   for item in frontier["outputs"]}
        witnesses = []
        for item in value["pareto_protocols"]:
            masks, columns = check_matrix(item["generator_rows"], 5)
            tensor = LogicalTensor(masks[:5])
            distance, coefficient = distance_up_to(columns, 5)
            key, = normalized_keys([item["output"]["canonical_key_words"]])
            identifier = outputs.get((5, key))
            if (not tensor.intrinsic() or distance != 3 or identifier is None
                    or (item["protocol_length_n"], item["space_footprint_S"], item["d_Z"], item["error_coefficient"])
                    != (len(columns), len(masks), distance, coefficient)):
                raise ValueError("Invalid finite witness")
            choices = [p for p in frontier["protocols"] if p["output_id"] == identifier and p["d_Z"] == distance
                       and p["n"] <= len(columns) and p["S"] <= len(masks)]
            if not choices:
                raise ValueError("Fresh finite witness has an uncovered Pareto metric")
            chosen = min(choices, key=lambda p: (p["n"], p["S"], p["index"]))
            basis = find_equivalence_basis(tensor, LogicalTensor(tuple(int(s, 2) for s in chosen["generator_matrix_rows"][:5])))
            if basis is None:
                raise ValueError("Incorrect output-equivalence key")
            witnesses.append(dict(output_id=identifier, n=len(columns), S=len(masks), d_Z=distance,
                error_coefficient=coefficient, generator_matrix_rows=item["generator_rows"],
                dominator=chosen["index"], output_basis=list(basis)))
        checks = dict(output_keys=[list(k) for k in sorted(keys)],
            marked_code_orbits=value["statistics"]["marked_code_orbits"],
            weighted_nondegenerate_subspaces=value["statistics"]["nondegenerate_subspaces"], witnesses=witnesses)
    return dict(schema="length52-large-quotient-case-replay-v1", status="pass", index=index, sector=sector,
        canonical_key_sha256=row["canonical_key_sha256"], complete_finite_case_freshly_recomputed=True,
        dependencies=bindings, checks=checks, is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=int, required=True)
    parser.add_argument("--sector", choices=("primitive", "nonprimitive"), required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = replay(root, args.case, args.sector, args.native.resolve(), args.work_directory.resolve(), args.workers)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("dependencies", "checks")}, indent=2))


if __name__ == "__main__":
    main()
