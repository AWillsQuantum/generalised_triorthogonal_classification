"""Verify all logical dimensions for a finite set of cubic parent spaces."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess

from logical_spaces import logical_label_space
from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from selected_spaces import read_selected
from space_codec import binary_rank, checked_support
from space_lifts import evaluation_rows
from verify_finite_sector_closure import verify as verify_q6_closure
from verify_low_q_pruning import verify as verify_low_q
from verify_primitive_restrictions import verify as primitive_restrictions
from verify_profile_exclusions import necessary_profiles
from verify_protocol_case_census import normalized_keys, write_pointed_cases
from verify_protocols import verify as verify_frontier
from verify_small_source_closure import all_pointings
from verify_tensor_authorities import standard_tensor


def check_raw_census(row, dimension):
    q, floor, stats = row["q"], row["minimum_distance"], row["statistics"]
    if (q not in (3, 5) or floor not in (3, 4) or row["enumeration_mode"] != "raw_isotropic_subspaces"
            or stats["source_spaces"] != 1 or stats["supports_processed"] != 1
            or stats["quotient_dimension_filtered_supports"] or stats["eligible_supports"] != int(dimension >= q)
            or stats["raw_enumerated_supports"] != stats["eligible_supports"] or stats["marked_orbit_enumerated_supports"]
            or stats["marked_code_orbits"] or stats["isotropic_subspace_statistics_exact"] is not True
            or stats["radical_statistics_exact"] is not True):
        raise ValueError("Incomplete raw finite census")
    radicals = stats["radical_dimension_counts"]
    if (any(not key.isdigit() or not 1 <= int(key) <= q or type(count) is not int or count < 0 for key, count in radicals.items())
            or any(type(stats[key]) is not int or stats[key] < 0 for key in
                   ("isotropic_subspaces", "nondegenerate_subspaces", "canonical_output_orbits"))
            or sum(radicals.values())+stats["nondegenerate_subspaces"] != stats["isotropic_subspaces"]):
        raise ValueError("Incorrect exact subspace or radical mass")
    keys = normalized_keys(row["output_keys"])
    if (len(keys) != len(row["output_keys"]) or len(keys) != stats["canonical_output_orbits"]
            or bool(keys) != bool(stats["nondegenerate_subspaces"])
            or keys != normalized_keys(r["output_key_words"] for r in row["witnesses"])):
        raise ValueError("Incomplete output profile or finite witnesses")
    return keys


def replay_case(case, census, native, work, workers):
    index, q, floor = case["index"], census["q"], census["minimum_distance"]
    manifest = work / f"case{index}_q{q}_d{floor}.utsp"
    output = manifest.with_suffix(".json")
    with manifest.open("wb") as stream:
        write_pointed_cases(stream, dict(maximum_protocol_length=52, cases=[dict(index=0,
            points=case["points"], ambient_dimension=case["ambient_dimension"])]))
    subprocess.run([str(native), "manifest-catalogue", "--input", str(manifest), "--output", str(output),
        "--q", str(q), "--minimum-distance", str(floor), "--workers", str(workers), "--emit-positive-supports"],
        capture_output=True, text=True, check=True)
    result = json.loads(output.read_bytes())
    if (result["status"] != "complete" or result["matrix_scope"] != "full_projective"
            or result["enumeration_mode"] != "raw_isotropic_subspaces" or result["logical_qubits"] != q
            or result["minimum_distance"] != floor
            or result["support_range"] != dict(start=0, count=1, end_exclusive=1)):
        raise ValueError("Incomplete fresh finite replay")
    for name in ("eligible_supports", "supports_processed", "raw_enumerated_supports", "marked_orbit_enumerated_supports",
                 "quotient_dimension_filtered_supports", "isotropic_subspaces", "nondegenerate_subspaces", "radical_dimension_counts"):
        if result["statistics"][name] != census["statistics"][name]:
            raise ValueError("Fresh finite replay differs on "+name)
    keys = normalized_keys(key["canonical_key_words"] for row in result["positive_supports"] for key in row["output_keys"])
    if keys != normalized_keys(census["output_keys"]):
        raise ValueError("Fresh complete output profile differs")
    return dict(case=index, q=q, minimum_distance=floor, complete_census_freshly_recomputed=True,
                isotropic_subspaces=result["statistics"]["isotropic_subspaces"],
                nondegenerate_subspaces=result["statistics"]["nondegenerate_subspaces"])


def verify(root, path, native, work, workers, replay_indices=()):
    bindings = {}
    def load(name, expected=None):
        raw = (root / name).read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if expected is not None and digest != expected:
            raise ValueError("Changed finite evidence binding")
        bindings[name] = digest
        return json.loads(raw)
    data = load(path.relative_to(root).as_posix())
    if (data["schema"] != "finite-cubic-parent-logical-censuses-v1"
            or data["protocol_length_interval"] != [51, 52] or data["minimum_distance"] != 3
            or data["matrix_scope"] != "full_projective"):
        raise ValueError("Wrong finite parent-sector scope")
    partition = load(data["source_domain"], data["source_domain_sha256"])
    if data["sources"] != partition["separate_sources"]:
        raise ValueError("The separate parent-source set is incomplete")
    source_keys = {(r["c"], r["m"], r["index"]) for r in data["sources"]}
    if len(source_keys) != len(data["sources"]) or any(c != 52 or m != 7 for c, m, _ in source_keys):
        raise ValueError("Invalid finite cubic sources")
    sources, shards = read_selected(root, source_keys)
    low_q = verify_low_q(root)
    frontier_name = "data/protocols/pareto_frontier.json"
    frontier = load(frontier_name)
    verify_frontier(root / frontier_name)
    tensors = {r["output_id"]: LogicalTensor(tuple(int(s, 2) for s in r["generator_matrix_rows"][:r["q"]]))
               for r in frontier["protocols"]}
    outputs = {(r["q"], tuple(int(w, 0) for w in r["tensor_key_words"])): r["output_id"] for r in frontier["outputs"]}
    restrictions = primitive_restrictions(load("certificates/tensors/q5_orbits.json"))
    primitive_keys = normalized_keys(r["output_key_words"] for r in restrictions["restrictions"])
    necessary = necessary_profiles(load("certificates/output_profiles/q6_necessary_profiles.json"))
    for q in (1, 2, 4):
        if any(standard_tensor(q, v << (q+q*(q-1)//2)).intrinsic() for v in range(1 << (q*(q-1)*(q-2)//6))):
            raise ValueError("Alternating small-tensor exclusion is false")
    cases, quotient_checks, candidate_indices, witnesses, replays = data["cases"], [], set(), [], []
    work.mkdir(parents=True, exist_ok=True)
    for i, case in enumerate(cases):
        points, h = checked_support(case["points"], case["ambient_dimension"]), case["ambient_dimension"]
        source = tuple(case["source"][k] for k in ("c", "m", "index"))
        if case["index"] != i or source not in source_keys or h not in (7, 8) or 0 in points:
            raise ValueError("Invalid finite stabiliser support")
        d3 = logical_label_space(points, h, 3)["quotient_dimension"]
        d4 = logical_label_space(points, h, 4)["quotient_dimension"]
        expected_branches = {(5, 3)} if h == 7 else {(3, 4)}
        if h == 7 and d4:
            raise ValueError("A nonconstant support has an unclosed higher-distance quotient")
        if h == 8:
            if (d3 != d4 or binary_rank([*evaluation_rows(points, h, 1)[1:], (1 << len(points))-1]) != h
                    or logical_label_space(points, h, 5)["quotient_dimension"]):
                raise ValueError("Constant-stabiliser or distance-five exclusion failed")
            q3, = [c for c in case["censuses"] if c["q"] == 3]
            if q3["statistics"]["nondegenerate_subspaces"]:
                expected_branches.add((5, 4))
        if {(c["q"], c["minimum_distance"]) for c in case["censuses"]} != expected_branches or len(case["censuses"]) != len(expected_branches):
            raise ValueError("Missing or repeated required finite census")
        for census in case["censuses"]:
            q, floor = census["q"], census["minimum_distance"]
            keys = check_raw_census(census, d3 if floor == 3 else d4)
            if q == 5:
                if keys & primitive_keys:
                    raise ValueError("A primitive higher-dimensional successor remains")
                if floor == 4 and keys:
                    raise ValueError("An intrinsic higher-dimensional hyperplane case remains")
                if floor == 3 and any(p.issubset(keys) for p in necessary):
                    candidate_indices.add(i)
            elif keys not in (set(), {(0x40,)}):
                raise ValueError("Incorrect alternating three-dimensional output")
            if h == 8 or i in replay_indices:
                replays.append(replay_case(case, census, native, work, workers))
            for row in census["witnesses"]:
                masks, columns = check_matrix(row["generator_matrix_rows"], q)
                tensor = LogicalTensor(masks[:q])
                distance, coefficient = distance_up_to(columns, q)
                key, = normalized_keys([row["output_key_words"]])
                identifier = outputs.get((q, key))
                if (identifier is None or not tensor.intrinsic() or distance != floor
                        or (row["n"], row["S"], row["d_Z"], row["error_coefficient"])
                        != (len(points), h+q, distance, coefficient) or len(columns) != len(points) or len(masks) != h+q):
                    raise ValueError("Incorrect finite protocol witness")
                basis = find_equivalence_basis(tensor, tensors[identifier])
                choices = [p["index"] for p in frontier["protocols"] if p["output_id"] == identifier and p["d_Z"] == distance
                           and p["n"] <= len(points) and p["S"] <= h+q]
                if basis is None or not choices:
                    raise ValueError("Incorrect output key or uncovered finite metric")
                witnesses.append(dict(case=i, q=q, output_id=identifier, dominator=min(choices), output_basis=list(basis)))
        quotient_checks.append(dict(case=i, d3=d3, d4=d4))
    q6data = load("data/protocol_sectors/length52_q6_large_quotient.json")
    q6cases = [r for r in q6data["cases"] if r["source_profile"].get("space_length") == 52]
    closure = verify_q6_closure(root)
    bindings.update(closure["dependencies"])
    canonical_inputs = [dict(index=i, points=r["points"], ambient_dimension=r["ambient_dimension"])
                        for i, r in enumerate(cases)]
    raw_cases = []
    for source, points in sorted(sources.items()):
        for row in all_pointings(points, source[1], 52):
            normal = [x for x in row["points"] if x]
            raw_cases.append(dict(source=source, n=len(row["points"]), h=row["ambient_dimension"]))
            canonical_inputs.append(dict(index=len(canonical_inputs), points=normal, ambient_dimension=row["ambient_dimension"]))
    for row in q6cases:
        canonical_inputs.append(dict(index=len(canonical_inputs), points=row["points"], ambient_dimension=row["ambient_dimension"]))
    manifest, output = work / "all_geometries.utsp", work / "all_canonical_keys.json"
    with manifest.open("wb") as stream:
        write_pointed_cases(stream, dict(maximum_protocol_length=52, cases=canonical_inputs))
    subprocess.run([str(native), "manifest-support-canonical-keys", "--input", str(manifest), "--output", str(output),
                    "--workers", str(workers)], capture_output=True, text=True, check=True)
    canonical = json.loads(output.read_bytes())
    if (canonical["status"] != "complete" or len(canonical["records"]) != len(canonical_inputs)
            or canonical["support_range"] != dict(start=0, count=len(canonical_inputs), end_exclusive=len(canonical_inputs))):
        raise ValueError("Incomplete canonical geometry comparison")
    for i, row in enumerate(canonical["records"]):
        if (row["support_index"] != i or row["source_index"] != i or row["original_record_index"] != i
                or row["support_length"] != len(canonical_inputs[i]["points"])
                or row["ambient_dimension"] != canonical_inputs[i]["ambient_dimension"]):
            raise ValueError("Incorrect canonical input identity")
    lookup = {}
    for case, row in zip(cases, canonical["records"][:len(cases)], strict=True):
        source = tuple(case["source"][k] for k in ("c", "m", "index"))
        key = source, tuple(row["key_words"])
        if key in lookup:
            raise ValueError("A finite source geometry occurs twice")
        lookup[key] = case["index"]
    covered, multiplicities = set(), Counter()
    for row, value in zip(raw_cases, canonical["records"][len(cases):len(cases)+len(raw_cases)], strict=True):
        key = row["source"], tuple(value["key_words"])
        if key not in lookup:
            raise ValueError("A raw source pointing has no complete logical census")
        index = lookup[key]
        if row["n"] < len(cases[index]["points"]) or row["h"] != cases[index]["ambient_dimension"]:
            raise ValueError("Zero-column normalisation changes the claimed geometry")
        covered.add(index)
        multiplicities[index] += 1
    if covered != set(range(len(cases))):
        raise ValueError("A finite census has no released parent pointing")
    q6bindings = []
    for row, value in zip(q6cases, canonical["records"][len(cases)+len(raw_cases):], strict=True):
        ref = row["source_profile"]
        key = (ref["space_length"], ref["affine_dimension"], ref["space_index"]), tuple(value["key_words"])
        index = lookup.get(key)
        if (index not in candidate_indices
                or normalized_keys(row["q5_output_keys"]) != normalized_keys(cases[index]["censuses"][0]["output_keys"])):
            raise ValueError("A required q6 successor is not bound to the complete predecessor profile")
        q6bindings.append(dict(case=index, q6_case_index=row["index"]))
    if {r["case"] for r in q6bindings} != candidate_indices or len(q6bindings) != len(candidate_indices):
        raise ValueError("A six-dimensional successor remains unclosed")
    return dict(schema="finite-cubic-parent-all-dimensions-verification-v1", status="pass", source_spaces=len(sources),
        raw_pointings=len(raw_cases), finite_normalised_cases=len(cases), q6_successor_bindings=q6bindings,
        all_raw_origins_reconstructed_from_compact_sources=True, all_input_equivalences_freshly_recomputed=True,
        all_distance_filtered_quotients_independently_checked=True,
        q5_h7_censuses_retained_complete_raw_results=True, higher_distance_censuses_freshly_recomputed=True,
        every_output_metric_covered_by_frontier=True, all_q_at_least7_excluded=True,
        quotient_checks=quotient_checks, raw_pointing_multiplicities=dict(multiplicities),
        witness_checks=witnesses, fresh_censuses=replays, source_shard_bindings=shards,
        dependencies=bindings, low_q_frontier_sha256=low_q["frontier_sha256"], is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--replay-case", type=int, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("Workers must be positive")
    result = verify(root, root / "data/protocol_sectors/length52_separate_cubic_censuses.json",
        args.native.resolve(), args.work_directory.resolve(), args.workers, args.replay_case)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("dependencies", "quotient_checks", "raw_pointing_multiplicities",
                                                               "source_shard_bindings", "witness_checks")}, indent=2))


if __name__ == "__main__":
    main()
