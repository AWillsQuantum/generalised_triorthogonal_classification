"""Bind the length-52 larger-quotient geometry and complete logical output profiles."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from logical_spaces import logical_label_space
from selected_spaces import read_selected
from space_codec import checked_support
from verify_finite_sector_closure import verify as verify_q6_closure
from verify_large_quotient_domain import canonical_check, pointed_support
from verify_primitive_restrictions import verify as primitive_restrictions
from verify_primitive_support_sector import support_key, validate_zero_result
from verify_profile_exclusions import necessary_profiles, verify as verify_profiles
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_frontier


def check_class(row, index, primitive, profile):
    points, h = row["points"], row["ambient_dimension"]
    checked_support(sorted(points), h)
    if (row["index"] != index or row["output_profile_index"] != index or 0 in points
            or row["support_length"] != len(points) or not 51 <= len(points) <= 52
            or support_key(points, h, row["key_words"]) != row["canonical_key_sha256"]
            or type(row["automorphism_group_order"]) is not int or row["automorphism_group_order"] <= 0):
        raise ValueError("Invalid canonical stabiliser geometry")
    if (profile["index"] != index or profile["canonical_geometry_sha256"] != row["canonical_key_sha256"]
            or type(profile["nondegenerate_subspaces"]) is not int or profile["nondegenerate_subspaces"] < 0
            or bool(profile["nondegenerate_subspaces"]) != bool(profile["output_keys"])
            or len(normalized_keys(profile["output_keys"])) != len(profile["output_keys"])):
        raise ValueError("The complete output profile is not bound to its support")
    d3 = logical_label_space(sorted(points), h, 3)["quotient_dimension"]
    d4 = logical_label_space(sorted(points), h, 4)["quotient_dimension"]
    if d3 < 15 or profile["output_keys"] and d4 >= 5:
        raise ValueError("The q5 profile has an unclosed higher-distance case")
    if ("primitive_union_case_index" in row) == ("primitive_target_census" in row):
        raise ValueError("Expected exactly one primitive-target census binding")
    if "primitive_union_case_index" in row:
        reference = row["primitive_union_case_index"]
        if not 0 <= reference < len(primitive["cases"]):
            raise ValueError("Invalid primitive-target class index")
        shared = primitive["cases"][reference]
        if any(row[key] != shared[key] for key in ("points", "ambient_dimension", "key_words", "canonical_key_sha256")):
            raise ValueError("Primitive-target census has different support coordinates")
        census = shared["result"]
        if census["subspace_orbits"] or census["weighted_subspace_count"]:
            raise ValueError("The primitive-target census is not empty")
    else:
        census = row["primitive_target_census"]
        validate_zero_result(census)
        if census["support_automorphism_group_order"] != row["automorphism_group_order"]:
            raise ValueError("Primitive-target support symmetry differs")
    if census["quotient_dimension"] != d3:
        raise ValueError("Primitive-target logical quotient dimension differs")
    return dict(index=index, d3=d3, d4=d4, q5_positive=bool(profile["output_keys"]))


def verify(root, path, native=None, work=None, workers=1):
    bindings = {}
    def load(name, expected=None):
        target = (root / name).resolve()
        if not target.is_relative_to(root):
            raise ValueError("A dependency escapes the public release")
        raw = target.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if expected is not None and expected != digest:
            raise ValueError("A finite dataset binding changed")
        bindings[name] = digest
        return json.loads(raw)
    data = load(path.relative_to(root).as_posix())
    if (data["schema"] != "finite-length52-large-quotient-domain-v1"
            or data["protocol_length_interval"] != [51, 52] or data["minimum_distance"] != 3
            or data["minimum_quotient_dimension"] != 15 or data["matrix_scope"] != "full_projective"):
        raise ValueError("Incorrect larger-quotient scope")
    domain = load(data["source_domain"], data["source_domain_sha256"])
    small = load("certificates/length52_small_quotient_blocks.json")
    if (small["status"] != "pass" or small["source_domain_sha256"] != data["source_domain_sha256"]
            or not small["complete_source_intervals_bound"]):
        raise ValueError("The complete source and low-q higher-distance cover is missing")
    block_data = load("data/protocol_sectors/length52_small_quotient_blocks.json", small["data_sha256"])
    blocks = {row["index"]: row for row in domain["main_blocks"]}
    expected_counts = Counter({i: row["larger_quotient_supports"] for i, row in blocks.items()
                               if row.get("larger_quotient_supports")})
    if any(blocks[i]["logical_cover"] != "raw_small_quotients" for i in expected_counts):
        raise ValueError("Large quotients must belong to the raw finite source cover")
    del block_data
    primitive = load(data["primitive_domain"], data["primitive_domain_sha256"])
    primitive_certificate = load("certificates/primitive_support_sector.json")
    if (primitive_certificate["status"] != "pass"
            or primitive_certificate["input_sha256"] != data["primitive_domain_sha256"]
            or not primitive_certificate["normalisation_equivalences_freshly_recomputed"]):
        raise ValueError("Missing verified shared primitive-target domain")
    profiles = load(data["profile_dataset"])
    if len(profiles["profiles"]) != len(data["classes"]):
        raise ValueError("Incomplete finite class-profile correspondence")
    necessary = necessary_profiles(load("certificates/output_profiles/q6_necessary_profiles.json"))
    profile_check = verify_profiles(profiles, necessary)
    restriction = primitive_restrictions(load("certificates/tensors/q5_orbits.json"))
    primitive_keys = normalized_keys(row["output_key_words"] for row in restriction["restrictions"])
    classes, checks, seen_keys, dimensions = data["classes"], [], set(), Counter()
    for index, row in enumerate(classes):
        profile = profiles["profiles"][index]
        value = check_class(row, index, primitive, profile)
        if row["canonical_key_sha256"] in seen_keys or normalized_keys(profile["output_keys"]) & primitive_keys:
            raise ValueError("Repeated class or contradiction with primitive-target absence")
        seen_keys.add(row["canonical_key_sha256"])
        checks.append(value)
        dimensions[value["d3"], value["d4"]] += 1
    keys = {(52, row["source"]["m"], row["source"]["index"]) for row in data["pointings"]}
    supports, shards = read_selected(root, keys)
    counts, used, seen, cache, canonical_inputs = Counter(), set(), set(), {}, []
    for index, row in enumerate(data["pointings"]):
        source = row["source"]
        key = source["c"], source["m"], source["index"]
        block = blocks[row["source_block_index"]]
        if (row["index"] != index or key[0] != 52 or key[1] != block["m"]
                or not block["source_interval"][0] <= key[2] < block["source_interval"][1]):
            raise ValueError("Pointing is outside its source interval")
        points, h = pointed_support(supports[key], key[1], row["parity_case"], row["origin"])
        identity = key, h, tuple(points)
        if (identity in seen or not 51 <= len(points) <= 52
                or (len(points), h) != (row["protocol_length"], row["ambient_dimension"])):
            raise ValueError("Repeated pointing or incorrect geometric parameters")
        seen.add(identity)
        normal = tuple(x for x in points if x)
        geometry = h, normal
        if geometry not in cache:
            cache[geometry] = logical_label_space(normal, h, 3)["quotient_dimension"]
        case_index = row["support_class_index"]
        case = classes[case_index]
        if (case["index"] != case_index or (len(normal), h) != (case["support_length"], case["ambient_dimension"])
                or cache[geometry] != row["quotient_dimension"] or cache[geometry] != checks[case_index]["d3"]):
            raise ValueError("Pointing and canonical class have different logical parameters")
        canonical_inputs.append(dict(index=index, points=list(normal), ambient_dimension=h, support_class_index=case_index))
        counts[row["source_block_index"]] += 1
        used.add(case_index)
    if counts != expected_counts or used != set(range(len(classes))):
        raise ValueError("The finite pointing domain does not cover every selected input and canonical class")
    candidate_indices = {row["index"] for row in profile_check["candidates"]}
    q6data = load("data/protocol_sectors/length52_q6_large_quotient.json")
    q6certificate = load("certificates/q6_large_quotient_census.json")
    if (q6certificate["status"] != "pass" or not q6certificate["every_support_output_profile_recomputed"]
            or q6certificate["input_sha256"] != bindings["data/protocol_sectors/length52_q6_large_quotient.json"]):
        raise ValueError("A complete fresh q6 census is missing")
    selected_q6 = {}
    for row in q6data["cases"]:
        reference = row["source_profile"]
        if reference.get("dataset") != profiles["id"]:
            continue
        index = reference["index"]
        if (index in selected_q6 or index not in candidate_indices
                or (row["points"], row["ambient_dimension"]) != (sorted(classes[index]["points"]), classes[index]["ambient_dimension"])
                or normalized_keys(row["q5_output_keys"]) != normalized_keys(profiles["profiles"][index]["output_keys"])):
            raise ValueError("A q6 successor is not bound to its predecessor geometry and profile")
        selected_q6[index] = row["index"]
    if set(selected_q6) != candidate_indices:
        raise ValueError("A required q6 successor has not been enumerated")
    closure = verify_q6_closure(root)
    bindings.update(closure["dependencies"])
    frontier_name = "data/protocols/pareto_frontier.json"
    frontier = load(frontier_name)
    verify_frontier(root / frontier_name)
    outputs = {(row["q"], tuple(int(w, 0) for w in row["tensor_key_words"])): row["output_id"] for row in frontier["outputs"]}
    metrics = []
    for row, profile in zip(classes, profiles["profiles"], strict=True):
        for key in sorted(normalized_keys(profile["output_keys"])):
            identifier = outputs.get((5, key))
            choices = [p["index"] for p in frontier["protocols"] if p["output_id"] == identifier and p["d_Z"] == 3
                       and p["n"] <= row["support_length"] and p["S"] <= 5+row["ambient_dimension"]]
            if not choices:
                raise ValueError("A complete q5 profile gives an uncovered exact-distance metric")
            metrics.append(dict(support_class_index=row["index"], output_id=identifier, d_Z=3,
                                n=row["support_length"], S=5+row["ambient_dimension"], dominator=min(choices)))
    canonical = dict(exact_support_equivalences_freshly_recomputed=False)
    if native is not None:
        work.mkdir(parents=True, exist_ok=True)
        canonical = canonical_check(dict(pointings=canonical_inputs, classes=classes), native, work, workers)
    return dict(schema="length52-large-quotient-sector-verification-v1", status="pass",
        protocol_length_interval=[51, 52], pointings=len(data["pointings"]), classes=len(classes), source_spaces=len(keys),
        primitive_shared_classes=sum("primitive_union_case_index" in row for row in classes),
        primitive_additional_classes=sum("primitive_target_census" in row for row in classes),
        q5_positive_profiles=sum(row["q5_positive"] for row in checks),
        q5_weighted_nondegenerate_subspaces=sum(row["nondegenerate_subspaces"] for row in profiles["profiles"]),
        distance_filtration=[dict(d3=a, d4=b, count=count) for (a, b), count in sorted(dimensions.items())],
        source_block_cardinalities_verified=True, every_pointing_reconstructed_from_compact_sources=True,
        every_quotient_dimension_independently_checked=True, all_positive_q5_profiles_have_exact_distance_three=True,
        q5_profiles_complete_finite_census_premise=True, q5_censuses_freshly_recomputed=False,
        primitive_zero_censuses_freshly_recomputed=False, q6_successor_case_indices=selected_q6,
        all_q_at_least7_excluded=True, q_at_most4_higher_distance_covered_by_bound_source_blocks=True,
        every_q5_output_metric_covered_by_frontier=True, quotient_checks=checks, q5_metric_checks=metrics,
        dependencies=bindings, source_shard_bindings=shards, **canonical, is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--work-directory", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1 or (args.native is None) != (args.work_directory is None):
        parser.error("Native equivalence checking needs an executable, work directory and positive workers")
    result = verify(root, root / "data/protocol_sectors/length52_large_quotient_domain.json",
                    args.native.resolve() if args.native else None,
                    args.work_directory.resolve() if args.work_directory else None, args.workers)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in
        ("dependencies", "source_shard_bindings", "quotient_checks", "q5_metric_checks")}, indent=2))


if __name__ == "__main__":
    main()
