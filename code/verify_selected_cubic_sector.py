"""Verify complete selected-cubic output profiles and their all-dimension closure."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from logical_spaces import logical_label_space
from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis, gate_tensor
from protocol_domain import open_domain, read_header, read_pointings, reconstruct, source_key
from selected_spaces import read_selected
from space_codec import binary_rank
from space_lifts import evaluation_rows
from verify_low_q_pruning import verify as verify_small_outputs
from verify_primitive_restrictions import verify as primitive_restrictions
from verify_profile_exclusions import necessary_profiles
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_frontier
from verify_tensor_authorities import standard_tensor


def check_profile(row, q, quotient_dimension):
    stats = row["statistics"]
    keys = normalized_keys(row["output_keys"])
    if (row["enumeration_mode"] != "hybrid_raw_and_marked_code_orbits"
            or stats["source_spaces"] != 1 or stats["supports_processed"] != 1
            or stats["quotient_dimension_filtered_supports"]
            or stats["eligible_supports"] != int(quotient_dimension >= q)
            or stats["raw_enumerated_supports"] or stats["marked_orbit_enumerated_supports"] != 1
            or stats["isotropic_subspace_statistics_exact"] is not False
            or stats["radical_statistics_exact"] is not False or stats["isotropic_subspaces"] is not None
            or stats["radical_dimension_counts"]):
        raise ValueError("Incomplete marked input domain or misleading raw-count claim")
    bits = q+q*(q-1)//2+q*(q-1)*(q-2)//6
    if (any(len(key) != 1 or not 0 < key[0] < 1 << bits for key in keys)
            or len(keys) != len(row["output_keys"]) or len(keys) != stats["canonical_output_orbits"]
            or bool(keys) != bool(stats["nondegenerate_subspaces"])
            or bool(keys) != bool(stats["marked_code_orbits"])
            or not 0 <= stats["marked_code_orbits"] <= stats["nondegenerate_subspaces"]):
        raise ValueError("Inconsistent complete output profile or marked subspace mass")
    return keys


def geometry(points, h):
    normal = tuple(x for x in points if x)
    d3 = logical_label_space(normal, h, 3)["quotient_dimension"]
    d4 = logical_label_space(normal, h, 4)["quotient_dimension"]
    stabilisers = evaluation_rows(normal, h, 1)[1:]
    constant = binary_rank([*stabilisers, (1 << len(normal))-1]) == h
    if d4 and not constant:
        raise ValueError("A positive distance-four quotient has no constant stabiliser")
    d5 = logical_label_space(normal, h, 5)["quotient_dimension"] if d4 else 0
    if d5:
        raise ValueError("Distance five or higher has not been excluded")
    return dict(d3=d3, d4=d4, d5=d5, constant_stabiliser=constant)


def verify(root, path):
    bindings = {}
    def read(name):
        file = root / name
        raw = file.read_bytes()
        bindings[name] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    name = path.relative_to(root).as_posix()
    data = read(name)
    if (data["schema"] != "finite-selected-cubic-protocol-censuses-v1"
            or data["protocol_length_interval"] != [51, 52] or data["matrix_scope"] != "full_projective"):
        raise ValueError("Incorrect finite selected-cubic scope")
    domain = read(data["source_domain"])
    if (bindings[data["source_domain"]] != data["source_domain_sha256"]
            or domain["schema"] != "finite-protocol-source-domain-v1"
            or domain["pointing_count"] != data["pointing_count"]
            or (domain["minimum_protocol_length"], domain["maximum_protocol_length"]) != (51, 52)):
        raise ValueError("Finite census is bound to a different pointed domain")
    partition = read(domain["source_partition"])
    if bindings[domain["source_partition"]] != domain["source_partition_sha256"]:
        raise ValueError("Changed full source partition")
    expected_sources = {row["index"]: row for row in partition["main_blocks"] if row["logical_cover"] == "selected_cubic_supports"}
    actual_sources = {}
    for index, source in enumerate(domain["sources"]):
        block_index = source["source_block_index"]
        if source["index"] != index or block_index in actual_sources or block_index not in expected_sources:
            raise ValueError("Unknown or repeated selected source block")
        block = expected_sources[block_index]
        if (source["space"] != dict(c=52, m=7, index=block["source_interval"][0])
                or block["source_interval"][1] != block["source_interval"][0]+1
                or domain["source_pointing_counts"][index] != block["pointed_supports"]):
            raise ValueError("Cubic source catalogue interval differs")
        actual_sources[block_index] = source
    if set(actual_sources) != set(expected_sources):
        raise ValueError("Selected cubic source set is incomplete")
    origin = read("certificates/length52_selected_cubic_origins.json")
    if (origin["status"] != "pass" or origin["source_domain_sha256"] != data["source_domain_sha256"]
            or origin["pointing_sha256"] != domain["pointing_sha256"]
            or origin["sources"] != len(domain["sources"]) or origin["pointed_supports"] != domain["pointing_count"]
            or not origin["all_origin_orbits_freshly_recomputed"] or not origin["every_finite_point_set_agrees"]):
        raise ValueError("A complete source origin-orbit certificate is missing")
    primitive = read(data["primitive_domain"])
    primitive_certificate = read("certificates/primitive_support_sector.json")
    if (bindings[data["primitive_domain"]] != data["primitive_domain_sha256"]
            or primitive_certificate["status"] != "pass"
            or primitive_certificate["input_sha256"] != data["primitive_domain_sha256"]
            or not primitive_certificate["normalisation_equivalences_freshly_recomputed"]):
        raise ValueError("Unbound primitive-target support equivalences")
    necessary = necessary_profiles(read("certificates/output_profiles/q6_necessary_profiles.json"))
    gates = read("certificates/output_profiles/certificate.json")
    for filename, digest in gates["files"].items():
        read("certificates/output_profiles/"+filename)
        if bindings["certificates/output_profiles/"+filename] != digest:
            raise ValueError("Changed finite tensor-profile gate")
    gate = read("certificates/output_profiles/q7_anchored_exclusions.json")
    if (gates["status"] != "pass" or gate["status"] != "pass"
            or gate["extensions_per_anchor"] != 1 << 22
            or gate["candidates_tested"] != gate["q6_target_key_count"]*(1 << 22)):
        raise ValueError("Incomplete seven-dimensional extension domain")
    excluded_q6 = [normalized_keys([[key] for key in row["q6_keys"]]) for row in gate["support_profiles"]
                   if row["surviving_anchored_extensions"] == 0]
    primitives = primitive_restrictions(read("certificates/tensors/q5_orbits.json"))
    primitive_keys = normalized_keys(row["output_key_words"] for row in primitives["restrictions"])
    for q in (1, 2, 4):
        shift = q+q*(q-1)//2
        if any(standard_tensor(q, value << shift).intrinsic() for value in range(1 << (q*(q-1)*(q-2)//6))):
            raise ValueError("Alternating small-output exclusion is false")
    frontier_name = "data/protocols/pareto_frontier.json"
    frontier = read(frontier_name)
    verify_frontier(root / frontier_name)
    verify_small_outputs(root)
    outputs = {(row["q"], tuple(int(word, 0) for word in row["tensor_key_words"])): row["output_id"]
               for row in frontier["outputs"]}
    tensors = {row["output_id"]: LogicalTensor(tuple(int(s, 2) for s in row["generator_matrix_rows"][:row["q"]]))
               for row in frontier["protocols"]}
    def dominator(q, keys, n, h, distance):
        records = []
        for key in sorted(keys):
            identifier = outputs.get((q, key))
            choices = [row["index"] for row in frontier["protocols"] if row["output_id"] == identifier
                       and row["d_Z"] == distance and row["n"] <= n and row["S"] <= q+h]
            if not choices:
                raise ValueError("A complete output profile gives an uncovered Pareto metric")
            records.append(dict(output_id=identifier, dominator=min(choices)))
        return records
    q6_records = {row["pointing_index"]: row for row in data["q6_profiles"]}
    if len(q6_records) != len(data["q6_profiles"]) or len(data["q5_profiles"]) != domain["pointing_count"]:
        raise ValueError("Repeated or incomplete selected profile set")
    supports, shard_bindings = read_selected(root, {source_key(row) for row in domain["sources"]})
    pointing_path = root / domain["pointing_data"]
    with pointing_path.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != domain["pointing_sha256"]:
            raise ValueError("Changed selected pointed domain")
    cache, geometry_rows, profile_checks, selected = {}, [], [], set()
    counts, d4_counts, sources_seen = Counter(), Counter(), Counter()
    if len(data["primitive_normalisation_indices"]) != domain["pointing_count"]:
        raise ValueError("Incomplete primitive input map")
    with open_domain(pointing_path) as stream:
        header = read_header(stream)
        if header != dict(source_count=len(domain["sources"]), count=domain["pointing_count"], maximum_protocol_length=52):
            raise ValueError("Pointing header disagrees with the finite source domain")
        for index, record in enumerate(read_pointings(stream, header)):
            points, h = reconstruct(domain["sources"], supports, record)
            if not 51 <= len(points) <= 52:
                raise ValueError("Pointing is outside the exact protocol interval")
            normal = tuple(x for x in points if x)
            key = h, normal
            if key not in cache:
                cache[key] = geometry(points, h)
            value = cache[key]
            sources_seen[record[0]] += 1
            d4_counts[value["d4"]] += 1
            reference = data["primitive_normalisation_indices"][index]
            if not 0 <= reference < len(primitive["normalisation_inputs"]):
                raise ValueError("Invalid primitive-target normalisation index")
            mapping = primitive["normalisation_inputs"][reference]
            if (mapping["index"] != reference or mapping["family"] != "c52_m07_pointings"
                    or (mapping["ambient_dimension"], tuple(mapping["points"])) != key):
                raise ValueError("Primitive-target domain omits this exact normalised input")
            primitive_result = primitive["cases"][mapping["case_index"]]["result"]
            if (primitive_result["subspace_orbits"] or primitive_result["weighted_subspace_count"]
                    or primitive_result["quotient_dimension"] != value["d3"]):
                raise ValueError("A primitive five-dimensional output remains possible")
            row = data["q5_profiles"][index]
            if row["pointing_index"] != index:
                raise ValueError("The complete q5 profiles are not indexed by their pointed inputs")
            keys = check_profile(row, 5, value["d3"])
            if keys & primitive_keys:
                raise ValueError("Q5 profile contradicts primitive-target absence")
            matches = [i for i, profile in enumerate(necessary) if profile.issubset(keys)]
            if matches:
                selected.add(index)
                if index not in q6_records:
                    raise ValueError("A required six-dimensional census is missing")
            if keys and value["d4"]:
                raise ValueError("A positive q5 support has an unclosed higher exact distance")
            current = dict(pointing_index=index, primitive_normalisation_index=reference,
                           q6_necessary_profile_indices=matches, q5_metrics=dominator(5, keys, len(points), h, 3))
            if index in q6_records:
                q6_keys = check_profile(q6_records[index], 6, value["d3"])
                if value["d4"] or not any(q6_keys.issubset(profile) for profile in excluded_q6):
                    raise ValueError("Higher-distance or higher-dimensional q6 successors are unclosed")
                current["q6_metrics"] = dominator(6, q6_keys, len(points), h, 3)
                counts["q6_nondegenerate_subspaces"] += q6_records[index]["statistics"]["nondegenerate_subspaces"]
            if value["d4"]:
                current["d4_CCZ_metric"] = dominator(3, {(0x40,)}, len(points), h, 4)
                if find_equivalence_basis(gate_tensor(3, "CCZ123"), tensors[outputs[3, (0x40,)]]) is None:
                    raise ValueError("Incorrect alternating three-dimensional output identification")
            counts["q5_nondegenerate_subspaces"] += row["statistics"]["nondegenerate_subspaces"]
            counts["q5_positive_supports"] += bool(keys)
            geometry_rows.append(dict(index=index, n=len(points), h=h, **value))
            profile_checks.append(current)
    if (selected != set(q6_records) or len(geometry_rows) != domain["pointing_count"]
            or [sources_seen[i] for i in range(len(domain["sources"]))] != domain["source_pointing_counts"]):
        raise ValueError("Incomplete selected successor or source-pointing domain")
    witnesses = []
    for index, row in enumerate(data["finite_witnesses"]):
        q = row["q"]
        masks, columns = check_matrix(row["generator_matrix_rows"], q)
        tensor = LogicalTensor(masks[:q])
        key, = normalized_keys([row["output_key_words"]])
        identifier = outputs.get((q, key))
        basis = find_equivalence_basis(tensor, tensors[identifier]) if identifier is not None else None
        distance, coefficient = distance_up_to(columns, q)
        if (not tensor.intrinsic() or basis is None or distance != 3 or q not in (5, 6)
                or (row["n"], row["S"], row["d_Z"], row["error_coefficient"])
                != (len(columns), len(masks), distance, coefficient)):
            raise ValueError("Invalid retained finite metric witness")
        witnesses.append(dict(index=index, output_basis=list(basis),
                             metric=dominator(q, {key}, len(columns), len(masks)-q, distance)))
    return dict(schema="selected-cubic-all-logical-dimensions-verification-v1", status="pass",
        protocol_length_interval=[51, 52], source_spaces=len(domain["sources"]), pointings=len(geometry_rows),
        q5_profiles=len(data["q5_profiles"]), q6_profiles=len(q6_records), counts=dict(counts),
        distance_four_quotient_histogram=dict(sorted(d4_counts.items())),
        all_distance_five_quotients_zero=True, all_positive_distance_four_quotients_have_constant_stabiliser=True,
        all_q_at_most4_higher_distance_candidates_pareto_dominated=True,
        all_q5_q6_positive_supports_have_exact_distance_three=True,
        all_q_at_least7_excluded=True, primitive_q5_absence_uses_bound_finite_target_censuses=True,
        finite_enumerations_freshly_recomputed=False, every_quotient_and_primitive_input_binding_recomputed=True,
        every_output_profile_metric_covered_by_frontier=True, finite_witnesses_checked=len(witnesses),
        quotient_checks=geometry_rows, profile_checks=profile_checks, witnesses=witnesses,
        source_shard_bindings=shard_bindings, dependencies=bindings, is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/length52_selected_cubic_censuses.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(root, args.input)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("quotient_checks", "profile_checks", "witnesses", "source_shard_bindings", "dependencies")}, indent=2))


if __name__ == "__main__":
    main()
