"""Compose finite protocol sectors, keeping parent-space completeness explicit."""

import argparse
import hashlib
import json
from pathlib import Path

from verify_cubic_parent_sector import verify as verify_cubic_parents
from verify_length52_large_quotients import verify as verify_large52
from verify_length52_small_blocks import verify as verify_small52
from verify_primitive_restrictions import verify as primitive_restrictions
from verify_profile_exclusions import necessary_profiles, verify as verify_profiles
from verify_protocol_case_census import normalized_keys
from verify_selected_cubic_sector import verify as verify_selected52
from verify_source_block_sector import verify as verify_blocks50
from verify_through48_sector import nondominated_metrics, verify as verify_through48


def require(value, message):
    if not value:
        raise ValueError(message)


class Evidence:
    """Hash each finite input once and reject inconsistent certificate bindings."""

    def __init__(self, root):
        self.root = Path(root).resolve()
        self.bindings = {}

    def path(self, name):
        path = (self.root / name).resolve()
        require(path.is_relative_to(self.root), "Evidence path escapes the release")
        return path

    def digest(self, name, expected=None):
        if name not in self.bindings:
            with self.path(name).open("rb") as stream:
                self.bindings[name] = hashlib.file_digest(stream, "sha256").hexdigest()
        value = self.bindings[name]
        require(expected is None or value == expected, "Changed evidence binding: "+name)
        return value

    def read(self, name, expected=None):
        self.digest(name, expected)
        return json.loads(self.path(name).read_bytes())

    def certificate(self, name, schema, flags=(), expected=None):
        data = self.read(name)
        require(data.get("schema") == schema and data.get("status") == "pass", "Incomplete certificate: "+name)
        require(all(data.get(flag) is True for flag in flags), "Missing certified implication: "+name)
        require(all(data.get(k) == v for k, v in (expected or {}).items()), "Mismatched certificate scope: "+name)
        for family in ("dependencies", "result_bindings", "supporting_data_sha256"):
            for path, digest in data.get(family, {}).items():
                self.digest(path, digest)
        return data


def check52_partition(domain, small, selected, large, separate):
    require(domain["source_spaces"] == small["totals"]["source_spaces"]+selected["source_spaces"]+separate["source_spaces"],
            "The length-52 source families do not exhaust the parent catalogue")
    require(small["source_blocks"]+selected["source_spaces"] == len(domain["main_blocks"]),
            "The main source-block partition has a gap or an overlap")
    require(small["totals"]["larger_quotient_supports"] == large["pointings"],
            "The deferred logical inputs do not match their source-block domain")
    require(small["selected_cubic_source_blocks"] == selected["source_spaces"]
            and small["separately_certified_sources"] == domain["separate_sources"]
            and len(domain["separate_sources"]) == separate["source_spaces"], "A separate source family is unclosed")
    require(selected["all_q_at_least7_excluded"] and large["all_q_at_least7_excluded"] and separate["all_q_at_least7_excluded"],
            "An arbitrary-logical-dimension sector remains open")
    return dict(source_spaces=domain["source_spaces"], main_source_blocks=len(domain["main_blocks"]),
        small_quotient_source_spaces=small["totals"]["source_spaces"], selected_cubic_sources=selected["source_spaces"],
        separate_cubic_sources=separate["source_spaces"],
        source_pointings=small["totals"]["pointed_supports"]+selected["pointings"]+separate["raw_pointings"],
        deferred_pointings_already_in_main_domain=large["pointings"])


def verify_through52(root, native, work, workers):
    evidence = Evidence(root)
    root = evidence.root
    frontier_name = "data/protocols/pareto_frontier.json"
    frontier = evidence.read(frontier_name)
    frontier_sha = evidence.digest(frontier_name)
    outputs = {(r["q"], tuple(int(w, 0) for w in r["tensor_key_words"])): r["output_id"] for r in frontier["outputs"]}
    def metric(row):
        key, = normalized_keys([row["output_key_words"]])
        identifier = outputs.get((row["q"], key))
        require(identifier is not None, "A finite metric has an unknown output")
        return identifier, row["d_Z"], row["n"], row["S"]
    source48 = evidence.read("data/protocol_sectors/through48_source_domain.json")
    source48_sha = evidence.digest("data/protocol_sectors/through48_source_domain.json")
    origin48 = evidence.certificate("certificates/through48_origin_replay.json", "finite-origin-quotient-replay-v1",
        ("all_origin_orbits_freshly_recomputed", "every_finite_point_set_agrees"),
        dict(source_domain_sha256=source48_sha, pointing_sha256=source48["pointing_sha256"],
             sources=len(source48["sources"]), pointed_supports=source48["pointing_count"]))
    evidence.digest(source48["pointing_data"], origin48["pointing_sha256"])
    evidence.certificate("certificates/through48_quotient_partition.json", "complete-quotient-domain-verification-v1",
        ("entire_quotient_partition_freshly_recomputed",), dict(source_domain_sha256=source48_sha,
        logical_census_sha256=evidence.digest("data/protocol_sectors/through48_logical_censuses.json"),
        pointing_supports=source48["pointing_count"], omitted_larger_quotient_inputs=0))
    small48 = evidence.certificate("certificates/through48_small_sources.json", "finite-small-source-closure-v1",
        ("all_origins_explicitly_enumerated", "all_required_logical_censuses_freshly_recomputed",
         "all_q_at_least5_excluded_at_distance3", "all_outputs_excluded_at_distance_at_least4", "exact_d3_q_at_most4_dominated"),
        dict(source_domain_sha256=source48_sha, maximum_protocol_length=48, low_q_frontier_sha256=frontier_sha))
    require({tuple(r["space"][k] for k in ("c", "m", "index")) for r in small48["sources"]}
            == {tuple(r[k] for k in ("c", "m", "index")) for r in source48["separate_sources"]},
            "The shorter-length separate source certificate covers different parents")
    previous48 = verify_through48(root)
    finite_name = "data/protocol_sectors/finite_shorter_source_blocks.json"
    evidence.digest(finite_name, previous48["finite_source_census_sha256"])
    evidence.certificate("certificates/finite_shorter_protocol_sources.json", "finite-pointed-source-closure-v1",
        ("all_intrinsic_logical_dimensions_exhausted", "all_five_pointing_cases_counted_independently",
         "finite_censuses_freshly_recomputed", "all_emitted_witnesses_strictly_dominated"),
        dict(data_sha256=previous48["finite_source_census_sha256"], frontier_sha256=frontier_sha))
    previous50 = verify_blocks50(root, root / "data/protocol_sectors/length50_source_domain.json",
                                root / "data/protocol_sectors/length50_blocks")
    evidence.digest("data/protocol_sectors/length50_source_domain.json", previous50["source_domain_sha256"])
    for name, digest in previous50["result_bindings"].items():
        evidence.digest(name, digest)
    require(previous48["exact_frontier_merge_agrees"] and previous50["cumulative_frontier_equals_release"],
            "The preceding cutoff has no complete finite frontier")
    require(previous48["primitive_q7_and_all_q8_excluded_given_complete_finite_censuses"]
            and previous50["all_q_at_least8_excluded_by_isotropic_subspace_containment"], "Unclosed preceding logical dimensions")
    domain_name = "data/protocol_sectors/length52_source_domain.json"
    domain = evidence.read(domain_name)
    small_name = "data/protocol_sectors/length52_small_quotient_blocks.json"
    small = verify_small52(root, root / small_name)
    evidence.digest(small_name, small["data_sha256"])
    evidence.digest(domain_name, small["source_domain_sha256"])
    selected_name = "data/protocol_sectors/length52_selected_cubic_censuses.json"
    selected = verify_selected52(root, root / selected_name)
    large = verify_large52(root, root / "data/protocol_sectors/length52_large_quotient_domain.json",
                           native, work / "large", workers)
    separate = verify_cubic_parents(root, root / "data/protocol_sectors/length52_separate_cubic_censuses.json",
                                    native, work / "cubic", workers)
    for report in (selected, large, separate):
        for name, digest in report["dependencies"].items():
            evidence.digest(name, digest)
    partition = check52_partition(domain, small, selected, large, separate)
    require(large["exact_support_equivalences_freshly_recomputed"] and separate["all_input_equivalences_freshly_recomputed"],
            "Some finite logical input equivalences are unverified")
    profile_name = "data/protocol_profiles/length52_small_quotient.json"
    profiles = evidence.read(profile_name, small["profile_sha256"])
    necessary = necessary_profiles(evidence.read("certificates/output_profiles/q6_necessary_profiles.json"))
    exclusion = verify_profiles(profiles, necessary)
    primitive = primitive_restrictions(evidence.read("certificates/tensors/q5_orbits.json"))
    primitive_keys = normalized_keys(row["output_key_words"] for row in primitive["restrictions"])
    require(exclusion["all_q6_absent_on_listed_profiles"] and not any(
        normalized_keys(row["output_keys"]) & primitive_keys for row in profiles["profiles"]),
        "The small-quotient higher-logical-dimension sector is unclosed")
    candidates = set()
    finite_small = evidence.read(small_name)
    for block in finite_small["blocks"]:
        for branch in block["branches"]:
            candidates.update(metric(row) for row in branch["witnesses"])
    del finite_small
    finite_selected = evidence.read(selected_name)
    for q in (5, 6):
        for row in finite_selected[f"q{q}_profiles"]:
            geometry = selected["quotient_checks"][row["pointing_index"]]
            for key in normalized_keys(row["output_keys"]):
                candidates.add((outputs[q, key], 3, geometry["n"], q+geometry["h"]))
    candidates.update(metric(row) for row in finite_selected["finite_witnesses"])
    candidates.update((row["output_id"], 3, row["n"], row["S"]) for row in large["q5_metric_checks"])
    finite_q6 = evidence.read("data/protocol_sectors/length52_q6_large_quotient.json")
    for case in finite_q6["cases"]:
        for key in normalized_keys(case["expected"]["output_keys"]):
            candidates.add((outputs[6, key], 3, len(case["points"]), 6+case["ambient_dimension"]))
    for case in evidence.read("data/protocol_sectors/length52_separate_cubic_censuses.json")["cases"]:
        for census in case["censuses"]:
            candidates.update(metric(row) for row in census["witnesses"])
    previous = {(r["output_id"], r["d_Z"], r["n"], r["S"]) for r in frontier["protocols"] if r["n"] <= 50}
    expected = {(r["output_id"], r["d_Z"], r["n"], r["S"]) for r in frontier["protocols"] if r["n"] <= 52}
    require(all(51 <= r[2] <= 52 for r in candidates), "A new finite metric is outside its length interval")
    merged = nondominated_metrics(previous | candidates)
    require(merged == expected, "The finite sector union does not reproduce the released through-52 frontier")
    for source in (previous48, previous50, selected, large, separate):
        for name, digest in source.get("source_shard_bindings", {}).items():
            evidence.digest(name, digest)
    return dict(schema="composed-protocol-coverage-v1", status="pass", maximum_protocol_length=52,
        minimum_exact_distance=3, matrix_scope="full_projective", output_equivalence="CNOT+S",
        logical_dimension_scope="all intrinsic dimensions", pareto_objectives=["n", "S"], exact_distance_is_part_of_key=True,
        preceding_frontier_counts={"48": previous48["merged_frontier_points"], "50": previous50["cumulative_frontier_points"]},
        length52_partition=partition, finite_candidate_metrics=len(candidates), cumulative_frontier_points=len(merged),
        cumulative_frontier_equals_release=True, all_finite_source_and_logical_sectors_composed=True,
        every_successor_obligation_closed=True, all_logical_dimensions_covered=True,
        small_quotient_q6_and_primitive_q7_exclusions_verified=True,
        zero_column_extensions_excluded_at_fixed_exact_distance=True,
        geometry_and_sector_consistency_checks_recomputed=True, all_logical_censuses_freshly_recomputed=False,
        conditional_premises=["Completeness of the affine unital-support catalogues through length 52",
                              "Completeness of each supplied finite raw or marked logical census on its explicit input domain"],
        parent_space_completeness_is_separate=True, is_global_completeness_certificate=False,
        frontier_metrics=[dict(output_id=o, d_Z=d, n=n, S=s) for o, d, n, s in sorted(merged)],
        dependencies=evidence.bindings)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(args.workers >= 1, "Workers must be positive")
    result = verify_through52(root, args.native.resolve(), args.work_directory.resolve(), args.workers)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("dependencies", "frontier_metrics")}, indent=2))


if __name__ == "__main__":
    main()
