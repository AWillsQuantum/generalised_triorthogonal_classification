"""Compose the complete protocol cover through length 54 from finite mathematical evidence."""

import argparse
import json
from pathlib import Path

from verify_chain_sector import verify as verify_chain
from verify_finite_sector_closure import verify as verify_finite_closure
from verify_high_dimensional_sector import verify as verify_high
from verify_large_quotient_domain import verify as verify_large, canonical_check as canonical_large
from verify_length54_source_partition import verify as verify_partition
from verify_low_dimensional_blocks import verify as verify_low
from verify_nonprimitive_subsector import verify as verify_marked
from verify_pointed_sources import verify as verify_pointed
from verify_primitive_restrictions import verify as primitive_restrictions
from verify_primitive_support_sector import verify as verify_primitive, canonical_check as canonical_primitive
from verify_profile_exclusions import necessary_profiles, verify as verify_profiles
from verify_profile_refinement import verify as verify_refinement
from verify_protocol_case_census import normalized_keys
from verify_protocol_cover import Evidence, require
from verify_protocols import verify as verify_frontier
from verify_through48_sector import nondominated_metrics


def verify(root, native, work, workers):
    evidence = Evidence(root)
    root = evidence.root
    frontier_name = "data/protocols/pareto_frontier.json"
    frontier = evidence.read(frontier_name)
    verify_frontier(root / frontier_name)
    previous = evidence.certificate("certificates/through52_protocol_cover.json", "composed-protocol-coverage-v1",
        ("cumulative_frontier_equals_release", "all_finite_source_and_logical_sectors_composed",
         "every_successor_obligation_closed", "all_logical_dimensions_covered", "exact_distance_is_part_of_key"),
        dict(maximum_protocol_length=52, minimum_exact_distance=3, matrix_scope="full_projective",
             output_equivalence="CNOT+S", pareto_objectives=["n", "S"], logical_dimension_scope="all intrinsic dimensions"))
    require(previous["dependencies"][frontier_name] == evidence.digest(frontier_name), "The preceding frontier differs")
    previous_metrics = {(r["output_id"], r["d_Z"], r["n"], r["S"]) for r in previous["frontier_metrics"]}
    require(previous_metrics == {(r["output_id"], r["d_Z"], r["n"], r["S"]) for r in frontier["protocols"] if r["n"] <= 52},
            "The preceding complete frontier has different metrics")
    low_name = "data/protocol_sectors/length54_low_dimensional_blocks.json"
    high_name = "data/protocol_sectors/length54_high_dimensional.json"
    low_report = verify_low(root, root / low_name)
    high_report = verify_high(root, root / high_name)
    low = evidence.read(low_name, low_report["data_sha256"])
    high = evidence.read(high_name, high_report["data_sha256"])
    partition = verify_partition(root)
    require(partition["source_intervals_disjoint_and_exhaustive"]
            and high_report["logical_dimensions_at_least_seven_excluded_by_common_isotropic_containment"],
            "An entire source dimension is unclosed")
    for name, digest in partition["supporting_files"].items():
        evidence.digest(name, digest)
    necessary = necessary_profiles(evidence.read("certificates/output_profiles/q6_necessary_profiles.json"))
    primitive_tensor = primitive_restrictions(evidence.read("certificates/tensors/q5_orbits.json"))
    primitive_keys = normalized_keys(row["output_key_words"] for row in primitive_tensor["restrictions"])
    profiles = evidence.read(low["profile_dataset"], low_report["output_profile_sha256"])
    small_exclusion = verify_profiles(profiles, necessary)
    require(not any(normalized_keys(r["output_keys"]) & primitive_keys for r in profiles["profiles"]),
            "Primitive seven-dimensional small-quotient outputs are not excluded")
    refined = verify_refinement(root, native, work / "refinement")
    for name, digest in refined["dependencies"].items():
        evidence.digest(name, digest)
    candidates = evidence.read("data/protocol_sectors/length54_q6_candidates.json")
    pointed = verify_pointed(root, candidates, low)
    require(pointed["exact_source_block_positions_checked"] and pointed["cases"] == refined["candidate_supports"],
            "The refined logical candidates are not tied to exact source-block positions")
    require(small_exclusion["excluded"]+refined["refined_blocks"] == len(low["blocks"]),
            "Some small-quotient source block has no higher-dimensional closure")
    finite = verify_finite_closure(root)
    for name, digest in finite["dependencies"].items():
        evidence.digest(name, digest)
    q6_name = "data/protocol_sectors/length54_q6_small_quotient.json"
    q6 = evidence.read(q6_name)
    require(len(q6["cases"]) == refined["canonical_support_classes"], "The complete q6 domain differs from its predecessor cover")

    domain_name = "data/protocol_sectors/large_quotient_domain.json"
    primitive_name = "data/protocol_sectors/primitive_support_union.json"
    domain = evidence.read(domain_name)
    primitive = evidence.read(primitive_name)
    primitive_report = verify_primitive(primitive)
    primitive_dir = work / "primitive"
    primitive_dir.mkdir(parents=True, exist_ok=True)
    primitive_report.update(canonical_primitive(primitive, native, primitive_dir, workers))
    large_report = verify_large(root, domain, primitive, low)
    large_dir = work / "large"
    large_dir.mkdir(parents=True, exist_ok=True)
    large_report.update(canonical_large(domain, native, large_dir, workers))
    require(large_report["selected_pointings"] == low_report["totals"]["larger_quotient_supports"],
            "The large-quotient domain does not close the source-block tails")
    marked_name = "data/protocol_sectors/nonprimitive_marked_subsector.json"
    marked = evidence.read(marked_name)
    marked_report = verify_marked(root, marked)
    evidence.certificate("certificates/nonprimitive_marked_subsector.json", "finite-nonprimitive-marked-sector-verification-v1",
        ("marked_enumerations_freshly_recomputed", "all_finite_candidates_covered_by_frontier"),
        dict(input_sha256=evidence.digest(marked_name), cases=len(marked["cases"]),
             orbits=marked_report["orbits"], frontier_sha256=evidence.digest(frontier_name)))
    marked_profiles = evidence.read(marked["profile_dataset"])
    marked_exclusion = verify_profiles(marked_profiles, necessary)
    require(marked_exclusion["all_q6_absent_on_listed_profiles"] and not any(
        normalized_keys(r["output_keys"]) & primitive_keys for r in marked_profiles["profiles"]),
        "Higher logical dimensions on the marked subsector remain open")
    chain_name = "data/protocol_sectors/nonprimitive_chain_sector.json"
    chain = evidence.read(chain_name)
    chain_report = verify_chain(root, root / chain_name)
    for name, digest in chain_report["dependencies"].items():
        evidence.digest(name, digest)
    require(chain_report["exact_input_partition_checked"] and chain_report["higher_dimension_profile_implication_checked"]
            and chain_report["q5_input_classes"]+chain_report["zero_column_dominated_classes"]+chain_report["separate_finite_classes"]
            == large_report["support_classes"], "The complete large-quotient logical partition has a gap or overlap")

    outputs = {(r["q"], tuple(int(w, 0) for w in r["tensor_key_words"])): r["output_id"] for r in frontier["outputs"]}
    metrics = set()
    def add(row):
        key, = normalized_keys([row["output_key_words"]])
        identifier = outputs.get((row["q"], key))
        require(identifier is not None and row["n"] <= 54 and row["d_Z"] >= 3, "Invalid terminal finite metric")
        metrics.add((identifier, row["d_Z"], row["n"], row["S"]))
    for block in low["blocks"]:
        for branch in block["branches"]:
            for row in branch["witnesses"]:
                add(row)
    for row in high["witnesses"]:
        add(row)
    for case in marked["cases"]:
        for row in case["orbits"]:
            add(row)
    for row in chain["witnesses"]:
        add(row)
    for case in q6["cases"]:
        for key in normalized_keys(case["expected"]["output_keys"]):
            metrics.add((outputs[6, key], 3, len(case["points"]), 6+case["ambient_dimension"]))
    merged = nondominated_metrics(previous_metrics | metrics)
    expected = {(r["output_id"], r["d_Z"], r["n"], r["S"]) for r in frontier["protocols"]}
    require(merged == expected, "The complete terminal finite metrics do not give the released frontier")
    for report in (pointed, large_report):
        for family in ("source_catalogue_sha256", "source_shard_bindings"):
            for name, digest in report.get(family, {}).items():
                evidence.digest(name, digest)
    for name, digest in marked_report["supporting_data_sha256"].items():
        evidence.digest(name, digest)
    return dict(schema="composed-protocol-coverage-v1", status="pass", maximum_protocol_length=54,
        minimum_exact_distance=3, matrix_scope="full_projective", output_equivalence="CNOT+S",
        logical_dimension_scope="all intrinsic dimensions", pareto_objectives=["n", "S"], exact_distance_is_part_of_key=True,
        preceding_frontier_counts={**previous["preceding_frontier_counts"], "52": previous["cumulative_frontier_points"]},
        length54_source_spaces=partition["source_spaces"], length54_source_blocks=partition["source_blocks"],
        lower_dimensional_pointings=low_report["totals"]["pointed_supports"],
        small_quotient_blocks=len(low["blocks"]), refined_blocks=refined["refined_blocks"],
        refined_positive_supports=refined["refined_positive_supports"], small_quotient_q6_cases=len(q6["cases"]),
        larger_quotient_pointings=large_report["selected_pointings"], larger_quotient_classes=large_report["support_classes"],
        primitive_target_classes=primitive_report["cases"], finite_marked_classes=marked_report["cases"],
        aggregate_chain_q5_classes=chain_report["q5_input_classes"], aggregate_chain_q6_classes=chain_report["q6_input_classes"],
        zero_column_dominated_large_quotient_classes=chain_report["zero_column_dominated_classes"],
        finite_candidate_metrics=len(metrics), cumulative_frontier_points=len(merged), output_classes=len({r[0] for r in merged}),
        cumulative_frontier_equals_release=True, all_finite_source_and_logical_sectors_composed=True,
        every_successor_obligation_closed=True, all_logical_dimensions_covered=True,
        all_refined_source_block_positions_checked=True, all_large_quotient_input_equivalences_recomputed=True,
        all_primitive_target_input_equivalences_recomputed=True,
        zero_column_extensions_excluded_at_fixed_exact_distance=True,
        geometry_and_sector_consistency_checks_recomputed=True, all_logical_censuses_freshly_recomputed=False,
        conditional_premises=["Completeness of the affine unital-support catalogues through length 54",
            "Completeness of each supplied finite raw or marked logical census on its explicit input domain",
            "Completeness of the aggregate chain censuses, including their complete successor selection"],
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
    result = verify(root, args.native.resolve(), args.work_directory.resolve(), args.workers)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("dependencies", "frontier_metrics")}, indent=2))


if __name__ == "__main__":
    main()
