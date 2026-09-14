"""Check finite space-sector composition without conflating it with a full rerun."""

import argparse
from collections import Counter
import json
from math import comb
from pathlib import Path
import re

from verify_protocol_cover import Evidence, require


def check9(data):
    routes = {str(row["multiplicity"]): row["classes"] for row in data["source_domains"]}
    source_counts = Counter()
    parts = data["generation_parts"]
    require(len(parts) == 3 and data["minimum_direction_filter"] is True, "Incomplete minimum-direction census")
    for part in parts:
        counts = part["source_counts_by_multiplicity"]
        require(sum(counts.values()) == part["source_count"]
                and part["raw_marked_extension_count"]-part["full_rank_marked_extension_count"] == counts.get("0", 0)
                and 0 <= part["retained_occurrence_count"] <= part["full_rank_marked_extension_count"]
                and part["within_batch_duplicate_count"] == 0
                and part["batch_distinct_candidate_record_count"] == part["retained_occurrence_count"],
                "Inconsistent source or minimum-direction accounting")
        if "minimum_multiplicity_counts" in part:
            require(sum(part["minimum_multiplicity_counts"].values()) == part["full_rank_marked_extension_count"],
                    "The minimum multiplicities do not partition the full lift domain")
        source_counts.update(counts)
    require(dict(source_counts) == routes and parts[2]["source_counts_by_multiplicity"] == {"2": routes["2"]}
            and parts[2]["retained_occurrence_count"] == 0, "The two-fibre route is not closed")
    local = data["local_exact_quotients"]
    require(len(local) == 2, "The local exact quotients do not cover both nonempty candidate parts")
    for part, row in zip(parts, local):
        require(row["source_candidate_count"] == part["retained_occurrence_count"]
                and row["singleton_class_count"]+row["repeated_exact_class_count"] == row["representative_count"]
                and row["negative_comparisons_independently_replayed"] is True,
                "Incomplete local exact quotient")
        for witness in (row["producer"], row["verifier"]):
            require(witness["source_member_count"] == row["source_candidate_count"]
                    and witness["affine_class_count"] == row["representative_count"]
                    and witness["signature_bucket_count"] == row["signature_bucket_count"], "The local quotient verifiers disagree")
    candidates = sum(row["source_candidate_count"] for row in local)
    representatives = sum(row["representative_count"] for row in local)
    for merge in (data["weighted_merge"], data["weighted_merge_verification"]):
        require(merge["input_record_count"] == merge["unique_representative_count"] == representatives
                and merge["input_member_count"] == merge["output_member_count"] == candidates
                and merge["literal_duplicate_count"] == 0
                and merge["singleton_signature_bucket_count"]+merge["repeated_representative_count"] == representatives,
                "The weighted merge loses candidates or representatives")
    final = data["global_exact_quotient"]
    require(final["source_candidate_count"] == candidates and final["source_local_representative_count"] == representatives
            and final["singleton_class_count"]+final["repeated_candidate_count"] == representatives
            and final["singleton_class_count"]+final["repeated_exact_class_count"] == final["representative_count"]
            and final["singleton_class_count"]+final["repeated_signature_bucket_count"] == final["signature_bucket_count"]
            and final["inequivalent_comparison_count"] == final["independently_replayed_negative_comparison_count"],
            "The final weighted exact quotient is inconsistent")
    require(final["producer"]["output_member_count"] == final["verifier"]["verified_member_count"] == candidates
            and final["producer"]["assignment_record_count"] == final["verifier"]["verified_assignment_count"] == final["repeated_candidate_count"]
            and final["producer"]["representative_count"] == final["verifier"]["representative_count"] == final["representative_count"],
            "The final weighted witness verification is incomplete")
    require(final["representative_count"] == data["output_domain"]["classes"], "Wrong final class count")
    return dict(full_rank_lifts_before_filter=sum(p["full_rank_marked_extension_count"] for p in parts),
                retained_candidates=candidates, local_representatives=representatives,
                affine_classes=final["representative_count"], minimum_two_candidates=0)


def check10(data):
    routes = {str(row["multiplicity"]): row["classes"] for row in data["source_domains"]}
    positions = Counter()
    for row in data["source_profile_intervals"]:
        key = str(row["multiplicity"])
        require(key in routes and row["first_source_index"] == positions[key] and row["source_count"] > 0,
                "The source profile has a gap or overlap")
        positions[key] += row["source_count"]
    require(dict(positions) == routes and set(data["source_profiles"]) == set(routes), "Incomplete source profile domain")
    total = 0
    for key, row in data["source_profiles"].items():
        multiplicity = int(key)
        require(row["source_count"] == routes[key] and row["core_length"] == 54-2*multiplicity,
                "The source profile uses a different parent sector")
        histogram = Counter()
        for profile, count in row["quadratic_profile_distribution"].items():
            match = re.fullmatch(r"rank(\d+)_lift(\d+)", profile)
            require(match is not None and type(count) is int and count >= 0, "Invalid rank profile")
            rank, lift = map(int, match.groups())
            require(rank+lift+10 == row["core_length"], "A rank-nullity profile fails")
            histogram[lift] += count
        require(sum(histogram.values()) == row["source_count"], "The rank distribution is incomplete")
        if multiplicity == 0:
            require(sum(((1 << lift)-1)*count for lift, count in histogram.items()) == row["full_rank_marked_extension_count"]
                    and row["contributing_source_count"] == row["source_count"]-histogram[0], "The graph-lift census disagrees with rank-nullity")
        else:
            require(row["contributing_source_count"]+row["zero_compatible_source_count"] == row["source_count"]
                    and row["full_rank_marked_extension_count"] >= row["compatible_fibre_count"], "The one-fibre profile is inconsistent")
        total += row["full_rank_marked_extension_count"]
    quotient = data["exact_quotient"]
    require(total == data["candidate_count"] == quotient["repeated_candidate_count"]
            and quotient["singleton_signature_bucket_count"] == 0
            and quotient["equivalent_comparison_count"] == total-quotient["affine_class_count"]
            and quotient["comparison_count"] == quotient["equivalent_comparison_count"]+quotient["inequivalent_comparison_count"]
            and quotient["inequivalent_comparison_count"] == quotient["independent_negative_replay_count"]
            and quotient["affine_class_count"] == quotient["final_affine_class_count"] == data["output_domain"]["classes"],
            "The complete exact quotient accounting fails")
    weighted = data["weighted_representatives"]
    require(weighted["count"] == quotient["affine_class_count"] and weighted["source_distinct_candidates"] == total,
            "The final representative ledger is inconsistent")
    for row in (weighted["producer"], weighted["verifier"]):
        require(row["source_member_count"] == total and row["affine_class_count"] == weighted["count"], "The representative verification is incomplete")
    return dict(full_rank_lifts=total, affine_classes=weighted["count"], source_profile_intervals=len(data["source_profile_intervals"]))


def check_upper(data, profile):
    source_count = data["source_domains"][0]["classes"]
    histogram = profile["lift_dimension_histogram"]
    require(profile["source_interval"] == [0, source_count] and profile["source_count"] == source_count
            and sum(histogram.values()) == source_count
            and sum(((1 << int(lift))-1)*count for lift, count in histogram.items()) == profile["nonaffine_graph_lifts"]
            and source_count-histogram.get("0", 0) == data["contributing_sources"], "The fresh source profile does not close the full domain")
    generation, quotient = data["generation"], data["exact_quotient"]
    count = profile["nonaffine_graph_lifts"]
    require(count == generation["candidate_occurrence_count"] == generation["input_distinct_candidate_count"]
            == generation["candidate_count"] == quotient["candidate_count"]
            and generation["within_batch_literal_duplicate_count"] == generation["global_literal_duplicate_count"] == 0,
            "The graph-lift census and complete source profiles disagree")
    first = classes = negatives = 0
    for row in data["exact_quotient_intervals"]:
        require(row["first_candidate"] == first and row["candidate_count"] >= row["affine_class_count"] > 0,
                "The exact quotient intervals have a gap or overlap")
        first += row["candidate_count"]
        classes += row["affine_class_count"]
        negatives += row["negative_comparison_count"]
    require(first == count and classes == quotient["affine_class_count"] == data["output_domain"]["classes"]
            and negatives == quotient["negative_comparison_count"] == quotient["independently_replayed_negative_comparisons"]
            and quotient["negative_replay_disagreements"] == 0, "The exact affine quotient is incomplete")
    for distribution, groups in ((generation["bucket_size_distribution"], quotient["signature_bucket_count"]),
                                  (data["representative_member_count_distribution"], classes)):
        require(sum(distribution.values()) == groups and sum(int(size)*number for size, number in distribution.items()) == count,
                "The candidate partition does not close")
    return dict(source_classes=source_count, contributing_sources=data["contributing_sources"],
                nonaffine_graph_lifts=count, affine_classes=classes, quotient_intervals=len(data["exact_quotient_intervals"]))


def verify(root, dimensions):
    evidence = Evidence(root)
    index = evidence.read("data/spaces/index.json")
    records = {}
    require(set(dimensions) <= {9, 10, 11, 12} and len(set(dimensions)) == len(dimensions), "Unsupported or repeated dimensions")
    for m in dimensions:
        data = evidence.read(f"data/space_contractions/middle_sectors/c54_m{m:02d}.json")
        require(data["schema"] == "finite-space-sector-census-v1" and (data["c"], data["m"]) == (54, m), "Different finite sector")
        bound = comb(54, 2)//((1 << m)-1)
        require([(row["multiplicity"], row["c"], row["m"]) for row in data["source_domains"]]
                == [(n, 54-2*n, m-1) for n in range(bound+1)], "The complete contraction cover is not present")
        for domain in [*data["source_domains"], data["output_domain"]]:
            row, = [s for length in index["lengths"] if length["c"] == domain["c"] for s in length["sectors"] if s["m"] == domain["m"]]
            require(row["data_complete"] and row["included_classes"] == row["expected_classes"] == domain["classes"]
                    and row["encoding_certificates"] == domain["encoding_certificates"], "A census uses a different catalogue sector")
            for item in domain["encoding_certificates"]:
                evidence.digest(item["path"], item["sha256"])
        for name, value in data["native_sources_sha256"].items():
            evidence.digest(name, value)
        if m in (9, 10):
            records[str(m)] = check9(data) if m == 9 else check10(data)
        else:
            profile = evidence.certificate(f"certificates/zero_fibre_native_c54_m{m:02d}_input.json", "native-zero-fibre-input-v1",
                ("complete_parent_sector_included", "bounded_memory_catalogue_stream"), dict(c=54, m=m))
            records[str(m)] = check_upper(data, profile)
    return dict(schema="middle-space-sector-census-composition-v1", status="pass", c=54, dimensions=list(dimensions), sectors=records,
                finite_domain_and_aggregate_accounting_verified=True, full_affine_quotients_freshly_recomputed=False,
                conditional_premises=["Every parent space catalogue is complete.",
                                      "The retained finite lift and exact affine quotient censuses are complete."],
                is_global_completeness_certificate=False, dependencies=evidence.bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--dimensions", nargs="+", type=int, default=[9, 10, 11, 12])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.dimensions)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "dependencies"}, indent=2))


if __name__ == "__main__":
    main()
