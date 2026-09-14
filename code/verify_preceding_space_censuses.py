"""Compose the length-50 and length-52 finite space censuses."""

import argparse
import json
from math import comb
from pathlib import Path

from prepare_space_contractions import task_text
from verify_protocol_cover import Evidence, require
from verify_zero_fibre_spaces import sector


def distribution(values, groups, mass):
    require(all(str(k).isdigit() and int(k) > 0 and type(v) is int and v >= 0
                for k, v in values.items())
            and sum(values.values()) == groups
            and sum(int(k)*v for k, v in values.items()) == mass,
            "A finite partition has the wrong size or mass")


def check_quotient(data):
    signature, quotient = data["signature_census"], data["exact_quotient"]
    count, classes = data["candidate_count"], data["class_count"]
    buckets = signature["signature_bucket_count"]
    negative = quotient["inequivalent_comparison_count"]
    require(type(count) is int and type(classes) is int and count >= classes >= buckets >= 0
            and signature["candidate_count"] == quotient["candidate_count"] == count
            and quotient["signature_bucket_count"] == buckets
            and quotient["affine_class_count"] == data["output_domain"]["classes"] == classes
            and quotient["equivalent_comparison_count"] == count-classes
            and type(negative) is int and negative >= 0
            and quotient["exact_transporter_comparison_count"] == count-classes+negative,
            "The exact quotient accounting does not close")
    distribution(signature["bucket_size_distribution"], buckets, count)
    if "member_count_distribution" in quotient:
        distribution(quotient["member_count_distribution"], classes, count)
    if "bucket_affine_class_count_distribution" in quotient:
        distribution(quotient["bucket_affine_class_count_distribution"], buckets, classes)
    intervals = data["exact_quotient_intervals"]
    if intervals:
        require([row["interval"] for row in intervals] == list(range(len(intervals)))
                and sum(row["candidate_count"] for row in intervals) == count
                and sum(row["affine_class_count"] for row in intervals) == classes,
                "The quotient intervals omit or repeat candidates")
        position = 0
        for row in intervals:
            require(row["candidate_count"] >= row["affine_class_count"] > 0,
                    "Invalid quotient interval counts")
            if "first_bucket" in row:
                require(row["first_bucket"] == position and row["last_bucket"] > position,
                        "The quotient bucket intervals have a gap or overlap")
                position = row["last_bucket"]
            else:
                require(0 < row["bucket_count"] <= row["affine_class_count"], "Invalid bucket count")
                position += row["bucket_count"]
            if "equivalent_comparison_count" in row:
                require(row["equivalent_comparison_count"] == row["candidate_count"]-row["affine_class_count"],
                        "The interval has inconsistent equivalent comparisons")
        require(position == buckets
                and sum(row.get("inequivalent_comparison_count", 0) for row in intervals) == negative,
                "The quotient intervals do not cover all buckets or comparisons")
    return dict(candidates=count, affine_classes=classes, signature_buckets=buckets,
                negative_comparisons=negative, quotient_intervals=len(intervals))


def check_generation(data):
    m, generation = data["m"], data["generation"]
    count = data["candidate_count"]
    sources = {str(row["multiplicity"]): row["classes"] for row in data["source_domains"]}
    if m == 8:
        require(generation["direction_marked_orbit_count"] == count,
                "The marked enumeration and affine quotient disagree")
    elif m == 9:
        full = generation["full_rank_marked_extension_count"]
        require(generation["raw_marked_extension_count"]-full == sources["0"]
                and generation["batch_distinct_candidate_record_count"] == count
                and generation["within_batch_duplicate_count"] == 0,
                "The complete lift and quotient domains disagree")
        if data["minimum_direction_filter"]:
            require(generation["source_counts_by_multiplicity"] == sources
                    and generation["source_count"] == sum(sources.values())
                    and sum(generation["minimum_multiplicity_counts"].values()) == full
                    and 0 <= generation["retained_occurrence_count"] == count <= full,
                    "The minimum-direction census is incomplete")
        else:
            require(count == full, "Unaccounted removal of unfiltered lifts")
    elif m == 10:
        require(data["minimum_direction_filter"] is False
                and generation["full_rank_occurrence_count"] == count
                and generation["within_batch_distinct_candidate_count"] == count
                and generation["within_batch_literal_duplicate_count"] == 0,
                "The one-fibre generation is incomplete")
    else:
        require(list(sources) == ["0"] and sources["0"] == data["zero_fibre_source_count"]
                and generation["candidate_count"] == count,
                "The zero-fibre source domain is incomplete")
        if sources["0"] == 0:
            require(count == 0, "An empty predecessor has nonempty graph extensions")
    if m in (9, 10):
        witness = data["independent_positive_replay"]
        signature, quotient = data["signature_census"], data["exact_quotient"]
        require(witness["verified_witness_count"] == witness["candidate_count"] == count
                and witness["class_count"] == data["class_count"]
                and witness["bucket_count"] == signature["signature_bucket_count"]
                and witness["negative_record_count"] == quotient["inequivalent_comparison_count"]
                and witness["ledger_sha256"] == signature["ledger_sha256"]
                and witness["bucket_index_sha256"] == signature["bucket_index_sha256"],
                "The retained affine-witness census uses a different domain")
        if quotient["inequivalent_comparison_count"]:
            negative = data["independent_negative_replay"]
            require(negative["negative_comparison_count"] == quotient["inequivalent_comparison_count"]
                    and negative["disagreement_count"] == 0, "Incomplete independent negative replay")


def bind_sector(evidence, index, domain):
    rows = [s for length in index["lengths"] if length["c"] == domain["c"]
            for s in length["sectors"] if s["m"] == domain["m"]]
    if not rows:
        require(domain["classes"] == 0 and domain["encoding_certificates"] == [],
                "An undeclared sector has nonempty data")
        return
    require(len(rows) == 1, "Repeated catalogue sector")
    row, = rows
    require(row["data_complete"] and row["included_classes"] == row["expected_classes"] == domain["classes"]
            and row["encoding_certificates"] == domain["encoding_certificates"],
            "The finite census is bound to a different catalogue")
    for item in domain["encoding_certificates"]:
        evidence.digest(item["path"], item["sha256"])


def check_core_partition(data, profiles, alternative):
    primary = {(row["c"], row["m"], row["index"]): row for row in data["primary_cores"]}
    other = {(row["c"], row["m"], row["index"]): row for row in data["alternative_cores"]}
    require(len(primary) == len(data["primary_cores"]) and len(other) == len(data["alternative_cores"])
            and primary.keys().isdisjoint(other) and primary.keys() | other.keys() == profiles.keys(),
            "The marked core families do not partition every contraction source")
    for key, row in primary.items():
        require(row["profile"] == profiles[key] and row["multiplicity"] == (data["c"]-key[0])//2
                and type(row["marked_orbits"]) is int and 0 <= row["marked_orbits"] <= profiles[key]["raw_marked_pairs"]
                and (row["marked_orbits"] == 0) == (profiles[key]["raw_marked_pairs"] == 0),
                "A primary orbit census uses a different source profile")
    for key, row in other.items():
        require(row["multiplicity"] == (data["c"]-key[0])//2, "Incorrect alternative multiplicity")
    require(sum(row["marked_orbits"] for row in primary.values()) == data["candidate_count"]
            and sum(profiles[key]["raw_marked_pairs"] for key in primary) == data["generation"]["raw_marked_pair_count"],
            "The primary core counts do not close")
    totals = alternative["totals"][str(data["c"])]
    require(totals["cases"] == len(other)
            and totals["raw_marked_pairs"] == sum(profiles[key]["raw_marked_pairs"] for key in other),
            "The alternative-direction cover leaves an unresolved source")
    return dict(primary_cores=len(primary), alternative_cores=len(other),
                raw_marked_pairs=sum(row["raw_marked_pairs"] for row in profiles.values()))


def verify(root):
    evidence = Evidence(root)
    index = evidence.read("data/spaces/index.json")
    inventory = evidence.read("data/space_contractions/preceding_sectors/index.json")
    expected = [(c, m) for c, last in ((50, 14), (52, 16)) for m in range(8, last+1)]
    require(inventory["schema"] == "preceding-space-sector-index-v1"
            and [(r["c"], r["m"]) for r in inventory["sectors"]] == expected,
            "The preceding-sector collection has a gap or duplicate")
    base = evidence.read("certificates/low_dimensional_spaces.json")
    require(base["status"] == "pass" and base["complete_and_pairwise_inequivalent"] is True
            and base["maximum_length"] == 54 and base["maximum_affine_dimension"] == 7,
            "The complete low-dimensional census is required")
    forms = {(r["c"], r["m"], r["index"]): r for r in base["canonical_maps"]}
    require(len(forms) == len(base["canonical_maps"]), "Repeated affine base record")
    for path, digest in base["catalogue_sha256"].items():
        evidence.digest(path, digest)
    alternative = evidence.certificate("certificates/alternative_contractions.json",
        "alternative-contraction-cover-verification-v1",
        ("marked_censuses_freshly_generated", "every_direction_freshly_replayed",
         "all_candidates_valid_and_contract_to_the_declared_source"), dict(unresolved_cases=0))
    results = []
    for item in inventory["sectors"]:
        data = evidence.read(item["path"], item["sha256"])
        c, m = item["c"], item["m"]
        bound = comb(c, 2)//((1 << m)-1)
        require(data["schema"] == "preceding-space-sector-census-v1" and (data["c"], data["m"]) == (c, m)
                and [(r["multiplicity"], r["c"], r["m"]) for r in data["source_domains"]]
                == [(n, c-2*n, m-1) for n in range(bound+1)]
                and item["candidate_count"] == data["candidate_count"]
                and item["class_count"] == data["class_count"], "The finite domain has a missing contraction route")
        for domain in [*data["source_domains"], data["output_domain"]]:
            bind_sector(evidence, index, domain)
        quotient = check_quotient(data)
        check_generation(data)
        result = dict(c=c, m=m, **quotient)
        if m == 8:
            profiles = {}
            for domain in data["source_domains"]:
                weight, n = domain["c"], domain["multiplicity"]
                if n > 3:
                    require(not any(k[0] == weight and k[1] < 7 for k in forms),
                            "A proper-span contraction core must also be included")
                points = sector(evidence, index, weight, 7)
                require(len(points) == domain["classes"], "The base domain has a different class count")
                for i, support in enumerate(points):
                    key = weight, 7, i
                    _, profiles[key] = task_text(f"c{weight}_m07_{i:09d}", support, forms[key], c)
            certificate = data["alternative_cover_certificate"]
            evidence.digest(certificate["path"], certificate["sha256"])
            result["marked_core_cover"] = check_core_partition(data, profiles, alternative)
        results.append(result)
    return dict(schema="preceding-space-census-composition-v1", status="pass", sectors=results,
                all_source_domains_bound_to_compact_catalogues=True,
                every_marked_core_profile_independently_recomputed=True,
                primary_and_alternative_core_domains_exhaustive=True,
                finite_domain_and_aggregate_accounting_verified=True,
                full_affine_quotients_freshly_recomputed=False,
                individual_assignment_streams_included=False,
                conditional_premises=["The predecessor space catalogues are complete.",
                    "The retained finite generation and exact affine quotient censuses are complete."],
                is_global_completeness_certificate=False, dependencies=evidence.bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "dependencies"}, indent=2))


if __name__ == "__main__":
    main()
