"""Compose the length-54 dimension-eight contraction domain and finite censuses."""

import argparse
import json
from pathlib import Path
from time import perf_counter

from prepare_space_contractions import prepare
from verify_alternative_contractions import CORES
from verify_protocol_cover import Evidence, require


def check_primary_partition(domain, census, alternative):
    records = {r["id"]: r for r in domain["records"]}
    require(len(records) == domain["cores"] == 321, "Repeated or absent core")
    excluded = {r["id"] for r in records.values() if (r["core_length"], 7, r["core_index"]) in CORES}
    primary = {r["core_id"] for r in census["primary_cores"]}
    require(len(primary) == len(census["primary_cores"]) == 317 and len(excluded) == 4
            and primary.isdisjoint(excluded) and primary | excluded == set(records),
            "The primary and alternative families do not partition the complete domain")
    for row in census["primary_cores"]:
        source = records[row["core_id"]]
        require(row["core_length"] == source["core_length"] and row["core_index"] == source["core_index"]
                and row["input_sha256"] == source["sha256"] and row["raw_marked_pairs"] == source["raw_marked_pairs"]
                and row["complete_marked_orbit_census"] is True and row["marked_representatives"] >= 0
                and (row["marked_representatives"] == 0) == (row["raw_marked_pairs"] == 0),
                "A primary census has the wrong input or is incomplete")
    require(sum(r["marked_representatives"] for r in census["primary_cores"]) == census["primary_marked_representatives"]
            and sum(r["raw_marked_pairs"] for r in census["primary_cores"]) == census["primary_raw_marked_pairs"],
            "The primary census totals do not close")
    totals = alternative["totals"]["54"]
    require(totals["cases"] == len(excluded)
            and totals["raw_marked_pairs"] == sum(records[key]["raw_marked_pairs"] for key in excluded)
            and census["primary_raw_marked_pairs"] + totals["raw_marked_pairs"] == domain["raw_marked_pairs"],
            "The alternative census does not cover the complementary input mass")
    return dict(primary_cores=len(primary), alternative_cores=len(excluded), all_cores=len(records),
                complete_raw_marked_pairs=domain["raw_marked_pairs"],
                primary_marked_representatives=census["primary_marked_representatives"],
                alternative_marked_representatives=totals["marked_representatives"])


def check_quotient(census):
    candidate = bucket = representative = equivalent = inequivalent = 0
    for row in census["exact_quotient_slices"]:
        require((row["first_candidate"], row["first_bucket"], row["first_representative"])
                == (candidate, bucket, representative), "The exact quotient has a gap or overlap")
        require(row["candidate_count"] >= row["affine_class_count"] >= row["bucket_count"] > 0
                and row["every_affine_witness_independently_verified"] is True
                and row["exact_transporter_tests_complete"] is True
                and row["equivalent_comparisons"] == row["candidate_count"]-row["affine_class_count"]
                and row["inequivalent_comparisons"] >= 0, "Incomplete exact affine quotient slice")
        candidate += row["candidate_count"]
        bucket += row["bucket_count"]
        representative += row["affine_class_count"]
        equivalent += row["equivalent_comparisons"]
        inequivalent += row["inequivalent_comparisons"]
    require(candidate == census["primary_marked_representatives"] == census["independently_verified_affine_witnesses"]
            and bucket == census["signature_bucket_count"] and representative == census["affine_class_count"]
            and equivalent == census["equivalent_comparisons"] and inequivalent == census["inequivalent_comparisons"],
            "The quotient totals disagree with its candidate domain")
    for distribution, target in ((census["signature_bucket_size_distribution"], candidate),
                                 (census["affine_classes_per_signature"], representative)):
        require(all(int(size) > 0 and type(count) is int and count >= 0 for size, count in distribution.items())
                and sum(distribution.values()) == bucket
                and sum(int(size)*count for size, count in distribution.items()) == target,
                "The signature or affine-class distribution does not close")
    require(census["literal_duplicate_candidates"] == 0, "Duplicate candidates need separate accounting")
    validity = census["validity"]
    require(validity["representative_count"] == representative and validity["source_member_count"] == candidate
            and all(validity[name] == 0 for name in ("invalid_weight_count", "invalid_rank_count", "invalid_moment_count", "zero_member_count")),
            "The representative validity census fails")
    return dict(candidates=candidate, signature_buckets=bucket, affine_classes=representative,
                quotient_slices=len(census["exact_quotient_slices"]),
                equivalent_comparisons=equivalent, inequivalent_comparisons=inequivalent)


def verify(root, work):
    started = perf_counter()
    evidence = Evidence(root)
    work = Path(work).resolve()
    require(not work.exists(), "Use a new input-reconstruction directory")
    census = evidence.read("data/space_contractions/c54_m08/primary_census.json")
    require(census["schema"] == "finite-space-contraction-census-v1" and (census["c"], census["m"]) == (54, 8), "Wrong space census")
    domain = evidence.read(census["input_domain"], census["input_domain_sha256"])
    require(prepare(evidence.root, work) == domain, "The complete contraction domain does not reproduce")
    evidence.digest("certificates/low_dimensional_spaces.json", domain["affine_census_sha256"])
    for name, value in domain["catalogue_sha256"].items():
        evidence.digest(name, value)
    for row in domain["records"]:
        evidence.digest("data/space_contractions/c54_m08/"+row["task"], row["sha256"])
    for name, value in census["native_sources_sha256"].items():
        evidence.digest(name, value)
    alternative = evidence.certificate("certificates/alternative_contractions.json",
        "alternative-contraction-cover-verification-v1",
        ("marked_censuses_freshly_generated", "every_direction_freshly_replayed",
         "all_candidates_valid_and_contract_to_the_declared_source"), dict(unresolved_cases=0))
    partition = check_primary_partition(domain, census, alternative)
    quotient = check_quotient(census)
    encoding = evidence.read(census["representatives_encoding"], census["representatives_encoding_sha256"])
    require(encoding["status"] == "pass" and (encoding["c"], encoding["m"], encoding["first_index"])
            == (54, 8, 0) and encoding["count"] == quotient["affine_classes"]
            and encoding["ordered_support_sha256"] == census["ordered_representative_support_sha256"]
            and all(encoding["checks"].values()), "The representative census uses a different compact catalogue")
    index = 0
    for row in encoding["shards"]:
        require(row["first_index"] == index and row["count"] > 0, "The representative shards have a gap")
        evidence.digest(row["path"], row["sha256"])
        index += row["count"]
    require(index == quotient["affine_classes"], "The representative domain is incomplete")
    return dict(schema="length54-m08-space-cover-composition-v1", status="pass",
        input_partition=partition, exact_quotient=quotient,
        every_input_reconstructed_from_the_complete_affine_census=True,
        primary_and_alternative_core_families_exhaust_the_domain=True,
        primary_aggregate_census_accounting_closes=True, primary_affine_quotient_freshly_recomputed=False,
        conditional_premises=["The primary finite marked-orbit census is complete.",
                              "The primary finite exact affine quotient census is complete."],
        is_global_completeness_certificate=False, dependencies=evidence.bindings,
        elapsed_seconds=round(perf_counter()-started, 6))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.work_directory)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "dependencies"}, indent=2))


if __name__ == "__main__":
    main()
