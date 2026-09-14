"""Compose the complete finite space induction, with explicit computational premises."""

import argparse
from collections import Counter
import json
from math import comb
from pathlib import Path

from verify_protocol_cover import Evidence, require


def dimension_bound(length):
    require(type(length) is int and length >= 16, "Invalid positive unital length")
    return (length+length//16)//3-1


def compose_sectors(counts, base_counts, routes):
    covered, records = {}, []
    for c in range(16, 55, 2):
        for m in range(4, dimension_bound(c)+1):
            key, count = (c, m), counts.get((c, m), 0)
            if m <= 7:
                require(count == base_counts.get(key, 0), "The complete affine base and catalogue disagree")
                reason = "complete_low_dimensional_affine_base"
            elif key in routes:
                require(routes[key]["classes"] == count, "A finite sector has a different output count")
                bound = comb(c, 2)//((1 << m)-1)
                require(2*bound < c, "The positive-core recurrence does not cover this sector")
                for n in range(bound+1):
                    source_length = c-2*n
                    for source_m in range(4, min(m-1, dimension_bound(source_length))+1):
                        require((source_length, source_m) in covered,
                                "A contraction source is not closed earlier in the induction")
                reason = routes[key]["proof"]
            else:
                require(comb(c, 2) < (1 << m)-1 and covered.get((c, m-1)) == 0 and count == 0,
                        "An unclosed finite space sector remains")
                reason = "empty_zero_fibre_successor"
            covered[key] = count
            records.append(dict(c=c, m=m, classes=count, proof=reason))
    require(all(count == 0 or key in covered for key, count in counts.items()),
            "A nonempty catalogue sector lies outside the proved dimension range")
    require(set(routes) <= set(covered), "A supplied finite census lies outside the induction")
    return records


def verify(root):
    evidence = Evidence(root)
    index = evidence.read("data/spaces/index.json")
    counts = {(length["c"], row["m"]): row["included_classes"]
              for length in index["lengths"] for row in length["sectors"]}
    require(len(counts) == sum(len(length["sectors"]) for length in index["lengths"])
            and all(row["data_complete"] and row["included_classes"] == row["expected_classes"]
                    for length in index["lengths"] for row in length["sectors"]), "Incomplete or duplicate space data")
    readback = evidence.certificate("certificates/all_spaces.json", "triorthogonal-included-space-verification-v1",
        ("validity_checked_for_all_included_representatives", "all_declared_catalogue_records_verified"),
        dict(catalogue_index_sha256=evidence.digest("data/spaces/index.json"), representatives=sum(counts.values())))
    base = evidence.read("certificates/low_dimensional_spaces.json")
    require(base["status"] == "pass" and base["complete_and_pairwise_inequivalent"] is True
            and base["maximum_length"] == 54 and base["maximum_affine_dimension"] == 7,
            "Incomplete low-dimensional orbit-mass census")
    for path, digest in base["catalogue_sha256"].items():
        evidence.digest(path, digest)
    base_counts = Counter((r["c"], r["m"]) for r in base["canonical_maps"])
    require(len({(r["c"], r["m"], r["index"]) for r in base["canonical_maps"]}) == sum(base_counts.values()),
            "Repeated low-dimensional affine class")
    finite = evidence.certificate("certificates/finite_space_censuses.json",
        "finite-affine-contraction-collection-verification-v1",
        ("all_catalogue_records_matched_by_affine_frames", "all_witness_streams_losslessly_verified",
         "all_candidates_and_positive_maps_freshly_replayed", "all_class_canonical_forms_freshly_recomputed"))
    routes = {(r["c"], r["m"]): dict(classes=r["class_count"], proof="fresh_finite_contraction_census")
              for r in finite["reports"]}
    require(len(routes) == finite["sectors"], "Repeated finite contraction sector")
    marked48 = evidence.certificate("certificates/c48_marked_space_census.json", "marked-c48-space-census-verification-v1",
        ("complete_core_domain_reconstructed", "proper_span_core_included", "all_marked_outputs_and_contractions_checked",
         "marked_enumeration_freshly_recomputed", "all_positive_affine_witnesses_freshly_replayed",
         "all_class_canonical_forms_freshly_recomputed", "complete_given_low_dimensional_base"), dict(c=48, m=8))
    routes[48, 8] = dict(classes=marked48["affine_classes"], proof="fresh_marked_contraction_census")
    preceding = evidence.certificate("certificates/preceding_space_censuses.json", "preceding-space-census-composition-v1",
        ("all_source_domains_bound_to_compact_catalogues", "every_marked_core_profile_independently_recomputed",
         "primary_and_alternative_core_domains_exhaustive", "finite_domain_and_aggregate_accounting_verified"))
    computational = []
    for row in preceding["sectors"]:
        key = row["c"], row["m"]
        require(key not in routes, "Repeated preceding finite sector")
        routes[key] = dict(classes=row["affine_classes"], proof="retained_finite_contraction_census")
        computational.append(dict(c=key[0], m=key[1],
                                  census=f"data/space_contractions/preceding_sectors/c{key[0]}_m{key[1]:02d}.json"))
    marked54 = evidence.certificate("certificates/c54_m08_space_cover.json", "length54-m08-space-cover-composition-v1",
        ("every_input_reconstructed_from_the_complete_affine_census",
         "primary_and_alternative_core_families_exhaust_the_domain", "primary_aggregate_census_accounting_closes"))
    routes[54, 8] = dict(classes=marked54["exact_quotient"]["affine_classes"], proof="retained_marked_contraction_census")
    computational.append(dict(c=54, m=8, census="data/space_contractions/c54_m08/primary_census.json"))
    middle = evidence.certificate("certificates/c54_middle_space_censuses.json", "middle-space-sector-census-composition-v1",
        ("finite_domain_and_aggregate_accounting_verified",), dict(c=54, dimensions=[9, 10, 11, 12]))
    for m in (9, 10, 11, 12):
        routes[54, m] = dict(classes=middle["sectors"][str(m)]["affine_classes"], proof="retained_finite_contraction_census")
        computational.append(dict(c=54, m=m, census=f"data/space_contractions/middle_sectors/c54_m{m:02d}.json"))
    for m in range(13, 18):
        name = f"certificates/zero_fibre_c54_m{m:02d}.json"
        if m <= 14:
            row = evidence.certificate(f"certificates/zero_fibre_c54_m{m:02d}_witness_replay.json",
                "separated-graph-witness-replay-v1",
                ("every_affine_shear_coset_independently_enumerated", "every_positive_affine_map_replayed",
                 "targets_separated_by_independent_exact_invariants", "complete_given_predecessor_catalogue"), dict(c=54, m=m))
            evidence.digest(name, row["certificate_sha256"])
        else:
            row = evidence.certificate(name, "complete-zero-fibre-space-recurrence-v1",
                ("complete_zero_fibre_cover_by_pair_averaging", "all_graph_lifts_freshly_enumerated",
                 "every_affine_map_replayed", "every_target_reached", "targets_pairwise_affine_inequivalent",
                 "complete_given_predecessor_catalogue"), dict(c=54, m=m))
        require(row["source_classes"] == counts[54, m-1], "A graph proof uses a different predecessor")
        routes[54, m] = dict(classes=row["target_classes"], proof="fresh_complete_graph_cover")
    records = compose_sectors(counts, base_counts, routes)
    totals = {str(c): sum(r["classes"] for r in records if r["c"] == c) for c in range(16, 55, 2)}
    require(sum(totals.values()) == readback["representatives"], "The induction does not cover every distributed record")
    return dict(schema="space-classification-induction-composition-v1", status="pass", maximum_length=54,
        representatives=sum(totals.values()), classes_by_length=totals, sectors=records,
        all_dimensions_covered_by_independent_bound=True, all_contraction_predecessors_closed_in_inductive_order=True,
        all_distributed_representatives_validated=True, fresh_finite_induction_complete_through_length=48,
        finite_censuses_for_independent_regeneration=computational,
        conditional_premises=["The listed retained finite marked-lift and exact affine quotient censuses are complete."],
        is_unconditional_machine_checked_completeness_certificate=False, dependencies=evidence.bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("sectors", "dependencies")}, indent=2))


if __name__ == "__main__":
    main()
