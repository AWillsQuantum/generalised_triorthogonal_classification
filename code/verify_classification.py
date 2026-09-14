"""Join the space induction, complete protocol domains and exact frontier."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

from verify_high_dimensional_sector import check_intervals
from verify_protocol_cover import Evidence, require
from verify_protocols import verify as verify_witnesses
from verify_space_classification_cover import verify as verify_spaces


def check_source_partition(index, domains):
    expected = {(row["c"], sector["m"]): sector["included_classes"]
                for row in index["lengths"] for sector in row["sectors"]}
    intervals = defaultdict(list)
    for c, m, first, end in domains:
        require((c, m) in expected, "Unknown protocol source sector")
        require(type(first) is int and type(end) is int and 0 <= first < end,
                "Invalid protocol source interval")
        intervals[c, m].append((first, end))
    for key, count in expected.items():
        check_intervals(intervals[key], count)
    return {str(c): sum(count for (length, _), count in expected.items() if length == c)
            for c in range(16, 55, 2)}


def check_frontier_scope(frontier, cover):
    require(frontier["maximum_protocol_length"] == cover["maximum_protocol_length"] == 54
            and frontier["minimum_distance"] == cover["minimum_exact_distance"] == 3,
            "Inconsistent protocol length or distance scope")
    require(frontier["equivalence"] == cover["output_equivalence"] == "CNOT+S"
            and frontier["intrinsic_outputs_only"] is True
            and cover["logical_dimension_scope"] == "all intrinsic dimensions",
            "Inconsistent output equivalence or logical scope")
    require(frontier["pareto_partition"] == ["output_id", "exact d_Z"]
            and cover["exact_distance_is_part_of_key"] is True
            and frontier["pareto_objectives"] == cover["pareto_objectives"] == ["n", "S"],
            "Inconsistent Pareto objectives")
    require(frontier["matrix_scope"] == "Full row rank with pairwise distinct complete columns; one zero column permitted"
            and cover["matrix_scope"] == "full_projective", "Inconsistent matrix scope")
    def metric(row):
        return tuple(row[k] for k in ("output_id", "d_Z", "n", "S"))
    metrics = [metric(row) for row in frontier["protocols"]]
    retained = [metric(row) for row in cover["frontier_metrics"]]
    require(len(set(metrics)) == len(metrics) and len(set(retained)) == len(retained)
            and set(metrics) == set(retained), "The witnesses and complete census have different frontiers")
    require(cover["cumulative_frontier_points"] == len(metrics)
            and cover["output_classes"] == len(frontier["outputs"]), "Incorrect frontier counts")
    for cutoff, count in cover["preceding_frontier_counts"].items():
        require(sum(row[2] <= int(cutoff) for row in metrics) == count,
                "A preceding cutoff has a different frontier")


def verify(root):
    evidence = Evidence(root)
    spaces = verify_spaces(evidence.root)
    for name, digest in spaces["dependencies"].items():
        evidence.digest(name, digest)
    index = evidence.read("data/spaces/index.json")
    frontier_name = "data/protocols/pareto_frontier.json"
    frontier = evidence.read(frontier_name)
    witnesses = verify_witnesses(evidence.path(frontier_name))
    cover = evidence.certificate("certificates/through54_protocol_cover.json", "composed-protocol-coverage-v1",
        ("cumulative_frontier_equals_release", "all_finite_source_and_logical_sectors_composed",
         "every_successor_obligation_closed", "all_logical_dimensions_covered",
         "all_refined_source_block_positions_checked", "all_large_quotient_input_equivalences_recomputed",
         "all_primitive_target_input_equivalences_recomputed", "zero_column_extensions_excluded_at_fixed_exact_distance",
         "geometry_and_sector_consistency_checks_recomputed"))
    check_frontier_scope(frontier, cover)
    domains = []
    def source(row):
        domains.append((row["c"], row["m"], row["index"], row["index"]+1))
    def block(row):
        domains.append((row["c"], row["m"], *row["source_interval"]))
    short = evidence.read("data/protocol_sectors/through48_source_domain.json")
    for row in short["sources"]:
        # The source-domain index is not the space index.
        source(row["space"])
    for row in short["separate_sources"]:
        source(row)
    for row in evidence.read("data/protocol_sectors/finite_shorter_source_blocks.json")["blocks"]:
        block(row)
    for row in evidence.read("data/protocol_sectors/length50_source_domain.json")["blocks"]:
        block(row)
    domain52 = evidence.read("data/protocol_sectors/length52_source_domain.json")
    for row in domain52["main_blocks"]:
        block(row)
    for row in domain52["separate_sources"]:
        source(row)
    for name in ("length54_low_dimensional_blocks", "length54_high_dimensional"):
        for row in evidence.read(f"data/protocol_sectors/{name}.json")["blocks"]:
            if row["c"] == 54:
                block(row)
    totals = check_source_partition(index, domains)
    require(totals == spaces["classes_by_length"] and totals["54"] == cover["length54_source_spaces"],
            "The protocol source partition and space induction disagree")
    for name in ("support_correspondence", "space_classification", "space_audit_composition",
                 "canonical_augmentation", "logical_dimension_closure", "output_profile_pruning",
                 "protocol_cover_composition", "witness_verification"):
        evidence.digest(f"theory/{name}.md")
    return dict(schema="complete-classification-composition-v1", status="pass",
        maximum_protocol_length=54, minimum_exact_distance=3, output_equivalence="CNOT+S",
        matrix_scope="full_projective", logical_dimension_scope="all intrinsic dimensions",
        pareto_objectives=["n", "S"], exact_distance_is_part_of_key=True,
        space_representatives=spaces["representatives"], space_classes_by_length=totals,
        protocol_frontier_points=witnesses["protocols"], output_classes=witnesses["outputs"],
        protocol_distance_counts=dict(sorted(Counter(row["d_Z"] for row in frontier["protocols"]).items())),
        every_protocol_parent_closed_by_space_induction=True,
        source_intervals_independently_composed=True, all_logical_successors_closed=True,
        exact_frontier_witnesses_independently_verified=True,
        retained_space_censuses=spaces["finite_censuses_for_independent_regeneration"],
        conditional_premises=spaces["conditional_premises"]+cover["conditional_premises"][1:],
        all_large_finite_censuses_freshly_recomputed=False,
        is_unconditional_machine_checked_completeness_certificate=False,
        dependencies=evidence.bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("dependencies", "retained_space_censuses")}, indent=2))


if __name__ == "__main__":
    main()
