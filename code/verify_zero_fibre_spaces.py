"""Exhaustively verify a zero-fibre space recurrence and its exact affine quotient."""

import argparse
from collections import Counter, defaultdict
import json
from math import comb
from pathlib import Path
import sys
from time import perf_counter

from selected_spaces import read_selected
from space_codec import binary_rank, validate_unital_support
from space_lifts import evaluation_rows, lift_family, reconstruct_lift
from verify_protocol_cover import Evidence, require

sys.path.insert(0, str(Path(__file__).resolve().parent / "space_algorithms"))
from exact_affine_transporter import prepare_support, transporter, witness_is_valid
from support_invariants import hyperplane_intersection_profile


def graph_labels(points, dimension):
    """Give the complete non-affine label quotient in deterministic Gray order."""
    validate_unital_support(points, dimension)
    particular, basis = lift_family(points, (), dimension)
    quadratic_rank = binary_rank(evaluation_rows(points, dimension, 2))
    require(particular == 0 and len(basis) == len(points)-quadratic_rank-dimension-1,
            "The graph quotient does not satisfy rank-nullity")
    return quadratic_rank, basis


def invariant(prepared):
    return tuple(sorted(Counter(prepared.difference_counts[1:]).items())), prepared.profile_multiset


def hyperplane_invariant(prepared):
    mask = sum(1 << x for x in prepared.points)
    return tuple(sorted(Counter(hyperplane_intersection_profile(
        prepared.dimension, mask, len(prepared.points))).items()))


def sector(evidence, index, length, dimension):
    rows = [s for row in index["lengths"] if row["c"] == length
            for s in row["sectors"] if s["m"] == dimension]
    require(len(rows) == 1, "The source index must explicitly declare this sector")
    row, = rows
    require(row["data_complete"] and row["included_classes"] == row["expected_classes"],
            "The compact sector is incomplete")
    for item in row["encoding_certificates"]:
        evidence.digest(item["path"], item["sha256"])
    keys = [(length, dimension, i) for i in range(row["included_classes"])]
    supports, bindings = read_selected(evidence.root, keys)
    for name, digest in bindings.items():
        evidence.digest(name, digest)
    return [tuple(supports[key]) for key in keys]


def verify(root, length, dimension, witness_data=None):
    started = perf_counter()
    require(16 <= length <= 54 and length % 2 == 0 and 2 <= dimension <= 18,
            "Invalid finite support sector")
    require((1 << dimension)-1 > comb(length, 2), "The complete cover is not a zero-fibre cover")
    evidence = Evidence(root)
    index = evidence.read("data/spaces/index.json")
    parents = sector(evidence, index, length, dimension-1)
    representatives = sector(evidence, index, length, dimension)
    targets = []
    buckets = defaultdict(list)
    hyperplanes = {}
    for i, points in enumerate(representatives):
        validate_unital_support(points, dimension)
        prepared = prepare_support(points, dimension)
        key = invariant(prepared)
        for other in buckets[key]:
            if i not in hyperplanes:
                hyperplanes[i] = hyperplane_invariant(prepared)
            if other not in hyperplanes:
                hyperplanes[other] = hyperplane_invariant(targets[other])
            if hyperplanes[i] != hyperplanes[other]:
                continue
            require(transporter(prepared, targets[other]) is None, "Repeated affine class in the target catalogue")
        buckets[key].append(i)
        targets.append(prepared)
    supplied = None
    if witness_data is not None:
        require((witness_data["c"], witness_data["m"]) == (length, dimension), "Different affine witness scope")
        supplied = iter(witness_data["assignments"])
    profiles, assignments, occurrences = [], [], Counter()
    for source_index, points in enumerate(parents):
        quadratic_rank, basis = graph_labels(points, dimension-1)
        profiles.append(dict(source_index=source_index, quadratic_rank=quadratic_rank,
                             lift_dimension=len(basis), nonaffine_lifts=(1 << len(basis))-1))
        labels = 0
        for ordinal in range(1, 1 << len(basis)):
            labels ^= basis[(ordinal & -ordinal).bit_length()-1]
            child = reconstruct_lift(points, (), labels, dimension-1)
            validate_unital_support(child, dimension)
            prepared = prepare_support(child, dimension)
            if supplied is not None:
                row = next(supplied, None)
                require(row is not None and (row["source_index"], row["lift_ordinal"]) == (source_index, ordinal),
                        "Missing or out-of-order graph-lift witness")
                target_index, witness = row["target_index"], row["affine_map"]
                require(type(target_index) is int and 0 <= target_index < len(targets), "Invalid target index")
            else:
                witness = None
                choices = buckets[invariant(prepared)]
                if len(choices) > 1:
                    refined = hyperplane_invariant(prepared)
                    choices = [i for i in choices if hyperplanes[i] == refined]
                for target_index in choices:
                    witness = transporter(prepared, targets[target_index])
                    if witness is not None:
                        break
                require(witness is not None,
                        f"Unassigned graph lift: c={length}, m={dimension}, source={source_index}, ordinal={ordinal}, support={child}")
            require(witness_is_valid(prepared, targets[target_index], tuple(witness)), "Invalid graph-lift affine map")
            occurrences[target_index] += 1
            assignments.append(dict(source_index=source_index, lift_ordinal=ordinal,
                                    target_index=target_index, affine_map=list(witness)))
    require(supplied is None or next(supplied, None) is None, "Extraneous graph-lift witnesses")
    require(set(occurrences) == set(range(len(targets))), "A target class is not generated by the source recurrence")
    return dict(schema="complete-zero-fibre-space-recurrence-v1", status="pass", c=length, m=dimension,
        source_classes=len(parents), target_classes=len(targets), nonaffine_graph_lifts=len(assignments),
        source_profile=profiles, source_lift_dimension_histogram=dict(sorted(Counter(r["lift_dimension"] for r in profiles).items())),
        target_multiplicities=[occurrences[i] for i in range(len(targets))], assignments=assignments,
        complete_zero_fibre_cover_by_pair_averaging=True, all_graph_lifts_freshly_enumerated=True,
        every_affine_map_replayed=True, every_target_reached=True, targets_pairwise_affine_inequivalent=True,
        complete_given_predecessor_catalogue=True, predecessor_completeness_is_separate=True,
        affine_maps_freshly_searched=witness_data is None, elapsed_seconds=perf_counter()-started,
        dependencies=evidence.bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--witnesses", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.length, args.dimension,
                    json.loads(args.witnesses.read_bytes()) if args.witnesses else None)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("source_profile", "target_multiplicities", "assignments", "dependencies")}, indent=2))


if __name__ == "__main__":
    main()
