"""Reproduce a conditional resource estimate, not a classification certificate."""

import argparse
from collections import Counter
from functools import reduce
from itertools import combinations
import json
import math
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "resource_estimates/length56_inputs.json"


def binary_basis(values):
    basis = {}
    for value in values:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return basis


def in_span(value, basis):
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            return False
        value ^= basis[pivot]
    return True


def action_orbits(points, maps):
    unseen = set(points)
    result = []
    while unseen:
        first = min(unseen)
        orbit, pending = {first}, [first]
        while pending:
            point = pending.pop()
            for permutation in maps:
                other = permutation[point]
                if other not in orbit:
                    orbit.add(other)
                    pending.append(other)
        if not orbit <= unseen:
            raise ValueError("The supplied domain is not invariant")
        unseen -= orbit
        result.append(orbit)
    return result


def support_profile(points):
    """Recheck a 56-point, eight-dimensional pilot support using igraph/Bliss."""
    import igraph

    points = sorted(points)
    if len(points) != 56 or len(set(points)) != 56 or not all(0 <= x < 256 for x in points):
        raise ValueError("Expected 56 distinct eight-bit points")
    columns = [1 | (x << 1) for x in points]
    if len(binary_basis(columns)) != 9:
        raise ValueError("The support must affinely span dimension eight")
    if any(sum(all(x >> j & 1 for j in term) for x in points) % 2
           for degree in range(4) for term in combinations(range(8), degree)):
        raise ValueError("Invalid unital triorthogonal support")

    # All codewords, not just generator rows, make coordinate automorphisms exact.
    words = [sum(((word & column).bit_count() % 2) << j
                 for j, column in enumerate(columns)) for word in range(512)]
    graph = igraph.Graph(n=568, edges=[(j, 56 + w) for w, row in enumerate(words)
                                     for j in range(56) if row >> j & 1])
    colours = [0] * 56 + [1] * 512
    order = int(graph.count_automorphisms(color=colours))
    point_maps, direction_maps = [], []
    support = set(points)
    for permutation in graph.automorphism_group(color=colours):
        images = [permutation[56 + (1 << i)] - 56 for i in range(9)]
        if images[0] != 1 or any(
            permutation[56 + word] - 56 != reduce(
                int.__xor__, (images[i] for i in range(9) if word >> i & 1), 0)
            for word in range(512)
        ):
            raise ValueError("Codeword permutation is not linear or does not fix one")
        point_map = [sum(((image & (1 | (x << 1))).bit_count() % 2) << i
                         for i, image in enumerate(images)) >> 1 for x in range(256)]
        direction_map = [point_map[x] ^ point_map[0] for x in range(256)]
        if {point_map[x] for x in points} != support:
            raise ValueError("Recovered affine map does not preserve the support")
        point_maps.append(point_map)
        direction_maps.append(direction_map)

    pair_counts = Counter(x ^ y for x, y in combinations(points, 2))
    admissible = [a for a in range(1, 256) if pair_counts[a] <= 6]
    for direction in admissible:
        derivative_support = support ^ {x ^ direction for x in support}
        if len(binary_basis(1 | (x << 1) for x in derivative_support)) != 9:
            raise ValueError("A contracted core lies outside the calibrated full-rank census")
    directions = action_orbits(admissible, direction_maps)
    marked_orbit = next(orbit for orbit in directions if 1 in orbit)
    if order % len(marked_orbit):
        raise ValueError("Invalid orbit-stabiliser quotient")

    terms = [()] + [(i,) for i in range(8)] + list(combinations(range(8), 2))
    features = [sum(all(x >> j & 1 for j in term) << i
                    for i, term in enumerate(terms)) for x in range(256)]
    quadratic_basis = binary_basis(features[x] for x in points)
    rank = len(quadratic_basis)
    histogram = Counter()
    for orbit in action_orbits(range(256), point_maps):
        origin = min(orbit)
        quotient_dimension = 48 - rank + int(
            origin not in support and in_span(features[origin], quadratic_basis))
        histogram[str(quotient_dimension)] += 1
    # Also retain the whole unital space as the stabiliser space.
    histogram[str(47 - rank)] += 1
    return dict(affine_automorphism_order=order,
                marked_automorphism_order=order // len(marked_orbit),
                admissible_direction_orbits=len(directions),
                quadratic_evaluation_rank=rank,
                quotient_dimension_histogram=dict(sorted(histogram.items())))


def estimate(data):
    by_stratum = {}
    core_forecasts = []
    for core in data["sampled_n4_cores"]:
        weight = sum(row["profile"]["marked_automorphism_order"] for row in core["samples"])
        counts = {}
        for dimension in (15, 16):
            counts[str(dimension)] = core["marked_orbit_count"] * sum(
                row["profile"]["marked_automorphism_order"]
                * row["profile"]["quotient_dimension_histogram"].get(str(dimension), 0)
                / row["profile"]["admissible_direction_orbits"]
                for row in core["samples"]
            ) / weight
        stratum = str(core["quadratic_nullity"])
        group = by_stratum.setdefault(stratum, dict(mass=0, counts=Counter()))
        group["mass"] += core["marked_orbit_count"]
        group["counts"].update(counts)
        core_forecasts.append(dict(core_id=core["core_id"], counts=counts))

    rates = {s: {k: row["counts"][k] / row["mass"] for k in ("15", "16")}
             for s, row in by_stratum.items()}
    # No positive observation is a forecast assumption, never an exclusion rule.
    proxies = {k: max(row[k] for row in rates.values()) for k in ("15", "16")}
    totals = Counter()
    strata = []
    for row in data["n4_population_strata"]:
        key = str(row["quadratic_nullity"])
        rate = rates.get(key, proxies)
        count = {k: row["marked_orbit_count"] * rate[k] for k in ("15", "16")}
        totals.update(count)
        strata.append(dict(**row, sampled=key in rates, estimated_hard_supports=count))

    timings = {k: mean(v) / 3600 for k, v in data["timing_core_seconds"].items()}
    tail = sum(totals[k] * timings[k] for k in ("15", "16"))
    background = sum(data["background_allowance_core_hours"].values())
    subtotal = tail + background
    allowance = subtotal * data["additional_work_allowance_fraction"]
    total = subtotal + allowance
    exponent = math.floor(math.log10(total) + 0.5)
    return dict(
        schema="length56-resource-estimate-v1", date=data["date"],
        status="conditional_engineering_forecast_not_a_bound_or_confidence_interval",
        estimate_core_hours=10 ** exponent, exponent=exponent,
        computational_scope=data["scope"],
        arithmetic=dict(mean_hard_key_core_hours=timings,
                        projected_hard_supports=dict(totals),
                        projected_hard_filter_core_hours=tail,
                        background_allowance_core_hours=background,
                        subtotal_core_hours=subtotal,
                        additional_work_allowance_core_hours=allowance,
                        unrounded_planning_total_core_hours=total),
        strata=strata, sampled_core_forecasts=core_forecasts,
        assumptions=data["assumptions"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--check-sample", action="store_true",
                        help="Recompute all retained support profiles; requires python-igraph")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    if args.check_sample:
        for core in data["sampled_n4_cores"]:
            for row in core["samples"]:
                if support_profile(row["support"]) != row["profile"]:
                    raise ValueError(f"Changed profile: {core['core_id']} / {row['sample_index']}")
    result = estimate(data)
    content = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(content, encoding="ascii")
    print(content, end="")


if __name__ == "__main__":
    main()
