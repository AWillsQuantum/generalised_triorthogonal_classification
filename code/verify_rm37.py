"""Certify all cubic Boolean-function orbits by exact graph forms and orbit mass."""

import argparse
from collections import Counter
import hashlib
import json
from math import comb
from pathlib import Path

from affine_graph import canonicalise_many, group_order, rank


def polynomial_support(monomials, m, degree):
    if len(set(monomials)) != len(monomials):
        raise ValueError("Duplicate polynomial monomial")
    if any(type(t) is not int or not 0 <= t < 1 << m or t.bit_count() > degree
           for t in monomials):
        raise ValueError("Polynomial lies outside the specified Reed-Muller code")
    return sum(1 << x for x in range(1 << m)
               if sum((x & t) == t for t in monomials) % 2)


def certify(data, command):
    m, degree = data["m"], data["maximum_degree"]
    if (m, degree) != (7, 3):
        raise ValueError("This source certificate concerns RM(3,7)")
    seeds = data["records"]
    if [r["index"] for r in seeds] != list(range(len(seeds))):
        raise ValueError("Incomplete or repeated source indices")
    masks = [polynomial_support(r["monomial_masks"], m, degree) for r in seeds]
    forms = canonicalise_many(command, [(m, mask) for mask in masks])
    canonical = [int(row["canonical_mask"], 16) for row in forms]
    if len(set(canonical)) != len(seeds):
        raise ValueError("Repeated affine class")
    ambient_order = group_order(m)
    mass = sum(ambient_order // row["stabilizer_order"] for row in forms)
    function_count = 1 << sum(comb(m, i) for i in range(degree + 1))
    if mass != function_count:
        raise ValueError("Affine orbit mass does not exhaust the Reed-Muller code")
    weights = Counter()
    full_weights = Counter()
    for seed, mask, row in zip(seeds, masks, forms, strict=True):
        support = tuple(x for x in range(1 << m) if mask >> x & 1)
        affine_rank = rank(x ^ support[0] for x in support) if support else -1
        row.update(index=seed["index"], weight=len(support), affine_rank=affine_rank,
                   orbit_size=ambient_order // row["stabilizer_order"])
        weights[len(support)] += 1
        if affine_rank == m:
            full_weights[len(support)] += 1
    summary = dict(schema="rm37-affine-orbit-mass-certificate-v1", status="pass",
                   m=m, degree=degree, classes=len(forms), affine_group_order=ambient_order,
                   exact_orbit_mass=mass, function_count=function_count,
                   canonical_supports_pairwise_distinct=True,
                   all_polynomials_and_affine_maps_replayed=True,
                   method="Point-hyperplane incidence graph; Bliss canonical forms and GMP exact stabiliser orders",
                   weight_counts=dict(sorted(weights.items())),
                   full_affine_rank_weight_counts=dict(sorted(full_weights.items())))
    return summary, forms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=Path,
                        default=Path(__file__).resolve().parents[1] / "data/seeds/rm37.json")
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    raw = args.seeds.read_bytes()
    summary, rows = certify(json.loads(raw), [args.native])
    args.output_directory.mkdir(parents=True, exist_ok=True)
    results = args.output_directory / "orbits.json"
    results.write_text(json.dumps(dict(m=7, records=rows), indent=2) + "\n")
    summary.update(seeds_sha256=hashlib.sha256(raw).hexdigest(),
                   orbits_sha256=hashlib.sha256(results.read_bytes()).hexdigest())
    (args.output_directory / "certificate.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
