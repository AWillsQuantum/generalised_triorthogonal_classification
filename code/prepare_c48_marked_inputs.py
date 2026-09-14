"""Prepare the complete even-fibre contraction cover at length 48, dimension 8."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from affine_graph import affine_images, group_order, transform
from prepare_space_contractions import CoreExtensionContext, count_zero_xor_index_subsets, count_zero_xor_index_subsets_walsh
from space_codec import validate_unital_support
from verify_protocol_cover import Evidence, require
from verify_zero_fibre_spaces import sector


def stabiliser_generators(support, dimension, form):
    """Extend the intrinsic stabiliser to the seven-dimensional ambient space."""
    require(dimension in (6, 7), "Unsupported core span")
    validate_unital_support(support, dimension)
    mask = sum(1 << x for x in support)
    require(form["m"] == dimension and form["c"] == len(support)
            and transform(int(form["canonical_mask"], 16), form["frame"], dimension) == mask,
            "The affine certificate describes a different core")
    order = form["stabilizer_order"]
    require(type(order) is int and order > 0 and group_order(dimension) % order == 0,
            "Invalid intrinsic stabiliser order")
    require(form["generators"] or order == 1, "A nontrivial intrinsic stabiliser has no generators")
    generators = set()
    for frame in form["generators"]:
        images = tuple(affine_images(frame, dimension))
        require({images[x] for x in support} == set(support), "A generator does not preserve the core")
        if dimension == 6:
            images = tuple(images[x & 63] | (x & 64) for x in range(128))
        generators.add(images)
    if dimension == 6:
        # These maps fix the core hyperplane pointwise and generate its kernel.
        for bit in range(6):
            generators.add(tuple(x ^ ((1 << bit) if x & 64 else 0) for x in range(128)))
        order *= 64
    if not generators:
        require(order == 1, "A nontrivial stabiliser has no generators")
        generators.add(tuple(range(128)))
    return sorted(generators), order


def prepare(root, output):
    evidence = Evidence(root)
    index = evidence.read("data/spaces/index.json")
    base = evidence.read("certificates/low_dimensional_spaces.json")
    require(base["status"] == "pass" and base["complete_and_pairwise_inequivalent"], "The complete affine base is required")
    forms = {(r["c"], r["m"], r["index"]):r for r in base["canonical_maps"]}
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for c, m in ((48, 7), (44, 7), (40, 6), (40, 7)):
        for i, support in enumerate(sector(evidence, index, c, m)):
            context = CoreExtensionContext.build(7, support, m+1)
            generators, order = stabiliser_generators(support, m, forms[c, m, i])
            n = (48-c)//2
            fibres = count_zero_xor_index_subsets(context.external_labels, n)
            require(fibres == count_zero_xor_index_subsets_walsh(context.external_labels, n),
                    "The independent compatible-fibre counts disagree")
            ell = len(context.lift_complement)
            require(0 <= ell <= 18, "Unsupported lift quotient")
            excluded = int(n == 0)
            raw = fibres*(1 << ell)-excluded
            identifier = f"c{c}_m{m:02d}_{i:09d}"
            lines = ["L48_M08_NATIVE_TASK_V1", identifier,
                     " ".join(map(str, (c, n, c, 128-c, m+1, ell, len(generators), order, fibres, raw, excluded)))]
            for values in (context.core_points, context.quadratic_columns, context.external_points,
                           context.external_labels, context.affine_rows, context.lift_complement, *generators):
                lines.append(" ".join(map(str, values)))
            payload = ("\n".join(lines)+"\n").encode("ascii")
            path = output / (identifier+".task")
            require(not path.exists() or path.read_bytes() == payload, "A different native input exists")
            path.write_bytes(payload)
            records.append(dict(c=c, m=m, source_index=i, multiplicity=n,
                                compatible_fibres=fibres, ell=ell, raw_pairs=raw,
                                proper_span_requires_output_rank_filter=m < 7,
                                stabiliser_order=order, task=path.name,
                                sha256=hashlib.sha256(payload).hexdigest()))
    counts = Counter((row["c"], row["m"]) for row in records)
    require(counts == {(48, 7):98, (44, 7):35, (40, 6):1, (40, 7):18}, "The contraction source cover is incomplete")
    result = dict(schema="c48-marked-contraction-domain-v1", c=48, m=8,
                  source_count=len(records), raw_pairs=sum(row["raw_pairs"] for row in records),
                  proper_span_cores_included=True, records=records, dependencies=evidence.bindings)
    (output / "domain.json").write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.root.resolve(), args.output.resolve())
    print(json.dumps({k:v for k,v in result.items() if k not in ("records", "dependencies")}))


if __name__ == "__main__":
    main()
