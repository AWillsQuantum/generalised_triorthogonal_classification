"""Reconstruct complete dimension-eight contraction inputs from public data."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

from affine_graph import affine_images, group_order, transform
from read_spaces import open_shard
from space_codec import read_header, read_records, validate_unital_support

sys.path.insert(0, str(Path(__file__).resolve().parent / "space_algorithms"))
from inverse_contraction_extensions import (CoreExtensionContext,
    count_zero_xor_index_subsets, count_zero_xor_index_subsets_walsh)


def task_text(identifier, support, form, target_length=54):
    """Return native input and exact counts for c=50,52,54 and m=8.

    The supplied stabiliser must come from a complete affine-graph census;
    checking its maps alone would not establish that it is the full group.
    """
    support = tuple(support)
    validate_unital_support(support, 7)
    weight = len(support)
    if weight not in (44, 48, 52) or not identifier or any(x.isspace() for x in identifier):
        raise ValueError("Invalid contraction core or identifier")
    if target_length not in (50, 52, 54) or weight > target_length:
        raise ValueError("Unsupported contraction target")
    mask = sum(1 << x for x in support)
    if (form["m"] != 7 or form["c"] != weight
            or transform(int(form["canonical_mask"], 16), form["frame"], 7) != mask):
        raise ValueError("The affine certificate describes another core")
    order = form["stabilizer_order"]
    if type(order) is not int or order <= 0 or group_order(7) % order:
        raise ValueError("Invalid stabiliser order")
    generators = set()
    for frame in form["generators"]:
        images = tuple(affine_images(frame, 7))
        if {images[x] for x in support} != set(support):
            raise ValueError("A generator does not preserve the core")
        generators.add(images)
    if not generators:
        if order != 1:
            raise ValueError("A nontrivial stabiliser has no generators")
        generators.add(tuple(range(128)))
    generators = sorted(generators)
    context = CoreExtensionContext.build(7, support, 8)
    multiplicity = (target_length-weight)//2
    fibres = count_zero_xor_index_subsets(context.external_labels, multiplicity)
    if fibres != count_zero_xor_index_subsets_walsh(context.external_labels, multiplicity):
        raise ValueError("Independent fibre counts disagree")
    dimension = len(context.lift_complement)
    if not 0 <= dimension <= 18:
        raise ValueError("Lift quotient exceeds the native kernel domain")
    excluded = int(multiplicity == 0)
    raw = fibres * (1 << dimension) - excluded
    lines = [f"L{target_length}_M08_NATIVE_TASK_V1", identifier,
             " ".join(map(str, (weight, multiplicity, weight, 128-weight, 8,
                                  dimension, len(generators), order, fibres, raw, excluded)))]
    for values in (context.core_points, context.quadratic_columns,
                   context.external_points, context.external_labels,
                   context.affine_rows, context.lift_complement, *generators):
        lines.append(" ".join(map(str, values)))
    profile = dict(core_length=weight, core_affine_dimension=7,
                   full_fibres=multiplicity, lift_quotient_dimension=dimension,
                   compatible_fibre_sets=fibres, raw_marked_pairs=raw,
                   rank_deficient_pairs=excluded, stabiliser_order=order)
    return "\n".join(lines) + "\n", profile


def prepare(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    certificate_bytes = (root / "certificates/low_dimensional_spaces.json").read_bytes()
    certificate = json.loads(certificate_bytes)
    if (certificate.get("status") != "pass"
            or certificate.get("complete_and_pairwise_inequivalent") is not True
            or certificate.get("maximum_length") != 54
            or certificate.get("maximum_affine_dimension") != 7):
        raise ValueError("The complete low-dimensional certificate is required")
    forms = {}
    for row in certificate["canonical_maps"]:
        key = row["c"], row["m"], row["index"]
        if key in forms:
            raise ValueError("Repeated catalogue index in affine certificate")
        forms[key] = row
    output.mkdir(parents=True, exist_ok=True)
    records, bindings = [], {}
    for name, digest in sorted(certificate["catalogue_sha256"].items()):
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError("External catalogue reference")
        with open_shard(path) as stream:
            header = read_header(stream)
            if header.m != 7 or header.c not in (44, 48, 52):
                continue
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise ValueError("Core catalogue digest mismatch")
            bindings[name] = digest
            for offset, support in enumerate(read_records(stream, header)):
                index = header.first_index + offset
                identifier = f"c{header.c}_m07_{index:09d}"
                text, profile = task_text(identifier, support, forms[header.c, 7, index])
                payload = text.encode("ascii")
                task = output / f"{identifier}.task"
                if task.exists() and task.read_bytes() != payload:
                    raise ValueError("Refusing to overwrite a different contraction input")
                task.write_bytes(payload)
                records.append(dict(id=identifier, core_index=index, **profile,
                                    task=task.name, sha256=hashlib.sha256(payload).hexdigest()))
    counts = Counter(row["core_length"] for row in records)
    if dict(counts) != {44: 35, 48: 98, 52: 188} or len({r["id"] for r in records}) != 321:
        raise ValueError("The complete contraction core domain was not reconstructed")
    report = dict(schema="space-contraction-domain-v1", target_length=54,
                  target_affine_dimension=8, input_domain_complete=True,
                  child_enumeration_performed=False, cores=len(records),
                  core_counts={str(c): n for c, n in sorted(counts.items())},
                  compatible_fibre_sets=sum(r["compatible_fibre_sets"] for r in records),
                  raw_marked_pairs=sum(r["raw_marked_pairs"] for r in records),
                  affine_census_sha256=hashlib.sha256(certificate_bytes).hexdigest(),
                  catalogue_sha256=bindings, records=records)
    (output / "domain.json").write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.root, args.output)
    print(json.dumps({k: v for k, v in report.items() if k not in ("records", "catalogue_sha256")}, indent=2))
