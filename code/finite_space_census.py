"""Finite inverse-contraction enumeration with an independently replayed affine quotient."""

import argparse
from collections import Counter
import json
from pathlib import Path
import struct
import subprocess
import sys
from time import perf_counter

import numpy as np

from affine_codeword_graph import canonicalise_many
from prepare_zero_fibre_input import file_digest
from space_batch import decode_bitmaps, validate_batch
from space_codec import binary_rank, validate_unital_support
from space_recursion import multiplicity_bound
from verify_protocol_cover import Evidence, require
from verify_zero_fibre_ledger import build_buckets
from verify_zero_fibre_spaces import sector

sys.path.insert(0, str(Path(__file__).resolve().parent / "space_algorithms"))
from inverse_contraction_extensions import CoreExtensionContext, count_zero_xor_index_subsets, xor_quadratic_evaluations


def points_of(mask):
    points = []
    while mask:
        bit = mask & -mask
        points.append(bit.bit_length()-1)
        mask ^= bit
    return tuple(points)


def generate_domain(root, length, dimension, replacement_files, output, maximum_candidates):
    evidence = Evidence(root)
    index = evidence.read("data/spaces/index.json")
    replacements = {}
    for path in replacement_files:
        data = json.loads(path.read_bytes())
        key = data["c"], data["m"]
        require(key not in replacements, "Repeated supplied predecessor sector")
        replacements[key] = data["class_supports"], file_digest(path)
    bound = multiplicity_bound(length, dimension)
    require(8 <= dimension <= 15 and 0 <= bound <= 5 and 2*bound < length, "Unsupported finite contraction domain")
    domains, contexts, estimated = [], [], 0
    for multiplicity in range(bound+1):
        c = length-2*multiplicity
        lower = dimension-1 if multiplicity <= 3 else max(4, (c-1).bit_length())
        upper = min(dimension-1, (c+c//16)//3-1)
        for m in range(lower, upper+1):
            key = c, m
            if key in replacements:
                supports, source_digest = replacements[key]
                supports = [tuple(p) for p in supports]
            else:
                supports, source_digest = sector(evidence, index, c, m), None
            domains.append(dict(c=c, m=m, multiplicity=multiplicity, classes=len(supports), supplied_sector_sha256=source_digest))
            for i, points in enumerate(supports):
                validate_unital_support(points, m)
                context = CoreExtensionContext.build(dimension-1, points, m+1)
                fibres = count_zero_xor_index_subsets(context.external_labels, multiplicity)
                ell = len(context.lift_complement)
                raw = fibres*(1 << ell)-(multiplicity == 0)
                require(raw >= 0, "Negative full-rank candidate count")
                estimated += raw
                require(estimated <= maximum_candidates, "The exact finite domain exceeds the requested candidate budget")
                contexts.append((context, dict(c=c, m=m, source_index=i, multiplicity=multiplicity,
                    quadratic_rank=context.solver.rank, ell=ell, compatible_fibres=fibres, candidate_upper_bound=raw)))
    width = (1 << dimension)//8
    count, profiles = 0, []
    with output.open("xb") as stream:
        for context, profile in contexts:
            n = profile["multiplicity"]
            def toggle(labels):
                return sum(3 << (2*x) for i, x in enumerate(context.core_points) if labels >> i & 1)
            base = sum(1 << (2*x) for x in context.core_points)
            toggles = [toggle(row) for row in context.lift_complement]
            generated = fibre_count = 0
            for fibres in context.iter_fiber_sets(n, 0):
                fibre_count += 1
                rhs = xor_quadratic_evaluations(fibres, dimension-1)
                particular = context.solver.solve(rhs)
                require(particular is not None, "A compatible fibre set has no solution")
                value = base ^ toggle(particular)
                value |= sum(3 << (2*x) for x in fibres)
                for ordinal in range(1 << len(toggles)):
                    if ordinal:
                        value ^= toggles[(ordinal & -ordinal).bit_length()-1]
                    if n == 0 and ordinal == 0:
                        continue
                    if profile["m"] != dimension-1:
                        points = points_of(value)
                        if binary_rank(x ^ points[0] for x in points) != dimension:
                            continue
                    stream.write(value.to_bytes(width, "little"))
                    generated += 1
            require(fibre_count == profile["compatible_fibres"], "The independent fixed-cardinality fibre counts disagree")
            require(generated <= profile["candidate_upper_bound"]
                    and (profile["m"] != dimension-1 or generated == profile["candidate_upper_bound"]),
                    "The complete graph-label enumeration does not close")
            count += generated
            profiles.append(dict(**profile, candidates=generated))
    require(output.stat().st_size == count*width, "The generated bitmap stream is incomplete")
    return dict(source_domains=domains, source_profiles=profiles, candidates=count,
                candidate_sha256=file_digest(output), dependencies=evidence.bindings,
                full_span_core_reduction_for_multiplicity_at_most_three=True)


def replay_quotient(ledger, bucket_path, class_path, assignment_path, length, dimension):
    width = (1 << dimension)//8
    class_width = 24+width
    raw_classes = class_path.read_bytes()
    require(len(raw_classes) % class_width == 0, "Truncated class stream")
    by_bucket, class_supports, member_counts = {}, [], []
    for offset in range(0, len(raw_classes), class_width):
        raw = raw_classes[offset:offset+class_width]
        bucket, local, representative = struct.unpack_from("<QII", raw)
        mask = int.from_bytes(raw[16:16+width], "little")
        points = points_of(mask)
        validate_unital_support(points, dimension)
        members, = struct.unpack_from("<Q", raw, 16+width)
        rows = by_bucket.setdefault(bucket, [])
        require(local == len(rows) and members > 0, "Invalid local class index or multiplicity")
        rows.append((points, members, representative))
        class_supports.append(points)
        member_counts.append(members)
    assignment_width = 4*((9+2*dimension)//4)
    dtype = np.dtype(dict(names=["class", "origin", "basis"], formats=["<u4", "<u2", ("<u2", dimension)],
                          offsets=[0, 4, 6], itemsize=assignment_width))
    count = buckets = 0
    with ledger.open("rb") as incoming, bucket_path.open("rb") as indices, assignment_path.open("rb") as assignments:
        while raw := indices.read(48+2*width):
            require(len(raw) == 48+2*width, "Truncated signature bucket")
            first, size = struct.unpack_from("<QQ", raw, 32)
            require(first == count and size > 0 and buckets in by_bucket, "Incomplete signature partition")
            classes = by_bucket[buckets]
            require(all(row[2] < size for row in classes), "A representative position lies outside its bucket")
            sources = np.asarray([row[0] for row in classes], dtype=np.uint32)
            observed = np.zeros(len(classes), dtype=np.int64)
            remaining = size
            position = 0
            while remaining:
                batch = min(4096, remaining, max(1, 32*1024*1024//(1 << dimension)))
                records = incoming.read(batch*(32+width))
                require(len(records) == batch*(32+width), "Truncated candidate ledger")
                records = np.frombuffer(records, dtype=np.uint8).reshape(batch, 32+width)
                require(np.all(records[:, :32] == np.frombuffer(raw[:32], dtype=np.uint8)), "Different candidate signature")
                targets = decode_bitmaps(records[:, 32:].copy().tobytes(), dimension, length)
                validate_batch(targets, dimension, length)
                encoded = assignments.read(batch*assignment_width)
                require(len(encoded) == batch*assignment_width, "Truncated affine assignment stream")
                maps = np.frombuffer(encoded, dtype=dtype)
                require(np.all(maps["class"] < len(classes)) and np.all(maps["origin"] < 1 << dimension)
                        and np.all(maps["basis"] < 1 << dimension), "Invalid affine map range")
                selected = sources[maps["class"]]
                images = np.broadcast_to(maps["origin"][:, None], selected.shape).astype(np.uint32).copy()
                for bit in range(dimension):
                    images ^= np.where((selected >> bit) & 1, maps["basis"][:, bit, None], 0).astype(np.uint32)
                require(np.array_equal(np.sort(images, axis=1), targets), "An affine assignment does not reproduce its candidate")
                # A linear map onto a full-affine-rank target is invertible.
                observed += np.bincount(maps["class"], minlength=len(classes))
                for local, (points, _, representative) in enumerate(classes):
                    if position <= representative < position+batch:
                        require(tuple(targets[representative-position]) == points, "A class representative is not its indicated candidate")
                remaining -= batch
                position += batch
            require(observed.tolist() == [row[1] for row in classes], "Affine-class member accounting fails")
            count += size
            buckets += 1
        require(not incoming.read(1) and not assignments.read(1) and buckets == len(by_bucket), "Trailing or missing quotient data")
    return dict(candidate_count=count, bucket_count=buckets, affine_classes=len(class_supports),
                class_supports=class_supports, member_counts=member_counts)


def run(root, c, m, replacements, native, graph_native, work, maximum_candidates, threads):
    started = perf_counter()
    work.mkdir(parents=True, exist_ok=True)
    masks, ledger, buckets, classes, assignments, negatives = [work / name for name in
        ("candidates.bin", "signatures.bin", "buckets.bin", "classes.bin", "assignments.bin", "negatives.bin")]
    require(not any(p.exists() for p in (masks, ledger, buckets, classes, assignments, negatives)), "A finite-census output exists")
    domain = generate_domain(root, c, m, replacements, masks, maximum_candidates)
    (work / "domain.json").write_text(json.dumps(domain, indent=2)+"\n", encoding="ascii")
    print(json.dumps(dict(c=c, m=m, stage="generated", candidates=domain["candidates"])), flush=True)
    return quotient_generated_domain(root, c, m, domain, native, graph_native, work, started, threads)


def quotient_generated_domain(root, c, m, domain, native, graph_native, work, started=None, threads=4):
    """Exactly quotient a complete supplied candidate stream, independently replaying it."""
    if started is None:
        started = perf_counter()
    masks, ledger, buckets, classes, assignments, negatives = [work / name for name in
        ("candidates.bin", "signatures.bin", "buckets.bin", "classes.bin", "assignments.bin", "negatives.bin")]
    require(not any(p.exists() for p in (ledger, buckets, classes, assignments, negatives)), "A quotient output exists")
    require(masks.stat().st_size == domain["candidates"]*((1 << m)//8)
            and file_digest(masks) == domain["candidate_sha256"], "The supplied candidate stream differs")
    if not domain["candidates"]:
        checked, quotient = dict(candidate_count=0, affine_classes=0, class_supports=[], member_counts=[]), {}
    else:
        result = subprocess.run([str(native/f"signature_c{c}_m{m:02d}"), "--input", str(masks),
            "--output", str(ledger), "--threads", str(threads)], capture_output=True, text=True, check=True)
        (work / "signature.json").write_text(result.stdout, encoding="ascii")
        indexed = build_buckets(ledger, buckets, m, c)
        require(indexed["candidate_count"] == domain["candidates"], "Signature generation lost candidates")
        result = subprocess.run([str(native/f"exact_quotient_c{c}_m{m:02d}"), "--ledger", str(ledger),
            "--bucket-index", str(buckets), "--first-bucket", "0", "--last-bucket", str(indexed["signature_bucket_count"]),
            "--assignments", str(assignments), "--classes", str(classes), "--negatives", str(negatives)],
            capture_output=True, text=True, check=True)
        quotient = json.loads(result.stdout)
        (work / "quotient.json").write_text(result.stdout, encoding="ascii")
        checked = replay_quotient(ledger, buckets, classes, assignments, c, m)
        require(checked["candidate_count"] == domain["candidates"] == quotient["candidate_count"]
                and checked["affine_classes"] == quotient["affine_class_count"], "The exact quotient census differs")
    evidence = Evidence(root)
    old = sector(evidence, evidence.read("data/spaces/index.json"), c, m)
    forms = canonicalise_many([graph_native], [(m, p) for p in (*checked["class_supports"], *old)]) if old or checked["class_supports"] else []
    new_keys = [tuple(row["canonical_points"]) for row in forms[:checked["affine_classes"]]]
    old_keys = [tuple(row["canonical_points"]) for row in forms[checked["affine_classes"]:]]
    new_key_set, old_key_set = set(new_keys), set(old_keys)
    require(len(new_key_set) == checked["affine_classes"], "The independent canonical graph finds equivalent quotient classes")
    result = dict(schema="complete-finite-contraction-census-v1", status="pass", c=c, m=m, **checked,
        domain=domain, exact_quotient=quotient, every_candidate_validated=True, every_positive_affine_witness_replayed=True,
        every_class_independently_separated=True, complete_given_predecessor_catalogues=True,
        predecessor_completeness_is_separate=True, old_target_classes=len(old),
        new_classes=len(new_key_set-old_key_set), unreached_old_targets=[i for i, key in enumerate(old_keys) if key not in new_key_set],
        graph_results=forms, elapsed_seconds=perf_counter()-started, comparison_dependencies=evidence.bindings)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--source-sector", type=Path, action="append", default=[])
    parser.add_argument("--native-directory", type=Path, required=True)
    parser.add_argument("--graph-native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--maximum-candidates", type=int, default=2_000_000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.root.resolve(), args.length, args.dimension, args.source_sector, args.native_directory.resolve(),
                 args.graph_native.resolve(), args.work_directory.resolve(), args.maximum_candidates, args.threads)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k in ("status", "c", "m", "candidate_count", "affine_classes", "old_target_classes", "new_classes", "unreached_old_targets", "elapsed_seconds")}, indent=2))


if __name__ == "__main__":
    main()
