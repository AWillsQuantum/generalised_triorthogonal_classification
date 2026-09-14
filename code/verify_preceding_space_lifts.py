"""Differential checks of configured native lifts against independent binary equations."""

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import subprocess

from enumerate_preceding_space_interval import generate
from finite_space_census import points_of
from prepare_contraction_interval import reference_profile
from prepare_space_contractions import CoreExtensionContext
from selected_spaces import iter_interval
from space_codec import validate_unital_support
from space_lifts import evaluation_rows, iter_lifts, reduce_vector
from verify_protocol_cover import Evidence, require


def minimum_multiplicity(points, dimension):
    counts = Counter(a ^ b for a, b in combinations(points, 2))
    return 0 if len(counts) < (1 << dimension)-1 else min(counts.values())


def shear_key(points, core, dimension, multiplicity, pivots):
    fibres = Counter(x >> 1 for x in points)
    full = tuple(sorted(x for x, count in fibres.items() if count == 2))
    require(tuple(sorted(x for x, count in fibres.items() if count == 1)) == core
            and len(full) == multiplicity and all(count in (1, 2) for count in fibres.values()),
            "The generated support contracts to a different core")
    labels = {x >> 1: x & 1 for x in points}
    vector = sum(labels[x] << i for i, x in enumerate(core))
    return full, reduce_vector(vector, pivots)


def verify(root, native, work, maximum_reference_lifts=262144):
    evidence = Evidence(root)
    native, work = Path(native).resolve(), Path(work).resolve()
    require(not work.exists(), "Use a new differential-check directory")
    work.mkdir(parents=True)
    evidence.digest("code/space_native/generic/signature.cpp")
    records = []
    for c in (50, 52):
        for m in (9, 10, 11):
            for n in range(max(1, 12-m)):
                bindings = {}
                candidates = []
                for i, points in iter_interval(evidence.root, c-2*n, m-1, 0, 8, bindings):
                    rank, ell, fibres = reference_profile(points, m, n)
                    count = fibres*(1 << ell)-int(n == 0)
                    if 0 < count <= maximum_reference_lifts:
                        candidates.append((i, points, count))
                for path, digest in bindings.items():
                    evidence.digest(path, digest)
                require(candidates, "No bounded positive source in the reference prefix")
                for i, core, count in sorted(candidates, key=lambda row: (row[2], row[0]))[:2]:
                    case = work / f"c{c}_m{m:02d}_n{n}_{i:09d}"
                    result = generate(root, c, m, n, i, 1, native, case, maximum_reference_lifts, 2)
                    for path, digest in result["dependencies"].items():
                        evidence.digest(path, digest)
                    require(result["full_rank_marked_lifts"] == count, "Independent raw lift counts disagree")
                    pivots = {}
                    for row in evaluation_rows(core, m-1, 1):
                        value = reduce_vector(row, pivots)
                        if value:
                            pivots[value.bit_length()-1] = value
                    expected, raw_count = set(), 0
                    context = CoreExtensionContext.build(m-1, core, m)
                    for fibres in context.iter_fiber_sets(n, 0):
                        for support in iter_lifts(core, fibres, m-1):
                            raw_count += 1
                            if not result["minimum_direction_filter"] or minimum_multiplicity(support, m) == n:
                                expected.add(shear_key(support, core, m, n, pivots))
                    require(raw_count == count, "Independent reference equations give a different lift count")
                    actual, width = set(), (1 << m)//8
                    with (case / "candidates.bin").open("rb") as stream:
                        while payload := stream.read(width):
                            require(len(payload) == width, "Truncated candidate mask")
                            support = points_of(int.from_bytes(payload, "little"))
                            validate_unital_support(support, m)
                            if result["minimum_direction_filter"]:
                                require(minimum_multiplicity(support, m) == n, "A retained support violates the minimum filter")
                            actual.add(shear_key(support, core, m, n, pivots))
                    require(actual == expected and len(actual) == result["candidates"],
                            "Native and independent affine-shear classes disagree")
                    if result["candidates"] and m != 9:
                        other = case / "independent_signature.bin"
                        subprocess.run([str(native / f"signature_c{c}_m{m:02d}"),
                            "--input", str(case / "candidates.bin"), "--output", str(other), "--threads", "2"],
                            text=True, capture_output=True, check=True)
                        require(other.read_bytes() == (case / "ledger.bin").read_bytes(),
                                "The fused and separate affine signatures disagree")
                    records.append(dict(c=c, m=m, multiplicity=n, source_index=i,
                                        full_rank_lifts=count, retained_lifts=result["candidates"]))
                    print(json.dumps(records[-1]), flush=True)
    return dict(schema="preceding-space-native-lift-verification-v1", status="pass", cases=records,
        all_native_outputs_independently_validated=True,
        all_families_equal_modulo_affine_shears=True,
        minimum_direction_filters_independently_replayed=True,
        fused_and_separate_signatures_agree=True,
        is_global_completeness_certificate=False, dependencies=evidence.bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--native-directory", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.native_directory, args.work_directory)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps(dict(status=result["status"], cases=len(result["cases"]))))


if __name__ == "__main__":
    main()
