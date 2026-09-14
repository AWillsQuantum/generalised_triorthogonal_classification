"""Cross-check native point-colour refinement and its affine covariance."""

import argparse
from bisect import bisect_left
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import subprocess

from space_codec import binary_rank
from verify_protocol_cover import require


def refine(points):
    points = tuple(sorted(points))
    differences = Counter(x ^ y for x in points for y in points if x != y)
    profiles = [tuple(sorted(differences[x ^ y] for y in points if y != x)) for x in points]
    palette = sorted(profiles)
    colours = [bisect_left(palette, p) for p in profiles]
    for _ in points:
        keys = [(colours[i], tuple(sorted((colours[j], differences[x ^ y])
                                         for j, y in enumerate(points) if i != j))) for i, x in enumerate(points)]
        palette = sorted(keys)
        refined = [bisect_left(palette, k) for k in keys]
        if refined == colours:
            return colours
        colours = refined
    raise ValueError("The point-colour partition did not stabilise")


def apply(point, frame):
    result = frame[0]
    for i, column in enumerate(frame[1:]):
        if point >> i & 1:
            result ^= column
    return result


def verify(native, root):
    rng = random.Random(1354)
    records, pairs = [], []
    for case in range(24):
        if case < 8:
            points = tuple(sorted(rng.sample(range(64), 54)))
        else:
            points = tuple(sorted(rng.sample(range(1 << 13), 54)))
        basis = [1 << i for i in range(13)]
        for _ in range(130):
            i, j = rng.sample(range(13), 2)
            basis[i] ^= basis[j]
        frame = [rng.randrange(1 << 13), *basis]
        require(binary_rank(basis) == 13, "The test frame is not invertible")
        transformed = tuple(sorted(apply(x, frame) for x in points))
        pairs.append((points, transformed, frame))
        records.extend((points, transformed))
    payload = "".join(" ".join(map(str, points))+"\n" for points in records)
    result = subprocess.run([str(native)], input=payload, text=True, capture_output=True, check=True)
    rows = [json.loads(line) for line in result.stdout.splitlines()]
    require(len(rows) == len(records), "The native colour output is incomplete")
    for points, colours in zip(records, rows, strict=True):
        require(colours == refine(points), "Native and reference point colours disagree")
    for i, (points, transformed, frame) in enumerate(pairs):
        target = dict(zip(transformed, rows[2*i+1], strict=True))
        require(all(rows[2*i][j] == target[apply(x, frame)] for j, x in enumerate(points)),
                "Point colours are not affine covariant")
    return dict(schema="point-colour-refinement-check-v1", status="pass", supports=len(records), affine_pairs=len(pairs),
        every_native_colour_equals_independent_reference=True, every_test_affine_image_preserves_colours=True,
        completeness_proof="Affine maps preserve difference-labelled neighbourhoods at every refinement step",
        dependencies={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in (
            "code/space_native/length54/m13_exact_quotient_kernel.cpp", "tests/native/m13_point_colours.cpp")})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.native.resolve(), Path(__file__).resolve().parents[1])
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
