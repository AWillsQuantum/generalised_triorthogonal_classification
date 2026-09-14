"""Exact codeword-incidence canonicalisation with independent affine-map replay."""

import json
import subprocess

from space_codec import binary_rank, checked_support


def canonicalise_many(command, records):
    records = [(m, checked_support(points, m)) for m, points in records]
    for m, points in records:
        if not 1 <= m <= 17 or not m+1 <= len(points) <= 63:
            raise ValueError("Invalid codeword-graph support parameters")
        if binary_rank(x ^ points[0] for x in points) != m:
            raise ValueError("The support is not full rank")
    payload = "".join(f"{m} {len(points)} "+" ".join(map(str, points))+"\n" for m, points in records)
    output = subprocess.run([str(x) for x in command], input=payload, text=True, capture_output=True, check=True)
    results = [json.loads(line) for line in output.stdout.splitlines()]
    if len(results) != len(records):
        raise ValueError("The canonical support stream is incomplete")
    for (m, points), result in zip(records, results, strict=True):
        canonical = checked_support(result["canonical_points"], m, len(points))
        frame = result["frame"]
        if (result["m"] != m or result["n"] != len(points) or len(frame) != m+1
                or any(type(x) is not int or not 0 <= x < 1 << m for x in frame)
                or binary_rank(frame[1:]) != m):
            raise ValueError("Invalid canonical affine frame")
        image = []
        for point in canonical:
            value = frame[0]
            for bit, column in enumerate(frame[1:]):
                if point >> bit & 1:
                    value ^= column
            image.append(value)
        if sorted(image) != list(points):
            raise ValueError("The canonical affine frame does not reproduce the input")
    return results
