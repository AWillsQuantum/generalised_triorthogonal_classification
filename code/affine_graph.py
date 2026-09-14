"""Exact affine graph canonicalisation, with independent map replay."""

import json
from math import prod
import subprocess


def rank(vectors):
    rows = {}
    for value in vectors:
        if not isinstance(value, int) or value < 0:
            raise ValueError("Vectors must be nonnegative integers")
        while value:
            pivot = value.bit_length() - 1
            if pivot not in rows:
                rows[pivot] = value
                break
            value ^= rows[pivot]
    return len(rows)


def affine_images(frame, m):
    size = 1 << m
    if (len(frame) != m + 1 or any(type(x) is not int or not 0 <= x < size for x in frame)
            or rank(frame[1:]) != m):
        raise ValueError("Invalid affine frame")
    images = [frame[0]]
    for column in frame[1:]:
        images += [x ^ column for x in images]
    if len(set(images)) != size:
        raise ValueError("Affine frame is not a bijection")
    return images


def transform(mask, frame, m):
    return sum(1 << y for x, y in enumerate(affine_images(frame, m)) if mask >> x & 1)


def group_order(m):
    size = 1 << m
    return size * prod(size - (1 << i) for i in range(m))


def canonicalise_many(command, records):
    """Each record is (ambient dimension, truth-table mask); returns exact forms."""
    records = list(records)
    for m, mask in records:
        if not 2 <= m <= 10 or type(mask) is not int or not 0 <= mask < 1 << (1 << m):
            raise ValueError("Invalid affine graph input")
    payload = "".join(f"{m} {mask:x}\n" for m, mask in records)
    result = subprocess.run([str(x) for x in command], input=payload, text=True,
                            capture_output=True, check=True)
    rows = [json.loads(line) for line in result.stdout.splitlines()]
    if len(rows) != len(records):
        raise ValueError("Incomplete graph result stream")
    for (m, mask), row in zip(records, rows, strict=True):
        canonical = int(row["canonical_mask"], 16)
        order = row["stabilizer_order"]
        if row["m"] != m or not 0 <= canonical < 1 << (1 << m):
            raise ValueError("Incorrect canonical support dimension")
        if type(order) is not int or order <= 0 or group_order(m) % order:
            raise ValueError("Invalid exact stabiliser order")
        if transform(canonical, row["frame"], m) != mask:
            raise ValueError("Canonical affine frame does not reproduce the input")
        for frame in row["generators"]:
            if transform(mask, frame, m) != mask:
                raise ValueError("Affine generator does not stabilise the support")
    return rows
