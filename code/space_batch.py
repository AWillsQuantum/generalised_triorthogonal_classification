"""Bounded-memory array operations for validating and encoding space supports."""

from itertools import combinations

import numpy as np

from space_codec import BITMAP, parameters, record_width


def checked_points(points, m, c):
    parameters(m, c)
    points = np.asarray(points)
    if points.ndim != 2 or points.shape[1] != c or points.dtype.kind not in "ui":
        raise ValueError("Expected a two-dimensional integer support array")
    if np.any(points < 0) or np.any(points >= 1 << m) or np.any(points[:, 1:] <= points[:, :-1]):
        raise ValueError("Support points are outside the domain or not strictly ordered")
    return points.astype(np.uint32, copy=False)


def ordered_point_bytes(points, m):
    width = (m + 7) // 8
    if width in (1, 2, 4):
        return points.astype(f"<u{width}", copy=False).tobytes()
    if width == 3:
        return points.astype("<u4").view(np.uint8).reshape(-1, 4)[:, :3].copy().tobytes()
    if width == 0 and np.all(points == 0):
        return b""
    raise ValueError("Unsupported point width")


def decode_text_supports(supports, matrices, m, c):
    """Check every supplied matrix against its binary-label support."""
    parameters(m, c)
    if len(supports) != len(matrices):
        raise ValueError("Support and matrix batch sizes differ")
    if any(len(s) != c or any(len(x) != m for x in s) for s in supports):
        raise ValueError("Incorrect binary point label lengths")
    if any(len(rows) != m + 1 or any(len(row) != c for row in rows) for rows in matrices):
        raise ValueError("Incorrect generator matrix shape")
    count = len(supports)
    sb = "".join("".join(s) for s in supports).encode("ascii")
    rb = "".join("".join(rows) for rows in matrices).encode("ascii")
    labels = np.frombuffer(sb, dtype=np.uint8).reshape(count, c, m)
    rows = np.frombuffer(rb, dtype=np.uint8).reshape(count, m + 1, c)
    if np.any((labels != 48) & (labels != 49)) or np.any(rows[:, 0, :] != 49):
        raise ValueError("Nonbinary point label or nonunital first row")
    if not np.array_equal(labels.transpose(0, 2, 1), rows[:, 1:, :]):
        raise ValueError("Generator matrix and support differ")
    weights = np.left_shift(np.uint32(1), np.arange(m - 1, -1, -1, dtype=np.int32))
    points = np.sum((labels - 48) * weights, axis=2, dtype=np.uint32)
    return checked_points(points, m, c)


def encode_bitmaps(points, m, c):
    points = checked_points(points, m, c)
    width = record_width(m, c, BITMAP)
    if len(points) * width > 128 * 1024 * 1024:
        raise ValueError("Bitmap batch exceeds the 128 MiB workspace limit")
    data = np.zeros((len(points), width), dtype=np.uint8)
    np.bitwise_or.at(data, (np.arange(len(points))[:, None], points >> 3),
                     np.left_shift(1, points & 7).astype(np.uint8))
    return data.tobytes()


def decode_bitmaps(data, m, c):
    width = record_width(m, c, BITMAP)
    if len(data) % width:
        raise ValueError("Truncated bitmap batch")
    count = len(data) // width
    if count * (1 << m) > 128 * 1024 * 1024:
        raise ValueError("Decoded bitmap batch exceeds the 128 MiB workspace limit")
    packed = np.frombuffer(data, dtype=np.uint8).reshape(count, width)
    bits = np.unpackbits(packed, axis=1, bitorder="little")
    if np.any(bits[:, 1 << m:]) or np.any(np.sum(bits, axis=1) != c):
        raise ValueError("Invalid bitmap weight or out-of-domain bit")
    return np.nonzero(bits)[1].astype(np.uint32).reshape(count, c)


def validate_batch(points, m, c):
    """Check affine spanning and all degree-at-most-three moments exactly."""
    points = checked_points(points, m, c)
    if not 0 < c <= 64 or c & 1:
        raise ValueError("This moment kernel requires positive even length at most 64")
    count = len(points)
    pivots = np.zeros((count, m), dtype=np.uint32)
    ranks = np.zeros(count, dtype=np.uint8)
    for column in range(1, c):
        value = points[:, column] ^ points[:, 0]
        for bit in range(m - 1, -1, -1):
            active = (value & (1 << bit)) != 0
            insert = active & (pivots[:, bit] == 0)
            pivots[insert, bit] = value[insert]
            ranks += insert
            value ^= np.where(active, pivots[:, bit], 0)
    if np.any(ranks != m):
        raise ValueError("A support does not affinely span its declared dimension")
    rows = np.zeros((count, m), dtype=np.uint64)
    weights = np.left_shift(np.uint64(1), np.arange(c, dtype=np.uint64))
    for bit in range(m):
        rows[:, bit] = np.sum(((points >> bit) & 1).astype(np.uint64) * weights, axis=1, dtype=np.uint64)
    for degree in (1, 2, 3):
        for indices in combinations(range(m), degree):
            overlap = rows[:, indices[0]].copy()
            for i in indices[1:]:
                overlap &= rows[:, i]
            if np.any(np.bitwise_count(overlap) & 1):
                raise ValueError("A support has an odd moment of degree at most three")
    return True
