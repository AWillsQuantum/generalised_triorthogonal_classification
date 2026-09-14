"""Lossless support encodings and independent expansion to binary matrices."""

from dataclasses import dataclass
from itertools import combinations
from math import comb
import struct


MAGIC = b"UTSPACE1"
HEADER = struct.Struct("<8sBBHQQI")
COLEX = 1
BITMAP = 2


def parameters(m, c):
    if not isinstance(m, int) or not 0 <= m <= 26:
        raise ValueError("Ambient dimension must be an integer in [0,26]")
    if not isinstance(c, int) or not 0 <= c <= min(65535, 1 << m):
        raise ValueError("Invalid support cardinality")
    return 1 << m


def record_width(m, c, codec=COLEX):
    size = parameters(m, c)
    if codec == COLEX:
        return ((comb(size, c) - 1).bit_length() + 7) // 8
    if codec == BITMAP:
        return (size + 7) // 8
    raise ValueError("Unknown support codec")


def checked_support(points, m, c=None):
    points = tuple(points)
    size = parameters(m, len(points) if c is None else c)
    if c is not None and len(points) != c:
        raise ValueError("Incorrect support cardinality")
    if any(not isinstance(x, int) or not 0 <= x < size for x in points):
        raise ValueError("Point outside the ambient space")
    if any(a >= b for a, b in zip(points, points[1:])):
        raise ValueError("Points must be distinct and strictly increasing")
    return points


def rank_support(points, m):
    points = checked_support(points, m)
    return sum(comb(x, i) for i, x in enumerate(points, 1))


def unrank_support(rank, m, c):
    size = parameters(m, c)
    if not isinstance(rank, int) or not 0 <= rank < comb(size, c):
        raise ValueError("Rank outside the fixed-weight domain")
    points = [0] * c
    upper = size - 1
    for i in range(c, 0, -1):
        low, high = i - 1, upper
        while low < high:
            middle = (low + high + 1) // 2
            if comb(middle, i) <= rank:
                low = middle
            else:
                high = middle - 1
        points[i - 1] = low
        rank -= comb(low, i)
        upper = low - 1
    if rank:
        raise ValueError("Nonzero remainder after rank decoding")
    return tuple(points)


def encode_support(points, m, codec=COLEX):
    points = checked_support(points, m)
    width = record_width(m, len(points), codec)
    value = rank_support(points, m) if codec == COLEX else sum(1 << x for x in points)
    return value.to_bytes(width, "little")


def decode_support(data, m, c, codec=COLEX):
    if len(data) != record_width(m, c, codec):
        raise ValueError("Incorrect record width")
    value = int.from_bytes(data, "little")
    if codec == COLEX:
        return unrank_support(value, m, c)
    if value >> (1 << m) or value.bit_count() != c:
        raise ValueError("Invalid bitmap domain or weight")
    points = []
    while value:
        bit = value & -value
        points.append(bit.bit_length() - 1)
        value ^= bit
    return tuple(points)


def generator_rows(points, m):
    points = checked_support(points, m)
    return ["1" * len(points)] + [
        "".join(str((x >> bit) & 1) for x in points)
        for bit in range(m - 1, -1, -1)
    ]


def binary_rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
    return len(pivots)


def validate_unital_support(points, m):
    points = checked_support(points, m)
    rows = generator_rows(points, m)
    row_masks = [int(row, 2) if row else 0 for row in rows]
    if binary_rank(row_masks) != m + 1:
        raise ValueError("Support does not affinely span its ambient space")
    for degree in (1, 2, 3):
        for selected in combinations(row_masks, degree):
            overlap = (1 << len(points)) - 1
            for row in selected:
                overlap &= row
            if overlap.bit_count() & 1:
                raise ValueError("Odd overlap of at most three rows")
    return True


def indicator_anf(points, m):
    """Return monomial masks; m remains explicit even for unused variables."""
    points = checked_support(points, m)
    if m > 20:
        raise ValueError("Explicit truth-table expansion is limited to m <= 20")
    coefficients = bytearray(1 << m)
    for point in points:
        coefficients[point] = 1
    for bit in range(m):
        step = 1 << bit
        for block in range(0, 1 << m, 2 * step):
            for offset in range(step):
                coefficients[block + step + offset] ^= coefficients[block + offset]
    return tuple(mask for mask, value in enumerate(coefficients) if value)


def polynomial_text(monomials, m):
    terms = []
    for mask in monomials:
        if not isinstance(mask, int) or not 0 <= mask < 1 << m:
            raise ValueError("Monomial outside the ambient variables")
        terms.append("*".join(f"x_{i+1}" for i in range(m)
                              if mask & (1 << (m-1-i))) or "1")
    return " + ".join(terms) or "0"


@dataclass(frozen=True)
class ShardHeader:
    codec: int
    m: int
    c: int
    first_index: int
    count: int

    @property
    def width(self):
        return record_width(self.m, self.c, self.codec)

    def pack(self):
        if not 0 <= self.first_index < 1 << 64 or not 0 <= self.count < 1 << 64:
            raise ValueError("Invalid shard index range")
        if self.first_index + self.count > 1 << 64:
            raise ValueError("Shard index range overflows uint64")
        return HEADER.pack(MAGIC, self.codec, self.m, self.c,
                           self.first_index, self.count, self.width)


def read_exact(stream, count):
    chunks = []
    remaining = count
    while remaining:
        part = stream.read(remaining)
        if not part:
            raise ValueError("Truncated support shard")
        chunks.append(part)
        remaining -= len(part)
    return b"".join(chunks)


def read_header(stream):
    magic, codec, m, c, first, count, width = HEADER.unpack(read_exact(stream, HEADER.size))
    if magic != MAGIC:
        raise ValueError("Incorrect support shard magic")
    header = ShardHeader(codec, m, c, first, count)
    if width != header.width or header.pack() != HEADER.pack(magic, codec, m, c, first, count, width):
        raise ValueError("Inconsistent support shard header")
    return header


def read_records(stream, header):
    for _ in range(header.count):
        yield decode_support(read_exact(stream, header.width), header.m, header.c, header.codec)
    if stream.read(1):
        raise ValueError("Trailing bytes after support records")


def write_shard(stream, header, supports):
    stream.write(header.pack())
    seen = 0
    for points in supports:
        if seen >= header.count:
            raise ValueError("Too many support records")
        points = checked_support(points, header.m, header.c)
        stream.write(encode_support(points, header.m, header.codec))
        seen += 1
    if seen != header.count:
        raise ValueError("Too few support records")
