"""Classify small binary kernel-pattern codes by exact column multiplicities."""

from __future__ import annotations

from functools import lru_cache
from itertools import product


def gf2_rank(vectors: tuple[int, ...]) -> int:
    basis: dict[int, int] = {}
    for vector in vectors:
        value = int(vector)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for remainder in compositions(total - first, parts - 1):
            yield (first, *remainder)


def apply_linear(vector: int, basis_images: tuple[int, ...]) -> int:
    image = 0
    for bit, basis_image in enumerate(basis_images):
        if (vector >> bit) & 1:
            image ^= basis_image
    return image


@lru_cache(maxsize=None)
def general_linear_maps(dimension: int) -> tuple[tuple[int, ...], ...]:
    if dimension == 0:
        return ((),)
    nonzero = range(1, 1 << dimension)
    return tuple(
        images
        for images in product(nonzero, repeat=dimension)
        if gf2_rank(tuple(images)) == dimension
    )


def transform_counts(
    counts: tuple[int, ...], basis_images: tuple[int, ...]
) -> tuple[int, ...]:
    transformed = [0] * len(counts)
    for vector, count in enumerate(counts):
        transformed[apply_linear(vector, basis_images)] += count
    return tuple(transformed)


def canonical_counts(counts: tuple[int, ...], dimension: int) -> tuple[int, ...]:
    return min(
        transform_counts(counts, linear_map)
        for linear_map in general_linear_maps(dimension)
    )


def codeword_weights(counts: tuple[int, ...], dimension: int) -> tuple[int, ...]:
    return tuple(
        sum(count for vector, count in enumerate(counts) if (functional & vector).bit_count() & 1)
        for functional in range(1, 1 << dimension)
    )


def is_self_orthogonal(counts: tuple[int, ...], dimension: int) -> bool:
    for left in range(dimension):
        for right in range(left, dimension):
            overlap = sum(
                count
                for vector, count in enumerate(counts)
                if ((vector >> left) & 1) and ((vector >> right) & 1)
            )
            if overlap % 2:
                return False
    return True


def classify_kernel_patterns(length: int, dimension: int) -> tuple[dict[str, object], ...]:
    if dimension < 0 or dimension > length:
        raise ValueError("Invalid kernel-pattern dimensions")
    if dimension == 0:
        return (
            {
                "canonical_column_multiplicities": [length],
                "nonzero_codeword_weights": [],
                "minimum_distance": None,
                "generator_rows": [],
            },
        )

    canonical = {}
    for counts in compositions(length, 1 << dimension):
        present = tuple(vector for vector, count in enumerate(counts) if count)
        if gf2_rank(present) != dimension:
            continue
        weights = codeword_weights(counts, dimension)
        if min(weights) < 4 or not is_self_orthogonal(counts, dimension):
            continue
        representative = canonical_counts(counts, dimension)
        canonical[representative] = codeword_weights(representative, dimension)

    result = []
    for counts, weights in sorted(canonical.items()):
        columns = tuple(
            vector for vector, multiplicity in enumerate(counts) for _ in range(multiplicity)
        )
        rows = [
            "".join(str((column >> bit) & 1) for column in columns)
            for bit in range(dimension)
        ]
        result.append(
            {
                "canonical_column_multiplicities": list(counts),
                "nonzero_codeword_weights": list(weights),
                "minimum_distance": min(weights),
                "generator_rows": rows,
            }
        )
    return tuple(result)
