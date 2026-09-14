"""Degree-two evaluation spans and contraction syndromes over F_2."""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from typing import Iterable


def degree_two_vector(point: int, m: int) -> int:
    coordinates = [
        (point >> (m - 1 - variable)) & 1 for variable in range(m)
    ]
    value = 1
    bit = 1
    for coordinate in coordinates:
        bit <<= 1
        if coordinate:
            value |= bit
    for left, right in combinations(range(m), 2):
        bit <<= 1
        if coordinates[left] & coordinates[right]:
            value |= bit
    return value


def degree_two_monomial_count(m: int) -> int:
    return 1 + m + m * (m - 1) // 2


class ReducedBasis:
    def __init__(self, vectors: Iterable[int] = ()) -> None:
        self._rows: dict[int, int] = {}
        for vector in vectors:
            self.add(vector)

    def reduce(self, vector: int) -> int:
        value = int(vector)
        for pivot in sorted(self._rows, reverse=True):
            if (value >> pivot) & 1:
                value ^= self._rows[pivot]
        return value

    def add(self, vector: int) -> bool:
        value = self.reduce(vector)
        if not value:
            return False
        pivot = value.bit_length() - 1
        for other_pivot, row in tuple(self._rows.items()):
            if (row >> pivot) & 1:
                self._rows[other_pivot] = row ^ value
        self._rows[pivot] = value
        return True

    @property
    def rank(self) -> int:
        return len(self._rows)


def core_profile(m: int, support: Iterable[int]) -> dict[str, object]:
    points = tuple(int(point) for point in support)
    basis = ReducedBasis(degree_two_vector(point, m) for point in points)
    syndromes = {
        point: basis.reduce(degree_two_vector(point, m))
        for point in range(1 << m)
    }
    support_set = set(points)
    complement_distribution = Counter(
        syndrome for point, syndrome in syndromes.items() if point not in support_set
    )
    return {
        "ambient_dimension": m,
        "monomial_count": degree_two_monomial_count(m),
        "evaluation_rank": basis.rank,
        "evaluation_kernel_dimension": len(points) - basis.rank,
        "syndromes": syndromes,
        "quadratic_closure_size": sum(syndrome == 0 for syndrome in syndromes.values()),
        "complement_syndrome_distribution": dict(sorted(complement_distribution.items())),
    }


def count_fixed_size_zero_xor_subsets(labels: Iterable[int], size: int) -> int:
    states: list[dict[int, int]] = [dict() for _ in range(size + 1)]
    states[0][0] = 1
    processed = 0
    for label in labels:
        processed += 1
        for selected in range(min(size, processed), 0, -1):
            destination = states[selected]
            for old_xor, count in states[selected - 1].items():
                new_xor = old_xor ^ int(label)
                destination[new_xor] = destination.get(new_xor, 0) + count
    return states[size].get(0, 0)
