"""Exact reconstruction primitives for inverse difference contractions."""

from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Iterable, Iterator

from quadratic_syndromes import degree_two_vector


class LinearBasis:
    def __init__(self) -> None:
        self.rows: dict[int, int] = {}

    def reduce(self, vector: int) -> int:
        value = int(vector)
        for pivot in sorted(self.rows, reverse=True):
            if value & (1 << pivot):
                value ^= self.rows[pivot]
        return value

    def add(self, vector: int) -> bool:
        value = self.reduce(vector)
        if not value:
            return False
        self.rows[value.bit_length() - 1] = value
        return True

    @property
    def rank(self) -> int:
        return len(self.rows)


class ColumnCombinationSolver:
    """Solve XOR combinations of columns and expose their relation space."""

    def __init__(self, columns: Iterable[int]) -> None:
        self.columns = tuple(int(column) for column in columns)
        self.rows: dict[int, tuple[int, int]] = {}
        relations = []
        for index, column in enumerate(self.columns):
            value = column
            coefficient_mask = 1 << index
            for pivot in sorted(self.rows, reverse=True):
                if value & (1 << pivot):
                    row, row_coefficients = self.rows[pivot]
                    value ^= row
                    coefficient_mask ^= row_coefficients
            if value:
                self.rows[value.bit_length() - 1] = (value, coefficient_mask)
            else:
                relations.append(coefficient_mask)
        self.relations = tuple(relations)

    @property
    def rank(self) -> int:
        return len(self.rows)

    def reduce(self, vector: int) -> int:
        value = int(vector)
        for pivot in sorted(self.rows, reverse=True):
            if value & (1 << pivot):
                value ^= self.rows[pivot][0]
        return value

    def solve(self, target: int) -> int | None:
        value = int(target)
        coefficients = 0
        for pivot in sorted(self.rows, reverse=True):
            if value & (1 << pivot):
                row, row_coefficients = self.rows[pivot]
                value ^= row
                coefficients ^= row_coefficients
        return None if value else coefficients

    def evaluate(self, coefficient_mask: int) -> int:
        value = 0
        remaining = int(coefficient_mask)
        while remaining:
            low_bit = remaining & -remaining
            value ^= self.columns[low_bit.bit_length() - 1]
            remaining ^= low_bit
        return value


def independent_masks(masks: Iterable[int]) -> tuple[int, ...]:
    basis = LinearBasis()
    result = []
    for mask in masks:
        if basis.add(int(mask)):
            result.append(int(mask))
    return tuple(result)


def affine_code_rows(points: tuple[int, ...], ambient_dimension: int) -> tuple[int, ...]:
    rows = [(1 << len(points)) - 1]
    for coordinate in range(ambient_dimension):
        rows.append(
            sum(
                ((point >> coordinate) & 1) << index
                for index, point in enumerate(points)
            )
        )
    return independent_masks(rows)


def iter_zero_xor_index_subsets(
    labels: Iterable[int], subset_size: int
) -> Iterator[tuple[int, ...]]:
    """Yield every fixed-size index subset with zero label XOR, once."""

    values = tuple(int(label) for label in labels)
    locations: dict[int, list[int]] = defaultdict(list)
    for index, label in enumerate(values):
        locations[label].append(index)

    if subset_size == 0:
        yield ()
        return
    if subset_size == 1:
        yield from ((index,) for index in locations.get(0, ()))
        return
    if subset_size == 2:
        for indices in locations.values():
            yield from combinations(indices, 2)
        return
    if subset_size == 3:
        for left in range(len(values) - 2):
            for middle in range(left + 1, len(values) - 1):
                matches = locations.get(values[left] ^ values[middle], ())
                start = bisect_right(matches, middle)
                for right in matches[start:]:
                    yield (left, middle, right)
        return
    if subset_size == 4:
        pairs_by_xor: dict[int, list[tuple[int, int]]] = defaultdict(list)
        for left in range(len(values) - 1):
            for right in range(left + 1, len(values)):
                pairs_by_xor[values[left] ^ values[right]].append((left, right))
        for pairs in pairs_by_xor.values():
            for pair_index, (left, middle) in enumerate(pairs[:-1]):
                for right, last in pairs[pair_index + 1 :]:
                    if middle < right:
                        yield (left, middle, right, last)
        return
    if subset_size == 5:
        pairs_by_xor: dict[int, list[tuple[int, int]]] = defaultdict(list)
        for left in range(len(values) - 1):
            for right in range(left + 1, len(values)):
                pairs_by_xor[values[left] ^ values[right]].append((left, right))
        for left in range(len(values) - 4):
            for middle in range(left + 1, len(values) - 3):
                for right in range(middle + 1, len(values) - 2):
                    matches = pairs_by_xor.get(
                        values[left] ^ values[middle] ^ values[right], ()
                    )
                    start = bisect_right(matches, (right, len(values)))
                    for fourth, last in matches[start:]:
                        yield (left, middle, right, fourth, last)
        return
    raise ValueError("The exact iterator supports subset sizes zero through five")


def count_zero_xor_index_subsets(labels: Iterable[int], subset_size: int) -> int:
    """Count fixed-size zero-XOR index subsets by exact dynamic programming."""

    if subset_size < 0:
        raise ValueError("The subset size must be nonnegative")
    counts: list[dict[int, int]] = [defaultdict(int) for _ in range(subset_size + 1)]
    counts[0][0] = 1
    seen = 0
    for original in labels:
        label = int(original)
        seen += 1
        for size in range(min(subset_size, seen), 0, -1):
            for value, count in tuple(counts[size - 1].items()):
                counts[size][value ^ label] += count
    return counts[subset_size].get(0, 0)


def count_zero_xor_index_subsets_walsh(
    labels: Iterable[int], subset_size: int
) -> int:
    """Count zero-XOR subsets by character averaging on the label span."""

    values = tuple(int(label) for label in labels)
    if not 0 <= subset_size <= len(values):
        return 0
    selection = LinearBasis()
    span_basis = []
    for value in values:
        if selection.add(value):
            span_basis.append(value)
    solver = ColumnCombinationSolver(span_basis)
    coordinates = tuple(solver.solve(value) for value in values)
    if any(value is None for value in coordinates):
        raise AssertionError("A label was not represented in its computed span")
    character_sum = 0
    for character in range(1 << len(span_basis)):
        negative = sum(
            ((character & int(coordinate)).bit_count() & 1)
            for coordinate in coordinates
        )
        positive = len(values) - negative
        coefficient = 0
        for selected_negative in range(
            max(0, subset_size - positive), min(subset_size, negative) + 1
        ):
            term = comb(negative, selected_negative) * comb(
                positive, subset_size - selected_negative
            )
            coefficient += -term if selected_negative & 1 else term
        character_sum += coefficient
    character_count = 1 << len(span_basis)
    if character_sum % character_count:
        raise AssertionError("Walsh character average is not integral")
    return character_sum // character_count


@dataclass(frozen=True)
class CoreExtensionContext:
    quotient_dimension: int
    core_points: tuple[int, ...]
    quadratic_columns: tuple[int, ...]
    external_points: tuple[int, ...]
    external_labels: tuple[int, ...]
    affine_rows: tuple[int, ...]
    lift_complement: tuple[int, ...]
    solver: ColumnCombinationSolver

    @classmethod
    def build(
        cls,
        quotient_dimension: int,
        core_points: Iterable[int],
        expected_core_code_dimension: int | None = None,
    ) -> "CoreExtensionContext":
        normalized = tuple(sorted(int(point) for point in core_points))
        if len(normalized) != len(set(normalized)):
            raise ValueError("The contraction core contains repeated points")
        if any(not 0 <= point < (1 << quotient_dimension) for point in normalized):
            raise ValueError("A contraction-core point lies outside the quotient")
        columns = tuple(
            degree_two_vector(point, quotient_dimension) for point in normalized
        )
        solver = ColumnCombinationSolver(columns)
        affine_rows = affine_code_rows(normalized, quotient_dimension)
        if expected_core_code_dimension is not None and len(affine_rows) != int(
            expected_core_code_dimension
        ):
            raise ValueError("The core affine-code dimension does not match its record")
        if any(solver.evaluate(row) for row in affine_rows):
            raise ValueError("An affine core row is not in the quadratic lift kernel")

        quotient_basis = LinearBasis()
        for row in affine_rows:
            if not quotient_basis.add(row):
                raise AssertionError("The retained affine rows are dependent")
        lift_complement = []
        for relation in solver.relations:
            if quotient_basis.add(relation):
                lift_complement.append(relation)
        expected_lift_dimension = len(normalized) - solver.rank - len(affine_rows)
        if len(lift_complement) != expected_lift_dimension:
            raise AssertionError("The lift-kernel quotient dimension is inconsistent")

        support_set = set(normalized)
        external_points = tuple(
            point
            for point in range(1 << quotient_dimension)
            if point not in support_set
        )
        external_labels = tuple(
            solver.reduce(degree_two_vector(point, quotient_dimension))
            for point in external_points
        )
        return cls(
            quotient_dimension=quotient_dimension,
            core_points=normalized,
            quadratic_columns=columns,
            external_points=external_points,
            external_labels=external_labels,
            affine_rows=affine_rows,
            lift_complement=tuple(lift_complement),
            solver=solver,
        )

    def lift_representatives(self, right_hand_side: int) -> tuple[int, ...]:
        particular = self.solver.solve(right_hand_side)
        if particular is None:
            return ()
        representatives = []
        for coefficients in range(1 << len(self.lift_complement)):
            lift = particular
            for index, relation in enumerate(self.lift_complement):
                if coefficients & (1 << index):
                    lift ^= relation
            if self.solver.evaluate(lift) != right_hand_side:
                raise AssertionError("A reconstructed lift does not solve the moment equation")
            representatives.append(lift)
        return tuple(representatives)

    def eligible_fiber_domain(
        self, multiplicity: int, kernel_dimension: int
    ) -> tuple[tuple[int, ...], tuple[int, ...]]:
        if kernel_dimension == 0:
            return self.external_points, self.external_labels
        if kernel_dimension == 1 and multiplicity == 4:
            selected = tuple(
                index
                for index, point in enumerate(self.external_points)
                if point & 1
            )
            return (
                tuple(self.external_points[index] for index in selected),
                tuple(self.external_labels[index] for index in selected),
            )
        raise ValueError("Unsupported nonzero-kernel contraction family")

    def iter_fiber_sets(
        self, multiplicity: int, kernel_dimension: int
    ) -> Iterator[tuple[int, ...]]:
        points, labels = self.eligible_fiber_domain(multiplicity, kernel_dimension)
        for indices in iter_zero_xor_index_subsets(labels, multiplicity):
            yield tuple(points[index] for index in indices)

    def reconstructed_support(
        self, fiber_points: tuple[int, ...], lift_mask: int
    ) -> tuple[int, ...]:
        singleton_points = tuple(
            (point << 1) | ((lift_mask >> index) & 1)
            for index, point in enumerate(self.core_points)
        )
        full_fibers = tuple(
            value
            for point in fiber_points
            for value in ((point << 1), (point << 1) | 1)
        )
        support = tuple(sorted((*singleton_points, *full_fibers)))
        if len(support) != len(self.core_points) + 2 * len(fiber_points):
            raise AssertionError("The reconstructed extension has the wrong length")
        if len(support) != len(set(support)):
            raise AssertionError("The reconstructed extension repeats a point")
        return support


def xor_quadratic_evaluations(
    points: Iterable[int], quotient_dimension: int
) -> int:
    result = 0
    for point in points:
        result ^= degree_two_vector(int(point), quotient_dimension)
    return result
