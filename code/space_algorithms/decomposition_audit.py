"""Independent exact decomposition checks for binary generator matrices."""

from __future__ import annotations


def gf2_rank(vectors: list[int] | tuple[int, ...]) -> int:
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


def gf2_rref(
    rows: list[int] | tuple[int, ...], width: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    work = [int(row) for row in rows if row]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(width):
        candidate = next(
            (
                index
                for index in range(pivot_row, len(work))
                if work[index] & (1 << column)
            ),
            None,
        )
        if candidate is None:
            continue
        work[pivot_row], work[candidate] = work[candidate], work[pivot_row]
        for index in range(len(work)):
            if index != pivot_row and work[index] & (1 << column):
                work[index] ^= work[pivot_row]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(work):
            break
    return tuple(work[:pivot_row]), tuple(pivots)


def gf2_nullspace(
    rows: list[int] | tuple[int, ...], width: int
) -> tuple[int, ...]:
    reduced, pivots = gf2_rref(rows, width)
    pivot_set = set(pivots)
    basis = []
    for free_column in range(width):
        if free_column in pivot_set:
            continue
        vector = 1 << free_column
        for row, pivot in zip(reduced, pivots, strict=True):
            if row & (1 << free_column):
                vector |= 1 << pivot
        basis.append(vector)
    return tuple(basis)


def column_vectors(rows: tuple[str, ...]) -> tuple[int, ...]:
    return tuple(
        sum(int(row[column]) << row_index for row_index, row in enumerate(rows))
        for column in range(len(rows[0]))
    )


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, value: int) -> int:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            parent = self.parent[value]
            self.parent[value] = root
            value = parent
        return root

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def matroid_components(rows: tuple[str, ...]) -> tuple[tuple[int, ...], ...]:
    columns = column_vectors(rows)
    target_rank = gf2_rank(columns)
    basis_indices = []
    basis_columns = []
    for index, column in enumerate(columns):
        if gf2_rank((*basis_columns, column)) > len(basis_columns):
            basis_indices.append(index)
            basis_columns.append(column)
        if len(basis_columns) == target_rank:
            break
    if len(basis_columns) != target_rank:
        raise ValueError("Could not select a full column basis")

    echelon: dict[int, tuple[int, int]] = {}
    for basis_position, column in enumerate(basis_columns):
        reduced = column
        coefficients = 1 << basis_position
        for pivot in sorted(echelon, reverse=True):
            if reduced & (1 << pivot):
                vector, vector_coefficients = echelon[pivot]
                reduced ^= vector
                coefficients ^= vector_coefficients
        if not reduced:
            raise ValueError("Selected column basis is dependent")
        echelon[reduced.bit_length() - 1] = (reduced, coefficients)

    basis_index_set = set(basis_indices)
    union_find = UnionFind(len(columns))
    for index, column in enumerate(columns):
        if index in basis_index_set:
            continue
        reduced = column
        coefficients = 0
        for pivot in sorted(echelon, reverse=True):
            if reduced & (1 << pivot):
                vector, vector_coefficients = echelon[pivot]
                reduced ^= vector
                coefficients ^= vector_coefficients
        if reduced:
            raise ValueError("Column is outside the selected basis span")
        for basis_position, basis_index in enumerate(basis_indices):
            if coefficients & (1 << basis_position):
                union_find.union(index, basis_index)
    groups: dict[int, list[int]] = {}
    for coordinate in range(len(columns)):
        groups.setdefault(union_find.find(coordinate), []).append(coordinate)
    return tuple(sorted((tuple(group) for group in groups.values()), key=lambda x: x[0]))


def row_coordinate_int(row: str) -> int:
    return sum((bit == "1") << index for index, bit in enumerate(row))


def stabilizer_components(
    rows: tuple[str, ...]
) -> tuple[tuple[tuple[int, ...], ...], int]:
    width = len(rows[0])
    code_basis, _ = gf2_rref(tuple(row_coordinate_int(row) for row in rows), width)
    dual_basis = gf2_nullspace(code_basis, width)
    equations = tuple(codeword & dual for codeword in code_basis for dual in dual_basis)
    stabilizer_basis = gf2_nullspace(equations, width)
    groups: dict[tuple[int, ...], list[int]] = {}
    for coordinate in range(width):
        signature = tuple(
            int(bool(vector & (1 << coordinate))) for vector in stabilizer_basis
        )
        groups.setdefault(signature, []).append(coordinate)
    components = tuple(
        sorted((tuple(group) for group in groups.values()), key=lambda x: x[0])
    )
    return components, len(stabilizer_basis)


def decomposition_profile(rows: tuple[str, ...]) -> dict[str, object]:
    matroid = matroid_components(rows)
    stabilizer, stabilizer_dimension = stabilizer_components(rows)
    return {
        "methods": [
            "binary_matroid_fundamental_circuits",
            "multiplicative_stabilizer_algebra",
        ],
        "matroid_components": [list(component) for component in matroid],
        "stabilizer_components": [list(component) for component in stabilizer],
        "matroid_component_lengths": [len(component) for component in matroid],
        "stabilizer_component_lengths": [len(component) for component in stabilizer],
        "stabilizer_dimension": stabilizer_dimension,
        "methods_agree": matroid == stabilizer,
        "indecomposable": len(matroid) == 1
        and len(stabilizer) == 1
        and stabilizer_dimension == 1,
    }
