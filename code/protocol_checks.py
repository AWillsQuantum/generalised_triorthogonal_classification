"""Independent exact checks of finite generalised triorthogonal witnesses."""

from collections import Counter
from itertools import combinations, combinations_with_replacement
from math import comb
import re

from space_codec import binary_rank


def check_matrix(rows, q):
    rows = tuple(rows)
    if not rows or not isinstance(q, int) or not 1 <= q <= len(rows):
        raise ValueError("Invalid logical dimension")
    n = len(rows[0])
    if not n or any(len(row) != n or set(row) - {"0", "1"} for row in rows):
        raise ValueError("Matrix must be nonempty, rectangular and binary")
    masks = tuple(int(row, 2) for row in rows)
    if binary_rank(masks) != len(rows):
        raise ValueError("Dependent matrix rows")
    columns = tuple(sum((row[j] == "1") << i for i, row in enumerate(rows))
                    for j in range(n))
    if len(set(columns)) != n:
        raise ValueError("Repeated complete columns")
    for degree in (1, 2, 3):
        for indices in combinations(range(len(rows)), degree):
            if indices[-1] < q:
                continue
            overlap = (1 << n) - 1
            for i in indices:
                overlap &= masks[i]
            if overlap.bit_count() & 1:
                raise ValueError("Odd overlap containing a stabiliser row")
    return masks, columns


def exact_distance_and_coefficient(columns, q):
    """Minimum-weight syndrome dynamic programming, counting distinct subsets."""
    columns = tuple(columns)
    if q < 1 or any(not isinstance(v, int) or v < 0 for v in columns):
        raise ValueError("Invalid columns or logical dimension")
    dimension = max(q, max(columns, default=0).bit_length())
    if dimension > 22:
        raise ValueError("This reference syndrome DP is limited to 22 matrix rows")
    unreachable = len(columns) + 1
    distances = [unreachable] * (1 << dimension)
    counts = [0] * (1 << dimension)
    distances[0], counts[0] = 0, 1
    # States paired by XOR with this column can be updated simultaneously.
    for column in columns:
        if not column:
            continue
        pivot = 1 << (column.bit_length() - 1)
        for start in range(0, len(distances), 2 * pivot):
            for left in range(start, start + pivot):
                right = left ^ column
                dl, dr = distances[left], distances[right]
                cl, cr = counts[left], counts[right]
                if dr + 1 < dl:
                    distances[left], counts[left] = dr + 1, cr
                elif dr + 1 == dl:
                    counts[left] += cr
                if dl + 1 < dr:
                    distances[right], counts[right] = dl + 1, cl
                elif dl + 1 == dr:
                    counts[right] += cl
    distance = min(distances[1:1 << q])
    if distance == unreachable:
        return None, 0
    return distance, sum(counts[value] for value in range(1, 1 << q)
                         if distances[value] == distance)


def distance_up_to(columns, q, maximum=5):
    """Exact small-weight search by disjoint subset matching.

    A result of (None, 0) means only that the distance exceeds the cutoff
    or that no logical error exists. No upper bound is asserted in that case.
    """
    columns = tuple(columns)
    if (type(q) is not int or q < 1 or type(maximum) is not int or maximum < 1
            or any(type(value) is not int or value < 0 for value in columns)):
        raise ValueError("Invalid small-weight distance parameters")
    logical_mask = (1 << q)-1
    subsets = {0: [(0, 0)]}

    def of_weight(weight):
        if weight not in subsets:
            values = []
            for indices in combinations(range(len(columns)), weight):
                syndrome, support = 0, 0
                for index in indices:
                    syndrome ^= columns[index]
                    support |= 1 << index
                values.append((syndrome, support))
            subsets[weight] = values
        return subsets[weight]

    for weight in range(1, min(maximum, len(columns))+1):
        left_weight = weight//2
        buckets = {}
        for syndrome, support in of_weight(left_weight):
            buckets.setdefault(syndrome >> q, []).append((syndrome, support))
        count = 0
        for syndrome, support in of_weight(weight-left_weight):
            for other_syndrome, other_support in buckets.get(syndrome >> q, ()):
                if not support & other_support and (syndrome ^ other_syndrome) & logical_mask:
                    count += 1
        multiplicity = comb(weight, left_weight)
        if count % multiplicity:
            raise AssertionError("Subset matching has inconsistent multiplicity")
        if count:
            return weight, count//multiplicity
    return None, 0


class LogicalTensor:
    """A symmetric trilinear tensor represented by binary logical rows."""

    def __init__(self, logical_rows):
        self.rows = tuple(logical_rows)
        self.q = len(self.rows)
        self.words = [0] * (1 << self.q)
        for vector in range(1, len(self.words)):
            bit = vector & -vector
            self.words[vector] = self.words[vector ^ bit] ^ self.rows[bit.bit_length() - 1]
        self.contraction_invariants = tuple(
            (self.value(v, v, v), binary_rank(
                sum(self.value(v, 1 << i, 1 << j) << j for j in range(self.q))
                for i in range(self.q)))
            for v in range(1 << self.q)
        )

    def value(self, left, middle, right):
        return (self.words[left] & self.words[middle] & self.words[right]).bit_count() & 1

    def intrinsic(self):
        return all(rank != 0 for _, rank in self.contraction_invariants[1:])

    def key_words(self):
        bits = []
        for right in range(self.q):
            e = 1 << right
            bits.append(self.value(e, e, e))
            bits.extend(self.value(1 << left, e, e) for left in range(right))
            bits.extend(self.value(1 << a, 1 << b, e)
                        for a, b in combinations(range(right), 2))
        result = [0] * ((len(bits) + 63) // 64)
        for i, bit in enumerate(bits):
            result[i // 64] |= bit << (i % 64)
        return [f"0x{word:016x}" for word in result]


def gate_tensor(q, gate):
    """Use the one-, three- and seven-parity decompositions of T, CS and CCZ."""
    matches = re.findall(r"(CCZ|CS|T)([1-9]+)", gate)
    if not matches or "".join(a + b for a, b in matches) != gate or not 1 <= q <= 9:
        raise ValueError("Invalid compact gate notation")
    columns = []
    for name, digits in matches:
        indices = tuple(int(digit) - 1 for digit in digits)
        if (len(indices) != {"T": 1, "CS": 2, "CCZ": 3}[name]
                or tuple(sorted(set(indices))) != indices or indices[-1] >= q):
            raise ValueError("Invalid gate indices")
        for mask in range(1, 1 << len(indices)):
            columns.append(sum(1 << index for i, index in enumerate(indices) if mask & (1 << i)))
    return LogicalTensor(tuple(sum(((value >> i) & 1) << j for j, value in enumerate(columns))
                               for i in range(q)))


def verify_basis(source, target, basis):
    if (source.q != target.q or len(basis) != source.q
            or any(not isinstance(v, int) or not 0 < v < 1 << source.q for v in basis)
            or binary_rank(basis) != source.q):
        raise ValueError("Invalid output equivalence basis")
    for i, j, k in combinations_with_replacement(range(source.q), 3):
        if source.value(basis[i], basis[j], basis[k]) != target.value(1 << i, 1 << j, 1 << k):
            raise ValueError("Output equivalence does not replay")
    return True


def find_equivalence_basis(source, target):
    """Exact small-q isomorphism search; returns target basis in source coordinates."""
    if source.q != target.q:
        return None
    if Counter(source.contraction_invariants) != Counter(target.contraction_invariants):
        return None
    q = source.q
    candidates = {i: [v for v in range(1, 1 << q)
                      if source.contraction_invariants[v] == target.contraction_invariants[1 << i]]
                  for i in range(q)}
    order = sorted(range(q), key=lambda i: len(candidates[i]))
    chosen = {}

    def extend(depth, span):
        if depth == q:
            return tuple(chosen[i] for i in range(q))
        i = order[depth]
        used = set(span.values())
        for vector in candidates[i]:
            if vector in used:
                continue
            if any(source.contraction_invariants[vector ^ value] != target.contraction_invariants[(1 << i) ^ key]
                   for key, value in span.items()):
                continue
            if any(source.value(vector, vector, value) != target.value(1 << i, 1 << i, 1 << j)
                   for j, value in chosen.items()):
                continue
            if any(source.value(vector, a, b) != target.value(1 << i, 1 << j, 1 << k)
                   for (j, a), (k, b) in combinations(chosen.items(), 2)):
                continue
            chosen[i] = vector
            enlarged = span | {key ^ (1 << i): value ^ vector for key, value in span.items()}
            result = extend(depth + 1, enlarged)
            if result is not None:
                return result
            del chosen[i]
        return None

    basis = extend(0, {0: 0})
    if basis is not None:
        verify_basis(source, target, basis)
    return basis
