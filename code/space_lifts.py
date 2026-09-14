"""Reference inverse-contraction equations and affine-shear quotient."""

from itertools import combinations

from space_codec import binary_rank, checked_support


def monomials(m, degree):
    return tuple(sum(1 << i for i in indices)
                 for size in range(degree + 1) for indices in combinations(range(m), size))


def evaluation_rows(points, m, degree):
    points = checked_support(points, m)
    return tuple(sum(((x & monomial) == monomial) << i for i, x in enumerate(points))
                 for monomial in monomials(m, degree))


def reduce_vector(value, pivots):
    for bit in sorted(pivots, reverse=True):
        if value & (1 << bit):
            value ^= pivots[bit]
    return value


def linear_solutions(rows, rhs, width):
    """Return one solution and a kernel basis, or None for an inconsistent system."""
    rows, rhs = tuple(rows), tuple(rhs)
    if (len(rows) != len(rhs) or width < 0
            or any(not 0 <= row < 1 << width for row in rows)
            or any(bit not in (0, 1) for bit in rhs)):
        raise ValueError("Invalid binary system")
    pivots = {}
    for row, bit in zip(rows, rhs):
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row, bit
                break
            other, value = pivots[pivot]
            row ^= other
            bit ^= value
        if not row and bit:
            return None

    def solve(value, homogeneous):
        for pivot, (row, bit) in sorted(pivots.items()):
            if ((value & row).bit_count() & 1) ^ (0 if homogeneous else bit):
                value ^= 1 << pivot
        return value

    particular = solve(0, False)
    kernel = tuple(solve(1 << i, True) for i in range(width) if i not in pivots)
    return particular, kernel


def lift_family(core, fibres, quotient_dimension):
    """All graph labels over fixed complete fibres, modulo affine shears.

    The core may lie in a proper affine subspace of the quotient domain.
    Full affine rank of the reconstructed target remains a separate test.
    """
    core = checked_support(core, quotient_dimension)
    fibres = checked_support(fibres, quotient_dimension)
    if set(core).intersection(fibres):
        raise ValueError("Singleton and complete fibres overlap")
    if any(row.bit_count() & 1 for row in evaluation_rows(core, quotient_dimension, 3)):
        raise ValueError("Core fails a moment of degree at most three")
    equations = evaluation_rows(core, quotient_dimension, 2)
    rhs = [row.bit_count() & 1 for row in evaluation_rows(fibres, quotient_dimension, 2)]
    result = linear_solutions(equations, rhs, len(core))
    if result is None:
        return None
    particular, kernel = result
    affine_pivots = {}
    for row in evaluation_rows(core, quotient_dimension, 1):
        value = reduce_vector(row, affine_pivots)
        if value:
            affine_pivots[value.bit_length() - 1] = value
    pivots = dict(affine_pivots)
    quotient = []
    for row in kernel:
        value = reduce_vector(row, pivots)
        if value:
            pivots[value.bit_length() - 1] = value
            quotient.append(value)
    return reduce_vector(particular, affine_pivots), tuple(quotient)


def reconstruct_lift(core, fibres, labels, quotient_dimension):
    core = checked_support(core, quotient_dimension)
    fibres = checked_support(fibres, quotient_dimension)
    if set(core).intersection(fibres) or not 0 <= labels < 1 << len(core):
        raise ValueError("Invalid inverse-contraction data")
    return tuple(sorted([(x << 1) | ((labels >> i) & 1) for i, x in enumerate(core)]
                        + [2*x+t for x in fibres for t in (0, 1)]))


def iter_lifts(core, fibres, quotient_dimension, full_rank=True):
    core, fibres = tuple(core), tuple(fibres)
    family = lift_family(core, fibres, quotient_dimension)
    if family is None:
        return
    labels, basis = family
    for index in range(1 << len(basis)):
        if index:
            labels ^= basis[(index & -index).bit_length() - 1]
        support = reconstruct_lift(core, fibres, labels, quotient_dimension)
        if full_rank and (not support or binary_rank(x ^ support[0] for x in support) != quotient_dimension + 1):
            continue
        yield support
