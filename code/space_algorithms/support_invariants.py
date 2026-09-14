"""Fast affine-invariant signatures for binary point supports.

These signatures are filters, not canonical labels.  Equal affine supports have
equal signatures, but unequal supports can also have equal signatures.
"""

from __future__ import annotations

from functools import lru_cache


def hyperplane_intersection_profile(m: int, support_mask: int, support_size: int) -> tuple[int, ...]:
    """Return the affine-hyperplane intersection profile.

    For every nonzero linear functional a, the two affine hyperplanes a.x=0
    and a.x=1 meet the support in k and support_size-k points.  Translation can
    swap the two sides, so the affine-invariant datum for a is min(k,n-k).
    The linear part of an affine transformation permutes nonzero a, so the
    sorted tuple over all a != 0 is an affine invariant.
    """

    points = points_from_mask(m, support_mask)
    _check_support_size(points, support_size)
    return _hyperplane_profile_from_points(m, points)


def derivative_weight_profile(m: int, support_mask: int, support_size: int) -> tuple[int, ...]:
    """Return the sorted derivative-weight profile over nonzero translations.

    For v != 0, the Boolean derivative f(x)+f(x+v) has weight
    2*(|C|-|C cap (C+v)|).  Translations are permuted by affine linear parts.
    """

    points = points_from_mask(m, support_mask)
    _check_support_size(points, support_size)
    counts = xor_difference_counts_from_points(m, points)
    return _derivative_profile_from_counts(support_size, counts)


def local_difference_count_profile(m: int, support_mask: int) -> tuple[tuple[int, ...], ...]:
    """Return a point-local affine invariant derived from difference counts.

    For every nonzero translation v, let N(v)=|C cap (C+v)|.  For each point
    p in C, record the sorted values N(p+q) as q ranges over C-{p}.  The final
    sorted list over p in C is invariant under the full affine group.
    """

    points = points_from_mask(m, support_mask)
    counts = xor_difference_counts_from_points(m, points)
    return _local_profile_from_counts(points, counts)


def combined_affine_profiles(
    m: int, support_mask: int
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[tuple[int, ...], ...]]:
    """Compute all three affine filters while sharing exact XOR counts."""

    points = points_from_mask(m, support_mask)
    counts = xor_difference_counts_from_points(m, points)
    return (
        _hyperplane_profile_from_points(m, points),
        _derivative_profile_from_counts(len(points), counts),
        _local_profile_from_counts(points, counts),
    )


def combined_basic_signature(m: int, support_mask: int, support_size: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    points = points_from_mask(m, support_mask)
    _check_support_size(points, support_size)
    counts = xor_difference_counts_from_points(m, points)
    return (
        _hyperplane_profile_from_points(m, points),
        _derivative_profile_from_counts(support_size, counts),
    )


def points_from_mask(m: int, support_mask: int) -> tuple[int, ...]:
    if m < 0 or support_mask < 0 or support_mask >> (1 << m):
        raise ValueError("Support mask lies outside the stated ambient space")
    points = []
    remaining = support_mask
    while remaining:
        least_bit = remaining & -remaining
        points.append(least_bit.bit_length() - 1)
        remaining ^= least_bit
    return tuple(points)


def xor_difference_counts_from_points(
    m: int, points: tuple[int, ...] | list[int]
) -> tuple[int, ...]:
    """Return |C intersect (C+v)| for v != 0, with entry zero kept at zero.

    Every unordered pair {p,q} contributes the two ordered pairs (p,q) and
    (q,p) to the nonzero difference p+q.  This is exact over F_2 and costs
    O(|C|^2 + 2^m), rather than constructing all 2^m translated masks.
    """

    normalized = tuple(int(point) for point in points)
    if len(normalized) != len(set(normalized)) or any(
        point < 0 or point >= (1 << m) for point in normalized
    ):
        raise ValueError("Points must be distinct members of the ambient space")
    counts = [0] * (1 << m)
    for index, left in enumerate(normalized):
        for right in normalized[index + 1 :]:
            counts[left ^ right] += 2
    return tuple(counts)


def _check_support_size(points: tuple[int, ...], support_size: int) -> None:
    if len(points) != support_size:
        raise ValueError("support_size does not match support_mask")


def _hyperplane_profile_from_points(
    m: int, points: tuple[int, ...]
) -> tuple[int, ...]:
    support_size = len(points)
    spectrum = [0] * (1 << m)
    for point in points:
        spectrum[point] = 1
    step = 1
    while step < (1 << m):
        for block in range(0, 1 << m, 2 * step):
            for offset in range(step):
                left = spectrum[block + offset]
                right = spectrum[block + step + offset]
                spectrum[block + offset] = left + right
                spectrum[block + step + offset] = left - right
        step *= 2
    intersections = []
    for coefficient in spectrum[1:]:
        if (support_size + coefficient) % 2:
            raise ValueError("Walsh coefficient has impossible parity")
        count = (support_size + coefficient) // 2
        intersections.append(min(count, support_size - count))
    return tuple(sorted(intersections))


def _derivative_profile_from_counts(
    support_size: int, counts: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(sorted(2 * (support_size - count) for count in counts[1:]))


def _local_profile_from_counts(
    points: tuple[int, ...], counts: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    profiles = tuple(
        sorted(
            tuple(sorted(counts[point ^ other] for other in points if other != point))
            for point in points
        )
    )
    if any(len(profile) != len(points) - 1 for profile in profiles):
        raise ValueError("Unexpected local profile length")
    return profiles


@lru_cache(maxsize=None)
def hyperplane_zero_masks(m: int) -> tuple[int, ...]:
    masks = []
    for functional in range(1, 1 << m):
        mask = 0
        for point in range(1 << m):
            if ((functional & point).bit_count() % 2) == 0:
                mask |= 1 << point
        masks.append(mask)
    return tuple(masks)


def translate_mask(m: int, support_mask: int, translation: int) -> int:
    translated = 0
    for point in range(1 << m):
        if (support_mask >> point) & 1:
            translated |= 1 << (point ^ translation)
    return translated


def mask_from_points(points: tuple[int, ...] | list[int]) -> int:
    mask = 0
    for point in points:
        mask |= 1 << point
    return mask
