"""Pruned exact affine-equivalence checks for binary supports."""

from __future__ import annotations

from collections import Counter
from functools import lru_cache

from affine_tools import coordinate_map, ordered_basis
from support_invariants import (
    mask_from_points,
    points_from_mask,
    xor_difference_counts_from_points,
)


def fast_affinely_equivalent(source: tuple[int, ...], target: tuple[int, ...], m: int) -> bool:
    """Return whether two full-rank supports are affinely equivalent.

    This is still an exact backtracking search.  The speedup over
    ``affine_tools.affinely_equivalent`` is that every partial linear map must
    preserve derivative counts on the prefix span and point-local profiles for
    any support point whose image is already determined.
    """

    source = tuple(sorted(source))
    target = tuple(sorted(target))
    if len(source) != len(target) or len(set(source)) != len(source) or len(set(target)) != len(target):
        return False

    source_mask = mask_from_points(source)
    target_mask = mask_from_points(target)
    source_difference_counts = difference_counts(m, source_mask)
    target_difference_counts = difference_counts(m, target_mask)
    if sorted(source_difference_counts[1:]) != sorted(target_difference_counts[1:]):
        return False

    source_point_profiles = point_profiles(source, source_difference_counts)
    target_point_profiles = point_profiles(target, target_difference_counts)
    if sorted(source_point_profiles.values()) != sorted(target_point_profiles.values()):
        return False

    source_profile_counts = Counter(source_point_profiles.values())
    source_origin = min(source, key=lambda point: (source_profile_counts[source_point_profiles[point]], point))
    source_basis = ordered_source_basis(source, source_origin, m, source_difference_counts)
    source_coordinate = coordinate_map(source_basis)
    source_support_masks = tuple(sorted(source_coordinate[point ^ source_origin] for point in source))
    source_vectors_by_coordinate = tuple(vector_from_coordinates(coordinate, source_basis) for coordinate in range(1 << m))
    source_prefix_spans = tuple(
        tuple(range(1, 1 << depth))
        for depth in range(m + 1)
    )

    target_set = set(target)
    origin_profile = source_point_profiles[source_origin]
    target_origins = [
        point for point in target
        if target_point_profiles[point] == origin_profile
    ]

    for target_origin in target_origins:
        target_vectors = tuple(sorted(point ^ target_origin for point in target if point != target_origin))
        candidate_images_by_depth = tuple(
            tuple(
                vector for vector in target_vectors
                if target_difference_counts[vector] == source_difference_counts[source_basis[depth]]
            )
            for depth in range(m)
        )
        if any(not candidates for candidates in candidate_images_by_depth):
            continue
        if _search(
            m=m,
            source_support_masks=source_support_masks,
            source_vectors_by_coordinate=source_vectors_by_coordinate,
            source_difference_counts=source_difference_counts,
            source_point_profiles=source_point_profiles,
            source_origin=source_origin,
            source_prefix_spans=source_prefix_spans,
            target_set=target_set,
            target_origin=target_origin,
            target_difference_counts=target_difference_counts,
            target_point_profiles=target_point_profiles,
            candidate_images_by_depth=candidate_images_by_depth,
            chosen=(),
        ):
            return True
    return False


def ordered_source_basis(
    source: tuple[int, ...],
    source_origin: int,
    m: int,
    source_difference_counts: tuple[int, ...],
) -> tuple[int, ...]:
    basis = ordered_basis(tuple(point ^ source_origin for point in source), m)
    count_frequencies = Counter(source_difference_counts[vector] for vector in range(1, 1 << m))
    return tuple(sorted(basis, key=lambda vector: (count_frequencies[source_difference_counts[vector]], vector)))


def _search(
    *,
    m: int,
    source_support_masks: tuple[int, ...],
    source_vectors_by_coordinate: tuple[int, ...],
    source_difference_counts: tuple[int, ...],
    source_point_profiles: dict[int, tuple[int, ...]],
    source_origin: int,
    source_prefix_spans: tuple[tuple[int, ...], ...],
    target_set: set[int],
    target_origin: int,
    target_difference_counts: tuple[int, ...],
    target_point_profiles: dict[int, tuple[int, ...]],
    candidate_images_by_depth: tuple[tuple[int, ...], ...],
    chosen: tuple[int, ...],
) -> bool:
    depth = len(chosen)
    if depth == m:
        return _full_image_matches(
            source_support_masks=source_support_masks,
            target_set=target_set,
            target_origin=target_origin,
            chosen=chosen,
        )

    for candidate in candidate_images_by_depth[depth]:
        if candidate in chosen:
            continue
        if rank((*chosen, candidate)) <= depth:
            continue
        next_chosen = (*chosen, candidate)
        if not _prefix_span_invariants_hold(
            source_vectors_by_coordinate=source_vectors_by_coordinate,
            source_difference_counts=source_difference_counts,
            target_difference_counts=target_difference_counts,
            prefix_coordinates=source_prefix_spans[depth + 1],
            chosen=next_chosen,
        ):
            continue
        if not _support_prefix_constraints_hold(
            source_support_masks=source_support_masks,
            source_vectors_by_coordinate=source_vectors_by_coordinate,
            source_point_profiles=source_point_profiles,
            source_origin=source_origin,
            target_set=target_set,
            target_origin=target_origin,
            target_point_profiles=target_point_profiles,
            chosen=next_chosen,
        ):
            continue
        if _search(
            m=m,
            source_support_masks=source_support_masks,
            source_vectors_by_coordinate=source_vectors_by_coordinate,
            source_difference_counts=source_difference_counts,
            source_point_profiles=source_point_profiles,
            source_origin=source_origin,
            source_prefix_spans=source_prefix_spans,
            target_set=target_set,
            target_origin=target_origin,
            target_difference_counts=target_difference_counts,
            target_point_profiles=target_point_profiles,
            candidate_images_by_depth=candidate_images_by_depth,
            chosen=next_chosen,
        ):
            return True
    return False


def _prefix_span_invariants_hold(
    *,
    source_vectors_by_coordinate: tuple[int, ...],
    source_difference_counts: tuple[int, ...],
    target_difference_counts: tuple[int, ...],
    prefix_coordinates: tuple[int, ...],
    chosen: tuple[int, ...],
) -> bool:
    for coordinate in prefix_coordinates:
        source_vector = source_vectors_by_coordinate[coordinate]
        target_vector = apply_linear_images(coordinate, chosen)
        if source_difference_counts[source_vector] != target_difference_counts[target_vector]:
            return False
    return True


def _support_prefix_constraints_hold(
    *,
    source_support_masks: tuple[int, ...],
    source_vectors_by_coordinate: tuple[int, ...],
    source_point_profiles: dict[int, tuple[int, ...]],
    source_origin: int,
    target_set: set[int],
    target_origin: int,
    target_point_profiles: dict[int, tuple[int, ...]],
    chosen: tuple[int, ...],
) -> bool:
    prefix_mask = (1 << len(chosen)) - 1
    for source_mask in source_support_masks:
        if source_mask & ~prefix_mask:
            continue
        source_point = source_origin ^ source_vectors_by_coordinate[source_mask]
        target_point = target_origin ^ apply_linear_images(source_mask, chosen)
        if target_point not in target_set:
            return False
        if source_point_profiles[source_point] != target_point_profiles[target_point]:
            return False
    return True


def _full_image_matches(
    *,
    source_support_masks: tuple[int, ...],
    target_set: set[int],
    target_origin: int,
    chosen: tuple[int, ...],
) -> bool:
    images = {
        target_origin ^ apply_linear_images(source_mask, chosen)
        for source_mask in source_support_masks
    }
    return images == target_set


def apply_linear_images(point: int, images: tuple[int, ...]) -> int:
    image = 0
    for index, basis_image in enumerate(images):
        if (point >> index) & 1:
            image ^= basis_image
    return image


def vector_from_coordinates(coordinates: int, basis: tuple[int, ...]) -> int:
    vector = 0
    for index, basis_vector in enumerate(basis):
        if (coordinates >> index) & 1:
            vector ^= basis_vector
    return vector


@lru_cache(maxsize=None)
def translate_mask(m: int, support_mask: int, translation: int) -> int:
    translated = 0
    for point in range(1 << m):
        if (support_mask >> point) & 1:
            translated |= 1 << (point ^ translation)
    return translated


def difference_counts(m: int, support_mask: int) -> tuple[int, ...]:
    return xor_difference_counts_from_points(m, points_from_mask(m, support_mask))


def point_profiles(points: tuple[int, ...], counts: tuple[int, ...]) -> dict[int, tuple[int, ...]]:
    return {
        point: tuple(sorted(counts[point ^ other] for other in points if other != point))
        for point in points
    }


def rank(vectors: tuple[int, ...]) -> int:
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
