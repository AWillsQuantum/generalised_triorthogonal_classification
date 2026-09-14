#pragma once

#include "utsp/native_core.hpp"

#include <cstdint>
#include <span>
#include <vector>

namespace utsp {

struct MarkedCodeCanonicalForm {
  // Exact canonical serialization of (G_0, <G_0,U>).  Equality of these
  // words is equivalent to pointed-support equivalence of U.
  std::vector<Mask> key_words;
  std::vector<std::vector<Mask>> quotient_automorphism_generators;
  std::uint64_t automorphism_group_order = 1;
  std::uint64_t canonical_search_nodes = 0;
};

struct IsotropicOrbitLevel {
  std::uint32_t dimension = 0;
  std::uint64_t parent_orbits = 0;
  std::uint64_t extension_vector_orbits = 0;
  std::uint64_t filtered_degenerate_children = 0;
  std::uint64_t filtered_nonprimitive_children = 0;
  std::uint64_t filtered_target_children = 0;
  std::uint64_t canonical_parent_rejections = 0;
  std::uint64_t duplicate_children = 0;
  std::uint64_t canonical_search_nodes = 0;
  std::uint64_t support_automorphism_group_order = 1;
  std::uint64_t weighted_subspace_count = 0;
  bool used_complete_primitive_target_recognizer = false;
  std::vector<std::vector<Mask>> representatives;
  std::vector<std::vector<Mask>> canonical_keys;
  std::vector<std::uint64_t> orbit_sizes;
};

struct Q3PredecessorExistenceResult {
  bool found = false;
  std::uint64_t cs_seed_orbits = 0;
  std::uint64_t extension_vector_orbits = 0;
  std::uint64_t primitive_target_orbits = 0;
  std::vector<Mask> witness_quotient_rows;
};

struct Q5OutputProfilePresenceResult {
  // Traversal counters through the canonical q=3 and q=4 parent levels.
  // The final q=5 extensions are tensor-classified without marked-code
  // canonicalization, so representatives and orbit sizes remain empty.
  IsotropicOrbitLevel traversal;
  std::vector<std::vector<Mask>> canonical_output_keys;
};

struct RadicalStratifiedQ8Result {
  bool has_isotropic_subspace = false;
  std::uint64_t seed_orbits = 0;
  std::uint64_t parent_orbits = 0;
  std::uint64_t extension_vector_orbits = 0;
  std::uint64_t rejected_radical_dimension = 0;
  std::uint64_t duplicate_children = 0;
  std::uint64_t canonical_search_nodes = 0;
  std::vector<std::uint64_t> seed_orbits_by_intrinsic_dimension;
  std::vector<bool> sector_witness_found;
  std::uint32_t witness_intrinsic_dimension = 0;
  std::vector<Mask> witness_quotient_rows;
};

[[nodiscard]] bool marked_code_canonicalizer_available();

// Canonicalize a quotient subspace U while marking the fixed stabilizer
// subcode G_0.  quotient_rows may be any basis of U.
[[nodiscard]] MarkedCodeCanonicalForm canonicalize_marked_code(
    const LabelSpace& label_space,
    std::span<const Mask> quotient_rows);

// Exact level-wise canonical augmentation under the full coordinate
// automorphism group of the pointed support.  The returned rows are RREF
// bases in V_X/L_X coordinates.
[[nodiscard]] IsotropicOrbitLevel enumerate_isotropic_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers = 1);

// Retain only nondegenerate tensors at every level.  This is the hereditary
// branch; primitive rank sectors must be unioned separately for completeness.
[[nodiscard]] IsotropicOrbitLevel
enumerate_hereditary_nondegenerate_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers = 1);

// Keep all subspaces below minimum_nondegenerate_dimension, then retain only
// nondegenerate tensors at that dimension and above.
[[nodiscard]] IsotropicOrbitLevel
enumerate_dimension_filtered_nondegenerate_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t minimum_nondegenerate_dimension,
    std::uint32_t workers = 1);

// Enumerate exactly the nondegenerate q-spaces in the supplied candidate
// sector having no nondegenerate hyperplane restriction.
[[nodiscard]] IsotropicOrbitLevel enumerate_primitive_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers = 1);

// Enumerate subspaces whose tensor at every partial dimension is isomorphic
// to a restriction of one of the supplied final target signatures.
[[nodiscard]] IsotropicOrbitLevel enumerate_target_tensor_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t search_dimension,
    std::uint32_t target_dimension,
    std::span<const std::uint64_t> standard_target_signatures,
    std::uint32_t nondegenerate_seed_dimension = 0,
    std::uint32_t workers = 1);

// Exact q=5 predecessor gate.  The finite q=5 restriction census proves that
// every target contains D3_03 or one of the two primitive q=3 tensors.
[[nodiscard]] IsotropicOrbitLevel
enumerate_q3_predecessor_target_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers = 1);

// Decide the same q=3 predecessor condition while retaining only one witness.
// D3_03 is searched by extending complete D2_02 seed orbits; the primitive
// D3_05/D3_06 branch is evaluated only when no D3_03 witness exists.
[[nodiscard]] Q3PredecessorExistenceResult
find_q3_predecessor_target_subspace(
    const LabelSpace& label_space,
    std::uint32_t workers = 1);

// Complete q=5 search from the certified D3_03/D3_05 through-q4 chain
// cover, together with the two primitive q=5 tensor orbits.
[[nodiscard]] IsotropicOrbitLevel
enumerate_q5_q3_chain_cover_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers = 1);

// Enumerate the exact set of intrinsic q=5 output classes present on a
// support.  This uses the certified q=3-through-q=4 chain cover and the two
// primitive q=5 branches, but deliberately skips final marked-code
// canonicalization.  It is an existence/profile gate, not an orbit census.
[[nodiscard]] Q5OutputProfilePresenceResult
enumerate_q5_q3_chain_output_profile_presence(
    const LabelSpace& label_space,
    std::uint32_t workers = 1);

// Complete q=5 search from the certified seven-class q=4 hitting set,
// together with the two primitive q=5 tensor orbits.
[[nodiscard]] IsotropicOrbitLevel
enumerate_q5_q4_hitting_set_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers = 1);

// Complete nondegenerate union seeded by every primitive tensor dimension
// through q=8: hereditary q=1, then primitive q=2,3,5,7 branches.  Exact
// marked-code keys remove overlaps between extension branches.
[[nodiscard]] IsotropicOrbitLevel
enumerate_complete_nondegenerate_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers = 1);

// Enumerate complete nondegenerate (q-1)-orbits, deduplicate that union, and
// extend each retained parent once. This is complete at q when the primitive
// q-sector has been excluded independently.
[[nodiscard]] IsotropicOrbitLevel
enumerate_nondegenerate_successor_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers = 1);

// For q=5, enumerate primitive outputs through their canonical zero-tensor
// four-dimensional hyperplane.  The primitive-tensor classification proves
// that this sector is exact.
[[nodiscard]] IsotropicOrbitLevel
enumerate_q5_zero_hyperplane_primitive_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers = 1);

// Complete q=5 union retaining the established hereditary and primitive q=2
// and q=3 branches, while replacing only the primitive q=5 target branch by
// the exact zero-hyperplane sector above.
[[nodiscard]] IsotropicOrbitLevel
enumerate_complete_q5_zero_hyperplane_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers = 1);

// Exact q=8 existence search stratified by the radical dimension of the
// output tensor. This uses the certified absence of primitive nondegenerate
// q=8 tensors, so the k=8 sector extends complete q=7 seeds.
[[nodiscard]] RadicalStratifiedQ8Result
search_q8_isotropic_subspaces_by_tensor_radical(
    const LabelSpace& label_space,
    std::uint32_t workers = 1);

}  // namespace utsp
