#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace utsp {

struct DirectTensorCanonicalForm {
  std::uint32_t logical_qubits = 0;
  std::uint64_t canonical_key = 0;
  std::vector<std::uint64_t> canonical_key_words;
  std::vector<std::uint32_t> new_basis_in_old_coordinates;
  std::uint64_t complete_bases_evaluated = 0;
  std::uint64_t prefix_branches_pruned = 0;
};

struct FiveQubitTensorOrbit {
  std::uint32_t canonical_key = 0;
  std::uint32_t canonical_standard_signature = 0;
  std::uint64_t orbit_size = 0;
  std::uint64_t stabilizer_size = 0;
  std::uint32_t radical_dimension = 0;
};

struct FiveQubitTensorOrbitCensus {
  std::uint64_t tensor_count = 0;
  std::uint64_t generated_linear_group_order = 0;
  std::uint64_t expected_linear_group_order = 0;
  std::uint64_t generator_transitions = 0;
  std::uint64_t direct_canonical_checks = 0;
  std::uint64_t direct_canonical_mismatches = 0;
  std::array<std::uint64_t, 6> orbit_counts_by_radical_dimension{};
  std::array<std::uint64_t, 6> tensor_counts_by_radical_dimension{};
  std::vector<FiveQubitTensorOrbit> orbits;
};

[[nodiscard]] std::uint32_t cubic_tensor_term_count(
    std::uint32_t logical_qubits);

[[nodiscard]] std::uint32_t cubic_tensor_repeated_term_count(
    std::uint32_t logical_qubits);

[[nodiscard]] std::uint32_t cubic_tensor_alternating_term_count(
    std::uint32_t logical_qubits);

[[nodiscard]] bool cubic_tensor_value(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature,
    std::uint32_t left,
    std::uint32_t middle,
    std::uint32_t right);

// Pull a tensor back along the basis whose entries are vectors in the old
// coordinates.  The result uses the standard singles, pairs, triples order.
[[nodiscard]] std::uint64_t transform_cubic_tensor_signature(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature,
    std::span<const std::uint32_t> new_basis_in_old_coordinates);

// Restrict a tensor on F_2^source_logical_qubits to the independent vectors
// in subspace_basis, returning a standard signature of the smaller tensor.
[[nodiscard]] std::uint64_t restrict_cubic_tensor_signature(
    std::uint32_t source_logical_qubits,
    std::uint64_t standard_tensor_signature,
    std::span<const std::uint32_t> subspace_basis);

[[nodiscard]] std::vector<std::uint32_t>
cubic_tensor_radical_basis_from_signature(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature);

[[nodiscard]] std::uint64_t
nondegenerate_tensor_hyperplane_count_from_signature(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature);

[[nodiscard]] DirectTensorCanonicalForm canonicalize_cubic_tensor_direct(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature);

[[nodiscard]] DirectTensorCanonicalForm canonicalize_cubic_tensor_from_labels(
    std::span<const std::uint64_t> label_rows);

[[nodiscard]] std::uint64_t interleaved_tensor_key(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature,
    const std::vector<std::uint32_t>& new_basis_in_old_coordinates);

[[nodiscard]] std::vector<std::uint64_t> interleaved_tensor_words_from_labels(
    std::span<const std::uint64_t> label_rows,
    const std::vector<std::uint32_t>& new_basis_in_old_coordinates);

// Exhaustively traverse all 2^25 binary symmetric cubic tensors under two
// generators of GL(5,2).  When requested, the separate prefix-pruned direct
// canonicalizer checks the canonical key of every resulting orbit.
[[nodiscard]] FiveQubitTensorOrbitCensus enumerate_five_qubit_tensor_orbits(
    bool check_with_direct_canonicalizer = true);

}  // namespace utsp
