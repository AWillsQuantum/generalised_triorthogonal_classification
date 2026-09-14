#pragma once

#include <cstdint>
#include <span>
#include <vector>

namespace utsp {

struct TensorExtensionAffineAction {
  std::uint32_t offset = 0;
  std::vector<std::uint32_t> columns;
};

struct TensorExtensionOrbit {
  std::uint32_t representative = 0;
  std::uint32_t orbit_size = 0;
};

struct TensorExtensionCensus {
  std::uint32_t parent_logical_qubits = 0;
  std::uint64_t parent_standard_signature = 0;
  std::uint32_t extension_bit_count = 0;
  std::uint32_t extension_state_count = 0;
  std::uint64_t parent_automorphism_group_order = 0;
  std::uint64_t marked_stabilizer_group_order = 0;
  std::uint64_t generator_transitions = 0;
  std::uint64_t transporter_replay_checks = 0;
  std::uint64_t transporter_replay_mismatches = 0;
  std::vector<TensorExtensionAffineAction> actions;
  std::vector<TensorExtensionOrbit> orbits;
  std::vector<std::uint32_t> orbit_index_by_state;
  std::vector<std::uint64_t> transport_to_representative_basis_by_state;
};

[[nodiscard]] std::uint32_t tensor_extension_bit_count(
    std::uint32_t parent_logical_qubits);

[[nodiscard]] std::uint64_t embed_tensor_parent_signature(
    std::uint32_t parent_logical_qubits,
    std::uint64_t parent_standard_signature);

[[nodiscard]] std::uint64_t combine_tensor_parent_and_extension(
    std::uint32_t parent_logical_qubits,
    std::uint64_t parent_standard_signature,
    std::uint32_t extension_state);

[[nodiscard]] std::uint64_t restrict_tensor_to_standard_hyperplane(
    std::uint32_t child_logical_qubits,
    std::uint64_t child_standard_signature);

[[nodiscard]] std::uint32_t tensor_extension_state(
    std::uint32_t parent_logical_qubits,
    std::uint64_t child_standard_signature);

[[nodiscard]] TensorExtensionCensus census_tensor_extensions(
    std::uint32_t parent_logical_qubits,
    std::uint64_t parent_standard_signature,
    std::span<const std::vector<std::uint32_t>>
        parent_automorphism_basis_generators,
    std::uint64_t parent_automorphism_group_order,
    bool retain_orbit_index_by_state = false,
    bool retain_transport_to_representative_basis = false);

}  // namespace utsp
