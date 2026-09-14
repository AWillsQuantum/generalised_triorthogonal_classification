#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace utsp {

struct TensorUnmarkClass {
  std::uint32_t class_index = 0;
  std::uint64_t representative_standard_signature = 0;
  std::uint64_t orbit_size = 0;
  std::uint64_t stabilizer_size = 0;
  std::uint32_t marked_node_count = 0;
};

struct TensorUnmarkChecks {
  bool authority_and_marked_ledger_headers_match = false;
  bool every_owner_and_transporter_map_header_matches = false;
  bool every_q6_marked_transport_replays = false;
  bool every_transition_affine_map_replays = false;
  bool sampled_transition_edges_replay_as_full_tensors = false;
  bool every_marked_parent_has_complete_extension_mass = false;
  bool component_count_matches_independent_burnside_count = false;
  bool unmarked_orbit_mass_is_two_to_63 = false;
  bool every_orbit_size_divides_gl7_order = false;
  bool class_map_written = false;
};

struct TensorUnmarkResult {
  bool pass = false;
  std::uint32_t q5_parent_count = 0;
  std::uint32_t q6_parent_count = 0;
  std::uint32_t q6_marked_orbit_count = 0;
  std::uint32_t transitions_per_marked_orbit = 0;
  std::uint64_t q7_marked_orbit_count = 0;
  std::uint64_t transition_edge_count = 0;
  std::uint32_t expected_unmarked_orbit_count = 0;
  std::uint64_t transition_affine_replay_checks = 0;
  std::uint64_t transition_edge_replay_checks = 0;
  std::uint64_t q7_tensor_mass = 0;
  TensorUnmarkChecks checks;
  std::vector<TensorUnmarkClass> classes;
  std::vector<std::uint32_t> class_by_marked_node;
};

[[nodiscard]] TensorUnmarkResult unmark_q7_tensor_orbits(
    const std::string& authority_path,
    const std::string& marked_ledger_path,
    const std::string& q6_map_directory,
    const std::string& q7_owner_map_directory,
    const std::string& class_map_output_path,
    std::uint32_t workers = 1);

}  // namespace utsp
