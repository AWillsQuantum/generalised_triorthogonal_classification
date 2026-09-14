#pragma once

#include "utsp/native_core.hpp"

#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <vector>

namespace utsp {

enum class ProtocolEnumerationMode : std::uint8_t {
  raw_subspaces,
  marked_code_orbits,
  final_dimension_nondegenerate_marked_code_orbits,
  q3_predecessor_target_marked_code_orbits,
  q3_predecessor_target_existence,
  q5_q3_chain_cover_marked_code_orbits,
  q5_q4_hitting_set_marked_code_orbits,
  q7_primitive_marked_code_orbits,
  hybrid,
};

struct ProtocolSupport {
  std::uint64_t record_index = 0;
  std::uint32_t source_index = 0;
  std::string source_id;
  std::string source_family;
  std::uint16_t source_space_length = 0;
  std::uint8_t ambient_dimension = 0;
  std::uint8_t parity_case = 0;
  std::int32_t origin = -1;
  std::vector<std::uint32_t> points;
};

struct OutputTensorKey {
  std::uint32_t logical_qubits = 0;
  std::vector<std::uint64_t> words;

  [[nodiscard]] bool operator<(const OutputTensorKey& other) const;
  [[nodiscard]] bool operator==(const OutputTensorKey& other) const = default;
};

struct ProtocolFrontierWitness {
  OutputTensorKey output;
  std::string known_class_id;
  std::uint32_t distance = 0;
  std::uint64_t error_coefficient = 0;
  std::uint32_t protocol_length = 0;
  std::uint32_t space_footprint = 0;
  std::uint64_t record_index = 0;
  std::uint32_t source_index = 0;
  std::string source_id;
  std::string source_family;
  std::uint16_t source_space_length = 0;
  std::uint8_t ambient_dimension = 0;
  std::uint8_t parity_case = 0;
  std::int32_t origin = -1;
  std::vector<std::uint32_t> points;
  std::vector<Mask> canonical_logical_rows;
};

struct ProtocolPositiveSupport {
  std::uint64_t support_offset = 0;
  std::uint64_t record_index = 0;
  std::vector<OutputTensorKey> output_keys;
};

struct ProtocolCatalogueResult {
  std::uint32_t logical_qubits = 0;
  std::uint32_t minimum_distance = 3;
  ProtocolEnumerationMode enumeration_mode =
      ProtocolEnumerationMode::raw_subspaces;
  std::uint32_t raw_maximum_quotient_dimension = 0;
  bool zero_isotropic_proof = false;
  std::uint32_t zero_isotropic_workers = 0;
  bool isotropic_subspace_statistics_exact = true;
  bool radical_statistics_exact = true;
  std::uint64_t supports_processed = 0;
  std::uint64_t quotient_dimension_filtered_supports = 0;
  std::uint64_t eligible_supports = 0;
  std::uint64_t raw_enumerated_supports = 0;
  std::uint64_t marked_orbit_enumerated_supports = 0;
  std::uint64_t isotropic_subspaces = 0;
  std::uint64_t nondegenerate_subspaces = 0;
  std::uint64_t marked_code_orbits = 0;
  std::uint64_t exact_distance_calls = 0;
  std::uint64_t raw_tensor_forms = 0;
  std::uint64_t canonical_output_orbits = 0;
  std::map<std::uint32_t, std::uint64_t> radical_dimension_counts;
  std::vector<std::vector<std::uint64_t>> raw_tensor_keys;
  std::vector<OutputTensorKey> canonical_output_keys;
  std::vector<ProtocolPositiveSupport> positive_supports;
  std::vector<ProtocolFrontierWitness> pareto_witnesses;
};

// Compact non-owning view of one cached quotient basis per support.  Offsets
// has support_count+1 entries into rows; an empty offsets span means no cache.
struct QuotientBasisCacheView {
  std::span<const std::uint64_t> offsets;
  std::span<const Mask> rows;

  [[nodiscard]] bool provided() const { return !offsets.empty(); }
  [[nodiscard]] std::size_t size() const;
  [[nodiscard]] std::span<const Mask> basis(std::size_t index) const;
  [[nodiscard]] QuotientBasisCacheView subview(
      std::size_t begin, std::size_t count) const;
};

// Classify every nondegenerate output tensor at one logical dimension. Tensors
// with radical are counted but omitted from the magic-output frontier because
// deleting their Clifford-trivial logical rows strictly improves footprint.
[[nodiscard]] ProtocolCatalogueResult classify_protocol_frontier(
    std::span<const ProtocolSupport> supports,
    std::uint32_t logical_qubits,
    std::uint32_t minimum_distance = 3,
    std::uint32_t maximum_distance_dp_dimension = 24,
    ProtocolEnumerationMode enumeration_mode =
        ProtocolEnumerationMode::raw_subspaces,
    std::uint32_t raw_maximum_quotient_dimension = 14,
    std::uint32_t orbit_workers = 1,
    std::uint32_t support_workers = 1,
    std::uint32_t zero_isotropic_workers = 0,
    std::uint32_t minimum_quotient_dimension = 0,
    std::uint32_t maximum_quotient_dimension =
        std::numeric_limits<std::uint32_t>::max(),
    bool collect_positive_supports = false,
    QuotientBasisCacheView cached_quotient_bases = {});

}  // namespace utsp
