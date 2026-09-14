#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace utsp {

using CentroidElement = std::array<std::uint64_t, 4>;

struct TensorIdempotentSplit {
  CentroidElement centroid_element{};
  std::vector<std::uint64_t> image_basis;
  std::vector<std::uint64_t> kernel_basis;
};

struct TensorCentroidAnalysis {
  std::uint32_t logical_qubits = 0;
  std::vector<CentroidElement> centroid_basis;
  std::optional<TensorIdempotentSplit> split;
  std::uint64_t centroid_elements_tested = 0;
  bool enumeration_limit_exceeded = false;
};

struct TensorCanonicalComponent {
  std::uint32_t dimension = 0;
  std::vector<std::uint64_t> canonical_key_words;
  std::vector<std::uint64_t> basis_in_original_coordinates;
};

struct DecomposedTensorCanonicalForm {
  std::uint32_t logical_qubits = 0;
  std::vector<std::uint64_t> radical_basis;
  std::vector<TensorCanonicalComponent> components;
  std::vector<std::uint64_t> canonical_basis_in_original_coordinates;
  std::uint64_t centroid_elements_tested = 0;
};

[[nodiscard]] TensorCentroidAnalysis analyze_tensor_centroid(
    std::span<const std::uint64_t> label_rows,
    std::uint32_t maximum_enumerated_dimension = 20);

[[nodiscard]] DecomposedTensorCanonicalForm canonicalize_tensor_by_decomposition(
    std::span<const std::uint64_t> label_rows);

}  // namespace utsp
