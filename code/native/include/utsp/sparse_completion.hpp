#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

namespace utsp {

struct CanonicalPointSet {
  std::uint32_t ambient_dimension = 0;
  std::uint32_t linear_rank = 0;
  std::vector<std::uint32_t> points;
  std::vector<std::uint32_t> source_basis;
};

struct SparseCompletion {
  std::uint32_t logical_qubits = 0;
  std::uint64_t tensor_signature = 0;
  std::vector<std::uint32_t> columns;
  CanonicalPointSet canonical_form;

  [[nodiscard]] std::uint32_t weight() const {
    return static_cast<std::uint32_t>(columns.size());
  }
};

class SparseCompletionDecoder {
 public:
  struct SubsetRecord {
    std::uint64_t signature = 0;
    std::uint64_t low = 0;
    std::uint64_t high = 0;
  };

  explicit SparseCompletionDecoder(std::uint32_t logical_qubits);

  [[nodiscard]] std::optional<SparseCompletion> decode(
      std::uint64_t tensor_signature,
      std::uint32_t maximum_weight = 7);

  [[nodiscard]] std::uint32_t logical_qubits() const { return q_; }
  [[nodiscard]] std::uint32_t tensor_bit_count() const { return tensor_bits_; }

 private:
  [[nodiscard]] const std::vector<SubsetRecord>& subset_table(
      std::uint32_t weight);

  std::uint32_t q_ = 0;
  std::uint32_t tensor_bits_ = 0;
  std::vector<std::uint64_t> column_signatures_;
  std::array<std::vector<SubsetRecord>, 5> subset_tables_;
  std::array<bool, 5> subset_tables_ready_{};
};

[[nodiscard]] std::uint64_t moment_signature_for_column(
    std::uint32_t column,
    std::uint32_t logical_qubits);

[[nodiscard]] CanonicalPointSet canonicalize_point_set(
    const std::vector<std::uint32_t>& points,
    std::uint32_t ambient_dimension);

}  // namespace utsp
