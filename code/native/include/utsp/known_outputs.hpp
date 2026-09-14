#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace utsp {

struct KnownOutputClass {
  std::string_view class_id;
  std::uint8_t intrinsic_logical_qubits;
  std::uint8_t t_count;
  std::uint16_t four_qubit_canonical_mask;
};

struct CompletionChoice {
  std::uint32_t parent_length = 0;
  std::vector<std::uint32_t> columns;
};

class KnownOutputClassifier {
 public:
  KnownOutputClassifier();

  [[nodiscard]] const KnownOutputClass* classify(
      std::uint32_t logical_qubits,
      std::uint64_t tensor_signature) const;

  [[nodiscard]] std::optional<CompletionChoice> minimum_feasible_completion(
      std::uint32_t logical_qubits,
      std::uint64_t tensor_signature,
      std::uint32_t protocol_length,
      bool support_contains_zero,
      std::uint32_t maximum_parent_length) const;

 private:
  std::array<std::vector<std::vector<std::uint16_t>>, 5> completions_;
  std::array<std::vector<const KnownOutputClass*>, 5> classes_;
};

[[nodiscard]] std::span<const KnownOutputClass> known_output_classes();

}  // namespace utsp
