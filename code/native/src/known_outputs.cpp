#include "utsp/known_outputs.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <deque>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace utsp {
namespace {

constexpr std::array<KnownOutputClass, 22> kKnownClasses{{
    {"D1_01", 1, 1, 1},
    {"D2_01", 2, 2, 3},
    {"D2_02", 2, 3, 7},
    {"D3_01", 3, 3, 11},
    {"D3_02", 3, 4, 15},
    {"D3_03", 3, 4, 30},
    {"D3_04", 3, 5, 31},
    {"D3_05", 3, 6, 63},
    {"D3_06", 3, 7, 127},
    {"D4_01", 4, 4, 139},
    {"D4_02", 4, 5, 143},
    {"D4_03", 4, 5, 158},
    {"D4_04", 4, 6, 159},
    {"D4_05", 4, 7, 191},
    {"D4_06", 4, 7, 255},
    {"D4_07", 4, 6, 414},
    {"D4_08", 4, 7, 415},
    {"D4_09", 4, 6, 427},
    {"D4_10", 4, 5, 428},
    {"D4_11", 4, 6, 429},
    {"D4_12", 4, 7, 431},
    {"D4_13", 4, 7, 446},
}};

[[nodiscard]] std::vector<std::uint32_t> monomial_masks(std::uint32_t q) {
  std::vector<std::uint32_t> result;
  for (std::uint32_t first = 0; first < q; ++first) {
    result.push_back(std::uint32_t{1} << first);
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      result.push_back((std::uint32_t{1} << first) |
                       (std::uint32_t{1} << second));
    }
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      for (std::uint32_t third = second + 1; third < q; ++third) {
        result.push_back((std::uint32_t{1} << first) |
                         (std::uint32_t{1} << second) |
                         (std::uint32_t{1} << third));
      }
    }
  }
  return result;
}

[[nodiscard]] std::uint32_t column_signature(
    std::uint32_t column,
    std::span<const std::uint32_t> monomials) {
  std::uint32_t signature = 0;
  for (std::size_t index = 0; index < monomials.size(); ++index) {
    if ((column & monomials[index]) == monomials[index]) {
      signature |= std::uint32_t{1} << index;
    }
  }
  return signature;
}

[[nodiscard]] std::vector<std::vector<std::uint16_t>> completion_table(
    std::uint32_t q) {
  const auto monomials = monomial_masks(q);
  const std::uint32_t point_count = (std::uint32_t{1} << q) - 1;
  std::vector<std::uint32_t> signatures_by_column;
  signatures_by_column.reserve(point_count);
  for (std::uint32_t column = 1; column <= point_count; ++column) {
    signatures_by_column.push_back(column_signature(column, monomials));
  }
  std::vector<std::vector<std::uint16_t>> result(
      std::size_t{1} << monomials.size());
  const std::uint32_t subset_count = std::uint32_t{1} << point_count;
  std::vector<std::uint16_t> signatures(subset_count, 0);
  for (std::uint32_t mask = 0; mask < subset_count; ++mask) {
    if (mask != 0) {
      const auto bit = static_cast<std::uint32_t>(std::countr_zero(mask));
      signatures[mask] = signatures[mask & (mask - 1)] ^
                         signatures_by_column[bit];
    }
    result[signatures[mask]].push_back(static_cast<std::uint16_t>(mask));
  }
  return result;
}

using PointPermutation = std::array<std::uint8_t, 16>;

[[nodiscard]] std::vector<PointPermutation> linear_generators() {
  std::vector<PointPermutation> generators;
  for (std::uint32_t left = 0; left < 3; ++left) {
    const std::uint32_t right = left + 1;
    PointPermutation permutation{};
    for (std::uint32_t vector = 1; vector < 16; ++vector) {
      std::uint32_t transformed = vector;
      const auto left_bit = (vector >> left) & 1U;
      const auto right_bit = (vector >> right) & 1U;
      transformed &= ~((std::uint32_t{1} << left) |
                       (std::uint32_t{1} << right));
      transformed |= left_bit << right;
      transformed |= right_bit << left;
      permutation[vector] = static_cast<std::uint8_t>(transformed);
    }
    generators.push_back(permutation);
  }
  for (std::uint32_t target = 0; target < 4; ++target) {
    for (std::uint32_t control = 0; control < 4; ++control) {
      if (target == control) continue;
      PointPermutation permutation{};
      for (std::uint32_t vector = 1; vector < 16; ++vector) {
        std::uint32_t transformed = vector;
        if ((vector >> control) & 1U) transformed ^= std::uint32_t{1} << target;
        permutation[vector] = static_cast<std::uint8_t>(transformed);
      }
      generators.push_back(permutation);
    }
  }
  return generators;
}

[[nodiscard]] std::uint16_t transform_completion(
    std::uint16_t mask,
    const PointPermutation& permutation) {
  std::uint16_t transformed = 0;
  while (mask != 0) {
    const auto bit = static_cast<std::uint32_t>(std::countr_zero(mask));
    transformed |= std::uint16_t{1} << (permutation[bit + 1] - 1);
    mask &= mask - 1;
  }
  return transformed;
}

[[nodiscard]] std::vector<std::uint16_t> four_qubit_canonical_map() {
  constexpr std::uint16_t full = (std::uint16_t{1} << 15) - 1;
  constexpr std::uint16_t unseen = std::numeric_limits<std::uint16_t>::max();
  const auto generators = linear_generators();
  std::vector<std::uint16_t> canonical(std::size_t{1} << 14, unseen);
  for (std::uint16_t start = 0; start < (std::uint16_t{1} << 14); ++start) {
    if (canonical[start] != unseen) continue;
    std::vector<bool> in_orbit(std::size_t{1} << 14, false);
    std::deque<std::uint16_t> queue;
    std::vector<std::uint16_t> orbit;
    queue.push_back(start);
    in_orbit[start] = true;
    while (!queue.empty()) {
      const auto current = queue.front();
      queue.pop_front();
      orbit.push_back(current);
      for (const auto& permutation : generators) {
        auto transformed = transform_completion(current, permutation);
        transformed = std::min<std::uint16_t>(transformed, full ^ transformed);
        if (in_orbit[transformed]) continue;
        in_orbit[transformed] = true;
        queue.push_back(transformed);
      }
    }
    const auto representative = *std::min_element(orbit.begin(), orbit.end());
    for (const auto member : orbit) canonical[member] = representative;
  }
  return canonical;
}

[[nodiscard]] std::uint16_t pad_signature_to_four(
    std::uint32_t q, std::uint64_t signature) {
  const auto source = monomial_masks(q);
  const auto target = monomial_masks(4);
  std::uint16_t padded = 0;
  for (std::size_t index = 0; index < source.size(); ++index) {
    if (((signature >> index) & 1U) == 0) continue;
    const auto position = std::find(target.begin(), target.end(), source[index]);
    if (position == target.end()) throw std::logic_error("signature padding failed");
    padded |= std::uint16_t{1} << std::distance(target.begin(), position);
  }
  return padded;
}

[[nodiscard]] std::vector<std::uint32_t> mask_columns(
    std::uint16_t mask, std::uint32_t q) {
  std::vector<std::uint32_t> columns;
  for (std::uint32_t column = 1; column < (std::uint32_t{1} << q); ++column) {
    if ((mask >> (column - 1)) & 1U) columns.push_back(column);
  }
  return columns;
}

}  // namespace

std::span<const KnownOutputClass> known_output_classes() {
  return kKnownClasses;
}

KnownOutputClassifier::KnownOutputClassifier() {
  for (std::uint32_t q = 1; q <= 4; ++q) {
    completions_[q] = completion_table(q);
  }
  const auto canonical = four_qubit_canonical_map();
  std::unordered_map<std::uint16_t, const KnownOutputClass*> by_canonical;
  for (const auto& output : kKnownClasses) {
    by_canonical.emplace(output.four_qubit_canonical_mask, &output);
  }
  constexpr std::uint16_t full = (std::uint16_t{1} << 15) - 1;
  for (std::uint32_t q = 1; q <= 4; ++q) {
    classes_[q].assign(completions_[q].size(), nullptr);
    for (std::size_t signature = 1; signature < classes_[q].size(); ++signature) {
      const auto padded = pad_signature_to_four(q, signature);
      const auto& masks = completions_[4][padded];
      if (masks.empty()) throw std::logic_error("four-qubit moment map is not onto");
      const auto completion = *std::min_element(masks.begin(), masks.end());
      const auto quotient = std::min<std::uint16_t>(completion, full ^ completion);
      const auto found = by_canonical.find(canonical[quotient]);
      if (found == by_canonical.end()) {
        throw std::logic_error("nonzero output orbit has no known class");
      }
      classes_[q][signature] = found->second;
    }
  }
}

const KnownOutputClass* KnownOutputClassifier::classify(
    std::uint32_t logical_qubits,
    std::uint64_t tensor_signature) const {
  if (logical_qubits < 1 || logical_qubits > 4) {
    throw std::invalid_argument("known-output classifier supports q=1,...,4");
  }
  if (tensor_signature >= classes_[logical_qubits].size()) {
    throw std::invalid_argument("tensor signature is outside its q space");
  }
  return classes_[logical_qubits][tensor_signature];
}

std::optional<CompletionChoice> KnownOutputClassifier::minimum_feasible_completion(
    std::uint32_t logical_qubits,
    std::uint64_t tensor_signature,
    std::uint32_t protocol_length,
    bool support_contains_zero,
    std::uint32_t maximum_parent_length) const {
  if (logical_qubits < 1 || logical_qubits > 4 ||
      tensor_signature >= completions_[logical_qubits].size()) {
    throw std::invalid_argument("completion query is outside q=1,...,4");
  }
  std::optional<CompletionChoice> best;
  for (const auto mask : completions_[logical_qubits][tensor_signature]) {
    auto columns = mask_columns(mask, logical_qubits);
    const auto lower_length =
        protocol_length + static_cast<std::uint32_t>(columns.size());
    const bool explicit_zero = (lower_length & 1U) != 0;
    if (explicit_zero && support_contains_zero) continue;
    const auto parent_length = lower_length + static_cast<std::uint32_t>(explicit_zero);
    if (parent_length > maximum_parent_length) continue;
    CompletionChoice candidate{parent_length, std::move(columns)};
    if (!best || candidate.parent_length < best->parent_length ||
        (candidate.parent_length == best->parent_length &&
         candidate.columns < best->columns)) {
      best = std::move(candidate);
    }
  }
  return best;
}

}  // namespace utsp
