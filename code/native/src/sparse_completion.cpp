#include "utsp/sparse_completion.hpp"

#include <algorithm>
#include <bit>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace utsp {
namespace {

[[nodiscard]] std::uint32_t tensor_term_count(std::uint32_t q) {
  return q + q * (q - 1) / 2 + q * (q - 1) * (q - 2) / 6;
}

[[nodiscard]] std::uint32_t rank_of_points(
    const std::vector<std::uint32_t>& points) {
  std::array<std::uint32_t, 32> basis{};
  std::uint32_t rank = 0;
  for (auto value : points) {
    while (value != 0) {
      const auto pivot = 31U - static_cast<std::uint32_t>(std::countl_zero(value));
      if (basis[pivot] != 0) {
        value ^= basis[pivot];
      } else {
        basis[pivot] = value;
        ++rank;
        break;
      }
    }
  }
  return rank;
}

[[nodiscard]] std::optional<std::uint32_t> coordinates_in_basis(
    std::uint32_t point,
    const std::vector<std::uint32_t>& basis) {
  const std::uint32_t combinations = std::uint32_t{1} << basis.size();
  std::uint32_t value = 0;
  std::uint32_t previous_gray = 0;
  if (point == 0) return 0;
  for (std::uint32_t index = 1; index < combinations; ++index) {
    const auto gray = index ^ (index >> 1);
    const auto changed = gray ^ previous_gray;
    value ^= basis[std::countr_zero(changed)];
    if (value == point) return gray;
    previous_gray = gray;
  }
  return std::nullopt;
}

[[nodiscard]] bool subset_disjoint(
    const SparseCompletionDecoder::SubsetRecord& left,
    const SparseCompletionDecoder::SubsetRecord& right) {
  return (left.low & right.low) == 0 && (left.high & right.high) == 0;
}

[[nodiscard]] std::vector<std::uint32_t> subset_columns(
    const SparseCompletionDecoder::SubsetRecord& left,
    const SparseCompletionDecoder::SubsetRecord& right,
    std::uint32_t point_count) {
  std::vector<std::uint32_t> result;
  for (std::uint32_t position = 0; position < point_count; ++position) {
    const bool selected =
        position < 64
            ? (((left.low | right.low) >> position) & 1U) != 0
            : (((left.high | right.high) >> (position - 64)) & 1U) != 0;
    if (selected) result.push_back(position + 1);
  }
  return result;
}

}  // namespace

std::uint64_t moment_signature_for_column(
    std::uint32_t column,
    std::uint32_t logical_qubits) {
  if (logical_qubits < 1 || logical_qubits > 7 || column == 0 ||
      column >= (std::uint32_t{1} << logical_qubits)) {
    throw std::invalid_argument("column is outside a q=1,...,7 logical space");
  }
  std::uint64_t signature = 0;
  std::uint32_t bit = 0;
  for (std::uint32_t first = 0; first < logical_qubits; ++first, ++bit) {
    if ((column >> first) & 1U) signature |= std::uint64_t{1} << bit;
  }
  for (std::uint32_t first = 0; first < logical_qubits; ++first) {
    for (std::uint32_t second = first + 1; second < logical_qubits;
         ++second, ++bit) {
      if (((column >> first) & 1U) && ((column >> second) & 1U)) {
        signature |= std::uint64_t{1} << bit;
      }
    }
  }
  for (std::uint32_t first = 0; first < logical_qubits; ++first) {
    for (std::uint32_t second = first + 1; second < logical_qubits; ++second) {
      for (std::uint32_t third = second + 1; third < logical_qubits;
           ++third, ++bit) {
        if (((column >> first) & 1U) && ((column >> second) & 1U) &&
            ((column >> third) & 1U)) {
          signature |= std::uint64_t{1} << bit;
        }
      }
    }
  }
  return signature;
}

CanonicalPointSet canonicalize_point_set(
    const std::vector<std::uint32_t>& points,
    std::uint32_t ambient_dimension) {
  if (ambient_dimension > 31) {
    throw std::invalid_argument("point-set canonicalizer supports q <= 31");
  }
  std::unordered_set<std::uint32_t> distinct;
  for (const auto point : points) {
    if (point == 0 || point >= (std::uint32_t{1} << ambient_dimension) ||
        !distinct.insert(point).second) {
      throw std::invalid_argument("canonical point set is invalid");
    }
  }
  const auto rank = rank_of_points(points);
  if (rank == 0) return CanonicalPointSet{ambient_dimension, 0, {}, {}};

  CanonicalPointSet best;
  best.ambient_dimension = ambient_dimension;
  best.linear_rank = rank;
  std::vector<std::uint32_t> basis;
  std::vector<bool> used(points.size(), false);
  std::function<void()> extend = [&]() {
    if (basis.size() == rank) {
      std::vector<std::uint32_t> transformed;
      transformed.reserve(points.size());
      for (const auto point : points) {
        const auto coordinates = coordinates_in_basis(point, basis);
        if (!coordinates) {
          throw std::logic_error("canonical basis does not span point set");
        }
        transformed.push_back(*coordinates);
      }
      std::sort(transformed.begin(), transformed.end());
      if (best.points.empty() || transformed < best.points) {
        best.points = std::move(transformed);
        best.source_basis = basis;
      }
      return;
    }
    for (std::size_t index = 0; index < points.size(); ++index) {
      if (used[index]) continue;
      basis.push_back(points[index]);
      if (rank_of_points(basis) == basis.size()) {
        used[index] = true;
        extend();
        used[index] = false;
      }
      basis.pop_back();
    }
  };
  extend();
  if (best.points.empty()) throw std::logic_error("point-set canonicalization failed");
  return best;
}

SparseCompletionDecoder::SparseCompletionDecoder(std::uint32_t logical_qubits)
    : q_(logical_qubits), tensor_bits_(tensor_term_count(logical_qubits)) {
  if (q_ < 1 || q_ > 7 || tensor_bits_ > 64) {
    throw std::invalid_argument("sparse decoder supports q=1,...,7");
  }
  const std::uint32_t point_count = (std::uint32_t{1} << q_) - 1;
  column_signatures_.reserve(point_count);
  for (std::uint32_t column = 1; column <= point_count; ++column) {
    column_signatures_.push_back(moment_signature_for_column(column, q_));
  }
}

const std::vector<SparseCompletionDecoder::SubsetRecord>&
SparseCompletionDecoder::subset_table(std::uint32_t weight) {
  if (weight > 4) throw std::invalid_argument("cached subset table stops at weight four");
  if (subset_tables_ready_[weight]) return subset_tables_[weight];

  auto& table = subset_tables_[weight];
  if (weight == 0) {
    table.push_back({0, 0, 0});
  } else {
    std::function<void(std::uint32_t, std::uint32_t, SubsetRecord)> extend;
    extend = [&](std::uint32_t next,
                 std::uint32_t remaining,
                 SubsetRecord record) {
      if (remaining == 0) {
        table.push_back(record);
        return;
      }
      const auto point_count = static_cast<std::uint32_t>(column_signatures_.size());
      for (std::uint32_t position = next;
           position + remaining <= point_count; ++position) {
        auto child = record;
        child.signature ^= column_signatures_[position];
        if (position < 64) {
          child.low |= std::uint64_t{1} << position;
        } else {
          child.high |= std::uint64_t{1} << (position - 64);
        }
        extend(position + 1, remaining - 1, child);
      }
    };
    extend(0, weight, {});
  }
  std::sort(
      table.begin(), table.end(),
      [](const SubsetRecord& left, const SubsetRecord& right) {
        return left.signature < right.signature;
      });
  for (std::size_t index = 1; index < table.size(); ++index) {
    if (table[index - 1].signature == table[index].signature) {
      throw std::logic_error(
          "low-weight moment signatures collided below kernel distance 15");
    }
  }
  subset_tables_ready_[weight] = true;
  return table;
}

std::optional<SparseCompletion> SparseCompletionDecoder::decode(
    std::uint64_t tensor_signature,
    std::uint32_t maximum_weight) {
  if (tensor_bits_ < 64 && tensor_signature >= (std::uint64_t{1} << tensor_bits_)) {
    throw std::invalid_argument("tensor signature lies outside decoder space");
  }
  maximum_weight = std::min<std::uint32_t>(maximum_weight, 7);
  for (std::uint32_t weight = 0; weight <= maximum_weight; ++weight) {
    const auto left_weight = weight / 2;
    const auto right_weight = weight - left_weight;
    const auto& left_table = subset_table(left_weight);
    const auto& right_table = subset_table(right_weight);
    for (const auto& left : left_table) {
      const auto sought = tensor_signature ^ left.signature;
      const auto found = std::lower_bound(
          right_table.begin(), right_table.end(), sought,
          [](const SubsetRecord& record, std::uint64_t signature) {
            return record.signature < signature;
          });
      if (found == right_table.end() || found->signature != sought ||
          !subset_disjoint(left, *found)) {
        continue;
      }
      auto columns = subset_columns(
          left, *found, static_cast<std::uint32_t>(column_signatures_.size()));
      if (columns.size() != weight) continue;
      return SparseCompletion{
          q_,
          tensor_signature,
          columns,
          canonicalize_point_set(columns, q_),
      };
    }
  }
  return std::nullopt;
}

}  // namespace utsp
