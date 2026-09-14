#include "utsp/tensor_decomposition.hpp"

#include "utsp/native_core.hpp"
#include "utsp/tensor_canonical.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>

namespace utsp {
namespace {

[[nodiscard]] bool wide_zero(const CentroidElement& value) {
  return std::all_of(
      value.begin(), value.end(), [](std::uint64_t word) { return word == 0; });
}

[[nodiscard]] bool wide_equal(
    const CentroidElement& left,
    const CentroidElement& right) {
  return left == right;
}

void wide_set(CentroidElement& value, std::uint32_t bit) {
  value[bit / 64] |= std::uint64_t{1} << (bit % 64);
}

void wide_toggle(CentroidElement& value, std::uint32_t bit) {
  value[bit / 64] ^= std::uint64_t{1} << (bit % 64);
}

[[nodiscard]] bool wide_test(const CentroidElement& value, std::uint32_t bit) {
  return ((value[bit / 64] >> (bit % 64)) & 1U) != 0;
}

void wide_xor(CentroidElement& target, const CentroidElement& source) {
  for (std::size_t word = 0; word < target.size(); ++word) {
    target[word] ^= source[word];
  }
}

[[nodiscard]] bool wide_dot(
    const CentroidElement& left,
    const CentroidElement& right) {
  unsigned parity = 0;
  for (std::size_t word = 0; word < left.size(); ++word) {
    parity ^= std::popcount(left[word] & right[word]) & 1U;
  }
  return parity != 0;
}

[[nodiscard]] std::uint32_t wide_pivot(const CentroidElement& value) {
  for (std::size_t word = value.size(); word-- > 0;) {
    if (value[word] != 0) {
      return static_cast<std::uint32_t>(
          64 * word + 63 - std::countl_zero(value[word]));
    }
  }
  throw std::invalid_argument("zero wide vector has no pivot");
}

[[nodiscard]] std::vector<CentroidElement> wide_rref(
    std::span<const CentroidElement> vectors) {
  std::map<std::uint32_t, CentroidElement, std::greater<>> basis;
  for (auto value : vectors) {
    for (const auto& [pivot, row] : basis) {
      if (wide_test(value, pivot)) wide_xor(value, row);
    }
    if (wide_zero(value)) continue;
    const auto pivot = wide_pivot(value);
    for (auto& [other_pivot, row] : basis) {
      (void)other_pivot;
      if (wide_test(row, pivot)) wide_xor(row, value);
    }
    basis.emplace(pivot, value);
  }
  std::vector<CentroidElement> result;
  result.reserve(basis.size());
  for (const auto& [pivot, row] : basis) {
    (void)pivot;
    result.push_back(row);
  }
  return result;
}

[[nodiscard]] std::vector<CentroidElement> wide_nullspace(
    std::span<const CentroidElement> rows,
    std::uint32_t width) {
  if (width > 256) throw std::invalid_argument("wide nullspace exceeds 256 bits");
  const auto reduced = wide_rref(rows);
  std::vector<bool> pivots(width, false);
  for (const auto& row : reduced) pivots[wide_pivot(row)] = true;
  std::vector<CentroidElement> basis;
  for (std::uint32_t free = 0; free < width; ++free) {
    if (pivots[free]) continue;
    CentroidElement vector{};
    wide_set(vector, free);
    for (const auto& row : reduced) {
      if (wide_dot(row, vector)) wide_set(vector, wide_pivot(row));
    }
    basis.push_back(vector);
  }
  return basis;
}

[[nodiscard]] CentroidElement wide_reduce(
    CentroidElement vector,
    std::span<const CentroidElement> basis) {
  const auto reduced = wide_rref(basis);
  for (const auto& row : reduced) {
    const auto pivot = wide_pivot(row);
    if (wide_test(vector, pivot)) wide_xor(vector, row);
  }
  return vector;
}

[[nodiscard]] CentroidElement wide_linear_combination(
    std::uint64_t coefficients,
    std::span<const CentroidElement> basis) {
  CentroidElement result{};
  while (coefficients != 0) {
    const auto bit = static_cast<std::size_t>(std::countr_zero(coefficients));
    wide_xor(result, basis[bit]);
    coefficients &= coefficients - 1;
  }
  return result;
}

class LabelTensorValues {
 public:
  explicit LabelTensorValues(std::span<const std::uint64_t> label_rows)
      : combined_rows_(std::size_t{1} << label_rows.size(), 0) {
    for (std::size_t vector = 1; vector < combined_rows_.size(); ++vector) {
      const auto bit = static_cast<std::size_t>(std::countr_zero(vector));
      combined_rows_[vector] =
          combined_rows_[vector & (vector - 1)] ^ label_rows[bit];
    }
  }

  [[nodiscard]] bool operator()(
      std::uint32_t left,
      std::uint32_t middle,
      std::uint32_t right) const {
    return (std::popcount(
                combined_rows_[left] & combined_rows_[middle] &
                combined_rows_[right]) &
            1U) != 0;
  }

 private:
  std::vector<std::uint64_t> combined_rows_;
};

[[nodiscard]] std::vector<std::uint64_t> matrix_columns(
    const CentroidElement& matrix,
    std::uint32_t q) {
  std::vector<std::uint64_t> columns(q, 0);
  for (std::uint32_t input = 0; input < q; ++input) {
    for (std::uint32_t output = 0; output < q; ++output) {
      if (wide_test(matrix, output * q + input)) {
        columns[input] |= std::uint64_t{1} << output;
      }
    }
  }
  return columns;
}

[[nodiscard]] std::uint64_t apply_matrix(
    std::span<const std::uint64_t> columns,
    std::uint64_t vector) {
  std::uint64_t result = 0;
  while (vector != 0) {
    const auto bit = static_cast<std::size_t>(std::countr_zero(vector));
    result ^= columns[bit];
    vector &= vector - 1;
  }
  return result;
}

[[nodiscard]] bool is_idempotent(
    const CentroidElement& matrix,
    std::uint32_t q) {
  const auto columns = matrix_columns(matrix, q);
  for (std::uint32_t input = 0; input < q; ++input) {
    if (apply_matrix(columns, columns[input]) != columns[input]) return false;
  }
  return true;
}

[[nodiscard]] TensorIdempotentSplit split_from_idempotent(
    const CentroidElement& matrix,
    std::uint32_t q) {
  const auto columns = matrix_columns(matrix, q);
  auto image = rref_basis(columns);
  std::vector<std::uint64_t> rows(q, 0);
  for (std::uint32_t output = 0; output < q; ++output) {
    for (std::uint32_t input = 0; input < q; ++input) {
      if ((columns[input] >> output) & 1U) {
        rows[output] |= std::uint64_t{1} << input;
      }
    }
  }
  auto kernel = nullspace_basis(rows, q);
  if (image.empty() || kernel.empty() || image.size() + kernel.size() != q) {
    throw std::logic_error("nontrivial idempotent gave an invalid direct split");
  }
  std::vector<std::uint64_t> combined = image;
  combined.insert(combined.end(), kernel.begin(), kernel.end());
  if (gf2_rank(combined) != q) {
    throw std::logic_error("idempotent image and kernel do not complement");
  }
  return TensorIdempotentSplit{matrix, std::move(image), std::move(kernel)};
}

[[nodiscard]] std::vector<std::uint64_t> combine_label_rows(
    std::span<const std::uint64_t> coefficients,
    std::span<const std::uint64_t> label_rows) {
  std::vector<std::uint64_t> result;
  result.reserve(coefficients.size());
  for (const auto vector : coefficients) {
    result.push_back(linear_combination(vector, label_rows));
  }
  return result;
}

[[nodiscard]] bool component_less(
    const TensorCanonicalComponent& left,
    const TensorCanonicalComponent& right) {
  return std::tie(left.dimension, left.canonical_key_words) <
         std::tie(right.dimension, right.canonical_key_words);
}

struct RecursiveCanonicalization {
  std::vector<std::uint64_t> radical_basis;
  std::vector<TensorCanonicalComponent> components;
  std::uint64_t centroid_elements_tested = 0;
};

[[nodiscard]] RecursiveCanonicalization canonicalize_recursive(
    std::span<const std::uint64_t> label_rows) {
  const auto q = static_cast<std::uint32_t>(label_rows.size());
  RecursiveCanonicalization result;
  const auto radical = tensor_radical_basis(label_rows);
  std::vector<std::uint64_t> span_basis = radical;
  std::vector<std::uint64_t> complement;
  for (std::uint32_t coordinate = 0; coordinate < q; ++coordinate) {
    const std::uint64_t unit = std::uint64_t{1} << coordinate;
    if (reduce_vector(unit, span_basis) == 0) continue;
    complement.push_back(unit);
    span_basis.push_back(unit);
    span_basis = rref_basis(span_basis);
  }
  if (span_basis.size() != q) {
    throw std::logic_error("failed to complement the tensor radical");
  }
  result.radical_basis = radical;
  if (complement.empty()) return result;

  const auto core_rows = combine_label_rows(complement, label_rows);
  if (!tensor_radical_basis(core_rows).empty()) {
    throw std::logic_error("tensor quotient by its radical remained degenerate");
  }
  const auto centroid = analyze_tensor_centroid(core_rows);
  result.centroid_elements_tested += centroid.centroid_elements_tested;
  if (centroid.enumeration_limit_exceeded) {
    throw std::runtime_error("centroid idempotent enumeration limit was exceeded");
  }
  if (!centroid.split) {
    if (core_rows.size() > 8) {
      throw std::runtime_error(
          "primitive tensor component exceeds the direct q=8 gate");
    }
    const auto direct = canonicalize_cubic_tensor_from_labels(core_rows);
    TensorCanonicalComponent component;
    component.dimension = static_cast<std::uint32_t>(core_rows.size());
    component.canonical_key_words = direct.canonical_key_words;
    for (const auto vector : direct.new_basis_in_old_coordinates) {
      component.basis_in_original_coordinates.push_back(
          linear_combination(vector, complement));
    }
    result.components.push_back(std::move(component));
    return result;
  }

  const std::array<std::vector<std::uint64_t>, 2> split_bases{
      centroid.split->image_basis,
      centroid.split->kernel_basis,
  };
  for (const auto& split_basis : split_bases) {
    const auto child_rows = combine_label_rows(split_basis, core_rows);
    auto child = canonicalize_recursive(child_rows);
    result.centroid_elements_tested += child.centroid_elements_tested;
    for (const auto child_radical : child.radical_basis) {
      const auto in_core = linear_combination(child_radical, split_basis);
      result.radical_basis.push_back(linear_combination(in_core, complement));
    }
    for (auto& component : child.components) {
      for (auto& vector : component.basis_in_original_coordinates) {
        const auto in_core = linear_combination(vector, split_basis);
        vector = linear_combination(in_core, complement);
      }
      result.components.push_back(std::move(component));
    }
  }
  std::sort(result.components.begin(), result.components.end(), component_less);
  return result;
}

}  // namespace

TensorCentroidAnalysis analyze_tensor_centroid(
    std::span<const std::uint64_t> label_rows,
    std::uint32_t maximum_enumerated_dimension) {
  const auto q = static_cast<std::uint32_t>(label_rows.size());
  if (q < 1 || q > 16) {
    throw std::invalid_argument("centroid analyzer currently supports q<=16");
  }
  const LabelTensorValues values(label_rows);
  std::vector<CentroidElement> constraints;
  constraints.reserve(static_cast<std::size_t>(q) * q * q);
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = 0; second < q; ++second) {
      for (std::uint32_t third = 0; third < q; ++third) {
        CentroidElement equation{};
        for (std::uint32_t image = 0; image < q; ++image) {
          if (values(1U << image, 1U << second, 1U << third)) {
            wide_toggle(equation, image * q + first);
          }
          if (values(1U << first, 1U << image, 1U << third)) {
            wide_toggle(equation, image * q + second);
          }
        }
        constraints.push_back(equation);
      }
    }
  }

  TensorCentroidAnalysis result;
  result.logical_qubits = q;
  result.centroid_basis = wide_nullspace(constraints, q * q);
  CentroidElement identity{};
  for (std::uint32_t index = 0; index < q; ++index) {
    wide_set(identity, index * q + index);
  }
  if (!wide_zero(wide_reduce(identity, result.centroid_basis))) {
    throw std::logic_error("tensor centroid does not contain the identity");
  }
  if (result.centroid_basis.size() > maximum_enumerated_dimension ||
      result.centroid_basis.size() >= 63) {
    result.enumeration_limit_exceeded = true;
    return result;
  }

  const std::uint64_t element_count =
      std::uint64_t{1} << result.centroid_basis.size();
  for (std::uint64_t coefficients = 1; coefficients < element_count; ++coefficients) {
    const auto element =
        wide_linear_combination(coefficients, result.centroid_basis);
    ++result.centroid_elements_tested;
    if (wide_equal(element, identity) || !is_idempotent(element, q)) continue;
    result.split = split_from_idempotent(element, q);
    break;
  }
  return result;
}

DecomposedTensorCanonicalForm canonicalize_tensor_by_decomposition(
    std::span<const std::uint64_t> label_rows) {
  if (label_rows.empty() || label_rows.size() > 16) {
    throw std::invalid_argument("decomposition canonicalizer supports q<=16");
  }
  auto recursive = canonicalize_recursive(label_rows);
  std::sort(recursive.components.begin(), recursive.components.end(), component_less);

  DecomposedTensorCanonicalForm result;
  result.logical_qubits = static_cast<std::uint32_t>(label_rows.size());
  result.radical_basis = std::move(recursive.radical_basis);
  result.components = std::move(recursive.components);
  result.centroid_elements_tested = recursive.centroid_elements_tested;
  for (const auto& component : result.components) {
    result.canonical_basis_in_original_coordinates.insert(
        result.canonical_basis_in_original_coordinates.end(),
        component.basis_in_original_coordinates.begin(),
        component.basis_in_original_coordinates.end());
  }
  result.canonical_basis_in_original_coordinates.insert(
      result.canonical_basis_in_original_coordinates.end(),
      result.radical_basis.begin(), result.radical_basis.end());
  if (gf2_rank(result.canonical_basis_in_original_coordinates) !=
      result.logical_qubits) {
    throw std::logic_error("decomposition transporter is not a logical basis");
  }
  return result;
}

}  // namespace utsp
