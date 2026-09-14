#include "utsp/tensor_canonical.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>

namespace utsp {

std::uint32_t cubic_tensor_term_count(std::uint32_t q) {
  return q + q * (q - 1) / 2 + q * (q - 1) * (q - 2) / 6;
}

std::uint32_t cubic_tensor_repeated_term_count(std::uint32_t q) {
  return q + q * (q - 1) / 2;
}

std::uint32_t cubic_tensor_alternating_term_count(std::uint32_t q) {
  return q * (q - 1) * (q - 2) / 6;
}

bool cubic_tensor_value(
    std::uint32_t q,
    std::uint64_t signature,
    std::uint32_t left,
    std::uint32_t middle,
    std::uint32_t right) {
  if (q < 1 || q > 7 ||
      left >= (std::uint32_t{1} << q) ||
      middle >= (std::uint32_t{1} << q) ||
      right >= (std::uint32_t{1} << q)) {
    throw std::invalid_argument("tensor evaluation requires q=1,...,7 vectors");
  }
  const auto terms = cubic_tensor_term_count(q);
  if (signature >= (std::uint64_t{1} << terms)) {
    throw std::invalid_argument("tensor signature lies outside its q space");
  }
  bool value = false;
  std::uint32_t bit = 0;
  for (std::uint32_t first = 0; first < q; ++first, ++bit) {
    if (((signature >> bit) & 1U) == 0) continue;
    value ^= ((left >> first) & 1U) && ((middle >> first) & 1U) &&
             ((right >> first) & 1U);
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second, ++bit) {
      if (((signature >> bit) & 1U) == 0) continue;
      const bool lf = (left >> first) & 1U;
      const bool ls = (left >> second) & 1U;
      const bool mf = (middle >> first) & 1U;
      const bool ms = (middle >> second) & 1U;
      const bool rf = (right >> first) & 1U;
      const bool rs = (right >> second) & 1U;
      value ^= lf && mf && rs;
      value ^= lf && ms && rf;
      value ^= ls && mf && rf;
      value ^= lf && ms && rs;
      value ^= ls && mf && rs;
      value ^= ls && ms && rf;
    }
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      for (std::uint32_t third = second + 1; third < q; ++third, ++bit) {
        if (((signature >> bit) & 1U) == 0) continue;
        for (const auto [left_index, middle_index, right_index] :
             std::array<std::array<std::uint32_t, 3>, 6>{{
                 {first, second, third},
                 {first, third, second},
                 {second, first, third},
                 {second, third, first},
                 {third, first, second},
                 {third, second, first},
             }}) {
          value ^= ((left >> left_index) & 1U) &&
                   ((middle >> middle_index) & 1U) &&
                   ((right >> right_index) & 1U);
        }
      }
    }
  }
  return value;
}

namespace {

class TensorValues {
 public:
  TensorValues(std::uint32_t q, std::uint64_t signature)
      : size_(std::uint32_t{1} << q),
        values_(static_cast<std::size_t>(size_) * size_ * size_, 0) {
    for (std::uint32_t left = 0; left < size_; ++left) {
      for (std::uint32_t middle = 0; middle < size_; ++middle) {
        for (std::uint32_t right = 0; right < size_; ++right) {
          values_[index(left, middle, right)] =
              cubic_tensor_value(q, signature, left, middle, right);
        }
      }
    }
  }

  [[nodiscard]] bool operator()(
      std::uint32_t left,
      std::uint32_t middle,
      std::uint32_t right) const {
    return values_[index(left, middle, right)] != 0;
  }

 private:
  [[nodiscard]] std::size_t index(
      std::uint32_t left,
      std::uint32_t middle,
      std::uint32_t right) const {
    return (static_cast<std::size_t>(left) * size_ + middle) * size_ + right;
  }

  std::uint32_t size_;
  std::vector<std::uint8_t> values_;
};

class PackedTensorContractions {
 public:
  PackedTensorContractions(std::uint32_t q, std::uint64_t signature)
      : q_(q),
        slices_(static_cast<std::size_t>(q) * q, 0) {
    if (q_ < 1 || q_ > 7) {
      throw std::invalid_argument(
          "packed tensor contractions support q=1,...,7");
    }
    const auto term_count = cubic_tensor_term_count(q_);
    if (signature >= (std::uint64_t{1} << term_count)) {
      throw std::invalid_argument("tensor signature lies outside its q space");
    }
    const auto set = [&](std::uint32_t contracted,
                         std::uint32_t left,
                         std::uint32_t right) {
      slices_[static_cast<std::size_t>(contracted) * q_ + left] ^=
          std::uint32_t{1} << right;
    };
    std::uint32_t bit = 0;
    for (std::uint32_t first = 0; first < q_; ++first, ++bit) {
      if ((signature >> bit) & 1U) set(first, first, first);
    }
    for (std::uint32_t first = 0; first < q_; ++first) {
      for (std::uint32_t second = first + 1; second < q_; ++second, ++bit) {
        if (((signature >> bit) & 1U) == 0) continue;
        set(first, first, second);
        set(first, second, first);
        set(first, second, second);
        set(second, first, first);
        set(second, first, second);
        set(second, second, first);
      }
    }
    for (std::uint32_t first = 0; first < q_; ++first) {
      for (std::uint32_t second = first + 1; second < q_; ++second) {
        for (std::uint32_t third = second + 1; third < q_; ++third, ++bit) {
          if (((signature >> bit) & 1U) == 0) continue;
          set(first, second, third);
          set(first, third, second);
          set(second, first, third);
          set(second, third, first);
          set(third, first, second);
          set(third, second, first);
        }
      }
    }
    if (bit != term_count) {
      throw std::logic_error("packed contraction decoded the wrong term count");
    }
  }

  [[nodiscard]] std::vector<std::uint32_t> matrix(
      std::uint32_t contracted_vector) const {
    std::vector<std::uint32_t> result(q_, 0);
    while (contracted_vector != 0) {
      const auto contracted = static_cast<std::uint32_t>(
          std::countr_zero(contracted_vector));
      for (std::uint32_t left = 0; left < q_; ++left) {
        result[left] ^=
            slices_[static_cast<std::size_t>(contracted) * q_ + left];
      }
      contracted_vector &= contracted_vector - 1;
    }
    return result;
  }

  [[nodiscard]] bool vanishes_on(
      std::uint32_t contracted_vector,
      std::span<const std::uint32_t> subspace_basis) const {
    const auto rows = matrix(contracted_vector);
    for (const auto left : subspace_basis) {
      std::uint32_t contraction = 0;
      auto coefficients = left;
      while (coefficients != 0) {
        const auto bit = static_cast<std::uint32_t>(
            std::countr_zero(coefficients));
        contraction ^= rows[bit];
        coefficients &= coefficients - 1;
      }
      for (const auto right : subspace_basis) {
        if (std::popcount(contraction & right) & 1U) return false;
      }
    }
    return true;
  }

 private:
  std::uint32_t q_;
  std::vector<std::uint32_t> slices_;
};

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

template <typename Values>
[[nodiscard]] std::vector<std::uint8_t> extension_block(
    const Values& values,
    const std::vector<std::uint32_t>& basis,
    std::uint32_t candidate) {
  std::vector<std::uint8_t> block;
  block.reserve(1 + basis.size() + basis.size() * (basis.size() - 1) / 2);
  block.push_back(values(candidate, candidate, candidate));
  for (const auto prior : basis) {
    block.push_back(values(prior, prior, candidate));
  }
  for (std::size_t first = 0; first < basis.size(); ++first) {
    for (std::size_t second = first + 1; second < basis.size(); ++second) {
      block.push_back(values(basis[first], basis[second], candidate));
    }
  }
  return block;
}

[[nodiscard]] int compare_prefix(
    const std::vector<std::uint8_t>& prefix,
    const std::vector<std::uint8_t>& best) {
  for (std::size_t index = 0; index < prefix.size(); ++index) {
    if (prefix[index] < best[index]) return -1;
    if (prefix[index] > best[index]) return 1;
  }
  return 0;
}

struct SpanMask {
  std::array<std::uint64_t, 4> words{};

  [[nodiscard]] bool contains(std::uint32_t point) const {
    return ((words[point / 64] >> (point % 64)) & 1U) != 0;
  }

  void insert(std::uint32_t point) {
    words[point / 64] |= std::uint64_t{1} << (point % 64);
  }
};

[[nodiscard]] SpanMask extend_span(
    SpanMask span,
    std::uint32_t vector,
    std::uint32_t vector_space_size) {
  SpanMask result = span;
  for (std::uint32_t point = 0; point < vector_space_size; ++point) {
    if (span.contains(point)) result.insert(point ^ vector);
  }
  return result;
}

[[nodiscard]] std::uint64_t pack_key(std::span<const std::uint8_t> bits) {
  std::uint64_t result = 0;
  for (std::size_t index = 0; index < bits.size(); ++index) {
    result |= std::uint64_t{bits[index]} << index;
  }
  return result;
}

[[nodiscard]] std::vector<std::uint64_t> pack_key_words(
    std::span<const std::uint8_t> bits) {
  std::vector<std::uint64_t> words((bits.size() + 63) / 64, 0);
  for (std::size_t index = 0; index < bits.size(); ++index) {
    words[index / 64] |= std::uint64_t{bits[index]} << (index % 64);
  }
  return words;
}

template <typename Values>
[[nodiscard]] DirectTensorCanonicalForm canonicalize_values(
    std::uint32_t logical_qubits,
    const Values& values) {
  const auto term_count = cubic_tensor_term_count(logical_qubits);
  const std::uint32_t vector_space_size = std::uint32_t{1} << logical_qubits;

  DirectTensorCanonicalForm result;
  result.logical_qubits = logical_qubits;
  std::vector<std::uint8_t> best;
  std::vector<std::uint8_t> prefix;
  std::vector<std::uint32_t> basis;
  std::function<void(SpanMask)> extend = [&](SpanMask span) {
    if (basis.size() == logical_qubits) {
      ++result.complete_bases_evaluated;
      if (best.empty() || prefix < best) {
        best = prefix;
        result.new_basis_in_old_coordinates = basis;
      }
      return;
    }
    for (std::uint32_t candidate = 1; candidate < vector_space_size; ++candidate) {
      if (span.contains(candidate)) continue;
      const auto block = extension_block(values, basis, candidate);
      prefix.insert(prefix.end(), block.begin(), block.end());
      if (!best.empty() && compare_prefix(prefix, best) > 0) {
        ++result.prefix_branches_pruned;
      } else {
        basis.push_back(candidate);
        extend(extend_span(span, candidate, vector_space_size));
        basis.pop_back();
      }
      prefix.resize(prefix.size() - block.size());
    }
  };
  SpanMask zero_span;
  zero_span.insert(0);
  extend(zero_span);
  if (best.size() != term_count) {
    throw std::logic_error("direct tensor canonicalization found no basis");
  }
  result.canonical_key_words = pack_key_words(best);
  result.canonical_key = result.canonical_key_words.empty()
                             ? 0
                             : result.canonical_key_words[0];
  return result;
}

[[nodiscard]] std::uint32_t rank32(
    std::span<const std::uint32_t> vectors) {
  std::array<std::uint32_t, 32> pivots{};
  std::uint32_t rank = 0;
  for (auto value : vectors) {
    while (value != 0) {
      const auto pivot =
          31U - static_cast<std::uint32_t>(std::countl_zero(value));
      if (pivots[pivot] != 0) {
        value ^= pivots[pivot];
      } else {
        pivots[pivot] = value;
        ++rank;
        break;
      }
    }
  }
  return rank;
}

[[nodiscard]] std::uint32_t linear_combination32(
    std::uint32_t coefficients,
    std::span<const std::uint32_t> vectors) {
  std::uint32_t result = 0;
  while (coefficients != 0) {
    const auto bit = static_cast<std::uint32_t>(
        std::countr_zero(coefficients));
    result ^= vectors[bit];
    coefficients &= coefficients - 1;
  }
  return result;
}

constexpr std::uint32_t kFiveQubitTensorBits = 25;
constexpr std::uint32_t kFiveQubitTensorCount =
    std::uint32_t{1} << kFiveQubitTensorBits;
constexpr std::uint64_t kGlFiveOrder =
    std::uint64_t{31} * 30 * 28 * 24 * 16;

using FiveBasis = std::array<std::uint32_t, 5>;

constexpr FiveBasis kCycleBasis{2, 4, 8, 16, 1};
constexpr FiveBasis kTransvectionBasis{3, 2, 4, 8, 16};

class PackedSignatureMap {
 public:
  explicit PackedSignatureMap(
      const std::array<std::uint32_t, kFiveQubitTensorBits>& images) {
    for (std::uint32_t chunk = 0; chunk < chunks_.size(); ++chunk) {
      for (std::uint32_t mask = 0; mask < chunks_[chunk].size(); ++mask) {
        std::uint32_t image = 0;
        for (std::uint32_t offset = 0; offset < 5; ++offset) {
          if ((mask >> offset) & 1U) image ^= images[5 * chunk + offset];
        }
        chunks_[chunk][mask] = image;
      }
    }
  }

  [[nodiscard]] std::uint32_t operator()(std::uint32_t signature) const {
    std::uint32_t result = 0;
    for (std::uint32_t chunk = 0; chunk < chunks_.size(); ++chunk) {
      result ^= chunks_[chunk][(signature >> (5 * chunk)) & 31U];
    }
    return result;
  }

 private:
  std::array<std::array<std::uint32_t, 32>, 5> chunks_{};
};

[[nodiscard]] PackedSignatureMap tensor_action(const FiveBasis& basis) {
  std::array<std::uint32_t, kFiveQubitTensorBits> images{};
  for (std::uint32_t bit = 0; bit < images.size(); ++bit) {
    images[bit] = static_cast<std::uint32_t>(
        transform_cubic_tensor_signature(5, std::uint64_t{1} << bit, basis));
  }
  return PackedSignatureMap(images);
}

[[nodiscard]] PackedSignatureMap standard_to_interleaved_map() {
  constexpr FiveBasis identity{1, 2, 4, 8, 16};
  std::array<std::uint32_t, kFiveQubitTensorBits> images{};
  const std::vector<std::uint32_t> basis(identity.begin(), identity.end());
  for (std::uint32_t bit = 0; bit < images.size(); ++bit) {
    images[bit] = static_cast<std::uint32_t>(
        interleaved_tensor_key(5, std::uint64_t{1} << bit, basis));
  }
  return PackedSignatureMap(images);
}

[[nodiscard]] bool mark_state(
    std::vector<std::uint64_t>& seen, std::uint32_t state) {
  auto& word = seen[state / 64];
  const auto mask = std::uint64_t{1} << (state % 64);
  if ((word & mask) != 0) return false;
  word |= mask;
  return true;
}

[[nodiscard]] bool interleaved_key_less(
    std::uint32_t left, std::uint32_t right) {
  const auto difference = left ^ right;
  if (difference == 0) return false;
  const auto first_difference = static_cast<std::uint32_t>(
      std::countr_zero(difference));
  return ((left >> first_difference) & 1U) == 0;
}

[[nodiscard]] std::uint32_t packed_basis(const FiveBasis& basis) {
  std::uint32_t result = 0;
  for (std::uint32_t index = 0; index < basis.size(); ++index) {
    result |= basis[index] << (5 * index);
  }
  return result;
}

[[nodiscard]] FiveBasis unpacked_basis(std::uint32_t packed) {
  FiveBasis basis{};
  for (std::uint32_t index = 0; index < basis.size(); ++index) {
    basis[index] = (packed >> (5 * index)) & 31U;
  }
  return basis;
}

[[nodiscard]] std::uint32_t compose_packed_basis(
    std::uint32_t packed, const FiveBasis& generator) {
  const auto basis = unpacked_basis(packed);
  FiveBasis transformed{};
  for (std::uint32_t target = 0; target < transformed.size(); ++target) {
    for (std::uint32_t source = 0; source < basis.size(); ++source) {
      if ((generator[target] >> source) & 1U) {
        transformed[target] ^= basis[source];
      }
    }
  }
  return packed_basis(transformed);
}

[[nodiscard]] std::uint64_t generated_five_qubit_linear_group_order() {
  constexpr FiveBasis identity{1, 2, 4, 8, 16};
  const std::array<FiveBasis, 2> generators{
      kCycleBasis, kTransvectionBasis};
  std::vector<std::uint64_t> seen(kFiveQubitTensorCount / 64, 0);
  std::vector<std::uint32_t> queue;
  queue.reserve(static_cast<std::size_t>(kGlFiveOrder));
  const auto start = packed_basis(identity);
  if (!mark_state(seen, start)) {
    throw std::logic_error("identity basis was unexpectedly already visited");
  }
  queue.push_back(start);
  for (std::size_t cursor = 0; cursor < queue.size(); ++cursor) {
    for (const auto& generator : generators) {
      const auto next = compose_packed_basis(queue[cursor], generator);
      if (mark_state(seen, next)) queue.push_back(next);
    }
  }
  return queue.size();
}

}  // namespace

std::uint64_t restrict_cubic_tensor_signature(
    std::uint32_t source_logical_qubits,
    std::uint64_t standard_tensor_signature,
    std::span<const std::uint32_t> basis) {
  if (source_logical_qubits < 1 || source_logical_qubits > 7 ||
      basis.empty() || basis.size() > source_logical_qubits ||
      rank32(basis) != basis.size()) {
    throw std::invalid_argument(
        "tensor restriction requires an independent nonempty basis");
  }
  const auto vector_space_size =
      std::uint32_t{1} << source_logical_qubits;
  if (std::any_of(
          basis.begin(), basis.end(),
          [vector_space_size](std::uint32_t vector) {
            return vector >= vector_space_size;
          })) {
    throw std::invalid_argument("tensor basis vector lies outside its q space");
  }
  const auto source_terms = cubic_tensor_term_count(source_logical_qubits);
  if (standard_tensor_signature >= (std::uint64_t{1} << source_terms)) {
    throw std::invalid_argument("tensor signature lies outside its q space");
  }
  const auto target_q = static_cast<std::uint32_t>(basis.size());
  const auto target_terms = cubic_tensor_term_count(target_q);
  std::uint64_t transformed = 0;
  std::uint32_t bit = 0;
  for (std::uint32_t first = 0; first < target_q; ++first, ++bit) {
    if (cubic_tensor_value(
            source_logical_qubits,
            standard_tensor_signature,
            basis[first],
            basis[first],
            basis[first])) {
      transformed |= std::uint64_t{1} << bit;
    }
  }
  for (std::uint32_t first = 0; first < target_q; ++first) {
    for (std::uint32_t second = first + 1; second < target_q;
         ++second, ++bit) {
      const bool left = cubic_tensor_value(
          source_logical_qubits,
          standard_tensor_signature,
          basis[first],
          basis[first],
          basis[second]);
      const bool right = cubic_tensor_value(
          source_logical_qubits,
          standard_tensor_signature,
          basis[first],
          basis[second],
          basis[second]);
      if (left != right) {
        throw std::logic_error(
            "basis change left the binary symmetric cubic-tensor space");
      }
      if (left) transformed |= std::uint64_t{1} << bit;
    }
  }
  for (std::uint32_t first = 0; first < target_q; ++first) {
    for (std::uint32_t second = first + 1; second < target_q; ++second) {
      for (std::uint32_t third = second + 1; third < target_q;
           ++third, ++bit) {
        if (cubic_tensor_value(
                source_logical_qubits,
                standard_tensor_signature,
                basis[first],
                basis[second],
                basis[third])) {
          transformed |= std::uint64_t{1} << bit;
        }
      }
    }
  }
  if (bit != target_terms) {
    throw std::logic_error("tensor transformation emitted the wrong bit count");
  }
  return transformed;
}

std::uint64_t transform_cubic_tensor_signature(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature,
    std::span<const std::uint32_t> basis) {
  if (basis.size() != logical_qubits) {
    throw std::invalid_argument(
        "tensor transformation requires a full basis");
  }
  return restrict_cubic_tensor_signature(
      logical_qubits, standard_tensor_signature, basis);
}

std::vector<std::uint32_t> cubic_tensor_radical_basis_from_signature(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature) {
  if (logical_qubits < 1 || logical_qubits > 7) {
    throw std::invalid_argument("packed tensor radical supports q=1,...,7");
  }
  const PackedTensorContractions contractions(
      logical_qubits, standard_tensor_signature);
  const auto vector_space_size = std::uint32_t{1} << logical_qubits;
  std::array<std::uint32_t, 32> pivots{};
  std::vector<std::uint32_t> basis;
  for (std::uint32_t candidate = 1; candidate < vector_space_size;
       ++candidate) {
    const auto matrix = contractions.matrix(candidate);
    if (std::any_of(matrix.begin(), matrix.end(),
                    [](std::uint32_t row) { return row != 0; })) {
      continue;
    }
    auto reduced = candidate;
    while (reduced != 0) {
      const auto pivot =
          31U - static_cast<std::uint32_t>(std::countl_zero(reduced));
      if (pivots[pivot] != 0) {
        reduced ^= pivots[pivot];
      } else {
        pivots[pivot] = reduced;
        basis.push_back(candidate);
        break;
      }
    }
  }
  return basis;
}

std::uint64_t nondegenerate_tensor_hyperplane_count_from_signature(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature) {
  if (logical_qubits < 2 || logical_qubits > 7) {
    throw std::invalid_argument(
        "packed tensor hyperplane count supports q=2,...,7");
  }
  const PackedTensorContractions contractions(
      logical_qubits, standard_tensor_signature);
  const auto vector_space_size = std::uint32_t{1} << logical_qubits;
  std::uint64_t nondegenerate = 0;
  for (std::uint32_t functional = 1; functional < vector_space_size;
       ++functional) {
    const auto pivot = static_cast<std::uint32_t>(
        std::countr_zero(functional));
    std::vector<std::uint32_t> kernel_basis;
    kernel_basis.reserve(logical_qubits - 1);
    for (std::uint32_t coordinate = 0; coordinate < logical_qubits;
         ++coordinate) {
      if (coordinate == pivot) continue;
      auto vector = std::uint32_t{1} << coordinate;
      if ((functional >> coordinate) & 1U) {
        vector |= std::uint32_t{1} << pivot;
      }
      kernel_basis.push_back(vector);
    }

    bool has_radical = false;
    const auto kernel_size = std::uint32_t{1} << kernel_basis.size();
    for (std::uint32_t coefficients = 1;
         coefficients < kernel_size && !has_radical; ++coefficients) {
      const auto candidate =
          linear_combination32(coefficients, kernel_basis);
      has_radical = contractions.vanishes_on(candidate, kernel_basis);
    }
    if (!has_radical) ++nondegenerate;
  }
  return nondegenerate;
}

DirectTensorCanonicalForm canonicalize_cubic_tensor_direct(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature) {
  if (logical_qubits < 1 || logical_qubits > 7) {
    throw std::invalid_argument("direct canonicalizer currently supports q<=7");
  }
  const auto term_count = cubic_tensor_term_count(logical_qubits);
  if (standard_tensor_signature >= (std::uint64_t{1} << term_count)) {
    throw std::invalid_argument("tensor signature lies outside its q space");
  }
  const TensorValues values(logical_qubits, standard_tensor_signature);
  return canonicalize_values(logical_qubits, values);
}

DirectTensorCanonicalForm canonicalize_cubic_tensor_from_labels(
    std::span<const std::uint64_t> label_rows) {
  if (label_rows.empty() || label_rows.size() > 8) {
    throw std::invalid_argument("label-row canonicalizer currently supports q<=8");
  }
  const LabelTensorValues values(label_rows);
  return canonicalize_values(
      static_cast<std::uint32_t>(label_rows.size()), values);
}

std::uint64_t interleaved_tensor_key(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature,
    const std::vector<std::uint32_t>& new_basis_in_old_coordinates) {
  if (new_basis_in_old_coordinates.size() != logical_qubits) {
    throw std::invalid_argument("tensor basis has the wrong length");
  }
  const TensorValues values(logical_qubits, standard_tensor_signature);
  std::vector<std::uint8_t> bits;
  std::vector<std::uint32_t> prefix;
  for (const auto vector : new_basis_in_old_coordinates) {
    const auto block = extension_block(values, prefix, vector);
    bits.insert(bits.end(), block.begin(), block.end());
    prefix.push_back(vector);
  }
  return pack_key(bits);
}

std::vector<std::uint64_t> interleaved_tensor_words_from_labels(
    std::span<const std::uint64_t> label_rows,
    const std::vector<std::uint32_t>& new_basis_in_old_coordinates) {
  if (new_basis_in_old_coordinates.size() != label_rows.size()) {
    throw std::invalid_argument("tensor basis has the wrong length");
  }
  const LabelTensorValues values(label_rows);
  std::vector<std::uint8_t> bits;
  std::vector<std::uint32_t> prefix;
  for (const auto vector : new_basis_in_old_coordinates) {
    const auto block = extension_block(values, prefix, vector);
    bits.insert(bits.end(), block.begin(), block.end());
    prefix.push_back(vector);
  }
  return pack_key_words(bits);
}

FiveQubitTensorOrbitCensus enumerate_five_qubit_tensor_orbits(
    bool check_with_direct_canonicalizer) {
  FiveQubitTensorOrbitCensus result;
  result.expected_linear_group_order = kGlFiveOrder;
  result.generated_linear_group_order =
      generated_five_qubit_linear_group_order();
  if (result.generated_linear_group_order != kGlFiveOrder) {
    throw std::logic_error(
        "five-qubit tensor generators did not generate GL(5,2)");
  }

  const std::array<PackedSignatureMap, 2> generators{
      tensor_action(kCycleBasis), tensor_action(kTransvectionBasis)};
  const auto interleaved = standard_to_interleaved_map();
  std::vector<std::uint64_t> seen(kFiveQubitTensorCount / 64, 0);
  std::vector<std::uint32_t> queue;
  queue.reserve(static_cast<std::size_t>(kGlFiveOrder));

  for (std::uint32_t start = 0; start < kFiveQubitTensorCount; ++start) {
    if (!mark_state(seen, start)) continue;
    queue.clear();
    queue.push_back(start);
    std::uint32_t canonical_key = 0;
    std::uint32_t canonical_signature = 0;
    for (std::size_t cursor = 0; cursor < queue.size(); ++cursor) {
      const auto current = queue[cursor];
      const auto key = interleaved(current);
      if (cursor == 0 || interleaved_key_less(key, canonical_key)) {
        canonical_key = key;
        canonical_signature = current;
      }
      for (const auto& generator : generators) {
        ++result.generator_transitions;
        const auto next = generator(current);
        if (mark_state(seen, next)) queue.push_back(next);
      }
    }

    const auto orbit_size = static_cast<std::uint64_t>(queue.size());
    if (orbit_size == 0 || kGlFiveOrder % orbit_size != 0) {
      throw std::logic_error(
          "five-qubit tensor orbit size does not divide |GL(5,2)|");
    }
    const auto radical_dimension = static_cast<std::uint32_t>(
        cubic_tensor_radical_basis_from_signature(
            5, canonical_signature).size());
    if (radical_dimension > 5) {
      throw std::logic_error("five-qubit tensor radical has invalid dimension");
    }
    if (check_with_direct_canonicalizer) {
      const auto direct =
          canonicalize_cubic_tensor_direct(5, canonical_signature);
      ++result.direct_canonical_checks;
      if (direct.canonical_key != canonical_key ||
          interleaved_tensor_key(
              5,
              canonical_signature,
              direct.new_basis_in_old_coordinates) != canonical_key) {
        ++result.direct_canonical_mismatches;
      }
    }

    ++result.orbit_counts_by_radical_dimension[radical_dimension];
    result.tensor_counts_by_radical_dimension[radical_dimension] += orbit_size;
    result.tensor_count += orbit_size;
    result.orbits.push_back(FiveQubitTensorOrbit{
        canonical_key,
        canonical_signature,
        orbit_size,
        kGlFiveOrder / orbit_size,
        radical_dimension,
    });
  }

  if (result.tensor_count != kFiveQubitTensorCount) {
    throw std::logic_error("five-qubit tensor orbit traversal was incomplete");
  }
  std::sort(
      result.orbits.begin(), result.orbits.end(),
      [](const FiveQubitTensorOrbit& left,
         const FiveQubitTensorOrbit& right) {
        return left.canonical_key < right.canonical_key;
      });
  return result;
}

}  // namespace utsp
