#include "utsp/tensor_primitives.hpp"

#include "utsp/tensor_canonical.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace utsp {
namespace {

struct AffineBitAction {
  std::uint64_t offset = 0;
  std::vector<std::uint64_t> columns;

  [[nodiscard]] std::uint64_t apply(std::uint64_t value) const {
    auto result = offset;
    while (value != 0) {
      const auto bit = static_cast<std::size_t>(std::countr_zero(value));
      result ^= columns[bit];
      value &= value - 1;
    }
    return result;
  }
};

struct Q7AlternatingAuthorityEntry {
  std::uint32_t id;
  std::vector<std::array<std::uint32_t, 3>> triples;
  std::uint64_t stabilizer_order;
};

struct Q8AlternatingAuthorityEntry {
  std::uint32_t id;
  std::vector<std::array<std::uint32_t, 3>> triples;
  std::uint64_t stabilizer_order;
};

[[nodiscard]] std::uint64_t low_mask(std::uint32_t width) {
  if (width > 63) {
    throw std::invalid_argument("packed tensor sector exceeds 63 bits");
  }
  return width == 0 ? 0 : ((std::uint64_t{1} << width) - 1);
}

[[nodiscard]] std::vector<std::uint32_t> identity_basis(std::uint32_t q) {
  std::vector<std::uint32_t> basis(q);
  for (std::uint32_t index = 0; index < q; ++index) {
    basis[index] = std::uint32_t{1} << index;
  }
  return basis;
}

[[nodiscard]] std::uint32_t pair_bit(
    std::uint32_t q,
    std::uint32_t left,
    std::uint32_t right) {
  if (left >= right || right >= q) {
    throw std::invalid_argument("invalid pair coefficient");
  }
  std::uint32_t bit = q;
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second, ++bit) {
      if (first == left && second == right) return bit;
    }
  }
  throw std::logic_error("pair coefficient was not found");
}

[[nodiscard]] std::uint32_t triple_index(
    std::uint32_t q,
    std::uint32_t first,
    std::uint32_t second,
    std::uint32_t third) {
  if (first >= second || second >= third || third >= q) {
    throw std::invalid_argument("invalid triple coefficient");
  }
  std::uint32_t index = 0;
  for (std::uint32_t left = 0; left < q; ++left) {
    for (std::uint32_t middle = left + 1; middle < q; ++middle) {
      for (std::uint32_t right = middle + 1; right < q; ++right, ++index) {
        if (left == first && middle == second && right == third) return index;
      }
    }
  }
  throw std::logic_error("triple coefficient was not found");
}

[[nodiscard]] std::uint32_t triple_bit(
    std::uint32_t q,
    std::uint32_t first,
    std::uint32_t second,
    std::uint32_t third) {
  return cubic_tensor_repeated_term_count(q) +
         triple_index(q, first, second, third);
}

[[nodiscard]] std::uint64_t alternating_signature_from_triples(
    std::uint32_t q,
    std::span<const std::array<std::uint32_t, 3>> triples) {
  std::uint64_t signature = 0;
  for (const auto& triple : triples) {
    signature |= std::uint64_t{1}
                 << triple_bit(q, triple[0], triple[1], triple[2]);
  }
  return signature;
}

[[nodiscard]] std::uint64_t alternating_state_from_triples(
    std::uint32_t q,
    std::span<const std::array<std::uint32_t, 3>> triples) {
  std::uint64_t state = 0;
  for (const auto& triple : triples) {
    state |= std::uint64_t{1}
             << triple_index(q, triple[0], triple[1], triple[2]);
  }
  return state;
}

[[nodiscard]] std::uint64_t rank_one_cube_signature(
    std::uint32_t q,
    std::uint32_t functional) {
  if (functional == 0 || functional >= (std::uint32_t{1} << q)) {
    throw std::invalid_argument("rank-one functional must be nonzero");
  }
  std::uint64_t signature = 0;
  for (std::uint32_t first = 0; first < q; ++first) {
    if ((functional >> first) & 1U) {
      signature |= std::uint64_t{1} << first;
    }
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      if (((functional >> first) & 1U) &&
          ((functional >> second) & 1U)) {
        signature |= std::uint64_t{1} << pair_bit(q, first, second);
      }
    }
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      for (std::uint32_t third = second + 1; third < q; ++third) {
        if (((functional >> first) & 1U) &&
            ((functional >> second) & 1U) &&
            ((functional >> third) & 1U)) {
          signature |= std::uint64_t{1}
                       << triple_bit(q, first, second, third);
        }
      }
    }
  }
  return signature;
}

[[nodiscard]] std::uint32_t binary_rank(
    std::vector<std::uint32_t> rows) {
  std::uint32_t rank = 0;
  for (std::uint32_t column = 0; column < 32; ++column) {
    const auto pivot = std::find_if(
        rows.begin() + rank,
        rows.end(),
        [column](std::uint32_t row) {
          return ((row >> column) & 1U) != 0;
        });
    if (pivot == rows.end()) continue;
    std::iter_swap(rows.begin() + rank, pivot);
    for (std::size_t row = 0; row < rows.size(); ++row) {
      if (row != rank && ((rows[row] >> column) & 1U)) {
        rows[row] ^= rows[rank];
      }
    }
    ++rank;
  }
  return rank;
}

[[nodiscard]] std::vector<std::uint32_t> functional_kernel_basis(
    std::uint32_t q,
    std::uint32_t functional) {
  const auto pivot = static_cast<std::uint32_t>(
      std::countr_zero(functional));
  std::vector<std::uint32_t> basis;
  basis.reserve(q - 1);
  for (std::uint32_t coordinate = 0; coordinate < q; ++coordinate) {
    if (coordinate == pivot) continue;
    auto vector = std::uint32_t{1} << coordinate;
    if ((functional >> coordinate) & 1U) {
      vector |= std::uint32_t{1} << pivot;
    }
    basis.push_back(vector);
  }
  return basis;
}

[[nodiscard]] bool alternating_form_is_q7_family_restriction(
    std::uint32_t q,
    std::uint64_t alternating_signature,
    std::uint32_t functional) {
  const auto repeated = cubic_tensor_repeated_term_count(q);
  if ((alternating_signature & low_mask(repeated)) != 0 || functional == 0) {
    return false;
  }
  const auto kernel = functional_kernel_basis(q, functional);
  for (std::size_t first = 0; first < kernel.size(); ++first) {
    for (std::size_t second = first + 1; second < kernel.size(); ++second) {
      for (std::size_t third = second + 1; third < kernel.size(); ++third) {
        if (cubic_tensor_value(
                q,
                alternating_signature,
                kernel[first],
                kernel[second],
                kernel[third])) {
          return false;
        }
      }
    }
  }

  const auto transverse = std::uint32_t{1}
                          << std::countr_zero(functional);
  std::vector<std::uint32_t> bilinear_rows(kernel.size(), 0);
  for (std::size_t first = 0; first < kernel.size(); ++first) {
    for (std::size_t second = first + 1; second < kernel.size(); ++second) {
      if (!cubic_tensor_value(
              q,
              alternating_signature,
              transverse,
              kernel[first],
              kernel[second])) {
        continue;
      }
      bilinear_rows[first] |= std::uint32_t{1} << second;
      bilinear_rows[second] |= std::uint32_t{1} << first;
    }
  }
  const auto rank = binary_rank(std::move(bilinear_rows));
  const auto kernel_dimension = q - 1;
  const auto radical_dimension = kernel_dimension - rank;
  return radical_dimension <= 6 - kernel_dimension;
}

[[nodiscard]] bool is_q7_primitive_family_restriction(
    std::uint32_t q,
    std::uint64_t signature) {
  if (q == 0 || q > 6) return false;
  const auto term_count = cubic_tensor_term_count(q);
  if (signature >= (std::uint64_t{1} << term_count)) return false;

  std::vector<std::uint32_t> gram_rows(q, 0);
  for (std::uint32_t coordinate = 0; coordinate < q; ++coordinate) {
    if ((signature >> coordinate) & 1U) {
      gram_rows[coordinate] |= std::uint32_t{1} << coordinate;
    }
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      if ((signature >> pair_bit(q, first, second)) & 1U) {
        gram_rows[first] |= std::uint32_t{1} << second;
        gram_rows[second] |= std::uint32_t{1} << first;
      }
    }
  }
  const auto gram_rank = binary_rank(gram_rows);
  if (gram_rank == 0) {
    if (signature == 0) return true;
    for (std::uint32_t functional = 1;
         functional < (std::uint32_t{1} << q); ++functional) {
      if (alternating_form_is_q7_family_restriction(
              q, signature, functional)) {
        return true;
      }
    }
    return false;
  }
  if (gram_rank != 1) return false;

  std::uint32_t functional = 0;
  for (std::uint32_t coordinate = 0; coordinate < q; ++coordinate) {
    if ((signature >> coordinate) & 1U) {
      functional |= std::uint32_t{1} << coordinate;
    }
  }
  if (functional == 0) return false;
  for (std::uint32_t coordinate = 0; coordinate < q; ++coordinate) {
    const auto expected = ((functional >> coordinate) & 1U)
                              ? functional
                              : std::uint32_t{0};
    if (gram_rows[coordinate] != expected) return false;
  }
  return alternating_form_is_q7_family_restriction(
      q,
      signature ^ rank_one_cube_signature(q, functional),
      functional);
}

[[nodiscard]] const std::vector<Q7AlternatingAuthorityEntry>&
q7_alternating_authority_entries() {
  // Coordinates are zero-based transcriptions of the published 1-based
  // representatives.  Entry zero is the zero form.
  static const std::vector<Q7AlternatingAuthorityEntry> entries{
      {0, {}, 163849992929280ULL},
      {1, {{{0, 1, 2}}}, 13872660480ULL},
      {2, {{{0, 1, 2}}, {{2, 3, 4}}}, 70778880ULL},
      {3, {{{0, 1, 2}}, {{3, 4, 5}}}, 3612672ULL},
      {4, {{{0, 1, 2}}, {{2, 3, 4}}, {{0, 4, 5}}}, 2752512ULL},
      {5,
       {{{0, 1, 2}}, {{1, 2, 3}}, {{2, 3, 4}}, {{1, 3, 5}}, {{0, 4, 5}}},
       7741440ULL},
      {6, {{{0, 1, 2}}, {{2, 3, 4}}, {{4, 5, 6}}}, 73728ULL},
      {7, {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 6}}, {{2, 4, 6}}}, 688128ULL},
      {8, {{{0, 1, 2}}, {{0, 5, 6}}, {{1, 3, 5}}, {{2, 4, 6}}}, 9216ULL},
      {9, {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 6}}}, 92897280ULL},
      {10,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 6}}, {{1, 3, 5}}, {{2, 4, 6}}},
       12096ULL},
      {11,
       {{{0, 1, 2}}, {{1, 2, 3}}, {{2, 3, 4}}, {{1, 3, 5}},
        {{0, 4, 5}}, {{2, 5, 6}}},
       368640ULL},
  };
  return entries;
}

[[nodiscard]] const std::vector<Q8AlternatingAuthorityEntry>&
q8_alternating_authority_entries() {
  // Hora--Pudlak numbering, independently tabulated with stabilizer orders by
  // O'Brien--Vojtechovsky.  Coordinates below are zero-based.
  static const std::vector<Q8AlternatingAuthorityEntry> entries{
      {0, {}, 5348063769211699200ULL},
      {1, {{{0, 1, 2}}}, 55046716784640ULL},
      {2, {{{0, 1, 2}}, {{2, 3, 4}}}, 63417876480ULL},
      {3, {{{0, 1, 2}}, {{3, 4, 5}}}, 1387266048ULL},
      {4, {{{0, 1, 2}}, {{2, 3, 4}}, {{0, 4, 5}}}, 1056964608ULL},
      {5,
       {{{0, 1, 2}}, {{1, 2, 3}}, {{2, 3, 4}}, {{1, 3, 5}}, {{0, 4, 5}}},
       2972712960ULL},
      {6, {{{0, 1, 2}}, {{2, 3, 4}}, {{4, 5, 6}}}, 9437184ULL},
      {7, {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 6}}, {{2, 4, 6}}}, 88080384ULL},
      {8, {{{0, 1, 2}}, {{0, 5, 6}}, {{1, 3, 5}}, {{2, 4, 6}}}, 1179648ULL},
      {9, {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 6}}}, 11890851840ULL},
      {10,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 6}}, {{1, 3, 5}}, {{2, 4, 6}}},
       1548288ULL},
      {11,
       {{{0, 1, 2}}, {{1, 2, 3}}, {{2, 3, 4}}, {{1, 3, 5}},
        {{0, 4, 5}}, {{2, 5, 6}}},
       47185920ULL},
      {12, {{{0, 1, 2}}, {{2, 3, 4}}, {{5, 6, 7}}}, 1935360ULL},
      {13, {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 6, 7}}, {{1, 3, 5}}}, 2359296ULL},
      {14, {{{0, 1, 2}}, {{0, 3, 4}}, {{1, 5, 7}}, {{2, 3, 6}}}, 18432ULL},
      {15,
       {{{0, 1, 2}}, {{2, 3, 4}}, {{4, 5, 6}}, {{0, 6, 7}}},
       294912ULL},
      {16,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 7}}, {{1, 3, 5}}, {{1, 4, 6}}},
       393216ULL},
      {17,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 7}}, {{2, 3, 6}}, {{1, 4, 5}}},
       3072ULL},
      {18,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 7}}, {{2, 3, 6}}, {{1, 5, 6}}},
       49152ULL},
      {19,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{1, 6, 7}}, {{2, 4, 5}}, {{3, 5, 6}}},
       2304ULL},
      {20,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 6, 7}}, {{1, 3, 5}},
        {{1, 4, 7}}, {{2, 3, 6}}},
       24576ULL},
      {21,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 7}}, {{2, 3, 6}},
        {{1, 4, 7}}, {{1, 5, 6}}},
       192ULL},
      {22,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{1, 4, 6}}, {{1, 6, 7}},
        {{2, 5, 7}}, {{3, 5, 6}}},
       336ULL},
      {23,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 7}}, {{1, 3, 5}},
        {{1, 4, 6}}, {{2, 4, 5}}, {{3, 4, 5}}},
       1474560ULL},
      {24,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{1, 4, 6}}, {{1, 4, 7}},
        {{1, 5, 7}}, {{2, 3, 7}}, {{3, 5, 6}}},
       96768ULL},
      {25,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 6}}, {{0, 6, 7}},
        {{1, 4, 7}}, {{1, 5, 6}}, {{2, 3, 6}}, {{2, 4, 5}}},
       10752ULL},
      {26,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 6, 7}}, {{1, 3, 5}},
        {{1, 4, 7}}, {{2, 3, 6}}, {{2, 4, 5}}, {{3, 4, 5}}},
       432ULL},
      {27,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{1, 4, 7}}, {{2, 4, 5}},
        {{3, 6, 7}}, {{4, 5, 6}}},
       768ULL},
      {28,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{0, 5, 6}}, {{1, 3, 5}},
        {{1, 4, 6}}, {{1, 5, 6}}, {{2, 5, 7}}},
       3072ULL},
      {29,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{1, 3, 5}}, {{3, 5, 7}}, {{4, 6, 7}}},
       3072ULL},
      {30,
       {{{0, 1, 2}}, {{0, 3, 4}}, {{1, 4, 7}}, {{2, 3, 6}},
        {{2, 5, 7}}, {{4, 5, 6}}},
       46080ULL},
      {31, {{{0, 1, 2}}, {{0, 3, 4}}, {{3, 4, 6}}, {{5, 6, 7}}}, 82944ULL},
  };
  return entries;
}

[[nodiscard]] bool q8_alternating_value(
    std::uint64_t alternating_state,
    std::uint32_t left,
    std::uint32_t middle,
    std::uint32_t right) {
  std::uint32_t bit = 0;
  bool value = false;
  for (std::uint32_t first = 0; first < 8; ++first) {
    for (std::uint32_t second = first + 1; second < 8; ++second) {
      for (std::uint32_t third = second + 1; third < 8; ++third, ++bit) {
        if (((alternating_state >> bit) & 1U) == 0) continue;
        value ^= ((left >> first) & 1U) &&
                 ((middle >> second) & 1U) &&
                 ((right >> third) & 1U);
        value ^= ((left >> first) & 1U) &&
                 ((middle >> third) & 1U) &&
                 ((right >> second) & 1U);
        value ^= ((left >> second) & 1U) &&
                 ((middle >> first) & 1U) &&
                 ((right >> third) & 1U);
        value ^= ((left >> second) & 1U) &&
                 ((middle >> third) & 1U) &&
                 ((right >> first) & 1U);
        value ^= ((left >> third) & 1U) &&
                 ((middle >> first) & 1U) &&
                 ((right >> second) & 1U);
        value ^= ((left >> third) & 1U) &&
                 ((middle >> second) & 1U) &&
                 ((right >> first) & 1U);
      }
    }
  }
  return value;
}

[[nodiscard]] bool q8_pair_tensor_value(
    std::uint32_t first,
    std::uint32_t second,
    std::uint32_t left,
    std::uint32_t middle,
    std::uint32_t right) {
  const bool lf = (left >> first) & 1U;
  const bool ls = (left >> second) & 1U;
  const bool mf = (middle >> first) & 1U;
  const bool ms = (middle >> second) & 1U;
  const bool rf = (right >> first) & 1U;
  const bool rs = (right >> second) & 1U;
  return (lf && mf && rs) ^ (lf && ms && rf) ^ (ls && mf && rf) ^
         (lf && ms && rs) ^ (ls && mf && rs) ^ (ls && ms && rf);
}

[[nodiscard]] bool q8_exceptional_tensor_value(
    PrimitiveGramSector sector,
    std::uint64_t alternating_state,
    std::uint32_t rank_one_functional,
    std::uint32_t left,
    std::uint32_t middle,
    std::uint32_t right) {
  bool value =
      q8_alternating_value(alternating_state, left, middle, right);
  if (sector == PrimitiveGramSector::rank_one) {
    value ^= (std::popcount(left & rank_one_functional) & 1U) &&
             (std::popcount(middle & rank_one_functional) & 1U) &&
             (std::popcount(right & rank_one_functional) & 1U);
  } else if (sector == PrimitiveGramSector::full_alternating) {
    for (std::uint32_t pair = 0; pair < 4; ++pair) {
      value ^= q8_pair_tensor_value(
          2 * pair, 2 * pair + 1, left, middle, right);
    }
  }
  return value;
}

class Q8TensorContractions {
 public:
  Q8TensorContractions(
      PrimitiveGramSector sector,
      std::uint64_t alternating_state,
      std::uint32_t rank_one_functional)
      : slices_(64, 0) {
    for (std::uint32_t contracted = 0; contracted < 8; ++contracted) {
      for (std::uint32_t left = 0; left < 8; ++left) {
        std::uint32_t row = 0;
        for (std::uint32_t right = 0; right < 8; ++right) {
          if (q8_exceptional_tensor_value(
                  sector,
                  alternating_state,
                  rank_one_functional,
                  std::uint32_t{1} << contracted,
                  std::uint32_t{1} << left,
                  std::uint32_t{1} << right)) {
            row |= std::uint32_t{1} << right;
          }
        }
        slices_[8 * contracted + left] = row;
      }
    }
  }

  [[nodiscard]] std::array<std::uint32_t, 8> matrix(
      std::uint32_t contracted) const {
    std::array<std::uint32_t, 8> result{};
    while (contracted != 0) {
      const auto bit = static_cast<std::uint32_t>(
          std::countr_zero(contracted));
      for (std::uint32_t row = 0; row < 8; ++row) {
        result[row] ^= slices_[8 * bit + row];
      }
      contracted &= contracted - 1;
    }
    return result;
  }

  [[nodiscard]] bool vanishes_on(
      std::uint32_t contracted,
      std::span<const std::uint32_t> basis) const {
    const auto rows = matrix(contracted);
    for (const auto left : basis) {
      std::uint32_t contraction = 0;
      auto coefficients = left;
      while (coefficients != 0) {
        const auto bit = static_cast<std::uint32_t>(
            std::countr_zero(coefficients));
        contraction ^= rows[bit];
        coefficients &= coefficients - 1;
      }
      for (const auto right : basis) {
        if (std::popcount(contraction & right) & 1U) return false;
      }
    }
    return true;
  }

 private:
  std::vector<std::uint32_t> slices_;
};

[[nodiscard]] std::uint32_t q8_tensor_radical_dimension(
    const Q8TensorContractions& contractions) {
  std::uint32_t radical_size = 1;
  for (std::uint32_t vector = 1; vector < 256; ++vector) {
    const auto matrix = contractions.matrix(vector);
    if (std::all_of(
            matrix.begin(), matrix.end(),
            [](std::uint32_t row) { return row == 0; })) {
      ++radical_size;
    }
  }
  if (!std::has_single_bit(radical_size)) {
    throw std::logic_error("q8 tensor radical is not a binary subspace");
  }
  return static_cast<std::uint32_t>(std::countr_zero(radical_size));
}

[[nodiscard]] bool q8_has_nondegenerate_hyperplane(
    const Q8TensorContractions& contractions) {
  for (std::uint32_t functional = 1; functional < 256; ++functional) {
    const auto basis = functional_kernel_basis(8, functional);
    bool has_radical = false;
    for (std::uint32_t coefficients = 1;
         coefficients < 128 && !has_radical; ++coefficients) {
      std::uint32_t vector = 0;
      auto remaining = coefficients;
      while (remaining != 0) {
        const auto bit = static_cast<std::uint32_t>(
            std::countr_zero(remaining));
        vector ^= basis[bit];
        remaining &= remaining - 1;
      }
      has_radical = contractions.vanishes_on(vector, basis);
    }
    if (!has_radical) return true;
  }
  return false;
}

[[nodiscard]] std::uint64_t repeated_signature(
    std::uint32_t q,
    PrimitiveGramSector sector) {
  switch (sector) {
    case PrimitiveGramSector::zero:
      return 0;
    case PrimitiveGramSector::rank_one:
      return 1;
    case PrimitiveGramSector::full_alternating: {
      if ((q & 1U) != 0) {
        throw std::invalid_argument(
            "a full alternating Gram form requires even dimension");
      }
      std::uint64_t result = 0;
      for (std::uint32_t pair = 0; pair < q / 2; ++pair) {
        result |= std::uint64_t{1}
                  << pair_bit(q, 2 * pair, 2 * pair + 1);
      }
      return result;
    }
  }
  throw std::logic_error("unknown primitive Gram sector");
}

[[nodiscard]] std::uint64_t general_linear_group_order(std::uint32_t q) {
  std::uint64_t result = 1;
  const auto space = std::uint64_t{1} << q;
  for (std::uint32_t dimension = 0; dimension < q; ++dimension) {
    const auto factor = space - (std::uint64_t{1} << dimension);
    if (result > std::numeric_limits<std::uint64_t>::max() / factor) {
      throw std::overflow_error("GL(q,2) order does not fit in 64 bits");
    }
    result *= factor;
  }
  return result;
}

[[nodiscard]] std::uint64_t symplectic_group_order(std::uint32_t q) {
  if ((q & 1U) != 0) {
    throw std::invalid_argument("symplectic group requires even dimension");
  }
  const auto half = q / 2;
  std::uint64_t result = std::uint64_t{1} << (half * half);
  for (std::uint32_t index = 1; index <= half; ++index) {
    const auto factor = (std::uint64_t{1} << (2 * index)) - 1;
    if (result > std::numeric_limits<std::uint64_t>::max() / factor) {
      throw std::overflow_error("Sp(q,2) order does not fit in 64 bits");
    }
    result *= factor;
  }
  return result;
}

[[nodiscard]] bool symplectic_pair(
    std::uint32_t left,
    std::uint32_t right,
    std::uint32_t q) {
  bool value = false;
  for (std::uint32_t pair = 0; pair < q / 2; ++pair) {
    value ^= ((left >> (2 * pair)) & 1U) &&
             ((right >> (2 * pair + 1)) & 1U);
    value ^= ((left >> (2 * pair + 1)) & 1U) &&
             ((right >> (2 * pair)) & 1U);
  }
  return value;
}

[[nodiscard]] std::vector<std::uint32_t> symplectic_transvection_basis(
    std::uint32_t q,
    std::uint32_t vector) {
  auto basis = identity_basis(q);
  for (std::uint32_t coordinate = 0; coordinate < q; ++coordinate) {
    if (symplectic_pair(
            std::uint32_t{1} << coordinate, vector, q)) {
      basis[coordinate] ^= vector;
    }
  }
  return basis;
}

[[nodiscard]] std::vector<std::uint32_t> symplectic_generator_labels(
    std::uint32_t q) {
  std::vector<std::uint32_t> labels;
  labels.reserve(q + q / 2 - 1);
  for (std::uint32_t coordinate = 0; coordinate < q; ++coordinate) {
    labels.push_back(std::uint32_t{1} << coordinate);
  }
  for (std::uint32_t pair = 0; pair + 1 < q / 2; ++pair) {
    labels.push_back(
        (std::uint32_t{1} << (2 * pair + 1)) |
        (std::uint32_t{1} << (2 * pair + 2)));
  }
  return labels;
}

[[nodiscard]] std::uint32_t generated_transvection_label_count(
    std::uint32_t q,
    std::span<const std::uint32_t> generators) {
  const auto space = std::uint32_t{1} << q;
  std::vector<bool> generated(space, false);
  for (const auto vector : generators) generated[vector] = true;
  bool changed = true;
  while (changed) {
    changed = false;
    for (std::uint32_t left = 1; left < space; ++left) {
      if (!generated[left]) continue;
      for (std::uint32_t right = left + 1; right < space; ++right) {
        if (!generated[right] || !symplectic_pair(left, right, q)) continue;
        const auto sum = left ^ right;
        if (!generated[sum]) {
          generated[sum] = true;
          changed = true;
        }
      }
    }
  }
  return static_cast<std::uint32_t>(
      std::count(generated.begin(), generated.end(), true));
}

[[nodiscard]] std::vector<std::vector<std::uint32_t>> gram_stabilizer_generators(
    std::uint32_t q,
    PrimitiveGramSector sector,
    std::uint32_t& transvection_label_count) {
  std::vector<std::vector<std::uint32_t>> result;
  transvection_label_count = 0;
  if (sector == PrimitiveGramSector::zero) {
    for (std::uint32_t index = 0; index + 1 < q; ++index) {
      auto forward = identity_basis(q);
      forward[index] ^= std::uint32_t{1} << (index + 1);
      result.push_back(std::move(forward));
      auto backward = identity_basis(q);
      backward[index + 1] ^= std::uint32_t{1} << index;
      result.push_back(std::move(backward));
    }
    return result;
  }
  if (sector == PrimitiveGramSector::rank_one) {
    if (q > 1) {
      auto shear = identity_basis(q);
      shear[0] ^= std::uint32_t{1} << 1;
      result.push_back(std::move(shear));
    }
    for (std::uint32_t index = 1; index + 1 < q; ++index) {
      auto forward = identity_basis(q);
      forward[index] ^= std::uint32_t{1} << (index + 1);
      result.push_back(std::move(forward));
      auto backward = identity_basis(q);
      backward[index + 1] ^= std::uint32_t{1} << index;
      result.push_back(std::move(backward));
    }
    return result;
  }
  const auto labels = symplectic_generator_labels(q);
  transvection_label_count =
      generated_transvection_label_count(q, labels);
  if (transvection_label_count != (std::uint32_t{1} << q) - 1) {
    throw std::logic_error(
        "symplectic generators do not generate every transvection label");
  }
  result.reserve(labels.size());
  for (const auto vector : labels) {
    result.push_back(symplectic_transvection_basis(q, vector));
  }
  return result;
}

[[nodiscard]] AffineBitAction induce_alternating_action(
    std::uint32_t q,
    std::uint64_t base_signature,
    std::span<const std::uint32_t> basis) {
  const auto repeated = cubic_tensor_repeated_term_count(q);
  const auto alternating = cubic_tensor_alternating_term_count(q);
  const auto repeated_mask = low_mask(repeated);
  const auto alternating_mask = low_mask(alternating);
  const auto base_repeated = base_signature & repeated_mask;
  const auto transformed_base =
      transform_cubic_tensor_signature(q, base_signature, basis);
  if ((transformed_base & repeated_mask) != base_repeated) {
    throw std::logic_error("Gram-stabilizer generator changed the Gram form");
  }
  AffineBitAction action;
  action.offset = (transformed_base >> repeated) & alternating_mask;
  action.columns.reserve(alternating);
  for (std::uint32_t bit = 0; bit < alternating; ++bit) {
    const auto signature =
        base_signature | (std::uint64_t{1} << (repeated + bit));
    const auto transformed =
        transform_cubic_tensor_signature(q, signature, basis);
    if ((transformed & repeated_mask) != base_repeated) {
      throw std::logic_error("alternating basis tensor changed the Gram form");
    }
    action.columns.push_back(
        ((transformed >> repeated) & alternating_mask) ^ action.offset);
  }
  return action;
}

[[nodiscard]] std::uint64_t gram_stabilizer_order(
    std::uint32_t q,
    PrimitiveGramSector sector) {
  switch (sector) {
    case PrimitiveGramSector::zero:
      return general_linear_group_order(q);
    case PrimitiveGramSector::rank_one:
      return (std::uint64_t{1} << (q - 1)) *
             general_linear_group_order(q - 1);
    case PrimitiveGramSector::full_alternating:
      return symplectic_group_order(q);
  }
  throw std::logic_error("unknown primitive Gram sector");
}

void enumerate_pivot_sets(
    std::uint32_t width,
    std::uint32_t dimension,
    std::uint32_t next,
    std::vector<std::uint32_t>& pivots,
    const std::function<void(std::span<const std::uint32_t>)>& visitor) {
  if (pivots.size() == dimension) {
    visitor(pivots);
    return;
  }
  const auto remaining = dimension - static_cast<std::uint32_t>(pivots.size());
  for (std::uint32_t pivot = next; pivot + remaining <= width; ++pivot) {
    pivots.push_back(pivot);
    enumerate_pivot_sets(
        width, dimension, pivot + 1, pivots, visitor);
    pivots.pop_back();
  }
}

void enumerate_rref_subspaces(
    std::uint32_t width,
    std::uint32_t dimension,
    const std::function<void(std::span<const std::uint32_t>)>& visitor) {
  std::vector<std::uint32_t> pivots;
  enumerate_pivot_sets(
      width,
      dimension,
      0,
      pivots,
      [&](std::span<const std::uint32_t> selected_pivots) {
        std::uint32_t pivot_mask = 0;
        std::vector<std::uint32_t> rows;
        rows.reserve(dimension);
        for (const auto pivot : selected_pivots) {
          pivot_mask |= std::uint32_t{1} << pivot;
          rows.push_back(std::uint32_t{1} << pivot);
        }
        std::vector<std::pair<std::uint32_t, std::uint32_t>> free_entries;
        for (std::uint32_t row = 0; row < dimension; ++row) {
          for (std::uint32_t column = selected_pivots[row] + 1;
               column < width; ++column) {
            if (((pivot_mask >> column) & 1U) == 0) {
              free_entries.emplace_back(row, column);
            }
          }
        }
        if (free_entries.size() >= 63) {
          throw std::invalid_argument("RREF subspace enumerator is too wide");
        }
        const auto assignments =
            std::uint64_t{1} << free_entries.size();
        for (std::uint64_t assignment = 0;
             assignment < assignments; ++assignment) {
          auto candidate = rows;
          for (std::size_t entry = 0; entry < free_entries.size(); ++entry) {
            if (((assignment >> entry) & 1U) == 0) continue;
            const auto [row, column] = free_entries[entry];
            candidate[row] |= std::uint32_t{1} << column;
          }
          visitor(candidate);
        }
      });
}

}  // namespace

std::string_view primitive_gram_sector_name(PrimitiveGramSector sector) {
  switch (sector) {
    case PrimitiveGramSector::zero:
      return "zero";
    case PrimitiveGramSector::rank_one:
      return "rank-one";
    case PrimitiveGramSector::full_alternating:
      return "full-alternating";
  }
  throw std::logic_error("unknown primitive Gram sector");
}

std::vector<PrimitiveGramSector> primitive_gram_sectors(std::uint32_t q) {
  if (q < 2 || q > 6) {
    throw std::invalid_argument("primitive tensor census supports q=2,...,6");
  }
  std::vector<PrimitiveGramSector> result{
      PrimitiveGramSector::zero, PrimitiveGramSector::rank_one};
  if ((q & 1U) == 0) {
    result.push_back(PrimitiveGramSector::full_alternating);
  }
  return result;
}

bool TensorRestrictionProfile::allows(
    std::uint32_t dimension,
    std::span<const std::uint64_t> canonical_key) const {
  if (dimension == 0 || dimension >= keys_by_dimension.size()) return false;
  return keys_by_dimension[dimension].contains(
      std::vector<std::uint64_t>(canonical_key.begin(), canonical_key.end()));
}

bool TensorRestrictionProfile::allows_standard_signature(
    std::uint32_t dimension,
    std::uint64_t signature) const {
  if (!uses_q7_primitive_family_recognizer) {
    throw std::logic_error(
        "standard-signature recognition is unavailable for this profile");
  }
  return is_q7_primitive_family_restriction(dimension, signature);
}

PrimitiveTensorSectorCensus census_primitive_tensor_sector(
    std::uint32_t q,
    PrimitiveGramSector sector) {
  const auto valid = primitive_gram_sectors(q);
  if (std::find(valid.begin(), valid.end(), sector) == valid.end()) {
    throw std::invalid_argument("Gram sector is not available at this q");
  }

  PrimitiveTensorSectorCensus result;
  result.logical_qubits = q;
  result.sector = sector;
  result.repeated_signature = repeated_signature(q, sector);
  result.alternating_dimension = cubic_tensor_alternating_term_count(q);
  result.tensor_count = std::uint64_t{1} << result.alternating_dimension;
  result.gram_stabilizer_order = gram_stabilizer_order(q, sector);

  std::uint32_t transvection_labels = 0;
  const auto basis_generators =
      gram_stabilizer_generators(q, sector, transvection_labels);
  result.generated_symplectic_transvection_labels = transvection_labels;
  result.action_generator_count =
      static_cast<std::uint32_t>(basis_generators.size());
  std::vector<AffineBitAction> actions;
  actions.reserve(basis_generators.size());
  for (const auto& basis : basis_generators) {
    actions.push_back(
        induce_alternating_action(q, result.repeated_signature, basis));
  }

  std::vector<std::uint8_t> visited(result.tensor_count, 0);
  std::vector<std::uint64_t> queue;
  std::map<std::uint64_t, std::pair<std::uint64_t, std::uint64_t>> histogram;
  const auto repeated = cubic_tensor_repeated_term_count(q);
  for (std::uint64_t representative = 0;
       representative < result.tensor_count; ++representative) {
    if (visited[representative]) continue;
    queue.clear();
    queue.push_back(representative);
    visited[representative] = 1;
    std::size_t head = 0;
    while (head < queue.size()) {
      const auto current = queue[head++];
      for (const auto& action : actions) {
        const auto image = action.apply(current);
        if (image >= result.tensor_count) {
          throw std::logic_error("Gram-stabilizer action left its state space");
        }
        if (!visited[image]) {
          visited[image] = 1;
          queue.push_back(image);
        }
      }
    }
    const auto orbit_size = static_cast<std::uint64_t>(queue.size());
    if (result.gram_stabilizer_order % orbit_size != 0) {
      throw std::logic_error("tensor orbit size does not divide group order");
    }
    ++result.tensor_orbit_count;
    const auto signature =
        result.repeated_signature | (representative << repeated);
    if (!cubic_tensor_radical_basis_from_signature(q, signature).empty()) {
      continue;
    }
    ++result.nondegenerate_tensor_orbit_count;
    result.nondegenerate_tensor_count += orbit_size;
    const auto hyperplanes =
        nondegenerate_tensor_hyperplane_count_from_signature(q, signature);
    auto& bucket = histogram[hyperplanes];
    ++bucket.first;
    bucket.second += orbit_size;
    if (hyperplanes != 0) continue;
    ++result.primitive_tensor_orbit_count;
    result.primitive_tensor_count += orbit_size;
    result.primitive_orbits.push_back(PrimitiveTensorOrbit{
        representative,
        signature,
        orbit_size,
        result.gram_stabilizer_order / orbit_size,
        hyperplanes,
    });
  }

  if (std::count(visited.begin(), visited.end(), std::uint8_t{1}) !=
      static_cast<std::ptrdiff_t>(result.tensor_count)) {
    throw std::logic_error("primitive tensor census did not cover its states");
  }
  for (const auto& [hyperplanes, counts] : histogram) {
    result.nondegenerate_hyperplane_histogram.push_back(
        TensorHyperplaneBucket{
            hyperplanes,
            counts.first,
            counts.second,
        });
  }
  return result;
}

Q7PrimitiveTensorAuthorityCensus
census_q7_primitive_tensors_from_authority(std::uint32_t workers) {
  if (workers == 0) {
    throw std::invalid_argument("q7 authority census requires workers >= 1");
  }
  constexpr std::uint32_t q = 7;
  Q7PrimitiveTensorAuthorityCensus result;
  result.general_linear_group_order = general_linear_group_order(q);
  result.alternating_tensor_count =
      std::uint64_t{1} << cubic_tensor_alternating_term_count(q);

  const auto& authority = q7_alternating_authority_entries();
  result.alternating_orbits.reserve(authority.size());
  std::vector<std::pair<std::uint32_t, std::uint64_t>> signatures;
  signatures.reserve(authority.size());
  for (const auto& entry : authority) {
    if (entry.stabilizer_order == 0 ||
        result.general_linear_group_order % entry.stabilizer_order != 0) {
      throw std::logic_error(
          "q7 alternating authority stabilizer does not divide GL(7,2)");
    }
    const auto signature =
        alternating_signature_from_triples(q, entry.triples);
    const auto orbit_size =
        result.general_linear_group_order / entry.stabilizer_order;
    result.alternating_orbit_mass += orbit_size;
    const auto radical =
        cubic_tensor_radical_basis_from_signature(q, signature);
    const auto hyperplanes = radical.empty()
                                 ? nondegenerate_tensor_hyperplane_count_from_signature(
                                       q, signature)
                                 : 0;
    result.alternating_orbits.push_back(Q7AlternatingFormOrbit{
        entry.id,
        signature,
        entry.stabilizer_order,
        orbit_size,
        static_cast<std::uint32_t>(radical.size()),
        hyperplanes,
    });
    signatures.emplace_back(entry.id, signature);
    if (radical.empty() && hyperplanes == 0) {
      result.primitive_orbits.push_back(Q7PrimitiveTensorOrbit{
          PrimitiveGramSector::zero,
          entry.id,
          0,
          signature,
          {},
          1,
      });
    }
  }
  if (result.alternating_orbit_mass != result.alternating_tensor_count) {
    throw std::logic_error(
        "q7 alternating authority does not cover all alternating forms");
  }

  struct RankOneCandidate {
    std::uint32_t authority_id;
    std::uint32_t functional;
    std::uint64_t signature;
  };
  std::vector<RankOneCandidate> rank_one_candidates;
  for (const auto& [authority_id, alternating_signature] : signatures) {
    for (std::uint32_t functional = 1;
         functional < (std::uint32_t{1} << q); ++functional) {
      ++result.rank_one_pairs_tested;
      const auto signature =
          alternating_signature ^ rank_one_cube_signature(q, functional);
      if (!cubic_tensor_radical_basis_from_signature(q, signature).empty()) {
        continue;
      }
      ++result.rank_one_nondegenerate_pairs;
      if (nondegenerate_tensor_hyperplane_count_from_signature(q, signature) !=
          0) {
        continue;
      }
      ++result.rank_one_primitive_pairs;
      rank_one_candidates.push_back(
          RankOneCandidate{authority_id, functional, signature});
    }
  }

  if (rank_one_candidates.size() == 1) {
    const auto& candidate = rank_one_candidates.front();
    result.primitive_orbits.push_back(Q7PrimitiveTensorOrbit{
        PrimitiveGramSector::rank_one,
        candidate.authority_id,
        candidate.functional,
        candidate.signature,
        {},
        1,
    });
    return result;
  }

  std::vector<std::vector<std::uint64_t>> canonical_keys(
      rank_one_candidates.size());
  std::atomic<std::size_t> next_candidate{0};
  const auto worker_count = std::min<std::size_t>(
      workers, std::max<std::size_t>(1, rank_one_candidates.size()));
  std::vector<std::thread> threads;
  threads.reserve(worker_count);
  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    threads.emplace_back([&] {
      while (true) {
        const auto index = next_candidate.fetch_add(1);
        if (index >= rank_one_candidates.size()) return;
        canonical_keys[index] =
            canonicalize_cubic_tensor_direct(
                q, rank_one_candidates[index].signature)
                .canonical_key_words;
      }
    });
  }
  for (auto& thread : threads) thread.join();

  std::map<std::vector<std::uint64_t>, std::size_t> rank_one_orbits;
  for (std::size_t index = 0; index < rank_one_candidates.size(); ++index) {
      const auto& candidate = rank_one_candidates[index];
      const auto [iterator, inserted] = rank_one_orbits.emplace(
          canonical_keys[index], result.primitive_orbits.size());
      if (inserted) {
        result.primitive_orbits.push_back(Q7PrimitiveTensorOrbit{
            PrimitiveGramSector::rank_one,
            candidate.authority_id,
            candidate.functional,
            candidate.signature,
            canonical_keys[index],
            1,
        });
      } else {
        ++result.primitive_orbits[iterator->second].represented_authority_pairs;
      }
  }
  return result;
}

Q8PrimitiveTensorAuthorityScreen
screen_q8_primitive_tensors_from_authority() {
  constexpr std::uint32_t q = 8;
  Q8PrimitiveTensorAuthorityScreen result;
  result.general_linear_group_order = general_linear_group_order(q);
  result.alternating_tensor_count = std::uint64_t{1} << 56;
  result.full_alternating_sector_theorem_excluded = true;

  const auto& authority = q8_alternating_authority_entries();
  result.alternating_orbits.reserve(authority.size());
  std::vector<std::pair<std::uint32_t, std::uint64_t>> states;
  states.reserve(authority.size());
  for (const auto& entry : authority) {
    if (entry.stabilizer_order == 0 ||
        result.general_linear_group_order % entry.stabilizer_order != 0) {
      throw std::logic_error(
          "q8 alternating authority stabilizer does not divide GL(8,2)");
    }
    const auto state = alternating_state_from_triples(q, entry.triples);
    const auto orbit_size =
        result.general_linear_group_order / entry.stabilizer_order;
    result.alternating_orbit_mass += orbit_size;
    const Q8TensorContractions contractions(
        PrimitiveGramSector::zero, state, 0);
    const auto radical_dimension =
        q8_tensor_radical_dimension(contractions);
    result.alternating_orbits.push_back(Q8AlternatingFormOrbit{
        entry.id,
        state,
        entry.stabilizer_order,
        orbit_size,
        radical_dimension,
    });
    states.emplace_back(entry.id, state);
    if (radical_dimension == 0 &&
        !q8_has_nondegenerate_hyperplane(contractions)) {
      result.primitive_candidates.push_back(Q8PrimitiveTensorCandidate{
          PrimitiveGramSector::zero,
          entry.id,
          0,
          state,
      });
    }
  }
  if (result.alternating_orbit_mass != result.alternating_tensor_count) {
    throw std::logic_error(
        "q8 alternating authority does not cover all alternating forms");
  }

  for (const auto& [authority_id, state] : states) {
    for (std::uint32_t functional = 1; functional < 256; ++functional) {
      ++result.rank_one_pairs_tested;
      const Q8TensorContractions contractions(
          PrimitiveGramSector::rank_one, state, functional);
      if (q8_tensor_radical_dimension(contractions) != 0) continue;
      ++result.rank_one_nondegenerate_pairs;
      if (q8_has_nondegenerate_hyperplane(contractions)) continue;
      ++result.rank_one_primitive_pairs;
      result.primitive_candidates.push_back(Q8PrimitiveTensorCandidate{
          PrimitiveGramSector::rank_one,
          authority_id,
          functional,
          state,
      });
    }
  }

  // The full-alternating proof leaves A=0 as its only possible candidate.
  // This explicit witness shows that candidate has a nondegenerate
  // hyperplane: x=e1+e3, y=e1, z=e2+e4 all lie as required in x^perp,
  // but T(x,y,z)=1, so x is not the restricted radical vector.
  const auto x = (std::uint32_t{1} << 0) | (std::uint32_t{1} << 2);
  const auto y = std::uint32_t{1} << 0;
  const auto z = (std::uint32_t{1} << 1) | (std::uint32_t{1} << 3);
  if (!q8_exceptional_tensor_value(
          PrimitiveGramSector::full_alternating, 0, 0, x, y, z)) {
    throw std::logic_error(
        "q8 full-alternating exclusion witness did not evaluate to one");
  }
  return result;
}

TensorRestrictionProfile build_tensor_restriction_profile(
    std::uint32_t q,
    std::span<const std::uint64_t> signatures,
    std::uint32_t minimum_dimension,
    std::uint32_t workers) {
  if (q < 1 || q > 7 || signatures.empty() || minimum_dimension == 0 ||
      minimum_dimension > q || workers == 0) {
    throw std::invalid_argument(
        "tensor restriction profile has invalid dimensions or no targets");
  }
  const auto term_count = cubic_tensor_term_count(q);
  const auto signature_limit = std::uint64_t{1} << term_count;
  for (const auto signature : signatures) {
    if (signature >= signature_limit) {
      throw std::invalid_argument(
          "target tensor signature lies outside its q space");
    }
  }

  TensorRestrictionProfile result;
  result.target_dimension = q;
  result.keys_by_dimension.resize(q + 1);
  auto maximum_profile_dimension = q;
  if (q == 7) {
    const auto authority = census_q7_primitive_tensors_from_authority();
    std::set<std::uint64_t> complete_primitive_signatures;
    for (const auto& orbit : authority.primitive_orbits) {
      complete_primitive_signatures.insert(orbit.standard_tensor_signature);
    }
    const std::set<std::uint64_t> supplied_signatures(
        signatures.begin(), signatures.end());
    if (supplied_signatures != complete_primitive_signatures) {
      throw std::invalid_argument(
          "q7 target profile requires the complete primitive authority set");
    }
    result.covers_all_primitive_orbits = true;
    result.uses_q7_primitive_family_recognizer = true;
    return result;
  }
  for (std::uint32_t dimension = minimum_dimension;
       dimension <= maximum_profile_dimension; ++dimension) {
    std::set<std::uint64_t> raw_restrictions;
    enumerate_rref_subspaces(
        q,
        dimension,
        [&](std::span<const std::uint32_t> basis) {
          for (const auto signature : signatures) {
            raw_restrictions.insert(
                restrict_cubic_tensor_signature(q, signature, basis));
          }
        });
    const std::vector<std::uint64_t> restrictions(
        raw_restrictions.begin(), raw_restrictions.end());
    std::vector<std::vector<std::uint64_t>> canonical_keys(
        restrictions.size());
    std::atomic<std::size_t> next_restriction{0};
    const auto worker_count = std::min<std::size_t>(
        workers, std::max<std::size_t>(1, restrictions.size()));
    std::vector<std::thread> threads;
    threads.reserve(worker_count);
    for (std::size_t worker = 0; worker < worker_count; ++worker) {
      threads.emplace_back([&] {
        while (true) {
          const auto index = next_restriction.fetch_add(1);
          if (index >= restrictions.size()) return;
          canonical_keys[index] =
              canonicalize_cubic_tensor_direct(
                  dimension, restrictions[index])
                  .canonical_key_words;
        }
      });
    }
    for (auto& thread : threads) thread.join();
    for (auto& key : canonical_keys) {
      result.keys_by_dimension[dimension].insert(std::move(key));
    }
  }
  if (q >= 2 && q <= 5) {
    std::set<std::vector<std::uint64_t>> primitive_keys;
    for (const auto sector : primitive_gram_sectors(q)) {
      const auto census = census_primitive_tensor_sector(q, sector);
      for (const auto& orbit : census.primitive_orbits) {
        primitive_keys.insert(
            canonicalize_cubic_tensor_direct(
                q, orbit.standard_tensor_signature)
                .canonical_key_words);
      }
    }
    result.covers_all_primitive_orbits =
        !primitive_keys.empty() &&
        result.keys_by_dimension[q] == primitive_keys;
  }
  return result;
}

bool has_nondegenerate_tensor_restriction(
    std::uint32_t q,
    std::uint64_t signature,
    std::uint32_t restriction_dimension) {
  if (q < 1 || q > 7 || restriction_dimension == 0 ||
      restriction_dimension > q) {
    throw std::invalid_argument(
        "nondegenerate tensor restriction has invalid dimensions");
  }
  const auto term_count = cubic_tensor_term_count(q);
  if (signature >= (std::uint64_t{1} << term_count)) {
    throw std::invalid_argument(
        "tensor signature lies outside its q space");
  }
  bool found = false;
  enumerate_rref_subspaces(
      q,
      restriction_dimension,
      [&](std::span<const std::uint32_t> basis) {
        if (found) return;
        const auto restricted = restrict_cubic_tensor_signature(
            q, signature, basis);
        found = cubic_tensor_radical_basis_from_signature(
                    restriction_dimension, restricted)
                    .empty();
      });
  return found;
}

}  // namespace utsp
